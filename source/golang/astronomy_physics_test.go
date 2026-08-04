package astronomy

import (
	"math"
	"testing"
)

func TestRefractionAngle(t *testing.T) {
	// Reference values from the Python implementation (identical algorithm).
	cases := []struct {
		alt      float64
		expected float64
	}{
		{-10, 0.581196063713},
		{-5, 0.617520817695},
		{-1, 0.646580620881},
		{0, 0.483032123074},
		{5, 0.161235443402},
		{10, 0.090128013386},
		{30, 0.029099750246},
		{60, 0.009752471283},
		{89, 0.000264253855},
	}
	for _, c := range cases {
		got := RefractionAngle(NormalRefraction, c.alt)
		if diff := absDiff(got, c.expected); diff > 1.0e-11 {
			t.Errorf("RefractionAngle(%g) = %0.12f, expected %0.12f (diff %g)", c.alt, got, c.expected, diff)
		}
	}
	// JplHorizons uses the same model for altitudes >= -1.
	if diff := absDiff(RefractionAngle(JplHorizonsRefraction, 0.0), 0.483032123074); diff > 1.0e-11 {
		t.Errorf("JPL refraction mismatch: %g", diff)
	}
	// Airless must give zero.
	if RefractionAngle(NoRefraction, 45.0) != 0.0 {
		t.Error("NoRefraction should give zero")
	}
	// Invalid altitudes are not corrected.
	if RefractionAngle(NormalRefraction, -91.0) != 0.0 || RefractionAngle(NormalRefraction, 91.0) != 0.0 {
		t.Error("invalid altitudes should not be corrected")
	}
}

func TestObserverGravity(t *testing.T) {
	cases := []struct {
		lat, height, expected float64
	}{
		{29.0, 0.0, 9.792473237342},
		{45.0, 1000.0, 9.803112946370},
		{90.0, 0.0, 9.832184937859},
	}
	for _, c := range cases {
		got := ObserverGravity(c.lat, c.height)
		if diff := absDiff(got, c.expected); diff > 1.0e-11 {
			t.Errorf("ObserverGravity(%g, %g) = %0.12f, expected %0.12f", c.lat, c.height, got, c.expected)
		}
	}
}

func TestObscuration(t *testing.T) {
	if got := obscuration(1.0, 1.0, 3.0); got != 0.0 {
		t.Errorf("non-overlapping discs should give 0, got %g", got)
	}
	if got := obscuration(1.0, 1.0, 0.0); got != 1.0 {
		t.Errorf("coincident equal discs should give 1, got %g", got)
	}
	if got := obscuration(2.0, 1.0, 0.0); math.Abs(got-0.25) > 1.0e-15 {
		t.Errorf("small disc inside big disc should give 0.25, got %g", got)
	}
	if got := obscuration(1.0, 2.0, 0.0); got != 1.0 {
		t.Errorf("big disc inside small disc should give 1, got %g", got)
	}
	// Tangent discs give zero overlap.
	if got := obscuration(1.0, 1.0, 2.0); got != 0.0 {
		t.Errorf("tangent discs should give 0, got %g", got)
	}
	// Half-overlap of two unit discs (centers 1 apart): each lens is two
	// 60-degree sectors minus triangles; total area = 2*(pi/3 - sqrt(3)/4)...
	// verify it is between 0 and 1 and symmetric.
	obs := obscuration(1.0, 1.0, 1.0)
	if obs <= 0.0 || obs >= 1.0 {
		t.Errorf("half-overlap should be in (0,1), got %g", obs)
	}
	// Invalid radii.
	if obscuration(0.0, 1.0, 0.5) != 0.0 {
		t.Error("invalid radius should give 0")
	}
}

func TestEclipseKindFromUmbra(t *testing.T) {
	if eclipseKindFromUmbra(0.02) != EclipseKindTotal {
		t.Error("positive umbra should be total")
	}
	if eclipseKindFromUmbra(-0.02) != EclipseKindAnnular {
		t.Error("negative umbra should be annular")
	}
}

func TestMaxAltitudeSlope(t *testing.T) {
	if _, err := maxAltitudeSlope(Earth, 0.0); err == nil {
		t.Error("Earth should not be allowed")
	}
	if _, err := maxAltitudeSlope(InvalidBody, 0.0); err == nil {
		t.Error("invalid body should fail")
	}
	moonSlope, err := maxAltitudeSlope(Moon, 0.0)
	if err != nil {
		t.Fatal(err)
	}
	// At the equator the Moon's slope is dominated by its RA rate of 4.5 deg/day.
	if diff := absDiff(moonSlope, 360.0/solarDaysPerSiderealDay-4.5); diff > 1.0e-12 {
		t.Errorf("Moon equator slope = %g, expected %g", moonSlope, 360.0/solarDaysPerSiderealDay-4.5)
	}
}

func TestFindAscent(t *testing.T) {
	// Mock altitude-difference function: rises from -2 to +2 across t=0..4.
	// f(t) = t - 2.
	altDiff := func(time AstroTime) (float64, error) {
		return time.Ut - 2.0, nil
	}
	t1 := TimeFromUniversalDays(0.0)
	t2 := TimeFromUniversalDays(4.0)
	ascent, ok, err := findAscent(0, 1.0, t1, t2, -2.0, 2.0, altDiff)
	if err != nil {
		t.Fatal(err)
	}
	if !ok {
		t.Fatal("expected to find an ascent")
	}
	if ascent.tx.Ut != 0.0 || ascent.ty.Ut != 4.0 {
		t.Errorf("ascent interval = [%g, %g], expected [0, 4]", ascent.tx.Ut, ascent.ty.Ut)
	}

	// A purely descending function must not produce an ascent.
	_, ok, err = findAscent(0, 1.0, t1, t2, +2.0, -2.0, altDiff)
	if err != nil {
		t.Fatal(err)
	}
	if ok {
		t.Error("descending function should not produce an ascent")
	}
}

func TestMoonMagnitude(t *testing.T) {
	// Reference: Python Illumination(Moon, 2019-06-24T15:45:37Z).
	got := moonMagnitude(81.608458773575, 1.016824802299, 0.002696881050)
	if diff := absDiff(got, -10.279489956918); diff > 1.0e-9 {
		t.Errorf("Moon magnitude = %g, expected -10.279489956918 (diff %g)", got, diff)
	}
}

func TestVisualMagnitude(t *testing.T) {
	// Reference: Python Illumination for Venus and Jupiter.
	venus, err := visualMagnitude(Venus, 19.684808500163, 0.721832629895, 1.666572915237)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(venus, -3.842863697359); diff > 1.0e-9 {
		t.Errorf("Venus magnitude = %g, expected -3.842863697359", venus)
	}
	jupiter, err := visualMagnitude(Jupiter, 2.874270604969, 5.294820010652, 4.307013639930)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(jupiter, -2.595491403719); diff > 1.0e-9 {
		t.Errorf("Jupiter magnitude = %g, expected -2.595491403719", jupiter)
	}
	if _, err := visualMagnitude(Earth, 10.0, 1.0, 1.0); err == nil {
		t.Error("Earth should not have a visual magnitude")
	}
}

func TestCalcShadow(t *testing.T) {
	time := TimeFromUniversalDays(0.0)
	// Target on the shadow axis: no lateral offset.
	target := AstroVector{X: 0.0, Y: 0.0, Z: 1.0, T: time}
	dir := AstroVector{X: 0.0, Y: 0.0, Z: 1.0, T: time}
	s := calcShadow(6378.1366, time, target, dir)
	if s.r != 0.0 {
		t.Errorf("on-axis target should have r = 0, got %g", s.r)
	}
	if math.Abs(s.u-1.0) > 1.0e-15 {
		t.Errorf("u should be 1.0, got %g", s.u)
	}

	// Target offset by 1 AU perpendicular to the axis.
	target2 := AstroVector{X: 1.0, Y: 0.0, Z: 1.0, T: time}
	s2 := calcShadow(6378.1366, time, target2, dir)
	if diff := absDiff(s2.r, KmPerAu); diff > 1.0e-6 {
		t.Errorf("perpendicular offset should give r = 1 AU in km, got %g", s2.r)
	}
}

func TestJupiterMoonElem2pv(t *testing.T) {
	// A nearly circular orbit in the equatorial plane.
	mu := 2.8248942843381399e-07
	elem := [6]float64{0.00282, 0.0, 0.0, 0.0, 0.0, 0.0}
	time := TimeFromUniversalDays(0.0)
	state := jupiterMoonElem2pv(time, mu, elem)
	// Position should start at (a, 0, 0) with velocity (0, n*a, 0).
	n := math.Sqrt(mu / (elem[0] * elem[0] * elem[0]))
	if diff := absDiff(state.X, elem[0]); diff > 1.0e-15 {
		t.Errorf("X = %g, expected %g", state.X, elem[0])
	}
	if diff := absDiff(state.Y, 0.0); diff > 1.0e-15 {
		t.Errorf("Y = %g, expected 0", state.Y)
	}
	if diff := absDiff(state.Vx, 0.0); diff > 1.0e-15 {
		t.Errorf("Vx = %g, expected 0", state.Vx)
	}
	if diff := absDiff(state.Vy, n*elem[0]); diff > 1.0e-15 {
		t.Errorf("Vy = %g, expected %g", state.Vy, n*elem[0])
	}
}
