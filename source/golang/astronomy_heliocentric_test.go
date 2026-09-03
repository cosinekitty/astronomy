package astronomy

import (
	"testing"
)

func TestHelioVector(t *testing.T) {
	// Reference values from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	// Pluto (via the gravity simulator).
	pluto, err := HelioVector(Pluto, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(pluto.X, 1.2407758350790576e+01); diff > 1.0e-12 {
		t.Errorf("Pluto x error = %g", diff)
	}
	// Sun is at the origin.
	sun, err := HelioVector(Sun, time)
	if err != nil {
		t.Fatal(err)
	}
	if sun.Length() != 0.0 {
		t.Error("Sun should be at origin")
	}
	// Moon is Earth + geocentric Moon.
	moon, err := HelioVector(Moon, time)
	if err != nil {
		t.Fatal(err)
	}
	earth := calcEarth(time)
	gm := GeoMoon(time)
	if diff := absDiff(moon.X, earth.X+gm.X); diff > 1.0e-15 {
		t.Errorf("Moon x mismatch: %g", diff)
	}
}

func TestHelioDistance(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	d, err := HelioDistance(Jupiter, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(d, 5.2948200106504437); diff > 1.0e-12 {
		t.Errorf("Jupiter helio distance error = %g", diff)
	}
	d, err = HelioDistance(Sun, time)
	if err != nil {
		t.Fatal(err)
	}
	if d != 0.0 {
		t.Errorf("Sun helio distance should be 0, got %g", d)
	}
}

func TestBaryState(t *testing.T) {
	// Reference values from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	state, err := BaryState(Mars, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(state.X, -1.0324815173507131e+00); diff > 1.0e-12 {
		t.Errorf("bary Mars x error = %g", diff)
	}
	if diff := absDiff(state.Y, 1.1641629491140129e+00); diff > 1.0e-12 {
		t.Errorf("bary Mars y error = %g", diff)
	}
	if diff := absDiff(state.Z, 5.6158377804275639e-01); diff > 1.0e-12 {
		t.Errorf("bary Mars z error = %g", diff)
	}
	if diff := absDiff(state.Vx, -1.0395333634235157e-02); diff > 1.0e-15 {
		t.Errorf("bary Mars vx error = %g", diff)
	}
	// SSB is at the origin.
	ssb, err := BaryState(Ssb, time)
	if err != nil {
		t.Fatal(err)
	}
	if ssb.X != 0.0 || ssb.Y != 0.0 || ssb.Z != 0.0 {
		t.Error("SSB should be at origin in barycentric coordinates")
	}
}

func TestHelioStateMoon(t *testing.T) {
	// Reference values from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	state, err := HelioState(Moon, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(state.X, 4.8676682902330903e-02); diff > 1.0e-12 {
		t.Errorf("helio Moon x error = %g", diff)
	}
	if diff := absDiff(state.Y, -9.3177576090776471e-01); diff > 1.0e-12 {
		t.Errorf("helio Moon y error = %g", diff)
	}
	if diff := absDiff(state.Vx, 1.6948727466270977e-02); diff > 1.0e-12 {
		t.Errorf("helio Moon vx error = %g", diff)
	}
}

func TestEclipticLongitude(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	lon, err := EclipticLongitude(Mars, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(lon, 129.014306945392); diff > 1.0e-9 {
		t.Errorf("Mars ecliptic longitude = %0.12f, expected 129.014306945392", lon)
	}
	if _, err := EclipticLongitude(Sun, time); err == nil {
		t.Error("expected error for Sun ecliptic longitude")
	}
}

func TestIllumination(t *testing.T) {
	// Reference values from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	illum, err := Illumination(Mars, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(illum.Mag, 1.805704385814); diff > 1.0e-9 {
		t.Errorf("Mars mag = %0.12f, expected 1.805704385814", illum.Mag)
	}
	if diff := absDiff(illum.PhaseAngle, 13.672425351972); diff > 1.0e-9 {
		t.Errorf("Mars phase = %0.12f, expected 13.672425351972", illum.PhaseAngle)
	}
	if diff := absDiff(illum.PhaseFraction, 0.985831495097); diff > 1.0e-9 {
		t.Errorf("Mars phase fraction = %0.12f, expected 0.985831495097", illum.PhaseFraction)
	}
	if diff := absDiff(illum.HelioDist, 1.646973604252); diff > 1.0e-9 {
		t.Errorf("Mars helio dist = %0.12f, expected 1.646973604252", illum.HelioDist)
	}

	sat, err := Illumination(Saturn, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(sat.Mag, 0.005171536245); diff > 1.0e-9 {
		t.Errorf("Saturn mag = %0.12f, expected 0.005171536245", sat.Mag)
	}
	if diff := absDiff(sat.RingTilt, -24.086756644497); diff > 1.0e-9 {
		t.Errorf("Saturn ring tilt = %0.12f, expected -24.086756644497", sat.RingTilt)
	}
	if diff := absDiff(sat.HelioDist, 10.051007587361); diff > 1.0e-9 {
		t.Errorf("Saturn helio dist = %0.12f, expected 10.051007587361", sat.HelioDist)
	}

	if _, err := Illumination(Earth, time); err == nil {
		t.Error("expected error for Earth illumination")
	}
}

func TestRotation_HOR_ECL(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	a := Rotation_HOR_ECL(time, observer)
	b := Rotation_ECL_HOR(time, observer)
	rt := CombineRotation(a, b)
	ident := IdentityMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(rt.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-14 {
				t.Errorf("HOR_ECL roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}
