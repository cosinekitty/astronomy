package astronomy

import (
	"math"
	"testing"
)

func TestEclipticGeoMoon(t *testing.T) {
	// Reference values from the Python implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	s := EclipticGeoMoon(time)
	if diff := absDiff(s.Lat, -4.851798346972); diff > 1.0e-11 {
		t.Errorf("ecliptic latitude = %0.12f, expected -4.851798346972", s.Lat)
	}
	if diff := absDiff(s.Lon, 354.590602869581); diff > 1.0e-9 {
		t.Errorf("ecliptic longitude = %0.12f, expected 354.590602869581", s.Lon)
	}
	if diff := absDiff(s.Dist, 0.002696881050); diff > 1.0e-12 {
		t.Errorf("ecliptic distance = %0.12f, expected 0.002696881050", s.Dist)
	}
}

func TestSeasons(t *testing.T) {
	// Reference values from the Python implementation.
	sea, err := Seasons(2023)
	if err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		label    string
		got      AstroTime
		expected string
	}{
		{"mar", sea.MarEquinox, "2023-03-20T21:24:23.770Z"},
		{"jun", sea.JunSolstice, "2023-06-21T14:57:22.411Z"},
		{"sep", sea.SepEquinox, "2023-09-23T06:49:59.173Z"},
		{"dec", sea.DecSolstice, "2023-12-22T03:27:34.200Z"},
	}
	for _, c := range cases {
		want := timeFromIso(c.expected, t)
		if diff := absDiff(c.got.Ut, want.Ut); diff > 1.2e-7 {
			t.Errorf("seasons %s = %0.9f, expected %s (diff %g days)", c.label, c.got.Ut, c.expected, diff)
		}
	}
}

func TestConstellation(t *testing.T) {
	// Reference values from the Python implementation.
	c, err := Constellation(15.0, 25.0)
	if err != nil {
		t.Fatal(err)
	}
	if c.Symbol != "Boo" {
		t.Errorf("symbol = %s, expected Boo", c.Symbol)
	}
	if c.Name != "Bootes" {
		t.Errorf("name = %s, expected Bootes", c.Name)
	}
	if diff := absDiff(c.Ra1875, 14.908526912852); diff > 1.0e-9 {
		t.Errorf("ra1875 = %0.12f, expected 14.908526912852", c.Ra1875)
	}
	if diff := absDiff(c.Dec1875, 25.497006480218); diff > 1.0e-9 {
		t.Errorf("dec1875 = %0.12f, expected 25.497006480218", c.Dec1875)
	}
	if _, err := Constellation(0.0, 91.0); err == nil {
		t.Error("expected error for invalid declination")
	}
}

func TestRotationAxis(t *testing.T) {
	// Reference values from the Python implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	ax, err := RotationAxis(Mars, &time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(ax.Ra, 21.177338213769); diff > 1.0e-9 {
		t.Errorf("Mars RA = %0.12f, expected 21.177338213769", ax.Ra)
	}
	if diff := absDiff(ax.Dec, 52.874406829529); diff > 1.0e-9 {
		t.Errorf("Mars Dec = %0.12f, expected 52.874406829529", ax.Dec)
	}
	if diff := absDiff(ax.Spin, 2496477.461263322271); diff > 1.0e-6 {
		t.Errorf("Mars spin = %0.12f, expected 2496477.461263322271", ax.Spin)
	}
	earth, err := RotationAxis(Earth, &time)
	if err != nil {
		t.Fatal(err)
	}
	// The Earth's north pole in J2000 equatorial coordinates is near
	// RA = 23.96 hours, Dec = 89.89 degrees (reference: Python).
	if diff := absDiff(earth.Dec, 89.893368008); diff > 1.0e-6 {
		t.Errorf("Earth axis Dec = %0.12f, expected near 89.89", earth.Dec)
	}
	if diff := absDiff(earth.Ra, 23.958926161); diff > 1.0e-6 {
		t.Errorf("Earth axis RA = %0.12f, expected near 23.96", earth.Ra)
	}
}

func TestGeoMoonState(t *testing.T) {
	// Reference values from the Python implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	s := GeoMoonState(time)
	if diff := absDiff(s.X, 0.002674037027); diff > 1.0e-11 {
		t.Errorf("x = %0.12f, expected 0.002674037027", s.X)
	}
	if diff := absDiff(s.Y, -0.000153161032); diff > 1.0e-11 {
		t.Errorf("y = %0.12f, expected -0.000153161032", s.Y)
	}
	if diff := absDiff(s.Vx, 0.000042675762); diff > 1.0e-12 {
		t.Errorf("vx = %0.12f, expected 0.000042675762", s.Vx)
	}
	if diff := absDiff(s.Vy, 0.000519190910); diff > 1.0e-12 {
		t.Errorf("vy = %0.12f, expected 0.000519190910", s.Vy)
	}
	if diff := absDiff(s.Vz, 0.000204155399); diff > 1.0e-12 {
		t.Errorf("vz = %0.12f, expected 0.000204155399", s.Vz)
	}
}

func TestGeoEmbState(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	s := GeoEmbState(time)
	// EMB is the Moon position divided by (1 + mass ratio).
	d := 1.0 + earthMoonMassRatio
	if diff := absDiff(s.X, 0.002674037027/d); diff > 1.0e-11 {
		t.Errorf("EMB x error = %g", diff)
	}
	if diff := absDiff(s.Y, -0.000153161032/d); diff > 1.0e-11 {
		t.Errorf("EMB y error = %g", diff)
	}
}

func TestEquatorialToEcliptic(t *testing.T) {
	// EQJ vector of the north ecliptic pole should map to elat=90 in ECT.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	// The ecliptic north pole in EQJ coordinates is (0, sin(23.44), cos(23.44)).
	obliquity := RadiansFromDegrees(23.4392911)
	// The ecliptic north pole in EQJ coordinates has RA = 18h (270 degrees),
	// so its y component is negative.
	pole := AstroVector{X: 0.0, Y: -math.Sin(obliquity), Z: math.Cos(obliquity), T: time}
	ecl := EquatorialToEcliptic(pole)
	if diff := absDiff(ecl.Elat, 90.0); diff > 0.02 {
		t.Errorf("ecliptic pole elat = %0.9f, expected 90", ecl.Elat)
	}

	// Reference from the C implementation: Astronomy_Ecliptic((1,0,0)) gives
	// elon = 0.267563546075, elat = 0.000229232305.
	ecl = EquatorialToEcliptic(AstroVector{X: 1.0, Y: 0.0, Z: 0.0, T: time})
	if diff := absDiff(ecl.Elon, 0.267563546075); diff > 1.0e-9 {
		t.Errorf("ecliptic elon = %0.12f, expected 0.267563546075", ecl.Elon)
	}
	if diff := absDiff(ecl.Elat, 0.000229232305); diff > 1.0e-9 {
		t.Errorf("ecliptic elat = %0.12f, expected 0.000229232305", ecl.Elat)
	}
}

func TestSunPosition(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	sp := SunPosition(time)
	if diff := absDiff(sp.Elon, 92.856080675903); diff > 1.0e-9 {
		t.Errorf("Sun ecliptic longitude = %0.12f, expected 92.856080675903", sp.Elon)
	}
	if diff := absDiff(sp.Elat, 0.000040994467); diff > 1.0e-9 {
		t.Errorf("Sun ecliptic latitude = %0.12f, expected 0.000040994467", sp.Elat)
	}
}

func TestRotation_EQJ_ECT(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	eqjEct := Rotation_EQJ_ECT(time)
	ectEqj := Rotation_ECT_EQJ(time)
	rt := CombineRotation(eqjEct, ectEqj)
	ident := IdentityMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(rt.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-14 {
				t.Errorf("EQJ_ECT roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}
