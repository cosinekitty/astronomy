package astronomy

import (
	"math"
	"testing"
)

func TestHorizon(t *testing.T) {
	// Reference values from the Python implementation.
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)

	hor := Horizon(time, observer, 18.0, 30.0, NoRefraction)
	if diff := absDiff(hor.Azimuth, 338.490014795096); diff > 1.0e-9 {
		t.Errorf("azimuth = %0.12f, expected 338.490014795096", hor.Azimuth)
	}
	if diff := absDiff(hor.Altitude, -27.350949611085); diff > 1.0e-9 {
		t.Errorf("altitude = %0.12f, expected -27.350949611085", hor.Altitude)
	}

	hor = Horizon(time, observer, 18.0, 30.0, NormalRefraction)
	if diff := absDiff(hor.Azimuth, 338.490014795096); diff > 1.0e-9 {
		t.Errorf("refracted azimuth = %0.12f", hor.Azimuth)
	}
	if diff := absDiff(hor.Altitude, -26.895807342568); diff > 1.0e-9 {
		t.Errorf("refracted altitude = %0.12f, expected -26.895807342568", hor.Altitude)
	}
	if diff := absDiff(hor.Ra, 18.013029960213); diff > 1.0e-9 {
		t.Errorf("refracted ra = %0.12f", hor.Ra)
	}
	if diff := absDiff(hor.Dec, 30.422642492188); diff > 1.0e-9 {
		t.Errorf("refracted dec = %0.12f", hor.Dec)
	}
}

func TestLibration(t *testing.T) {
	// Reference values from the Python implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	lib := Libration(time)
	if diff := absDiff(lib.Elon, -2.747257046332); diff > 1.0e-9 {
		t.Errorf("elon = %0.12f, expected -2.747257046332", lib.Elon)
	}
	if diff := absDiff(lib.Elat, 6.311759826781); diff > 1.0e-9 {
		t.Errorf("elat = %0.12f, expected 6.311759826781", lib.Elat)
	}
	if diff := absDiff(lib.Mlon, 354.595129819365); diff > 1.0e-9 {
		t.Errorf("mlon = %0.12f, expected 354.595129819365", lib.Mlon)
	}
	if diff := absDiff(lib.Mlat, -4.851798346972); diff > 1.0e-9 {
		t.Errorf("mlat = %0.12f, expected -4.851798346972", lib.Mlat)
	}
	if diff := absDiff(lib.DistKm, 403447.662575781520); diff > 1.0e-6 {
		t.Errorf("dist_km = %0.12f, expected 403447.662575781520", lib.DistKm)
	}
	if diff := absDiff(lib.DiamDeg, 0.493476622822); diff > 1.0e-9 {
		t.Errorf("diam_deg = %0.12f, expected 0.493476622822", lib.DiamDeg)
	}
}

func TestMoonEclipticLatitudeDegrees(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	lat := MoonEclipticLatitudeDegrees(time)
	// Cross-check with GeoMoon: the ecliptic latitude is close to the Moon's
	// declination rotated back to the ecliptic; just check it is a small angle.
	if math.Abs(lat) > 20.0 {
		t.Errorf("implausible Moon ecliptic latitude = %g", lat)
	}
}

func TestTerraRoundTrip(t *testing.T) {
	// Convert an observer to a geocentric vector and back.
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	st := 12.5
	pos, _ := terra(observer, st)
	obs2 := inverseTerra(pos, st)
	if diff := absDiff(obs2.Latitude, observer.Latitude); diff > 1.0e-6 {
		t.Errorf("latitude round-trip error = %g", diff)
	}
	if diff := absDiff(obs2.Longitude, observer.Longitude); diff > 1.0e-6 {
		t.Errorf("longitude round-trip error = %g", diff)
	}
	if diff := absDiff(obs2.Height, observer.Height); diff > 1.0e-4 {
		t.Errorf("height round-trip error = %g", diff)
	}
}

func TestNutationRot(t *testing.T) {
	// Reference: Python reports dpsi = -16.297 arcsec, deps = -3.3007 arcsec
	// at this time; the rotation matrix derived from those angles should be
	// almost an identity (small angles).
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	tilt := etilt(&time)
	if diff := absDiff(tilt.Dpsi, -16.297019219816); diff > 1.0e-9 {
		t.Errorf("dpsi = %0.12f, expected -16.297019219816", tilt.Dpsi)
	}
	if diff := absDiff(tilt.Deps, -3.300723406454); diff > 1.0e-9 {
		t.Errorf("deps = %0.12f, expected -3.300723406454", tilt.Deps)
	}
	rot := nutationRot(&time, from2000)
	ident := IdentityMatrix()
	roundTrip := CombineRotation(rot, InverseRotation(rot))
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(roundTrip.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("nutation round-trip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}

func TestRotation_EQD_HOR(t *testing.T) {
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	rot := Rotation_EQD_HOR(time, observer)
	ident := IdentityMatrix()
	// A rotation matrix must have orthonormal columns.
	for j := 0; j < 3; j++ {
		col := [3]float64{rot.Rot[0][j], rot.Rot[1][j], rot.Rot[2][j]}
		norm := math.Sqrt(col[0]*col[0] + col[1]*col[1] + col[2]*col[2])
		if diff := absDiff(norm, 1.0); diff > 1.0e-12 {
			t.Errorf("column %d not unit length: %g", j, norm)
		}
	}
	// The zenith vector in EQD should rotate to (0,0,1) in HOR.
	eqdZenith := AstroVector{X: rot.Rot[0][2], Y: rot.Rot[1][2], Z: rot.Rot[2][2], T: time}
	horZenith := RotateVector(rot, eqdZenith)
	if diff := absDiff(horZenith.X, 0.0); diff > 1.0e-12 {
		t.Errorf("zenith x = %g", horZenith.X)
	}
	if diff := absDiff(horZenith.Z, 1.0); diff > 1.0e-12 {
		t.Errorf("zenith z = %g", horZenith.Z)
	}
	_ = ident
}

func TestSearch(t *testing.T) {
	// f(t) = t - 2 crosses zero at t = 2 (UT days since J2000).
	f := func(t AstroTime) (float64, error) {
		return t.Ut - 2.0, nil
	}
	t1 := TimeFromUniversalDays(0.0)
	t2 := TimeFromUniversalDays(4.0)
	root, err := Search(f, t1, t2, 1.0e-6)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(root.Ut, 2.0); diff > 1.0e-6 {
		t.Errorf("search root error = %g", diff)
	}

	// A function that never crosses zero must fail.
	g := func(t AstroTime) (float64, error) {
		return 1.0, nil
	}
	if _, err := Search(g, t1, t2, 1.0e-6); err == nil {
		t.Error("expected search failure for constant positive function")
	}
}

func TestInverseRefractionRoundTrip(t *testing.T) {
	for alt := -89.0; alt <= 89.0; alt += 10.0 {
		refr := RefractionAngle(NormalRefraction, alt)
		corrected := alt + refr
		inv := InverseRefractionAngle(NormalRefraction, corrected)
		check := corrected + inv
		if diff := absDiff(check, alt); diff > 2.0e-14 {
			t.Errorf("refraction round-trip error at alt=%g: %g", alt, diff)
		}
	}
}

func TestHorizonFromVectorRoundTrip(t *testing.T) {
	time := TimeFromUniversalDays(0.0)
	// A vector straight up in the horizontal system.
	vec := AstroVector{X: 0.0, Y: 0.0, Z: 1.0, T: time}
	sphere, err := HorizonFromVector(vec, NoRefraction)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(sphere.Lat, 90.0); diff > 1.0e-12 {
		t.Errorf("zenith altitude = %g", sphere.Lat)
	}
	// Vector toward the north.
	vec = AstroVector{X: 1.0, Y: 0.0, Z: 0.0, T: time}
	sphere, err = HorizonFromVector(vec, NoRefraction)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(sphere.Lon, 0.0); diff > 1.0e-12 {
		t.Errorf("north azimuth = %g", sphere.Lon)
	}
	// VectorFromHorizon must invert HorizonFromVector.
	back := VectorFromHorizon(sphere, time, NoRefraction)
	if diff := absDiff(back.X, 1.0); diff > 1.0e-12 {
		t.Errorf("back X = %g", back.X)
	}
}
