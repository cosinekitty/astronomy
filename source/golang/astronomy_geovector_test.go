package astronomy

import (
	"testing"
)

func TestGeoVector(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	gv, err := GeoVector(Mars, time, WithAberration)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(gv.X, -1.0758899615720006e+00); diff > 1.0e-12 {
		t.Errorf("GeoVector Mars x error = %g", diff)
	}
	if diff := absDiff(gv.Y, 2.0888980769720087e+00); diff > 1.0e-12 {
		t.Errorf("GeoVector Mars y error = %g", diff)
	}
	if diff := absDiff(gv.Z, 9.6247792683541467e-01); diff > 1.0e-12 {
		t.Errorf("GeoVector Mars z error = %g", diff)
	}

	// Earth is at the origin; Moon is geocentric.
	earth, err := GeoVector(Earth, time, WithAberration)
	if err != nil {
		t.Fatal(err)
	}
	if earth.Length() != 0.0 {
		t.Error("GeoVector(Earth) should be zero")
	}
	moon, err := GeoVector(Moon, time, WithAberration)
	if err != nil {
		t.Fatal(err)
	}
	gm := GeoMoon(time)
	if diff := absDiff(moon.X, gm.X); diff > 1.0e-15 {
		t.Errorf("GeoVector(Moon) x mismatch: %g", diff)
	}
}

func TestBackdatePositionStar(t *testing.T) {
	// Define a star and verify its backdated position works.
	if err := DefineStar(Star1, 3.5, 45.0, 100.0); err != nil {
		t.Fatal(err)
	}
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	pos, err := BackdatePosition(time, Earth, Star1, WithAberration)
	if err != nil {
		t.Fatal(err)
	}
	if !isfinite(pos.X) || !isfinite(pos.Y) || !isfinite(pos.Z) {
		t.Error("star backdated position is not finite")
	}
}

func TestLagrangePoint(t *testing.T) {
	// Reference from the C implementation: Sun-Earth L4.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	lp, err := LagrangePoint(4, time, Sun, Earth)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(lp.X, 9.0235761392391545e-01); diff > 1.0e-12 {
		t.Errorf("L4 x error = %g", diff)
	}
	if diff := absDiff(lp.Y, -4.2925741091267255e-01); diff > 1.0e-12 {
		t.Errorf("L4 y error = %g", diff)
	}
	if diff := absDiff(lp.Z, -1.8608699723263533e-01); diff > 1.0e-12 {
		t.Errorf("L4 z error = %g", diff)
	}
	if _, err := LagrangePoint(1, time, InvalidBody, Earth); err == nil {
		t.Error("expected error for invalid major body")
	}
}

func TestSearchRelativeLongitude(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	root, err := SearchRelativeLongitude(Mars, 45.0, time)
	if err != nil {
		t.Fatal(err)
	}
	want := TimeFromCalendar(2021, 1, 21, 13, 30, 48.307)
	if diff := absDiff(root.Ut, want.Ut); diff > 1.0e-5 {
		t.Errorf("relative longitude time error = %g days", diff)
	}
	if _, err := SearchRelativeLongitude(Earth, 45.0, time); err == nil {
		t.Error("expected error for Earth")
	}
	if _, err := SearchRelativeLongitude(Sun, 45.0, time); err == nil {
		t.Error("expected error for Sun")
	}
}
