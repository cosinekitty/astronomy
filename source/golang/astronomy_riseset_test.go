package astronomy

import (
	"math"
	"testing"
)

func TestEquator(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	eq, err := Equator(Mars, &time, observer, EquatorJ2000, WithAberration)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(eq.Ra, 7.816766033729); diff > 1.0e-9 {
		t.Errorf("Mars RA = %0.12f, expected 7.816766033729", eq.Ra)
	}
	if diff := absDiff(eq.Dec, 22.274737976068); diff > 1.0e-9 {
		t.Errorf("Mars Dec = %0.12f, expected 22.274737976068", eq.Dec)
	}
	if diff := absDiff(eq.Dist, 2.539142577487); diff > 1.0e-9 {
		t.Errorf("Mars dist = %0.12f, expected 2.539142577487", eq.Dist)
	}
}

func TestAngleFromSunAndElongation(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	angle, err := AngleFromSun(Mars, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(angle, 22.515990265309); diff > 1.0e-9 {
		t.Errorf("Mars angle from Sun = %0.12f, expected 22.515990265309", angle)
	}
	elong, err := Elongation(Mars, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(elong.Elongation, 22.515990265309); diff > 1.0e-9 {
		t.Errorf("Mars elongation = %0.12f", elong.Elongation)
	}
	if diff := absDiff(elong.EclipticSeparation, 22.486283399513); diff > 1.0e-6 {
		t.Errorf("Mars ecliptic separation = %0.12f (diff %g)", elong.EclipticSeparation, diff)
	}
	if elong.Visibility != VisibilityEvening {
		t.Errorf("Mars visibility = %d, expected evening", elong.Visibility)
	}
}

func TestPairLongitudeAndMoonPhase(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	pl, err := PairLongitude(Mars, Sun, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(pl, 22.486283399513); diff > 1.0e-6 {
		t.Errorf("pair longitude = %0.12f, expected 22.486283399513 (diff %g)", pl, diff)
	}
	phase, err := MoonPhase(time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(phase, 261.729013321958); diff > 1.0e-6 {
		t.Errorf("moon phase = %0.12f, expected 261.729013321958 (diff %g)", phase, diff)
	}
}

func TestHourAngle(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	ha, err := HourAngle(Mars, &time, observer)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(ha, 20.691665395070); diff > 1.0e-9 {
		t.Errorf("Mars hour angle = %0.12f, expected 20.691665395070", ha)
	}
}

func TestSearchMoonQuarter(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	mq, err := SearchMoonQuarter(time)
	if err != nil {
		t.Fatal(err)
	}
	if mq.Quarter != 3 {
		t.Errorf("moon quarter = %d, expected 3", mq.Quarter)
	}
	want := TimeFromCalendar(2019, 6, 25, 9, 47, 5.867)
	if diff := absDiff(mq.Time.Ut, want.Ut); diff > 0.01 {
		t.Errorf("moon quarter time error = %g days", diff)
	}
	next, err := NextMoonQuarter(mq)
	if err != nil {
		t.Fatal(err)
	}
	if next.Quarter != 0 {
		t.Errorf("next moon quarter = %d, expected 0", next.Quarter)
	}
	wantNext := TimeFromCalendar(2019, 7, 2, 19, 16, 46.600)
	if diff := absDiff(next.Time.Ut, wantNext.Ut); diff > 0.01 {
		t.Errorf("next moon quarter time error = %g days", diff)
	}
}

func TestSearchPlanetApsis(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	apsis, err := SearchPlanetApsis(Mars, time)
	if err != nil {
		t.Fatal(err)
	}
	if apsis.Kind != ApsisKindApocenter {
		t.Errorf("expected aphelion, got kind %d", apsis.Kind)
	}
	want := TimeFromCalendar(2019, 8, 26, 1, 1, 26.732)
	if diff := absDiff(apsis.Time.Ut, want.Ut); diff > 1.0e-5 {
		t.Errorf("Mars aphelion time error = %g days", diff)
	}
	if diff := absDiff(apsis.DistAu, 1.666072695); diff > 1.0e-6 {
		t.Errorf("Mars aphelion distance error = %g", diff)
	}
	next, err := NextPlanetApsis(Mars, apsis)
	if err != nil {
		t.Fatal(err)
	}
	if next.Kind != ApsisKindPericenter {
		t.Errorf("expected perihelion after aphelion, got %d", next.Kind)
	}
	wantNext := TimeFromCalendar(2020, 8, 3, 9, 50, 49.517)
	if diff := absDiff(next.Time.Ut, wantNext.Ut); diff > 1.0e-5 {
		t.Errorf("Mars perihelion time error = %g days", diff)
	}
}

func TestSearchRiseSet(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	sunrise, err := SearchRiseSet(Sun, observer, +1, time, 1.0, 0.0)
	if err != nil {
		t.Fatal(err)
	}
	want := TimeFromCalendar(2019, 6, 25, 10, 26, 51.668)
	if diff := absDiff(sunrise.Ut, want.Ut); diff > 0.001 {
		t.Errorf("sunrise error = %g days", diff)
	}
	moonrise, err := SearchRiseSet(Moon, observer, +1, time, 2.0, 0.0)
	if err != nil {
		t.Fatal(err)
	}
	wantMoon := TimeFromCalendar(2019, 6, 25, 5, 34, 37.690)
	if diff := absDiff(moonrise.Ut, wantMoon.Ut); diff > 0.001 {
		t.Errorf("moonrise error = %g days", diff)
	}
}

func TestSearchHourAngle(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	hor, root, err := SearchHourAngle(Mars, observer, 4.5, time, +1)
	if err != nil {
		t.Fatal(err)
	}
	want := TimeFromCalendar(2019, 6, 24, 23, 33, 42.395)
	if diff := absDiff(root.Ut, want.Ut); diff > 0.001 {
		t.Errorf("search hour angle error = %g days", diff)
	}
	_ = hor
}

func TestSearchAltitude(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	root, err := SearchAltitude(Sun, observer, +1, time, 1.0, 10.0)
	if err != nil {
		t.Fatal(err)
	}
	want := TimeFromCalendar(2019, 6, 25, 11, 21, 0.0)
	if diff := absDiff(root.Ut, want.Ut); diff > 0.001 {
		t.Errorf("sun altitude search error = %g days", diff)
	}
}

func TestSunAltitude(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	alt, err := SunAltitude(time, observer)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(alt) > 90.0 {
		t.Errorf("implausible sun altitude = %g", alt)
	}
}
