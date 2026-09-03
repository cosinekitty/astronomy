package astronomy

import (
	"testing"
)

func TestSearchLunarEclipse(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	eclipse, err := SearchLunarEclipse(time)
	if err != nil {
		t.Fatal(err)
	}
	if eclipse.Kind != EclipseKindPartial {
		t.Errorf("lunar eclipse kind = %d, expected partial (2)", eclipse.Kind)
	}
	if diff := absDiff(eclipse.Obscuration, 0.661139131); diff > 1.0e-6 {
		t.Errorf("lunar obscuration = %0.9f, expected 0.661139131", eclipse.Obscuration)
	}
	if diff := absDiff(eclipse.SdPenum, 167.235248); diff > 1.0e-4 {
		t.Errorf("penumbral semi-duration = %0.6f, expected 167.235248", eclipse.SdPenum)
	}
	if diff := absDiff(eclipse.SdPartial, 89.409181); diff > 1.0e-4 {
		t.Errorf("partial semi-duration = %0.6f, expected 89.409181", eclipse.SdPartial)
	}
	want := TimeFromCalendar(2019, 7, 16, 21, 30, 41.806)
	if diff := absDiff(eclipse.Peak.Ut, want.Ut); diff > 1.0e-4 {
		t.Errorf("lunar eclipse peak error = %g days", diff)
	}
}

func TestNextLunarEclipse(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	first, err := SearchLunarEclipse(time)
	if err != nil {
		t.Fatal(err)
	}
	second, err := NextLunarEclipse(first.Peak)
	if err != nil {
		t.Fatal(err)
	}
	if second.Peak.Ut <= first.Peak.Ut {
		t.Error("next lunar eclipse should be after the first")
	}
}

func TestSearchGlobalSolarEclipse(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	eclipse, err := SearchGlobalSolarEclipse(time)
	if err != nil {
		t.Fatal(err)
	}
	if eclipse.Kind != EclipseKindTotal {
		t.Errorf("global solar eclipse kind = %d, expected total (4)", eclipse.Kind)
	}
	if eclipse.Obscuration != 1.0 {
		t.Errorf("total eclipse obscuration should be 1.0, got %g", eclipse.Obscuration)
	}
	if diff := absDiff(eclipse.Distance, 4123.831689); diff > 1.0e-3 {
		t.Errorf("eclipse distance = %0.6f, expected 4123.831689", eclipse.Distance)
	}
	if diff := absDiff(eclipse.Latitude, -17.388447); diff > 1.0e-4 {
		t.Errorf("eclipse latitude = %0.6f, expected -17.388447", eclipse.Latitude)
	}
	if diff := absDiff(eclipse.Longitude, -108.990531); diff > 1.0e-4 {
		t.Errorf("eclipse longitude = %0.6f, expected -108.990531", eclipse.Longitude)
	}
	want := TimeFromCalendar(2019, 7, 2, 19, 22, 56.168)
	if diff := absDiff(eclipse.Peak.Ut, want.Ut); diff > 1.0e-4 {
		t.Errorf("global solar eclipse peak error = %g days", diff)
	}
}

func TestSearchLocalSolarEclipse(t *testing.T) {
	// Reference from the C implementation (the next visible local eclipse at
	// the given observer is the 2023 eclipse, seen as partial from Florida).
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	eclipse, err := SearchLocalSolarEclipse(time, observer)
	if err != nil {
		t.Fatal(err)
	}
	if eclipse.Kind != EclipseKindPartial {
		t.Errorf("local solar eclipse kind = %d, expected partial (2)", eclipse.Kind)
	}
	if diff := absDiff(eclipse.Obscuration, 0.535128698); diff > 1.0e-6 {
		t.Errorf("local obscuration = %0.9f, expected 0.535128698", eclipse.Obscuration)
	}
	check := func(label string, got AstroTime, want float64) {
		if diff := absDiff(got.Ut, want); diff > 1.0e-4 {
			t.Errorf("%s error = %g days (got %0.9f, want %0.9f)", label, diff, got.Ut, want)
		}
	}
	check("partial_begin", eclipse.PartialBegin.Time, 8687.161661996)
	check("peak", eclipse.Peak.Time, 8687.226469609)
	check("partial_end", eclipse.PartialEnd.Time, 8687.292867498)
}

func TestSearchTransit(t *testing.T) {
	// Reference from the C implementation: Mercury transit of 2019-11-11.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	transit, err := SearchTransit(Mercury, time)
	if err != nil {
		t.Fatal(err)
	}
	check := func(label string, got AstroTime, want float64) {
		if diff := absDiff(got.Ut, want); diff > 1.0e-4 {
			t.Errorf("%s error = %g days", label, diff)
		}
	}
	check("transit start", transit.Start, 7254.024088962)
	check("transit peak", transit.Peak, 7254.138119051)
	check("transit finish", transit.Finish, 7254.252211598)
	if diff := absDiff(transit.Separation, 1.319527); diff > 1.0e-4 {
		t.Errorf("transit separation = %0.6f, expected 1.319527", transit.Separation)
	}

	next, err := NextTransit(Mercury, transit.Peak)
	if err != nil {
		t.Fatal(err)
	}
	if next.Peak.Ut <= transit.Peak.Ut {
		t.Error("next transit should be after the first")
	}

	if _, err := SearchTransit(Mars, time); err == nil {
		t.Error("expected error for non-Mercury/Venus transit")
	}
}
