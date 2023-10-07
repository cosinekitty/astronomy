package astronomy

import (
	"math"
	"testing"
)

func TestTerrestrialTime(t *testing.T) {
	ut1 := 8673.385964336088
	tt1 := 8673.386817342714
	time1 := TimeFromUniversalDays(ut1)
	diff1 := math.Abs(time1.Tt - tt1)
	if diff1 > 1.0e-16 {
		t.Errorf("Excessive Terrestrial Time Error = %g for time1", diff1)
	}

	ut2 := 8773.385964336088
	tt2 := 8773.386819208281
	time2 := time1.AddDays(100.0)
	diff2 := math.Abs(time2.Tt - tt2)
	if diff2 > 1.0e-16 {
		t.Errorf("Excessive Terrestrial Time Error = %g for time2", diff2)
	}

	utdiff2 := math.Abs((time2.Ut - ut2))
	if utdiff2 > 1.0e-16 {
		t.Errorf("Excessive Universal Time Error = %g for time2", utdiff2)
	}
}

func TestAngleBetween(t *testing.T) {
	time := TimeFromUniversalDays(0.0)
	a := AstroVector{1.0, 0.0, 0.0, time}
	b := AstroVector{1.0, 1.0, 0.0, time}
	angle := AngleBetween(a, b)
	diff := math.Abs(angle - 45.0)
	if diff > 1.0e-14 {
		t.Errorf("Excessive angle error = %g", diff)
	}
}

func TestCalendar(t *testing.T) {
	ut := 8679.201044872223
	cal, err := CalendarFromDays(ut)
	if err != nil {
		t.Error(err)
	}
	// (2023, 10, 6, 16, 49, 30.276975631713867)
	if cal.Year != 2023 {
		t.Errorf("Expected year = 2023 but found %d", cal.Year)
	}
	if cal.Month != 10 {
		t.Errorf("Expected month = 10 but found %d", cal.Month)
	}
	if cal.Day != 6 {
		t.Errorf("Expected day = 6 but found %d.", cal.Day)
	}
	if cal.Hour != 16 {
		t.Errorf("Expected hour = 16 but found %d.", cal.Hour)
	}
	if cal.Minute != 49 {
		t.Errorf("Expected minute = 49 but found %d.", cal.Minute)
	}
	diff := math.Abs(cal.Second - 30.276975631713867)
	if diff > 1.0e-16 {
		t.Errorf("Excessive error calculating calendar seconds: %0.16g", diff)
	}
}

func TestMoon(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	utdiff := math.Abs(time.Ut - 7114.156678240741)
	if utdiff > 1.0e-16 {
		t.Errorf("Excessive ut diff = %g", utdiff)
	}
	ttdiff := math.Abs(time.Tt - 7114.157503413434)
	if ttdiff > 1.0e-16 {
		t.Errorf("Excessive tt diff = %g", ttdiff)
	}

	// GeoMoon should return:
	// Vector(0.002674037026701135, -0.0001531610316600666, -0.00031501599270694294, Time('2019-06-24T15:45:37.000Z'))
}
