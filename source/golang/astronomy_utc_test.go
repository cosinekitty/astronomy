package astronomy

import (
	"testing"
)

func TestUtcRoundTrip(t *testing.T) {
	// A UTC time should round-trip through AstroTime.
	utc := UtcTime{Year: 2023, Month: 10, Day: 6, Hour: 16, Minute: 49, Second: 30.5}
	time := TimeFromUtc(utc)
	back := UtcFromTime(time)
	if back.Year != utc.Year || back.Month != utc.Month || back.Day != utc.Day {
		t.Errorf("date mismatch: %+v vs %+v", back, utc)
	}
	if back.Hour != utc.Hour || back.Minute != utc.Minute {
		t.Errorf("time mismatch: %+v vs %+v", back, utc)
	}
	if diff := absDiff(back.Second, utc.Second); diff > 1e-4 {
		t.Errorf("second mismatch: %g", diff)
	}
}

func TestFormatTime(t *testing.T) {
	time := TimeFromCalendar(2023, 10, 6, 16, 49, 30.5)
	if got := FormatTime(time, TimeFormatDay); got != "2023-10-06" {
		t.Errorf("day format = %q", got)
	}
	if got := FormatTime(time, TimeFormatMinute); got != "2023-10-06T16:49" {
		t.Errorf("minute format = %q", got)
	}
	if got := FormatTime(time, TimeFormatSecond); got != "2023-10-06T16:49:30" {
		t.Errorf("second format = %q", got)
	}
	if got := FormatTime(time, TimeFormatMilli); got != "2023-10-06T16:49:30.500" {
		t.Errorf("milli format = %q", got)
	}
}

func TestReset(t *testing.T) {
	// Reset must clear the star definitions and allow them to be redefined.
	if err := DefineStar(Star1, 3.0, 40.0, 50.0); err != nil {
		t.Fatal(err)
	}
	Reset()
	if userDefinedStar(Star1) != nil {
		t.Error("Reset should clear user-defined stars")
	}
	// Reset must not break subsequent calculations.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	if _, err := HelioVector(Mars, time); err != nil {
		t.Fatal(err)
	}
}
