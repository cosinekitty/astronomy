package astronomy

import (
	"testing"
)

func TestSearchMaxElongation(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	el, err := SearchMaxElongation(Venus, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(el.Time.Ut, 7388.425960236); diff > 1.0e-6 {
		t.Errorf("max elongation time error = %g days", diff)
	}
	if diff := absDiff(el.Elongation, 46.077543); diff > 1.0e-6 {
		t.Errorf("max elongation = %0.6f, expected 46.077543", el.Elongation)
	}
	if el.Visibility != VisibilityEvening {
		t.Errorf("visibility = %d, expected evening (1)", el.Visibility)
	}
	if _, err := SearchMaxElongation(Mars, time); err == nil {
		t.Error("expected error for non-Mercury/Venus")
	}
}

func TestSearchPeakMagnitude(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	illum, err := SearchPeakMagnitude(Venus, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(illum.Mag, -4.731874); diff > 1.0e-6 {
		t.Errorf("Venus peak magnitude = %0.6f, expected -4.731874", illum.Mag)
	}
	if diff := absDiff(illum.PhaseAngle, 118.995482); diff > 1.0e-6 {
		t.Errorf("Venus peak phase = %0.6f, expected 118.995482", illum.PhaseAngle)
	}
	if _, err := SearchPeakMagnitude(Mars, time); err == nil {
		t.Error("expected error for non-Venus")
	}
}

func TestSetDeltaTFunction(t *testing.T) {
	// Replace Delta T with a simple constant model and verify the effect on
	// the UT/TT conversion.
	SetDeltaTFunction(func(ut float64) float64 { return 60.0 })
	defer SetDeltaTFunction(nil) // restore the default
	time := TimeFromUniversalDays(0.0)
	if diff := absDiff(time.Tt, 60.0/86400.0); diff > 1.0e-15 {
		t.Errorf("TT with 60s Delta T = %0.15f, expected %0.15f", time.Tt, 60.0/86400.0)
	}
	SetDeltaTFunction(nil)
	// After restoring, the standard Espenak-Meeus value applies.
	time2 := TimeFromUniversalDays(0.0)
	if diff := absDiff(time2.Tt, 63.8285/86400.0); diff > 1.0e-6 {
		t.Errorf("TT with default Delta T = %0.9f", time2.Tt)
	}
}
