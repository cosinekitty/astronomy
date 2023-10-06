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
