package astronomy

import (
	"math"
	"testing"
)

func TestTerrestrialTime(t *testing.T) {
	ut := 8673.385964336088
	tt := TerrestrialTime(ut)
	diff := math.Abs(tt - 8673.386817342714)
	if diff > 1.0e-16 {
		t.Errorf("Excessive Terrestrial Time Error = %g", diff)
	}
}
