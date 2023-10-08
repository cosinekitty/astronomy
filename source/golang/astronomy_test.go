package astronomy

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
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

	vec := GeoMoon(time)

	// GeoMoon should return:
	// Vector(0.002674037026701135, -0.0001531610316600666, -0.00031501599270694294, Time('2019-06-24T15:45:37.000Z'))
	dx := (vec.X - (0.002674037026701135))
	dy := (vec.Y - (-0.0001531610316600666))
	dz := (vec.Z - (-0.00031501599270694294))
	diff := math.Sqrt(dx*dx + dy*dy + dz*dz)
	if diff > 4.34e-19 {
		t.Errorf("Excessive position error for Moon: %g, x=%g, y=%g, z=%g", diff, vec.X, vec.Y, vec.Z)
	}
}

func float(s string) float64 {
	x, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic("Cannot parse float")
	}
	return x
}

func TestAtmosphere(t *testing.T) {
	const filename = "../../generate/riseset/atmosphere.csv"
	const tolerance = 8.8e-11
	file, err := os.Open(filename)
	if err != nil {
		t.Errorf("Cannot open file: %s", filename)
	}
	lnum := 0
	scanner := bufio.NewScanner(file)
	maxdiff := 0.0
	ncases := 0
	for scanner.Scan() {
		lnum += 1
		line := scanner.Text()
		if lnum > 1 {
			tokens := strings.Split(line, ",")
			if len(tokens) != 5 {
				t.Errorf("%s line %d: expected 5 numeric tokens, found %d", filename, lnum, len(tokens))
			}
			elevation := float(tokens[0])
			temperature := float(tokens[1])
			pressure := float(tokens[2])
			// ignore tokens[3] = absolute density
			relativeDensity := float(tokens[4])
			atmos, err := Atmosphere(elevation)
			if err != nil {
				t.Error(err)
			}

			var diff float64
			diff = math.Abs(atmos.Temperature - temperature)
			maxdiff = math.Max(maxdiff, diff)
			if diff > tolerance {
				t.Errorf("Excessive temperature diff = %g", diff)
			}

			diff = math.Abs(atmos.Pressure - pressure)
			maxdiff = math.Max(maxdiff, diff)
			if diff > tolerance {
				t.Errorf("Excessive pressure diff = %g", diff)
			}

			diff = math.Abs(atmos.Density - relativeDensity)
			maxdiff = math.Max(maxdiff, diff)
			if diff > tolerance {
				t.Errorf("Excessive density diff = %g", diff)
			}

			ncases += 1
		}
	}
	t.Logf("maxdiff = %g", maxdiff)
	if ncases != 34 {
		t.Errorf("Expected 34 test cases but processed %d", ncases)
	}
}
