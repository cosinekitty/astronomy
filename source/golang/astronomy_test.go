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

func TestSiderealTime(t *testing.T) {
	const correct = 9.3983699280076483
	var diff float64
	time := TimeFromCalendar(2022, 3, 15, 21, 50, 0.0)
	t.Logf("time.tt = %0.16f, time.ut = %0.16f", time.Tt, time.Ut)
	et := etilt(&time)
	t.Logf("etilt tt=%0.16f, dpsi=%0.16f, deps=%0.16f, ee=%0.16f, mobl=%0.16f, tobl=%0.16f", et.Tt, et.Dpsi, et.Deps, et.Ee, et.Mobl, et.Tobl)
	diff = math.Abs(et.Tt - 8109.4105648890845)
	if diff > 1.0e-16 {
		t.Errorf("et.Tt error = %g", diff)
	}
	diff = math.Abs(et.Dpsi - (-13.38138780522237))
	if diff > 1.0e-16 {
		t.Errorf("et.Dpsi error = %g", diff)
	}
	diff = math.Abs(et.Deps - 5.8448952225938449)
	if diff > 1.0e-16 {
		t.Errorf("et.Deps error = %g", diff)
	}
	diff = math.Abs(et.Ee - (-0.8184968463212278))
	if diff > 1.0e-16 {
		t.Errorf("et.Ee error = %g", diff)
	}
	diff = math.Abs(et.Mobl - 23.436390874072426)
	if diff > 1.0e-16 {
		t.Errorf("et.Mobl error = %g", diff)
	}
	diff = math.Abs(et.Tobl - 23.4380144560787)
	if diff > 1.0e-16 {
		t.Errorf("et.Tobl error = %g", diff)
	}
	gast := SiderealTime(&time)
	diff = math.Abs(gast - correct)
	t.Logf("gast=%0.16f, correct=%0.16f, diff=%g", gast, correct, diff)
	if diff > 1.0e-15 {
		t.Error("Excessive sidereal time error")
	}
}

func jmverify(t *testing.T, m StateVector, x, y, z, vx, vy, vz float64) {
	const posTolerance = 9.0e-4
	const velTolerance = 9.0e-4
	var dx, dy, dz, mag float64
	dx = x - m.X
	dy = y - m.Y
	dz = z - m.Z
	mag = math.Sqrt(x*x + y*y + z*z)
	posDiff := math.Sqrt(dx*dx+dy*dy+dz*dz) / mag
	if posDiff > posTolerance {
		t.Errorf("jmverify: excessive position error %g", posDiff)
	}
	dx = vx - m.Vx
	dy = vy - m.Vy
	dz = vz - m.Vz
	mag = math.Sqrt(vx*vx + vy*vy + vz*vz)
	velDiff := math.Sqrt(dx*dx+dy*dy+dz*dz) / mag
	if velDiff > velTolerance {
		t.Errorf("jmverify: excessive velocity error %g", velDiff)
	}
	t.Logf("jmverify: posdiff=%g, veldiff=%g", posDiff, velDiff)
}

func jmtest(t *testing.T, jd float64,
	x0, y0, z0 float64,
	vx0, vy0, vz0 float64,
	x1, y1, z1 float64,
	vx1, vy1, vz1 float64,
	x2, y2, z2 float64,
	vx2, vy2, vz2 float64,
	x3, y3, z3 float64,
	vx3, vy3, vz3 float64) {
	tt := jd - 2451545.0 // convert JD to J2000 TT

	time := TimeFromTerrestrialDays(tt)
	jm := JupiterMoons(time)
	jmverify(t, jm.Io, x0, y0, z0, vx0, vy0, vz0)
	jmverify(t, jm.Europa, x1, y1, z1, vx1, vy1, vz1)
	jmverify(t, jm.Ganymede, x2, y2, z2, vx2, vy2, vz2)
	jmverify(t, jm.Callisto, x3, y3, z3, vx3, vy3, vz3)
}

func TestJupiterMoons(t *testing.T) {
	// 2426545.000000000 = A.D. 1931-Jul-22 12:00:00.0000 TDB
	jmtest(t, 2426545.0,
		1.982926146239590e-03, -1.832295139436290e-03, -8.421212456711890e-04,
		7.134149248554827e-03, 6.252992190334442e-03, 3.104307708231514e-03,
		-2.001703395523389e-03, -3.616838721655651e-03, -1.778587701398564e-03,
		7.114838436292209e-03, -3.188565831252877e-03, -1.343652914203527e-03,
		-3.569347026083638e-03, 5.624051966893631e-03, 2.602293356154182e-03,
		-5.449491834975323e-03, -2.794160033857133e-03, -1.420730759322332e-03,
		1.557331942227016e-03, 1.133602597831958e-02, 5.408245936181074e-03,
		-4.672704130230531e-03, 5.635487435093144e-04, 2.184495893761995e-04)

	//	2441335.000000000 = A.D. 1972-Jan-18 12:00:00.0000 TDB
	jmtest(t, 2441335.0,
		2.007699805595389e-03, -1.797454980051448e-03, -8.231484302662159e-04,
		6.998721305470736e-03, 6.429630159128681e-03, 3.171049780949413e-03,
		-4.443100267506271e-03, 2.925401293242437e-04, 8.409570852858875e-05,
		-4.902927233366754e-04, -7.226402913797873e-03, -3.382285097453537e-03,
		-6.829630727289741e-03, -1.889906681645492e-03, -9.820646931516248e-04,
		1.867675646626441e-03, -5.432667808431406e-03, -2.549409251738509e-03,
		-1.179742330021488e-02, 4.258274263479233e-03, 1.861994234230584e-03,
		-1.718979470855996e-03, -3.945479662311145e-03, -1.896520523122100e-03)

	// 2460555.000000000 = A.D. 2024-Sep-01 12:00:00.0000 TDB
	jmtest(t, 2460555.0,
		2.540448557411970e-03, -1.127014656496033e-03, -4.961264601034464e-04,
		4.394105339296697e-03, 8.082707005382311e-03, 3.917932947639376e-03,
		6.442853937260593e-04, -4.041048874596572e-03, -1.897349473425028e-03,
		7.801383991098988e-03, 1.000491184714297e-03, 6.643020423094265e-04,
		6.782582099995973e-03, -2.062347066586191e-03, -8.823648306081377e-04,
		1.978362595879211e-03, 5.375818246356209e-03, 2.609500123996761e-03,
		1.025314438251558e-02, 6.413374763992526e-03, 3.174110303332572e-03,
		-2.716265475976206e-03, 3.562173193821677e-03, 1.638391747868563e-03)
}

func gravityCheck(t *testing.T, latitude, height, expected float64) {
	actual := ObserverGravity(latitude, height)
	diff := math.Abs(actual - expected)
	if diff > 5.0e-7 {
		t.Errorf("Excessive gravity error = %g for latitude = %f, height = %f", diff, latitude, height)
	}
}

func TestObserverGravity(t *testing.T) {
	gravityCheck(t, 0.0000, 0.0, 9.780325)
	gravityCheck(t, 1.0000, 0.0, 9.780341)
	gravityCheck(t, 2.0000, 0.0, 9.780388)
	gravityCheck(t, 3.0000, 0.0, 9.780467)
	gravityCheck(t, 4.0000, 0.0, 9.780577)
	gravityCheck(t, 5.0000, 0.0, 9.780718)
	gravityCheck(t, 6.0000, 0.0, 9.780889)
	gravityCheck(t, 7.0000, 0.0, 9.781092)
	gravityCheck(t, 8.0000, 0.0, 9.781325)
	gravityCheck(t, 9.0000, 0.0, 9.781589)
	gravityCheck(t, 10.0000, 0.0, 9.781882)
	gravityCheck(t, 11.0000, 0.0, 9.782205)
	gravityCheck(t, 12.0000, 0.0, 9.782558)
	gravityCheck(t, 13.0000, 0.0, 9.782939)
	gravityCheck(t, 14.0000, 0.0, 9.783348)
	gravityCheck(t, 15.0000, 0.0, 9.783785)
	gravityCheck(t, 16.0000, 0.0, 9.784249)
	gravityCheck(t, 17.0000, 0.0, 9.784740)
	gravityCheck(t, 18.0000, 0.0, 9.785258)
	gravityCheck(t, 19.0000, 0.0, 9.785800)
	gravityCheck(t, 20.0000, 0.0, 9.786368)
	gravityCheck(t, 21.0000, 0.0, 9.786960)
	gravityCheck(t, 22.0000, 0.0, 9.787575)
	gravityCheck(t, 23.0000, 0.0, 9.788213)
	gravityCheck(t, 24.0000, 0.0, 9.788873)
	gravityCheck(t, 25.0000, 0.0, 9.789554)
	gravityCheck(t, 26.0000, 0.0, 9.790256)
	gravityCheck(t, 27.0000, 0.0, 9.790976)
	gravityCheck(t, 28.0000, 0.0, 9.791716)
	gravityCheck(t, 29.0000, 0.0, 9.792473)
	gravityCheck(t, 30.0000, 0.0, 9.793247)
	gravityCheck(t, 31.0000, 0.0, 9.794037)
	gravityCheck(t, 32.0000, 0.0, 9.794842)
	gravityCheck(t, 33.0000, 0.0, 9.795661)
	gravityCheck(t, 34.0000, 0.0, 9.796492)
	gravityCheck(t, 35.0000, 0.0, 9.797336)
	gravityCheck(t, 36.0000, 0.0, 9.798191)
	gravityCheck(t, 37.0000, 0.0, 9.799055)
	gravityCheck(t, 38.0000, 0.0, 9.799928)
	gravityCheck(t, 39.0000, 0.0, 9.800809)
	gravityCheck(t, 40.0000, 0.0, 9.801697)
	gravityCheck(t, 41.0000, 0.0, 9.802590)
	gravityCheck(t, 42.0000, 0.0, 9.803488)
	gravityCheck(t, 43.0000, 0.0, 9.804389)
	gravityCheck(t, 44.0000, 0.0, 9.805293)
	gravityCheck(t, 45.0000, 0.0, 9.806198)
	gravityCheck(t, 46.0000, 0.0, 9.807103)
	gravityCheck(t, 47.0000, 0.0, 9.808007)
	gravityCheck(t, 48.0000, 0.0, 9.808909)
	gravityCheck(t, 49.0000, 0.0, 9.809808)
	gravityCheck(t, 50.0000, 0.0, 9.810702)
	gravityCheck(t, 51.0000, 0.0, 9.811591)
	gravityCheck(t, 52.0000, 0.0, 9.812474)
	gravityCheck(t, 53.0000, 0.0, 9.813349)
	gravityCheck(t, 54.0000, 0.0, 9.814216)
	gravityCheck(t, 55.0000, 0.0, 9.815073)
	gravityCheck(t, 56.0000, 0.0, 9.815919)
	gravityCheck(t, 57.0000, 0.0, 9.816754)
	gravityCheck(t, 58.0000, 0.0, 9.817576)
	gravityCheck(t, 59.0000, 0.0, 9.818384)
	gravityCheck(t, 60.0000, 0.0, 9.819177)
	gravityCheck(t, 61.0000, 0.0, 9.819955)
	gravityCheck(t, 62.0000, 0.0, 9.820715)
	gravityCheck(t, 63.0000, 0.0, 9.821459)
	gravityCheck(t, 64.0000, 0.0, 9.822183)
	gravityCheck(t, 65.0000, 0.0, 9.822889)
	gravityCheck(t, 66.0000, 0.0, 9.823574)
	gravityCheck(t, 67.0000, 0.0, 9.824238)
	gravityCheck(t, 68.0000, 0.0, 9.824880)
	gravityCheck(t, 69.0000, 0.0, 9.825499)
	gravityCheck(t, 70.0000, 0.0, 9.826095)
	gravityCheck(t, 71.0000, 0.0, 9.826666)
	gravityCheck(t, 72.0000, 0.0, 9.827213)
	gravityCheck(t, 73.0000, 0.0, 9.827734)
	gravityCheck(t, 74.0000, 0.0, 9.828229)
	gravityCheck(t, 75.0000, 0.0, 9.828697)
	gravityCheck(t, 76.0000, 0.0, 9.829137)
	gravityCheck(t, 77.0000, 0.0, 9.829550)
	gravityCheck(t, 78.0000, 0.0, 9.829934)
	gravityCheck(t, 79.0000, 0.0, 9.830289)
	gravityCheck(t, 80.0000, 0.0, 9.830614)
	gravityCheck(t, 81.0000, 0.0, 9.830910)
	gravityCheck(t, 82.0000, 0.0, 9.831176)
	gravityCheck(t, 83.0000, 0.0, 9.831411)
	gravityCheck(t, 84.0000, 0.0, 9.831616)
	gravityCheck(t, 85.0000, 0.0, 9.831789)
	gravityCheck(t, 86.0000, 0.0, 9.831931)
	gravityCheck(t, 87.0000, 0.0, 9.832042)
	gravityCheck(t, 88.0000, 0.0, 9.832121)
	gravityCheck(t, 89.0000, 0.0, 9.832169)
	gravityCheck(t, 90.0000, 0.0, 9.832185)
}

func refractionCheck(t *testing.T, refraction Refraction, altitude, expected float64) {
	actual := RefractionAngle(refraction, altitude)
	diff := math.Abs(expected - actual)
	if diff > 1.0e-16 {
		t.Errorf("Refraction discrepancy = %g", diff)
	}

	check := InverseRefractionAngle(refraction, altitude+actual)
	diff = math.Abs(expected + check)
	if diff > 2.6e-15 {
		t.Errorf("Inverse refraction discrepancy = %g, expected = %g, check = %g", diff, expected, check)
	}
}

func TestRefraction(t *testing.T) {
	refractionCheck(t, NormalRefraction, 0.0, 0.4830321230741662)
	refractionCheck(t, NoRefraction, 0.0, 0.0)
	refractionCheck(t, NormalRefraction, 10.0, 0.09012801338558875)
	refractionCheck(t, NormalRefraction, 20.0, 0.04568673807086863)
	refractionCheck(t, NormalRefraction, -80.0, 0.0726495079641601)
	refractionCheck(t, NormalRefraction, -90.0, 0.0)
}

func TestHorizon(t *testing.T) {
	// Test case generated using the Python version of Astronomy Engine.
	// 2023-11-01T18:07:13.348Z
	time := TimeFromCalendar(2023, 11, 1, 18, 7, 13.348)
	ttExpected := 8705.255869201963
	if diff := math.Abs(ttExpected - time.Tt); diff > 1.0e-16 {
		t.Errorf("excessive tt diff = %g", diff)
	}
	utExpected := 8705.255015601851
	if diff := math.Abs(utExpected - time.Ut); diff > 1.0e-16 {
		t.Errorf("excessive ut diff = %g", diff)
	}

	observer := Observer{Latitude: 35.0, Longitude: -87.0, Height: 0.0}
	ra := 12.0
	dec := 5.0
	topo, err := Horizon(time, observer, ra, dec, NormalRefraction)
	if err != nil {
		t.Fatal(err)
	} else {
		azimuthExpected := 245.2118420327564
		altitudeExpected := 38.42409530728063
		raExpected := 12.001062402159066
		decExpected := 5.014149039382782

		if diff := math.Abs(topo.Azimuth - azimuthExpected); diff > 1.0e-13 {
			t.Errorf("excessive azimuth diff = %g, calculated = %g", diff, topo.Azimuth)
		}

		if diff := math.Abs(topo.Altitude - altitudeExpected); diff > 3.0e-14 {
			t.Errorf("excessive altitude diff = %g, calculated = %g", diff, topo.Altitude)
		}

		if diff := math.Abs(topo.Ra - raExpected); diff > 1.0e-16 {
			t.Errorf("excessive RA diff = %g, calculated = %g", diff, topo.Ra)
		}

		if diff := math.Abs(topo.Dec - decExpected); diff > 1.0e-16 {
			t.Errorf("excessive DEC diff = %g, calculated = %g", diff, topo.Dec)
		}
	}
}
