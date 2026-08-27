package main

import (
	"fmt"
	"math"
	"strconv"

	"github.com/cosinekitty/astronomy/source/golang"
)

// printTime formats an AstroTime as a UTC string with whole seconds and a
// trailing 'Z', rounding to the nearest second, to match the C demo PrintTime.
func printTime(t astronomy.AstroTime) string {
	rounded := t.AddDays(0.5 / 86400.0)
	cal, err := astronomy.CalendarFromDays(rounded.Ut)
	if err != nil {
		return "<error>"
	}
	return fmt.Sprintf("%04d-%02d-%02dT%02d:%02d:%02dZ",
		cal.Year, cal.Month, cal.Day, cal.Hour, cal.Minute, int(cal.Second))
}

// quarterName maps a Moon quarter index to its name, matching the C demo.
func quarterName(quarter int) string {
	switch quarter {
	case 0:
		return "New Moon"
	case 1:
		return "First Quarter"
	case 2:
		return "Full Moon"
	case 3:
		return "Third Quarter"
	default:
		return "INVALID QUARTER"
	}
}

// demoPositions displays equatorial and horizontal coordinates of the Sun,
// Moon, and planets.
func demoPositions(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	fmt.Printf("UTC date = ")
	fmt.Print(printTime(time))
	fmt.Printf("\n")
	fmt.Printf("BODY           RA      DEC       AZ      ALT\n")
	bodies := []astronomy.Body{
		astronomy.Sun, astronomy.Moon, astronomy.Mercury, astronomy.Venus, astronomy.Mars,
		astronomy.Jupiter, astronomy.Saturn, astronomy.Uranus, astronomy.Neptune, astronomy.Pluto,
	}
	for _, body := range bodies {
		equ2000, err := astronomy.Equator(body, &time, observer, astronomy.EquatorJ2000, astronomy.WithAberration)
		if err != nil {
			return err
		}
		equOfDate, err := astronomy.Equator(body, &time, observer, astronomy.EquatorOfDate, astronomy.WithAberration)
		if err != nil {
			return err
		}
		hor := astronomy.Horizon(time, observer, equOfDate.Ra, equOfDate.Dec, astronomy.NormalRefraction)
		fmt.Printf("%-8s %8.2f %8.2f %8.2f %8.2f\n", astronomy.BodyName(body), equ2000.Ra, equ2000.Dec, hor.Azimuth, hor.Altitude)
	}
	return nil
}

// demoRiseSet displays sunrise, sunset, moonrise, and moonset times.
func demoRiseSet(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	fmt.Printf("search   : ")
	fmt.Print(printTime(time))
	fmt.Printf("\n")
	events := []struct {
		name      string
		body      astronomy.Body
		direction int
	}{
		{"sunrise", astronomy.Sun, +1},
		{"sunset", astronomy.Sun, -1},
		{"moonrise", astronomy.Moon, +1},
		{"moonset", astronomy.Moon, -1},
	}
	for _, ev := range events {
		fmt.Printf("%-8s : ", ev.name)
		root, err := astronomy.SearchRiseSet(ev.body, observer, ev.direction, time, 300.0, 0.0)
		if err != nil {
			fmt.Printf("ERROR %v\n", err)
			continue
		}
		fmt.Print(printTime(root))
		fmt.Printf("\n")
	}
	return nil
}

// demoMoonPhase shows the Moon's phase and the next 10 quarter phases.
func demoMoonPhase(args []string) error {
	var time astronomy.AstroTime
	var err error
	if len(args) >= 1 {
		time, err = parseTime(args[0])
	} else {
		time = astronomy.CurrentTime()
	}
	if err != nil {
		return err
	}

	phase, err := astronomy.MoonPhase(time)
	if err != nil {
		return err
	}
	fmt.Print(printTime(time))
	fmt.Printf(" : Moon's ecliptic phase angle = %.3f degrees.\n", phase)

	illum, err := astronomy.Illumination(astronomy.Moon, time)
	if err != nil {
		return err
	}
	fmt.Print(printTime(time))
	fmt.Printf(" : Moon's illuminated fraction = %.2f%%.", 100.0*illum.PhaseFraction)

	fmt.Printf("\nThe next 10 lunar quarters are:\n")
	mq, err := astronomy.SearchMoonQuarter(time)
	if err != nil {
		return err
	}
	for i := 0; i < 10; i++ {
		fmt.Print(printTime(mq.Time))
		fmt.Printf(" : %s\n", quarterName(mq.Quarter))
		mq, err = astronomy.NextMoonQuarter(mq)
		if err != nil {
			return err
		}
	}
	return nil
}

// demoSeasons displays the equinoxes and solstices of a calendar year.
func demoSeasons(args []string) error {
	if len(args) != 1 {
		return fmt.Errorf("USAGE: seasons year")
	}
	year, err := strconv.Atoi(args[0])
	if err != nil {
		return err
	}
	sea, err := astronomy.Seasons(year)
	if err != nil {
		return err
	}
	rows := []struct {
		name string
		time astronomy.AstroTime
	}{
		{"March equinox", sea.MarEquinox},
		{"June solstice", sea.JunSolstice},
		{"September equinox", sea.SepEquinox},
		{"December solstice", sea.DecSolstice},
	}
	for _, r := range rows {
		fmt.Printf("%-17s : ", r.name)
		fmt.Print(printTime(r.time))
		fmt.Printf("\n")
	}
	return nil
}

// demoCulminate finds the next culmination (meridian crossing) of the Sun,
// Moon, and planets as seen by an observer.
func demoCulminate(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	fmt.Printf("search   : ")
	fmt.Print(printTime(time))
	fmt.Printf("\n")
	bodies := []astronomy.Body{
		astronomy.Sun, astronomy.Moon, astronomy.Mercury, astronomy.Venus, astronomy.Mars,
		astronomy.Jupiter, astronomy.Saturn, astronomy.Uranus, astronomy.Neptune, astronomy.Pluto,
	}
	for _, body := range bodies {
		fmt.Printf("%-8s : ", astronomy.BodyName(body))
		hor, root, err := astronomy.SearchHourAngle(body, observer, 0.0, time, +1)
		if err != nil {
			fmt.Printf("ERROR %v\n", err)
			continue
		}
		fmt.Print(printTime(root))
		fmt.Printf("  altitude=%6.2f  azimuth=%6.2f\n", hor.Altitude, hor.Azimuth)
	}
	return nil
}

// horizontalCoords converts an ecliptic longitude to horizontal coordinates
// using a precomputed ecliptic->horizontal rotation matrix.
func horizontalCoords(eclipticLongitude float64, time astronomy.AstroTime, rot astronomy.RotationMatrix) (astronomy.Spherical, error) {
	eclip := astronomy.Spherical{Lat: 0.0, Lon: eclipticLongitude, Dist: 1.0}
	eclVec := astronomy.VectorFromSphere(eclip, time)
	horVec := astronomy.RotateVector(rot, eclVec)
	return astronomy.HorizonFromVector(horVec, astronomy.NormalRefraction)
}

// searchForCrossing uses a binary search on ecliptic longitude to find where
// the horizontal altitude crosses zero between e1 and e2.
func searchForCrossing(e1, e2 float64, time astronomy.AstroTime, rot astronomy.RotationMatrix) (float64, astronomy.Spherical, error) {
	const tolerance = 1.0e-6 // one-millionth of a degree is close enough
	for {
		e3 := (e1 + e2) / 2.0
		h3, err := horizontalCoords(e3, time, rot)
		if err != nil {
			return 0.0, astronomy.Spherical{}, err
		}
		if math.Abs(e2-e1) < tolerance {
			return e3, h3, nil
		}
		if h3.Lat < 0.0 {
			e1 = e3
		} else {
			e2 = e3
		}
	}
}

// demoHorizon finds the two azimuths where the ecliptic intersects the
// observer's horizon, matching the C horizon demo.
func demoHorizon(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	const numSamples = 4
	ecliplon := func(i int) float64 { return (360.0 * float64(i)) / numSamples }

	// Rotation matrix that converts J2000 ecliptic vectors to horizontal vectors.
	rot := astronomy.Rotation_ECL_HOR(time, observer)

	// Sample several points around the ecliptic.
	hor := make([]astronomy.Spherical, numSamples)
	for i := 0; i < numSamples; i++ {
		hor[i], err = horizontalCoords(ecliplon(i), time, rot)
		if err != nil {
			return err
		}
	}

	for i := 0; i < numSamples; i++ {
		a1 := hor[i].Lat
		a2 := hor[(i+1)%numSamples].Lat
		e1 := ecliplon(i)
		e2 := ecliplon(i + 1)
		if a1*a2 <= 0.0 {
			var ex float64
			var hx astronomy.Spherical
			if a2 > a1 {
				ex, hx, err = searchForCrossing(e1, e2, time, rot)
			} else {
				ex, hx, err = searchForCrossing(e2, e1, time, rot)
			}
			if err != nil {
				return err
			}
			direction := "descends"
			if hx.Lon > 0.0 && hx.Lon < 180.0 {
				direction = "ascends "
			}
			fmt.Printf("Ecliptic longitude %9.4f %s through horizon az %9.4f, alt %12.6g\n",
				ex, direction, hx.Lon, hx.Lat)
		}
	}
	return nil
}

// demoLunarEclipse displays the next 10 partial/total lunar eclipses.
func demoLunarEclipse(args []string) error {
	var time astronomy.AstroTime
	var err error
	if len(args) >= 1 {
		time, err = parseTime(args[0])
	} else {
		time = astronomy.CurrentTime()
	}
	if err != nil {
		return err
	}

	const minutesPerDay = 24.0 * 60.0
	eclipse, err := astronomy.SearchLunarEclipse(time)
	if err != nil {
		return err
	}
	count := 0
	for {
		if eclipse.Kind != astronomy.EclipseKindPenumbral {
			fmt.Print(printTime(eclipse.Peak.AddDays(-eclipse.SdPartial / minutesPerDay)))
			fmt.Printf(" - Partial eclipse begins.\n")
			if eclipse.SdTotal > 0.0 {
				fmt.Print(printTime(eclipse.Peak.AddDays(-eclipse.SdTotal / minutesPerDay)))
				fmt.Printf(" - Total eclipse begins.\n")
			}
			kindName := "partial"
			if eclipse.Kind == astronomy.EclipseKindTotal {
				kindName = "total"
			}
			fmt.Print(printTime(eclipse.Peak))
			fmt.Printf(" - Peak of %s eclipse.\n", kindName)
			if eclipse.SdTotal > 0.0 {
				fmt.Print(printTime(eclipse.Peak.AddDays(+eclipse.SdTotal / minutesPerDay)))
				fmt.Printf(" - Total eclipse ends.\n")
			}
			fmt.Print(printTime(eclipse.Peak.AddDays(+eclipse.SdPartial / minutesPerDay)))
			fmt.Printf(" - Partial eclipse ends.\n\n")
			count++
			if count == 10 {
				return nil
			}
		}
		eclipse, err = astronomy.NextLunarEclipse(eclipse.Peak)
		if err != nil {
			return err
		}
	}
}

// demoSolarTime displays the true solar time for the observer.
func demoSolarTime(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	ha, err := astronomy.HourAngle(astronomy.Sun, &time, observer)
	if err != nil {
		return err
	}
	solarTimeHours := math.Mod(ha+12.0, 24.0)

	milli := int(math.Round(solarTimeHours * 3600000.0))
	second := milli / 1000
	milli %= 1000
	minute := second / 60
	second %= 60
	hour := minute / 60
	minute %= 60
	hour %= 24

	fmt.Printf("True solar time = %.4f hours (%02d:%02d:%02d.%03d)\n",
		solarTimeHours, hour, minute, second, milli)
	return nil
}

// demoGravity displays the acceleration of gravity at a latitude and height.
func demoGravity(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("USAGE: gravity latitude height")
	}
	lat, err := strconv.ParseFloat(args[0], 64)
	if err != nil {
		return err
	}
	height := 0.0
	if len(args) >= 2 {
		height, err = strconv.ParseFloat(args[1], 64)
		if err != nil {
			return err
		}
	}
	fmt.Printf("latitude = %8.4f,  height = %6.0f,  gravity = %8.6f\n",
		lat, height, astronomy.ObserverGravity(lat, height))
	return nil
}

// demoGalactic converts galactic coordinates to horizontal coordinates.
func demoGalactic(args []string) error {
	if len(args) < 4 {
		return fmt.Errorf("USAGE: galactic olat olon glat glon [yyyy-mm-ddThh:mm:ssZ]")
	}
	lat, err := strconv.ParseFloat(args[0], 64)
	if err != nil {
		return err
	}
	lon, err := strconv.ParseFloat(args[1], 64)
	if err != nil {
		return err
	}
	glat, err := strconv.ParseFloat(args[2], 64)
	if err != nil {
		return err
	}
	glon, err := strconv.ParseFloat(args[3], 64)
	if err != nil {
		return err
	}
	observer := astronomy.Observer{Latitude: lat, Longitude: lon}
	var time astronomy.AstroTime
	if len(args) >= 5 {
		time, err = parseTime(args[4])
		if err != nil {
			return err
		}
	} else {
		time = astronomy.CurrentTime()
	}

	// Rotation matrix that converts galactic coordinates to horizontal coordinates.
	rot := astronomy.CombineRotation(astronomy.Rotation_GAL_EQJ(), astronomy.Rotation_EQJ_HOR(time, observer))

	gsphere := astronomy.Spherical{Lat: glat, Lon: glon, Dist: 1.0}
	gvec := astronomy.VectorFromSphere(gsphere, time)
	hvec := astronomy.RotateVector(rot, gvec)
	hsphere, err := astronomy.HorizonFromVector(hvec, astronomy.NoRefraction)
	if err != nil {
		return err
	}
	fmt.Printf("altitude = %10.3f, azimuth = %10.3f\n", hsphere.Lat, hsphere.Lon)
	return nil
}
