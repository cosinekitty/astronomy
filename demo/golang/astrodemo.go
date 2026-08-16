package main

import (
	"fmt"
	"strconv"

	"github.com/cosinekitty/astronomy/source/golang"
)

// demoPositions displays equatorial and horizontal coordinates of the Sun,
// Moon, and planets.
func demoPositions(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	fmt.Printf("UTC date = %s\n\n", formatTime(time))
	fmt.Println("BODY           RA      DEC       AZ      ALT")
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
	fmt.Printf("search   : %s\n", formatTime(time))
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
		root, err := astronomy.SearchRiseSet(ev.body, observer, ev.direction, time, 300.0, 0.0)
		if err != nil {
			fmt.Printf("%-8s : ERROR %v\n", ev.name, err)
			continue
		}
		fmt.Printf("%-8s : %s\n", ev.name, formatTime(root))
	}
	return nil
}

// demoMoonPhase shows the current Moon phase and the next several quarters.
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
	names := []string{"new", "first quarter", "full", "last quarter"}
	fmt.Printf("Moon phase at %s = %.2f degrees\n\n", formatTime(time), phase)
	fmt.Printf("Next 8 quarter phases:\n")
	mq, err := astronomy.SearchMoonQuarter(time)
	if err != nil {
		return err
	}
	for i := 0; i < 8; i++ {
		fmt.Printf("  %-12s %s\n", names[mq.Quarter], formatTime(mq.Time))
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
		return fmt.Errorf("expected <year>")
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
		fmt.Printf("%-18s : %s\n", r.name, formatTime(r.time))
	}
	return nil
}

// demoCulminate finds the next culmination (highest point) of the Sun, Moon,
// and planets as seen by an observer.
func demoCulminate(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	fmt.Printf("search   : %s\n\n", formatTime(time))
	bodies := []astronomy.Body{
		astronomy.Sun, astronomy.Moon, astronomy.Mercury, astronomy.Venus, astronomy.Mars,
		astronomy.Jupiter, astronomy.Saturn, astronomy.Uranus, astronomy.Neptune, astronomy.Pluto,
	}
	for _, body := range bodies {
		_, root, err := astronomy.SearchHourAngle(body, observer, 0.0, time, +1)
		if err != nil {
			fmt.Printf("%-8s : ERROR %v\n", astronomy.BodyName(body), err)
			continue
		}
		fmt.Printf("%-8s : %s\n", astronomy.BodyName(body), formatTime(root))
	}
	return nil
}

// demoLunarEclipse displays the next 10 lunar eclipses after the given time.
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
	eclipse, err := astronomy.SearchLunarEclipse(time)
	if err != nil {
		return err
	}
	names := []string{"none", "penumbral", "partial", "annular", "total"}
	for i := 0; i < 10; i++ {
		fmt.Printf("%2d. %s : %-10s sdPenum=%.1f sdPartial=%.1f sdTotal=%.1f\n",
			i+1, formatTime(eclipse.Peak), names[eclipse.Kind], eclipse.SdPenum, eclipse.SdPartial, eclipse.SdTotal)
		eclipse, err = astronomy.NextLunarEclipse(eclipse.Peak)
		if err != nil {
			return err
		}
	}
	return nil
}

// demoSolarTime demonstrates using the Sun's hour angle to calculate true
// solar time.
func demoSolarTime(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	fmt.Printf("At %s\n", formatTime(time))
	// The true solar time is found from the Sun's hour angle: noon is 12:00.
	_, root, err := astronomy.SearchHourAngle(astronomy.Sun, observer, 12.0, time, +1)
	if err != nil {
		return err
	}
	fmt.Printf("The next time the Sun reaches true solar noon is %s\n", formatTime(root))
	return nil
}

// demoGravity displays the acceleration of gravity at a series of latitudes.
func demoGravity(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("expected <latitude> [height]")
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
	fmt.Printf("latitude=%g height=%g gravity=%0.9f m/s^2\n", lat, height, astronomy.ObserverGravity(lat, height))
	return nil
}

// demoGalactic converts galactic coordinates to horizontal coordinates.
func demoGalactic(args []string) error {
	if len(args) < 4 {
		return fmt.Errorf("expected <lat> <lon> <glat> <glon> [dateTime]")
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
	// Convert galactic coordinates to a vector in the galactic frame, then
	// rotate to the J2000 equatorial frame, then to the horizontal frame.
	galVec := astronomy.VectorFromSphere(astronomy.Spherical{Lat: glat, Lon: glon, Dist: 1.0}, time)
	eqjVec := astronomy.RotateVector(astronomy.Rotation_GAL_EQJ(), galVec)
	horVec := astronomy.RotateVector(astronomy.Rotation_EQJ_HOR(time, observer), eqjVec)
	hor, err := astronomy.HorizonFromVector(horVec, astronomy.NormalRefraction)
	if err != nil {
		return err
	}
	fmt.Printf("galactic (%.3f, %.3f) -> horizontal az=%.3f alt=%.3f\n", glon, glat, hor.Lon, hor.Lat)
	return nil
}

// demoHorizon finds where the ecliptic crosses the observer's horizon.
func demoHorizon(args []string) error {
	observer, time, err := parseObserverAndTime(args)
	if err != nil {
		return err
	}
	// The ecliptic plane is defined by its normal vector in EQJ coordinates.
	eclipticNormal := astronomy.RotateVector(astronomy.Rotation_EQJ_ECL(),
		astronomy.AstroVector{X: 0.0, Y: 0.0, Z: 1.0, T: time})
	// Rotate the ecliptic normal into horizontal coordinates.
	horNormal := astronomy.RotateVector(astronomy.Rotation_EQJ_HOR(time, observer), eclipticNormal)
	fmt.Printf("ecliptic normal in horizontal coordinates: (%.3f, %.3f, %.3f)\n",
		horNormal.X, horNormal.Y, horNormal.Z)
	return nil
}
