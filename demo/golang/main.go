// Package main contains example programs for Astronomy Engine in Go.
package main

import (
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/cosinekitty/astronomy/source/golang"
)

const usage = `Usage: %s <demo> [args...]

Available demos:
  positions  <lat> <lon> [dateTime]   - equatorial and horizontal coordinates
  riseset    <lat> <lon> [dateTime]   - sunrise, sunset, moonrise, moonset
  moonphase  [dateTime]               - Moon phase and next quarters
  seasons    <year>                   - equinoxes and solstices for a year
  culminate  <lat> <lon> [dateTime]   - next culmination of Sun, Moon, planets
  lunar_eclipse [dateTime]            - next 10 lunar eclipses
  solar_time <lat> <lon> [dateTime]   - true solar time
  gravity    <latitude> <height>      - acceleration of gravity
  galactic   <lat> <lon> <glat> <glon> [dateTime] - convert galactic to horizontal
  horizon    <lat> <lon> [dateTime]   - where the ecliptic crosses the horizon
`

// parseObserverAndTime parses a latitude, longitude, and optional UTC time
// from the command line arguments.
func parseObserverAndTime(args []string) (astronomy.Observer, astronomy.AstroTime, error) {
	if len(args) < 2 || len(args) > 3 {
		return astronomy.Observer{}, astronomy.AstroTime{}, errors.New("expected <lat> <lon> [dateTime]")
	}
	lat, err := strconv.ParseFloat(args[0], 64)
	if err != nil {
		return astronomy.Observer{}, astronomy.AstroTime{}, err
	}
	lon, err := strconv.ParseFloat(args[1], 64)
	if err != nil {
		return astronomy.Observer{}, astronomy.AstroTime{}, err
	}
	observer := astronomy.Observer{Latitude: lat, Longitude: lon, Height: 0.0}
	var time astronomy.AstroTime
	if len(args) == 3 {
		time, err = parseTime(args[2])
		if err != nil {
			return astronomy.Observer{}, astronomy.AstroTime{}, err
		}
	} else {
		time = astronomy.CurrentTime()
	}
	return observer, time, nil
}

// parseTime parses a UTC timestamp in the form YYYY-MM-DDTHH:MM:SSZ.
func parseTime(text string) (astronomy.AstroTime, error) {
	fields := strings.SplitN(text, "T", 2)
	if len(fields) != 2 {
		return astronomy.AstroTime{}, errors.New("expected dateTime in the form YYYY-MM-DDTHH:MM:SSZ")
	}
	date := strings.Split(fields[0], "-")
	if len(date) != 3 {
		return astronomy.AstroTime{}, errors.New("invalid date portion")
	}
	clock := strings.TrimSuffix(fields[1], "Z")
	parts := strings.Split(clock, ":")
	if len(parts) != 3 {
		return astronomy.AstroTime{}, errors.New("invalid time portion")
	}
	atoi := func(s string) (int, error) {
		return strconv.Atoi(s)
	}
	year, err := atoi(date[0])
	if err != nil {
		return astronomy.AstroTime{}, err
	}
	month, err := atoi(date[1])
	if err != nil {
		return astronomy.AstroTime{}, err
	}
	day, err := atoi(date[2])
	if err != nil {
		return astronomy.AstroTime{}, err
	}
	hour, err := atoi(parts[0])
	if err != nil {
		return astronomy.AstroTime{}, err
	}
	minute, err := atoi(parts[1])
	if err != nil {
		return astronomy.AstroTime{}, err
	}
	second, err := strconv.ParseFloat(parts[2], 64)
	if err != nil {
		return astronomy.AstroTime{}, err
	}
	return astronomy.TimeFromCalendar(year, month, day, hour, minute, second), nil
}

func main() {
	prog := os.Args[0]
	if len(os.Args) < 2 {
		fmt.Printf(usage, prog)
		os.Exit(1)
	}
	cmd := os.Args[1]
	args := os.Args[2:]
	var err error
	switch cmd {
	case "positions":
		err = demoPositions(args)
	case "riseset":
		err = demoRiseSet(args)
	case "moonphase":
		err = demoMoonPhase(args)
	case "seasons":
		err = demoSeasons(args)
	case "culminate":
		err = demoCulminate(args)
	case "lunar_eclipse":
		err = demoLunarEclipse(args)
	case "solar_time":
		err = demoSolarTime(args)
	case "gravity":
		err = demoGravity(args)
	case "galactic":
		err = demoGalactic(args)
	case "horizon":
		err = demoHorizon(args)
	default:
		fmt.Printf(usage, prog)
		os.Exit(1)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
		os.Exit(1)
	}
}
