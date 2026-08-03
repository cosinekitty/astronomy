// Package astronomy implements [Astronomy Engine] for the Go programming language.
//
// It provides a suite of well-tested functions for calculating positions of the
// Sun, Moon, and planets, along with many practical phenomena visible from observers
// on the Earth, such as sunrise, sunset, seasons, eclipses, transits, and lunar phases.
//
// [Astronomy Engine]: https://github.com/cosinekitty/astronomy
package astronomy

import (
	"errors"
	"math"
)

const (
	DaysPerTropicalYear       = 365.24217                                 // the number of days in one tropical year
	SpeedOfLightAuPerDay      = 173.1446326846693                         // the speed of light in vacuum expressed in astronomical units per day
	KmPerAu                   = 1.4959787069098932e+8                     // the number of kilometers in one astronomical unit
	AuPerLightYear            = 63241.07708807546                         // the number of astronomical units in one light year
	SunRadiusKm               = 695700.0                                  // the radius of the Sun in kilometers
	MercuryEquatorialRadiusKm = 2440.5                                    // the equatorial radius of Mercury in kilometers
	MercuryPolarRadiusKm      = 2438.3                                    // the polar radius of Mercury in kilometers
	VenusRadiusKm             = 6051.8                                    // the radius of Venus in kilometers
	EarthEquatorialRadiusKm   = 6378.1366                                 // the equatorial radius of the Earth in kilometers
	EarthEquatorialRadiusAu   = EarthEquatorialRadiusKm / KmPerAu         //the equatorial radius of the Earth in astronomical units
	EarthFlattening           = 0.996647180302104                         // the ratio of the Earth's polar radius to its equatorial radius
	EarthPolarRadiusKm        = EarthEquatorialRadiusKm * EarthFlattening // the polar radius of the Earth in kilometers
	MoonEquatorialRadiusKm    = 1738.1                                    // the Moon's equatorial radius in kilometers
	MoonPolarRadiusKm         = 1736.0                                    // the Moon's polar radius in kilometers
	MarsEquatorialRadiusKm    = 3396.2                                    // the equatorial radius of Mars in kilometers
	MarsPolarRadiusKm         = 3376.2                                    // the polar radius of Mars in kilometers
	JupiterEquatorialRadiusKm = 71492.0                                   // the equatorial radius of Jupiter in kilometers
	JupiterPolarRadiusKm      = 66854.0                                   // the polar radius of Jupiter in kilometers
	JupiterMeanRadiusKm       = 69911.0                                   // the volumetric mean radius of Jupiter in kilometers
	IoRadiusKm                = 1821.6                                    // the radius of Jupiter's moon Io in kilometers
	EuropaRadiusKm            = 1560.8                                    // the radius of Jupiter's moon Europa in kilometers
	GanymedeRadiusKm          = 2631.2                                    // the radius of Jupiter's moon Ganymede in kilometers
	CallistoRadiusKm          = 2410.3                                    // the radius of Jupiter's moon Callisto in kilometers
	SaturnEquatorialRadiusKm  = 60268.0                                   // the equatorial radius of Saturn in kilometers
	SaturnPolarRadiusKm       = 54364.0                                   // the polar radius of Saturn in kilometers
	UranusEquatorialRadiusKm  = 25559.0                                   // the equatorial radius of Uranus in kilometers
	UranusPolarRadiusKm       = 24973.0                                   // the polar radius of Uranus in kilometers
	NeptuneEquatorialRadiusKm = 24764.0                                   // the equatorial radius of Neptune in kilometers
	NeptunePolarRadiusKm      = 24341.0                                   // the polar radius of Neptune in kilometers
	PlutoRadiusKm             = 1188.3                                    // the radius of Pluto in kilometers
)

const (
	arc       = 3600.0 * 180.0 / math.Pi // arcseconds per radian
	asecToRad = 1.0 / arc                // radians per arcsecond
)

const (
	cosObliq2000       = 0.9174821430670688 // cosine of the mean obliquity of the J2000 ecliptic
	sinObliq2000       = 0.3977769691083922 // sine of the mean obliquity of the J2000 ecliptic
	earthMoonMassRatio = 81.30056           // ratio of Earth mass to Moon mass
	sunGM              = 0.2959122082855911e-03
	mercuryGM          = 0.4912547451450812e-10
	venusGM            = 0.7243452486162703e-09
	earthGM            = 0.8887692390113509e-09
	marsGM             = 0.9549535105779258e-10
	jupiterGM          = 0.2825345909524226e-06
	saturnGM           = 0.8459715185680659e-07
	uranusGM           = 0.1292024916781969e-07
	neptuneGM          = 0.1524358900784276e-07
	plutoGM            = 0.2188699765425970e-11
	moonGM             = earthGM / earthMoonMassRatio
)

const (
	earthOrbitalPeriod   = 365.256 // days for Earth to orbit the Sun
	neptuneOrbitalPeriod = 60189.0 // days for Neptune to orbit the Sun
)

func isfinite(x float64) bool {
	return !math.IsInf(x, 0) && !math.IsNaN(x)
}

type AstroTime struct {
	// Ut holds the floating point number of days of Universal Time since noon UTC January 1, 2000.
	// Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
	// not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds.
	// This approximation is sufficient for the accuracy requirements of Astronomy Engine.
	//
	// Universal Time Coordinate (UTC) is the international standard for legal and civil
	// timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
	// UTC is kept in sync with unpredictable observed changes in the Earth's rotation
	// by occasionally adding leap seconds as needed.
	//
	// UT1 is an idealized time scale based on observed rotation of the Earth, which
	// gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
	// large scale weather events like hurricanes, and internal seismic and convection effects.
	// Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
	// is adjusted by a scheduled whole number of leap seconds as needed.
	//
	// The value in Ut is appropriate for any calculation involving the Earth's rotation,
	// such as calculating rise/set times, culumination, and anything involving apparent
	// sidereal time.
	//
	// Before the era of atomic timekeeping, days based on the Earth's rotation
	// were often known as ``mean solar days''.
	Ut float64

	// Tt holds a Terrestrial Time value expressed in days.
	// Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
	// In this system, days are not based on Earth rotations, but instead by
	// the number of elapsed [SI seconds] divided by 86400.
	// Unlike Ut, Tt increases uniformly without adjustments for changes in the Earth's rotation.
	//
	// The value in Tt is used for calculations of movements not involving the Earth's rotation,
	// such as the orbits of planets around the Sun, or the Moon around the Earth.
	//
	// Historically, Terrestrial Time has also been known by the term ``Ephemeris Time'' (ET).
	//
	// [SI seconds]: https://physics.nist.gov/cuu/Units/second.html
	Tt float64

	psi float64 // For internal use only. Used to optimize Earth tilt calculations.
	eps float64 // For internal use only.  Used to optimize Earth tilt calculations.
	st  float64 // For internal use only.  Lazy-caches sidereal time (Earth rotation).
}

func makeTime(ut float64, tt float64) AstroTime {
	return AstroTime{
		Ut:  ut,
		Tt:  tt,
		psi: math.NaN(),
		eps: math.NaN(),
		st:  math.NaN(),
	}
}

func terrestrialTime(ut float64) float64 {
	return ut + deltaTime(ut)/86400.0
}

func universalTime(tt float64) float64 {
	// This is the inverse function of terrestrialTime.
	// This is an iterative numerical solver, but because
	// the relationship between UT and TT is almost perfectly linear,
	// it converges extremely fast (never more than 3 iterations).
	var dt = terrestrialTime(tt) - tt
	for {
		var ut = tt - dt
		var ttCheck = terrestrialTime(ut)
		var err = ttCheck - tt
		if math.Abs(err) < 1.0e-12 {
			return ut
		}
		dt += err
	}
}

func deltaTime(ut float64) float64 {
	// Fred Espenak writes about Delta-T generically here:
	// https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
	// https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html
	//
	// He provides polynomial approximations for distant years here:
	// https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
	//
	// They start with a year value 'y' such that y=2000 corresponds
	// to the UTC Date 15-January-2000. Convert difference in days
	// to mean tropical years.

	var y, u, u2, u3, u4, u5, u6, u7 float64

	y = 2000 + ((ut - 14) / DaysPerTropicalYear)

	if y < -500 {
		u = (y - 1820) / 100
		return -20 + (32 * u * u)
	}
	if y < 500 {
		u = y / 100
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		u5 = u2 * u3
		u6 = u3 * u3
		return 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3 - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6
	}
	if y < 1600 {
		u = (y - 1000) / 100
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		u5 = u2 * u3
		u6 = u3 * u3
		return 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3 - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6
	}
	if y < 1700 {
		u = y - 1600
		u2 = u * u
		u3 = u * u2
		return 120 - 0.9808*u - 0.01532*u2 + u3/7129.0
	}
	if y < 1800 {
		u = y - 1700
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		return 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000
	}
	if y < 1860 {
		u = y - 1800
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		u5 = u2 * u3
		u6 = u3 * u3
		u7 = u3 * u4
		return 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4 + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7
	}
	if y < 1900 {
		u = y - 1860
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		u5 = u2 * u3
		return 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3 - 0.0004473624*u4 + u5/233174
	}
	if y < 1920 {
		u = y - 1900
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		return -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4
	}
	if y < 1941 {
		u = y - 1920
		u2 = u * u
		u3 = u * u2
		return 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3
	}
	if y < 1961 {
		u = y - 1950
		u2 = u * u
		u3 = u * u2
		return 29.07 + 0.407*u - u2/233 + u3/2547
	}
	if y < 1986 {
		u = y - 1975
		u2 = u * u
		u3 = u * u2
		return 45.45 + 1.067*u - u2/260 - u3/718
	}
	if y < 2005 {
		u = y - 2000
		u2 = u * u
		u3 = u * u2
		u4 = u2 * u2
		u5 = u2 * u3
		return 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4 + 0.00002373599*u5
	}
	if y < 2050 {
		u = y - 2000
		return 62.92 + 0.32217*u + 0.005589*u*u
	}
	if y < 2150 {
		u = (y - 1820) / 100
		return -20 + 32*u*u - 0.5628*(2150-y)
	}

	// all years after 2150
	u = (y - 1820) / 100
	return -20 + (32 * u * u)
}

// TimeFromUniversalDays converts a UTC number of days since January 1, 2000
// into an AstroTime value that can be used for astronomy calculations.
func TimeFromUniversalDays(ut float64) AstroTime {
	return makeTime(ut, terrestrialTime(ut))
}

// TimeFromTerrestrialDays converts a Terrestrial Time (TT) day value,
// also known as ephemeris days, to an AstroTime value that can be used for astronomy calculations.
func TimeFromTerrestrialDays(tt float64) AstroTime {
	return makeTime(universalTime(tt), tt)
}

// Given an AstroTime value, creates a new AstroTime value that is the
// specified number of UT days in the future (positive) or past (negative).
func (time AstroTime) AddDays(days float64) AstroTime {
	return TimeFromUniversalDays(time.Ut + days)
}

// CalendarDateTime represents a Gregorian calendar date and time within
// plus or minus 1 million years from the year 0.
type CalendarDateTime struct {
	Year   int     // The year value in the range -999999 to +999999.
	Month  int     // The calendar month in the range 1..12.
	Day    int     // The day of the month in the range 1..31.
	Hour   int     // The hour in the range 0..23.
	Minute int     // The minute in the range 0..59.
	Second float64 // The real-valued second in the half-open range [0, 60).
}

// CalendarFromDays converts a J2000 day value to a Gregorian calendar date and time.
func CalendarFromDays(ut float64) (*CalendarDateTime, error) {
	if !isfinite(ut) {
		return nil, errors.New("Non-finite value of ut")
	}

	// Adapted from the NOVAS C 3.1 function cal_date().
	// Convert fractional days since J2000 into Gregorian calendar date/time.
	djd := ut + 2451545.5
	jd := int64(math.Floor(djd))
	var x float64
	x = 24.0 * math.Mod(djd, 1.0)
	if x < 0.0 {
		x += 24.0
	}
	hour := int(x)
	x = 60.0 * math.Mod(x, 1.0)
	minute := int(x)
	second := 60.0 * math.Mod(x, 1.0)

	// This is my own adjustment to the NOVAS cal_date logic
	// so that it can handle dates much farther back in the past.
	// I add c*400 years worth of days at the front,
	// then subtract c*400 years at the back,
	// which avoids negative values in the formulas that mess up
	// the calendar date calculations.
	// Any multiple of 400 years has the same number of days,
	// because it eliminates all the special cases for leap years.
	const c = 2500

	var k, n, m int64
	k = jd + (68569 + c*146097)
	n = (4 * k) / 146097
	k = k - (146097*n+3)/4
	m = (4000 * (k + 1)) / 1461001
	k = k - (1461*m)/4 + 31

	month := (int)((80 * k) / 2447)
	day := (int)(k - int64(2447*month)/80)
	k = int64(month / 11)

	month = int(int64(month+2) - 12*k)
	year := int(100*(n-49) + m + k - 400*c)

	if year < -999999 || year > +999999 {
		return nil, errors.New("The supplied time is too far from the year 2000 to be represented.")
	}

	if month < 1 || month > 12 || day < 1 || day > 31 {
		return nil, errors.New("Internal error: invalid calendar date calculated.")
	}

	return &CalendarDateTime{year, month, day, hour, minute, second}, nil
}

func DaysFromCalendar(year, month, day, hour, minute int, second float64) float64 {
	// This formula is adapted from NOVAS C 3.1 function julian_date(),
	// which in turn comes from Henry F. Fliegel & Thomas C. Van Flendern:
	// Communications of the ACM, Vol 11, No 10, October 1968, p. 657.
	// See: https://dl.acm.org/doi/pdf/10.1145/364096.364097
	//
	// [Don Cross - 2023-02-25] I modified the formula so that it will
	// work correctly with years as far back as -999999.
	y := int64(year)
	m := int64(month)
	d := int64(day)
	f := (14 - m) / 12
	y2000 := (d - 365972956) + (1461*(y+1000000-f))/4 + (367*(m-2+12*f))/12 - (3*((y+1000100-f)/100))/4
	ut := (float64(y2000) - 0.5) + (float64(hour) / 24.0) + (float64(minute) / 1440.0) + (second / 86400.0)
	return ut
}

// TimeFromCalendar returns an AstroTime value for a date and time expressed in civil UTC.
func TimeFromCalendar(year, month, day, hour, minute int, second float64) AstroTime {
	ut := DaysFromCalendar(year, month, day, hour, minute, second)
	return TimeFromUniversalDays(ut)
}

// Atmosphere calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.
// elevationMeters is the elevation above sea level at which to calculate atmospheric variables.
// It must be in the range -500 to +100000 or an error will occur.
// 1. COESA, U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, DC, 1976.
// 2. Jursa, A. S., Ed., Handbook of Geophysics and the Space Environment, Air Force Geophysics Laboratory, 1985.
// See:
// https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml
// https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
// https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf
func Atmosphere(elevationMeters float64) (AtmosphereInfo, error) {
	if !isfinite(elevationMeters) || elevationMeters < -500.0 || elevationMeters > 100000.0 {
		return AtmosphereInfo{}, errors.New("Invalid elevation")
	}
	const P0 = 101325.0 // pressure at sea level [pascals]
	const T0 = 288.15   // temperature at sea level [kelvins]
	const T1 = 216.65   // temperature between 20 km and 32 km [kelvins]

	var temperature, pressure float64
	switch {
	case elevationMeters <= 11000.0:
		{
			temperature = T0 - 0.0065*elevationMeters
			pressure = P0 * math.Pow(T0/temperature, -5.25577)
		}
	case elevationMeters <= 20000.0:
		{
			temperature = T1
			pressure = 22632.0 * math.Exp(-0.00015768832*(elevationMeters-11000.0))
		}
	default:
		{
			temperature = T1 + 0.001*(elevationMeters-20000.0)
			pressure = 5474.87 * math.Pow(T1/temperature, 34.16319)
		}
	}
	// The density is calculated relative to the sea level value.
	// Using the ideal gas law PV=nRT, we deduce that density is proportional to P/T.
	density := (pressure / temperature) / (P0 / T0)
	return AtmosphereInfo{pressure, temperature, density}, nil
}

// AstroVector represents a position in 3D space at a given time.
// Usually the distance components are expressed in astronomical units (AU).
// The origin and orientation system depends on context.
// Occasionally AstroVector is used to represent a velocity vector,
// in which case the component units are astronomical units per day.
type AstroVector struct {
	X float64
	Y float64
	Z float64
	T AstroTime
}

// Returns the length of vec expressed in the same distance units as vec's components.
func (vec AstroVector) Length() float64 {
	return math.Sqrt(vec.X*vec.X + vec.Y*vec.Y + vec.Z*vec.Z)
}

// StateVector represents the combined position and velocity of a body at a given moment of time.
type StateVector struct {
	X  float64
	Y  float64
	Z  float64
	Vx float64
	Vy float64
	Vz float64
	T  AstroTime
}

// Position returns the position vector inside a state vector.
func (state StateVector) Position() AstroVector {
	return AstroVector{state.X, state.Y, state.Z, state.T}
}

// Position returns the velocity vector inside a state vector.
func (state StateVector) Velocity() AstroVector {
	return AstroVector{state.Vx, state.Vy, state.Vz, state.T}
}

// Spherical coordinates for a body in space
type Spherical struct {
	Lat  float64 // a latitude-like angle expressed in degrees
	Lon  float64 // a longitude-like angle expressed in degrees
	Dist float64 // a distance expressed in astronomical units
}

type Body int

const (
	InvalidBody Body = -1
	Mercury     Body = 0   // The planet Mercury
	Venus            = 1   // The planet Venus
	Earth            = 2   // The planet Earth
	Mars             = 3   // The planet Mars
	Jupiter          = 4   // The planet Jupiter
	Saturn           = 5   // The planet Saturn
	Uranus           = 6   // The planet Uranus
	Neptune          = 7   // The planet Neptune
	Pluto            = 8   // The dwarf planet Pluto
	Sun              = 9   // The Sun
	Moon             = 10  // The Earth's Moon
	Emb              = 11  // The Earth/Moon Barycenter
	Ssb              = 12  // The Solar System Barycenter
	Star1            = 101 // User-defined star #1
	Star2            = 102 // User-defined star #2
	Star3            = 103 // User-defined star #3
	Star4            = 104 // User-defined star #4
	Star5            = 105 // User-defined star #5
	Star6            = 106 // User-defined star #6
	Star7            = 107 // User-defined star #7
	Star8            = 108 // User-defined star #8
)

// stardef stores the equatorial coordinates and distance of a user-defined star.
type stardef struct {
	ra   float64 // right ascension in sidereal hours (J2000 equatorial coordinates)
	dec  float64 // declination in degrees (J2000 equatorial coordinates)
	dist float64 // distance from the Sun in astronomical units
}

// starTable holds the definitions of the eight user-defined stars (Star1..Star8).
var starTable [8]stardef

// getStarPointer returns the definition slot for a user-defined star, or nil
// if the body is not one of the eight user-defined stars.
func getStarPointer(body Body) *stardef {
	if body >= Star1 && body <= Star8 {
		return &starTable[body-Star1]
	}
	return nil
}

// userDefinedStar returns the definition of a user-defined star, or nil if
// the body is not a star or has not been defined yet.
func userDefinedStar(body Body) *stardef {
	star := getStarPointer(body)
	if star != nil && star.dist > 0.0 {
		return star
	}
	return nil
}

// DefineStar assigns equatorial coordinates to a user-defined star.
// Some Astronomy Engine functions allow their body parameter to be a
// user-defined fixed point in the sky, loosely called a "star".
// This function assigns a right ascension, declination, and distance to one of
// the eight user-defined stars Star1 .. Star8. Stars are not valid until defined.
func DefineStar(body Body, ra float64, dec float64, distanceLightYears float64) error {
	star := getStarPointer(body)
	if star == nil {
		return errors.New("Invalid star body")
	}
	if !isfinite(ra) || ra < 0.0 || ra >= 24.0 {
		return errors.New("Invalid right ascension for star")
	}
	if !isfinite(dec) || dec < -90.0 || dec > +90.0 {
		return errors.New("Invalid declination for star")
	}
	if !isfinite(distanceLightYears) || distanceLightYears < 1.0 {
		return errors.New("Invalid distance for star")
	}
	star.ra = ra
	star.dec = dec
	star.dist = distanceLightYears * AuPerLightYear
	return nil
}

// BodyName returns the English name of the given body, or "" if invalid.
func BodyName(body Body) string {
	switch body {
	case Mercury:
		return "Mercury"
	case Venus:
		return "Venus"
	case Earth:
		return "Earth"
	case Mars:
		return "Mars"
	case Jupiter:
		return "Jupiter"
	case Saturn:
		return "Saturn"
	case Uranus:
		return "Uranus"
	case Neptune:
		return "Neptune"
	case Pluto:
		return "Pluto"
	case Sun:
		return "Sun"
	case Moon:
		return "Moon"
	case Emb:
		return "EMB"
	case Ssb:
		return "SSB"
	default:
		return ""
	}
}

// BodyCode returns the Body value corresponding to the given English name
// (case-sensitive), or InvalidBody if the name is not recognized.
func BodyCode(name string) Body {
	switch name {
	case "Mercury":
		return Mercury
	case "Venus":
		return Venus
	case "Earth":
		return Earth
	case "Mars":
		return Mars
	case "Jupiter":
		return Jupiter
	case "Saturn":
		return Saturn
	case "Uranus":
		return Uranus
	case "Neptune":
		return Neptune
	case "Pluto":
		return Pluto
	case "Sun":
		return Sun
	case "Moon":
		return Moon
	case "EMB":
		return Emb
	case "SSB":
		return Ssb
	default:
		return InvalidBody
	}
}

// isSuperiorPlanet returns true for planets that are farther from the Sun
// than the Earth is.
func isSuperiorPlanet(body Body) bool {
	switch body {
	case Mars, Jupiter, Saturn, Uranus, Neptune, Pluto:
		return true
	default:
		return false
	}
}

// PlanetOrbitalPeriod returns the average number of days it takes for a planet
// to orbit the Sun. Returns 0.0 if body is not one of the nine planets.
func PlanetOrbitalPeriod(body Body) float64 {
	switch body {
	case Mercury:
		return 87.969
	case Venus:
		return 224.701
	case Earth:
		return earthOrbitalPeriod
	case Mars:
		return 686.980
	case Jupiter:
		return 4332.589
	case Saturn:
		return 10759.22
	case Uranus:
		return 30685.4
	case Neptune:
		return neptuneOrbitalPeriod
	case Pluto:
		return 90560.0
	default:
		return 0.0 // invalid body
	}
}

// MassProduct returns the product of the body's mass and the gravitational
// constant G, expressed in astronomical units, days, and solar masses.
// Returns 0.0 for invalid bodies and for the Earth/Moon barycenter's
// component bodies.
func MassProduct(body Body) float64 {
	switch body {
	case Sun:
		return sunGM
	case Mercury:
		return mercuryGM
	case Venus:
		return venusGM
	case Earth:
		return earthGM
	case Moon:
		return moonGM
	case Emb:
		return earthGM + moonGM
	case Mars:
		return marsGM
	case Jupiter:
		return jupiterGM
	case Saturn:
		return saturnGM
	case Uranus:
		return uranusGM
	case Neptune:
		return neptuneGM
	case Pluto:
		return plutoGM
	default:
		return 0.0 // invalid body
	}
}

// bodyState holds the position and velocity of a body at a moment of time.
type bodyState struct {
	tt float64     // Terrestrial Time in J2000 days
	r  terseVector // position [au]
	v  terseVector // velocity [au/day]
}

// exportState converts a terse bodyState into a full StateVector at the given time.
func exportState(terse bodyState, time AstroTime) StateVector {
	return StateVector{
		X:  terse.r.X,
		Y:  terse.r.Y,
		Z:  terse.r.Z,
		Vx: terse.v.X,
		Vy: terse.v.Y,
		Vz: terse.v.Z,
		T:  time,
	}
}

// updatePosition advances a position by dt days using the current velocity and
// acceleration, integrating acceleration by the midpoint method.
func updatePosition(dt float64, r terseVector, v terseVector, a terseVector) terseVector {
	r.X += (v.X + a.X*dt/2) * dt
	r.Y += (v.Y + a.Y*dt/2) * dt
	r.Z += (v.Z + a.Z*dt/2) * dt
	return r
}

// updateVelocity advances a velocity by dt days using the given acceleration.
func updateVelocity(dt float64, v terseVector, a terseVector) terseVector {
	v.X += dt * a.X
	v.Y += dt * a.Y
	v.Z += dt * a.Z
	return v
}

// Rotation_EQJ_ECL returns a rotation matrix that converts EQJ (J2000 mean
// equator) coordinates to ECL (J2000 mean ecliptic) coordinates.
func Rotation_EQJ_ECL() RotationMatrix {
	var r RotationMatrix
	r.Rot[0][0] = 1.0
	r.Rot[0][1] = 0.0
	r.Rot[0][2] = 0.0
	r.Rot[1][0] = 0.0
	r.Rot[1][1] = +cosObliq2000
	r.Rot[1][2] = -sinObliq2000
	r.Rot[2][0] = 0.0
	r.Rot[2][1] = +sinObliq2000
	r.Rot[2][2] = +cosObliq2000
	return r
}

// Rotation_ECL_EQJ returns a rotation matrix that converts ECL (J2000 mean
// ecliptic) coordinates to EQJ (J2000 mean equator) coordinates.
func Rotation_ECL_EQJ() RotationMatrix {
	var r RotationMatrix
	r.Rot[0][0] = 1.0
	r.Rot[0][1] = 0.0
	r.Rot[0][2] = 0.0
	r.Rot[1][0] = 0.0
	r.Rot[1][1] = +cosObliq2000
	r.Rot[1][2] = +sinObliq2000
	r.Rot[2][0] = 0.0
	r.Rot[2][1] = -sinObliq2000
	r.Rot[2][2] = +cosObliq2000
	return r
}

// Rotation_EQJ_GAL returns a rotation matrix that converts EQJ (J2000 mean
// equator) coordinates to GAL (galactic, IAU 1958) coordinates.
func Rotation_EQJ_GAL() RotationMatrix {
	var rot RotationMatrix
	rot.Rot[0][0] = -0.0548624779711344
	rot.Rot[0][1] = +0.4941095946388765
	rot.Rot[0][2] = -0.8676668813529025
	rot.Rot[1][0] = -0.8734572784246782
	rot.Rot[1][1] = -0.4447938112296831
	rot.Rot[1][2] = -0.1980677870294097
	rot.Rot[2][0] = -0.4838000529948520
	rot.Rot[2][1] = +0.7470034631630423
	rot.Rot[2][2] = +0.4559861124470794
	return rot
}

// Rotation_GAL_EQJ returns a rotation matrix that converts GAL (galactic,
// IAU 1958) coordinates to EQJ (J2000 mean equator) coordinates.
func Rotation_GAL_EQJ() RotationMatrix {
	var rot RotationMatrix
	rot.Rot[0][0] = -0.0548624779711344
	rot.Rot[0][1] = -0.8734572784246782
	rot.Rot[0][2] = -0.4838000529948520
	rot.Rot[1][0] = +0.4941095946388765
	rot.Rot[1][1] = -0.4447938112296831
	rot.Rot[1][2] = +0.7470034631630423
	rot.Rot[2][0] = -0.8676668813529025
	rot.Rot[2][1] = -0.1980677870294097
	rot.Rot[2][2] = +0.4559861124470794
	return rot
}

// The location of a point on or near the surface of the Earth
type Observer struct {
	Latitude  float64 // Latitude degrees north (positive) or south (negative) of the equator
	Longitude float64 // Longitude east (positive) or west (negative) of the prime meridian passing through Greenwich, England
	Height    float64 // Height above mean sea level in meters
}

// A location of a body expressed in angular coordinates relative to the Earth's equator
type Equatorial struct {
	Ra   float64     // Right Ascension in sidereal hours, in the half-open range [0, 24)
	Dec  float64     // Declination in degrees north (positive) or south (negative) of the celestial equator, in the closed range [-90, +90].
	Dist float64     // Distance of an object in astronomical units [AU]
	Vec  AstroVector // The position expressed as a Cartesian vector in the same equatorial orientation system
}

// A location of a body expressed in angular coordinates relative to the plane of the Earth's orbit around the Sun
type Ecliptic struct {
	Vec  AstroVector // The object position expressed as a Cartesian vector in the same ecliptic orientation system
	Elat float64     // Eclilptic latitude in degrees north (positive) or south (negative) with respect to the ecliptic plane, in the closed range [-90, +90].
	Elon float64     // Ecliptic longitude in degrees, in the half-open range [0, 360).
}

// A location of a body as seen from an observer's point of view on or near the surface of the Earth.
// The topocentric position can optionally be corrected for atmospheric refraction.
type Topocentric struct {
	Azimuth  float64 // The compass direction of the object, where north = 0, east = 90, south = 180, and west = 270.
	Altitude float64 // The angle above (positive) or below (negative) the observer's horizon in degrees, in the range [-90, +90].
	Ra       float64 // The body's apparent right ascension
	Dec      float64 // The body's apparent declination
}

// RotationMatrix is a 3x3 matrix used to convert a vector from one orientation system to another.
type RotationMatrix struct {
	Rot [3][3]float64
}

type Refraction int

const (
	NoRefraction Refraction = iota
	NormalRefraction
	JplHorizonsRefraction
)

// AtmosphereInfo contains information about idealized atmospheric variables at a given elevation.
type AtmosphereInfo struct {
	Pressure    float64 // atmospheric pressure in pascals
	Temperature float64 // atmospheric temperature in kelvins
	Density     float64 // atmospheric density relative to sea level
}

type SeasonsInfo struct {
	MarEquinox  AstroTime
	JunSolstice AstroTime
	SepEquinox  AstroTime
	DecSolstice AstroTime
}

type AstroMoonQuarter struct {
	Quarter int
	Time    AstroTime
}

type TimeFormat int

const (
	TimeFormatDay TimeFormat = iota
	TimeFormatMinute
	TimeFormatSecond
	TimeFormatMilli
)

type AstroSearchFunc func(context interface{}, time AstroTime) float64

type DeltaTimeFunc func(ut float64) float64

type LibrationInfo struct {
	Elat    float64
	Elon    float64
	Mlat    float64
	Mlon    float64
	DistKm  float64
	DiamDeg float64
}

type AxisInfo struct {
	Ra    float64
	Dec   float64
	Spin  float64
	North AstroVector
}

type NodeEventKind int

const (
	InvalidNode    NodeEventKind = 0
	AscendingNode                = 1
	DescendingNode               = -1
)

type NodeEventInfo struct {
	Time AstroTime
	Kind NodeEventKind
}

type terseVector struct {
	X float64
	Y float64
	Z float64
}

type bodyGravCalc struct {
	Tt float64     // J2000 terrestrial time [days]
	R  terseVector // Position [au]
	V  terseVector // Velocity [au/day]
	A  terseVector // Acceleration [au/day^2]
}

type JupiterMoonsInfo struct {
	Io       StateVector
	Europa   StateVector
	Ganymede StateVector
	Callisto StateVector
}

type constelInfo struct {
	symbol string
	name   string
}

type constelBoundary struct {
	index int
	raLo  float64
	raHi  float64
	decLo float64
}

// DegreesFromRadians converts an angle expressed in radians to an angle expressed in degrees.
func DegreesFromRadians(radians float64) float64 {
	return radians * (180.0 / math.Pi)
}

// RadiansFromDegrees converts an angle expressed in degrees to an angle expressed in radians.
func RadiansFromDegrees(degrees float64) float64 {
	return degrees * (math.Pi / 180.0)
}

func dcos(degrees float64) float64 {
	return math.Cos(RadiansFromDegrees(degrees))
}

func dsin(degrees float64) float64 {
	return math.Sin(RadiansFromDegrees(degrees))
}

// AngleBetween calculates the angle in degrees between two vectors.
// Given a pair of vectors avec and bvec, this function returns the
// angle in degrees between the vectors in 3D space.
// The angle is measured in the plane that contains both vectors.
// The returned value is in the closed range [0, 180].
func AngleBetween(avec AstroVector, bvec AstroVector) float64 {
	r := avec.Length() * bvec.Length()
	if r < 1.0e-8 {
		panic("Cannot find angle between vectors because they are too short.")
	}
	dot := (avec.X*bvec.X + avec.Y*bvec.Y + avec.Z*bvec.Z) / r
	if dot <= -1.0 {
		return 180.0
	}
	if dot >= +1.0 {
		return 0.0
	}
	return DegreesFromRadians(math.Acos(dot))
}

func eclOblToEquVec(ecl AstroVector, obl float64) AstroVector {
	c := math.Cos(obl)
	s := math.Sin(obl)
	return AstroVector{
		ecl.X,
		ecl.Y*c - ecl.Z*s,
		ecl.Y*s + ecl.Z*c,
		ecl.T,
	}
}

func meanObliq(tt float64) float64 {
	t := tt / 36525.0
	asec := ((((-0.0000000434*t-0.000000576)*t+0.00200340)*t-0.0001831)*t-46.836769)*t + 84381.406
	return asec / 3600.0
}

func iau2000b(time *AstroTime) {
	// Earth axis nutation function.
	// Adapted from the NOVAS C 3.1 function of the same name.
	// Lazy-evaluates nutation angles only if not yet calculated.
	if math.IsNaN(time.psi) {
		t := time.Tt / 36525.0
		const asec = 360.0 * 3600.0
		elp := math.Mod(1287104.79305+t*129596581.0481, asec) * asecToRad
		f := math.Mod(335779.526232+t*1739527262.8478, asec) * asecToRad
		d := math.Mod(1072260.70369+t*1602961601.2090, asec) * asecToRad
		om := math.Mod(450160.398036-t*6962890.5431, asec) * asecToRad

		sarg := math.Sin(om)
		carg := math.Cos(om)
		dp := (-172064161.0-174666.0*t)*sarg + 33386.0*carg
		de := (92052331.0+9086.0*t)*carg + 15377.0*sarg

		arg := 2.0 * (f - d + om)
		sarg = math.Sin(arg)
		carg = math.Cos(arg)
		dp += (-13170906.0-1675.0*t)*sarg - 13696.0*carg
		de += (5730336.0-3015.0*t)*carg - 4587.0*sarg

		arg = 2.0 * (f + om)
		sarg = math.Sin(arg)
		carg = math.Cos(arg)
		dp += (-2276413.0-234.0*t)*sarg + 2796.0*carg
		de += (978459.0-485.0*t)*carg + 1374.0*sarg

		arg = 2.0 * om
		sarg = math.Sin(arg)
		carg = math.Cos(arg)
		dp += (2074554.0+207.0*t)*sarg - 698.0*carg
		de += (-897492.0+470.0*t)*carg - 291.0*sarg

		sarg = math.Sin(elp)
		carg = math.Cos(elp)
		dp += (1475877.0-3633.0*t)*sarg + 11817.0*carg
		de += (73871.0-184.0*t)*carg - 1924.0*sarg

		time.psi = -0.000135 + (dp * 1.0e-7)
		time.eps = +0.000388 + (de * 1.0e-7)
	}
}

type earthTiltInfo struct {
	Tt   float64
	Dpsi float64
	Deps float64
	Ee   float64
	Mobl float64
	Tobl float64
}

func etilt(time *AstroTime) earthTiltInfo {
	iau2000b(time)
	mobl := meanObliq(time.Tt)
	tobl := mobl + (time.eps / 3600.0)
	ee := time.psi * dcos(mobl) / 15.0
	return earthTiltInfo{
		time.Tt,
		time.psi,
		time.eps,
		ee,
		mobl,
		tobl,
	}
}

// Earth Rotation Angle
func era(ut float64) float64 {
	thet1 := 0.7790572732640 + 0.00273781191135448*ut
	thet3 := math.Mod(ut, 1.0)
	theta := 360.0 * math.Mod(thet1+thet3, 1.0)
	if theta < 0.0 {
		theta += 360.0
	}
	return theta
}

func eclToEquVec(ecl AstroVector) AstroVector {
	oblDeg := meanObliq(ecl.T.Tt)
	oblRad := RadiansFromDegrees(oblDeg)
	return eclOblToEquVec(ecl, oblRad)
}

// Given a date and time, SiderealTime calculates the rotation of the
// Earth, represented by the equatorial angle of the Greenwich prime meridian
// with respect to distant stars (not the Sun, which moves relative to background
// stars by almost one degree per day).
// This angle is called Greenwich Apparent Sidereal Time (GAST).
// GAST is measured in sidereal hours in the half-open range [0, 24).
// When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
// corrected at that time for precession and nutation of the Earth's axis.
// In this context, the "equinox" is the direction in space where the Earth's
// orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
// at the location on the Earth's orbit of the (seasonal) March equinox.
// As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
// then starts over at 0.
// To convert to degrees, multiply the return value by 15.
// As an optimization, this function caches the sidereal time value in the time parameter.
// The value is reused later as needed, to avoid redundant calculations.
func SiderealTime(time *AstroTime) float64 {
	if math.IsNaN(time.st) {
		t := time.Tt / 36525.0
		eqeq := 15.0 * etilt(time).Ee // Replace with eqeq=0 to get GMST instead of GAST (if we ever need it)
		theta := era(time.Ut)
		st := (eqeq + 0.014506 + ((((-0.0000000368*t-0.000029956)*t-0.00000044)*t+1.3915817)*t+4612.156534)*t)
		gst := math.Mod(st/3600.0+theta, 360.0) / 15.0
		if gst < 0.0 {
			gst += 24.0
		}
		time.st = gst
	}
	return time.st
}

type precessDirection int

const (
	from2000 precessDirection = 0
	into2000
)

func precessionRot(time AstroTime, dir precessDirection) RotationMatrix {
	t := time.Tt / 36525.0

	var eps0, psia, omegaa, chia float64

	eps0 = 84381.406

	psia = (((((-0.0000000951*t+0.000132851)*t-0.00114045)*t-1.0790069)*t + 5038.481507) * t)

	omegaa = (((((+0.0000003337*t-0.000000467)*t-0.00772503)*t+0.0512623)*t-0.025754)*t + eps0)

	chia = (((((-0.0000000560*t+0.000170663)*t-0.00121197)*t-2.3814292)*t + 10.556403) * t)

	eps0 *= asecToRad
	psia *= asecToRad
	omegaa *= asecToRad
	chia *= asecToRad

	sa := math.Sin(eps0)
	ca := math.Cos(eps0)
	sb := math.Sin(-psia)
	cb := math.Cos(-psia)
	sc := math.Sin(-omegaa)
	cc := math.Cos(-omegaa)
	sd := math.Sin(chia)
	cd := math.Cos(chia)

	xx := cd*cb - sb*sd*cc
	yx := cd*sb*ca + sd*cc*cb*ca - sa*sd*sc
	zx := cd*sb*sa + sd*cc*cb*sa + ca*sd*sc
	xy := -sd*cb - sb*cd*cc
	yy := -sd*sb*ca + cd*cc*cb*ca - sa*cd*sc
	zy := -sd*sb*sa + cd*cc*cb*sa + ca*cd*sc
	xz := sb * sc
	yz := -sc*cb*ca - sa*cc
	zz := -sc*cb*sa + cc*ca

	rot := RotationMatrix{}
	switch {
	case dir == into2000:
		{
			// Perform rotation from other epoch to J2000.0.
			rot.Rot[0][0] = xx
			rot.Rot[0][1] = yx
			rot.Rot[0][2] = zx
			rot.Rot[1][0] = xy
			rot.Rot[1][1] = yy
			rot.Rot[1][2] = zy
			rot.Rot[2][0] = xz
			rot.Rot[2][1] = yz
			rot.Rot[2][2] = zz
		}
	case dir == from2000:
		{
			// Perform rotation from J2000.0 to other epoch.
			rot.Rot[0][0] = xx
			rot.Rot[0][1] = xy
			rot.Rot[0][2] = xz
			rot.Rot[1][0] = yx
			rot.Rot[1][1] = yy
			rot.Rot[1][2] = yz
			rot.Rot[2][0] = zx
			rot.Rot[2][1] = zy
			rot.Rot[2][2] = zz

		}
	default:
		panic("Unsupported precession direction")
	}
	return rot
}

func precession(pos AstroVector, dir precessDirection) AstroVector {
	r := precessionRot(pos.T, dir)
	return RotateVector(r, pos)
}

// RotateVector applies a rotation to a vector, yielding a vector in another orientation system.
func RotateVector(rotation RotationMatrix, vector AstroVector) AstroVector {
	return AstroVector{
		rotation.Rot[0][0]*vector.X + rotation.Rot[1][0]*vector.Y + rotation.Rot[2][0]*vector.Z,
		rotation.Rot[0][1]*vector.X + rotation.Rot[1][1]*vector.Y + rotation.Rot[2][1]*vector.Z,
		rotation.Rot[0][2]*vector.X + rotation.Rot[1][2]*vector.Y + rotation.Rot[2][2]*vector.Z,
		vector.T,
	}
}

// RotateState applies a rotation to a state vector, rotating both its position
// and velocity components into another orientation system.
func RotateState(rotation RotationMatrix, state StateVector) StateVector {
	return StateVector{
		X:  rotation.Rot[0][0]*state.X + rotation.Rot[1][0]*state.Y + rotation.Rot[2][0]*state.Z,
		Y:  rotation.Rot[0][1]*state.X + rotation.Rot[1][1]*state.Y + rotation.Rot[2][1]*state.Z,
		Z:  rotation.Rot[0][2]*state.X + rotation.Rot[1][2]*state.Y + rotation.Rot[2][2]*state.Z,
		Vx: rotation.Rot[0][0]*state.Vx + rotation.Rot[1][0]*state.Vy + rotation.Rot[2][0]*state.Vz,
		Vy: rotation.Rot[0][1]*state.Vx + rotation.Rot[1][1]*state.Vy + rotation.Rot[2][1]*state.Vz,
		Vz: rotation.Rot[0][2]*state.Vx + rotation.Rot[1][2]*state.Vy + rotation.Rot[2][2]*state.Vz,
		T:  state.T,
	}
}

// IdentityMatrix returns a rotation matrix that has no effect on orientation.
// It can be the starting point for building custom rotations with Pivot.
func IdentityMatrix() RotationMatrix {
	var r RotationMatrix
	r.Rot[0][0] = 1.0
	r.Rot[1][1] = 1.0
	r.Rot[2][2] = 1.0
	return r
}

// InverseRotation returns a rotation matrix that undoes another rotation.
// Multiplying a vector by the inverse matrix is equivalent to rotating it
// back through the opposite orientation.
func InverseRotation(rotation RotationMatrix) RotationMatrix {
	var inverse RotationMatrix
	inverse.Rot[0][0] = rotation.Rot[0][0]
	inverse.Rot[0][1] = rotation.Rot[1][0]
	inverse.Rot[0][2] = rotation.Rot[2][0]
	inverse.Rot[1][0] = rotation.Rot[0][1]
	inverse.Rot[1][1] = rotation.Rot[1][1]
	inverse.Rot[1][2] = rotation.Rot[2][1]
	inverse.Rot[2][0] = rotation.Rot[0][2]
	inverse.Rot[2][1] = rotation.Rot[1][2]
	inverse.Rot[2][2] = rotation.Rot[2][2]
	return inverse
}

// CombineRotation creates a rotation matrix that is equivalent to applying
// rotation a first, followed by rotation b.
func CombineRotation(a RotationMatrix, b RotationMatrix) RotationMatrix {
	var c RotationMatrix
	c.Rot[0][0] = b.Rot[0][0]*a.Rot[0][0] + b.Rot[1][0]*a.Rot[0][1] + b.Rot[2][0]*a.Rot[0][2]
	c.Rot[1][0] = b.Rot[0][0]*a.Rot[1][0] + b.Rot[1][0]*a.Rot[1][1] + b.Rot[2][0]*a.Rot[1][2]
	c.Rot[2][0] = b.Rot[0][0]*a.Rot[2][0] + b.Rot[1][0]*a.Rot[2][1] + b.Rot[2][0]*a.Rot[2][2]
	c.Rot[0][1] = b.Rot[0][1]*a.Rot[0][0] + b.Rot[1][1]*a.Rot[0][1] + b.Rot[2][1]*a.Rot[0][2]
	c.Rot[1][1] = b.Rot[0][1]*a.Rot[1][0] + b.Rot[1][1]*a.Rot[1][1] + b.Rot[2][1]*a.Rot[1][2]
	c.Rot[2][1] = b.Rot[0][1]*a.Rot[2][0] + b.Rot[1][1]*a.Rot[2][1] + b.Rot[2][1]*a.Rot[2][2]
	c.Rot[0][2] = b.Rot[0][2]*a.Rot[0][0] + b.Rot[1][2]*a.Rot[0][1] + b.Rot[2][2]*a.Rot[0][2]
	c.Rot[1][2] = b.Rot[0][2]*a.Rot[1][0] + b.Rot[1][2]*a.Rot[1][1] + b.Rot[2][2]*a.Rot[1][2]
	c.Rot[2][2] = b.Rot[0][2]*a.Rot[2][0] + b.Rot[1][2]*a.Rot[2][1] + b.Rot[2][2]*a.Rot[2][2]
	return c
}

// Pivot re-orients a rotation matrix by rotating it by the given angle in
// degrees around one of its coordinate axes (0=x, 1=y, 2=z).
// Positive angles are counterclockwise as seen from the positive axis direction.
func Pivot(rotation RotationMatrix, axis int, angle float64) (RotationMatrix, error) {
	if axis < 0 || axis > 2 {
		return RotationMatrix{}, errors.New("Invalid coordinate axis")
	}
	if !isfinite(angle) {
		return RotationMatrix{}, errors.New("Invalid angle")
	}
	radians := RadiansFromDegrees(angle)
	c := math.Cos(radians)
	s := math.Sin(radians)

	// Maintain the right-hand rule: pick the (i,j,k) axis order such that i x j = k.
	i := (axis + 1) % 3
	j := (axis + 2) % 3
	k := axis

	var p RotationMatrix
	p.Rot[i][i] = c*rotation.Rot[i][i] - s*rotation.Rot[i][j]
	p.Rot[i][j] = s*rotation.Rot[i][i] + c*rotation.Rot[i][j]
	p.Rot[i][k] = rotation.Rot[i][k]

	p.Rot[j][i] = c*rotation.Rot[j][i] - s*rotation.Rot[j][j]
	p.Rot[j][j] = s*rotation.Rot[j][i] + c*rotation.Rot[j][j]
	p.Rot[j][k] = rotation.Rot[j][k]

	p.Rot[k][i] = c*rotation.Rot[k][i] - s*rotation.Rot[k][j]
	p.Rot[k][j] = s*rotation.Rot[k][i] + c*rotation.Rot[k][j]
	p.Rot[k][k] = rotation.Rot[k][k]

	return p, nil
}

// VectorFromSphere converts spherical coordinates (latitude, longitude,
// distance in AU) into a Cartesian vector at the given time.
func VectorFromSphere(sphere Spherical, time AstroTime) AstroVector {
	radlat := RadiansFromDegrees(sphere.Lat)
	radlon := RadiansFromDegrees(sphere.Lon)
	rcoslat := sphere.Dist * math.Cos(radlat)
	return AstroVector{
		X: rcoslat * math.Cos(radlon),
		Y: rcoslat * math.Sin(radlon),
		Z: sphere.Dist * math.Sin(radlat),
		T: time,
	}
}

// SphereFromVector converts a Cartesian vector into spherical coordinates.
// Returns an error if the vector has zero length, which would make the
// latitude and longitude indeterminate.
func SphereFromVector(vector AstroVector) (Spherical, error) {
	xyproj := vector.X*vector.X + vector.Y*vector.Y
	dist := math.Sqrt(xyproj + vector.Z*vector.Z)
	if xyproj == 0.0 {
		if vector.Z == 0.0 {
			return Spherical{}, errors.New("Zero-length vector has no spherical coordinates")
		}
		lat := -90.0
		if vector.Z > 0.0 {
			lat = +90.0
		}
		return Spherical{Lat: lat, Lon: 0.0, Dist: dist}, nil
	}
	lon := DegreesFromRadians(math.Atan2(vector.Y, vector.X))
	if lon < 0.0 {
		lon += 360.0
	}
	lat := DegreesFromRadians(math.Atan2(vector.Z, math.Sqrt(xyproj)))
	return Spherical{Lat: lat, Lon: lon, Dist: dist}, nil
}

// EquatorFromVector converts an equatorial Cartesian vector into right
// ascension (in sidereal hours), declination, and distance.
func EquatorFromVector(vector AstroVector) (Equatorial, error) {
	sphere, err := SphereFromVector(vector)
	if err != nil {
		return Equatorial{}, err
	}
	return Equatorial{
		Ra:   sphere.Lon / 15.0, // convert degrees to sidereal hours
		Dec:  sphere.Lat,
		Dist: sphere.Dist,
		Vec:  vector,
	}, nil
}

// cubeRoot returns the real cube root of x.
func cubeRoot(x float64) float64 {
	return math.Cbrt(x)
}

// quadInterp fits a parabola through three equally spaced samples and finds
// the root (t) nearest the middle sample, along with the derivative df/dt.
// Returns ok=false when a unique in-bounds root cannot be found.
func quadInterp(tm, dt, fa, fm, fb float64) (t float64, dfDt float64, ok bool) {
	Q := (fb+fa)/2.0 - fm
	R := (fb - fa) / 2.0
	S := fm

	var x float64
	if Q == 0.0 {
		// This is a line, not a parabola.
		if R == 0.0 {
			return 0.0, 0.0, false // a horizontal line; can't make progress
		}
		x = -S / R
		if x < -1.0 || x > +1.0 {
			return 0.0, 0.0, false // out of bounds
		}
	} else {
		// This really is a parabola. Find roots x1, x2.
		u := R*R - 4*Q*S
		if u <= 0.0 {
			return 0.0, 0.0, false // no real roots, or vertex tangent to axis
		}
		ru := math.Sqrt(u)
		x1 := (-R + ru) / (2.0 * Q)
		x2 := (-R - ru) / (2.0 * Q)
		if -1.0 <= x1 && x1 <= +1.0 {
			if -1.0 <= x2 && x2 <= +1.0 {
				return 0.0, 0.0, false // two roots in bounds; not a unique crossing
			}
			x = x1
		} else if -1.0 <= x2 && x2 <= +1.0 {
			x = x2
		} else {
			return 0.0, 0.0, false // neither root is within bounds
		}
	}

	return tm + x*dt, (2*Q*x + R) / dt, true
}

// longitudeOffset returns the value diff reduced into the half-open range
// (-180, +180] degrees.
func longitudeOffset(diff float64) float64 {
	for diff <= -180.0 {
		diff += 360.0
	}
	for diff > 180.0 {
		diff -= 360.0
	}
	return diff
}

// normalizeLongitude returns the value lon reduced into the half-open range
// [0, 360) degrees.
func normalizeLongitude(lon float64) float64 {
	for lon < 0.0 {
		lon += 360.0
	}
	for lon >= 360.0 {
		lon -= 360.0
	}
	return lon
}

// toggleAzimuthDirection converts an azimuth measured counterclockwise from
// north into one measured clockwise from north, and vice versa.
func toggleAzimuthDirection(az float64) float64 {
	az = 360.0 - az
	if az >= 360.0 {
		az -= 360.0
	} else if az < 0.0 {
		az += 360.0
	}
	return az
}

// rotateEquatorialToEcliptic rotates a Cartesian position from the equatorial
// frame to the ecliptic frame, using the given obliquity angle in radians.
func rotateEquatorialToEcliptic(pos [3]float64, obliqRadians float64, time AstroTime) Ecliptic {
	cosOb := math.Cos(obliqRadians)
	sinOb := math.Sin(obliqRadians)

	vec := AstroVector{
		X: pos[0],
		Y: pos[1]*cosOb + pos[2]*sinOb,
		Z: -pos[1]*sinOb + pos[2]*cosOb,
		T: time,
	}

	xyproj := math.Hypot(vec.X, vec.Y)
	elon := 0.0
	if xyproj > 0.0 {
		elon = DegreesFromRadians(math.Atan2(vec.Y, vec.X))
		if elon < 0.0 {
			elon += 360.0
		}
	}
	return Ecliptic{
		Vec:  vec,
		Elat: DegreesFromRadians(math.Atan2(vec.Z, xyproj)),
		Elon: elon,
	}
}

// vsopTerm is a single trigonometric term in a VSOP87 series.
type vsopTerm struct {
	amp float64 // amplitude
	pha float64 // phase
	fre float64 // frequency
}

// vsopSeries is one power-of-t series of vsopTerm values.
type vsopSeries struct {
	nterms int
	term   []vsopTerm
}

// vsopFormula holds all the power-of-t series for one spherical coordinate.
type vsopFormula struct {
	nseries int
	series  []vsopSeries
}

// vsopModel holds the longitude, latitude, and radius formulas for one body.
type vsopModel struct {
	formula [3]vsopFormula
}

const (
	vsopLonIndex = 0 // index of the longitude formula inside vsopModel
	vsopLatIndex = 1 // index of the latitude formula inside vsopModel
	vsopRadIndex = 2 // index of the radius formula inside vsopModel
)

const daysPerMillennium = 365250.0

// vsopFormulaCalc evaluates the VSOP87B trigonometric series for one body and
// returns ecliptic spherical coordinates: [longitude, latitude, radius] in
// radians and astronomical units. The time t is expressed in thousands of
// Julian years (millennia) since J2000.
func vsopFormulaCalc(model *vsopModel, t float64) [3]float64 {
	var sphere [3]float64
	for k := 0; k < 3; k++ {
		tpower := 1.0
		formula := &model.formula[k]
		sphere[k] = 0.0
		for s := 0; s < formula.nseries; s++ {
			sum := 0.0
			series := &formula.series[s]
			for i := 0; i < series.nterms; i++ {
				term := &series.term[i]
				sum += term.amp * math.Cos(term.pha+(t*term.fre))
			}
			incr := tpower * sum
			if k == vsopLonIndex {
				incr = math.Mod(incr, 2.0*math.Pi) // improve precision for longitudes, which can be hundreds of radians
			}
			sphere[k] += incr
			tpower *= t
		}
	}
	return sphere
}

// vsopDerivCalc evaluates the time derivatives of the VSOP87B series and
// returns d(lon)/dt, d(lat)/dt, dr/dt in units per millennium.
func vsopDerivCalc(model *vsopModel, t float64) [3]float64 {
	var deriv [3]float64
	for k := 0; k < 3; k++ {
		tpower := 1.0
		dpower := 0.0
		formula := &model.formula[k]
		deriv[k] = 0.0
		for s := 0; s < formula.nseries; s++ {
			sinSum := 0.0
			cosSum := 0.0
			series := &formula.series[s]
			for i := 0; i < series.nterms; i++ {
				term := &series.term[i]
				angle := term.pha + (t * term.fre)
				sinSum += term.amp * term.fre * math.Sin(angle)
				if s > 0 {
					cosSum += term.amp * math.Cos(angle)
				}
			}
			deriv[k] += (float64(s) * dpower * cosSum) - (tpower * sinSum)
			dpower = tpower
			tpower *= t
		}
	}
	return deriv
}

// vsopRotate converts a VSOP87A ecliptic Cartesian vector into an equatorial
// Cartesian vector, applying the FK5 correction to the J2000 mean equator.
func vsopRotate(ecl [3]float64) terseVector {
	//     X        +1.000000000000  +0.000000440360  -0.000000190919   X
	//     Y     =  -0.000000479966  +0.917482137087  -0.397776982902   Y
	//     Z FK5     0.000000000000  +0.397776982902  +0.917482137087   Z VSOP87A
	return terseVector{
		X: ecl[0] + 0.000000440360*ecl[1] - 0.000000190919*ecl[2],
		Y: -0.000000479966*ecl[0] + 0.917482137087*ecl[1] - 0.397776982902*ecl[2],
		Z: 0.397776982902*ecl[1] + 0.917482137087*ecl[2],
	}
}

// vsopSphereToRect converts ecliptic spherical coordinates (longitude and
// latitude in radians, radius in AU) into ecliptic Cartesian coordinates.
func vsopSphereToRect(lon float64, lat float64, radius float64) [3]float64 {
	rCosLat := radius * math.Cos(lat)
	return [3]float64{
		rCosLat * math.Cos(lon),
		rCosLat * math.Sin(lon),
		radius * math.Sin(lat),
	}
}

type addSolTerm struct {
	coeffl float64
	coeffs float64
	coeffg float64
	coeffp float64
	p      int
	q      int
	r      int
	s      int
}

type pascalArray struct {
	array [13][4]float64
}

func (p *pascalArray) access(x, y int) *float64 {
	return &p.array[x+6][y-1]
}

type moonContext struct {
	T                        float64
	DGAM                     float64
	DLAM, N, GAM1C, SINPI    float64
	L0, L, LS, F, D, S       float64
	DL0, DL, DLS, DF, DD, DS float64
	CO                       pascalArray
	SI                       pascalArray
}

func addThe(c1, s1, c2, s2 float64, c, s *float64) {
	*c = c1*c2 - s1*s2
	*s = s1*c2 + c1*s2
}

func sineRev(phi float64) float64 {
	// sine of phi in revolutions, not radians
	return math.Sin(2.0 * math.Pi * phi)
}

func (mc *moonContext) longPeriodic() {
	S1 := sineRev(0.19833 + 0.05611*mc.T)
	S2 := sineRev(0.27869 + 0.04508*mc.T)
	S3 := sineRev(0.16827 - 0.36903*mc.T)
	S4 := sineRev(0.34734 - 5.37261*mc.T)
	S5 := sineRev(0.10498 - 5.37899*mc.T)
	S6 := sineRev(0.42681 - 0.41855*mc.T)
	S7 := sineRev(0.14943 - 5.37511*mc.T)

	mc.DL0 = 0.84*S1 + 0.31*S2 + 14.27*S3 + 7.26*S4 + 0.28*S5 + 0.24*S6
	mc.DL = 2.94*S1 + 0.31*S2 + 14.27*S3 + 9.34*S4 + 1.12*S5 + 0.83*S6
	mc.DLS = -6.40*S1 - 1.89*S6
	mc.DF = 0.21*S1 + 0.31*S2 + 14.27*S3 - 88.70*S4 - 15.30*S5 + 0.24*S6 - 1.86*S7
	mc.DD = mc.DL0 - mc.DLS
	mc.DGAM = -3332e-9*sineRev(0.59734-5.37261*mc.T) - 539e-9*sineRev(0.35498-5.37899*mc.T) - 64e-9*sineRev(0.39943-5.37511*mc.T)
}

func (mc *moonContext) term(p, q, r, s int, x, y *float64) {
	*x = 1.0
	*y = 0.0
	if p != 0 {
		addThe(*x, *y, *mc.CO.access(p, 1), *mc.SI.access(p, 1), x, y)
	}
	if q != 0 {
		addThe(*x, *y, *mc.CO.access(q, 2), *mc.SI.access(q, 2), x, y)
	}
	if r != 0 {
		addThe(*x, *y, *mc.CO.access(r, 3), *mc.SI.access(r, 3), x, y)
	}
	if s != 0 {
		addThe(*x, *y, *mc.CO.access(s, 4), *mc.SI.access(s, 4), x, y)
	}
}

func (mc *moonContext) addSol(coeffl, coeffs, coeffg, coeffp float64, p, q, r, s int) {
	var x, y float64
	mc.term(p, q, r, s, &x, &y)
	mc.DLAM += coeffl * y
	mc.DS += coeffs * y
	mc.GAM1C += coeffg * x
	mc.SINPI += coeffp * x
}

func (mc *moonContext) addn(coeffn float64, p, q, r, s int) {
	var x, y float64
	mc.term(p, q, r, s, &x, &y)
	mc.N += coeffn * y
}

func (mc *moonContext) solarn() {
	mc.N = 0.0
	mc.addn(-526.069, 0, 0, 1, -2)
	mc.addn(-3.352, 0, 0, 1, -4)
	mc.addn(+44.297, +1, 0, 1, -2)
	mc.addn(-6.000, +1, 0, 1, -4)
	mc.addn(+20.599, -1, 0, 1, 0)
	mc.addn(-30.598, -1, 0, 1, -2)
	mc.addn(-24.649, -2, 0, 1, 0)
	mc.addn(-2.000, -2, 0, 1, -2)
	mc.addn(-22.571, 0, +1, 1, -2)
	mc.addn(+10.985, 0, -1, 1, -2)
}

func (mc *moonContext) planetary() {
	mc.DLAM += +0.82*sineRev(0.7736-62.5512*mc.T) + 0.31*sineRev(0.0466-125.1025*mc.T)
	mc.DLAM += +0.35*sineRev(0.5785-25.1042*mc.T) + 0.66*sineRev(0.4591+1335.8075*mc.T)
	mc.DLAM += +0.64*sineRev(0.3130-91.5680*mc.T) + 1.14*sineRev(0.1480+1331.2898*mc.T)
	mc.DLAM += +0.21*sineRev(0.5918+1056.5859*mc.T) + 0.44*sineRev(0.5784+1322.8595*mc.T)
	mc.DLAM += +0.24*sineRev(0.2275-5.7374*mc.T) + 0.28*sineRev(0.2965+2.6929*mc.T)
	mc.DLAM += +0.33 * sineRev(0.3132+6.3368*mc.T)
}

func twoPiFrac(x float64) float64 {
	return (2.0 * math.Pi) * (x - math.Floor(x))
}

func makeMoonContext(centuriesSinceJ2000 float64) *moonContext {
	var I, J, MAX int
	var T2, ARG, FAC float64
	var c, s float64

	mc := moonContext{}

	mc.T = centuriesSinceJ2000
	T2 = mc.T * mc.T
	mc.DLAM = 0.0
	mc.DS = 0.0
	mc.GAM1C = 0
	mc.SINPI = 3422.7
	mc.longPeriodic()
	mc.L0 = twoPiFrac(0.60643382+1336.85522467*mc.T-0.00000313*T2) + mc.DL0/arc
	mc.L = twoPiFrac(0.37489701+1325.55240982*mc.T+0.00002565*T2) + mc.DL/arc
	mc.LS = twoPiFrac(0.99312619+99.99735956*mc.T-0.00000044*T2) + mc.DLS/arc
	mc.F = twoPiFrac(0.25909118+1342.22782980*mc.T-0.00000892*T2) + mc.DF/arc
	mc.D = twoPiFrac(0.82736186+1236.85308708*mc.T-0.00000397*T2) + mc.DD/arc

	for I = 1; I <= 4; I += 1 {
		switch I {
		case 1:
			{
				ARG = mc.L
				MAX = 4
				FAC = 1.000002208
			}
		case 2:
			{
				ARG = mc.LS
				MAX = 3
				FAC = 0.997504612 - 0.002495388*mc.T
			}
		case 3:
			{
				ARG = mc.F
				MAX = 4
				FAC = 1.000002708 + 139.978*mc.DGAM
			}
		default:
			{
				ARG = mc.D
				MAX = 6
				FAC = 1.0
			}
		}
		*mc.CO.access(0, I) = 1.0
		*mc.CO.access(1, I) = math.Cos(ARG) * FAC
		*mc.SI.access(0, I) = 0.0
		*mc.SI.access(1, I) = math.Sin(ARG) * FAC

		for J = 2; J <= MAX; J += 1 {
			addThe(*mc.CO.access(J-1, I), *mc.SI.access(J-1, I), *mc.CO.access(1, I), *mc.SI.access(1, I), &c, &s)
			*mc.CO.access(J, I) = c
			*mc.SI.access(J, I) = s
		}

		for J = 1; J <= MAX; J += 1 {
			*mc.CO.access(-J, I) = +*mc.CO.access(J, I)
			*mc.SI.access(-J, I) = -*mc.SI.access(J, I)
		}
	}

	return &mc
}

type moonResult struct {
	geoEclipLon float64
	geoEclipLat float64
	distanceAu  float64
}

func (mc *moonContext) calcMoon() moonResult {
	for _, t := range moonAddSolTerms {
		mc.addSol(t.coeffl, t.coeffs, t.coeffg, t.coeffp, t.p, t.q, t.r, t.s)
	}
	mc.solarn()
	mc.planetary()
	mc.S = mc.F + mc.DS/arc
	latSeconds := (1.000002708+139.978*mc.DGAM)*(18518.511+1.189+mc.GAM1C)*math.Sin(mc.S) - 6.24*math.Sin(3*mc.S) + mc.N
	return moonResult{
		twoPiFrac((mc.L0 + mc.DLAM/arc) / (2.0 * math.Pi)),
		latSeconds / arc,
		(arc * EarthEquatorialRadiusAu) / (0.999953253 * mc.SINPI),
	}
}

// GeoMoon calculates the equatorial geocentric position of the Moon at a given time.
// The returned vector indicates the Moon's center relative to the Earth's center.
// The vector components are expressed in AU (astronomical units).
// The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
// In Astronomy Engine, this orientation is called EQJ.
func GeoMoon(time AstroTime) AstroVector {
	mc := makeMoonContext(time.Tt / 36525.0)
	moon := mc.calcMoon()

	// Convert geocentric ecliptic spherical coordinates to Cartesian coordinates.
	distCosLat := moon.distanceAu * math.Cos(moon.geoEclipLat)

	gepos := AstroVector{
		distCosLat * math.Cos(moon.geoEclipLon),
		distCosLat * math.Sin(moon.geoEclipLon),
		moon.distanceAu * math.Sin(moon.geoEclipLat),
		time,
	}

	// Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
	mpos1 := eclToEquVec(gepos)

	// Convert from mean equinox of date to J2000.
	mpos2 := precession(mpos1, into2000)

	return mpos2
}

//--- Generated code begins here ------------------------------------------------------------------

var vsopLonMercury_0 = []vsopTerm{
	{4.40250710144, 0.00000000000, 0.00000000000},      // 0
	{0.40989414977, 1.48302034195, 26087.90314157420},  // 1
	{0.05046294200, 4.47785489551, 52175.80628314840},  // 2
	{0.00855346844, 1.16520322459, 78263.70942472259},  // 3
	{0.00165590362, 4.11969163423, 104351.61256629678}, // 4
	{0.00034561897, 0.77930768443, 130439.51570787099}, // 5
	{0.00007583476, 3.71348404924, 156527.41884944518}, // 6
}

var vsopLonMercury_1 = []vsopTerm{
	{26087.90313685529, 0.00000000000, 0.00000000000},  // 0
	{0.01131199811, 6.21874197797, 26087.90314157420},  // 1
	{0.00292242298, 3.04449355541, 52175.80628314840},  // 2
	{0.00075775081, 6.08568821653, 78263.70942472259},  // 3
	{0.00019676525, 2.80965111777, 104351.61256629678}, // 4
}

var vsopLonMercury = []vsopSeries{
	{7, vsopLonMercury_0},
	{5, vsopLonMercury_1},
}

var vsopLatMercury_0 = []vsopTerm{
	{0.11737528961, 1.98357498767, 26087.90314157420},  // 0
	{0.02388076996, 5.03738959686, 52175.80628314840},  // 1
	{0.01222839532, 3.14159265359, 0.00000000000},      // 2
	{0.00543251810, 1.79644363964, 78263.70942472259},  // 3
	{0.00129778770, 4.83232503958, 104351.61256629678}, // 4
	{0.00031866927, 1.58088495658, 130439.51570787099}, // 5
	{0.00007963301, 4.60972126127, 156527.41884944518}, // 6
}

var vsopLatMercury_1 = []vsopTerm{
	{0.00274646065, 3.95008450011, 26087.90314157420}, // 0
	{0.00099737713, 3.14159265359, 0.00000000000},     // 1
}

var vsopLatMercury = []vsopSeries{
	{7, vsopLatMercury_0},
	{2, vsopLatMercury_1},
}

var vsopRadMercury_0 = []vsopTerm{
	{0.39528271651, 0.00000000000, 0.00000000000},      // 0
	{0.07834131818, 6.19233722598, 26087.90314157420},  // 1
	{0.00795525558, 2.95989690104, 52175.80628314840},  // 2
	{0.00121281764, 6.01064153797, 78263.70942472259},  // 3
	{0.00021921969, 2.77820093972, 104351.61256629678}, // 4
	{0.00004354065, 5.82894543774, 130439.51570787099}, // 5
}

var vsopRadMercury_1 = []vsopTerm{
	{0.00217347740, 4.65617158665, 26087.90314157420}, // 0
	{0.00044141826, 1.42385544001, 52175.80628314840}, // 1
}

var vsopRadMercury = []vsopSeries{
	{6, vsopRadMercury_0},
	{2, vsopRadMercury_1},
}

var vsopLonVenus_0 = []vsopTerm{
	{3.17614666774, 0.00000000000, 0.00000000000},     // 0
	{0.01353968419, 5.59313319619, 10213.28554621100}, // 1
	{0.00089891645, 5.30650047764, 20426.57109242200}, // 2
	{0.00005477194, 4.41630661466, 7860.41939243920},  // 3
	{0.00003455741, 2.69964447820, 11790.62908865880}, // 4
	{0.00002372061, 2.99377542079, 3930.20969621960},  // 5
	{0.00001317168, 5.18668228402, 26.29831979980},    // 6
	{0.00001664146, 4.25018630147, 1577.34354244780},  // 7
	{0.00001438387, 4.15745084182, 9683.59458111640},  // 8
	{0.00001200521, 6.15357116043, 30639.85663863300}, // 9
}

var vsopLonVenus_1 = []vsopTerm{
	{10213.28554621638, 0.00000000000, 0.00000000000}, // 0
	{0.00095617813, 2.46406511110, 10213.28554621100}, // 1
	{0.00007787201, 0.62478482220, 20426.57109242200}, // 2
}

var vsopLonVenus = []vsopSeries{
	{10, vsopLonVenus_0},
	{3, vsopLonVenus_1},
}

var vsopLatVenus_0 = []vsopTerm{
	{0.05923638472, 0.26702775812, 10213.28554621100}, // 0
	{0.00040107978, 1.14737178112, 20426.57109242200}, // 1
	{0.00032814918, 3.14159265359, 0.00000000000},     // 2
}

var vsopLatVenus_1 = []vsopTerm{
	{0.00287821243, 1.88964962838, 10213.28554621100}, // 0
}

var vsopLatVenus = []vsopSeries{
	{3, vsopLatVenus_0},
	{1, vsopLatVenus_1},
}

var vsopRadVenus_0 = []vsopTerm{
	{0.72334820891, 0.00000000000, 0.00000000000},     // 0
	{0.00489824182, 4.02151831717, 10213.28554621100}, // 1
	{0.00001658058, 4.90206728031, 20426.57109242200}, // 2
	{0.00001378043, 1.12846591367, 11790.62908865880}, // 3
	{0.00001632096, 2.84548795207, 7860.41939243920},  // 4
	{0.00000498395, 2.58682193892, 9683.59458111640},  // 5
	{0.00000221985, 2.01346696541, 19367.18916223280}, // 6
	{0.00000237454, 2.55136053886, 15720.83878487840}, // 7
}

var vsopRadVenus_1 = []vsopTerm{
	{0.00034551041, 0.89198706276, 10213.28554621100}, // 0
}

var vsopRadVenus = []vsopSeries{
	{8, vsopRadVenus_0},
	{1, vsopRadVenus_1},
}

var vsopLonEarth_0 = []vsopTerm{
	{1.75347045673, 0.00000000000, 0.00000000000},     // 0
	{0.03341656453, 4.66925680415, 6283.07584999140},  // 1
	{0.00034894275, 4.62610242189, 12566.15169998280}, // 2
	{0.00003417572, 2.82886579754, 3.52311834900},     // 3
	{0.00003497056, 2.74411783405, 5753.38488489680},  // 4
	{0.00003135899, 3.62767041756, 77713.77146812050}, // 5
	{0.00002676218, 4.41808345438, 7860.41939243920},  // 6
	{0.00002342691, 6.13516214446, 3930.20969621960},  // 7
	{0.00001273165, 2.03709657878, 529.69096509460},   // 8
	{0.00001324294, 0.74246341673, 11506.76976979360}, // 9
	{0.00000901854, 2.04505446477, 26.29831979980},    // 10
	{0.00001199167, 1.10962946234, 1577.34354244780},  // 11
	{0.00000857223, 3.50849152283, 398.14900340820},   // 12
	{0.00000779786, 1.17882681962, 5223.69391980220},  // 13
	{0.00000990250, 5.23268072088, 5884.92684658320},  // 14
	{0.00000753141, 2.53339052847, 5507.55323866740},  // 15
	{0.00000505267, 4.58292599973, 18849.22754997420}, // 16
	{0.00000492392, 4.20505711826, 775.52261132400},   // 17
	{0.00000356672, 2.91954114478, 0.06731030280},     // 18
	{0.00000284125, 1.89869240932, 796.29800681640},   // 19
	{0.00000242879, 0.34481445893, 5486.77784317500},  // 20
	{0.00000317087, 5.84901948512, 11790.62908865880}, // 21
	{0.00000271112, 0.31486255375, 10977.07880469900}, // 22
	{0.00000206217, 4.80646631478, 2544.31441988340},  // 23
	{0.00000205478, 1.86953770281, 5573.14280143310},  // 24
	{0.00000202318, 2.45767790232, 6069.77675455340},  // 25
	{0.00000126225, 1.08295459501, 20.77539549240},    // 26
	{0.00000155516, 0.83306084617, 213.29909543800},   // 27
}

var vsopLonEarth_1 = []vsopTerm{
	{6283.07584999140, 0.00000000000, 0.00000000000},  // 0
	{0.00206058863, 2.67823455808, 6283.07584999140},  // 1
	{0.00004303419, 2.63512233481, 12566.15169998280}, // 2
}

var vsopLonEarth_2 = []vsopTerm{
	{0.00008721859, 1.07253635559, 6283.07584999140}, // 0
}

var vsopLonEarth = []vsopSeries{
	{28, vsopLonEarth_0},
	{3, vsopLonEarth_1},
	{1, vsopLonEarth_2},
}

var vsopLatEarth_1 = []vsopTerm{
	{0.00227777722, 3.41376620530, 6283.07584999140},  // 0
	{0.00003805678, 3.37063423795, 12566.15169998280}, // 1
}

var vsopLatEarth = []vsopSeries{
	{0, nil},
	{2, vsopLatEarth_1},
}

var vsopRadEarth_0 = []vsopTerm{
	{1.00013988784, 0.00000000000, 0.00000000000},      // 0
	{0.01670699632, 3.09846350258, 6283.07584999140},   // 1
	{0.00013956024, 3.05524609456, 12566.15169998280},  // 2
	{0.00003083720, 5.19846674381, 77713.77146812050},  // 3
	{0.00001628463, 1.17387558054, 5753.38488489680},   // 4
	{0.00001575572, 2.84685214877, 7860.41939243920},   // 5
	{0.00000924799, 5.45292236722, 11506.76976979360},  // 6
	{0.00000542439, 4.56409151453, 3930.20969621960},   // 7
	{0.00000472110, 3.66100022149, 5884.92684658320},   // 8
	{0.00000085831, 1.27079125277, 161000.68573767410}, // 9
	{0.00000057056, 2.01374292245, 83996.84731811189},  // 10
	{0.00000055736, 5.24159799170, 71430.69561812909},  // 11
	{0.00000174844, 3.01193636733, 18849.22754997420},  // 12
	{0.00000243181, 4.27349530790, 11790.62908865880},  // 13
}

var vsopRadEarth_1 = []vsopTerm{
	{0.00103018607, 1.10748968172, 6283.07584999140},  // 0
	{0.00001721238, 1.06442300386, 12566.15169998280}, // 1
}

var vsopRadEarth_2 = []vsopTerm{
	{0.00004359385, 5.78455133808, 6283.07584999140}, // 0
}

var vsopRadEarth = []vsopSeries{
	{14, vsopRadEarth_0},
	{2, vsopRadEarth_1},
	{1, vsopRadEarth_2},
}

var vsopLonMars_0 = []vsopTerm{
	{6.20347711581, 0.00000000000, 0.00000000000},     // 0
	{0.18656368093, 5.05037100270, 3340.61242669980},  // 1
	{0.01108216816, 5.40099836344, 6681.22485339960},  // 2
	{0.00091798406, 5.75478744667, 10021.83728009940}, // 3
	{0.00027744987, 5.97049513147, 3.52311834900},     // 4
	{0.00010610235, 2.93958560338, 2281.23049651060},  // 5
	{0.00012315897, 0.84956094002, 2810.92146160520},  // 6
	{0.00008926784, 4.15697846427, 0.01725365220},     // 7
	{0.00008715691, 6.11005153139, 13362.44970679920}, // 8
	{0.00006797556, 0.36462229657, 398.14900340820},   // 9
	{0.00007774872, 3.33968761376, 5621.84292321040},  // 10
	{0.00003575078, 1.66186505710, 2544.31441988340},  // 11
	{0.00004161108, 0.22814971327, 2942.46342329160},  // 12
	{0.00003075252, 0.85696614132, 191.44826611160},   // 13
	{0.00002628117, 0.64806124465, 3337.08930835080},  // 14
	{0.00002937546, 6.07893711402, 0.06731030280},     // 15
	{0.00002389414, 5.03896442664, 796.29800681640},   // 16
	{0.00002579844, 0.02996736156, 3344.13554504880},  // 17
	{0.00001528141, 1.14979301996, 6151.53388830500},  // 18
	{0.00001798806, 0.65634057445, 529.69096509460},   // 19
	{0.00001264357, 3.62275122593, 5092.15195811580},  // 20
	{0.00001286228, 3.06796065034, 2146.16541647520},  // 21
	{0.00001546404, 2.91579701718, 1751.53953141600},  // 22
	{0.00001024902, 3.69334099279, 8962.45534991020},  // 23
	{0.00000891566, 0.18293837498, 16703.06213349900}, // 24
	{0.00000858759, 2.40093811940, 2914.01423582380},  // 25
	{0.00000832715, 2.46418619474, 3340.59517304760},  // 26
	{0.00000832720, 4.49495782139, 3340.62968035200},  // 27
	{0.00000712902, 3.66335473479, 1059.38193018920},  // 28
	{0.00000748723, 3.82248614017, 155.42039943420},   // 29
	{0.00000723861, 0.67497311481, 3738.76143010800},  // 30
	{0.00000635548, 2.92182225127, 8432.76438481560},  // 31
	{0.00000655162, 0.48864064125, 3127.31333126180},  // 32
	{0.00000550474, 3.81001042328, 0.98032106820},     // 33
	{0.00000552750, 4.47479317037, 1748.01641306700},  // 34
	{0.00000425966, 0.55364317304, 6283.07584999140},  // 35
	{0.00000415131, 0.49662285038, 213.29909543800},   // 36
	{0.00000472167, 3.62547124025, 1194.44701022460},  // 37
	{0.00000306551, 0.38052848348, 6684.74797174860},  // 38
	{0.00000312141, 0.99853944405, 6677.70173505060},  // 39
	{0.00000293198, 4.22131299634, 20.77539549240},    // 40
	{0.00000302375, 4.48618007156, 3532.06069281140},  // 41
	{0.00000274027, 0.54222167059, 3340.54511639700},  // 42
	{0.00000281079, 5.88163521788, 1349.86740965880},  // 43
	{0.00000231183, 1.28242156993, 3870.30339179440},  // 44
	{0.00000283602, 5.76885434940, 3149.16416058820},  // 45
	{0.00000236117, 5.75503217933, 3333.49887969900},  // 46
	{0.00000274033, 0.13372524985, 3340.67973700260},  // 47
	{0.00000299395, 2.78323740866, 6254.62666252360},  // 48
}

var vsopLonMars_1 = []vsopTerm{
	{3340.61242700512, 0.00000000000, 0.00000000000},  // 0
	{0.01457554523, 3.60433733236, 3340.61242669980},  // 1
	{0.00168414711, 3.92318567804, 6681.22485339960},  // 2
	{0.00020622975, 4.26108844583, 10021.83728009940}, // 3
	{0.00003452392, 4.73210393190, 3.52311834900},     // 4
	{0.00002586332, 4.60670058555, 13362.44970679920}, // 5
	{0.00000841535, 4.45864030426, 2281.23049651060},  // 6
}

var vsopLonMars_2 = []vsopTerm{
	{0.00058152577, 2.04961712429, 3340.61242669980}, // 0
	{0.00013459579, 2.45738706163, 6681.22485339960}, // 1
}

var vsopLonMars = []vsopSeries{
	{49, vsopLonMars_0},
	{7, vsopLonMars_1},
	{2, vsopLonMars_2},
}

var vsopLatMars_0 = []vsopTerm{
	{0.03197134986, 3.76832042431, 3340.61242669980},  // 0
	{0.00298033234, 4.10616996305, 6681.22485339960},  // 1
	{0.00289104742, 0.00000000000, 0.00000000000},     // 2
	{0.00031365539, 4.44651053090, 10021.83728009940}, // 3
	{0.00003484100, 4.78812549260, 13362.44970679920}, // 4
}

var vsopLatMars_1 = []vsopTerm{
	{0.00217310991, 6.04472194776, 3340.61242669980}, // 0
	{0.00020976948, 3.14159265359, 0.00000000000},    // 1
	{0.00012834709, 1.60810667915, 6681.22485339960}, // 2
}

var vsopLatMars = []vsopSeries{
	{5, vsopLatMars_0},
	{3, vsopLatMars_1},
}

var vsopRadMars_0 = []vsopTerm{
	{1.53033488271, 0.00000000000, 0.00000000000},     // 0
	{0.14184953160, 3.47971283528, 3340.61242669980},  // 1
	{0.00660776362, 3.81783443019, 6681.22485339960},  // 2
	{0.00046179117, 4.15595316782, 10021.83728009940}, // 3
	{0.00008109733, 5.55958416318, 2810.92146160520},  // 4
	{0.00007485318, 1.77239078402, 5621.84292321040},  // 5
	{0.00005523191, 1.36436303770, 2281.23049651060},  // 6
	{0.00003825160, 4.49407183687, 13362.44970679920}, // 7
	{0.00002306537, 0.09081579001, 2544.31441988340},  // 8
	{0.00001999396, 5.36059617709, 3337.08930835080},  // 9
	{0.00002484394, 4.92545639920, 2942.46342329160},  // 10
	{0.00001960195, 4.74249437639, 3344.13554504880},  // 11
	{0.00001167119, 2.11260868341, 5092.15195811580},  // 12
	{0.00001102816, 5.00908403998, 398.14900340820},   // 13
	{0.00000899066, 4.40791133207, 529.69096509460},   // 14
	{0.00000992252, 5.83861961952, 6151.53388830500},  // 15
	{0.00000807354, 2.10217065501, 1059.38193018920},  // 16
	{0.00000797915, 3.44839203899, 796.29800681640},   // 17
	{0.00000740975, 1.49906336885, 2146.16541647520},  // 18
}

var vsopRadMars_1 = []vsopTerm{
	{0.01107433345, 2.03250524857, 3340.61242669980},  // 0
	{0.00103175887, 2.37071847807, 6681.22485339960},  // 1
	{0.00012877200, 0.00000000000, 0.00000000000},     // 2
	{0.00010815880, 2.70888095665, 10021.83728009940}, // 3
}

var vsopRadMars_2 = []vsopTerm{
	{0.00044242249, 0.47930604954, 3340.61242669980}, // 0
	{0.00008138042, 0.86998389204, 6681.22485339960}, // 1
}

var vsopRadMars = []vsopSeries{
	{19, vsopRadMars_0},
	{4, vsopRadMars_1},
	{2, vsopRadMars_2},
}

var vsopLonJupiter_0 = []vsopTerm{
	{0.59954691494, 0.00000000000, 0.00000000000},    // 0
	{0.09695898719, 5.06191793158, 529.69096509460},  // 1
	{0.00573610142, 1.44406205629, 7.11354700080},    // 2
	{0.00306389205, 5.41734730184, 1059.38193018920}, // 3
	{0.00097178296, 4.14264726552, 632.78373931320},  // 4
	{0.00072903078, 3.64042916389, 522.57741809380},  // 5
	{0.00064263975, 3.41145165351, 103.09277421860},  // 6
	{0.00039806064, 2.29376740788, 419.48464387520},  // 7
	{0.00038857767, 1.27231755835, 316.39186965660},  // 8
	{0.00027964629, 1.78454591820, 536.80451209540},  // 9
	{0.00013589730, 5.77481040790, 1589.07289528380}, // 10
	{0.00008246349, 3.58227925840, 206.18554843720},  // 11
	{0.00008768704, 3.63000308199, 949.17560896980},  // 12
	{0.00007368042, 5.08101194270, 735.87651353180},  // 13
	{0.00006263150, 0.02497628807, 213.29909543800},  // 14
	{0.00006114062, 4.51319998626, 1162.47470440780}, // 15
	{0.00004905396, 1.32084470588, 110.20632121940},  // 16
	{0.00005305285, 1.30671216791, 14.22709400160},   // 17
	{0.00005305441, 4.18625634012, 1052.26838318840}, // 18
	{0.00004647248, 4.69958103684, 3.93215326310},    // 19
	{0.00003045023, 4.31676431084, 426.59819087600},  // 20
	{0.00002609999, 1.56667394063, 846.08283475120},  // 21
	{0.00002028191, 1.06376530715, 3.18139373770},    // 22
	{0.00001764763, 2.14148655117, 1066.49547719000}, // 23
	{0.00001722972, 3.88036268267, 1265.56747862640}, // 24
	{0.00001920945, 0.97168196472, 639.89728631400},  // 25
	{0.00001633223, 3.58201833555, 515.46387109300},  // 26
	{0.00001431999, 4.29685556046, 625.67019231240},  // 27
	{0.00000973272, 4.09764549134, 95.97922721780},   // 28
}

var vsopLonJupiter_1 = []vsopTerm{
	{529.69096508814, 0.00000000000, 0.00000000000},  // 0
	{0.00489503243, 4.22082939470, 529.69096509460},  // 1
	{0.00228917222, 6.02646855621, 7.11354700080},    // 2
	{0.00030099479, 4.54540782858, 1059.38193018920}, // 3
	{0.00020720920, 5.45943156902, 522.57741809380},  // 4
	{0.00012103653, 0.16994816098, 536.80451209540},  // 5
	{0.00006067987, 4.42422292017, 103.09277421860},  // 6
	{0.00005433968, 3.98480737746, 419.48464387520},  // 7
	{0.00004237744, 5.89008707199, 14.22709400160},   // 8
}

var vsopLonJupiter_2 = []vsopTerm{
	{0.00047233601, 4.32148536482, 7.11354700080},   // 0
	{0.00030649436, 2.92977788700, 529.69096509460}, // 1
	{0.00014837605, 3.14159265359, 0.00000000000},   // 2
}

var vsopLonJupiter = []vsopSeries{
	{29, vsopLonJupiter_0},
	{9, vsopLonJupiter_1},
	{3, vsopLonJupiter_2},
}

var vsopLatJupiter_0 = []vsopTerm{
	{0.02268615702, 3.55852606721, 529.69096509460},  // 0
	{0.00109971634, 3.90809347197, 1059.38193018920}, // 1
	{0.00110090358, 0.00000000000, 0.00000000000},    // 2
	{0.00008101428, 3.60509572885, 522.57741809380},  // 3
	{0.00006043996, 4.25883108339, 1589.07289528380}, // 4
	{0.00006437782, 0.30627119215, 536.80451209540},  // 5
}

var vsopLatJupiter_1 = []vsopTerm{
	{0.00078203446, 1.52377859742, 529.69096509460}, // 0
}

var vsopLatJupiter = []vsopSeries{
	{6, vsopLatJupiter_0},
	{1, vsopLatJupiter_1},
}

var vsopRadJupiter_0 = []vsopTerm{
	{5.20887429326, 0.00000000000, 0.00000000000},    // 0
	{0.25209327119, 3.49108639871, 529.69096509460},  // 1
	{0.00610599976, 3.84115365948, 1059.38193018920}, // 2
	{0.00282029458, 2.57419881293, 632.78373931320},  // 3
	{0.00187647346, 2.07590383214, 522.57741809380},  // 4
	{0.00086792905, 0.71001145545, 419.48464387520},  // 5
	{0.00072062974, 0.21465724607, 536.80451209540},  // 6
	{0.00065517248, 5.97995884790, 316.39186965660},  // 7
	{0.00029134542, 1.67759379655, 103.09277421860},  // 8
	{0.00030135335, 2.16132003734, 949.17560896980},  // 9
	{0.00023453271, 3.54023522184, 735.87651353180},  // 10
	{0.00022283743, 4.19362594399, 1589.07289528380}, // 11
	{0.00023947298, 0.27458037480, 7.11354700080},    // 12
	{0.00013032614, 2.96042965363, 1162.47470440780}, // 13
	{0.00009703360, 1.90669633585, 206.18554843720},  // 14
	{0.00012749023, 2.71550286592, 1052.26838318840}, // 15
	{0.00007057931, 2.18184839926, 1265.56747862640}, // 16
	{0.00006137703, 6.26418240033, 846.08283475120},  // 17
	{0.00002616976, 2.00994012876, 1581.95934828300}, // 18
}

var vsopRadJupiter_1 = []vsopTerm{
	{0.01271801520, 2.64937512894, 529.69096509460},  // 0
	{0.00061661816, 3.00076460387, 1059.38193018920}, // 1
	{0.00053443713, 3.89717383175, 522.57741809380},  // 2
	{0.00031185171, 4.88276958012, 536.80451209540},  // 3
	{0.00041390269, 0.00000000000, 0.00000000000},    // 4
}

var vsopRadJupiter = []vsopSeries{
	{19, vsopRadJupiter_0},
	{5, vsopRadJupiter_1},
}

var vsopLonSaturn_0 = []vsopTerm{
	{0.87401354025, 0.00000000000, 0.00000000000},   // 0
	{0.11107659762, 3.96205090159, 213.29909543800}, // 1
	{0.01414150957, 4.58581516874, 7.11354700080},   // 2
	{0.00398379389, 0.52112032699, 206.18554843720}, // 3
	{0.00350769243, 3.30329907896, 426.59819087600}, // 4
	{0.00206816305, 0.24658372002, 103.09277421860}, // 5
	{0.00079271300, 3.84007056878, 220.41264243880}, // 6
	{0.00023990355, 4.66976924553, 110.20632121940}, // 7
	{0.00016573588, 0.43719228296, 419.48464387520}, // 8
	{0.00014906995, 5.76903183869, 316.39186965660}, // 9
	{0.00015820290, 0.93809155235, 632.78373931320}, // 10
	{0.00014609559, 1.56518472000, 3.93215326310},   // 11
	{0.00013160301, 4.44891291899, 14.22709400160},  // 12
	{0.00015053543, 2.71669915667, 639.89728631400}, // 13
	{0.00013005299, 5.98119023644, 11.04570026390},  // 14
	{0.00010725067, 3.12939523827, 202.25339517410}, // 15
	{0.00005863206, 0.23656938524, 529.69096509460}, // 16
	{0.00005227757, 4.20783365759, 3.18139373770},   // 17
	{0.00006126317, 1.76328667907, 277.03499374140}, // 18
	{0.00005019687, 3.17787728405, 433.71173787680}, // 19
	{0.00004592550, 0.61977744975, 199.07200143640}, // 20
	{0.00004005867, 2.24479718502, 63.73589830340},  // 21
	{0.00002953796, 0.98280366998, 95.97922721780},  // 22
	{0.00003873670, 3.22283226966, 138.51749687070}, // 23
	{0.00002461186, 2.03163875071, 735.87651353180}, // 24
	{0.00003269484, 0.77492638211, 949.17560896980}, // 25
	{0.00001758145, 3.26580109940, 522.57741809380}, // 26
	{0.00001640172, 5.50504453050, 846.08283475120}, // 27
	{0.00001391327, 4.02333150505, 323.50541665740}, // 28
	{0.00001580648, 4.37265307169, 309.27832265580}, // 29
	{0.00001123498, 2.83726798446, 415.55249061210}, // 30
	{0.00001017275, 3.71700135395, 227.52618943960}, // 31
	{0.00000848642, 3.19150170830, 209.36694217490}, // 32
}

var vsopLonSaturn_1 = []vsopTerm{
	{213.29909521690, 0.00000000000, 0.00000000000}, // 0
	{0.01297370862, 1.82834923978, 213.29909543800}, // 1
	{0.00564345393, 2.88499717272, 7.11354700080},   // 2
	{0.00093734369, 1.06311793502, 426.59819087600}, // 3
	{0.00107674962, 2.27769131009, 206.18554843720}, // 4
	{0.00040244455, 2.04108104671, 220.41264243880}, // 5
	{0.00019941774, 1.27954390470, 103.09277421860}, // 6
	{0.00010511678, 2.74880342130, 14.22709400160},  // 7
	{0.00006416106, 0.38238295041, 639.89728631400}, // 8
	{0.00004848994, 2.43037610229, 419.48464387520}, // 9
	{0.00004056892, 2.92133209468, 110.20632121940}, // 10
	{0.00003768635, 3.64965330780, 3.93215326310},   // 11
}

var vsopLonSaturn_2 = []vsopTerm{
	{0.00116441330, 1.17988132879, 7.11354700080},   // 0
	{0.00091841837, 0.07325195840, 213.29909543800}, // 1
	{0.00036661728, 0.00000000000, 0.00000000000},   // 2
	{0.00015274496, 4.06493179167, 206.18554843720}, // 3
}

var vsopLonSaturn = []vsopSeries{
	{33, vsopLonSaturn_0},
	{12, vsopLonSaturn_1},
	{4, vsopLonSaturn_2},
}

var vsopLatSaturn_0 = []vsopTerm{
	{0.04330678039, 3.60284428399, 213.29909543800}, // 0
	{0.00240348302, 2.85238489373, 426.59819087600}, // 1
	{0.00084745939, 0.00000000000, 0.00000000000},   // 2
	{0.00030863357, 3.48441504555, 220.41264243880}, // 3
	{0.00034116062, 0.57297307557, 206.18554843720}, // 4
	{0.00014734070, 2.11846596715, 639.89728631400}, // 5
	{0.00009916667, 5.79003188904, 419.48464387520}, // 6
	{0.00006993564, 4.73604689720, 7.11354700080},   // 7
	{0.00004807588, 5.43305312061, 316.39186965660}, // 8
}

var vsopLatSaturn_1 = []vsopTerm{
	{0.00198927992, 4.93901017903, 213.29909543800}, // 0
	{0.00036947916, 3.14159265359, 0.00000000000},   // 1
	{0.00017966989, 0.51979431110, 426.59819087600}, // 2
}

var vsopLatSaturn = []vsopSeries{
	{9, vsopLatSaturn_0},
	{3, vsopLatSaturn_1},
}

var vsopRadSaturn_0 = []vsopTerm{
	{9.55758135486, 0.00000000000, 0.00000000000},    // 0
	{0.52921382865, 2.39226219573, 213.29909543800},  // 1
	{0.01873679867, 5.23549604660, 206.18554843720},  // 2
	{0.01464663929, 1.64763042902, 426.59819087600},  // 3
	{0.00821891141, 5.93520042303, 316.39186965660},  // 4
	{0.00547506923, 5.01532618980, 103.09277421860},  // 5
	{0.00371684650, 2.27114821115, 220.41264243880},  // 6
	{0.00361778765, 3.13904301847, 7.11354700080},    // 7
	{0.00140617506, 5.70406606781, 632.78373931320},  // 8
	{0.00108974848, 3.29313390175, 110.20632121940},  // 9
	{0.00069006962, 5.94099540992, 419.48464387520},  // 10
	{0.00061053367, 0.94037691801, 639.89728631400},  // 11
	{0.00048913294, 1.55733638681, 202.25339517410},  // 12
	{0.00034143772, 0.19519102597, 277.03499374140},  // 13
	{0.00032401773, 5.47084567016, 949.17560896980},  // 14
	{0.00020936596, 0.46349251129, 735.87651353180},  // 15
	{0.00009796004, 5.20477537945, 1265.56747862640}, // 16
	{0.00011993338, 5.98050967385, 846.08283475120},  // 17
	{0.00020839300, 1.52102476129, 433.71173787680},  // 18
	{0.00015298404, 3.05943814940, 529.69096509460},  // 19
	{0.00006465823, 0.17732249942, 1052.26838318840}, // 20
	{0.00011380257, 1.73105427040, 522.57741809380},  // 21
	{0.00003419618, 4.94550542171, 1581.95934828300}, // 22
}

var vsopRadSaturn_1 = []vsopTerm{
	{0.06182981340, 0.25843511480, 213.29909543800}, // 0
	{0.00506577242, 0.71114625261, 206.18554843720}, // 1
	{0.00341394029, 5.79635741658, 426.59819087600}, // 2
	{0.00188491195, 0.47215589652, 220.41264243880}, // 3
	{0.00186261486, 3.14159265359, 0.00000000000},   // 4
	{0.00143891146, 1.40744822888, 7.11354700080},   // 5
}

var vsopRadSaturn_2 = []vsopTerm{
	{0.00436902572, 4.78671677509, 213.29909543800}, // 0
}

var vsopRadSaturn = []vsopSeries{
	{23, vsopRadSaturn_0},
	{6, vsopRadSaturn_1},
	{1, vsopRadSaturn_2},
}

var vsopLonUranus_0 = []vsopTerm{
	{5.48129294297, 0.00000000000, 0.00000000000},   // 0
	{0.09260408234, 0.89106421507, 74.78159856730},  // 1
	{0.01504247898, 3.62719260920, 1.48447270830},   // 2
	{0.00365981674, 1.89962179044, 73.29712585900},  // 3
	{0.00272328168, 3.35823706307, 149.56319713460}, // 4
	{0.00070328461, 5.39254450063, 63.73589830340},  // 5
	{0.00068892678, 6.09292483287, 76.26607127560},  // 6
	{0.00061998615, 2.26952066061, 2.96894541660},   // 7
	{0.00061950719, 2.85098872691, 11.04570026390},  // 8
	{0.00026468770, 3.14152083966, 71.81265315070},  // 9
	{0.00025710476, 6.11379840493, 454.90936652730}, // 10
	{0.00021078850, 4.36059339067, 148.07872442630}, // 11
	{0.00017818647, 1.74436930289, 36.64856292950},  // 12
	{0.00014613507, 4.73732166022, 3.93215326310},   // 13
	{0.00011162509, 5.82681796350, 224.34479570190}, // 14
	{0.00010997910, 0.48865004018, 138.51749687070}, // 15
	{0.00009527478, 2.95516862826, 35.16409022120},  // 16
	{0.00007545601, 5.23626582400, 109.94568878850}, // 17
	{0.00004220241, 3.23328220918, 70.84944530420},  // 18
	{0.00004051900, 2.27755017300, 151.04766984290}, // 19
	{0.00003354596, 1.06549007380, 4.45341812490},   // 20
	{0.00002926718, 4.62903718891, 9.56122755560},   // 21
	{0.00003490340, 5.48306144511, 146.59425171800}, // 22
	{0.00003144069, 4.75199570434, 77.75054398390},  // 23
	{0.00002922333, 5.35235361027, 85.82729883120},  // 24
	{0.00002272788, 4.36600400036, 70.32818044240},  // 25
	{0.00002051219, 1.51773566586, 0.11187458460},   // 26
	{0.00002148602, 0.60745949945, 38.13303563780},  // 27
	{0.00001991643, 4.92437588682, 277.03499374140}, // 28
	{0.00001376226, 2.04283539351, 65.22037101170},  // 29
	{0.00001666902, 3.62744066769, 380.12776796000}, // 30
	{0.00001284107, 3.11347961505, 202.25339517410}, // 31
	{0.00001150429, 0.93343589092, 3.18139373770},   // 32
	{0.00001533221, 2.58594681212, 52.69019803950},  // 33
	{0.00001281604, 0.54271272721, 222.86032299360}, // 34
	{0.00001372139, 4.19641530878, 111.43016149680}, // 35
	{0.00001221029, 0.19900650030, 108.46121608020}, // 36
	{0.00000946181, 1.19253165736, 127.47179660680}, // 37
	{0.00001150989, 4.17898916639, 33.67961751290},  // 38
}

var vsopLonUranus_1 = []vsopTerm{
	{74.78159860910, 0.00000000000, 0.00000000000},  // 0
	{0.00154332863, 5.24158770553, 74.78159856730},  // 1
	{0.00024456474, 1.71260334156, 1.48447270830},   // 2
	{0.00009258442, 0.42829732350, 11.04570026390},  // 3
	{0.00008265977, 1.50218091379, 63.73589830340},  // 4
	{0.00009150160, 1.41213765216, 149.56319713460}, // 5
}

var vsopLonUranus = []vsopSeries{
	{39, vsopLonUranus_0},
	{6, vsopLonUranus_1},
}

var vsopLatUranus_0 = []vsopTerm{
	{0.01346277648, 2.61877810547, 74.78159856730},  // 0
	{0.00062341400, 5.08111189648, 149.56319713460}, // 1
	{0.00061601196, 3.14159265359, 0.00000000000},   // 2
	{0.00009963722, 1.61603805646, 76.26607127560},  // 3
	{0.00009926160, 0.57630380333, 73.29712585900},  // 4
}

var vsopLatUranus_1 = []vsopTerm{
	{0.00034101978, 0.01321929936, 74.78159856730}, // 0
}

var vsopLatUranus = []vsopSeries{
	{5, vsopLatUranus_0},
	{1, vsopLatUranus_1},
}

var vsopRadUranus_0 = []vsopTerm{
	{19.21264847206, 0.00000000000, 0.00000000000},  // 0
	{0.88784984413, 5.60377527014, 74.78159856730},  // 1
	{0.03440836062, 0.32836099706, 73.29712585900},  // 2
	{0.02055653860, 1.78295159330, 149.56319713460}, // 3
	{0.00649322410, 4.52247285911, 76.26607127560},  // 4
	{0.00602247865, 3.86003823674, 63.73589830340},  // 5
	{0.00496404167, 1.40139935333, 454.90936652730}, // 6
	{0.00338525369, 1.58002770318, 138.51749687070}, // 7
	{0.00243509114, 1.57086606044, 71.81265315070},  // 8
	{0.00190522303, 1.99809394714, 1.48447270830},   // 9
	{0.00161858838, 2.79137786799, 148.07872442630}, // 10
	{0.00143706183, 1.38368544947, 11.04570026390},  // 11
	{0.00093192405, 0.17437220467, 36.64856292950},  // 12
	{0.00071424548, 4.24509236074, 224.34479570190}, // 13
	{0.00089806014, 3.66105364565, 109.94568878850}, // 14
	{0.00039009723, 1.66971401684, 70.84944530420},  // 15
	{0.00046677296, 1.39976401694, 35.16409022120},  // 16
	{0.00039025624, 3.36234773834, 277.03499374140}, // 17
	{0.00036755274, 3.88649278513, 146.59425171800}, // 18
	{0.00030348723, 0.70100838798, 151.04766984290}, // 19
	{0.00029156413, 3.18056336700, 77.75054398390},  // 20
	{0.00022637073, 0.72518687029, 529.69096509460}, // 21
	{0.00011959076, 1.75043392140, 984.60033162190}, // 22
	{0.00025620756, 5.25656086672, 380.12776796000}, // 23
}

var vsopRadUranus_1 = []vsopTerm{
	{0.01479896629, 3.67205697578, 74.78159856730}, // 0
}

var vsopRadUranus = []vsopSeries{
	{24, vsopRadUranus_0},
	{1, vsopRadUranus_1},
}

var vsopLonNeptune_0 = []vsopTerm{
	{5.31188633046, 0.00000000000, 0.00000000000},   // 0
	{0.01798475530, 2.90101273890, 38.13303563780},  // 1
	{0.01019727652, 0.48580922867, 1.48447270830},   // 2
	{0.00124531845, 4.83008090676, 36.64856292950},  // 3
	{0.00042064466, 5.41054993053, 2.96894541660},   // 4
	{0.00037714584, 6.09221808686, 35.16409022120},  // 5
	{0.00033784738, 1.24488874087, 76.26607127560},  // 6
	{0.00016482741, 0.00007727998, 491.55792945680}, // 7
	{0.00009198584, 4.93747051954, 39.61750834610},  // 8
	{0.00008994250, 0.27462171806, 175.16605980020}, // 9
}

var vsopLonNeptune_1 = []vsopTerm{
	{38.13303563957, 0.00000000000, 0.00000000000}, // 0
	{0.00016604172, 4.86323329249, 1.48447270830},  // 1
	{0.00015744045, 2.27887427527, 38.13303563780}, // 2
}

var vsopLonNeptune = []vsopSeries{
	{10, vsopLonNeptune_0},
	{3, vsopLonNeptune_1},
}

var vsopLatNeptune_0 = []vsopTerm{
	{0.03088622933, 1.44104372644, 38.13303563780}, // 0
	{0.00027780087, 5.91271884599, 76.26607127560}, // 1
	{0.00027623609, 0.00000000000, 0.00000000000},  // 2
	{0.00015355489, 2.52123799551, 36.64856292950}, // 3
	{0.00015448133, 3.50877079215, 39.61750834610}, // 4
}

var vsopLatNeptune = []vsopSeries{
	{5, vsopLatNeptune_0},
}

var vsopRadNeptune_0 = []vsopTerm{
	{30.07013205828, 0.00000000000, 0.00000000000},   // 0
	{0.27062259632, 1.32999459377, 38.13303563780},   // 1
	{0.01691764014, 3.25186135653, 36.64856292950},   // 2
	{0.00807830553, 5.18592878704, 1.48447270830},    // 3
	{0.00537760510, 4.52113935896, 35.16409022120},   // 4
	{0.00495725141, 1.57105641650, 491.55792945680},  // 5
	{0.00274571975, 1.84552258866, 175.16605980020},  // 6
	{0.00012012320, 1.92059384991, 1021.24889455140}, // 7
	{0.00121801746, 5.79754470298, 76.26607127560},   // 8
	{0.00100896068, 0.37702724930, 73.29712585900},   // 9
	{0.00135134092, 3.37220609835, 39.61750834610},   // 10
	{0.00007571796, 1.07149207335, 388.46515523820},  // 11
}

var vsopRadNeptune = []vsopSeries{
	{12, vsopRadNeptune_0},
}

var vsop = []vsopModel{
	{[3]vsopFormula{{len(vsopLonMercury), vsopLonMercury}, {len(vsopLatMercury), vsopLatMercury}, {len(vsopRadMercury), vsopRadMercury}}},
	{[3]vsopFormula{{len(vsopLonVenus), vsopLonVenus}, {len(vsopLatVenus), vsopLatVenus}, {len(vsopRadVenus), vsopRadVenus}}},
	{[3]vsopFormula{{len(vsopLonEarth), vsopLonEarth}, {len(vsopLatEarth), vsopLatEarth}, {len(vsopRadEarth), vsopRadEarth}}},
	{[3]vsopFormula{{len(vsopLonMars), vsopLonMars}, {len(vsopLatMars), vsopLatMars}, {len(vsopRadMars), vsopRadMars}}},
	{[3]vsopFormula{{len(vsopLonJupiter), vsopLonJupiter}, {len(vsopLatJupiter), vsopLatJupiter}, {len(vsopRadJupiter), vsopRadJupiter}}},
	{[3]vsopFormula{{len(vsopLonSaturn), vsopLonSaturn}, {len(vsopLatSaturn), vsopLatSaturn}, {len(vsopRadSaturn), vsopRadSaturn}}},
	{[3]vsopFormula{{len(vsopLonUranus), vsopLonUranus}, {len(vsopLatUranus), vsopLatUranus}, {len(vsopRadUranus), vsopRadUranus}}},
	{[3]vsopFormula{{len(vsopLonNeptune), vsopLonNeptune}, {len(vsopLatNeptune), vsopLatNeptune}, {len(vsopRadNeptune), vsopRadNeptune}}},
}

var constelNames = [...]constelInfo{
	{"And", "Andromeda"},
	{"Ant", "Antila"},
	{"Aps", "Apus"},
	{"Aql", "Aquila"},
	{"Aqr", "Aquarius"},
	{"Ara", "Ara"},
	{"Ari", "Aries"},
	{"Aur", "Auriga"},
	{"Boo", "Bootes"},
	{"Cae", "Caelum"},
	{"Cam", "Camelopardis"},
	{"Cap", "Capricornus"},
	{"Car", "Carina"},
	{"Cas", "Cassiopeia"},
	{"Cen", "Centaurus"},
	{"Cep", "Cepheus"},
	{"Cet", "Cetus"},
	{"Cha", "Chamaeleon"},
	{"Cir", "Circinus"},
	{"CMa", "Canis Major"},
	{"CMi", "Canis Minor"},
	{"Cnc", "Cancer"},
	{"Col", "Columba"},
	{"Com", "Coma Berenices"},
	{"CrA", "Corona Australis"},
	{"CrB", "Corona Borealis"},
	{"Crt", "Crater"},
	{"Cru", "Crux"},
	{"Crv", "Corvus"},
	{"CVn", "Canes Venatici"},
	{"Cyg", "Cygnus"},
	{"Del", "Delphinus"},
	{"Dor", "Dorado"},
	{"Dra", "Draco"},
	{"Equ", "Equuleus"},
	{"Eri", "Eridanus"},
	{"For", "Fornax"},
	{"Gem", "Gemini"},
	{"Gru", "Grus"},
	{"Her", "Hercules"},
	{"Hor", "Horologium"},
	{"Hya", "Hydra"},
	{"Hyi", "Hydrus"},
	{"Ind", "Indus"},
	{"Lac", "Lacerta"},
	{"Leo", "Leo"},
	{"Lep", "Lepus"},
	{"Lib", "Libra"},
	{"LMi", "Leo Minor"},
	{"Lup", "Lupus"},
	{"Lyn", "Lynx"},
	{"Lyr", "Lyra"},
	{"Men", "Mensa"},
	{"Mic", "Microscopium"},
	{"Mon", "Monoceros"},
	{"Mus", "Musca"},
	{"Nor", "Norma"},
	{"Oct", "Octans"},
	{"Oph", "Ophiuchus"},
	{"Ori", "Orion"},
	{"Pav", "Pavo"},
	{"Peg", "Pegasus"},
	{"Per", "Perseus"},
	{"Phe", "Phoenix"},
	{"Pic", "Pictor"},
	{"PsA", "Pisces Austrinus"},
	{"Psc", "Pisces"},
	{"Pup", "Puppis"},
	{"Pyx", "Pyxis"},
	{"Ret", "Reticulum"},
	{"Scl", "Sculptor"},
	{"Sco", "Scorpius"},
	{"Sct", "Scutum"},
	{"Ser", "Serpens"},
	{"Sex", "Sextans"},
	{"Sge", "Sagitta"},
	{"Sgr", "Sagittarius"},
	{"Tau", "Taurus"},
	{"Tel", "Telescopium"},
	{"TrA", "Triangulum Australe"},
	{"Tri", "Triangulum"},
	{"Tuc", "Tucana"},
	{"UMa", "Ursa Major"},
	{"UMi", "Ursa Minor"},
	{"Vel", "Vela"},
	{"Vir", "Virgo"},
	{"Vol", "Volans"},
	{"Vul", "Vulpecula"},
}

var constelBounds = [...]constelBoundary{
	{83, 0, 8640, 2112},       // UMi
	{83, 2880, 5220, 2076},    // UMi
	{83, 7560, 8280, 2068},    // UMi
	{83, 6480, 7560, 2064},    // UMi
	{15, 0, 2880, 2040},       // Cep
	{10, 3300, 3840, 1968},    // Cam
	{15, 0, 1800, 1920},       // Cep
	{10, 3840, 5220, 1920},    // Cam
	{83, 6300, 6480, 1920},    // UMi
	{33, 7260, 7560, 1920},    // Dra
	{15, 0, 1263, 1848},       // Cep
	{10, 4140, 4890, 1848},    // Cam
	{83, 5952, 6300, 1800},    // UMi
	{15, 7260, 7440, 1800},    // Cep
	{10, 2868, 3300, 1764},    // Cam
	{33, 3300, 4080, 1764},    // Dra
	{83, 4680, 5952, 1680},    // UMi
	{13, 1116, 1230, 1632},    // Cas
	{33, 7350, 7440, 1608},    // Dra
	{33, 4080, 4320, 1596},    // Dra
	{15, 0, 120, 1584},        // Cep
	{83, 5040, 5640, 1584},    // UMi
	{15, 8490, 8640, 1584},    // Cep
	{33, 4320, 4860, 1536},    // Dra
	{33, 4860, 5190, 1512},    // Dra
	{15, 8340, 8490, 1512},    // Cep
	{10, 2196, 2520, 1488},    // Cam
	{33, 7200, 7350, 1476},    // Dra
	{15, 7393.2, 7416, 1462},  // Cep
	{10, 2520, 2868, 1440},    // Cam
	{82, 2868, 3030, 1440},    // UMa
	{33, 7116, 7200, 1428},    // Dra
	{15, 7200, 7393.2, 1428},  // Cep
	{15, 8232, 8340, 1418},    // Cep
	{13, 0, 876, 1404},        // Cas
	{33, 6990, 7116, 1392},    // Dra
	{13, 612, 687, 1380},      // Cas
	{13, 876, 1116, 1368},     // Cas
	{10, 1116, 1140, 1368},    // Cam
	{15, 8034, 8232, 1350},    // Cep
	{10, 1800, 2196, 1344},    // Cam
	{82, 5052, 5190, 1332},    // UMa
	{33, 5190, 6990, 1332},    // Dra
	{10, 1140, 1200, 1320},    // Cam
	{15, 7968, 8034, 1320},    // Cep
	{15, 7416, 7908, 1316},    // Cep
	{13, 0, 612, 1296},        // Cas
	{50, 2196, 2340, 1296},    // Lyn
	{82, 4350, 4860, 1272},    // UMa
	{33, 5490, 5670, 1272},    // Dra
	{15, 7908, 7968, 1266},    // Cep
	{10, 1200, 1800, 1260},    // Cam
	{13, 8232, 8400, 1260},    // Cas
	{33, 5670, 6120, 1236},    // Dra
	{62, 735, 906, 1212},      // Per
	{33, 6120, 6564, 1212},    // Dra
	{13, 0, 492, 1200},        // Cas
	{62, 492, 600, 1200},      // Per
	{50, 2340, 2448, 1200},    // Lyn
	{13, 8400, 8640, 1200},    // Cas
	{82, 4860, 5052, 1164},    // UMa
	{13, 0, 402, 1152},        // Cas
	{13, 8490, 8640, 1152},    // Cas
	{39, 6543, 6564, 1140},    // Her
	{33, 6564, 6870, 1140},    // Dra
	{30, 6870, 6900, 1140},    // Cyg
	{62, 600, 735, 1128},      // Per
	{82, 3030, 3300, 1128},    // UMa
	{13, 60, 312, 1104},       // Cas
	{82, 4320, 4350, 1080},    // UMa
	{50, 2448, 2652, 1068},    // Lyn
	{30, 7887, 7908, 1056},    // Cyg
	{30, 7875, 7887, 1050},    // Cyg
	{30, 6900, 6984, 1044},    // Cyg
	{82, 3300, 3660, 1008},    // UMa
	{82, 3660, 3882, 960},     // UMa
	{8, 5556, 5670, 960},      // Boo
	{39, 5670, 5880, 960},     // Her
	{50, 3330, 3450, 954},     // Lyn
	{0, 0, 906, 882},          // And
	{62, 906, 924, 882},       // Per
	{51, 6969, 6984, 876},     // Lyr
	{62, 1620, 1689, 864},     // Per
	{30, 7824, 7875, 864},     // Cyg
	{44, 7875, 7920, 864},     // Lac
	{7, 2352, 2652, 852},      // Aur
	{50, 2652, 2790, 852},     // Lyn
	{0, 0, 720, 840},          // And
	{44, 7920, 8214, 840},     // Lac
	{44, 8214, 8232, 828},     // Lac
	{0, 8232, 8460, 828},      // And
	{62, 924, 978, 816},       // Per
	{82, 3882, 3960, 816},     // UMa
	{29, 4320, 4440, 816},     // CVn
	{50, 2790, 3330, 804},     // Lyn
	{48, 3330, 3558, 804},     // LMi
	{0, 258, 507, 792},        // And
	{8, 5466, 5556, 792},      // Boo
	{0, 8460, 8550, 770},      // And
	{29, 4440, 4770, 768},     // CVn
	{0, 8550, 8640, 752},      // And
	{29, 5025, 5052, 738},     // CVn
	{80, 870, 978, 736},       // Tri
	{62, 978, 1620, 736},      // Per
	{7, 1620, 1710, 720},      // Aur
	{51, 6543, 6969, 720},     // Lyr
	{82, 3960, 4320, 696},     // UMa
	{30, 7080, 7530, 696},     // Cyg
	{7, 1710, 2118, 684},      // Aur
	{48, 3558, 3780, 684},     // LMi
	{29, 4770, 5025, 684},     // CVn
	{0, 0, 24, 672},           // And
	{80, 507, 600, 672},       // Tri
	{7, 2118, 2352, 672},      // Aur
	{37, 2838, 2880, 672},     // Gem
	{30, 7530, 7824, 672},     // Cyg
	{30, 6933, 7080, 660},     // Cyg
	{80, 690, 870, 654},       // Tri
	{25, 5820, 5880, 648},     // CrB
	{8, 5430, 5466, 624},      // Boo
	{25, 5466, 5820, 624},     // CrB
	{51, 6612, 6792, 624},     // Lyr
	{48, 3870, 3960, 612},     // LMi
	{51, 6792, 6933, 612},     // Lyr
	{80, 600, 690, 600},       // Tri
	{66, 258, 306, 570},       // Psc
	{48, 3780, 3870, 564},     // LMi
	{87, 7650, 7710, 564},     // Vul
	{77, 2052, 2118, 548},     // Tau
	{0, 24, 51, 528},          // And
	{73, 5730, 5772, 528},     // Ser
	{37, 2118, 2238, 516},     // Gem
	{87, 7140, 7290, 510},     // Vul
	{87, 6792, 6930, 506},     // Vul
	{0, 51, 306, 504},         // And
	{87, 7290, 7404, 492},     // Vul
	{37, 2811, 2838, 480},     // Gem
	{87, 7404, 7650, 468},     // Vul
	{87, 6930, 7140, 460},     // Vul
	{6, 1182, 1212, 456},      // Ari
	{75, 6792, 6840, 444},     // Sge
	{59, 2052, 2076, 432},     // Ori
	{37, 2238, 2271, 420},     // Gem
	{75, 6840, 7140, 388},     // Sge
	{77, 1788, 1920, 384},     // Tau
	{39, 5730, 5790, 384},     // Her
	{75, 7140, 7290, 378},     // Sge
	{77, 1662, 1788, 372},     // Tau
	{77, 1920, 2016, 372},     // Tau
	{23, 4620, 4860, 360},     // Com
	{39, 6210, 6570, 344},     // Her
	{23, 4272, 4620, 336},     // Com
	{37, 2700, 2811, 324},     // Gem
	{39, 6030, 6210, 308},     // Her
	{61, 0, 51, 300},          // Peg
	{77, 2016, 2076, 300},     // Tau
	{37, 2520, 2700, 300},     // Gem
	{61, 7602, 7680, 300},     // Peg
	{37, 2271, 2496, 288},     // Gem
	{39, 6570, 6792, 288},     // Her
	{31, 7515, 7578, 284},     // Del
	{61, 7578, 7602, 284},     // Peg
	{45, 4146, 4272, 264},     // Leo
	{59, 2247, 2271, 240},     // Ori
	{37, 2496, 2520, 240},     // Gem
	{21, 2811, 2853, 240},     // Cnc
	{61, 8580, 8640, 240},     // Peg
	{6, 600, 1182, 238},       // Ari
	{31, 7251, 7308, 204},     // Del
	{8, 4860, 5430, 192},      // Boo
	{61, 8190, 8580, 180},     // Peg
	{21, 2853, 3330, 168},     // Cnc
	{45, 3330, 3870, 168},     // Leo
	{58, 6570, 6718.4, 150},   // Oph
	{3, 6718.4, 6792, 150},    // Aql
	{31, 7500, 7515, 144},     // Del
	{20, 2520, 2526, 132},     // CMi
	{73, 6570, 6633, 108},     // Ser
	{39, 5790, 6030, 96},      // Her
	{58, 6570, 6633, 72},      // Oph
	{61, 7728, 7800, 66},      // Peg
	{66, 0, 720, 48},          // Psc
	{73, 6690, 6792, 48},      // Ser
	{31, 7308, 7500, 48},      // Del
	{34, 7500, 7680, 48},      // Equ
	{61, 7680, 7728, 48},      // Peg
	{61, 7920, 8190, 48},      // Peg
	{61, 7800, 7920, 42},      // Peg
	{20, 2526, 2592, 36},      // CMi
	{77, 1290, 1662, 0},       // Tau
	{59, 1662, 1680, 0},       // Ori
	{20, 2592, 2910, 0},       // CMi
	{85, 5280, 5430, 0},       // Vir
	{58, 6420, 6570, 0},       // Oph
	{16, 954, 1182, -42},      // Cet
	{77, 1182, 1290, -42},     // Tau
	{73, 5430, 5856, -78},     // Ser
	{59, 1680, 1830, -96},     // Ori
	{59, 2100, 2247, -96},     // Ori
	{73, 6420, 6468, -96},     // Ser
	{73, 6570, 6690, -96},     // Ser
	{3, 6690, 6792, -96},      // Aql
	{66, 8190, 8580, -96},     // Psc
	{45, 3870, 4146, -144},    // Leo
	{85, 4146, 4260, -144},    // Vir
	{66, 0, 120, -168},        // Psc
	{66, 8580, 8640, -168},    // Psc
	{85, 5130, 5280, -192},    // Vir
	{58, 5730, 5856, -192},    // Oph
	{3, 7200, 7392, -216},     // Aql
	{4, 7680, 7872, -216},     // Aqr
	{58, 6180, 6468, -240},    // Oph
	{54, 2100, 2910, -264},    // Mon
	{35, 1770, 1830, -264},    // Eri
	{59, 1830, 2100, -264},    // Ori
	{41, 2910, 3012, -264},    // Hya
	{74, 3450, 3870, -264},    // Sex
	{85, 4260, 4620, -264},    // Vir
	{58, 6330, 6360, -280},    // Oph
	{3, 6792, 7200, -288.8},   // Aql
	{35, 1740, 1770, -348},    // Eri
	{4, 7392, 7680, -360},     // Aqr
	{73, 6180, 6570, -384},    // Ser
	{72, 6570, 6792, -384},    // Sct
	{41, 3012, 3090, -408},    // Hya
	{58, 5856, 5895, -438},    // Oph
	{41, 3090, 3270, -456},    // Hya
	{26, 3870, 3900, -456},    // Crt
	{71, 5856, 5895, -462},    // Sco
	{47, 5640, 5730, -480},    // Lib
	{28, 4530, 4620, -528},    // Crv
	{85, 4620, 5130, -528},    // Vir
	{41, 3270, 3510, -576},    // Hya
	{16, 600, 954, -585.2},    // Cet
	{35, 954, 1350, -585.2},   // Eri
	{26, 3900, 4260, -588},    // Crt
	{28, 4260, 4530, -588},    // Crv
	{47, 5130, 5370, -588},    // Lib
	{58, 5856, 6030, -590},    // Oph
	{16, 0, 600, -612},        // Cet
	{11, 7680, 7872, -612},    // Cap
	{4, 7872, 8580, -612},     // Aqr
	{16, 8580, 8640, -612},    // Cet
	{41, 3510, 3690, -636},    // Hya
	{35, 1692, 1740, -654},    // Eri
	{46, 1740, 2202, -654},    // Lep
	{11, 7200, 7680, -672},    // Cap
	{41, 3690, 3810, -700},    // Hya
	{41, 4530, 5370, -708},    // Hya
	{47, 5370, 5640, -708},    // Lib
	{71, 5640, 5760, -708},    // Sco
	{35, 1650, 1692, -720},    // Eri
	{58, 6030, 6336, -720},    // Oph
	{76, 6336, 6420, -720},    // Sgr
	{41, 3810, 3900, -748},    // Hya
	{19, 2202, 2652, -792},    // CMa
	{41, 4410, 4530, -792},    // Hya
	{41, 3900, 4410, -840},    // Hya
	{36, 1260, 1350, -864},    // For
	{68, 3012, 3372, -882},    // Pyx
	{35, 1536, 1650, -888},    // Eri
	{76, 6420, 6900, -888},    // Sgr
	{65, 7680, 8280, -888},    // PsA
	{70, 8280, 8400, -888},    // Scl
	{36, 1080, 1260, -950},    // For
	{1, 3372, 3960, -954},     // Ant
	{70, 0, 600, -960},        // Scl
	{36, 600, 1080, -960},     // For
	{35, 1392, 1536, -960},    // Eri
	{70, 8400, 8640, -960},    // Scl
	{14, 5100, 5370, -1008},   // Cen
	{49, 5640, 5760, -1008},   // Lup
	{71, 5760, 5911.5, -1008}, // Sco
	{9, 1740, 1800, -1032},    // Cae
	{22, 1800, 2370, -1032},   // Col
	{67, 2880, 3012, -1032},   // Pup
	{35, 1230, 1392, -1056},   // Eri
	{71, 5911.5, 6420, -1092}, // Sco
	{24, 6420, 6900, -1092},   // CrA
	{76, 6900, 7320, -1092},   // Sgr
	{53, 7320, 7680, -1092},   // Mic
	{35, 1080, 1230, -1104},   // Eri
	{9, 1620, 1740, -1116},    // Cae
	{49, 5520, 5640, -1152},   // Lup
	{63, 0, 840, -1156},       // Phe
	{35, 960, 1080, -1176},    // Eri
	{40, 1470, 1536, -1176},   // Hor
	{9, 1536, 1620, -1176},    // Cae
	{38, 7680, 7920, -1200},   // Gru
	{67, 2160, 2880, -1218},   // Pup
	{84, 2880, 2940, -1218},   // Vel
	{35, 870, 960, -1224},     // Eri
	{40, 1380, 1470, -1224},   // Hor
	{63, 0, 660, -1236},       // Phe
	{12, 2160, 2220, -1260},   // Car
	{84, 2940, 3042, -1272},   // Vel
	{40, 1260, 1380, -1276},   // Hor
	{32, 1380, 1440, -1276},   // Dor
	{63, 0, 570, -1284},       // Phe
	{35, 780, 870, -1296},     // Eri
	{64, 1620, 1800, -1296},   // Pic
	{49, 5418, 5520, -1296},   // Lup
	{84, 3042, 3180, -1308},   // Vel
	{12, 2220, 2340, -1320},   // Car
	{14, 4260, 4620, -1320},   // Cen
	{49, 5100, 5418, -1320},   // Lup
	{56, 5418, 5520, -1320},   // Nor
	{32, 1440, 1560, -1356},   // Dor
	{84, 3180, 3960, -1356},   // Vel
	{14, 3960, 4050, -1356},   // Cen
	{5, 6300, 6480, -1368},    // Ara
	{78, 6480, 7320, -1368},   // Tel
	{38, 7920, 8400, -1368},   // Gru
	{40, 1152, 1260, -1380},   // Hor
	{64, 1800, 1980, -1380},   // Pic
	{12, 2340, 2460, -1392},   // Car
	{63, 0, 480, -1404},       // Phe
	{35, 480, 780, -1404},     // Eri
	{63, 8400, 8640, -1404},   // Phe
	{32, 1560, 1650, -1416},   // Dor
	{56, 5520, 5911.5, -1440}, // Nor
	{43, 7320, 7680, -1440},   // Ind
	{64, 1980, 2160, -1464},   // Pic
	{18, 5460, 5520, -1464},   // Cir
	{5, 5911.5, 5970, -1464},  // Ara
	{18, 5370, 5460, -1526},   // Cir
	{5, 5970, 6030, -1526},    // Ara
	{64, 2160, 2460, -1536},   // Pic
	{12, 2460, 3252, -1536},   // Car
	{14, 4050, 4260, -1536},   // Cen
	{27, 4260, 4620, -1536},   // Cru
	{14, 4620, 5232, -1536},   // Cen
	{18, 4860, 4920, -1560},   // Cir
	{5, 6030, 6060, -1560},    // Ara
	{40, 780, 1152, -1620},    // Hor
	{69, 1152, 1650, -1620},   // Ret
	{18, 5310, 5370, -1620},   // Cir
	{5, 6060, 6300, -1620},    // Ara
	{60, 6300, 6480, -1620},   // Pav
	{81, 7920, 8400, -1620},   // Tuc
	{32, 1650, 2370, -1680},   // Dor
	{18, 4920, 5310, -1680},   // Cir
	{79, 5310, 6120, -1680},   // TrA
	{81, 0, 480, -1800},       // Tuc
	{42, 1260, 1650, -1800},   // Hyi
	{86, 2370, 3252, -1800},   // Vol
	{12, 3252, 4050, -1800},   // Car
	{55, 4050, 4920, -1800},   // Mus
	{60, 6480, 7680, -1800},   // Pav
	{43, 7680, 8400, -1800},   // Ind
	{81, 8400, 8640, -1800},   // Tuc
	{81, 270, 480, -1824},     // Tuc
	{42, 0, 1260, -1980},      // Hyi
	{17, 2760, 4920, -1980},   // Cha
	{2, 4920, 6480, -1980},    // Aps
	{52, 1260, 2760, -2040},   // Men
	{57, 0, 8640, -2160},      // Oct
}

var moonAddSolTerms = [...]addSolTerm{

	{13.9020, 14.0600, -0.0010, 0.2607, 0, 0, 0, 4},
	{0.4030, -4.0100, 0.3940, 0.0023, 0, 0, 0, 3},
	{2369.9120, 2373.3600, 0.6010, 28.2333, 0, 0, 0, 2},
	{-125.1540, -112.7900, -0.7250, -0.9781, 0, 0, 0, 1},
	{1.9790, 6.9800, -0.4450, 0.0433, 1, 0, 0, 4},
	{191.9530, 192.7200, 0.0290, 3.0861, 1, 0, 0, 2},
	{-8.4660, -13.5100, 0.4550, -0.1093, 1, 0, 0, 1},
	{22639.5000, 22609.0700, 0.0790, 186.5398, 1, 0, 0, 0},
	{18.6090, 3.5900, -0.0940, 0.0118, 1, 0, 0, -1},
	{-4586.4650, -4578.1300, -0.0770, 34.3117, 1, 0, 0, -2},
	{3.2150, 5.4400, 0.1920, -0.0386, 1, 0, 0, -3},
	{-38.4280, -38.6400, 0.0010, 0.6008, 1, 0, 0, -4},
	{-0.3930, -1.4300, -0.0920, 0.0086, 1, 0, 0, -6},
	{-0.2890, -1.5900, 0.1230, -0.0053, 0, 1, 0, 4},
	{-24.4200, -25.1000, 0.0400, -0.3000, 0, 1, 0, 2},
	{18.0230, 17.9300, 0.0070, 0.1494, 0, 1, 0, 1},
	{-668.1460, -126.9800, -1.3020, -0.3997, 0, 1, 0, 0},
	{0.5600, 0.3200, -0.0010, -0.0037, 0, 1, 0, -1},
	{-165.1450, -165.0600, 0.0540, 1.9178, 0, 1, 0, -2},
	{-1.8770, -6.4600, -0.4160, 0.0339, 0, 1, 0, -4},
	{0.2130, 1.0200, -0.0740, 0.0054, 2, 0, 0, 4},
	{14.3870, 14.7800, -0.0170, 0.2833, 2, 0, 0, 2},
	{-0.5860, -1.2000, 0.0540, -0.0100, 2, 0, 0, 1},
	{769.0160, 767.9600, 0.1070, 10.1657, 2, 0, 0, 0},
	{1.7500, 2.0100, -0.0180, 0.0155, 2, 0, 0, -1},
	{-211.6560, -152.5300, 5.6790, -0.3039, 2, 0, 0, -2},
	{1.2250, 0.9100, -0.0300, -0.0088, 2, 0, 0, -3},
	{-30.7730, -34.0700, -0.3080, 0.3722, 2, 0, 0, -4},
	{-0.5700, -1.4000, -0.0740, 0.0109, 2, 0, 0, -6},
	{-2.9210, -11.7500, 0.7870, -0.0484, 1, 1, 0, 2},
	{1.2670, 1.5200, -0.0220, 0.0164, 1, 1, 0, 1},
	{-109.6730, -115.1800, 0.4610, -0.9490, 1, 1, 0, 0},
	{-205.9620, -182.3600, 2.0560, 1.4437, 1, 1, 0, -2},
	{0.2330, 0.3600, 0.0120, -0.0025, 1, 1, 0, -3},
	{-4.3910, -9.6600, -0.4710, 0.0673, 1, 1, 0, -4},
	{0.2830, 1.5300, -0.1110, 0.0060, 1, -1, 0, 4},
	{14.5770, 31.7000, -1.5400, 0.2302, 1, -1, 0, 2},
	{147.6870, 138.7600, 0.6790, 1.1528, 1, -1, 0, 0},
	{-1.0890, 0.5500, 0.0210, 0.0000, 1, -1, 0, -1},
	{28.4750, 23.5900, -0.4430, -0.2257, 1, -1, 0, -2},
	{-0.2760, -0.3800, -0.0060, -0.0036, 1, -1, 0, -3},
	{0.6360, 2.2700, 0.1460, -0.0102, 1, -1, 0, -4},
	{-0.1890, -1.6800, 0.1310, -0.0028, 0, 2, 0, 2},
	{-7.4860, -0.6600, -0.0370, -0.0086, 0, 2, 0, 0},
	{-8.0960, -16.3500, -0.7400, 0.0918, 0, 2, 0, -2},
	{-5.7410, -0.0400, 0.0000, -0.0009, 0, 0, 2, 2},
	{0.2550, 0.0000, 0.0000, 0.0000, 0, 0, 2, 1},
	{-411.6080, -0.2000, 0.0000, -0.0124, 0, 0, 2, 0},
	{0.5840, 0.8400, 0.0000, 0.0071, 0, 0, 2, -1},
	{-55.1730, -52.1400, 0.0000, -0.1052, 0, 0, 2, -2},
	{0.2540, 0.2500, 0.0000, -0.0017, 0, 0, 2, -3},
	{0.0250, -1.6700, 0.0000, 0.0031, 0, 0, 2, -4},
	{1.0600, 2.9600, -0.1660, 0.0243, 3, 0, 0, 2},
	{36.1240, 50.6400, -1.3000, 0.6215, 3, 0, 0, 0},
	{-13.1930, -16.4000, 0.2580, -0.1187, 3, 0, 0, -2},
	{-1.1870, -0.7400, 0.0420, 0.0074, 3, 0, 0, -4},
	{-0.2930, -0.3100, -0.0020, 0.0046, 3, 0, 0, -6},
	{-0.2900, -1.4500, 0.1160, -0.0051, 2, 1, 0, 2},
	{-7.6490, -10.5600, 0.2590, -0.1038, 2, 1, 0, 0},
	{-8.6270, -7.5900, 0.0780, -0.0192, 2, 1, 0, -2},
	{-2.7400, -2.5400, 0.0220, 0.0324, 2, 1, 0, -4},
	{1.1810, 3.3200, -0.2120, 0.0213, 2, -1, 0, 2},
	{9.7030, 11.6700, -0.1510, 0.1268, 2, -1, 0, 0},
	{-0.3520, -0.3700, 0.0010, -0.0028, 2, -1, 0, -1},
	{-2.4940, -1.1700, -0.0030, -0.0017, 2, -1, 0, -2},
	{0.3600, 0.2000, -0.0120, -0.0043, 2, -1, 0, -4},
	{-1.1670, -1.2500, 0.0080, -0.0106, 1, 2, 0, 0},
	{-7.4120, -6.1200, 0.1170, 0.0484, 1, 2, 0, -2},
	{-0.3110, -0.6500, -0.0320, 0.0044, 1, 2, 0, -4},
	{0.7570, 1.8200, -0.1050, 0.0112, 1, -2, 0, 2},
	{2.5800, 2.3200, 0.0270, 0.0196, 1, -2, 0, 0},
	{2.5330, 2.4000, -0.0140, -0.0212, 1, -2, 0, -2},
	{-0.3440, -0.5700, -0.0250, 0.0036, 0, 3, 0, -2},
	{-0.9920, -0.0200, 0.0000, 0.0000, 1, 0, 2, 2},
	{-45.0990, -0.0200, 0.0000, -0.0010, 1, 0, 2, 0},
	{-0.1790, -9.5200, 0.0000, -0.0833, 1, 0, 2, -2},
	{-0.3010, -0.3300, 0.0000, 0.0014, 1, 0, 2, -4},
	{-6.3820, -3.3700, 0.0000, -0.0481, 1, 0, -2, 2},
	{39.5280, 85.1300, 0.0000, -0.7136, 1, 0, -2, 0},
	{9.3660, 0.7100, 0.0000, -0.0112, 1, 0, -2, -2},
	{0.2020, 0.0200, 0.0000, 0.0000, 1, 0, -2, -4},
	{0.4150, 0.1000, 0.0000, 0.0013, 0, 1, 2, 0},
	{-2.1520, -2.2600, 0.0000, -0.0066, 0, 1, 2, -2},
	{-1.4400, -1.3000, 0.0000, 0.0014, 0, 1, -2, 2},
	{0.3840, -0.0400, 0.0000, 0.0000, 0, 1, -2, -2},
	{1.9380, 3.6000, -0.1450, 0.0401, 4, 0, 0, 0},
	{-0.9520, -1.5800, 0.0520, -0.0130, 4, 0, 0, -2},
	{-0.5510, -0.9400, 0.0320, -0.0097, 3, 1, 0, 0},
	{-0.4820, -0.5700, 0.0050, -0.0045, 3, 1, 0, -2},
	{0.6810, 0.9600, -0.0260, 0.0115, 3, -1, 0, 0},
	{-0.2970, -0.2700, 0.0020, -0.0009, 2, 2, 0, -2},
	{0.2540, 0.2100, -0.0030, 0.0000, 2, -2, 0, -2},
	{-0.2500, -0.2200, 0.0040, 0.0014, 1, 3, 0, -2},
	{-3.9960, 0.0000, 0.0000, 0.0004, 2, 0, 2, 0},
	{0.5570, -0.7500, 0.0000, -0.0090, 2, 0, 2, -2},
	{-0.4590, -0.3800, 0.0000, -0.0053, 2, 0, -2, 2},
	{-1.2980, 0.7400, 0.0000, 0.0004, 2, 0, -2, 0},
	{0.5380, 1.1400, 0.0000, -0.0141, 2, 0, -2, -2},
	{0.2630, 0.0200, 0.0000, 0.0000, 1, 1, 2, 0},
	{0.4260, 0.0700, 0.0000, -0.0006, 1, 1, -2, -2},
	{-0.3040, 0.0300, 0.0000, 0.0003, 1, -1, 2, 0},
	{-0.3720, -0.1900, 0.0000, -0.0027, 1, -1, -2, 2},
	{0.4180, 0.0000, 0.0000, 0.0000, 0, 0, 4, 0},
	{-0.3300, -0.0400, 0.0000, 0.0000, 3, 0, 2, 0},
}
