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
	cosObliq2000 = 0.9174821430670688 // cosine of the mean obliquity of the J2000 ecliptic
	sinObliq2000 = 0.3977769691083922 // sine of the mean obliquity of the J2000 ecliptic
	earthMoonMassRatio = 81.30056     // ratio of Earth mass to Moon mass
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
	earthOrbitalPeriod   = 365.256   // days for Earth to orbit the Sun
	neptuneOrbitalPeriod = 60189.0   // days for Neptune to orbit the Sun
	solarDaysPerSiderealDay = 0.9972695717592592
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
	InvalidBody Body  = -1
	Mercury     Body  = 0 // The planet Mercury
	Venus             = 1 // The planet Venus
	Earth             = 2 // The planet Earth
	Mars              = 3 // The planet Mars
	Jupiter           = 4 // The planet Jupiter
	Saturn            = 5 // The planet Saturn
	Uranus            = 6 // The planet Uranus
	Neptune           = 7 // The planet Neptune
	Pluto             = 8 // The dwarf planet Pluto
	Sun               = 9 // The Sun
	Moon              = 10 // The Earth's Moon
	Emb               = 11 // The Earth/Moon Barycenter
	Ssb               = 12 // The Solar System Barycenter
	Star1             = 101 // User-defined star #1
	Star2             = 102 // User-defined star #2
	Star3             = 103 // User-defined star #3
	Star4             = 104 // User-defined star #4
	Star5             = 105 // User-defined star #5
	Star6             = 106 // User-defined star #6
	Star7             = 107 // User-defined star #7
	Star8             = 108 // User-defined star #8
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
	tt float64    // Terrestrial Time in J2000 days
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
	into2000 precessDirection = 1
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
	R := (fb-fa)/2.0
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

// RefractionAngle returns the angular adjustment in degrees to add to the given
// altitude angle to correct for atmospheric refraction. The refraction model
// is selected by the refraction parameter: NoRefraction (or an invalid value)
// returns 0, while NormalRefraction and JplHorizonsRefraction apply the
// Saemundsson model. Altitudes outside [-90, +90] are not corrected.
func RefractionAngle(refraction Refraction, altitude float64) float64 {
	var refr, hd float64
	if altitude < -90.0 || altitude > +90.0 {
		return 0.0 // no attempt to correct an invalid altitude
	}
	if refraction == NormalRefraction || refraction == JplHorizonsRefraction {
		hd = altitude
		if hd < -1.0 {
			hd = -1.0
		}
		refr = (1.02 / math.Tan((hd+10.3/(hd+5.11))*RadiansFromDegrees(1.0))) / 60.0
		if refraction == NormalRefraction && altitude < -1.0 {
			// In "normal" mode we gradually reduce refraction toward the nadir
			// so that we never get an altitude angle less than -90 degrees.
			refr *= (altitude + 90.0) / 89.0
		}
	} else {
		refr = 0.0
	}
	return refr
}

// eclipseKind classifies an eclipse as total, annular, partial, or penumbral.
type EclipseKind int

const (
	EclipseKindTotal     EclipseKind = iota // the Moon completely blocks the Sun (total solar eclipse)
	EclipseKindAnnular                      // the Moon is entirely inside the Sun's disc (annular solar eclipse)
	EclipseKindPartial                      // the Moon only partially covers the Sun
	EclipseKindPenumbral                    // the Moon passes only through the Earth's penumbra
)

// eclipseKindFromUmbra classifies a solar eclipse by the size of the umbra
// radius: positive means the umbra reaches the ground (total), otherwise
// the eclipse is annular.
func eclipseKindFromUmbra(k float64) EclipseKind {
	if k > 0.014 {
		return EclipseKindTotal
	}
	return EclipseKindAnnular
}

// obscuration returns the fractional area of the first disc that is covered
// by the second disc. a and b are the radii of the two discs, and c is the
// distance between their centers, all expressed in the same units.
func obscuration(a, b, c float64) float64 {
	if a <= 0.0 || b <= 0.0 {
		return 0.0 // invalid radius
	}
	if c < 0.0 {
		return 0.0 // invalid distance between centers
	}
	if c >= a+b {
		return 0.0 // the discs are too far apart to overlap
	}
	if c == 0.0 {
		// The discs have a common center: one disc is inside the other.
		if a <= b {
			return 1.0
		}
		return (b * b) / (a * a)
	}
	x := (a*a - b*b + c*c) / (2 * c)
	radicand := a*a - x*x
	if radicand <= 0.0 {
		// The circumferences do not intersect or are tangent; one disc is inside the other.
		if a <= b {
			return 1.0
		}
		return (b * b) / (a * a)
	}
	y := math.Sqrt(radicand)
	// The discs overlap in a pair of lens-shaped areas.
	lens1 := a*a*math.Acos(x/a) - x*y
	lens2 := b*b*math.Acos((c-x)/b) - (c-x)*y
	return (lens1 + lens2) / (math.Pi * a * a)
}

// horizonDipAngle calculates the angle in degrees that the visible horizon
// dips below the observer's horizontal plane, due to the curvature of the
// Earth and the observer's height above the ground.
func horizonDipAngle(observer Observer, metersAboveGround float64) float64 {
	phi := RadiansFromDegrees(observer.Latitude)
	sinphi := math.Sin(phi)
	cosphi := math.Cos(phi)
	c := 1.0 / math.Hypot(cosphi, sinphi*EarthFlattening)
	s := c * (EarthFlattening * EarthFlattening)
	htKm := (observer.Height - metersAboveGround) / 1000.0
	ach := EarthEquatorialRadiusKm*c + htKm
	ash := EarthEquatorialRadiusKm*s + htKm
	radiusM := 1000.0 * math.Hypot(ach*cosphi, ash*sinphi)

	// Correct for refraction of a ray of light tangent to the Earth's surface.
	k := 0.175 * math.Pow(1.0-(6.5e-3/283.15)*(observer.Height-(2.0/3.0)*metersAboveGround), 3.256)
	return DegreesFromRadians(-(math.Sqrt(2*(1-k)*metersAboveGround/radiusM) / (1 - k)))
}

// ObserverGravity returns the acceleration of gravity in m/s^2 at the given
// latitude in degrees and height above sea level in meters.
func ObserverGravity(latitude float64, height float64) float64 {
	s := math.Sin(RadiansFromDegrees(latitude))
	s2 := s * s
	g0 := 9.7803253359 * (1.0 + 0.00193185265241*s2) / math.Sqrt(1.0-0.00669437999013*s2)
	return g0 * (1.0 - (3.15704e-07-2.10269e-09*s2)*height + 7.37452e-14*height*height)
}

// maxAltitudeSlope returns the maximum possible rate of change in degrees per
// day of the body's altitude as seen by an observer at the given latitude.
func maxAltitudeSlope(body Body, latitude float64) (float64, error) {
	if !isfinite(latitude) || latitude < -90.0 || latitude > +90.0 {
		return math.NaN(), errors.New("Invalid latitude")
	}
	var derivRa, derivDec float64
	switch body {
	case Moon:
		derivRa = +4.5
		derivDec = +8.2
	case Sun:
		derivRa = +0.8
		derivDec = +0.5
	case Mercury:
		derivRa = -1.6
		derivDec = +1.0
	case Venus:
		derivRa = -0.8
		derivDec = +0.6
	case Mars:
		derivRa = -0.5
		derivDec = +0.4
	case Jupiter, Saturn, Uranus, Neptune, Pluto:
		derivRa = -0.2
		derivDec = +0.2
	case Earth:
		return math.NaN(), errors.New("Earth not allowed for altitude search")
	default:
		if userDefinedStar(body) != nil {
			// Parallax and aberration are tiny for a star; use a generous buffer.
			derivRa = -0.008
			derivDec = +0.008
			break
		}
		return math.NaN(), errors.New("Invalid body for altitude search")
	}
	latrad := RadiansFromDegrees(latitude)
	value := math.Abs(((360.0 / solarDaysPerSiderealDay) - derivRa)*math.Cos(latrad)) + math.Abs(derivDec*math.Sin(latrad))
	return value, nil
}

// altitudeContext holds the state needed to evaluate the altitude of a body
// relative to a target altitude, used by the rise/set search algorithms.
type altitudeContext struct {
	body           Body
	direction      int
	observer       Observer
	bodyRadiusAu   float64
	targetAltitude float64
}

// ascent holds a time interval [tx, ty] within which the altitude-difference
// function rises from non-positive to positive.
type ascent struct {
	tx AstroTime
	ty AstroTime
	ax float64
	ay float64
}

// findAscent searches the interval [t1, t2] for a sub-interval where the
// altitude-difference function (evaluated by altDiff) rises from non-positive
// to positive. Returns ok=false when no such ascent exists, and a non-nil
// error for actual failures.
func findAscent(depth int, maxDerivAlt float64, t1, t2 AstroTime, a1, a2 float64, altDiff func(time AstroTime) (float64, error)) (ascent, bool, error) {
	if a1 < 0.0 && a2 >= 0.0 {
		// Trivial success: the endpoints already rise through zero.
		return ascent{tx: t1, ty: t2, ax: a1, ay: a2}, true, nil
	}
	if a1 >= 0.0 && a2 < 0.0 {
		// Trivial failure: assume the Nyquist condition prevents an ascent.
		return ascent{}, false, nil
	}
	if depth > 17 {
		// Safety valve: do not allow unlimited recursion.
		return ascent{}, false, errors.New("Maximum recursion depth exceeded in rise/set search")
	}
	dt := (t2.Ut - t1.Ut) / 2
	if dt*86400.0 < 1.0 {
		return ascent{}, false, nil
	}
	absA1 := math.Abs(a1)
	absA2 := math.Abs(a2)
	da := absA1
	if absA2 < absA1 {
		da = absA2
	}
	if da > maxDerivAlt*(dt/2) {
		// Prune: the altitude cannot change fast enough to reach zero.
		return ascent{}, false, nil
	}
	tm := TimeFromUniversalDays((t1.Ut + t2.Ut) / 2)
	altValue, err := altDiff(tm)
	if err != nil {
		return ascent{}, false, err
	}
	// Recurse into the left interval, then the right.
	left, ok, err := findAscent(depth+1, maxDerivAlt, t1, tm, a1, altValue, altDiff)
	if err != nil {
		return ascent{}, false, err
	}
	if !ok {
		return findAscent(depth+1, maxDerivAlt, tm, t2, altValue, a2, altDiff)
	}
	return left, true, nil
}

// moonMagnitude estimates the visual magnitude of the Moon from its phase
// angle and its heliocentric and geocentric distances.
func moonMagnitude(phase float64, helioDist float64, geoDist float64) float64 {
	rad := RadiansFromDegrees(phase)
	rad2 := rad * rad
	rad4 := rad2 * rad2
	mag := -12.717 + 1.49*math.Abs(rad) + 0.0431*rad4
	moonMeanDistanceAu := 385000.6 / KmPerAu
	geoAu := geoDist / moonMeanDistanceAu
	mag += 5 * math.Log10(helioDist*geoAu)
	return mag
}

// visualMagnitude estimates the visual magnitude of a planet from its phase
// angle and its heliocentric and geocentric distances.
func visualMagnitude(body Body, phase float64, helioDist float64, geoDist float64) (float64, error) {
	// For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212
	var c0, c1, c2, c3 float64
	switch body {
	case Mercury:
		c0, c1, c2, c3 = -0.60, +4.98, -4.88, +3.02
	case Venus:
		if phase < 163.6 {
			c0, c1, c2, c3 = -4.47, +1.03, +0.57, +0.13
		} else {
			c0, c1 = 0.98, -1.02
		}
	case Mars:
		c0, c1 = -1.52, +1.60
	case Jupiter:
		c0, c1 = -9.40, +0.50
	case Uranus:
		c0, c1 = -7.19, +0.25
	case Neptune:
		c0 = -6.87
	case Pluto:
		c0, c1 = -1.00, +4.00
	default:
		return 0.0, errors.New("Invalid body for visual magnitude")
	}
	x := phase / 100
	mag := c0 + x*(c1+x*(c2+x*c3))
	mag += 5.0 * math.Log10(helioDist*geoDist)
	return mag, nil
}

// jupiterMoonElem2pv converts a set of orbital elements of a Galilean moon
// into a Cartesian state vector, translated from the FORTRAN routine ELEM2PV
// of the L1.2 ephemeris model.
func jupiterMoonElem2pv(time AstroTime, mu float64, elem [6]float64) StateVector {
	a := elem[0]
	al := elem[1]
	k := elem[2]
	h := elem[3]
	q := elem[4]
	p := elem[5]

	an := math.Sqrt(mu / (a * a * a))

	ee := al + k*math.Sin(al) - h*math.Cos(al)
	for {
		ce := math.Cos(ee)
		se := math.Sin(ee)
		de := (al - ee + k*se - h*ce) / (1.0 - k*ce - h*se)
		ee += de
		if math.Abs(de) < 1.0e-12 {
			break
		}
	}

	ce := math.Cos(ee)
	se := math.Sin(ee)
	dle := h*ce - k*se
	rsam1 := -k*ce - h*se
	asr := 1.0 / (1.0 + rsam1)
	phi := math.Sqrt(1.0 - k*k - h*h)
	psi := 1.0 / (1.0 + phi)
	x1 := a * (ce - k - psi*h*dle)
	y1 := a * (se - h + psi*k*dle)
	vx1 := an * asr * a * (-se - psi*h*rsam1)
	vy1 := an * asr * a * (+ce + psi*k*rsam1)
	f2 := 2.0 * math.Sqrt(1.0-q*q-p*p)
	p2 := 1.0 - 2.0*p*p
	q2 := 1.0 - 2.0*q*q
	pq := 2.0 * p * q

	return StateVector{
		X:  x1*p2 + y1*pq,
		Y:  x1*pq + y1*q2,
		Z:  (q*y1 - x1*p) * f2,
		Vx: vx1*p2 + vy1*pq,
		Vy: vx1*pq + vy1*q2,
		Vz: (q*vy1 - vx1*p) * f2,
		T:  time,
	}
}

// shadow describes the alignment of a shadow-casting body with a target body,
// used for finding eclipses.
type shadow struct {
	time   AstroTime
	u      float64 // dot product of the shadow direction and target position
	r      float64 // km distance between the target and the shadow axis
	k      float64 // umbra radius in km at the shadow plane
	p      float64 // penumbra radius in km at the shadow plane
	target AstroVector
	dir    AstroVector
}

// calcShadow calculates the geometry of the shadow cast by one body onto
// another, given the casting body's radius and the target and direction
// vectors expressed relative to the casting body.
func calcShadow(bodyRadiusKm float64, time AstroTime, target AstroVector, dir AstroVector) shadow {
	s := shadow{
		time:   time,
		target: target,
		dir:    dir,
	}
	s.u = (dir.X*target.X + dir.Y*target.Y + dir.Z*target.Z) / (dir.X*dir.X + dir.Y*dir.Y + dir.Z*dir.Z)
	dx := (s.u * dir.X) - target.X
	dy := (s.u * dir.Y) - target.Y
	dz := (s.u * dir.Z) - target.Z
	s.r = KmPerAu * math.Sqrt(dx*dx+dy*dy+dz*dz)
	s.k = +SunRadiusKm - (1.0+s.u)*(SunRadiusKm-bodyRadiusKm)
	s.p = -SunRadiusKm + (1.0+s.u)*(SunRadiusKm+bodyRadiusKm)
	return s
}

// calcVsop evaluates the VSOP87 model for a body at the given time, returning
// the heliocentric equatorial position vector.
func calcVsop(model *vsopModel, time AstroTime) AstroVector {
	t := time.Tt / daysPerMillennium
	sphere := vsopFormulaCalc(model, t)
	eclip := vsopSphereToRect(sphere[vsopLonIndex], sphere[vsopLatIndex], sphere[vsopRadIndex])
	pos := vsopRotate(eclip)
	return AstroVector{X: pos.X, Y: pos.Y, Z: pos.Z, T: time}
}

// calcVsopPosVel evaluates the VSOP87 model for a body at the given Terrestrial
// Time in days since J2000, returning both position and velocity.
func calcVsopPosVel(model *vsopModel, tt float64) bodyState {
	t := tt / daysPerMillennium
	sphere := vsopFormulaCalc(model, t)
	deriv := vsopDerivCalc(model, t)
	eclip := vsopSphereToRect(sphere[vsopLonIndex], sphere[vsopLatIndex], sphere[vsopRadIndex])

	state := bodyState{tt: tt}
	state.r = vsopRotate(eclip)

	// Use spherical coordinates and spherical derivatives to calculate the
	// velocity vector in rectangular coordinates.
	coslon := math.Cos(sphere[vsopLonIndex])
	sinlon := math.Sin(sphere[vsopLonIndex])
	coslat := math.Cos(sphere[vsopLatIndex])
	sinlat := math.Sin(sphere[vsopLatIndex])
	r := sphere[vsopRadIndex]
	dlonDt := deriv[vsopLonIndex]
	dlatDt := deriv[vsopLatIndex]
	drDt := deriv[vsopRadIndex]

	eclip[0] = (drDt*coslat*coslon) - (r*sinlat*coslon*dlatDt) - (r*coslat*sinlon*dlonDt)
	eclip[1] = (drDt*coslat*sinlon) - (r*sinlat*sinlon*dlatDt) + (r*coslat*coslon*dlonDt)
	eclip[2] = (drDt * sinlat) + (r * coslat * dlatDt)

	state.v = vsopRotate(eclip)
	// Convert speed units from [AU/millennium] to [AU/day].
	state.v.X /= daysPerMillennium
	state.v.Y /= daysPerMillennium
	state.v.Z /= daysPerMillennium
	return state
}

// DeltaT returns the difference in seconds between Terrestrial Time (TT) and
// Universal Time (UT), using the default Delta T model by Espenak and Meeus.
func DeltaT(ut float64) float64 {
	return deltaTime(ut)
}

// synodicPeriod returns the average number of days it takes for a planet to
// return to the same relative longitude with respect to the Earth, as seen
// from the Sun.
func synodicPeriod(body Body) (float64, error) {
	if body == Earth {
		return 0.0, errors.New("The Earth does not have a synodic period as seen from itself")
	}
	if body == Moon {
		return 29.530588, nil
	}
	tp := PlanetOrbitalPeriod(body)
	if tp <= 0.0 {
		return 0.0, errors.New("Invalid body for synodic period")
	}
	return math.Abs(earthOrbitalPeriod / (earthOrbitalPeriod/tp - 1.0)), nil
}

// clampIndex returns the nearest valid index to index, in the range [0, n).
func clampIndex(index, n int) int {
	if index < 0 {
		return 0
	}
	if index >= n {
		return n - 1
	}
	return index
}

// correctLightTravel solves for the position of a body at the earlier time
// when the light we see now left the body. The calc function must return the
// body's position at a given time. The returned vector is the light-travel
// time corrected position.
func correctLightTravel(calc func(time AstroTime) (AstroVector, error), time AstroTime) (AstroVector, error) {
	ltime := time
	var pos AstroVector
	var err error
	for iter := 0; iter < 10; iter++ {
		pos, err = calc(ltime)
		if err != nil {
			return AstroVector{}, err
		}
		distance := pos.Length()
		if distance > SpeedOfLightAuPerDay {
			return AstroVector{}, errors.New("Cannot correct light travel time for distance greater than one light-day")
		}
		ltime2 := time.AddDays(-distance / SpeedOfLightAuPerDay)
		dt := math.Abs(ltime2.Tt - ltime.Tt)
		if dt < 1.0e-9 {
			return pos, nil // 86.4 microseconds
		}
		ltime = ltime2
	}
	return AstroVector{}, errors.New("Light travel time solver did not converge")
}

// spin rotates a vector by the given angle in degrees about the z-axis.
func spin(angle float64, pos [3]float64) [3]float64 {
	angr := RadiansFromDegrees(angle)
	cosang := math.Cos(angr)
	sinang := math.Sin(angr)
	return [3]float64{
		+cosang*pos[0] + sinang*pos[1],
		-sinang*pos[0] + cosang*pos[1],
		pos[2],
	}
}

// terra converts an observer's geographic position and sidereal time into
// geocentric equatorial position and velocity vectors (both in AU and AU/day).
func terra(observer Observer, st float64) (terseVector, terseVector) {
	const angVel = 7.2921150e-5
	phi := RadiansFromDegrees(observer.Latitude)
	sinphi := math.Sin(phi)
	cosphi := math.Cos(phi)
	c := 1.0 / math.Hypot(cosphi, sinphi*EarthFlattening)
	s := c * (EarthFlattening * EarthFlattening)
	htKm := observer.Height / 1000.0
	ach := EarthEquatorialRadiusKm*c + htKm
	ash := EarthEquatorialRadiusKm*s + htKm
	stlocl := RadiansFromDegrees(15.0*st + observer.Longitude)
	sinst := math.Sin(stlocl)
	cosst := math.Cos(stlocl)

	pos := terseVector{
		X: ach * cosphi * cosst / KmPerAu,
		Y: ach * cosphi * sinst / KmPerAu,
		Z: ash * sinphi / KmPerAu,
	}
	vel := terseVector{
		X: -(angVel * 86400.0 / KmPerAu) * ach * cosphi * sinst,
		Y: +(angVel * 86400.0 / KmPerAu) * ach * cosphi * cosst,
		Z: 0.0,
	}
	return pos, vel
}

// inverseTerra converts a geocentric equatorial position vector (in AU) at the
// given sidereal time into an observer's geographic position.
func inverseTerra(ovec terseVector, st float64) Observer {
	x := ovec.X * KmPerAu
	y := ovec.Y * KmPerAu
	z := ovec.Z * KmPerAu
	p := math.Hypot(x, y)

	var lonDeg, latDeg, heightKm float64
	if p < 1.0e-6 {
		// Special case: within 1 millimeter of a pole.
		lonDeg = 0.0
		if z > 0.0 {
			latDeg = +90.0
		} else {
			latDeg = -90.0
		}
		heightKm = math.Abs(z) - EarthPolarRadiusKm
	} else {
		stlocl := math.Atan2(y, x)
		lonDeg = DegreesFromRadians(stlocl) - (15.0 * st)
		for lonDeg <= -180.0 {
			lonDeg += 360.0
		}
		for lonDeg > +180.0 {
			lonDeg -= 360.0
		}
		// Numerically solve for exact latitude using Newton's Method.
		f := EarthFlattening * EarthFlattening
		lat := math.Atan2(z, p)
		for count := 0; ; count++ {
			if count > 10 {
				panic("inverseTerra did not converge")
			}
			// Calculate the error function W(lat) and its derivative D.
			c := math.Cos(lat)
			s := math.Sin(lat)
			factor := (f - 1) * EarthEquatorialRadiusKm
			c2 := c * c
			s2 := s * s
			radicand := c2 + f*s2
			denom := math.Sqrt(radicand)
			w := (factor*s*c)/denom - z*c + p*s
			if math.Abs(w) < 1.0e-8 {
				break
			}
			d := factor*((c2-s2)/denom-s2*c2*(f-1)/(factor*radicand)) + z*s + p*c
			lat -= w / d
		}
		latDeg = DegreesFromRadians(lat)
		// Solve for exact height in meters using the less risky denominator.
		adjust := EarthEquatorialRadiusKm / math.Sqrt(math.Cos(lat)*math.Cos(lat)+f*math.Sin(lat)*math.Sin(lat))
		if math.Abs(math.Sin(lat)) > math.Abs(math.Cos(lat)) {
			heightKm = z/math.Sin(lat) - f*adjust
		} else {
			heightKm = p/math.Cos(lat) - adjust
		}
	}
	return Observer{Latitude: latDeg, Longitude: lonDeg, Height: 1000.0 * heightKm}
}

// nutationRot creates a rotation matrix that adds or removes nutation from an
// equatorial vector of date.
func nutationRot(time *AstroTime, dir precessDirection) RotationMatrix {
	tilt := etilt(time)
	oblm := RadiansFromDegrees(tilt.Mobl)
	oblt := RadiansFromDegrees(tilt.Tobl)
	psi := tilt.Dpsi * asecToRad
	cobm := math.Cos(oblm)
	sobm := math.Sin(oblm)
	cobt := math.Cos(oblt)
	sobt := math.Sin(oblt)
	cpsi := math.Cos(psi)
	spsi := math.Sin(psi)

	xx := cpsi
	yx := -spsi * cobm
	zx := -spsi * sobm
	xy := spsi * cobt
	yy := cpsi*cobm*cobt + sobm*sobt
	zy := cpsi*sobm*cobt - cobm*sobt
	xz := spsi * sobt
	yz := cpsi*cobm*sobt - sobm*cobt
	zz := cpsi*sobm*sobt + cobm*cobt

	var rotation RotationMatrix
	if dir == from2000 {
		// convert J2000 to of-date
		rotation.Rot[0][0] = xx
		rotation.Rot[0][1] = xy
		rotation.Rot[0][2] = xz
		rotation.Rot[1][0] = yx
		rotation.Rot[1][1] = yy
		rotation.Rot[1][2] = yz
		rotation.Rot[2][0] = zx
		rotation.Rot[2][1] = zy
		rotation.Rot[2][2] = zz
	} else {
		// convert of-date to J2000
		rotation.Rot[0][0] = xx
		rotation.Rot[0][1] = yx
		rotation.Rot[0][2] = zx
		rotation.Rot[1][0] = xy
		rotation.Rot[1][1] = yy
		rotation.Rot[1][2] = zy
		rotation.Rot[2][0] = xz
		rotation.Rot[2][1] = yz
		rotation.Rot[2][2] = zz
	}
	return rotation
}

// nutationPosVel applies the nutation rotation to both the position and
// velocity components of a state vector.
func nutationPosVel(state StateVector, dir precessDirection) StateVector {
	return RotateState(nutationRot(&state.T, dir), state)
}

// precessionPosVel applies the precession rotation to both the position and
// velocity components of a state vector.
func precessionPosVel(state StateVector, dir precessDirection) StateVector {
	return RotateState(precessionRot(state.T, dir), state)
}

// Rotation_EQD_HOR returns a rotation matrix that converts equatorial
// of-date (EQD) coordinates to horizontal (HOR) coordinates for the given
// observer.
func Rotation_EQD_HOR(time AstroTime, observer Observer) RotationMatrix {
	uze := [3]float64{
		math.Cos(RadiansFromDegrees(observer.Latitude)) * math.Cos(RadiansFromDegrees(observer.Longitude)),
		math.Cos(RadiansFromDegrees(observer.Latitude)) * math.Sin(RadiansFromDegrees(observer.Longitude)),
		math.Sin(RadiansFromDegrees(observer.Latitude)),
	}
	une := [3]float64{
		-math.Sin(RadiansFromDegrees(observer.Latitude)) * math.Cos(RadiansFromDegrees(observer.Longitude)),
		-math.Sin(RadiansFromDegrees(observer.Latitude)) * math.Sin(RadiansFromDegrees(observer.Longitude)),
		math.Cos(RadiansFromDegrees(observer.Latitude)),
	}
	uwe := [3]float64{
		math.Sin(RadiansFromDegrees(observer.Longitude)),
		-math.Cos(RadiansFromDegrees(observer.Longitude)),
		0.0,
	}
	spinAngle := -15.0 * SiderealTime(&time)
	uz := spin(spinAngle, uze)
	un := spin(spinAngle, une)
	uw := spin(spinAngle, uwe)

	var rot RotationMatrix
	rot.Rot[0][0] = un[0]
	rot.Rot[1][0] = un[1]
	rot.Rot[2][0] = un[2]
	rot.Rot[0][1] = uw[0]
	rot.Rot[1][1] = uw[1]
	rot.Rot[2][1] = uw[2]
	rot.Rot[0][2] = uz[0]
	rot.Rot[1][2] = uz[1]
	rot.Rot[2][2] = uz[2]
	return rot
}

// Horizon calculates the azimuth and altitude of a celestial body as seen by
// an observer at the given time. The input right ascension and declination
// must be equator-of-date coordinates. The returned structure optionally
// corrects for atmospheric refraction.
func Horizon(time AstroTime, observer Observer, ra float64, dec float64, refraction Refraction) Topocentric {
	latrad := RadiansFromDegrees(observer.Latitude)
	lonrad := RadiansFromDegrees(observer.Longitude)
	decrad := RadiansFromDegrees(dec)
	rarad := ra * 15.0 * (math.Pi / 180.0)

	sinlat := math.Sin(latrad)
	coslat := math.Cos(latrad)
	sinlon := math.Sin(lonrad)
	coslon := math.Cos(lonrad)
	sindc := math.Sin(decrad)
	cosdc := math.Cos(decrad)
	sinra := math.Sin(rarad)
	cosra := math.Cos(rarad)

	// Uncorrected unit vectors for zenith, north, and west.
	uze := [3]float64{coslat * coslon, coslat * sinlon, sinlat}
	une := [3]float64{-sinlat * coslon, -sinlat * sinlon, coslat}
	uwe := [3]float64{sinlon, -coslon, 0.0}

	spinAngle := -15.0 * SiderealTime(&time)
	uz := spin(spinAngle, uze)
	un := spin(spinAngle, une)
	uw := spin(spinAngle, uwe)

	p := [3]float64{cosdc * cosra, cosdc * sinra, sindc}
	pn := p[0]*un[0] + p[1]*un[1] + p[2]*un[2]
	pw := p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2]
	pz := p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2]

	proj := math.Hypot(pn, pw)
	var az float64
	if proj > 0.0 {
		az = -DegreesFromRadians(math.Atan2(pw, pn))
		if az < 0.0 {
			az += 360.0
		}
	} else {
		az = 0.0
	}
	zd := DegreesFromRadians(math.Atan2(proj, pz))
	raOut := ra
	decOut := dec

	if refraction == NormalRefraction || refraction == JplHorizonsRefraction {
		zd0 := zd
		refr := RefractionAngle(refraction, 90.0-zd)
		zd -= refr
		if refr > 0.0 && zd > 3.0e-4 {
			sinzd := math.Sin(RadiansFromDegrees(zd))
			coszd := math.Cos(RadiansFromDegrees(zd))
			sinzd0 := math.Sin(RadiansFromDegrees(zd0))
			coszd0 := math.Cos(RadiansFromDegrees(zd0))
			var pr [3]float64
			for j := 0; j < 3; j++ {
				pr[j] = ((p[j] - coszd0*uz[j]) / sinzd0 * sinzd) + uz[j]*coszd
			}
			proj = math.Hypot(pr[0], pr[1])
			if proj > 0 {
				raOut = (DegreesFromRadians(math.Atan2(pr[1], pr[0]))) / 15.0
				if raOut < 0.0 {
					raOut += 24.0
				}
			} else {
				raOut = 0.0
			}
			decOut = DegreesFromRadians(math.Atan2(pr[2], proj))
		}
	}

	return Topocentric{
		Azimuth:  az,
		Altitude: 90.0 - zd,
		Ra:       raOut,
		Dec:      decOut,
	}
}

// HorizonFromVector converts a horizontal Cartesian vector into azimuth and
// altitude angles, optionally correcting for atmospheric refraction.
func HorizonFromVector(vector AstroVector, refraction Refraction) (Spherical, error) {
	sphere, err := SphereFromVector(vector)
	if err != nil {
		return Spherical{}, err
	}
	// Convert azimuth from counterclockwise-from-north to clockwise-from-north.
	sphere.Lon = toggleAzimuthDirection(sphere.Lon)
	sphere.Lat += RefractionAngle(refraction, sphere.Lat)
	return sphere, nil
}

// VectorFromHorizon converts apparent horizontal angular coordinates into a
// horizontal Cartesian vector, optionally removing atmospheric refraction.
func VectorFromHorizon(sphere Spherical, time AstroTime, refraction Refraction) AstroVector {
	// Convert azimuth from clockwise-from-north to counterclockwise-from-north.
	sphere.Lon = toggleAzimuthDirection(sphere.Lon)
	// Reverse any applied refraction.
	sphere.Lat += InverseRefractionAngle(refraction, sphere.Lat)
	return VectorFromSphere(sphere, time)
}

// InverseRefractionAngle returns the negative angular correction in degrees
// to subtract from an observed (refracted) altitude to obtain the true,
// unrefracted altitude.
func InverseRefractionAngle(refraction Refraction, bentAltitude float64) float64 {
	if bentAltitude < -90.0 || bentAltitude > +90.0 {
		return 0.0 // no attempt to correct an invalid altitude
	}
	altitude := bentAltitude - RefractionAngle(refraction, bentAltitude)
	for {
		// See how close we got.
		diff := (altitude + RefractionAngle(refraction, altitude)) - bentAltitude
		if math.Abs(diff) < 1.0e-14 {
			return altitude - bentAltitude
		}
		altitude -= diff
	}
}

// Search finds the time when a continuous function f crosses zero from
// negative to positive within the time range [t1, t2]. The dtToleranceSeconds
// value sets the maximum error in the returned time.
func Search(f func(t AstroTime) (float64, error), t1 AstroTime, t2 AstroTime, dtToleranceSeconds float64) (AstroTime, error) {
	dtDays := math.Abs(dtToleranceSeconds / 86400.0)
	f1, err := f(t1)
	if err != nil {
		return AstroTime{}, err
	}
	f2, err := f(t2)
	if err != nil {
		return AstroTime{}, err
	}
	const iterLimit = 20
	iter := 0
	calcFmid := true
	var fmid float64

	for {
		iter++
		if iter > iterLimit {
			return AstroTime{}, errors.New("Search did not converge")
		}
		dt := (t2.Tt - t1.Tt) / 2.0
		tmid := t1.AddDays(dt)
		if math.Abs(dt) < dtDays {
			return tmid, nil // close enough to the event
		}
		if calcFmid {
			fmid, err = f(tmid)
			if err != nil {
				return AstroTime{}, err
			}
		} else {
			calcFmid = true // we already have the correct value from the previous loop
		}

		// Quadratic interpolation: fit a parabola through (t1,f1), (tmid,fmid), (t2,f2).
		if qUt, qDfDt, ok := quadInterp(tmid.Ut, t2.Ut-tmid.Ut, f1, fmid, f2); ok {
			tq := TimeFromUniversalDays(qUt)
			fq, err := f(tq)
			if err != nil {
				return AstroTime{}, err
			}
			if qDfDt != 0.0 {
				dtGuess := math.Abs(fq / qDfDt)
				if dtGuess < dtDays {
					return tq, nil // the estimated time error is small enough
				}
				// Try guessing a tighter boundary with the interpolated root at the center.
				dtGuess *= 1.2
				if dtGuess < dt/10.0 {
					tleft := tq.AddDays(-dtGuess)
					tright := tq.AddDays(+dtGuess)
					if (tleft.Ut-t1.Ut)*(tleft.Ut-t2.Ut) < 0 {
						if (tright.Ut-t1.Ut)*(tright.Ut-t2.Ut) < 0 {
							fleft, err := f(tleft)
							if err != nil {
								return AstroTime{}, err
							}
							fright, err := f(tright)
							if err != nil {
								return AstroTime{}, err
							}
							if fleft < 0.0 && fright >= 0.0 {
								f1 = fleft
								f2 = fright
								t1 = tleft
								t2 = tright
								fmid = fq
								calcFmid = false // no need to recalculate fmid next loop
								continue
							}
						}
					}
				}
			}
		}

		// After quadratic interpolation, divide the region and pick the half that contains a root.
		if f1 < 0.0 && fmid >= 0.0 {
			t2 = tmid
			f2 = fmid
			continue
		}
		if fmid < 0.0 && f2 >= 0.0 {
			t1 = tmid
			f1 = fmid
			continue
		}
		return AstroTime{}, errors.New("Search failure: no ascending zero-crossing in the search window")
	}
}

// MoonEclipticLatitudeDegrees returns the Moon's ecliptic latitude in degrees
// at the given time.
func MoonEclipticLatitudeDegrees(time AstroTime) float64 {
	mc := makeMoonContext(time.Tt / 36525.0)
	moon := mc.calcMoon()
	return DegreesFromRadians(moon.geoEclipLat)
}

// Libration calculates the Moon's libration angles: the optical and physical
// librations of the Moon as seen from the Earth.
func Libration(time AstroTime) LibrationInfo {
	sinI := math.Sin(RadiansFromDegrees(1.543)) // sin of the Moon's axis inclination
	cosI := math.Cos(RadiansFromDegrees(1.543)) // cos of the Moon's axis inclination
	const moonMeanRadiusKm = 1737.4

	t := time.Tt / 36525.0
	t2 := t * t
	t3 := t2 * t
	t4 := t2 * t2

	mc := makeMoonContext(time.Tt / 36525.0)
	moon := mc.calcMoon()
	mlon := moon.geoEclipLon
	mlat := moon.geoEclipLat
	distKm := moon.distanceAu * KmPerAu

	var lib LibrationInfo
	lib.Mlon = DegreesFromRadians(mlon)
	lib.Mlat = DegreesFromRadians(mlat)
	lib.DistKm = distKm
	lib.DiamDeg = 2.0 * DegreesFromRadians(math.Atan(moonMeanRadiusKm/math.Sqrt(distKm*distKm-moonMeanRadiusKm*moonMeanRadiusKm)))

	// Moon's argument of latitude in radians.
	f := RadiansFromDegrees(normalizeLongitude(93.2720950 + 483202.0175233*t - 0.0036539*t2 - t3/3526000 + t4/863310000))
	// Moon's ascending node's mean longitude in radians.
	omega := RadiansFromDegrees(normalizeLongitude(125.0445479 - 1934.1362891*t + 0.0020754*t2 + t3/467441 - t4/60616000))
	// Sun's mean anomaly.
	m := RadiansFromDegrees(normalizeLongitude(357.5291092 + 35999.0502909*t - 0.0001536*t2 + t3/24490000))
	// Moon's mean anomaly.
	mdash := RadiansFromDegrees(normalizeLongitude(134.9633964 + 477198.8675055*t + 0.0087414*t2 + t3/69699 - t4/14712000))
	// Moon's mean elongation.
	d := RadiansFromDegrees(normalizeLongitude(297.8501921 + 445267.1114034*t - 0.0018819*t2 + t3/545868 - t4/113065000))
	// Eccentricity of the Earth's orbit.
	e := 1.0 - 0.002516*t - 0.0000074*t2

	// Optical librations.
	w := mlon - omega
	a := math.Atan2(math.Sin(w)*math.Cos(mlat)*cosI-math.Sin(mlat)*sinI, math.Cos(w)*math.Cos(mlat))
	ldash := longitudeOffset(DegreesFromRadians(a - f))
	bdash := math.Asin(-math.Sin(w)*math.Cos(mlat)*sinI - math.Sin(mlat)*cosI)

	// Physical librations.
	k1 := RadiansFromDegrees(119.75 + 131.849*t)
	k2 := RadiansFromDegrees(72.56 + 20.186*t)

	rho := (-0.02752*math.Cos(mdash) +
		-0.02245*math.Sin(f) +
		+0.00684*math.Cos(mdash-2*f) +
		-0.00293*math.Cos(2*f) +
		-0.00085*math.Cos(2*f-2*d) +
		-0.00054*math.Cos(mdash-2*d) +
		-0.00020*math.Sin(mdash+f) +
		-0.00020*math.Cos(mdash+2*f) +
		-0.00020*math.Cos(mdash-f) +
		+0.00014*math.Cos(mdash+2*f-2*d))

	sigma := (-0.02816*math.Sin(mdash) +
		+0.02244*math.Cos(f) +
		-0.00682*math.Sin(mdash-2*f) +
		-0.00279*math.Sin(2*f) +
		-0.00083*math.Sin(2*f-2*d) +
		+0.00069*math.Sin(mdash-2*d) +
		+0.00040*math.Cos(mdash+f) +
		-0.00025*math.Sin(2*mdash) +
		-0.00023*math.Sin(mdash+2*f) +
		+0.00020*math.Cos(mdash-f) +
		+0.00019*math.Sin(mdash-f) +
		+0.00013*math.Sin(mdash+2*f-2*d) +
		-0.00010*math.Cos(mdash-3*f))

	tau := (+0.02520*e*math.Sin(m) +
		+0.00473*math.Sin(2*mdash-2*f) +
		-0.00467*math.Sin(mdash) +
		+0.00396*math.Sin(k1) +
		+0.00276*math.Sin(2*mdash-2*d) +
		+0.00196*math.Sin(omega) +
		-0.00183*math.Cos(mdash-f) +
		+0.00115*math.Sin(mdash-2*d) +
		-0.00096*math.Sin(mdash-d) +
		+0.00046*math.Sin(2*f-2*d) +
		-0.00039*math.Sin(mdash-f) +
		-0.00032*math.Sin(mdash-m-d) +
		+0.00027*math.Sin(2*mdash-m-2*d) +
		+0.00023*math.Sin(k2) +
		-0.00014*math.Sin(2*d) +
		+0.00014*math.Cos(2*mdash-2*f) +
		-0.00012*math.Sin(mdash-2*f) +
		-0.00012*math.Sin(2*mdash) +
		+0.00011*math.Sin(2*mdash-2*m-2*d))

	ldash2 := -tau + (rho*math.Cos(a)+sigma*math.Sin(a))*math.Tan(bdash)
	bdashDeg := DegreesFromRadians(bdash)
	bdash2 := sigma*math.Cos(a) - rho*math.Sin(a)

	lib.Elon = ldash + ldash2
	lib.Elat = bdashDeg + bdash2
	return lib
}

// SolarEclipseObscuration calculates the fraction of the Sun's disc covered
// by the Moon, as seen from the given lunacentric observer position. hm is the
// heliocentric Moon vector and lo is the lunacentric observer vector.
func SolarEclipseObscuration(hm AstroVector, lo AstroVector) (float64, error) {
	// Find the heliocentric observer position.
	ho := AstroVector{
		X: hm.X + lo.X,
		Y: hm.Y + lo.Y,
		Z: hm.Z + lo.Z,
		T: lo.T,
	}
	sunAu := ho.Length()
	// Calculate the apparent angular radius of the Sun for the observer.
	sunRadius := math.Asin((SunRadiusKm / KmPerAu) / sunAu)
	// Calculate the apparent angular radius of the Moon for the observer.
	moonRadius := math.Asin((MoonPolarRadiusKm / KmPerAu) / lo.Length())
	// Calculate the apparent angular separation between the Sun and Moon centers.
	sepAngle := AngleBetween(lo, ho)
	obscuration := obscuration(sunRadius, moonRadius, RadiansFromDegrees(sepAngle))
	// HACK: in marginal cases, clamp obscuration to less than 1.0.
	if obscuration > 0.9999 {
		obscuration = 0.9999
	}
	return obscuration, nil
}

// jupiterMoonElem holds the orbital element series for one Galilean moon.
type jupiterMoonElem struct {
	mu   float64
	al   [2]float64
	a    vsopSeries
	l    vsopSeries
	z    vsopSeries
	zeta vsopSeries
}

// CalcEarth calculates the heliocentric equatorial position of the Earth at
// the given time using the VSOP87 model.
func calcEarth(time AstroTime) AstroVector {
	return calcVsop(&vsop[Earth], time)
}

// majorBodies holds the barycentric positions and velocities of the Sun and
// the four giant planets, used by the gravity simulator and Pluto calculations.
type majorBodies struct {
	Sun     bodyState
	Jupiter bodyState
	Saturn  bodyState
	Uranus  bodyState
	Neptune bodyState
}

// calcSolarSystemBarycenter calculates the heliocentric position of the Solar
// System Barycenter, considering the gravitational pulls of the four giant
// planets.
func calcSolarSystemBarycenter(time AstroTime) AstroVector {
	var ssb AstroVector
	ssb.T = time

	// Adjust the SSB for each giant planet's pull on the Sun.
	var bary majorBodies
	majorBodyBary(&bary, time.Tt)

	ssb.X = -bary.Sun.r.X
	ssb.Y = -bary.Sun.r.Y
	ssb.Z = -bary.Sun.r.Z
	return ssb
}

// adjustBarycenterPosVel adjusts the given SSB state by the effect of one
// major body, and returns that body's heliocentric position and velocity.
func adjustBarycenterPosVel(ssb *bodyState, tt float64, body Body, planetGM float64) bodyState {
	shift := planetGM / (planetGM + sunGM)
	planet := calcVsopPosVel(&vsop[body], tt)
	ssb.r.X += shift * planet.r.X
	ssb.r.Y += shift * planet.r.Y
	ssb.r.Z += shift * planet.r.Z
	ssb.v.X += shift * planet.v.X
	ssb.v.Y += shift * planet.v.Y
	ssb.v.Z += shift * planet.v.Z
	return planet
}

// majorBodyBary computes the barycentric state of the Sun and the four giant
// planets at the given Terrestrial Time.
func majorBodyBary(bary *majorBodies, tt float64) {
	bary.Sun.tt = tt
	bary.Sun.r = terseVector{}
	bary.Sun.v = terseVector{}

	bary.Jupiter = adjustBarycenterPosVel(&bary.Sun, tt, Jupiter, jupiterGM)
	bary.Saturn = adjustBarycenterPosVel(&bary.Sun, tt, Saturn, saturnGM)
	bary.Uranus = adjustBarycenterPosVel(&bary.Sun, tt, Uranus, uranusGM)
	bary.Neptune = adjustBarycenterPosVel(&bary.Sun, tt, Neptune, neptuneGM)

	// Convert planet states from heliocentric to barycentric.
	bary.Jupiter.r.X -= bary.Sun.r.X
	bary.Jupiter.r.Y -= bary.Sun.r.Y
	bary.Jupiter.r.Z -= bary.Sun.r.Z
	bary.Jupiter.v.X -= bary.Sun.v.X
	bary.Jupiter.v.Y -= bary.Sun.v.Y
	bary.Jupiter.v.Z -= bary.Sun.v.Z

	bary.Saturn.r.X -= bary.Sun.r.X
	bary.Saturn.r.Y -= bary.Sun.r.Y
	bary.Saturn.r.Z -= bary.Sun.r.Z
	bary.Saturn.v.X -= bary.Sun.v.X
	bary.Saturn.v.Y -= bary.Sun.v.Y
	bary.Saturn.v.Z -= bary.Sun.v.Z

	bary.Uranus.r.X -= bary.Sun.r.X
	bary.Uranus.r.Y -= bary.Sun.r.Y
	bary.Uranus.r.Z -= bary.Sun.r.Z
	bary.Uranus.v.X -= bary.Sun.v.X
	bary.Uranus.v.Y -= bary.Sun.v.Y
	bary.Uranus.v.Z -= bary.Sun.v.Z

	bary.Neptune.r.X -= bary.Sun.r.X
	bary.Neptune.r.Y -= bary.Sun.r.Y
	bary.Neptune.r.Z -= bary.Sun.r.Z
	bary.Neptune.v.X -= bary.Sun.v.X
	bary.Neptune.v.Y -= bary.Sun.v.Y
	bary.Neptune.v.Z -= bary.Sun.v.Z

	// Convert the heliocentric SSB to a barycentric Sun.
	bary.Sun.r.X = -bary.Sun.r.X
	bary.Sun.r.Y = -bary.Sun.r.Y
	bary.Sun.r.Z = -bary.Sun.r.Z
	bary.Sun.v.X = -bary.Sun.v.X
	bary.Sun.v.Y = -bary.Sun.v.Y
	bary.Sun.v.Z = -bary.Sun.v.Z
}

// smallBodyAcceleration calculates the acceleration experienced by a small
// body at the given barycentric position, due to the Sun and the four giant
// planets.
func smallBodyAcceleration(smallPos terseVector, bary *majorBodies) terseVector {
	var acc terseVector
	addAcceleration(&acc, smallPos, sunGM, bary.Sun.r)
	addAcceleration(&acc, smallPos, jupiterGM, bary.Jupiter.r)
	addAcceleration(&acc, smallPos, saturnGM, bary.Saturn.r)
	addAcceleration(&acc, smallPos, uranusGM, bary.Uranus.r)
	addAcceleration(&acc, smallPos, neptuneGM, bary.Neptune.r)
	return acc
}

// addAcceleration accumulates the gravitational pull of a major body on a
// small body at the given position.
func addAcceleration(acc *terseVector, smallPos terseVector, gm float64, majorPos terseVector) {
	dx := majorPos.X - smallPos.X
	dy := majorPos.Y - smallPos.Y
	dz := majorPos.Z - smallPos.Z
	r2 := dx*dx + dy*dy + dz*dz
	pull := gm / (r2 * math.Sqrt(r2))
	acc.X += dx * pull
	acc.Y += dy * pull
	acc.Z += dz * pull
}

// gravFromState converts a heliocentric position/velocity state into a
// barycentric gravity calculation including the acceleration at that state.
func gravFromState(bary *majorBodies, state *bodyState) bodyGravCalc {
	majorBodyBary(bary, state.tt)
	calc := bodyGravCalc{
		Tt: state.tt,
		R:  terseVector{X: state.r.X + bary.Sun.r.X, Y: state.r.Y + bary.Sun.r.Y, Z: state.r.Z + bary.Sun.r.Z},
		V:  terseVector{X: state.v.X + bary.Sun.v.X, Y: state.v.Y + bary.Sun.v.Y, Z: state.v.Z + bary.Sun.v.Z},
	}
	calc.A = smallBodyAcceleration(calc.R, bary)
	return calc
}

// gravSim advances a gravity calculation from time tt1 to tt2 using a
// predictor-corrector integration of the small body's motion.
func gravSim(bary2 *majorBodies, tt2 float64, calc1 *bodyGravCalc) bodyGravCalc {
	dt := tt2 - calc1.Tt

	// Calculate where the major bodies (Sun, Jupiter...Neptune) will be at the
	// next time step.
	majorBodyBary(bary2, tt2)

	// Estimate the position of the small body as if the current acceleration
	// applies across the whole time interval.
	approxPos := updatePosition(dt, calc1.R, calc1.V, calc1.A)

	// Calculate the acceleration at the approximate next location.
	acc := smallBodyAcceleration(approxPos, bary2)

	// Average the endpoint accelerations as the mean effective acceleration.
	acc = terseVector{
		X: (acc.X + calc1.A.X) / 2,
		Y: (acc.Y + calc1.A.Y) / 2,
		Z: (acc.Z + calc1.A.Z) / 2,
	}

	// Refine the estimates using the mean acceleration.
	var calc2 bodyGravCalc
	calc2.R = updatePosition(dt, calc1.R, calc1.V, acc)
	calc2.V = updateVelocity(dt, calc1.V, acc)
	calc2.A = smallBodyAcceleration(calc2.R, bary2)
	calc2.Tt = tt2
	return calc2
}

// plutoCache caches pre-calculated Pluto gravity segments so that repeated
// Pluto calculations at nearby times are fast.
var plutoCache []*bodySegment

// bodySegment is a series of gravity states interpolating Pluto's motion
// between two entries of the plutoStateTable.
type bodySegment struct {
	step [plutoNsteps]bodyGravCalc
}

// getSegment finds or calculates a cached gravity segment that contains the
// given Terrestrial Time. Returns -1 for the segment index when the time is
// outside the covered range.
func getSegment(cache []*bodySegment, tt float64) (int, *bodySegment, error) {
	if tt < plutoStateTable[0].tt || tt > plutoStateTable[plutoNumStates-1].tt {
		return -1, nil, nil
	}

	segIndex := clampIndex(int(math.Floor((tt-plutoStateTable[0].tt)/plutoTimeStep)), plutoNumStates-1)
	seg := cache[segIndex]
	if seg == nil {
		// Calculate the segment: pick the pair of bracketing body states.
		seg = &bodySegment{}
		cache[segIndex] = seg
		var bary majorBodies
		seg.step[0] = gravFromState(&bary, &plutoStateTable[segIndex])
		seg.step[plutoNsteps-1] = gravFromState(&bary, &plutoStateTable[segIndex+1])

		// Simulate forwards from the lower time bound.
		stepTT := seg.step[0].Tt
		for i := 1; i < plutoNsteps-1; i++ {
			stepTT += plutoDt
			seg.step[i] = gravSim(&bary, stepTT, &seg.step[i-1])
		}

		// Simulate backwards from the upper time bound.
		var reverse bodySegment
		reverse.step[plutoNsteps-1] = seg.step[plutoNsteps-1]
		stepTT = reverse.step[plutoNsteps-1].Tt
		for i := plutoNsteps - 2; i > 0; i-- {
			stepTT -= plutoDt
			reverse.step[i] = gravSim(&bary, stepTT, &reverse.step[i+1])
		}

		// Fade-mix the two series so that there are no discontinuities.
		for i := plutoNsteps - 2; i > 0; i-- {
			ramp := float64(i) / float64(plutoNsteps-1)
			seg.step[i].R = vecRamp(seg.step[i].R, reverse.step[i].R, ramp)
			seg.step[i].V = vecRamp(seg.step[i].V, reverse.step[i].V, ramp)
			seg.step[i].A = vecRamp(seg.step[i].A, reverse.step[i].A, ramp)
		}
	}
	return segIndex, seg, nil
}

// vecRamp linearly interpolates between two terse vectors.
func vecRamp(a, b terseVector, ramp float64) terseVector {
	return terseVector{
		X: (1-ramp)*a.X + ramp*b.X,
		Y: (1-ramp)*a.Y + ramp*b.Y,
		Z: (1-ramp)*a.Z + ramp*b.Z,
	}
}

// calcPlutoOneWay computes Pluto's barycentric state at targetTt by stepping
// forward (or backward) from a known initial state.
func calcPlutoOneWay(bary *majorBodies, initState *bodyState, targetTt float64, dt float64) bodyGravCalc {
	calc := gravFromState(bary, initState)
	n := int(math.Ceil((targetTt - calc.Tt) / dt))
	for i := 0; i < n; i++ {
		next := targetTt
		if i+1 != n {
			next = calc.Tt + dt
		}
		calc = gravSim(bary, next, &calc)
	}
	return calc
}

// calcPluto calculates the position and velocity of Pluto at the given time.
// When helio is true, the result is converted from barycentric to heliocentric.
func calcPluto(time AstroTime, helio bool) bodyState {
	if plutoCache == nil {
		plutoCache = make([]*bodySegment, plutoNumStates-1)
	}

	var bstate bodyState
	bstate.tt = time.Tt

	segIndex, seg, err := getSegment(plutoCache, time.Tt)
	if err != nil {
		return bstate
	}

	var bary majorBodies
	if segIndex < 0 {
		// The target time is outside the covered range; crawl from the nearest
		// endpoint.
		if time.Tt < plutoStateTable[0].tt {
			calc := calcPlutoOneWay(&bary, &plutoStateTable[0], time.Tt, -plutoDt)
			bstate.r = calc.R
			bstate.v = calc.V
		} else {
			calc := calcPlutoOneWay(&bary, &plutoStateTable[plutoNumStates-1], time.Tt, +plutoDt)
			bstate.r = calc.R
			bstate.v = calc.V
		}
	} else {
		left := clampIndex(int(math.Floor((time.Tt-seg.step[0].Tt)/plutoDt)), plutoNsteps-1)
		s1 := &seg.step[left]
		s2 := &seg.step[left+1]

		// Find the mean acceleration vector over the interval.
		acc := terseVector{
			X: (s1.A.X + s2.A.X) / 2,
			Y: (s1.A.Y + s2.A.Y) / 2,
			Z: (s1.A.Z + s2.A.Z) / 2,
		}

		// Extrapolate away from t1 in the positive time direction.
		ra := updatePosition(time.Tt-s1.Tt, s1.R, s1.V, acc)
		va := updateVelocity(time.Tt-s1.Tt, s1.V, acc)

		// Extrapolate away from t2 in the negative time direction.
		rb := updatePosition(time.Tt-s2.Tt, s2.R, s2.V, acc)
		vb := updateVelocity(time.Tt-s2.Tt, s2.V, acc)

		// Blend the two position estimates with a fade.
		ramp := (time.Tt - s1.Tt) / plutoDt
		bstate.r = vecRamp(ra, rb, ramp)
		bstate.v = vecRamp(va, vb, ramp)

		if helio {
			majorBodyBary(&bary, time.Tt)
		}
	}

	if helio {
		// Convert barycentric coordinates back to heliocentric coordinates.
		bstate.r.X -= bary.Sun.r.X
		bstate.r.Y -= bary.Sun.r.Y
		bstate.r.Z -= bary.Sun.r.Z
		bstate.v.X -= bary.Sun.v.X
		bstate.v.Y -= bary.Sun.v.Y
		bstate.v.Z -= bary.Sun.v.Z
	}
	return bstate
}

// EquatorDate selects the epoch of the Earth's equator for equatorial
// coordinates.
type EquatorDate int

const (
	EquatorOfDate EquatorDate = 0 // the Earth's equator at the given date/time
	EquatorJ2000  EquatorDate = 1 // the Earth's equator at the J2000 epoch
)

// geoPos converts an observer's geographic position to a geocentric equatorial
// position vector in J2000 coordinates.
func geoPos(time *AstroTime, observer Observer) [3]float64 {
	gast := SiderealTime(time)
	pos, _ := terra(observer, gast)
	vec := gyration(AstroVector{X: pos.X, Y: pos.Y, Z: pos.Z, T: *time}, into2000)
	return [3]float64{vec.X, vec.Y, vec.Z}
}

// ObserverVector calculates the geocentric equatorial position vector of an
// observer on the surface of the Earth at the given time.
func ObserverVector(time *AstroTime, observer Observer, equdate EquatorDate) (AstroVector, error) {
	gast := SiderealTime(time)
	pos, _ := terra(observer, gast)
	var vec AstroVector
	switch equdate {
	case EquatorOfDate:
		vec = AstroVector{X: pos.X, Y: pos.Y, Z: pos.Z, T: *time}
	case EquatorJ2000:
		vec = gyration(AstroVector{X: pos.X, Y: pos.Y, Z: pos.Z, T: *time}, into2000)
	default:
		return AstroVector{}, errors.New("Invalid equinox of date")
	}
	return vec, nil
}

// ObserverState calculates the geocentric equatorial position and velocity of
// an observer on the surface of the Earth at the given time.
func ObserverState(time *AstroTime, observer Observer, equdate EquatorDate) (StateVector, error) {
	gast := SiderealTime(time)
	pos, vel := terra(observer, gast)
	var state StateVector
	switch equdate {
	case EquatorOfDate:
		state = StateVector{X: pos.X, Y: pos.Y, Z: pos.Z, Vx: vel.X, Vy: vel.Y, Vz: vel.Z, T: *time}
	case EquatorJ2000:
		state = gyrationPosVel(StateVector{X: pos.X, Y: pos.Y, Z: pos.Z, Vx: vel.X, Vy: vel.Y, Vz: vel.Z, T: *time}, into2000)
	default:
		return StateVector{}, errors.New("Invalid equinox of date")
	}
	return state, nil
}

// VectorObserver is the inverse of ObserverVector: it converts a geocentric
// equatorial position vector into the corresponding geographic location.
func VectorObserver(vector *AstroVector, equdate EquatorDate) Observer {
	pos1 := [3]float64{vector.X, vector.Y, vector.Z}
	if equdate == EquatorJ2000 {
		mean := precession(AstroVector{X: pos1[0], Y: pos1[1], Z: pos1[2], T: vector.T}, from2000)
		eqd := nutation(mean, from2000)
		pos1 = [3]float64{eqd.X, eqd.Y, eqd.Z}
	}
	gast := SiderealTime(&vector.T)
	return inverseTerra(terseVector{X: pos1[0], Y: pos1[1], Z: pos1[2]}, gast)
}

// saturnMagnitude estimates the visual magnitude of Saturn, accounting for the
// tilt of its rings.
func saturnMagnitude(phase float64, helioDist float64, geoDist float64, gc AstroVector, time AstroTime) (float64, float64, error) {
	// Based on formulas by Paul Schlyter.
	// We must handle Saturn's rings as a major component of its visual magnitude.
	eclip := EquatorialToEcliptic(gc)
	ir := RadiansFromDegrees(28.06)                       // tilt of Saturn's rings to the ecliptic
	nr := RadiansFromDegrees(169.51 + 3.82e-5*time.Tt)   // ascending node of Saturn's rings

	lat := RadiansFromDegrees(eclip.Elat)
	lon := RadiansFromDegrees(eclip.Elon)
	tilt := math.Asin(math.Sin(lat)*math.Cos(ir) - math.Cos(lat)*math.Sin(ir)*math.Sin(lon-nr))
	sinTilt := math.Sin(math.Abs(tilt))

	mag := -9.0 + 0.044*phase
	mag += sinTilt * (-2.6 + 1.2*sinTilt)
	mag += 5.0 * math.Log10(helioDist*geoDist)
	return mag, DegreesFromRadians(tilt), nil
}

// moonNodeSearchFunc evaluates the Moon's ecliptic latitude, signed by the
// node kind, used for finding where the Moon's orbit crosses the ecliptic.
func moonNodeSearchFunc(kind NodeEventKind, time AstroTime) (float64, error) {
	eclip := EclipticGeoMoon(time)
	return eclip.Lat * float64(kind), nil
}

// SearchMoonNode searches for the first ascending or descending node of the
// Moon after the given start time. An ascending node is where the Moon's orbit
// crosses the ecliptic from south to north; a descending node is the reverse.
func SearchMoonNode(startTime AstroTime) (NodeEventInfo, error) {
	const stepDays = 10.0
	time1 := startTime
	eclip1 := EclipticGeoMoon(time1)
	for {
		time2 := time1.AddDays(stepDays)
		eclip2 := EclipticGeoMoon(time2)
		if eclip1.Lat*eclip2.Lat <= 0.0 {
			// There is a node somewhere inside this closed time interval.
			kind := DescendingNode
			if eclip2.Lat > eclip1.Lat {
				kind = AscendingNode
			}
			root, err := Search(func(t AstroTime) (float64, error) {
				return moonNodeSearchFunc(kind, t)
			}, time1, time2, 1.0)
			if err != nil {
				return NodeEventInfo{}, err
			}
			return NodeEventInfo{Time: root, Kind: kind}, nil
		}
		time1 = time2
		eclip1 = eclip2
	}
}

// NextMoonNode finds the next lunar node event after the given previous node.
func NextMoonNode(prevNode NodeEventInfo) (NodeEventInfo, error) {
	time := prevNode.Time.AddDays(10.0)
	// Search for the next node of the opposite kind.
	opposite := DescendingNode
	if prevNode.Kind == DescendingNode {
		opposite = AscendingNode
	}
	node, err := SearchMoonNode(time)
	if err != nil {
		return NodeEventInfo{}, err
	}
	return node, nil
}

// Rotation_EQJ_HOR returns a rotation matrix that converts J2000 mean equator
// (EQJ) coordinates to horizontal (HOR) coordinates.
func Rotation_EQJ_HOR(time AstroTime, observer Observer) RotationMatrix {
	return InverseRotation(Rotation_HOR_EQJ(time, observer))
}

// Rotation_ECL_EQD returns a rotation matrix that converts J2000 mean ecliptic
// (ECL) coordinates to equatorial of-date (EQD) coordinates.
func Rotation_ECL_EQD(time AstroTime) RotationMatrix {
	return InverseRotation(Rotation_EQD_ECL(time))
}

// Rotation_ECL_HOR returns a rotation matrix that converts J2000 mean ecliptic
// (ECL) coordinates to horizontal (HOR) coordinates.
func Rotation_ECL_HOR(time AstroTime, observer Observer) RotationMatrix {
	return CombineRotation(Rotation_ECL_EQD(time), Rotation_EQD_HOR(time, observer))
}

// calcJupiterMoon calculates the jovicentric position and velocity of one of
// Jupiter's four Galilean moons (Io, Europa, Ganymede, Callisto) at the given
// time, expressed in the J2000 equatorial (EQJ) frame.
func calcJupiterMoon(time AstroTime, mindex int) StateVector {
	m := &jupiterMoonModel[mindex]
	t := time.Tt + 18262.5 // time since 1950-01-01T00:00:00Z

	// Calculate 6 orbital elements at the given time t.
	var elem [6]float64

	elem[0] = 0.0
	for k := 0; k < m.a.nterms; k++ {
		arg := m.a.term[k].pha + (t * m.a.term[k].fre)
		elem[0] += m.a.term[k].amp * math.Cos(arg)
	}

	elem[1] = m.al[0] + (t * m.al[1])
	for k := 0; k < m.l.nterms; k++ {
		arg := m.l.term[k].pha + (t * m.l.term[k].fre)
		elem[1] += m.l.term[k].amp * math.Sin(arg)
	}
	elem[1] = math.Mod(elem[1], 2.0*math.Pi)
	if elem[1] < 0.0 {
		elem[1] += 2.0 * math.Pi
	}

	for k := 0; k < m.z.nterms; k++ {
		arg := m.z.term[k].pha + (t * m.z.term[k].fre)
		elem[2] += m.z.term[k].amp * math.Cos(arg)
		elem[3] += m.z.term[k].amp * math.Sin(arg)
	}

	for k := 0; k < m.zeta.nterms; k++ {
		arg := m.zeta.term[k].pha + (t * m.zeta.term[k].fre)
		elem[4] += m.zeta.term[k].amp * math.Cos(arg)
		elem[5] += m.zeta.term[k].amp * math.Sin(arg)
	}

	// Convert the orbital elements into a position vector in the Jupiter
	// equatorial system (JUP), then re-orient to the J2000 equatorial frame.
	state := jupiterMoonElem2pv(time, m.mu, elem)
	return RotateState(rotationJupEqj, state)
}

// JupiterMoons calculates the jovicentric positions and velocities of
// Jupiter's four largest moons (Io, Europa, Ganymede, and Callisto) at the
// given time, expressed in the J2000 equatorial (EQJ) frame.
func JupiterMoons(time AstroTime) JupiterMoonsInfo {
	return JupiterMoonsInfo{
		Io:       calcJupiterMoon(time, 0),
		Europa:   calcJupiterMoon(time, 1),
		Ganymede: calcJupiterMoon(time, 2),
		Callisto: calcJupiterMoon(time, 3),
	}
}

// nutation applies the nutation rotation matrix to a vector.
func nutation(pos AstroVector, dir precessDirection) AstroVector {
	return RotateVector(nutationRot(&pos.T, dir), pos)
}

// Rotation_EQJ_EQD returns a rotation matrix that converts J2000 mean equator
// (EQJ) coordinates to equatorial of-date (EQD) coordinates.
func Rotation_EQJ_EQD(time AstroTime) RotationMatrix {
	prec := precessionRot(time, from2000)
	nut := nutationRot(&time, from2000)
	return CombineRotation(prec, nut)
}

// Rotation_EQD_EQJ returns a rotation matrix that converts equatorial of-date
// (EQD) coordinates to J2000 mean equator (EQJ) coordinates.
func Rotation_EQD_EQJ(time AstroTime) RotationMatrix {
	nut := nutationRot(&time, into2000)
	prec := precessionRot(time, into2000)
	return CombineRotation(nut, prec)
}

// Rotation_HOR_EQD returns a rotation matrix that converts horizontal (HOR)
// coordinates to equatorial of-date (EQD) coordinates.
func Rotation_HOR_EQD(time AstroTime, observer Observer) RotationMatrix {
	return InverseRotation(Rotation_EQD_HOR(time, observer))
}

// Rotation_ECT_EQD returns a rotation matrix that converts true ecliptic of
// date (ECT) coordinates to equatorial of-date (EQD) coordinates.
func Rotation_ECT_EQD(time AstroTime) RotationMatrix {
	tilt := etilt(&time)
	tobl := RadiansFromDegrees(tilt.Tobl)
	cosTobl := math.Cos(tobl)
	sinTobl := math.Sin(tobl)

	var m RotationMatrix
	// EQD.x = ECT.x
	m.Rot[0][0] = 1.0
	m.Rot[1][0] = 0.0
	m.Rot[2][0] = 0.0
	// EQD.y = +cos*ECT.y - sin*ECT.z
	m.Rot[0][1] = 0.0
	m.Rot[1][1] = +cosTobl
	m.Rot[2][1] = -sinTobl
	// EQD.z = +sin*ECT.y + cos*ECT.z
	m.Rot[0][2] = 0.0
	m.Rot[1][2] = +sinTobl
	m.Rot[2][2] = +cosTobl
	return m
}

// Rotation_EQD_ECT returns a rotation matrix that converts equatorial of-date
// (EQD) coordinates to true ecliptic of date (ECT) coordinates.
func Rotation_EQD_ECT(time AstroTime) RotationMatrix {
	return InverseRotation(Rotation_ECT_EQD(time))
}

// ApsisKind indicates whether an apsis event is a closest approach
// (pericenter) or a farthest distance (apocenter).
type ApsisKind int

const (
	ApsisKindPericenter ApsisKind = iota // the body is at its closest approach
	ApsisKindApocenter                   // the body is at its farthest distance
	ApsisKindInvalid                     // undefined or invalid apsis
)

// ApsisInfo describes a closest-approach or farthest-distance event of an
// orbiting body, such as the Moon's perigee and apogee.
type ApsisInfo struct {
	Time    AstroTime // the date and time of the apsis
	Kind    ApsisKind // whether this is a pericenter or apocenter event
	DistAu  float64   // the distance between the centers of the bodies in AU
	DistKm  float64   // the distance between the centers of the bodies in kilometers
}

// moonDistance returns the geocentric distance to the Moon in AU.
func moonDistance(t AstroTime) float64 {
	mc := makeMoonContext(t.Tt / 36525.0)
	moon := mc.calcMoon()
	return moon.distanceAu
}

// moonDistanceSlope evaluates the rate of change of the Moon's distance.
func moonDistanceSlope(time AstroTime, direction float64) float64 {
	dt := 0.001
	t1 := time.AddDays(-dt / 2.0)
	t2 := time.AddDays(+dt / 2.0)
	return direction * (moonDistance(t2) - moonDistance(t1)) / dt
}

// SearchLunarApsis finds the next perigee (closest approach) or apogee
// (farthest distance) of the Moon after the given start time.
func SearchLunarApsis(startTime AstroTime) (ApsisInfo, error) {
	const increment = 5.0 // number of days to skip in each iteration
	t1 := startTime
	m1 := moonDistanceSlope(t1, +1.0)

	for iter := 0.0; iter < 2.0*29.530588; iter += increment {
		t2 := t1.AddDays(increment)
		m2 := moonDistanceSlope(t2, +1.0)

		if m1*m2 <= 0.0 {
			// There is a change of slope polarity within [t1, t2], so this
			// range contains an apsis. Figure out whether it is perigee or apogee.
			var kind ApsisKind
			var root AstroTime
			var err error
			if m1 < 0.0 || m2 > 0.0 {
				// Minimum-distance event: perigee.
				kind = ApsisKindPericenter
				root, err = Search(func(t AstroTime) (float64, error) {
					return moonDistanceSlope(t, +1.0), nil
				}, t1, t2, 1.0)
			} else if m1 > 0.0 || m2 < 0.0 {
				// Maximum-distance event: apogee.
				kind = ApsisKindApocenter
				root, err = Search(func(t AstroTime) (float64, error) {
					return moonDistanceSlope(t, -1.0), nil
				}, t1, t2, 1.0)
			} else {
				return ApsisInfo{}, errors.New("Both Moon distance slopes are zero")
			}
			if err != nil {
				return ApsisInfo{}, err
			}
			distAu := moonDistance(root)
			return ApsisInfo{Time: root, Kind: kind, DistAu: distAu, DistKm: distAu * KmPerAu}, nil
		}

		// Not yet found a slope polarity change; keep searching.
		t1 = t2
		m1 = m2
	}
	return ApsisInfo{}, errors.New("Could not find a lunar apsis within two synodic months")
}

// NextLunarApsis finds the next lunar perigee or apogee event after the given
// apsis event. Given an apogee event, it finds the next perigee, and vice versa.
func NextLunarApsis(apsis ApsisInfo) (ApsisInfo, error) {
	if apsis.Kind != ApsisKindApocenter && apsis.Kind != ApsisKindPericenter {
		return ApsisInfo{}, errors.New("Invalid apsis kind")
	}
	const skip = 11.0 // days to skip to start looking for the next apsis event
	next, err := SearchLunarApsis(apsis.Time.AddDays(skip))
	if err != nil {
		return ApsisInfo{}, err
	}
	// Verify that we found the opposite apsis from the previous one.
	if next.Kind+apsis.Kind != 1 {
		return ApsisInfo{}, errors.New("Found same apsis kind instead of opposite")
	}
	return next, nil
}

// gyration combines nutation and precession into a single operation.
func gyration(pos AstroVector, dir precessDirection) AstroVector {
	if dir == into2000 {
		return precession(nutation(pos, dir), dir)
	}
	return nutation(precession(pos, dir), dir)
}

// gyrationPosVel combines nutation and precession for a state vector.
func gyrationPosVel(state StateVector, dir precessDirection) StateVector {
	if dir == into2000 {
		return precessionPosVel(nutationPosVel(state, dir), dir)
	}
	return nutationPosVel(precessionPosVel(state, dir), dir)
}

// EclipticGeoMoon calculates the geocentric ecliptic coordinates of the Moon
// at the given time, using the true equinox of date.
func EclipticGeoMoon(time AstroTime) Spherical {
	mc := makeMoonContext(time.Tt / 36525.0)
	moon := mc.calcMoon()

	// CalcMoon produces ecliptic coordinates in mean equinox of date (ECM).
	distCosLat := moon.distanceAu * math.Cos(moon.geoEclipLat)
	ecm := [3]float64{
		distCosLat * math.Cos(moon.geoEclipLon),
		distCosLat * math.Sin(moon.geoEclipLon),
		moon.distanceAu * math.Sin(moon.geoEclipLat),
	}

	// Obtain true and mean obliquity angles for the given time.
	// This also pre-calculates and caches the nutation in `time`.
	tilt := etilt(&time)

	// Convert ecliptic coordinates (ECM) to equatorial coordinates (EQM),
	// both in mean equinox of date.
	eqm := eclOblToEquVec(AstroVector{X: ecm[0], Y: ecm[1], Z: ecm[2], T: time}, RadiansFromDegrees(tilt.Mobl))

	// Add nutation to convert to true equatorial coordinates of date (EQD).
	eqd := nutation(eqm, from2000)

	// Convert back to ecliptic using the true equinox of date (ECT).
	eclip := rotateEquatorialToEcliptic([3]float64{eqd.X, eqd.Y, eqd.Z}, RadiansFromDegrees(tilt.Tobl), time)
	return Spherical{Lat: eclip.Elat, Lon: eclip.Elon, Dist: moon.distanceAu}
}

// EquatorialToEcliptic converts a J2000 mean equator (EQJ) vector into true
// ecliptic of date (ECT) coordinates, returning both Cartesian and spherical
// results.
func EquatorialToEcliptic(eqj AstroVector) Ecliptic {
	// Calculate nutation and obliquity for this time.
	tilt := etilt(&eqj.T)

	// Convert mean J2000 equator (EQJ) to true equator of date (EQD).
	mean := precession(eqj, from2000)
	eqd := nutation(mean, from2000)

	// Rotate from EQD to true ecliptic of date (ECT).
	return rotateEquatorialToEcliptic([3]float64{eqd.X, eqd.Y, eqd.Z}, RadiansFromDegrees(tilt.Tobl), eqj.T)
}

// GeoMoonState calculates the geocentric position and velocity of the Moon at
// the given time, using numerical differentiation of GeoMoon.
func GeoMoonState(time AstroTime) StateVector {
	const dt = 1.0e-5 // 0.864 seconds
	t1 := time.AddDays(-dt)
	t2 := time.AddDays(+dt)
	r1 := GeoMoon(t1)
	r2 := GeoMoon(t2)
	return StateVector{
		X:  (r1.X + r2.X) / 2,
		Y:  (r1.Y + r2.Y) / 2,
		Z:  (r1.Z + r2.Z) / 2,
		Vx: (r2.X - r1.X) / (2 * dt),
		Vy: (r2.Y - r1.Y) / (2 * dt),
		Vz: (r2.Z - r1.Z) / (2 * dt),
		T:  time,
	}
}

// GeoEmbState calculates the geocentric position and velocity of the
// Earth/Moon barycenter at the given time.
func GeoEmbState(time AstroTime) StateVector {
	s := GeoMoonState(time)
	const d = 1.0 + earthMoonMassRatio
	return StateVector{
		X:  s.X / d,
		Y:  s.Y / d,
		Z:  s.Z / d,
		Vx: s.Vx / d,
		Vy: s.Vy / d,
		Vz: s.Vz / d,
		T:  time,
	}
}

// EarthRotationAxis calculates the orientation of the Earth's rotation axis at
// the given time, returning the north pole vector and the rotation angle of
// the prime meridian.
func EarthRotationAxis(time *AstroTime) AxisInfo {
	// Start with a north pole vector in equator-of-date coordinates: (0,0,1),
	// then convert to J2000 coordinates.
	north := gyration(AstroVector{X: 0, Y: 0, Z: 1, T: *time}, into2000)

	equ, err := EquatorFromVector(north)
	if err != nil {
		panic(err)
	}
	// Use a modified version of the era() function that does not trim to
	// 0..360 degrees, corrected to give the correct angle at J2000.
	spin := 190.41375788700253 + (360.9856122880876 * time.Ut)
	return AxisInfo{
		Ra:    equ.Ra,
		Dec:   equ.Dec,
		Spin:  spin,
		North: north,
	}
}

// RotationAxis calculates the orientation of a body's rotation axis and the
// rotation angle of its prime meridian at the given time, using the IAU
// Working Group on Cartographics and Rotational Elements 2015 report formulas.
func RotationAxis(body Body, time *AstroTime) (AxisInfo, error) {
	var ra, dec, w float64
	const d = 36525.0
	tt := time.Tt
	T := tt / d

	switch body {
	case Sun:
		ra = 286.13
		dec = 63.87
		w = 84.176 + (14.1844 * tt)
	case Mercury:
		ra = 281.0103 - (0.0328 * T)
		dec = 61.4155 - (0.0049 * T)
		w = 329.5988 +
			(6.1385108 * tt) +
			(0.01067257 * math.Sin(RadiansFromDegrees(174.7910857+4.092335*tt))) -
			(0.00112309 * math.Sin(RadiansFromDegrees(349.5821714+8.184670*tt))) -
			(0.00011040 * math.Sin(RadiansFromDegrees(164.3732571+12.277005*tt))) -
			(0.00002539 * math.Sin(RadiansFromDegrees(339.1643429+16.369340*tt))) -
			(0.00000571 * math.Sin(RadiansFromDegrees(153.9554286+20.461675*tt)))
	case Venus:
		ra = 272.76
		dec = 67.16
		w = 160.20 - (1.4813688 * tt)
	case Earth:
		return EarthRotationAxis(time), nil
	case Moon:
		e1 := RadiansFromDegrees(125.045 - 0.0529921*tt)
		e2 := RadiansFromDegrees(250.089 - 0.1059842*tt)
		e3 := RadiansFromDegrees(260.008 + 13.0120009*tt)
		e4 := RadiansFromDegrees(176.625 + 13.3407154*tt)
		e5 := RadiansFromDegrees(357.529 + 0.9856003*tt)
		e6 := RadiansFromDegrees(311.589 + 26.4057084*tt)
		e7 := RadiansFromDegrees(134.963 + 13.0649930*tt)
		e8 := RadiansFromDegrees(276.617 + 0.3287146*tt)
		e9 := RadiansFromDegrees(34.226 + 1.7484877*tt)
		e10 := RadiansFromDegrees(15.134 - 0.1589763*tt)
		e11 := RadiansFromDegrees(119.743 + 0.0036096*tt)
		e12 := RadiansFromDegrees(239.961 + 0.1643573*tt)
		e13 := RadiansFromDegrees(25.053 + 12.9590088*tt)

		ra = 269.9949 + 0.0031*T -
			3.8787*math.Sin(e1) -
			0.1204*math.Sin(e2) +
			0.0700*math.Sin(e3) -
			0.0172*math.Sin(e4) +
			0.0072*math.Sin(e6) -
			0.0052*math.Sin(e10) +
			0.0043*math.Sin(e13)

		dec = 66.5392 + 0.0130*T +
			1.5419*math.Cos(e1) +
			0.0239*math.Cos(e2) -
			0.0278*math.Cos(e3) +
			0.0068*math.Cos(e4) -
			0.0029*math.Cos(e6) +
			0.0009*math.Cos(e7) +
			0.0008*math.Cos(e10) -
			0.0009*math.Cos(e13)

		w = 38.3213 + (13.17635815-1.4e-12*tt)*tt +
			3.5610*math.Sin(e1) +
			0.1208*math.Sin(e2) -
			0.0642*math.Sin(e3) +
			0.0158*math.Sin(e4) +
			0.0252*math.Sin(e5) -
			0.0066*math.Sin(e6) -
			0.0047*math.Sin(e7) -
			0.0046*math.Sin(e8) +
			0.0028*math.Sin(e9) +
			0.0052*math.Sin(e10) +
			0.0040*math.Sin(e11) +
			0.0019*math.Sin(e12) -
			0.0044*math.Sin(e13)
	case Mars:
		ra = 317.269202 - 0.10927547*T +
			0.000068*math.Sin(RadiansFromDegrees(198.991226+19139.4819985*T)) +
			0.000238*math.Sin(RadiansFromDegrees(226.292679+38280.8511281*T)) +
			0.000052*math.Sin(RadiansFromDegrees(249.663391+57420.7251593*T)) +
			0.000009*math.Sin(RadiansFromDegrees(266.183510+76560.6367950*T)) +
			0.419057*math.Sin(RadiansFromDegrees(79.398797+0.5042615*T))
		dec = 54.432516 - 0.05827105*T +
			0.000051*math.Cos(RadiansFromDegrees(122.433576+19139.9407476*T)) +
			0.000141*math.Cos(RadiansFromDegrees(43.058401+38280.8753272*T)) +
			0.000031*math.Cos(RadiansFromDegrees(57.663379+57420.7517205*T)) +
			0.000005*math.Cos(RadiansFromDegrees(79.476401+76560.6495004*T)) +
			1.591274*math.Cos(RadiansFromDegrees(166.325722+0.5042615*T))
		w = 176.049863 + 350.891982443297*tt +
			0.000145*math.Sin(RadiansFromDegrees(129.071773+19140.0328244*T)) +
			0.000157*math.Sin(RadiansFromDegrees(36.352167+38281.0473591*T)) +
			0.000040*math.Sin(RadiansFromDegrees(56.668646+57420.9295360*T)) +
			0.000001*math.Sin(RadiansFromDegrees(67.364003+76560.2552215*T)) +
			0.000001*math.Sin(RadiansFromDegrees(104.792680+95700.4387578*T)) +
			0.584542*math.Sin(RadiansFromDegrees(95.391654+0.5042615*T))
	case Jupiter:
		ja := RadiansFromDegrees(99.360714 + 4850.4046*T)
		jb := RadiansFromDegrees(175.895369 + 1191.9605*T)
		jc := RadiansFromDegrees(300.323162 + 262.5475*T)
		jd := RadiansFromDegrees(114.012305 + 6070.2476*T)
		je := RadiansFromDegrees(49.511251 + 64.3000*T)
		ra = 268.056595 - 0.006499*T +
			0.000117*math.Sin(ja) +
			0.000938*math.Sin(jb) +
			0.001432*math.Sin(jc) +
			0.000030*math.Sin(jd) +
			0.002150*math.Sin(je)
		dec = 64.495303 + 0.002413*T +
			0.000050*math.Cos(ja) +
			0.000404*math.Cos(jb) +
			0.000617*math.Cos(jc) -
			0.000013*math.Cos(jd) +
			0.000926*math.Cos(je)
		w = 284.95 + 870.536*tt
	case Saturn:
		ra = 40.589 - 0.036*T
		dec = 83.537 - 0.004*T
		w = 38.90 + 810.7939024*tt
	case Uranus:
		ra = 257.311
		dec = -15.175
		w = 203.81 - 501.1600928*tt
	case Neptune:
		n := RadiansFromDegrees(357.85 + 52.316*T)
		ra = 299.36 + 0.70*math.Sin(n)
		dec = 43.46 - 0.51*math.Cos(n)
		w = 249.978 + 541.1397757*tt - 0.48*math.Sin(n)
	case Pluto:
		ra = 132.993
		dec = -6.163
		w = 302.695 + 56.3625225*tt
	default:
		return AxisInfo{}, errors.New("Invalid body for rotation axis")
	}

	radlat := RadiansFromDegrees(dec)
	radlon := RadiansFromDegrees(ra)
	rcoslat := math.Cos(radlat)
	axis := AxisInfo{
		Ra:   ra / 15.0, // convert degrees to sidereal hours
		Dec:  dec,
		Spin: w,
		North: AstroVector{
			X: rcoslat * math.Cos(radlon),
			Y: rcoslat * math.Sin(radlon),
			Z: math.Sin(radlat),
			T: *time,
		},
	}
	return axis, nil
}

// ConstellationInfo identifies the constellation that contains a point in the
// sky, based on the International Astronomical Union's 1922 constellation
// boundaries in the B1875 epoch.
type ConstellationInfo struct {
	Symbol  string  // the 3-letter abbreviation of the constellation
	Name    string  // the full name of the constellation
	Ra1875  float64 // right ascension in B1875 coordinates (sidereal hours)
	Dec1875 float64 // declination in B1875 coordinates (degrees)
}

// constellationEpochTime is the time of the B1875 epoch, used for converting
// coordinates into the constellation boundary frame.
var constellationEpochTime = TimeFromUniversalDays(-45655.74141261017)

// Constellation finds the constellation that contains the given J2000
// equatorial coordinates (right ascension in sidereal hours, declination in
// degrees).
func Constellation(ra float64, dec float64) (ConstellationInfo, error) {
	if dec < -90.0 || dec > +90.0 {
		return ConstellationInfo{}, errors.New("Invalid declination for constellation")
	}
	// Allow right ascension to wrap around; clamp to [0, 24) sidereal hours.
	ra = math.Mod(ra, 24.0)
	if ra < 0.0 {
		ra += 24.0
	}

	// Convert coordinates from J2000 to the B1875 epoch.
	vec2000 := VectorFromSphere(Spherical{Lat: dec, Lon: ra * 15.0, Dist: 1.0}, constellationEpochTime)
	rot := Rotation_EQJ_EQD(constellationEpochTime)
	vec1875 := RotateVector(rot, vec2000)
	b1875, err := EquatorFromVector(vec1875)
	if err != nil {
		return ConstellationInfo{}, err
	}

	// Convert DEC from degrees and RA from hours to the compact angle units
	// used in the constellation boundary table.
	xRa := (24.0 * 15.0) * b1875.Ra
	xDec := 24.0 * b1875.Dec

	// Search for the constellation using the B1875 coordinates.
	index := -1
	for i := range constelBounds {
		b := &constelBounds[i]
		if b.decLo <= xDec && b.raHi > xRa && b.raLo <= xRa {
			index = b.index
			break
		}
	}
	if index < 0 || index >= len(constelNames) {
		return ConstellationInfo{}, errors.New("Constellation not found")
	}
	return ConstellationInfo{
		Symbol:  constelNames[index].symbol,
		Name:    constelNames[index].name,
		Ra1875:  b1875.Ra,
		Dec1875: b1875.Dec,
	}, nil
}

// findSeasonChange searches for the moment when the Sun reaches the given
// ecliptic longitude, starting from a date in the given year.
func findSeasonChange(targetLon float64, year, month, day int) (AstroTime, error) {
	startTime := TimeFromCalendar(year, month, day, 0, 0, 0.0)
	return SearchSunLongitude(targetLon, startTime, 20.0)
}

// SunPosition calculates the geocentric ecliptic coordinates of the Sun at the
// given time, using the true equinox of date.
func SunPosition(time AstroTime) Ecliptic {
	// Correct for light travel time from the Sun, otherwise season
	// calculations (equinox, solstice) will all be early by about 8 minutes.
	adjustedTime := time.AddDays(-1.0 / SpeedOfLightAuPerDay)

	// Find the heliocentric equatorial coordinates of Earth at the J2000 epoch,
	// then negate to get the geocentric position of the Sun.
	earth2000 := calcEarth(adjustedTime)
	sun2000 := AstroVector{X: -earth2000.X, Y: -earth2000.Y, Z: -earth2000.Z, T: adjustedTime}

	// Convert from J2000 mean equator (EQJ) to mean equator of date (EQM).
	tilt := etilt(&adjustedTime)
	meanEquinox := precession(sun2000, from2000)

	// Add nutation to convert to true equator of date (EQD).
	eqd := nutation(meanEquinox, from2000)

	// Convert EQD to true ecliptic of date (ECT).
	return rotateEquatorialToEcliptic([3]float64{eqd.X, eqd.Y, eqd.Z}, RadiansFromDegrees(tilt.Tobl), time)
}

// sunOffset returns the difference between the Sun's ecliptic longitude and a
// target longitude, reduced into the range (-180, +180].
func sunOffset(targetLon float64, time AstroTime) (float64, error) {
	ecl := SunPosition(time)
	return longitudeOffset(ecl.Elon - targetLon), nil
}

// SearchSunLongitude finds the next time the Sun reaches the given apparent
// ecliptic longitude in degrees, after startTime.
func SearchSunLongitude(targetLon float64, startTime AstroTime, limitDays float64) (AstroTime, error) {
	t2 := startTime.AddDays(limitDays)
	return Search(func(t AstroTime) (float64, error) {
		return sunOffset(targetLon, t)
	}, startTime, t2, 0.01)
}

// Seasons calculates the March and September equinoxes and the June and
// December solstices for the given calendar year.
func Seasons(year int) (SeasonsInfo, error) {
	marEquinox, err := findSeasonChange(0, year, 3, 10)
	if err != nil {
		return SeasonsInfo{}, err
	}
	junSolstice, err := findSeasonChange(90, year, 6, 10)
	if err != nil {
		return SeasonsInfo{}, err
	}
	sepEquinox, err := findSeasonChange(180, year, 9, 10)
	if err != nil {
		return SeasonsInfo{}, err
	}
	decSolstice, err := findSeasonChange(270, year, 12, 10)
	if err != nil {
		return SeasonsInfo{}, err
	}
	return SeasonsInfo{
		MarEquinox:  marEquinox,
		JunSolstice: junSolstice,
		SepEquinox:  sepEquinox,
		DecSolstice: decSolstice,
	}, nil
}

// Rotation_EQJ_ECT returns a rotation matrix that converts J2000 mean equator
// (EQJ) coordinates to true ecliptic of date (ECT) coordinates.
func Rotation_EQJ_ECT(time AstroTime) RotationMatrix {
	return CombineRotation(Rotation_EQJ_EQD(time), Rotation_EQD_ECT(time))
}

// Rotation_ECT_EQJ returns a rotation matrix that converts true ecliptic of
// date (ECT) coordinates to J2000 mean equator (EQJ) coordinates.
func Rotation_ECT_EQJ(time AstroTime) RotationMatrix {
	return InverseRotation(Rotation_EQJ_ECT(time))
}

// Rotation_HOR_EQJ returns a rotation matrix that converts horizontal (HOR)
// coordinates to J2000 mean equator (EQJ) coordinates.
func Rotation_HOR_EQJ(time AstroTime, observer Observer) RotationMatrix {
	return CombineRotation(Rotation_HOR_EQD(time, observer), Rotation_EQD_EQJ(time))
}

// Rotation_EQD_ECL returns a rotation matrix that converts equatorial of-date
// (EQD) coordinates to J2000 mean ecliptic (ECL) coordinates.
func Rotation_EQD_ECL(time AstroTime) RotationMatrix {
	return CombineRotation(Rotation_EQD_EQJ(time), Rotation_EQJ_ECL())
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

//$ASTRO_GO_PLUTO_CONST()

//$ASTRO_GO_PLUTO_TABLE()

//$ASTRO_GO_VSOP(Mercury)
//$ASTRO_GO_VSOP(Venus)
//$ASTRO_GO_VSOP(Earth)
//$ASTRO_GO_VSOP(Mars)
//$ASTRO_GO_VSOP(Jupiter)
//$ASTRO_GO_VSOP(Saturn)
//$ASTRO_GO_VSOP(Uranus)
//$ASTRO_GO_VSOP(Neptune)

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

//$ASTRO_GO_JUPITER_MOONS()

//$ASTRO_CONSTEL()

var moonAddSolTerms = [...]addSolTerm{
	//$ASTRO_ADDSOL()
}
