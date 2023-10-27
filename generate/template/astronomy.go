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
	DaysPerTropicalYear       = 365.24217 // the number of days in one tropical year
	SecondsPerDay             = 86400.0   // the number of seconds in one day
	SolarDaysPerSiderealDay   = 0.9972695717592592
	SpeedOfLightAuPerDay      = 173.1446326846693                         // the speed of light in vacuum expressed in astronomical units per day
	KmPerAu                   = 1.4959787069098932e+8                     // the number of kilometers in one astronomical unit
	AuPerLightYear            = 63241.07708807546                         // the number of astronomical units in one light year
	AsecToRad                 = 4.848136811095359935899141e-6             // factor to convert arcseconds to radians
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
	//  Masses of the Sun and outer planets, used for:
	//  (1) Calculating the Solar System Barycenter
	//  (2) Integrating the movement of Pluto
	//
	//  https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
	//
	//  Page 10 in the above document describes the constants used in the DE405 ephemeris.
	//  The following are G*M values (gravity constant * mass) in [au^3 / day^2].
	//  This side-steps issues of not knowing the exact values of G and masses M[i];
	//  the products GM[i] are known extremely accurately.

	earthMoonMassRatio = 81.30056

	gmSun     = 0.2959122082855911e-03
	gmMercury = 0.4912547451450812e-10
	gmVenus   = 0.7243452486162703e-09
	gmEarth   = 0.8887692390113509e-09
	gmMoon    = gmEarth / earthMoonMassRatio
	gmMars    = 0.9549535105779258e-10
	gmJupiter = 0.2825345909524226e-06
	gmSaturn  = 0.8459715185680659e-07
	gmUranus  = 0.1292024916781969e-07
	gmNeptune = 0.1524358900784276e-07
	gmPluto   = 0.2188699765425970e-11
)

const (
	arc       = 3600.0 * 180.0 / math.Pi // arcseconds per radian
	asecToRad = 1.0 / arc                // radians per arcsecond
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

// Returns the scalar dot product of two vectors.
func Dot(a, b AstroVector) float64 {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z
}

// Returns the length of vec expressed in the same distance units as vec's components.
func (vec AstroVector) Length() float64 {
	return math.Sqrt(Dot(vec, vec))
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

type EclipseKind int

// The different kinds of lunar/solar eclipses.
const (
	NoEclipse        = iota // no eclipse found
	PenumbralEclipse        // A penumbral lunar eclipse. (Never used for a solar eclipse.)
	PartialEclipse          // A partial lunar/solar eclipse.
	AnnularEclipse          // An annular solar eclipse. (Never used for a lunar eclipse.)
	TotalEclipse            // A total lunar/solar eclipse.
)

type Body int

const (
	InvalidBody Body  = -1
	Mercury           = iota // The planet Mercury
	Venus                    // The planet Venus
	Earth                    // The planet Earth
	Mars                     // The planet Mars
	Jupiter                  // The planet Jupiter
	Saturn                   // The planet Saturn
	Uranus                   // The planet Uranus
	Neptune                  // The planet Neptune
	Pluto                    // The dwarf planet Pluto
	Sun                      // The Sun
	Moon                     // The Earth's Moon
	Emb                      // The Earth/Moon Barycenter
	Ssb                      // The Solar System Barycenter
	Star1       = 101        // User-defined star #1
	Star2       = 102        // User-defined star #2
	Star3       = 103        // User-defined star #3
	Star4       = 104        // User-defined star #4
	Star5       = 105        // User-defined star #5
	Star6       = 106        // User-defined star #6
	Star7       = 107        // User-defined star #7
	Star8       = 108        // User-defined star #8
)

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

// Creates a rotation matrix that represents no rotation at all.
func IdentityMatrix() RotationMatrix {
	r := RotationMatrix{}
	r.Rot[0][0] = 1.0
	r.Rot[0][1] = 0.0
	r.Rot[0][2] = 0.0
	r.Rot[1][0] = 0.0
	r.Rot[1][1] = 1.0
	r.Rot[1][2] = 0.0
	r.Rot[2][0] = 0.0
	r.Rot[2][1] = 0.0
	r.Rot[2][2] = 1.0
	return r
}

// Calculates the inverse of a rotation matrix.
// Given a rotation matrix that performs some coordinate transform,
// this function returns the matrix that reverses that transform.
func InverseRotation(rotation RotationMatrix) RotationMatrix {
	inverse := RotationMatrix{}
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

func (tv *terseVector) increment(other terseVector) {
	tv.X += other.X
	tv.Y += other.Y
	tv.Z += other.Z
}

func (left terseVector) sub(right terseVector) terseVector {
	return terseVector{left.X - right.X, left.Y - right.Y, left.Z - right.Z}
}

func (tv terseVector) quadrature() float64 {
	return tv.X*tv.X + tv.Y*tv.Y + tv.Z*tv.Z
}

func (tv terseVector) timesScalar(k float64) terseVector {
	return terseVector{k * tv.X, k * tv.Y, k * tv.Z}
}

type bodyState struct {
	Tt float64     // Terrestrial Time in J2000 days
	R  terseVector // position [au]
	V  terseVector // velocity [au/day]
}

type majorBodies struct {
	Sun     bodyState
	Jupiter bodyState
	Saturn  bodyState
	Uranus  bodyState
	Neptune bodyState
}

func (m *majorBodies) accelerationIncrement(smallBodyPos terseVector, gm float64, majorBodyPos terseVector) terseVector {
	var delta terseVector
	delta = majorBodyPos.sub(smallBodyPos)
	r2 := delta.quadrature()
	return delta.timesScalar(gm / (r2 * math.Sqrt(r2)))
}

func (m *majorBodies) majorBodiesAcceleration(smallBodyPos terseVector) terseVector {
	// Use barycentric coordinates of the Sun and major planets to calculate
	// the gravitational acceleration vector experienced by a small body at location 'smallBodyPos'.
	var accel terseVector
	accel = m.accelerationIncrement(smallBodyPos, gmSun, m.Sun.R)
	accel.increment(m.accelerationIncrement(smallBodyPos, gmJupiter, m.Jupiter.R))
	accel.increment(m.accelerationIncrement(smallBodyPos, gmSaturn, m.Saturn.R))
	accel.increment(m.accelerationIncrement(smallBodyPos, gmUranus, m.Uranus.R))
	accel.increment(m.accelerationIncrement(smallBodyPos, gmNeptune, m.Neptune.R))
	return accel
}

func exportState(terse bodyState, time AstroTime) StateVector {
	return StateVector{
		terse.R.X, terse.R.Y, terse.R.Z,
		terse.V.X, terse.V.Y, terse.V.Z,
		time,
	}
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

type vsopTerm struct {
	amplitude float64
	phase     float64
	frequency float64
}

type vsopSeries struct {
	term []vsopTerm
}

type vsopFormula struct {
	series []vsopSeries
}

type vsopModel struct {
	lon vsopFormula
	lat vsopFormula
	rad vsopFormula
}

type jupiterMoon struct {
	mu   float64
	al0  float64
	al1  float64
	a    []vsopTerm
	l    []vsopTerm
	z    []vsopTerm
	zeta []vsopTerm
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

// Returns the product of mass and universal gravitational constant of a Solar System body.
// For problems involving the gravitational interactions of Solar System bodies,
// it is helpful to know the product GM, where G = the universal gravitational constant
// and M = the mass of the body. In practice, GM is known to a higher precision than
// either G or M alone, and thus using the product results in the most accurate results.
// This function returns the product GM in the units au^3/day^2.
// The values come from page 10 of a JPL memorandum regarding the DE405/LE405 ephemeris:
// https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
func MassProduct(body Body) float64 {
	switch body {
	case Sun:
		return gmSun
	case Mercury:
		return gmMercury
	case Venus:
		return gmVenus
	case Earth:
		return gmEarth
	case Moon:
		return gmMoon
	case Emb:
		return gmEarth + gmMoon
	case Mars:
		return gmMars
	case Jupiter:
		return gmJupiter
	case Saturn:
		return gmSaturn
	case Uranus:
		return gmUranus
	case Neptune:
		return gmNeptune
	case Pluto:
		return gmPluto
	default:
		return -1.0 // invalid body
	}
}

func longitudeOffset(diff float64) float64 {
	offset := diff
	for offset <= -180.0 {
		offset += 360.0
	}
	for offset > 180.0 {
		offset -= 360.0
	}
	return offset
}

func normalizeLongitude(lon float64) float64 {
	for lon < 0.0 {
		lon += 360.0
	}
	for lon >= 360.0 {
		lon -= 360.0
	}
	return lon
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

type starDef struct {
	ra   float64
	dec  float64
	dist float64
}

var starTable [8]starDef

func getStar(body Body) *starDef {
	if body >= Star1 && body <= Star8 {
		return &starTable[body-Star1]
	}
	return nil
}

func userDefinedStar(body Body) *starDef {
	star := getStar(body)
	if star != nil && star.dist > 0.0 {
		return star
	}
	return nil
}

func DefineStar(body Body, ra, dec, distanceLightYears float64) error {
	star := getStar(body)
	if star == nil {
		return errors.New("Invalid body value for a user-defined star.")
	}
	if !isfinite(ra) || ra < 0.0 || ra >= 24.0 {
		return errors.New("Invalid right ascension for user-defined star.")
	}
	if !isfinite(dec) || dec < -90.0 || dec > +90.0 {
		return errors.New("Invalid declination for user-defined star.")
	}
	if !isfinite(distanceLightYears) || distanceLightYears < 1.0 {
		return errors.New("Invalid heliocentric distance for user-defined star. Must be at least 1 light year.")
	}
	star.ra = ra
	star.dec = dec
	star.dist = distanceLightYears * AuPerLightYear
	return nil // success
}

func isSuperiorPlanet(body Body) bool {
	switch body {
	case Mars:
	case Jupiter:
	case Saturn:
	case Uranus:
	case Neptune:
	case Pluto:
		return true
	}
	return false
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

func nutationRot(time *AstroTime, dir precessDirection) RotationMatrix {
	tilt := etilt(time)
	oblm := RadiansFromDegrees(tilt.Mobl)
	oblt := RadiansFromDegrees(tilt.Tobl)
	psi := tilt.Dpsi * AsecToRad

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

	r := RotationMatrix{}

	switch {
	case dir == from2000:
		{
			// convert J2000 to of-date
			r.Rot[0][0] = xx
			r.Rot[0][1] = xy
			r.Rot[0][2] = xz
			r.Rot[1][0] = yx
			r.Rot[1][1] = yy
			r.Rot[1][2] = yz
			r.Rot[2][0] = zx
			r.Rot[2][1] = zy
			r.Rot[2][2] = zz
		}
	case dir == into2000:
		{
			// convert of-date to J2000
			r.Rot[0][0] = xx
			r.Rot[0][1] = yx
			r.Rot[0][2] = zx
			r.Rot[1][0] = xy
			r.Rot[1][1] = yy
			r.Rot[1][2] = zy
			r.Rot[2][0] = xz
			r.Rot[2][1] = yz
			r.Rot[2][2] = zz
		}
	default:
		panic("Invalid precess direction.")
	}

	return r
}

func nutation(pos AstroVector, dir precessDirection) AstroVector {
	rot := nutationRot(&pos.T, dir)
	return RotateVector(rot, pos)
}

func nutationPosVel(state StateVector, dir precessDirection) StateVector {
	rot := nutationRot(&state.T, dir)
	return RotateState(rot, state)
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

func RotateState(rotation RotationMatrix, state StateVector) StateVector {
	return StateVector{
		rotation.Rot[0][0]*state.X + rotation.Rot[1][0]*state.Y + rotation.Rot[2][0]*state.Z,
		rotation.Rot[0][1]*state.X + rotation.Rot[1][1]*state.Y + rotation.Rot[2][1]*state.Z,
		rotation.Rot[0][2]*state.X + rotation.Rot[1][2]*state.Y + rotation.Rot[2][2]*state.Z,
		rotation.Rot[0][0]*state.Vx + rotation.Rot[1][0]*state.Vy + rotation.Rot[2][0]*state.Vz,
		rotation.Rot[0][1]*state.Vx + rotation.Rot[1][1]*state.Vy + rotation.Rot[2][1]*state.Vz,
		rotation.Rot[0][2]*state.Vx + rotation.Rot[1][2]*state.Vy + rotation.Rot[2][2]*state.Vz,
		state.T,
	}
}

// CombineRotation combines the effects of two consecutive rotation matrices into a single rotation matrix.
func CombineRotation(a, b RotationMatrix) RotationMatrix {

	// Use matrix multiplication: c = b*a.
	// We put 'b' on the left and 'a' on the right because,
	// just like when you use a matrix M to rotate a vector V,
	// you put the M on the left in the product M*V.
	// We can think of this as 'b' rotating all the 3 column vectors in 'a'.

	c := RotationMatrix{}
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

// Calculates a rotation matrix that converts equator-of-date (EQD) to J2000 mean equator (EQJ).
func RotationEqdEqj(time *AstroTime) RotationMatrix {
	prec := precessionRot(*time, from2000)
	nut := nutationRot(time, from2000)
	return CombineRotation(prec, nut)
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

//--- Eclipse/transit code begins

type shadowInfo struct {
	time   AstroTime
	u      float64     // dot product of (heliocentric earth) and (geocentric moon): defines the shadow plane where the Moon is
	r      float64     // km distance between center of Moon and the line passing through the centers of the Sun and Earth.
	k      float64     // umbra radius in km, at the shadow plane
	p      float64     // penumbra radius in km, at the shadow plane
	target AstroVector // coordinates of target body relative to shadow-casting body at 'time'
	dir    AstroVector // heliocentric coordinates of shadow-casting body at 'time'
}

func calcShadow(bodyRadiusKm float64, time AstroTime, target AstroVector, dir AstroVector) shadowInfo {
	u := Dot(dir, target) / Dot(dir, dir)
	dx := (u * dir.X) - target.X
	dy := (u * dir.Y) - target.Y
	dz := (u * dir.Z) - target.Z
	r := KmPerAu * math.Sqrt(dx*dx+dy*dy+dz*dz)
	k := +SunRadiusKm - (1.0+u)*(SunRadiusKm-bodyRadiusKm)
	p := -SunRadiusKm + (1.0+u)*(SunRadiusKm-bodyRadiusKm)
	return shadowInfo{time, u, r, k, p, target, dir}
}

func eclipseKindFromUmbra(k float64) EclipseKind {
	// The umbra radius tells us what kind of eclipse the observer sees.
	// If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular.
	// HACK: I added a tiny bias (14 meters) to match Espenak test data.
	if k > 0.014 {
		return TotalEclipse
	}
	return AnnularEclipse
}

//--- Eclipse/transit code ends

func jupiterMoonElemToPv(time AstroTime, mu, A, AL, K, H, Q, P float64) StateVector {
	// Translation of FORTRAN subroutine ELEM2PV from:
	// https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

	AN := math.Sqrt(mu / (A * A * A))

	var CE, SE, DE float64
	EE := AL + K*math.Sin(AL) - H*math.Cos(AL)
	for {
		CE = math.Cos(EE)
		SE = math.Sin(EE)
		DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE)
		EE += DE
		if math.Abs(DE) < 1.0e-12 {
			break
		}
	}

	CE = math.Cos(EE)
	SE = math.Sin(EE)
	DLE := H*CE - K*SE
	RSAM1 := -K*CE - H*SE
	ASR := 1.0 / (1.0 + RSAM1)
	PHI := math.Sqrt(1.0 - K*K - H*H)
	PSI := 1.0 / (1.0 + PHI)
	X1 := A * (CE - K - PSI*H*DLE)
	Y1 := A * (SE - H + PSI*K*DLE)
	VX1 := AN * ASR * A * (-SE - PSI*H*RSAM1)
	VY1 := AN * ASR * A * (+CE + PSI*K*RSAM1)
	F2 := 2.0 * math.Sqrt(1.0-Q*Q-P*P)
	P2 := 1.0 - 2.0*P*P
	Q2 := 1.0 - 2.0*Q*Q
	PQ := 2.0 * P * Q

	return StateVector{
		X1*P2 + Y1*PQ,
		X1*PQ + Y1*Q2,
		(Q*Y1 - X1*P) * F2,
		VX1*P2 + VY1*PQ,
		VX1*PQ + VY1*Q2,
		(Q*VY1 - VX1*P) * F2,
		time,
	}
}

func calcJupiterMoon(time AstroTime, m *jupiterMoon) StateVector {
	// This is a translation of FORTRAN code by Duriez, Lainey, and Vienne:
	// https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

	t := time.Tt + 18262.5 // number of days since 1950-01-01T00:00:00Z

	// Calculate 6 orbital elements at the given time t.
	elem0 := 0.0
	for _, term := range m.a {
		elem0 += term.amplitude * math.Cos(term.phase+(t*term.frequency))
	}

	elem1 := m.al0 + (t * m.al1)
	for _, term := range m.l {
		elem1 += term.amplitude * math.Sin(term.phase+(t*term.frequency))
	}

	elem1 = math.Mod(elem1, 2.0*math.Pi)
	if elem1 < 0 {
		elem1 += 2.0 * math.Pi
	}

	elem2 := 0.0
	elem3 := 0.0
	for _, term := range m.z {
		arg := term.phase + (t * term.frequency)
		elem2 += term.amplitude * math.Cos(arg)
		elem3 += term.amplitude * math.Sin(arg)
	}

	elem4 := 0.0
	elem5 := 0.0
	for _, term := range m.zeta {
		arg := term.phase + (t * term.frequency)
		elem4 += term.amplitude * math.Cos(arg)
		elem5 += term.amplitude * math.Sin(arg)
	}

	// Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP).
	state := jupiterMoonElemToPv(time, m.mu, elem0, elem1, elem2, elem3, elem4, elem5)

	// Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ).
	return RotateState(rotationJupEqj, state)
}

// Calculates Jovicentric positoins and velocities of Jupiter's largest 4 moons.
// Calculates position and velocity vectors for Jupiter's moons
// Io, Europa, Ganymede, and Callisto, at the given date and time.
// The vectors are jovicentric (relative to the center of Jupiter).
// Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
// The position components are expressed in astronomical units (AU), and
// the velocity components are in AU/day.
// To convert to heliocentric position vectors, call HelioVector with
// Jupiter as the body to get Jupiter's heliocentric position,
// then add the jovicentric moon positions.
// Likewise, you can call #Astronomy.GeoVector
// to convert to geocentric positions; however, you will have to manually
// correct for light travel time from the Jupiter system to Earth to
// figure out what time to pass to `JupiterMoons` to get an accurate picture
// of how Jupiter and its moons look from Earth.
func JupiterMoons(time AstroTime) JupiterMoonsInfo {
	return JupiterMoonsInfo{
		Io:       calcJupiterMoon(time, &jupiterMoonModel[0]),
		Europa:   calcJupiterMoon(time, &jupiterMoonModel[1]),
		Ganymede: calcJupiterMoon(time, &jupiterMoonModel[2]),
		Callisto: calcJupiterMoon(time, &jupiterMoonModel[3]),
	}
}

//--- Pluto calc begins

func clampIndex(frac float64, nsteps int) int {
	var index int
	index = int(math.Floor(frac))
	if index < 0 {
		return 0
	}
	if index >= nsteps {
		return nsteps - 1
	}
	return index
}

//--- Pluto calc ends

type ascentInfo struct {
	valid bool
	tx    AstroTime
	ty    AstroTime
	ax    float64
	ay    float64
}

type SearchContext interface {
	Eval(time AstroTime) (float64, error)
}

func findAscent(depth int, context SearchContext, maxDerivAlt float64, t1, t2 AstroTime, a1, a2 float64) ascentInfo {
	// See if we can find any time interval where the altitude-diff function
	// rises from non-positive to positive.

	if a1 < 0.0 && a2 >= 0.0 {
		// Trivial success case: the endpoints already rise through zero.
		return ascentInfo{valid: true, tx: t1, ty: t2, ax: a1, ay: a2}
	}

	if a1 >= 0.0 && a2 < 0.0 {
		// Trivial failure case: Assume Nyquist condition prevents an ascent.
		return ascentInfo{valid: false}
	}

	if depth > 17 {
		// Safety valve: do not allow unlimited recursion.
		// This should never happen if the rest of the logic is working correctly,
		// so fail the whole search if it does happen. It's a bug!
		panic("Excessive recursion in rise/set ascent search.")
	}

	// Both altitudes are on the same side of zero: both are negative, or both are non-negative.
	// There could be a convex "hill" or a concave "valley" that passes through zero.
	// In polar regions sometimes there is a rise/set or set/rise pair within minutes of each other.
	// For example, the Moon can be below the horizon, then the very top of it becomes
	// visible (moonrise) for a few minutes, then it moves sideways and down below
	// the horizon again (moonset). We want to catch these cases.
	// However, for efficiency and practicality concerns, because the rise/set search itself
	// has a 0.1 second threshold, we do not worry about rise/set pairs that are less than
	// one second apart. These are marginal cases that are rendered highly uncertain
	// anyway, due to unpredictable atmospheric refraction conditions (air temperature and pressure).

	dt := t2.Ut - t1.Ut
	if dt*SecondsPerDay < 1.0 {
		return ascentInfo{valid: false}
	}

	// Is it possible to reach zero from the altitude that is closer to zero?
	da := math.Min(math.Abs(a1), math.Abs(a2))

	// Without loss of generality, assume |a1| <= |a2|.
	// (Reverse the argument in the case |a2| < |a1|.)
	// Imagine you have to "drive" from a1 to 0, then back to a2.
	// You can't go faster than max_deriv_alt. If you can't reach 0 in half the time,
	// you certainly don't have time to reach 0, turn around, and still make your way
	// back up to a2 (which is at least as far from 0 than a1 is) in the time interval dt.
	// Therefore, the time threshold is half the time interval, or dt/2.
	if da > maxDerivAlt*(dt/2) {
		// Prune: the altitude cannot change fast enough to reach zero.
		return ascentInfo{valid: false}
	}

	// Bisect the time interval and evaluate the altitude at the midpoint.
	tmid := TimeFromUniversalDays((t1.Ut + t2.Ut) / 2.0)
	amid, err := context.Eval(tmid)
	if err != nil {
		panic("Altitude context should not have failed to evaluate.")
	}

	// Recurse to the left interval.
	ascent := findAscent(1+depth, context, maxDerivAlt, t1, tmid, a1, amid)
	if !ascent.valid {
		// Recurse to the right interval.
		ascent = findAscent(1+depth, context, maxDerivAlt, tmid, t2, amid, a2)
	}

	return ascent
}

func horizonDipAngle(observer Observer, metersAboveGround float64) float64 {
	// Calculate the effective radius of the Earth at ground level below the observer.
	// Correct for the Earth's oblateness.
	phi := RadiansFromDegrees(observer.Latitude)
	sinphi := math.Sin(phi)
	cosphi := math.Cos(phi)
	c := 1.0 / math.Hypot(cosphi, sinphi*EarthFlattening)
	s := c * (EarthFlattening * EarthFlattening)
	htkm := (observer.Height - metersAboveGround) / 1000.0 // height of ground above sea level
	ach := EarthEquatorialRadiusKm*c + htkm
	ash := EarthEquatorialRadiusKm*s + htkm
	radiusMeters := 1000.0 * math.Hypot(ach*cosphi, ash*sinphi)

	// Correct refraction of a ray of light traveling tangent to the Earth's surface.
	// Based on: https://www.largeformatphotography.info/sunmooncalc/SMCalc.js
	// which in turn derives from:
	// Sweer, John. 1938.  The Path of a Ray of Light Tangent to the Surface of the Earth.
	// Journal of the Optical Society of America 28 (September):327-329.

	// k = refraction index
	k := 0.175 * math.Pow(1.0-(6.5e-3/283.15)*(observer.Height-(2.0/3.0)*metersAboveGround), 3.256)

	// Calculate how far below the observer's horizontal plane the observed horizon dips.
	return DegreesFromRadians(-(math.Sqrt(2*(1-k)*metersAboveGround/radiusMeters) / (1 - k)))
}

func maxAltitudeSlope(body Body, latitude float64) float64 {
	// Calculate the maximum possible rate that this body's altitude
	// could change [degrees/day] as seen by this observer.
	// First use experimentally determined extreme bounds for this body
	// of how much topocentric RA and DEC can ever change per rate of time.
	// We need minimum possible d(RA)/dt, and maximum possible magnitude of d(DEC)/dt.
	// Conservatively, we round d(RA)/dt down, d(DEC)/dt up.
	// Then calculate the resulting maximum possible altitude change rate.

	var derivRa, derivDec float64

	switch body {
	case Moon:
		derivRa = +4.5
		derivDec = +8.2
		break

	case Sun:
		derivRa = +0.8
		derivDec = +0.5
		break

	case Mercury:
		derivRa = -1.6
		derivDec = +1.0
		break

	case Venus:
		derivRa = -0.8
		derivDec = +0.6
		break

	case Mars:
		derivRa = -0.5
		derivDec = +0.4
		break

	case Jupiter:
	case Saturn:
	case Uranus:
	case Neptune:
	case Pluto:
		derivRa = -0.2
		derivDec = +0.2
		break

	case Star1:
	case Star2:
	case Star3:
	case Star4:
	case Star5:
	case Star6:
	case Star7:
	case Star8:
		// The minimum allowed heliocentric distance of a user-defined star
		// is one light-year. This can cause a tiny amount of parallax (about 0.001 degrees).
		// Also, including stellar aberration (22 arcsec = 0.006 degrees), we provide a
		// generous safety buffer of 0.008 degrees.
		derivRa = -0.008
		derivDec = +0.008
		break

	default:
		return -1.0 // invalid body
	}

	latrad := RadiansFromDegrees(latitude)
	return math.Abs(((360.0/SolarDaysPerSiderealDay)-derivRa)*math.Cos(latrad)) + math.Abs(derivDec*math.Sin(latrad))
}

func moonMagnitude(phase float64, helioDist float64, geoDist float64) float64 {
	// https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
	rad := RadiansFromDegrees(phase)
	rad2 := rad * rad
	rad4 := rad2 * rad2
	mag := -12.717 + 1.49*math.Abs(rad) + 0.0431*rad4
	moonMeanDistanceAu := 385000.6 / KmPerAu
	geo_au := geoDist / moonMeanDistanceAu
	mag += 5.0 * math.Log10(helioDist*geo_au)
	return mag
}

//--- Generated code begins here ------------------------------------------------------------------

//$ASTRO_CONSTEL()

var moonAddSolTerms = [...]addSolTerm{
	//$ASTRO_ADDSOL()
}

var jupiterMoonModel = [...]jupiterMoon{
	//$ASTRO_JUPITER_MOONS()
}

var rotationJupEqj = RotationMatrix{
	[3][3]float64{
		//$ASTRO_JUP_EQJ_ROT()
	},
}
