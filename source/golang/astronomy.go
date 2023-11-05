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
	arc                  = 3600.0 * 180.0 / math.Pi // arcseconds per radian
	asecToRad            = 1.0 / arc                // radians per arcsecond
	earthOrbitalPeriod   = 365.256
	neptuneOrbitalPeriod = 60189.0
	daysPerMillennium    = 365250.0
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
		return nil, errors.New("non-finite value of ut")
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
		return nil, errors.New("the supplied time is too far from the year 2000 to be represented.")
	}

	if month < 1 || month > 12 || day < 1 || day > 31 {
		return nil, errors.New("internal error: invalid calendar date calculated.")
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
		return AtmosphereInfo{}, errors.New("invalid elevation")
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

// Converts Cartesian coordinates to spherical coordinates.
func SphereFromVector(vector AstroVector) Spherical {
	xyproj := vector.X*vector.X + vector.Y*vector.Y
	dist := math.Sqrt(xyproj + vector.Z*vector.Z)
	var lat, lon float64
	if xyproj == 0.0 {
		if vector.Z == 0.0 {
			// Indeterminate coordinates; pos vector has zero length.
			return Spherical{0.0, 0.0, 0.0}
		}
		lon = 0.0
		if vector.Z < 0.0 {
			lat = -90.0
		} else {
			lat = +90.0
		}
	} else {
		lon = datan2(vector.Y, vector.X)
		if lon < 0.0 {
			lon += 360.0
		}
		lat = datan2(vector.Z, math.Sqrt(xyproj))
	}
	return Spherical{lat, lon, dist}
}

// Given an equatorial vector, calculates equatorial angular coordinates.
func EquatorFromVector(vector AstroVector) Equatorial {
	sphere := SphereFromVector(vector)
	return Equatorial{
		Ra:   sphere.Lon / 15.0, // convert degrees to sidereal hours
		Dec:  sphere.Lat,
		Dist: sphere.Dist,
		Vec:  vector,
	}
}

// Converts cartesian coordinates to horizontal coordinates.
// Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
//
// IMPORTANT: This function differs from SphereFromVector in two ways:
//   - SphereFromVector returns a `Lon` value that represents azimuth defined counterclockwise
//     from north (e.g., west = +90), but this function represents a clockwise rotation
//     (e.g., east = +90). The difference is because `SphereFromVector` is intended
//     to preserve the vector "right-hand rule", while this function defines azimuth in a more
//     traditional way as used in navigation and cartography.
//   - This function optionally corrects for atmospheric refraction, while `SphereFromVector`
//     does not.
//
// The returned structure contains the azimuth in `Lon`.
// It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
//
// The altitude is stored in `Lat`.
//
// The distance to the observed object is stored in `Dist`,
// and is expressed in astronomical units (AU).
func HorizonFromVector(vector AstroVector, refraction Refraction) Spherical {
	sphere := SphereFromVector(vector)
	return Spherical{
		sphere.Lat + RefractionAngle(refraction, sphere.Lat),
		toggleAzimuthDirection(sphere.Lon),
		sphere.Dist,
	}
}

// Converts spherical coordinates to Cartesian coordinates.
func VectorFromSphere(sphere Spherical, time AstroTime) AstroVector {
	radlat := RadiansFromDegrees(sphere.Lat)
	radlon := RadiansFromDegrees(sphere.Lon)
	rcoslat := sphere.Dist * math.Cos(radlat)
	return AstroVector{
		rcoslat * math.Cos(radlon),
		rcoslat * math.Sin(radlon),
		sphere.Dist * math.Sin(radlat),
		time,
	}
}

func toggleAzimuthDirection(az float64) float64 {
	az = 360.0 - az
	if az >= 360.0 {
		az -= 360.0
	} else if az < 0.0 {
		az += 360.0
	}
	return az
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

// RefractionAngle calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
// Given an altitude angle and a refraction option, calculates
// the amount of "lift" caused by atmospheric refraction.
// This is the number of degrees higher in the sky an object appears
// due to the lensing of the Earth's atmosphere.
// This function works best near sea level.
// To correct for higher elevations, call Atmosphere for that
// elevation and multiply the refraction angle by the resulting relative density.
// The refraction parameter specifies which refraction correction to use.
// If set to NormalRefraction, uses a well-behaved refraction model that works well for
// all valid values (-90 to +90) of altitude.
// If set to JplHorizonsRefraction, this function returns a value compatible with the JPL Horizons tool.
// This is provided for internal unit tests that compare against JPL Horizons data.
// Any other value, including NoRefraction, causes this function to return 0.0.
// The return value is a non-negative value expressed in degrees of refraction above the horizontal.
func RefractionAngle(refraction Refraction, altitude float64) float64 {
	if altitude < -90.0 || altitude > +90.0 {
		return 0.0 // no attempt to correct an invalid altitude
	}
	var refr float64
	if refraction == NormalRefraction || refraction == JplHorizonsRefraction {
		// http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
		// JPL Horizons says it uses refraction algorithm from
		// Meeus "Astronomical Algorithms", 1991, p. 101-102.
		// I found the following Go implementation:
		// https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
		// This is a translation from the function "Saemundsson" there.
		// I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
		// This is important because the 'refr' formula below goes crazy near hd = -5.11.
		hd := altitude
		if hd < -1.0 {
			hd = -1.0
		}

		refr = (1.02 / dtan(hd+10.3/(hd+5.11))) / 60.0

		if refraction == NormalRefraction && altitude < -1.0 {
			// In "normal" mode we gradually reduce refraction toward the nadir
			// so that we never get an altitude angle less than -90 degrees.
			// When horizon angle is -1 degrees, the factor is exactly 1.
			// As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
			refr *= (altitude + 90.0) / 89.0
		}
	} else {
		// No refraction, or the refraction option is invalid.
		refr = 0.0
	}
	return refr
}

// Calculates the inverse of an atmospheric refraction angle.
// Given an observed altitude angle that includes atmospheric refraction,
// calculates the negative angular correction to obtain the unrefracted
// altitude. This is useful for cases where observed horizontal
// coordinates are to be converted to another orientation system,
// but refraction first must be removed from the observed position.
func InverseRefractionAngle(refraction Refraction, bentAltitude float64) float64 {
	if bentAltitude < -90.0 || bentAltitude > +90.0 {
		return 0.0 // no attempt to correct an invalid altitude
	}

	// Find the pre-adjusted altitude whose refraction correction leads to 'altitude'.
	altitude := bentAltitude - RefractionAngle(refraction, bentAltitude)
	for {
		diff := (altitude + RefractionAngle(refraction, altitude)) - bentAltitude
		if math.Abs(diff) < 1.0e-14 {
			return altitude - bentAltitude
		}
		altitude -= diff
	}
}

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

func (tv terseVector) toAstroVector(time AstroTime) AstroVector {
	return AstroVector{tv.X, tv.Y, tv.Z, time}
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

func updatePosition(dt float64, r, v, a terseVector) terseVector {
	return terseVector{
		r.X + dt*(v.X+dt*a.X/2.0),
		r.Y + dt*(v.Y+dt*a.Y/2.0),
		r.Z + dt*(v.Z+dt*a.Z/2.0),
	}
}

func updateVelocity(dt float64, v, a terseVector) terseVector {
	return terseVector{
		v.X + dt*a.X,
		v.Y + dt*a.Y,
		v.Z + dt*a.Z,
	}
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

func vsopFormulaCalc(formula *vsopFormula, t float64, clampAngle bool) float64 {
	coord := 0.0
	tpower := 1.0
	for _, series := range formula.series {
		sum := 0.0
		for _, term := range series.term {
			sum += term.amplitude * math.Cos(term.phase+(t*term.frequency))
		}
		incr := tpower * sum
		if clampAngle {
			// improve precision: longitude angles can be hundreds of radians
			incr = math.Mod(incr, 2.0*math.Pi)
		}
		coord += incr
		tpower *= t
	}
	return coord
}

func vsopDerivCalc(formula *vsopFormula, t float64) float64 {
	tpower := 1.0 // t^2
	dpower := 0.0 // t^(s-1)
	deriv := 0.0
	for s, series := range formula.series {
		sinSum := 0.0
		cosSum := 0.0
		for _, term := range series.term {
			angle := term.phase + (t * term.frequency)
			sinSum += term.amplitude * term.frequency * math.Sin(angle)
			if s > 0 {
				cosSum += term.amplitude * math.Cos(angle)
			}
		}
		deriv += (float64(s) * dpower * cosSum) - (tpower * sinSum)
		dpower = tpower
		tpower *= t
	}
	return deriv
}

func vsopRotate(eclip terseVector) terseVector {
	return terseVector{
		eclip.X + 0.000000440360*eclip.Y - 0.000000190919*eclip.Z,
		-0.000000479966*eclip.X + 0.917482137087*eclip.Y - 0.397776982902*eclip.Z,
		0.397776982902*eclip.Y + 0.917482137087*eclip.Z,
	}
}

func vsopSphereToRect(lon, lat, radius float64) terseVector {
	rcoslat := radius * math.Cos(lat)
	return terseVector{
		rcoslat * math.Cos(lon),
		rcoslat * math.Sin(lon),
		radius * math.Sin(lat),
	}
}

func calcVsop(model *vsopModel, time AstroTime) AstroVector {
	t := time.Tt / daysPerMillennium // millennia since 2000

	// Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
	lon := vsopFormulaCalc(&model.lon, t, true)
	lat := vsopFormulaCalc(&model.lat, t, false)
	rad := vsopFormulaCalc(&model.rad, t, false)

	// Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates.
	eclip := vsopSphereToRect(lon, lat, rad)

	// Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates.
	return vsopRotate(eclip).toAstroVector(time)
}

func calcVsopPosVel(model *vsopModel, tt float64) bodyState {
	t := tt / daysPerMillennium // millennia since 2000

	// Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
	lon := vsopFormulaCalc(&model.lon, t, true)
	lat := vsopFormulaCalc(&model.lat, t, false)
	rad := vsopFormulaCalc(&model.rad, t, false)

	eclipPos := vsopSphereToRect(lon, lat, rad)

	dlonDt := vsopDerivCalc(&model.lon, t)
	dlatDt := vsopDerivCalc(&model.lat, t)
	dradDt := vsopDerivCalc(&model.rad, t)

	// Use spherical coords and spherical derivatives to calculate
	// the velocity vector in rectangular coordinates.

	coslon := math.Cos(lon)
	sinlon := math.Sin(lon)
	coslat := math.Cos(lat)
	sinlat := math.Sin(lat)

	vx := +(dradDt * coslat * coslon) - (rad * sinlat * coslon * dlatDt) - (rad * coslat * sinlon * dlonDt)
	vy := +(dradDt * coslat * sinlon) - (rad * sinlat * sinlon * dlatDt) + (rad * coslat * coslon * dlonDt)
	vz := +(dradDt * sinlat) + (rad * coslat * dlatDt)

	// Convert speed units from [AU/millennium] to [AU/day].
	eclipVel := terseVector{
		vx / daysPerMillennium,
		vy / daysPerMillennium,
		vz / daysPerMillennium,
	}

	// Rotate the vectors from ecliptic to equatorial coordinates.
	equPos := vsopRotate(eclipPos)
	equVel := vsopRotate(eclipVel)
	return bodyState{tt, equPos, equVel}
}

// The function CorrectLightTravel solves a generalized problem of
// deducing how far in the past light must have left a target object
// to be seen by an observer at a specified time.
// This interface expresses an arbitrary position function as
// a function of time that is passed to CorrectLightTravel.
type PositionFunction interface {
	Position(time AstroTime) (AstroVector, error)
}

// Solves for light travel time, given a vector function representing distance as a function of time.
//
// When observing a distant object, for example Jupiter as seen from Earth,
// the amount of time it takes for light to travel from the object to the
// observer can significantly affect the object's apparent position.
// This function is a generic solver that figures out how long in the
// past light must have left the observed object to reach the observer
// at the specified observation time. It uses PositionFunction
// to express an arbitrary position vector as a function of time.
//
// This function repeatedly calls `func.Position`, passing a series of time
// estimates in the past. Then `func.Position` must return a relative state vector between
// the observer and the target. `CorrectLightTravel` keeps calling
// `func.Position` with more and more refined estimates of the time light must have
// left the target to arrive at the observer.
//
// For common use cases, it is simpler to use BackdatePosition
// for calculating the light travel time correction of one body observing another body.
//
// For geocentric calculations, GeoVector also backdates the returned
// position vector for light travel time, only it returns the observation time in
// the returned vector's `t` field rather than the backdated time.
func CorrectLightTravel(posFunc PositionFunction, time AstroTime) (*AstroVector, error) {
	ltime := time
	for iter := 0; iter < 10; iter++ {
		pos, err := posFunc.Position(ltime)
		if err != nil {
			return nil, err
		}
		// This solver does not support more than one light-day of distance,
		// because that would cause convergence problems and inaccurate
		// values for stellar aberration angles.
		lt := pos.Length() / SpeedOfLightAuPerDay
		if lt > 1.0 {
			return nil, errors.New("object is too distant for light-travel solver.")
		}
		ltime2 := time.AddDays(-lt)
		dt := math.Abs(ltime2.Tt - ltime.Tt)
		if dt < 1.0e-9 { // 86.4 microseconds
			return &pos, nil
		}
		ltime = ltime2
	}
	return nil, errors.New("light travel time correction did not converge.")
}

func Horizon(time AstroTime, observer Observer, ra, dec float64, refraction Refraction) (*Topocentric, error) {
	sinlat := dsin(observer.Latitude)
	coslat := dcos(observer.Latitude)
	sinlon := dsin(observer.Longitude)
	coslon := dcos(observer.Longitude)
	sindc := dsin(dec)
	cosdc := dcos(dec)
	sinra := dsin(ra * 15.0)
	cosra := dcos(ra * 15.0)

	// Calculate three mutually perpendicular unit vectors
	// in equatorial coordinates: uze, une, uwe.
	//
	// uze = The direction of the observer's local zenith (straight up).
	// une = The direction toward due north on the observer's horizon.
	// uwe = The direction toward due west on the observer's horizon.
	//
	// HOWEVER, these are uncorrected for the Earth's rotation due to the time of day.
	//
	// The components of these 3 vectors are as follows:
	// x = direction from center of Earth toward 0 degrees longitude (the prime meridian) on equator.
	// y = direction from center of Earth toward 90 degrees west longitude on equator.
	// z = direction from center of Earth toward the north pole.
	uze := AstroVector{coslat * coslon, coslat * sinlon, sinlat, time}
	une := AstroVector{-sinlat * coslon, -sinlat * sinlon, coslat, time}
	uwe := AstroVector{sinlon, -coslon, 0.0, time}

	// Correct the vectors uze, une, uwe for the Earth's rotation by calculating
	// sidereal time. Call spin() for each uncorrected vector to rotate about
	// the Earth's axis to yield corrected unit vectors uz, un, uw.
	// Multiply sidereal hours by -15 to convert to degrees and flip eastward
	// rotation of the Earth to westward apparent movement of objects with time.
	angle := -15.0 * SiderealTime(&time)
	uz := spin(angle, uze)
	un := spin(angle, une)
	uw := spin(angle, uwe)

	// Convert angular equatorial coordinates (RA, DEC) to
	// cartesian equatorial coordinates in 'p', using the
	// same orientation system as uze, une, uwe.
	p := AstroVector{cosdc * cosra, cosdc * sinra, sindc, time}

	// Use dot products of p with the zenith, north, and west
	// vectors to obtain the cartesian coordinates of the body in
	// the observer's horizontal orientation system.
	// pz = zenith component [-1, +1]
	// pn = north  component [-1, +1]
	// pw = west   component [-1, +1]
	pz := p.X*uz.X + p.Y*uz.Y + p.Z*uz.Z
	pn := p.X*un.X + p.Y*un.Y + p.Z*un.Z
	pw := p.X*uw.X + p.Y*uw.Y + p.Z*uw.Z

	// proj is the "shadow" of the body vector along the observer's flat ground.
	proj := math.Hypot(pn, pw)

	// Calculate az = azimuth (compass direction clockwise from East.)
	var az float64
	if proj > 0.0 {
		// If the body is not exactly straight up/down, it has an azimuth.
		// Invert the angle to produce degrees eastward from north.
		az = -datan2(pw, pn)
		if az < 0.0 {
			az += 360.0
		}
	} else {
		// The body is straight up/down, so it does not have an azimuth.
		// Report an arbitrary but reasonable value.
		az = 0.0
	}

	// zd = the angle of the body away from the observer's zenith, in degrees.
	zd := datan2(proj, pz)
	horRa := ra
	horDec := dec

	if refraction == NormalRefraction || refraction == JplHorizonsRefraction {
		zd0 := zd
		refr := RefractionAngle(refraction, 90.0-zd)
		zd -= refr

		if refr > 0.0 && zd > 3.0e-4 {
			sinzd := dsin(zd)
			coszd := dcos(zd)
			sinzd0 := dsin(zd0)
			coszd0 := dcos(zd0)

			prx := ((p.X-coszd0*uz.X)/sinzd0)*sinzd + uz.X*coszd
			pry := ((p.Y-coszd0*uz.Y)/sinzd0)*sinzd + uz.Y*coszd
			prz := ((p.Z-coszd0*uz.Z)/sinzd0)*sinzd + uz.Z*coszd

			proj = math.Hypot(prx, pry)
			if proj > 0.0 {
				horRa = datan2(pry, prx) / 15.0
				if horRa < 0.0 {
					horRa += 24.0
				}
			} else {
				horRa = 0.0
			}
			horDec = datan2(prz, proj)
		}
	} else if refraction != NoRefraction {
		return nil, errors.New("unsupported refraction option.")
	}

	return &Topocentric{az, 90.0 - zd, horRa, horDec}, nil
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

// PlanetOrbitalPeriod returns the average number of days it takes for a planet to orbit the Sun.
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
		return -1.0 // invalid body
	}
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

func dtan(degrees float64) float64 {
	return math.Tan(RadiansFromDegrees(degrees))
}

func datan2(y, x float64) float64 {
	return DegreesFromRadians(math.Atan2(y, x))
}

func obscuration(a, b, c float64) float64 {
	if a <= 0.0 {
		panic("Radius of first disc must be positive.")
	}
	if b <= 0.0 {
		panic("Radius of second disc must be positive.")
	}
	if c < 0.0 {
		panic("Distance between discs is not allowed to be negative.")
	}
	if c >= a+b {
		// The discs are too far apart to have any overlapping area.
		return 0.0
	}

	if c == 0.0 {
		// The discs have a common center. Therefore, one disc is inside the other.
		if a <= b {
			return 1.0
		}
		return (b * b) / (a * a)
	}

	x := (a*a - b*b + c*c) / (2 * c)
	radicand := a*a - x*x
	if radicand <= 0.0 {
		// The circumferences do not intersect, or are tangent.
		// We already ruled out the case of non-overlapping discs.
		// Therefore, one disc is inside the other.
		if a <= b {
			return 1.0
		}
		return (b * b) / (a * a)
	}

	// The discs overlap fractionally in a pair of lens-shaped areas.

	y := math.Sqrt(radicand)

	// Return the overlapping fractional area.
	// There are two lens-shaped areas, one to the left of x, the other to the right of x.
	// Each part is calculated by subtracting a triangular area from a sector's area.
	lens1 := a*a*math.Acos(x/a) - x*y
	lens2 := b*b*math.Acos((c-x)/b) - (c-x)*y

	// Find the fractional area with respect to the first disc.
	return (lens1 + lens2) / (math.Pi * a * a)
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

// Calculates the gravitational acceleration experienced by an observer on the Earth.
// This function implements the WGS 84 Ellipsoidal Gravity Formula.
// The result is a combination of inward gravitational acceleration
// with outward centrifugal acceleration, as experienced by an observer
// in the Earth's rotating frame of reference.
// The resulting value increases toward the Earth's poles and decreases
// toward the equator, consistent with changes of the weight measured
// by a spring scale of a fixed mass moved to different latitudes and heights
// on the Earth.
// The latitude is of the observer in degrees north or south of the equator.
// By formula symmetry, positive latitudes give the same answer as negative
// latitudes, so the sign does not matter.
// The height is specified above the sea level geoid in meters.
// No range checking is done; however, accuracy is only valid in the
// range 0 to 100000 meters.
// The return value is the gravitational acceleration expressed in meters per second squared.
func ObserverGravity(latitude, height float64) float64 {
	s := dsin(latitude)
	s2 := s * s
	g0 := 9.7803253359 * (1.0 + 0.00193185265241*s2) / math.Sqrt(1.0-0.00669437999013*s2)
	return g0 * (1.0 - (3.15704e-07-2.10269e-09*s2)*height + 7.37452e-14*height*height)
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
		return errors.New("invalid body value for a user-defined star.")
	}
	if !isfinite(ra) || ra < 0.0 || ra >= 24.0 {
		return errors.New("invalid right ascension for user-defined star.")
	}
	if !isfinite(dec) || dec < -90.0 || dec > +90.0 {
		return errors.New("invalid declination for user-defined star.")
	}
	if !isfinite(distanceLightYears) || distanceLightYears < 1.0 {
		return errors.New("invalid heliocentric distance for user-defined star. Must be at least 1 light year.")
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

func rotateEquatorialToEcliptic(pos AstroVector, obliqRadians float64) Ecliptic {
	cosOb := math.Cos(obliqRadians)
	sinOb := math.Sin(obliqRadians)

	ex := +pos.X
	ey := +pos.Y*cosOb + pos.Z*sinOb
	ez := -pos.Y*sinOb + pos.Z*cosOb

	xyproj := math.Hypot(ex, ey)
	elon := 0.0
	if xyproj > 0.0 {
		elon = datan2(ey, ex)
		if elon < 0.0 {
			elon += 360.0
		}
	}

	elat := datan2(ez, xyproj)
	vec := AstroVector{ex, ey, ez, pos.T}
	return Ecliptic{vec, elat, elon}
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

func spin(angle float64, pos AstroVector) AstroVector {
	cosang := dcos(angle)
	sinang := dsin(angle)
	return AstroVector{
		+cosang*pos.X + sinang*pos.Y,
		-sinang*pos.X + cosang*pos.Y,
		pos.Z,
		pos.T,
	}
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

// Pivot re-orients a rotation matrix by pivoting it by an angle around one of its axes.
// Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
// this function pivots the rotation matrix by that angle around that coordinate axis.
// For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
// to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
// of a telescope camera pointed at a given body, you can use Pivot twice:
// (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
// western axis by the body's altitude angle. The resulting rotation matrix will then
// reorient ECL coordinates to the orientation of your telescope camera.
// The axis parameter is an integer that selects which axis to pivot about: 0=x, 1=y, 2=z.
// The angle parameter is an angle in degrees indicating the amount of rotation around the specified axis.
// Positive angles indicate rotation counterclockwise as seen from the positive
// direction along that axis, looking towards the origin point of the orientation system.
// Any finite number of degrees is allowed, but best precision will result from keeping angle
// in the range [-360, +360].
func Pivot(rotation RotationMatrix, axis int, angle float64) (*RotationMatrix, error) {
	if axis < 0 || axis > 2 {
		return nil, errors.New("invalid coordinate axis for Pivot. Must be 0..2.")
	}
	if !isfinite(angle) {
		return nil, errors.New("angle is not a finite number. Not valid for Pivot.")
	}
	c := dcos(angle)
	s := dsin(angle)

	// We need to maintain the "right-hand" rule, no matter which
	// axis was selected. That means we pick (i, j, k) axis order
	// such that the following vector cross product is satisfied:
	// i x j = k
	i := (axis + 1) % 3
	j := (axis + 2) % 3
	k := axis

	r := RotationMatrix{}

	r.Rot[i][i] = c*rotation.Rot[i][i] - s*rotation.Rot[i][j]
	r.Rot[i][j] = s*rotation.Rot[i][i] + c*rotation.Rot[i][j]
	r.Rot[i][k] = rotation.Rot[i][k]

	r.Rot[j][i] = c*rotation.Rot[j][i] - s*rotation.Rot[j][j]
	r.Rot[j][j] = s*rotation.Rot[j][i] + c*rotation.Rot[j][j]
	r.Rot[j][k] = rotation.Rot[j][k]

	r.Rot[k][i] = c*rotation.Rot[k][i] - s*rotation.Rot[k][j]
	r.Rot[k][j] = s*rotation.Rot[k][i] + c*rotation.Rot[k][j]
	r.Rot[k][k] = rotation.Rot[k][k]

	return &r, nil
}

// Calculates a rotation matrix that converts equator-of-date (EQD) to J2000 mean equator (EQJ).
func RotationEqdEqj(time *AstroTime) RotationMatrix {
	prec := precessionRot(*time, from2000)
	nut := nutationRot(time, from2000)
	return CombineRotation(prec, nut)
}

// Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ).
func RotationEclEqj() RotationMatrix {
	// ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
	const c = 0.9174821430670688 // cos(ob)
	const s = 0.3977769691083922 // sin(ob)

	r := RotationMatrix{}
	r.Rot[0][0] = 1.0
	r.Rot[1][0] = 0.0
	r.Rot[2][0] = 0.0
	r.Rot[0][1] = 0.0
	r.Rot[1][1] = +c
	r.Rot[2][1] = -s
	r.Rot[0][2] = 0.0
	r.Rot[1][2] = +s
	r.Rot[2][2] = +c

	return r
}

// Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL).
func RotationEqjEcl() RotationMatrix {
	// ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
	const c = 0.9174821430670688 // cos(ob)
	const s = 0.3977769691083922 // sin(ob)

	r := RotationMatrix{}
	r.Rot[0][0] = 1.0
	r.Rot[1][0] = 0.0
	r.Rot[2][0] = 0.0
	r.Rot[0][1] = 0.0
	r.Rot[1][1] = +c
	r.Rot[2][1] = +s
	r.Rot[0][2] = 0.0
	r.Rot[1][2] = -s
	r.Rot[2][2] = +c

	return r
}

// Calculates a rotation matrix from J2000 mean equator (EQJ) to galactic (GAL).
func RotationEqjGal() RotationMatrix {
	r := RotationMatrix{}

	// This rotation matrix was calculated by the following script
	// in this same source code repository:
	// demo/python/galeqj_matrix.py

	r.Rot[0][0] = -0.0548624779711344
	r.Rot[0][1] = +0.4941095946388765
	r.Rot[0][2] = -0.8676668813529025
	r.Rot[1][0] = -0.8734572784246782
	r.Rot[1][1] = -0.4447938112296831
	r.Rot[1][2] = -0.1980677870294097
	r.Rot[2][0] = -0.4838000529948520
	r.Rot[2][1] = +0.7470034631630423
	r.Rot[2][2] = +0.4559861124470794

	return r
}

func RotationGalEqj() RotationMatrix {
	r := RotationMatrix{}

	// This rotation matrix was calculated by the following script
	// in this same source code repository:
	// demo/python/galeqj_matrix.py

	r.Rot[0][0] = -0.0548624779711344
	r.Rot[0][1] = -0.8734572784246782
	r.Rot[0][2] = -0.4838000529948520
	r.Rot[1][0] = +0.4941095946388765
	r.Rot[1][1] = -0.4447938112296831
	r.Rot[1][2] = +0.7470034631630423
	r.Rot[2][0] = -0.8676668813529025
	r.Rot[2][1] = -0.1980677870294097
	r.Rot[2][2] = +0.4559861124470794

	return r
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

func visualMagnitude(body Body, phase, helioDist, geoDist float64) (float64, error) {
	// For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212
	c0 := 0.0
	c1 := 0.0
	c2 := 0.0
	c3 := 0.0
	switch body {
	case Mercury:
		c0 = -0.60
		c1 = +4.98
		c2 = -4.88
		c3 = +3.02
		break
	case Venus:
		if phase < 163.6 {
			c0 = -4.47
			c1 = +1.03
			c2 = +0.57
			c3 = +0.13
		} else {
			c0 = 0.98
			c1 = -1.02
		}
		break
	case Mars:
		c0 = -1.52
		c1 = +1.60
		break
	case Jupiter:
		c0 = -9.40
		c1 = +0.50
		break
	case Uranus:
		c0 = -7.19
		c1 = +0.25
		break
	case Neptune:
		c0 = -6.87
		break
	case Pluto:
		c0 = -1.00
		c1 = +4.00
		break
	default:
		return 0.0, errors.New("invalid body for visual magnitude")
	}
	x := phase / 100.0
	mag := c0 + x*(c1+x*(c2+x*c3))
	mag += 5.0 * math.Log10(helioDist*geoDist)
	return mag, nil
}

//--- Search code begins here ---------------------------------------------------------------------

func quadInterp(tm, dt, fa, fm, fb float64) (bool, float64, float64) {
	var Q, R, S, u, ru, x, x1, x2, outT, outDfDt float64

	Q = (fb+fa)/2.0 - fm
	R = (fb - fa) / 2.0
	S = fm

	if Q == 0.0 {
		// This is a line, not a parabola.
		if R == 0.0 {
			return false, 0.0, 0.0 // This is a HORIZONTAL line... can't make progress!
		}
		x = -S / R
		if x < -1.0 || x > +1.0 {
			return false, 0.0, 0.0 // out of bounds
		}
	} else {
		// This really is a parabola. Find roots x1, x2.
		u = R*R - 4*Q*S
		if u <= 0.0 {
			return false, 0.0, 0.0 // can't solve if imaginary, or if vertex of parabola is tangent.
		}

		ru = math.Sqrt(u)
		x1 = (-R + ru) / (2.0 * Q)
		x2 = (-R - ru) / (2.0 * Q)
		if -1.0 <= x1 && x1 <= +1.0 {
			if -1.0 <= x2 && x2 <= +1.0 {
				return false, 0.0, 0.0 // two roots are within bounds; we require a unique zero-crossing.
			}
			x = x1
		} else if -1.0 <= x2 && x2 <= +1.0 {
			x = x2
		} else {
			return false, 0.0, 0.0 // neither root is within bounds
		}
	}

	outT = tm + x*dt
	outDfDt = (2*Q*x + R) / dt
	return true, outT, outDfDt
}

//--- Search code ends here ---------------------------------------------------------------------

// Calculates one of the 5 Lagrange points from body masses and state vectors.
// Given a more massive "major" body and a much less massive "minor" body,
// calculates one of the five Lagrange points in relation to the minor body's
// orbit around the major body. The parameter `point` is an integer that
// selects the Lagrange point as follows:
//
// 1 = the Lagrange point between the major body and minor body.
// 2 = the Lagrange point on the far side of the minor body.
// 3 = the Lagrange point on the far side of the major body.
// 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
// 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
//
// The caller passes in the state vector and mass for both bodies.
// The state vectors can be in any orientation and frame of reference.
// The body masses are expressed as GM products, where G = the universal
// gravitation constant and M = the body's mass. Thus the units for
// majorMass and minorMass must be au^3/day^2.
// Use MassProduct to obtain GM values for various solar system bodies.
//
// The function returns the state vector for the selected Lagrange point
// using the same orientation as the state vector parameters majorState and minorState,
// and the position and velocity components are with respect to the major body's center.
//
// Consider calling LagrangePoint, instead of this function, for simpler usage in most cases.
func LagrangePointFast(point int, majorState StateVector, majorMass float64, minorState StateVector, minorMass float64) (*StateVector, error) {
	const cos60 = 0.5
	const sin60 = 0.8660254037844386 // sqrt(3) / 2

	if point < 1 || point > 5 {
		return nil, errors.New("Invalid lagrange point. Must be integer 1..5.")
	}

	if !isfinite(majorMass) || majorMass <= 0.0 {
		return nil, errors.New("Major mass must be a positive number.")
	}

	if !isfinite(minorMass) || minorMass <= 0.0 {
		return nil, errors.New("Minor mass must be a positive number.")
	}

	// Find the relative position vector <dx, dy, dz>.
	dx := minorState.X - majorState.X
	dy := minorState.Y - majorState.Y
	dz := minorState.Z - majorState.Z
	R2 := dx*dx + dy*dy + dz*dz

	// R = Total distance between the bodies.
	R := math.Sqrt(R2)

	// Find the velocity vector <vx, vy, vz>.
	vx := minorState.Vx - majorState.Vx
	vy := minorState.Vy - majorState.Vy
	vz := minorState.Vz - majorState.Vz

	var p StateVector
	if point == 4 || point == 5 {
		// For L4 and L5, we need to find points 60 degrees away from the
		// line connecting the two bodies and in the instantaneous orbital plane.
		// Define the instantaneous orbital plane as the unique plane that contains
		// both the relative position vector and the relative velocity vector.

		// Take the cross product of position and velocity to find a normal vector <nx, ny, nz>.
		nx := dy*vz - dz*vy
		ny := dz*vx - dx*vz
		nz := dx*vy - dy*vx

		// Take the cross product normal*position to get a tangential vector <ux, uy, uz>.
		ux := ny*dz - nz*dy
		uy := nz*dx - nx*dz
		uz := nx*dy - ny*dx

		// Convert the tangential direction vector to a unit vector.
		U := math.Sqrt(ux*ux + uy*uy + uz*uz)
		ux /= U
		uy /= U
		uz /= U

		// Convert the relative position vector into a unit vector.
		dx /= R
		dy /= R
		dz /= R

		// Now we have two perpendicular unit vectors in the orbital plane: 'd' and 'u'.

		// Create new unit vectors rotated (+/-)60 degrees from the radius/tangent directions.
		var vert float64
		if point == 4 {
			vert = +sin60
		} else {
			vert = -sin60
		}

		// Rotated radial vector
		Dx := cos60*dx + vert*ux
		Dy := cos60*dy + vert*uy
		Dz := cos60*dz + vert*uz

		// Rotated tangent vector
		Ux := cos60*ux - vert*dx
		Uy := cos60*uy - vert*dy
		Uz := cos60*uz - vert*dz

		// Calculate L4/L5 positions relative to the major body.
		p.X = R * Dx
		p.Y = R * Dy
		p.Z = R * Dz

		// Use dot products to find radial and tangential components of the relative velocity.
		vrad := vx*dx + vy*dy + vz*dz
		vtan := vx*ux + vy*uy + vz*uz

		// Calculate L4/L5 velocities.
		p.Vx = vrad*Dx + vtan*Ux
		p.Vy = vrad*Dy + vtan*Uy
		p.Vz = vrad*Dz + vtan*Uz
	} else {
		// Calculate the distances of each body from their mutual barycenter.
		// r1 = negative distance of major mass from barycenter (e.g. Sun to the left of barycenter)
		// r2 = positive distance of minor mass from barycenter (e.g. Earth to the right of barycenter)
		r1 := -R * (minorMass / (majorMass + minorMass))
		r2 := +R * (majorMass / (majorMass + minorMass))

		// Calculate the square of the angular orbital speed in [rad^2 / day^2].
		omega2 := (majorMass + minorMass) / (R2 * R)

		// Use Newton's Method to numerically solve for the location where
		// outward centrifugal acceleration in the rotating frame of reference
		// is equal to net inward gravitational acceleration.
		// First derive a good initial guess based on approximate analysis.
		var scale, numer1, numer2 float64
		if point == 1 || point == 2 {
			scale = (majorMass / (majorMass + minorMass)) * math.Cbrt(minorMass/(3.0*majorMass))
			numer1 = -majorMass // The major mass is to the left of L1 and L2
			if point == 1 {
				scale = 1.0 - scale
				numer2 = +minorMass // The minor mass is to the right of L1.
			} else {
				scale = 1.0 + scale
				numer2 = -minorMass // The minor mass is to the left of L2.
			}
		} else { // point == 3
			scale = ((7.0/12.0)*minorMass - majorMass) / (minorMass + majorMass)
			numer1 = +majorMass // major mass is to the right of L3.
			numer2 = +minorMass // minor mass is to the right of L3.
		}

		// Iterate Newton's Method until it converges.
		x := R*scale - r1
		var deltax float64
		for {
			dr1 := x - r1
			dr2 := x - r2
			accel := omega2*x + numer1/(dr1*dr1) + numer2/(dr2*dr2)
			deriv := omega2 - 2*numer1/(dr1*dr1*dr1) - 2*numer2/(dr2*dr2*dr2)
			deltax = accel / deriv
			x -= deltax
			if math.Abs(deltax/R) < 1.0e-14 {
				break
			}
		}

		scale = (x - r1) / R

		p.X = scale * dx
		p.Y = scale * dy
		p.Z = scale * dz
		p.Vx = scale * vx
		p.Vy = scale * vy
		p.Vz = scale * vz
	}
	p.T = majorState.T
	return &p, nil
}

//--- Generated code begins here ------------------------------------------------------------------

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

var jupiterMoonModel = [...]jupiterMoon{
	// [0] Io
	{
		mu:  2.8248942843381399e-07,
		al0: 1.4462132960212239e+00,
		al1: 3.5515522861824000e+00,
		a: []vsopTerm{
			{0.0028210960212903, 0.0000000000000000e+00, 0.0000000000000000e+00},
		},
		l: []vsopTerm{
			{-0.0001925258348666, 4.9369589722644998e+00, 1.3584836583050000e-02},
			{-0.0000970803596076, 4.3188796477322002e+00, 1.3034138432430000e-02},
			{-0.0000898817416500, 1.9080016428616999e+00, 3.0506486715799999e-03},
			{-0.0000553101050262, 1.4936156681568999e+00, 1.2938928911549999e-02},
		},
		z: []vsopTerm{
			{0.0041510849668155, 4.0899396355450000e+00, -1.2906864146660001e-02},
			{0.0006260521444113, 1.4461888986270000e+00, 3.5515522949801999e+00},
			{0.0000352747346169, 2.1256287034577999e+00, 1.2727416566999999e-04},
		},
		zeta: []vsopTerm{
			{0.0003142172466014, 2.7964219722923001e+00, -2.3150960980000000e-03},
			{0.0000904169207946, 1.0477061879627001e+00, -5.6920638196000003e-04},
		},
	},
	// [1] Europa
	{
		mu:  2.8248327439289299e-07,
		al0: -3.7352634374713622e-01,
		al1: 1.7693227111234699e+00,
		a: []vsopTerm{
			{0.0044871037804314, 0.0000000000000000e+00, 0.0000000000000000e+00},
			{0.0000004324367498, 1.8196456062910000e+00, 1.7822295777568000e+00},
		},
		l: []vsopTerm{
			{0.0008576433172936, 4.3188693178264002e+00, 1.3034138308049999e-02},
			{0.0004549582875086, 1.4936531751079001e+00, 1.2938928819619999e-02},
			{0.0003248939825174, 1.8196494533458001e+00, 1.7822295777568000e+00},
			{-0.0003074250079334, 4.9377037005910998e+00, 1.3584832867240000e-02},
			{0.0001982386144784, 1.9079869054759999e+00, 3.0510121286900001e-03},
			{0.0001834063551804, 2.1402853388529000e+00, 1.4500978933800000e-03},
			{-0.0001434383188452, 5.6222140366630002e+00, 8.9111478887838003e-01},
			{-0.0000771939140944, 4.3002724372349999e+00, 2.6733443704265998e+00},
		},
		z: []vsopTerm{
			{-0.0093589104136341, 4.0899396509038999e+00, -1.2906864146660001e-02},
			{0.0002988994545555, 5.9097265185595003e+00, 1.7693227079461999e+00},
			{0.0002139036390350, 2.1256289300016000e+00, 1.2727418406999999e-04},
			{0.0001980963564781, 2.7435168292649998e+00, 6.7797343008999997e-04},
			{0.0001210388158965, 5.5839943711203004e+00, 3.2056614899999997e-05},
			{0.0000837042048393, 1.6094538368039000e+00, -9.0402165808846002e-01},
			{0.0000823525166369, 1.4461887708689001e+00, 3.5515522949801999e+00},
		},
		zeta: []vsopTerm{
			{0.0040404917832303, 1.0477063169425000e+00, -5.6920640539999997e-04},
			{0.0002200421034564, 3.3368857864364001e+00, -1.2491307306999999e-04},
			{0.0001662544744719, 2.4134862374710999e+00, 0.0000000000000000e+00},
			{0.0000590282470983, 5.9719930968366004e+00, -3.0561602250000000e-05},
		},
	},
	// [2] Ganymede
	{
		mu:  2.8249818418472298e-07,
		al0: 2.8740893911433479e-01,
		al1: 8.7820792358932798e-01,
		a: []vsopTerm{
			{0.0071566594572575, 0.0000000000000000e+00, 0.0000000000000000e+00},
			{0.0000013930299110, 1.1586745884981000e+00, 2.6733443704265998e+00},
		},
		l: []vsopTerm{
			{0.0002310797886226, 2.1402987195941998e+00, 1.4500978438400001e-03},
			{-0.0001828635964118, 4.3188672736968003e+00, 1.3034138282630000e-02},
			{0.0001512378778204, 4.9373102372298003e+00, 1.3584834812520000e-02},
			{-0.0001163720969778, 4.3002659861490002e+00, 2.6733443704265998e+00},
			{-0.0000955478069846, 1.4936612842567001e+00, 1.2938928798570001e-02},
			{0.0000815246854464, 5.6222137132535002e+00, 8.9111478887838003e-01},
			{-0.0000801219679602, 1.2995922951532000e+00, 1.0034433456728999e+00},
			{-0.0000607017260182, 6.4978769669238001e-01, 5.0172167043264004e-01},
		},
		z: []vsopTerm{
			{0.0014289811307319, 2.1256295942738999e+00, 1.2727413029000001e-04},
			{0.0007710931226760, 5.5836330003496002e+00, 3.2064341100000001e-05},
			{0.0005925911780766, 4.0899396636447998e+00, -1.2906864146660001e-02},
			{0.0002045597496146, 5.2713683670371996e+00, -1.2523544076106000e-01},
			{0.0001785118648258, 2.8743156721063001e-01, 8.7820792442520001e-01},
			{0.0001131999784893, 1.4462127277818000e+00, 3.5515522949801999e+00},
			{-0.0000658778169210, 2.2702423990985001e+00, -1.7951364394536999e+00},
			{0.0000497058888328, 5.9096792204858000e+00, 1.7693227129285001e+00},
		},
		zeta: []vsopTerm{
			{0.0015932721570848, 3.3368862796665000e+00, -1.2491307058000000e-04},
			{0.0008533093128905, 2.4133881688166001e+00, 0.0000000000000000e+00},
			{0.0003513347911037, 5.9720789850126996e+00, -3.0561017709999999e-05},
			{-0.0001441929255483, 1.0477061764435001e+00, -5.6920632124000004e-04},
		},
	},
	// [3] Callisto
	{
		mu:  2.8249214488990899e-07,
		al0: -3.6203412913757038e-01,
		al1: 3.7648623343382798e-01,
		a: []vsopTerm{
			{0.0125879701715314, 0.0000000000000000e+00, 0.0000000000000000e+00},
			{0.0000035952049470, 6.4965776007116005e-01, 5.0172168165034003e-01},
			{0.0000027580210652, 1.8084235781510001e+00, 3.1750660413359002e+00},
		},
		l: []vsopTerm{
			{0.0005586040123824, 2.1404207189814999e+00, 1.4500979323100001e-03},
			{-0.0003805813868176, 2.7358844897852999e+00, 2.9729650620000000e-05},
			{0.0002205152863262, 6.4979652596399995e-01, 5.0172167243580001e-01},
			{0.0001877895151158, 1.8084787604004999e+00, 3.1750660413359002e+00},
			{0.0000766916975242, 6.2720114319754998e+00, 1.3928364636651001e+00},
			{0.0000747056855106, 1.2995916202344000e+00, 1.0034433456728999e+00},
		},
		z: []vsopTerm{
			{0.0073755808467977, 5.5836071576083999e+00, 3.2065099140000001e-05},
			{0.0002065924169942, 5.9209831565786004e+00, 3.7648624194703001e-01},
			{0.0001589869764021, 2.8744006242622999e-01, 8.7820792442520001e-01},
			{-0.0001561131605348, 2.1257397865089001e+00, 1.2727441285000001e-04},
			{0.0001486043380971, 1.4462134301023000e+00, 3.5515522949801999e+00},
			{0.0000635073108731, 5.9096803285953996e+00, 1.7693227129285001e+00},
			{0.0000599351698525, 4.1125517584797997e+00, -2.7985797954588998e+00},
			{0.0000540660842731, 5.5390350845569003e+00, 2.8683408228299999e-03},
			{-0.0000489596900866, 4.6218149483337996e+00, -6.2695712529518999e-01},
		},
		zeta: []vsopTerm{
			{0.0038422977898495, 2.4133922085556998e+00, 0.0000000000000000e+00},
			{0.0022453891791894, 5.9721736773277003e+00, -3.0561255249999997e-05},
			{-0.0002604479450559, 3.3368746306408998e+00, -1.2491309972000001e-04},
			{0.0000332112143230, 5.5604137742336999e+00, 2.9003768850700000e-03},
		},
	},
}

var rotationJupEqj = RotationMatrix{
	[3][3]float64{
		{9.99432765338654e-01, -3.36771074697641e-02, 0.00000000000000e+00},
		{3.03959428906285e-02, 9.02057912352809e-01, 4.30543388542295e-01},
		{-1.44994559663353e-02, -4.30299169409101e-01, 9.02569881273754e-01},
	},
}
