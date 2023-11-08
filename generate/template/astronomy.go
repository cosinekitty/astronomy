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
	MeanSynodicMonth          = 29.530588                     // average number of days for Moon to return to the same phase
	SpeedOfLightAuPerDay      = 173.1446326846693             // the speed of light in vacuum expressed in astronomical units per day
	KmPerAu                   = 1.4959787069098932e+8         // the number of kilometers in one astronomical unit
	AuPerLightYear            = 63241.07708807546             // the number of astronomical units in one light year
	AsecToRad                 = 4.848136811095359935899141e-6 // factor to convert arcseconds to radians
	SunRadiusKm               = 695700.0                      // the radius of the Sun in kilometers
	SunRadiusAu               = SunRadiusKm / KmPerAu
	MercuryEquatorialRadiusKm = 2440.5                                    // the equatorial radius of Mercury in kilometers
	MercuryPolarRadiusKm      = 2438.3                                    // the polar radius of Mercury in kilometers
	VenusRadiusKm             = 6051.8                                    // the radius of Venus in kilometers
	EarthEquatorialRadiusKm   = 6378.1366                                 // the equatorial radius of the Earth in kilometers
	EarthEquatorialRadiusAu   = EarthEquatorialRadiusKm / KmPerAu         //the equatorial radius of the Earth in astronomical units
	EarthFlattening           = 0.996647180302104                         // the ratio of the Earth's polar radius to its equatorial radius
	EarthPolarRadiusKm        = EarthEquatorialRadiusKm * EarthFlattening // the polar radius of the Earth in kilometers
	MoonEquatorialRadiusKm    = 1738.1                                    // the Moon's equatorial radius in kilometers
	MoonPolarRadiusKm         = 1736.0                                    // the Moon's polar radius in kilometers
	MoonMeanRadiusKm          = 1737.4
	MoonPolarRadiusAu         = MoonPolarRadiusKm / KmPerAu
	MarsEquatorialRadiusKm    = 3396.2  // the equatorial radius of Mars in kilometers
	MarsPolarRadiusKm         = 3376.2  // the polar radius of Mars in kilometers
	JupiterEquatorialRadiusKm = 71492.0 // the equatorial radius of Jupiter in kilometers
	JupiterPolarRadiusKm      = 66854.0 // the polar radius of Jupiter in kilometers
	JupiterMeanRadiusKm       = 69911.0 // the volumetric mean radius of Jupiter in kilometers
	IoRadiusKm                = 1821.6  // the radius of Jupiter's moon Io in kilometers
	EuropaRadiusKm            = 1560.8  // the radius of Jupiter's moon Europa in kilometers
	GanymedeRadiusKm          = 2631.2  // the radius of Jupiter's moon Ganymede in kilometers
	CallistoRadiusKm          = 2410.3  // the radius of Jupiter's moon Callisto in kilometers
	SaturnEquatorialRadiusKm  = 60268.0 // the equatorial radius of Saturn in kilometers
	SaturnPolarRadiusKm       = 54364.0 // the polar radius of Saturn in kilometers
	UranusEquatorialRadiusKm  = 25559.0 // the equatorial radius of Uranus in kilometers
	UranusPolarRadiusKm       = 24973.0 // the polar radius of Uranus in kilometers
	NeptuneEquatorialRadiusKm = 24764.0 // the equatorial radius of Neptune in kilometers
	NeptunePolarRadiusKm      = 24341.0 // the polar radius of Neptune in kilometers
	PlutoRadiusKm             = 1188.3  // the radius of Pluto in kilometers
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
	angVel               = 7.2921150e-5
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

func (a AstroVector) Add(b AstroVector) AstroVector {
	return AstroVector{
		a.X + b.X,
		a.Y + b.Y,
		a.Z + b.Z,
		a.T,
	}
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

type Direction int

const (
	Rise = +1
	Set  = -1
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

// Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
func RotationEqdHor(time AstroTime, observer Observer) RotationMatrix {
	sinlat := dsin(observer.Latitude)
	coslat := dcos(observer.Latitude)
	sinlon := dsin(observer.Longitude)
	coslon := dcos(observer.Longitude)

	uze := AstroVector{coslat * coslon, coslat * sinlon, sinlat, time}
	une := AstroVector{-sinlat * coslon, -sinlat * sinlon, coslat, time}
	uwe := AstroVector{sinlon, -coslon, 0.0, time}

	// Multiply sidereal hours by -15 to convert to degrees and flip eastward
	// rotation of the Earth to westward apparent movement of objects with time.
	angle := -15.0 * SiderealTime(&time)
	uz := spin(angle, uze)
	un := spin(angle, une)
	uw := spin(angle, uwe)

	r := RotationMatrix{}
	r.Rot[0][0] = un.X
	r.Rot[1][0] = un.Y
	r.Rot[2][0] = un.X
	r.Rot[0][1] = uw.X
	r.Rot[1][1] = uw.Y
	r.Rot[2][1] = uw.Z
	r.Rot[0][2] = uz.X
	r.Rot[1][2] = uz.Y
	r.Rot[2][2] = uz.Z

	return r
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

func synodicPeriod(body Body) float64 {
	// The Earth does not have a synodic period as seen from itself.
	if body == Earth {
		return 0.0
	}
	if body == Moon {
		return MeanSynodicMonth
	}
	tp := PlanetOrbitalPeriod(body)
	return math.Abs(earthOrbitalPeriod / (earthOrbitalPeriod/tp - 1.0))
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

func datan(x float64) float64 {
	return DegreesFromRadians(math.Atan(x))
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

func precessionPosVel(state StateVector, dir precessDirection) StateVector {
	rot := precessionRot(state.T, dir)
	return RotateState(rot, state)
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

func terra(observer Observer, time *AstroTime) StateVector {
	st := SiderealTime(time)
	phi := RadiansFromDegrees(observer.Latitude)
	sinphi := math.Sin(phi)
	cosphi := math.Cos(phi)
	c := 1.0 / math.Hypot(cosphi, EarthFlattening*sinphi)
	s := (EarthFlattening * EarthFlattening) * c
	htKm := observer.Height / 1000.0
	ach := EarthEquatorialRadiusKm*c + htKm
	ash := EarthEquatorialRadiusKm*s + htKm
	stlocl := RadiansFromDegrees(15.0*st + observer.Longitude)
	sinst := math.Sin(stlocl)
	cosst := math.Cos(stlocl)
	return StateVector{
		ach * cosphi * cosst / KmPerAu,
		ach * cosphi * sinst / KmPerAu,
		ash * sinphi / KmPerAu,
		-(angVel * 86400.0 / KmPerAu) * ach * cosphi * sinst,
		+(angVel * 86400.0 / KmPerAu) * ach * cosphi * cosst,
		0.0,
		*time,
	}
}

func inverseTerra(ovec AstroVector) Observer {
	var lonDeg, latDeg, heightKm float64

	// Convert from AU to kilometers.
	x := ovec.X * KmPerAu
	y := ovec.Y * KmPerAu
	z := ovec.Z * KmPerAu
	p := math.Hypot(x, y)
	if p < 1.0e-6 {
		// Special case: within 1 millimeter of a pole!
		// Use arbitrary longitude, and latitude determined by polarity of z.
		lonDeg = 0.0
		if z > 0.0 {
			latDeg = +90.0
		} else {
			latDeg = -90.0
		}
		// Elevation is calculated directly from z
		heightKm = math.Abs(z) - EarthPolarRadiusKm
	} else {
		stlocl := math.Atan2(y, x)
		st := SiderealTime(&ovec.T)
		// Calculate exact longitude.
		lonDeg = DegreesFromRadians(stlocl - (15.0 * st))
		// Normalize longitude to the range (-180, +180].
		for lonDeg <= -180.0 {
			lonDeg += 360.0
		}
		for lonDeg > +180.0 {
			lonDeg -= 360.0
		}
		// Numerically solve for exact latitude, using Newton's Method.
		F := EarthFlattening * EarthFlattening
		// Start with initial latitude estimate, based on a spherical Earth.
		lat := math.Atan2(z, p)
		var c, s, denom float64
		count := 0
		for {
			count++
			if count > 10 {
				panic("inverse_terra solver failed to converge.")
			}
			// Calculate the error function W(lat).
			// We try to find the root of W, meaning where the error is 0.
			c = math.Cos(lat)
			s = math.Sin(lat)
			factor := (F - 1) * EarthEquatorialRadiusKm
			c2 := c * c
			s2 := s * s
			radicand := c2 + F*s2
			denom = math.Sqrt(radicand)
			W := (factor*s*c)/denom - z*c + p*s
			if math.Abs(W) < 1.0e-8 {
				break // The error is now negligible.
			}
			// Error is still too large. Find the next estimate.
			// Calculate D = the derivative of W with respect to lat.
			D := factor*((c2-s2)/denom-s2*c2*(F-1)/(factor*radicand)) + z*s + p*c
			lat -= W / D
		}
		// We now have a solution for the latitude in radians.
		latDeg = DegreesFromRadians(lat)
		// Solve for exact height in kilometers.
		// There are two formulas I can use. Use whichever has the less risky denominator.
		adjust := EarthEquatorialRadiusKm / denom
		if math.Abs(s) > math.Abs(c) {
			heightKm = z/s - F*adjust
		} else {
			heightKm = p/c - adjust
		}
	}
	return Observer{latDeg, lonDeg, 1000.0 * heightKm}
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

func makeMoonContextFromTime(time AstroTime) *moonContext {
	return makeMoonContext(time.Tt / 36525.0)
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

func moonEclipticLatitudeDegrees(time AstroTime) float64 {
	context := makeMoonContextFromTime(time)
	moon := context.calcMoon()
	return DegreesFromRadians(moon.geoEclipLat)
}

// Calculates spherical ecliptic geocentric position of the Moon.
// Given a time of observation, calculates the Moon's geocentric position
// in ecliptic spherical coordinates. Provides the ecliptic latitude and
// longitude in degrees, and the geocentric distance in astronomical units (AU).
//
// The ecliptic angles are measured in "ECT": relative to the true ecliptic plane and
// equatorial plane at the specified time. This means the Earth's equator
// is corrected for precession and nutation, and the plane of the Earth's
// orbit is corrected for gradual obliquity drift.
//
// This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
// which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
// It is adapted from Turbo Pascal code from the book
// Astronomy on the Personal Computer by Montenbruck and Pfleger:
// https://www.springer.com/us/book/9783540672210
//
// To calculate a J2000 mean equator vector instead, use GeoMoon.
func EclipticGeoMoon(time AstroTime) Spherical {
	context := makeMoonContextFromTime(time)
	moon := context.calcMoon()

	// Convert spherical coordinates to a vector.
	// The moonResult angles are already expressed in radians.
	distCosLat := moon.distanceAu * math.Cos(moon.geoEclipLat)
	ecm := AstroVector{
		distCosLat * math.Cos(moon.geoEclipLon),
		distCosLat * math.Sin(moon.geoEclipLon),
		moon.distanceAu * math.Sin(moon.geoEclipLat),
		time,
	}

	// Obtain true and mean obliquity angles for the given time.
	// This serves to pre-calculate the nutation also, and cache it in `time`.
	et := etilt(&time)

	// Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
	eqm := eclToEquVec(ecm)

	// Add nutation to convert ECM to true equatorial coordinates of date (EQD).
	eqd := nutation(eqm, from2000)

	// Convert back to ecliptic, this time in true equinox of date (ECT).
	eclip := rotateEquatorialToEcliptic(eqd, RadiansFromDegrees(et.Tobl))

	return Spherical{eclip.Elat, eclip.Elon, moon.distanceAu}
}

// GeoMoon calculates the equatorial geocentric position of the Moon at a given time.
// The returned vector indicates the Moon's center relative to the Earth's center.
// The vector components are expressed in AU (astronomical units).
// The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
// In Astronomy Engine, this orientation is called EQJ.
func GeoMoon(time AstroTime) AstroVector {
	mc := makeMoonContextFromTime(time)
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

// Calculates the Moon's libration angles at a given moment in time.
// Libration is an observed back-and-forth wobble of the portion of the
// Moon visible from the Earth. It is caused by the imperfect tidal locking
// of the Moon's fixed rotation rate, compared to its variable angular speed
// of orbit around the Earth.
// This function calculates a pair of perpendicular libration angles,
// one representing rotation of the Moon in ecliptic longitude `Elon`, the other
// in ecliptic latitude `Elat`, both relative to the Moon's mean Earth-facing position.
// This function also returns the geocentric position of the Moon
// expressed in ecliptic longitude `Mlon`, ecliptic latitude `Mlat`, the
// distance `DistKm` between the centers of the Earth and Moon expressed in kilometers,
// and the apparent angular diameter of the Moon `DiamDeg`.
func Libration(time AstroTime) LibrationInfo {
	t := time.Tt / 36525.0
	t2 := t * t
	t3 := t2 * t
	t4 := t2 * t2

	context := makeMoonContext(t)
	moon := context.calcMoon()

	lib := LibrationInfo{}
	lib.Mlon = DegreesFromRadians(moon.geoEclipLon)
	lib.Mlat = DegreesFromRadians(moon.geoEclipLat)
	lib.DistKm = moon.distanceAu * KmPerAu
	lib.DiamDeg = 2.0 * datan(MoonMeanRadiusKm/math.Sqrt(lib.DistKm*lib.DistKm-MoonMeanRadiusKm*MoonMeanRadiusKm))

	// Inclination angle
	I := RadiansFromDegrees(1.543)

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

	// Optical librations
	w := moon.geoEclipLon - omega
	a := math.Atan2(math.Sin(w)*math.Cos(moon.geoEclipLat)*math.Cos(I)-math.Sin(moon.geoEclipLat)*math.Sin(I), math.Cos(w)*math.Cos(moon.geoEclipLat))
	ldash := longitudeOffset(DegreesFromRadians(a - f))
	bdash := math.Asin(-math.Sin(w)*math.Cos(moon.geoEclipLat)*math.Sin(I) - math.Sin(moon.geoEclipLat)*math.Cos(I))

	// Physical librations
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
	bdash = DegreesFromRadians(bdash)
	bdash2 := sigma*math.Cos(a) - rho*math.Sin(a)

	lib.Elon = ldash + ldash2
	lib.Elat = bdash + bdash2
	return lib
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

// Searches for a time at which a function's value increases through zero.
// Certain astronomy calculations involve finding a time when an event occurs.
// Often such events can be defined as the root of a function:
// the time at which the function's value becomes zero.
//
// Search finds the *ascending root* of a function: the time at which
// the function's value becomes zero while having a positive slope. That is, as time increases,
// the function transitions from a negative value, through zero at a specific moment,
// to a positive value later. The goal of the search is to find that specific moment.
//
// Search uses a combination of bisection and quadratic interpolation
// to minimize the number of function calls. However, it is critical that the
// supplied time window be small enough that there cannot be more than one root
// (ascending or descending) within it; otherwise the search can fail.
// Beyond that, it helps to make the time window as small as possible, ideally
// such that the function itself resembles a smooth parabolic curve within that window.
//
// If an ascending root is not found, or more than one root
// (ascending and/or descending) exists within the window `t1`..`t2`,
// the search will return `null`.
//
// If the search does not converge within 20 iterations, it will return an error.
func Search(context SearchContext, t1, t2 AstroTime, dtToleranceSeconds float64) (*AstroTime, error) {
	const iterLimit = 20
	dtDays := math.Abs(dtToleranceSeconds / SecondsPerDay)
	f1, e1 := context.Eval(t1)
	if e1 != nil {
		return nil, e1
	}
	f2, e2 := context.Eval(t2)
	if e2 != nil {
		return nil, e2
	}
	iter := 0
	calcFmid := true
	fmid := 0.0
	for {
		iter++
		if iter > iterLimit {
			return nil, errors.New("Search did not converge within 20 iterations")
		}

		dt := (t2.Tt - t1.Tt) / 2.0
		tmid := t1.AddDays(dt)
		if math.Abs(dt) < dtDays {
			// We are close enough to the event to stop the search.
			return &tmid, nil
		}

		if calcFmid {
			if f3, e3 := context.Eval(tmid); e3 != nil {
				return nil, e3
			} else {
				fmid = f3
			}
		} else {
			calcFmid = true // we already have the correct value of fmid from the previous loop
		}

		// Quadratic interpolation:
		// Try to find a parabola that passes through the 3 points we have sampled:
		// (t1,f1), (tmid,fmid), (t2,f2)

		if success, qUt, qDfDt := quadInterp(tmid.Ut, t2.Ut-tmid.Ut, f1, fmid, f2); success {
			tq := TimeFromUniversalDays(qUt)
			fq, eq := context.Eval(tq)
			if eq != nil {
				return nil, eq
			}
			if qDfDt != 0.0 {
				dtGuess := math.Abs(fq / qDfDt)
				if dtGuess < dtDays {
					// The estimated time error is small enough that we can quit now.
					return &tq, nil
				}

				// Try guessing a tighter boundary with the interpolated root at the center.
				dtGuess *= 1.2
				if dtGuess < dt/10.0 {
					tleft := tq.AddDays(-dtGuess)
					tright := tq.AddDays(+dtGuess)
					if (tleft.Ut-t1.Ut)*(tleft.Ut-t2.Ut) < 0.0 {
						if (tright.Ut-t1.Ut)*(tright.Ut-t2.Ut) < 0.0 {
							fleft, eleft := context.Eval(tleft)
							if eleft != nil {
								return nil, eleft
							}
							fright, eright := context.Eval(tright)
							if eright != nil {
								return nil, eright
							}
							if fleft < 0.0 && fright >= 0.0 {
								f1 = fleft
								f2 = fright
								t1 = tleft
								t2 = tright
								fmid = fq
								calcFmid = false // save a little work -- no need to re-calculate fmid next time around the loop
								continue
							}
						}
					}
				}
			}
		}

		// After quadratic interpolation attempt.
		// Now just divide the region in two parts and pick whichever one appears to contain a root.
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

		// Either there is no ascending zero-crossing in this range
		// or the search window is too wide (more than one zero-crossing).
		return nil, nil
	}
}

const riseSetDt = 0.42 // 10.08 hours: Nyquist-safe for 22-hour period.

type searchContextAltitude struct {
	body           Body
	direction      Direction
	observer       Observer
	bodyRadiusAu   float64
	targetAltitude float64
}

/*
func (context *searchContextAltitude) Eval(time AstroTime) (float64, error) {
	ofdate := Equator(context.body, time, context.observer, EquatorEpoch.OfDate, Aberration.Corrected)
}

func internalSearchAltitude(
	body Body,
	observer Observer,
	direction Direction,
	startTime AstroTime,
	limitDays, bodyRadiusAu, targetAltitude float64) (*AstroTime, error) {
	max_deriv_alt := maxAltitudeSlope(body, observer.Latitude)
	context := searchContextAltitude{body, direction, observer, bodyRadiusAu, targetAltitude}

	// We allow searching forward or backward in time.
	// But we want to keep t1 < t2, so we need a few if/else statements.
	t1 := startTime
	t2 := t1
	a1 := context.Eval(t1)
	a2 := a1

	for {
		if limitDays < 0.0 {
			t1 = t2.AddDays(-riseSetDt)
			a1 = context.Eval(t1)
		} else {
			t2 = t1.AddDays(+riseSetDt)
			a2 = context.Eval(t2)
		}

		ascent := findAscent(0, context, max_deriv_alt, t1, t2, a1, a2)
		if ascent.valid {
			// We found a time interval [t1, t2] that contains an alt-diff
			// rising from negative a1 to non-negative a2.
			// Search for the time where the root occurs.
			time, error := Search(context, ascent.tx, ascent.ty, 0.1)
			if time != nil {
				// Now that we have a solution, we have to check whether it goes outside the time bounds.
				if limitDays < 0.0 {
					if time.ut < startTime.ut+limitDays {
						return nil, nil
					}
				} else {
					if time.ut > startTime.ut+limitDays {
						return nil, nil
					}
				}
				return time, nil // success
			}

			if error != nil {
				return nil, error
			}

			// The search should have succeeded. Something is wrong with the ascent finder!
			panic("Rise/set search failed after finding an ascent.")
		}

		// There is no ascent in this interval, so keep searching.
		if limitDays < 0.0 {
			if t1.Ut < startTime.Ut+limitDays {
				return nil, nil
			}
			t2 = t1
			a2 = a1
		} else {
			if t2.Ut > startTime.Ut+limitDays {
				return nil, nil
			}
			t1 = t2
			a1 = a2
		}
	}
}
*/

func solarEclipseObscuration(hm, lo AstroVector) float64 {
	// hm = heliocentric moon
	// lo = lunacentric observer
	// Find heliocentric observer.
	ho := hm.Add(lo)

	// Calculate the apparent angular radius of the Sun for the observer.
	sunRadius := math.Asin(SunRadiusAu / ho.Length())

	// Calculate the apparent angular radius of the Moon for the observer.
	moonRadius := math.Asin(MoonPolarRadiusAu / lo.Length())

	// Calculate the apparent angular separation between the Sun's center and the Moon's center.
	sunMoonSeparation := AngleBetween(lo, ho)

	// Find the fraction of the Sun's apparent disc area that is covered by the Moon.
	obscuration := obscuration(sunRadius, moonRadius, RadiansFromDegrees(sunMoonSeparation))

	// HACK: In marginal cases, we need to clamp obscuration to less than 1.0.
	// This function is never called for total eclipses, so it should never return 1.0.
	return math.Min(0.9999, obscuration)
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
