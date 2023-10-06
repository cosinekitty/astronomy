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

func TimeFromUniversalDays(ut float64) AstroTime {
	return makeTime(ut, terrestrialTime(ut))
}

func TimeFromTerrestrialDays(tt float64) AstroTime {
	return makeTime(universalTime(tt), tt)
}

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

type Spherical struct {
	Lat  float64
	Lon  float64
	Dist float64
}

type Body int

const (
	InvalidBody Body = -1
	Mercury          = 0
	Venus
	Earth
	Mars
	Jupiter
	Saturn
	Uranus
	Neptune
	Pluto
	Sun
	Moon
	Emb
	Ssb
	Star1 = 101
	Star2
	Star3
	Star4
	Star5
	Star6
	Star7
	Star8
)

type Observer struct {
	Latitude  float64
	Longitude float64
	Height    float64
}

type Equatorial struct {
	Ra   float64
	Dec  float64
	Dist float64
	Vec  AstroVector
}

type Ecliptic struct {
	Vec  AstroVector
	Elat float64
	Elon float64
}

type Topocentric struct {
	Azimuth  float64
	Altitude float64
	Ra       float64
	Dec      float64
}

type RotationMatrix struct {
	Rot [3][3]float64
}

type Refraction int

const (
	NoRefraction Refraction = iota
	NormalRefraction
	JplHorizonsRefraction
)

type AtmosphereInfo struct {
	Pressure    float64
	Temperature float64
	Density     float64
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

//$ASTRO_CONSTEL()
