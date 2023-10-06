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
