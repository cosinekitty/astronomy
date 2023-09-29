package main

import (
	"fmt"
)

const (
	C_AUDAY                      = 173.1446326846693
	KM_PER_AU                    = 1.4959787069098932e+8
	AU_PER_LY                    = 63241.07708807546
	DEG2RAD                      = 0.017453292519943296
	HOUR2RAD                     = 0.2617993877991494365
	RAD2DEG                      = 57.295779513082321
	RAD2HOUR                     = 3.819718634205488
	SUN_RADIUS_KM                = 695700.0
	MERCURY_EQUATORIAL_RADIUS_KM = 2440.5
	MERCURY_POLAR_RADIUS_KM      = 2438.3
	VENUS_RADIUS_KM              = 6051.8
	EARTH_EQUATORIAL_RADIUS_KM   = 6378.1366
	EARTH_FLATTENING             = 0.996647180302104
	EARTH_POLAR_RADIUS_KM        = (EARTH_EQUATORIAL_RADIUS_KM * EARTH_FLATTENING)
	MOON_EQUATORIAL_RADIUS_KM    = 1738.1
	MOON_POLAR_RADIUS_KM         = 1736.0
	MARS_EQUATORIAL_RADIUS_KM    = 3396.2
	MARS_POLAR_RADIUS_KM         = 3376.2
	JUPITER_EQUATORIAL_RADIUS_KM = 71492.0
	JUPITER_POLAR_RADIUS_KM      = 66854.0
	JUPITER_MEAN_RADIUS_KM       = 69911.0
	IO_RADIUS_KM                 = 1821.6
	EUROPA_RADIUS_KM             = 1560.8
	GANYMEDE_RADIUS_KM           = 2631.2
	CALLISTO_RADIUS_KM           = 2410.3
	SATURN_EQUATORIAL_RADIUS_KM  = 60268.0
	SATURN_POLAR_RADIUS_KM       = 54364.0
	URANUS_EQUATORIAL_RADIUS_KM  = 25559.0
	URANUS_POLAR_RADIUS_KM       = 24973.0
	NEPTUNE_EQUATORIAL_RADIUS_KM = 24764.0
	NEPTUNE_POLAR_RADIUS_KM      = 24341.0
	PLUTO_RADIUS_KM              = 1188.3
)

type astroStatus int

const (
	ASTRO_SUCCESS = iota
	ASTRO_NOT_INITIALIZED
	ASTRO_INVALID_BODY
	ASTRO_NO_CONVERGE
	ASTRO_BAD_TIME
	ASTRO_BAD_VECTOR
	ASTRO_SEARCH_FAILURE
	ASTRO_EARTH_NOT_ALLOWED
	ASTRO_NO_MOON_QUARTER
	ASTRO_WRONG_MOON_QUARTER
	ASTRO_INTERNAL_ERROR
	ASTRO_INVALID_PARAMETER
	ASTRO_FAIL_APSIS
	ASTRO_BUFFER_TOO_SMALL
	ASTRO_OUT_OF_MEMORY
	ASTRO_INCONSISTENT_TIMES
)

type astroTime struct {

	/**
	* @brief   UT1/UTC number of days since noon on January 1, 2000.
	*
	* The floating point number of days of Universal Time since noon UTC January 1, 2000.
	* Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
	* not exactly equivalent; UTC and UT1 can disagree by up to &plusmn;0.9 seconds.
	* This approximation is sufficient for the accuracy requirements of Astronomy Engine.
	*
	* Universal Time Coordinate (UTC) is the international standard for legal and civil
	* timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
	* UTC is kept in sync with unpredictable observed changes in the Earth's rotation
	* by occasionally adding leap seconds as needed.
	*
	* UT1 is an idealized time scale based on observed rotation of the Earth, which
	* gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
	* large scale weather events like hurricanes, and internal seismic and convection effects.
	* Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
	* is adjusted by a scheduled whole number of leap seconds as needed.
	*
	* The value in `ut` is appropriate for any calculation involving the Earth's rotation,
	* such as calculating rise/set times, culumination, and anything involving apparent
	* sidereal time.
	*
	* Before the era of atomic timekeeping, days based on the Earth's rotation
	* were often known as *mean solar days*.
	*/
	UT  float64
	/**
     *    Terrestrial Time days since noon on January 1, 2000.
     *
     * Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
     * In this system, days are not based on Earth rotations, but instead by
     * the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
     * divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
     * for changes in the Earth's rotation.
     *
     * The value in `tt` is used for calculations of movements not involving the Earth's rotation,
     * such as the orbits of planets around the Sun, or the Moon around the Earth.
     *
     * Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
     */
	TT  float64
	Psi float64 //For internal use only. Used to optimize Earth tilt calculations.
	Eps float64 //For internal use only.  Used to optimize Earth tilt calculations.
	St  float64 //For internal use only.  Lazy-caches sidereal time (Earth rotation).
}

type astroUTC struct {
	Year   int
	Month  int
	Day    int
	Hour   int
	Minute int
	Second float64
}

type astroVector struct {
	Status astroStatus
	X      float64
	Y      float64
	Z      float64
	T      astroTime
}

type astroStateVector struct {
	Status astroStatus
	X      float64
	Y      float64
	Z      float64
	Vx     float64
	Vy     float64
	Vz     float64
	T      astroTime
}

type astroSpherical struct {
	Status astroStatus
	Lat    float64
	Lon    float64
	Dist   float64
}

type astroAngleResult struct {
	Status astroStatus
	Angle  float64
}

type astroBody int

const (
	BODY_INVALID = -1
	BODY_MERCURY
	BODY_VENUS
	BODY_EARTH
	BODY_MARS
	BODY_JUPITER
	BODY_SATURN
	BODY_URANUS
	BODY_NEPTUNE
	BODY_PLUTO
	BODY_SUN
	BODY_MOON
	BODY_EMB
	BODY_SSB
	BODY_STAR1 = 101
	BODY_STAR2
	BODY_STAR3
	BODY_STAR4
	BODY_STAR5
	BODY_STAR6
	BODY_STAR7
	BODY_STAR8
)

type astroObserver struct {
	Latitude  float64
	Longitude float64
	Height    float64
}

type astroEquatorial struct {
	Status astroStatus
	RA     float64
	Dec    float64
	Dist   float64
	Vec    astroVector
}

type astroEcliptic struct {
	Status astroStatus
	Vec    astroVector
	Elat   float64
	Elon   float64
}

type astroHorizon struct {
	Azimuth  float64
	Altitude float64
	RA       float64
	Dec      float64
}

type astroRotation struct {
	Status astroStatus
	Rot    [3][3]float64
}

type astroRefraction int

const (
	REFRACTION_NONE = iota
	REFRACTION_NORMAL
	REFRACTION_JPLHOR
)

type astroAtmosphere struct {
	Status      astroStatus
	Pressure    float64
	Temperature float64
	Density     float64
}

type astroSearchResult struct {
	Status astroStatus
	Time   astroTime
}

type astroSeasons struct {
	Status      astroStatus
	MarEquinox  astroTime
	JunSolstice astroTime
	SepEquinox  astroTime
	DecSolstice astroTime
}

type astroMoonQuarter struct {
	Status  astroStatus
	Quarter int
	Time    astroTime
}

type astroFuncResult struct {
	Status astroStatus
	Value  float64
}

const (
	TimeFormatDay = iota
	TimeFormatMinute
	TimeFormatSecond
	TimeFormatMilli
)

type astroSearchFunc func(context interface{}, time astroTime) astroFuncResult

type astroDeltaTFunc func(ut float64) float64

type astroLibration struct {
	Elat    float64
	Elon    float64
	Mlat    float64
	Mlon    float64
	DistKm  float64
	DiamDeg float64
}

type astroAxis struct {
	Status astroStatus
	RA     float64
	Dec    float64
	Spin   float64
	North  astroVector
}

const TIME_TEXT_BYTES = 28

type astroNodeKind int

const (
	INVALID_NODE    = 0
	ASCENDING_NODE  = 1
	DESCENDING_NODE = -1
)

type astroNodeEvent struct {
	Status astroStatus
	Time   astroTime
	Kind   astroNodeKind
}

type TerseVector struct {
	X float64
	Y float64
	Z float64
}

// AstroTime represents a time in astronomy.

// BodyGravCalc represents gravitational calculations for a celestial body.
type BodyGravCalc struct {
	TT float64     // J2000 terrestrial time [days]
	R  TerseVector // Position [au]
	V  TerseVector // Velocity [au/day]
	A  TerseVector // Acceleration [au/day^2]
}

// GravSimEndpoint represents an endpoint in a gravitational simulation.

type astroJupiterMoons struct {
	Io       astroStateVector
	Europa   astroStateVector
	Ganymede astroStateVector
	Callisto astroStateVector
}
type GravSimEndpoint struct {
	Time        AstroTime
	Gravitators [1 + BODY_SUN]BodyGravCalc
	Bodies      *BodyGravCalc
}
type astroGravSim struct {
	OriginBody astroBody
	NumBodies  int
	Endpoint   [2]GravSimEndpoint
	Prev       *GravSimEndpoint
	Curr       *GravSimEndpoint
}

type astro_grav_sim_t = astroGravSim

type AstroStatus int

const (
	AstroSuccess AstroStatus = iota
	AstroError
)

type AstroFuncResult struct {
	Status AstroStatus // Status of the result.
	Value  float64     // The value returned by the function.
}

type AstroTime struct {
}

type AstroSearchResult struct {
}

type astroDeltatFunc func(ut float64) float64

var (
	daysPerTropicalYear                 = 365.24217
	deltaTFunc          astroDeltatFunc = AstronomyDeltaTEspenakMeeus
)

func AstronomySetDeltaTFunction(funcPtr astroDeltatFunc) {
	deltaTFunc = funcPtr
}

func TerrestrialTime(ut float64) float64 {
	return ut + deltaTFunc(ut)/86400.0
}

func AstronomyDeltaTEspenakMeeus(ut float64) float64 {

	/*
		Fred Espenak writes about Delta-T generically here:
		https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
		https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html

		He provides polynomial approximations for distant years here:
		https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html

		They start with a year value 'y' such that y=2000 corresponds
		to the UTC Date 15-January-2000. Convert difference in days
		to mean tropical years.
	*/

	var y, u, u2, u3, u4, u5, u6, u7 float64

	y = 2000 + ((ut - 14) / daysPerTropicalYear)

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

	/* all years after 2150 */
	u = (y - 1820) / 100
	return -20 + (32 * u * u)
}

func main() {
	d_TerrestrialTime := TerrestrialTime(235)
	fmt.Printf("Value of TerrestrialTime = %f\n", d_TerrestrialTime)
}

// test it -> go run astronomy.go
