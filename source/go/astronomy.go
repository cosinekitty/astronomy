package main



const C_AUDAY = 173.1446326846693
const KM_PER_AU = 1.4959787069098932e+08
const AU_PER_LY = 63241.07708807546
const DEG2RAD = 0.017453292519943295
const HOUR2RAD = 0.26179938779914946
const RAD2DEG = 57.29577951308232
const RAD2HOUR = 3.819718634205488
const SUN_RADIUS_KM = 695700.0
const MERCURY_EQUATORIAL_RADIUS_KM = 2440.5
const MERCURY_POLAR_RADIUS_KM = 2438.3
const VENUS_RADIUS_KM = 6051.8
const EARTH_EQUATORIAL_RADIUS_KM = 6378.1366
const EARTH_FLATTENING = 0.996647180302104
const EARTH_POLAR_RADIUS_KM = 0
const MOON_EQUATORIAL_RADIUS_KM = 1738.1
const MOON_POLAR_RADIUS_KM = 1736.0
const MARS_EQUATORIAL_RADIUS_KM = 3396.2
const MARS_POLAR_RADIUS_KM = 3376.2
const JUPITER_EQUATORIAL_RADIUS_KM = 71492.0
const JUPITER_POLAR_RADIUS_KM = 66854.0
const JUPITER_MEAN_RADIUS_KM = 69911.0
const IO_RADIUS_KM = 1821.6
const EUROPA_RADIUS_KM = 1560.8
const GANYMEDE_RADIUS_KM = 2631.2
const CALLISTO_RADIUS_KM = 2410.3
const SATURN_EQUATORIAL_RADIUS_KM = 60268.0
const SATURN_POLAR_RADIUS_KM = 54364.0
const URANUS_EQUATORIAL_RADIUS_KM = 25559.0
const URANUS_POLAR_RADIUS_KM = 24973.0
const NEPTUNE_EQUATORIAL_RADIUS_KM = 24764.0
const NEPTUNE_POLAR_RADIUS_KM = 24341.0
const PLUTO_RADIUS_KM = 1188.3
const TIME_TEXT_BYTES = 28

type astro_status_t int

const (
	ASTRO_SUCCESS = astro_status_t(iota)
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

type astro_time_t struct {
	Ut  float64
	Tt  float64
	Psi float64
	Eps float64
	St  float64
}
type astro_utc_t struct {
	Year   int
	Month  int
	Day    int
	Hour   int
	Minute int
	Second float64
}
type astro_vector_t struct {
	Status astro_status_t
	X      float64
	Y      float64
	Z      float64
	T      astro_time_t
}
type astro_state_vector_t struct {
	Status astro_status_t
	X      float64
	Y      float64
	Z      float64
	Vx     float64
	Vy     float64
	Vz     float64
	T      astro_time_t
}
type astro_spherical_t struct {
	Status astro_status_t
	Lat    float64
	Lon    float64
	Dist   float64
}
type astro_angle_result_t struct {
	Status astro_status_t
	Angle  float64
}
type astro_body_t int

const (
	BODY_INVALID astro_body_t = -1
	BODY_MERCURY astro_body_t = 0
	BODY_VENUS   astro_body_t = 1
	BODY_EARTH   astro_body_t = 2
	BODY_MARS    astro_body_t = 3
	BODY_JUPITER astro_body_t = 4
	BODY_SATURN  astro_body_t = 5
	BODY_URANUS  astro_body_t = 6
	BODY_NEPTUNE astro_body_t = 7
	BODY_PLUTO   astro_body_t = 8
	BODY_SUN     astro_body_t = 9
	BODY_MOON    astro_body_t = 10
	BODY_EMB     astro_body_t = 11
	BODY_SSB     astro_body_t = 12
	BODY_STAR1   astro_body_t = 101
	BODY_STAR2   astro_body_t = 102
	BODY_STAR3   astro_body_t = 103
	BODY_STAR4   astro_body_t = 104
	BODY_STAR5   astro_body_t = 105
	BODY_STAR6   astro_body_t = 106
	BODY_STAR7   astro_body_t = 107
	BODY_STAR8   astro_body_t = 108
)

type astro_observer_t struct {
	Latitude  float64
	Longitude float64
	Height    float64
}
type astro_equatorial_t struct {
	Status astro_status_t
	Ra     float64
	Dec    float64
	Dist   float64
	Vec    astro_vector_t
}
type astro_ecliptic_t struct {
	Status astro_status_t
	Vec    astro_vector_t
	Elat   float64
	Elon   float64
}
type astro_horizon_t struct {
	Azimuth  float64
	Altitude float64
	Ra       float64
	Dec      float64
}
type astro_rotation_t struct {
	Status astro_status_t
	Rot    [3][3]float64
}
type astro_refraction_t int

const (
	REFRACTION_NONE = astro_refraction_t(iota)
	REFRACTION_NORMAL
	REFRACTION_JPLHOR
)

type astro_atmosphere_t struct {
	Status      astro_status_t
	Pressure    float64
	Temperature float64
	Density     float64
}
type astro_search_result_t struct {
	Status astro_status_t
	Time   astro_time_t
}
type astro_seasons_t struct {
	Status       astro_status_t
	Mar_equinox  astro_time_t
	Jun_solstice astro_time_t
	Sep_equinox  astro_time_t
	Dec_solstice astro_time_t
}
type astro_moon_quarter_t struct {
	Status  astro_status_t
	Quarter int
	Time    astro_time_t
}
type astro_func_result_t struct {
	Status astro_status_t
	Value  float64
}
type astro_search_func_t func(context unsafe.Pointer, time astro_time_t) astro_func_result_t
type astro_deltat_func func(ut float64) float64
type astro_visibility_t int

const (
	VISIBLE_MORNING = astro_visibility_t(iota)
	VISIBLE_EVENING
)

type astro_elongation_t struct {
	Status              astro_status_t
	Time                astro_time_t
	Visibility          astro_visibility_t
	Elongation          float64
	Ecliptic_separation float64
}
type astro_hour_angle_t struct {
	Status astro_status_t
	Time   astro_time_t
	Hor    astro_horizon_t
}
type astro_illum_t struct {
	Status         astro_status_t
	Time           astro_time_t
	Mag            float64
	Phase_angle    float64
	Phase_fraction float64
	Helio_dist     float64
	Ring_tilt      float64
}
type astro_apsis_kind_t int

const (
	APSIS_PERICENTER = astro_apsis_kind_t(iota)
	APSIS_APOCENTER
	APSIS_INVALID
)

type astro_apsis_t struct {
	Status  astro_status_t
	Time    astro_time_t
	Kind    astro_apsis_kind_t
	Dist_au float64
	Dist_km float64
}
type astro_eclipse_kind_t int

const (
	ECLIPSE_NONE = astro_eclipse_kind_t(iota)
	ECLIPSE_PENUMBRAL
	ECLIPSE_PARTIAL
	ECLIPSE_ANNULAR
	ECLIPSE_TOTAL
)

type astro_lunar_eclipse_t struct {
	Status      astro_status_t
	Kind        astro_eclipse_kind_t
	Obscuration float64
	Peak        astro_time_t
	Sd_penum    float64
	Sd_partial  float64
	Sd_total    float64
}
type astro_global_solar_eclipse_t struct {
	Status      astro_status_t
	Kind        astro_eclipse_kind_t
	Obscuration float64
	Peak        astro_time_t
	Distance    float64
	Latitude    float64
	Longitude   float64
}
type astro_eclipse_event_t struct {
	Time     astro_time_t
	Altitude float64
}
type astro_local_solar_eclipse_t struct {
	Status        astro_status_t
	Kind          astro_eclipse_kind_t
	Obscuration   float64
	Partial_begin astro_eclipse_event_t
	Total_begin   astro_eclipse_event_t
	Peak          astro_eclipse_event_t
	Total_end     astro_eclipse_event_t
	Partial_end   astro_eclipse_event_t
}
type astro_transit_t struct {
	Status     astro_status_t
	Start      astro_time_t
	Peak       astro_time_t
	Finish     astro_time_t
	Separation float64
}
type astro_aberration_t int

const (
	ABERRATION = astro_aberration_t(iota)
	NO_ABERRATION
)

type astro_equator_date_t int

const (
	EQUATOR_J2000 = astro_equator_date_t(iota)
	EQUATOR_OF_DATE
)

type astro_direction_t int

const (
	DIRECTION_RISE astro_direction_t = +1
	DIRECTION_SET  astro_direction_t = -1
)

type astro_constellation_t struct {
	Status   astro_status_t
	Symbol   *byte
	Name     *byte
	Ra_1875  float64
	Dec_1875 float64
}
type astro_time_format_t int

const (
	TIME_FORMAT_DAY = astro_time_format_t(iota)
	TIME_FORMAT_MINUTE
	TIME_FORMAT_SECOND
	TIME_FORMAT_MILLI
)

type astro_libration_t struct {
	Elat     float64
	Elon     float64
	Mlat     float64
	Mlon     float64
	Dist_km  float64
	Diam_deg float64
}
type astro_axis_t struct {
	Status astro_status_t
	Ra     float64
	Dec    float64
	Spin   float64
	North  astro_vector_t
}
type astro_jupiter_moons_t struct {
	Io       astro_state_vector_t
	Europa   astro_state_vector_t
	Ganymede astro_state_vector_t
	Callisto astro_state_vector_t
}
type astro_node_kind_t int

const (
	INVALID_NODE    astro_node_kind_t = 0
	ASCENDING_NODE  astro_node_kind_t = +1
	DESCENDING_NODE astro_node_kind_t = -1
)

type astro_node_event_t struct {
	Status astro_status_t
	Time   astro_time_t
	Kind   astro_node_kind_t
}
type astro_grav_sim_s struct {
}
type astro_grav_sim_t astro_grav_sim_s
