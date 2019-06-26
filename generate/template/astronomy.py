#!/usr/bin/env python3

import math
import datetime

_PI2 = 2.0 * math.pi
_EPOCH = datetime.datetime(2000, 1, 1, 12)
_T0 = 2451545.0
_MJD_BASIS = 2400000.5
_Y2000_IN_MJD = _T0 - _MJD_BASIS
_DEG2RAD = 0.017453292519943296
_RAD2DEG = 57.295779513082321
_ASEC360 = 1296000.0
_ASEC2RAD = 4.848136811095359935899141e-6
_ARC = 3600.0 * 180.0 / math.pi     # arcseconds per radian
_C_AUDAY = 173.1446326846693        # speed of light in AU/day
_ERAD = 6378136.6                   # mean earth radius in meters
_AU = 1.4959787069098932e+11        # astronomical unit in meters
_KM_PER_AU = 1.4959787069098932e+8
_ANGVEL = 7.2921150e-5
_SECONDS_PER_DAY = 24.0 * 3600.0
_SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592
_MEAN_SYNODIC_MONTH = 29.530588
_EARTH_ORBITAL_PERIOD = 365.256
_REFRACTION_NEAR_HORIZON = 34.0 / 60.0
_SUN_RADIUS_AU  = 4.6505e-3
_MOON_RADIUS_AU = 1.15717e-5
_ASEC180 = 180.0 * 60.0 * 60.0

def _LongitudeOffset(diff):
    offset = diff
    while offset <= -180.0:
        offset += 360.0
    while offset > 180.0:
        offset -= 360.0
    return offset

def _NormalizeLongitude(lon):
    while lon < 0.0:
        lon += 360.0
    while lon >= 360.0:
        lon -= 360.0
    return lon

class Vector:
    def __init__(self, x, y, z, t):
        self.x = x
        self.y = y
        self.z = z
        self.t = t

    def Length(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

BODY_INVALID = -1
BODY_MERCURY = 0
BODY_VENUS = 1
BODY_EARTH = 2
BODY_MARS = 3
BODY_JUPITER = 4
BODY_SATURN = 5
BODY_URANUS = 6
BODY_NEPTUNE = 7
BODY_PLUTO = 8
BODY_SUN = 9
BODY_MOON = 10

BodyName = [
    'Mercury',
    'Venus',
    'Earth',
    'Mars',
    'Jupiter',
    'Saturn',
    'Uranus',
    'Neptune',
    'Pluto',
    'Sun',
    'Moon',
]

def BodyCode(name):
    return BodyName.index(name)

def _IsSuperiorPlanet(body):
    return body in [BODY_MARS, BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO]

_PlanetOrbitalPeriod = [
    87.969,
    224.701,
    _EARTH_ORBITAL_PERIOD,
    686.980,
    4332.589,
    10759.22,
    30685.4,
    60189.0,
    90560.0
]

class Error(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

class EarthNotAllowedError(Error):
    def __init__(self):
        Error.__init__(self, 'The Earth is not allowed as the body.')

class InvalidBodyError(Error):
    def __init__(self):
        Error.__init__(self, 'Invalid astronomical body.')

class BadVectorError(Error):
    def __init__(self):
        Error.__init__(self, 'Vector is too small to have a direction.')

def _SynodicPeriod(body):
    if body == BODY_EARTH:
        raise EarthNotAllowedError()
    if body < 0 or body >= len(_PlanetOrbitalPeriod):
        raise InvalidBodyError()
    if body == BODY_MOON:
        return _MEAN_SYNODIC_MONTH
    return abs(_EARTH_ORBITAL_PERIOD / (_EARTH_ORBITAL_PERIOD/_PlanetOrbitalPeriod[body] - 1.0))

def _AngleBetween(a, b):
    r = a.Length() * b.Length()
    if r < 1.0e-8:
        return BadVectorError()
    dot = (a.x*b.x + a.y*b.y + a.z*b.z) / r
    if dot <= -1.0:
        return 180.0
    if dot >= +1.0:
        return 0.0
    return _RAD2DEG * math.acos(dot)

class _delta_t_entry_t:
    def __init__(self, mjd, dt):
        self.mjd = mjd
        self.dt = dt

_DT = $ASTRO_DELTA_T()

def _DeltaT(mjd):
    if mjd <= _DT[0].mjd:
        return _DT[0].dt
    if mjd >= _DT[-1].mjd:
        return _DT[-1].dt
    # Do a binary search to find the pair of indexes this mjd lies between.
    lo = 0
    hi = len(_DT) - 2   # Make sure there is always an array element after the one we are looking at.
    while True:
        if lo > hi:
            # This should never happen unless there is a bug in the binary search.
            raise Error('Could not find delta-t value.')
        c = (lo + hi) // 2
        if mjd < _DT[c].mjd:
            hi = c-1
        elif mjd > _DT[c+1].mjd:
            lo = c+1
        else:
            frac = (mjd - _DT[c].mjd) / (_DT[c+1].mjd - _DT[c].mjd)
            return _DT[c].dt + frac*(_DT[c+1].dt - _DT[c].dt)

def _TerrestrialTime(ut):
    return ut + _DeltaT(ut + _Y2000_IN_MJD) / 86400.0

class Time:
    def __init__(self, ut):
        self.ut = ut
        self.tt = _TerrestrialTime(ut)

    @staticmethod
    def Make(year, month, day, hour, minute, second):
        micro = round((second % 1) * 1000000)
        second = math.floor(second - micro/1000000)
        d = datetime.datetime(year, month, day, hour, minute, second, micro)
        ut = (d - _EPOCH).total_seconds() / 86400
        return Time(ut)

    @staticmethod
    def Now():
        ut = (datetime.datetime.utcnow() - _EPOCH).total_seconds() / 86400.0
        return Time(ut)

    def AddDays(self, days):
        return Time(self.ut + days)

    def __str__(self):
        millis = round(self.ut * 86400000.0)
        n = _EPOCH + datetime.timedelta(milliseconds=millis)
        return '{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}Z'.format(n.year, n.month, n.day, n.hour, n.minute, n.second, math.floor(n.microsecond / 1000))

    def Utc(self):
        return _EPOCH + datetime.timedelta(days=self.ut)


class Observer:
    def __init__(self, latitude, longitude, height=0):
        self.latitude = latitude
        self.longitude = longitude
        self.height = height

_nals_t = [
    [ 0,    0,    0,    0,    1],
    [ 0,    0,    2,   -2,    2],
    [ 0,    0,    2,    0,    2],
    [ 0,    0,    0,    0,    2],
    [ 0,    1,    0,    0,    0],
    [ 0,    1,    2,   -2,    2],
    [ 1,    0,    0,    0,    0],
    [ 0,    0,    2,    0,    1],
    [ 1,    0,    2,    0,    2],
    [ 0,   -1,    2,   -2,    2],
    [ 0,    0,    2,   -2,    1],
    [-1,    0,    2,    0,    2],
    [-1,    0,    0,    2,    0],
    [ 1,    0,    0,    0,    1],
    [-1,    0,    0,    0,    1],
    [-1,    0,    2,    2,    2],
    [ 1,    0,    2,    0,    1],
    [-2,    0,    2,    0,    1],
    [ 0,    0,    0,    2,    0],
    [ 0,    0,    2,    2,    2],
    [ 0,   -2,    2,   -2,    2],
    [-2,    0,    0,    2,    0],
    [ 2,    0,    2,    0,    2],
    [ 1,    0,    2,   -2,    2],
    [-1,    0,    2,    0,    1],
    [ 2,    0,    0,    0,    0],
    [ 0,    0,    2,    0,    0],
    [ 0,    1,    0,    0,    1],
    [-1,    0,    0,    2,    1],
    [ 0,    2,    2,   -2,    2],
    [ 0,    0,   -2,    2,    0],
    [ 1,    0,    0,   -2,    1],
    [ 0,   -1,    0,    0,    1],
    [-1,    0,    2,    2,    1],
    [ 0,    2,    0,    0,    0],
    [ 1,    0,    2,    2,    2],
    [-2,    0,    2,    0,    0],
    [ 0,    1,    2,    0,    2],
    [ 0,    0,    2,    2,    1],
    [ 0,   -1,    2,    0,    2],
    [ 0,    0,    0,    2,    1],
    [ 1,    0,    2,   -2,    1],
    [ 2,    0,    2,   -2,    2],
    [-2,    0,    0,    2,    1],
    [ 2,    0,    2,    0,    1],
    [ 0,   -1,    2,   -2,    1],
    [ 0,    0,    0,   -2,    1],
    [-1,   -1,    0,    2,    0],
    [ 2,    0,    0,   -2,    1],
    [ 1,    0,    0,    2,    0],
    [ 0,    1,    2,   -2,    1],
    [ 1,   -1,    0,    0,    0],
    [-2,    0,    2,    0,    2],
    [ 3,    0,    2,    0,    2],
    [ 0,   -1,    0,    2,    0],
    [ 1,   -1,    2,    0,    2],
    [ 0,    0,    0,    1,    0],
    [-1,   -1,    2,    2,    2],
    [-1,    0,    2,    0,    0],
    [ 0,   -1,    2,    2,    2],
    [-2,    0,    0,    0,    1],
    [ 1,    1,    2,    0,    2],
    [ 2,    0,    0,    0,    1],
    [-1,    1,    0,    1,    0],
    [ 1,    1,    0,    0,    0],
    [ 1,    0,    2,    0,    0],
    [-1,    0,    2,   -2,    1],
    [ 1,    0,    0,    0,    2],
    [-1,    0,    0,    1,    0],
    [ 0,    0,    2,    1,    2],
    [-1,    0,    2,    4,    2],
    [-1,    1,    0,    1,    1],
    [ 0,   -2,    2,   -2,    1],
    [ 1,    0,    2,    2,    1],
    [-2,    0,    2,    2,    2],
    [-1,    0,    0,    0,    2],
    [ 1,    1,    2,   -2,    2]
]

_cls_t = [
    [-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0],
    [ -13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0],
    [  -2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0],
    [   2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0],
    [   1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0],
    [   -516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0],
    [    711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0],
    [   -387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0],
    [   -301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0],
    [    215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0],
    [    128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0],
    [    123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0],
    [    156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0],
    [     63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0],
    [    -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0],
    [    -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0],
    [    -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0],
    [     45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0],
    [     63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0],
    [    -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0],
    [     32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0],
    [    -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0],
    [    -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0],
    [     28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0],
    [     20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0],
    [     29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0],
    [     25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0],
    [    -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0],
    [     15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0],
    [    -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0],
    [     21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0],
    [    -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0],
    [    -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0],
    [    -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0],
    [     16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0],
    [     -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0],
    [    -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0],
    [      7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0],
    [     -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0],
    [     -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0],
    [     -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0],
    [      5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0],
    [      6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0],
    [     -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0],
    [     -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0],
    [     -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0],
    [     -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0],
    [      7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0],
    [      4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0],
    [      6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0],
    [      3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0],
    [      4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0],
    [     -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0],
    [     -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0],
    [      4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0],
    [     -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0],
    [     -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0],
    [     -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0],
    [     -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0],
    [     -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0],
    [     -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0],
    [      2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0],
    [      2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0],
    [      3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0],
    [     -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0],
    [      3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0],
    [     -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0],
    [     -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0],
    [      4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0],
    [      1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0],
    [     -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0],
    [      1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0],
    [     -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0],
    [     -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0],
    [      1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0],
    [      1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0],
    [      1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0]
]

class _iau2000b:
    def __init__(self, time):
        t = time.tt / 36525
        el  = ((485868.249036 + t * 1717915923.2178) % ASEC360) * ASEC2RAD
        elp = ((1287104.79305 + t * 129596581.0481)  % ASEC360) * ASEC2RAD
        f   = ((335779.526232 + t * 1739527262.8478) % ASEC360) * ASEC2RAD
        d   = ((1072260.70369 + t * 1602961601.2090) % ASEC360) * ASEC2RAD
        om  = ((450160.398036 - t * 6962890.5431)    % ASEC360) * ASEC2RAD
        dp = 0
        de = 0
        i = 76
        while i >= 0:
            arg = (_nals_t[i][0]*el + _nals_t[i][1]*elp + _nals_t[i][2]*f + _nals_t[i][3]*d + _nals_t[i][4]*om) % _PI2
            sarg = math.sin(arg)
            carg = math.cos(arg)
            dp += (_cls_t[i][0] + _cls_t[i][1] * t)*sarg + _cls_t[i][2]*carg
            de += (_cls_t[i][3] + _cls_t[i][4] * t)*carg + _cls_t[i][5]*sarg
            i -= 1
        self.dpsi = -0.000135 + (dp * 1.0e-7)
        self.deps = +0.000388 + (de * 1.0e-7)

def _mean_obliq(tt):
    t = tt / 36525
    asec = (
        (((( -  0.0000000434   * t
             -  0.000000576  ) * t
             +  0.00200340   ) * t
             -  0.0001831    ) * t
             - 46.836769     ) * t + 84381.406
    )
    return asec / 3600.0

class _e_tilt:
    def __init__(self, time):
        e = _iau2000b(time)
        self.dpsi = e.dpsi
        self.deps = e.deps
        self.mobl = _mean_obliq(time.tt)
        self.tobl = self.mobl + (e.deps / 3600.0)
        self.tt = time.tt
        self.ee = e.dpsi * math.cos(self.mobl * _DEG2RAD) / 15.0

def _ecl2equ_vec(time, ecl):
    obl = _mean_obliq(time.tt) * _DEG2RAD
    cos_obl = math.cos(obl)
    sin_obl = math.sin(obl)
    return [
        ecl[0],
        ecl[1]*cos_obl - ecl[2]*sin_obl,
        ecl[1]*sin_obl + ecl[2]*cos_obl
    ]

def _precession(tt1, pos1, tt2):
    eps0 = 84381.406
    if tt1 != 0 and tt2 != 0:
        raise Error('One of (tt1, tt2) must be zero.')
    t = (tt2 - tt1) / 36525
    if tt2 == 0:
        t = -t

    psia  = (((((-    0.0000000951  * t
                 +    0.000132851 ) * t
                 -    0.00114045  ) * t
                 -    1.0790069   ) * t
                 + 5038.481507    ) * t)

    omegaa = (((((+   0.0000003337  * t
                 -    0.000000467 ) * t
                 -    0.00772503  ) * t
                 +    0.0512623   ) * t
                 -    0.025754    ) * t + eps0)

    chia  = (((((-    0.0000000560  * t
                 +    0.000170663 ) * t
                 -    0.00121197  ) * t
                 -    2.3814292   ) * t
                 +   10.556403    ) * t)

    eps0 *= _ASEC2RAD
    psia *= _ASEC2RAD
    omegaa *= _ASEC2RAD
    chia *= _ASEC2RAD

    sa = math.sin(eps0)
    ca = math.cos(eps0)
    sb = math.sin(-psia)
    cb = math.cos(-psia)
    sc = math.sin(-omegaa)
    cc = math.cos(-omegaa)
    sd = math.sin(chia)
    cd = math.cos(chia)

    xx =  cd * cb - sb * sd * cc
    yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc
    zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc
    xy = -sd * cb - sb * cd * cc
    yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc
    zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc
    xz =  sb * sc
    yz = -sc * cb * ca - sa * cc
    zz = -sc * cb * sa + cc * ca
    if tt2 == 0.0:
        # Perform rotation from other epoch to J2000.0.
        return [
            xx * pos1[0] + xy * pos1[1] + xz * pos1[2],
            yx * pos1[0] + yy * pos1[1] + yz * pos1[2],
            zx * pos1[0] + zy * pos1[1] + zz * pos1[2]
        ]

    # Perform rotation from J2000.0 to other epoch.
    return [
        xx * pos1[0] + yx * pos1[1] + zx * pos1[2],
        xy * pos1[0] + yy * pos1[1] + zy * pos1[2],
        xz * pos1[0] + yz * pos1[1] + zz * pos1[2]
    ]

class Equatorial:
    def __init__(self, ra, dec, dist):
        self.ra = ra
        self.dec = dec
        self.dist = dist

def _vector2radec(pos):
    xyproj = pos[0]*pos[0] + pos[1]*pos[1]
    dist = math.sqrt(xyproj + pos[2]*pos[2])
    if xyproj == 0.0:
        if pos[2] == 0.0:
            # Indeterminate coordinates: pos vector has zero length.
            raise Error('Cannot convert vector to polar coordinates')
        ra = 0.0
        if pos[2] < 0.0:
            dec = -90.0
        else:
            dec = +90.0
    else:
        ra = math.atan2(pos[1], pos[0]) / (_DEG2RAD * 15)
        if ra < 0:
            ra += 24
        dec = _RAD2DEG * math.atan2(pos[2], math.sqrt(xyproj))
    return Equatorial(ra, dec, dist)


def _nutation(time, direction, inpos):
    tilt = _e_tilt(time)
    oblm = tilt.mobl * _DEG2RAD
    oblt = tilt.tobl * _DEG2RAD
    psi = tilt.dpsi * _ASEC2RAD
    cobm = math.cos(oblm)
    sobm = math.sin(oblm)
    cobt = math.cos(oblt)
    sobt = math.sin(oblt)
    cpsi = math.cos(psi)
    spsi = math.sin(psi)

    xx = cpsi
    yx = -spsi * cobm
    zx = -spsi * sobm
    xy = spsi * cobt
    yy = cpsi * cobm * cobt + sobm * sobt
    zy = cpsi * sobm * cobt - cobm * sobt
    xz = spsi * sobt
    yz = cpsi * cobm * sobt - sobm * cobt
    zz = cpsi * sobm * sobt + cobm * cobt

    if direction == 0:
        # forward rotation
        return [
            xx * inpos[0] + yx * inpos[1] + zx * inpos[2],
            xy * inpos[0] + yy * inpos[1] + zy * inpos[2],
            xz * inpos[0] + yz * inpos[1] + zz * inpos[2]
        ]

    # inverse rotation
    return [
        xx * inpos[0] + xy * inpos[1] + xz * inpos[2],
        yx * inpos[0] + yy * inpos[1] + yz * inpos[2],
        zx * inpos[0] + zy * inpos[1] + zz * inpos[2]
    ]

def _era(time):        # Earth Rotation Angle
    thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut
    thet3 = time.ut % 1.0
    theta = 360.0 * ((thet1 + thet3) % 1.0)
    if theta < 0.0:
        theta += 360.0
    return theta


def _sidereal_time(time):
    t = time.tt / 36525.0
    eqeq = 15.0 * _e_tilt(time).ee    # Replace with eqeq=0 to get GMST instead of GAST (if we ever need it)
    theta = _era(time)
    st = (eqeq + 0.014506 +
        (((( -    0.0000000368   * t
            -    0.000029956  ) * t
            -    0.00000044   ) * t
            +    1.3915817    ) * t
            + 4612.156534     ) * t)
    gst = ((st/3600.0 + theta) % 360.0) / 15.0
    if gst < 0.0:
        gst += 24.0
    return gst


def _terra(observer, st):
    erad_km = _ERAD / 1000.0
    df = 1.0 - 0.003352819697896    # flattening of the Earth
    df2 = df * df
    phi = observer.latitude * _DEG2RAD
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    c = 1.0 / math.sqrt(cosphi*cosphi + df2*sinphi*sinphi)
    s = df2 * c
    ht_km = observer.height / 1000.0
    ach = erad_km*c + ht_km
    ash = erad_km*s + ht_km
    stlocl = (15.0*st + observer.longitude) * _DEG2RAD
    sinst = math.sin(stlocl)
    cosst = math.cos(stlocl)
    return [
        ach * cosphi * cosst / _KM_PER_AU,
        ach * cosphi * sinst / _KM_PER_AU,
        ash * sinphi / _KM_PER_AU
    ]

def _geo_pos(time, observer):
    gast = _sidereal_time(time)
    pos1 = _terra(observer, gast)
    pos2 = _nutation(time, -1, pos1)
    outpos = _precession(time.tt, pos2, 0.0)
    return outpos

def _spin(angle, pos1):
    angr = angle * _DEG2RAD
    cosang = math.cos(angr)
    sinang = math.sin(angr)
    return [
        +cosang*pos1[0] + sinang*pos1[1],
        -sinang*pos1[0] + cosang*pos1[1],
        pos1[2]
    ]

def _ter2cel(time, vec1):
    gast = _sidereal_time(time)
    return _spin(-15.0 * gast, vec1)

#----------------------------------------------------------------------------
# BEGIN CalcMoon

class _Array1:
    def __init__(self, xmin, xmax):
        self.min = xmin
        self.array = [0] * (xmax - xmin + 1)

    def __getitem__(self, key):
        return self.array[key - self.min]

    def __setitem__(self, key, value):
        self.array[key - self.min] = value

class _Array2:
    def __init__(self, xmin, xmax, ymin, ymax):
        self.min = xmin
        self.array = [_Array1(ymin, ymax) for i in range(xmax - xmin + 1)]

    def __getitem__(self, key):
        return self.array[key - self.min]

    def __setitem__(self, key, value):
        self.array[key - self.min] = value

class _moonpos:
    def __init__(self, lon, lat, dist):
        self.geo_eclip_lon = lon
        self.geo_eclip_lat = lat
        self.distance_au = dist

def _CalcMoon(time):
    T = time.tt / 36525
    co = _Array2(-6, 6, 1, 4)
    si = _Array2(-6, 6, 1, 4)

    def AddThe(c1, s1, c2, s2):
        return (c1*c2 - s1*s2, s1*c2 + c1*s2)

    def Sine(phi):
        return math.sin(_PI2 * phi)

    def Frac(x):
        return x - math.floor(x)

    T2 = T*T
    DLAM = 0
    DS = 0
    GAM1C = 0
    SINPI = 3422.7000
    S1 = Sine(0.19833+0.05611*T)
    S2 = Sine(0.27869+0.04508*T)
    S3 = Sine(0.16827-0.36903*T)
    S4 = Sine(0.34734-5.37261*T)
    S5 = Sine(0.10498-5.37899*T)
    S6 = Sine(0.42681-0.41855*T)
    S7 = Sine(0.14943-5.37511*T)
    DL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6
    DL  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6
    DLS =-6.40*S1                                   -1.89*S6
    DF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7
    DD  = DL0-DLS
    DGAM  = ((-3332E-9 * Sine(0.59734-5.37261*T)
               -539E-9 * Sine(0.35498-5.37899*T)
                -64E-9 * Sine(0.39943-5.37511*T)))

    L0 = _PI2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + DL0/_ARC
    L  = _PI2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + DL /_ARC
    LS = _PI2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + DLS/_ARC
    F  = _PI2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + DF /_ARC
    D  = _PI2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + DD /_ARC

    I = 1
    while I <= 4:
        if I == 1:
            ARG=L; MAX=4; FAC=1.000002208
        elif I == 2:
            ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T
        elif I == 3:
            ARG=F; MAX=4; FAC=1.000002708+139.978*DGAM
        else:
            ARG=D; MAX=6; FAC=1.0

        co[0][I] = 1
        co[1][I] = math.cos(ARG) * FAC
        si[0][I] = 0
        si[1][I] = math.sin(ARG) * FAC

        J = 2
        while J <= MAX:
            co[J][I], si[J][I] = AddThe(co[J-1][I], si[J-1][I], co[1][I], si[1][I])
            J += 1

        J = 1
        while J <= MAX:
            co[-J][I] = +co[J][I]
            si[-J][I] = -si[J][I]
            J += 1

        I += 1

    def Term(p, q, r, s):
        result = (1, 0)
        I = [None, p, q, r, s]
        k = 1
        while k <= 4:
            if I[k] != 0:
                result = AddThe(result[0], result[1], co[I[k]][k], si[I[k]][k])
            k += 1
        return result

    def AddSol(coeffl, coeffs, coeffg, coeffp, p, q, r, s):
        nonlocal DLAM, DS, GAM1C, SINPI
        result = Term(p, q, r, s)
        DLAM += coeffl * result[1]
        DS += coeffs * result[1]
        GAM1C += coeffg * result[0]
        SINPI += coeffp * result[0]

    AddSol(    13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4)
    AddSol(     0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3)
    AddSol(  2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2)
    AddSol(  -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1)
    AddSol(     1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4)
    AddSol(   191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2)
    AddSol(    -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1)
    AddSol( 22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0)
    AddSol(    18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1)
    AddSol( -4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2)
    AddSol(    +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3)
    AddSol(   -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4)
    AddSol(    -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6)
    AddSol(    -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4)
    AddSol(   -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2)
    AddSol(    18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1)
    AddSol(  -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0)
    AddSol(     0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1)
    AddSol(  -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2)
    AddSol(    -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4)
    AddSol(     0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4)
    AddSol(    14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2)
    AddSol(    -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1)
    AddSol(   769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0)
    AddSol(    +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1)
    AddSol(  -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2)
    AddSol(    +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3)
    AddSol(   -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4)
    AddSol(    -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6)
    AddSol(    -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2)
    AddSol(    +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1)
    AddSol(  -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0)
    AddSol(  -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2)
    AddSol(     0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3)
    AddSol(    -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4)

    AddSol(     0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4)
    AddSol(    14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2)
    AddSol(   147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0)
    AddSol(    -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1)
    AddSol(    28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2)
    AddSol(    -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3)
    AddSol(     0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4)
    AddSol(    -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2)
    AddSol(    -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0)
    AddSol(    -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2)
    AddSol(    -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2)
    AddSol(     0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1)
    AddSol(  -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0)
    AddSol(     0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1)
    AddSol(   -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2)
    AddSol(     0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3)
    AddSol(    +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4)
    AddSol(     1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2)
    AddSol(    36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0)
    AddSol(   -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2)
    AddSol(    -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4)
    AddSol(    -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6)
    AddSol(    -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2)
    AddSol(    -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0)
    AddSol(    -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2)
    AddSol(    -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4)
    AddSol(     1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2)
    AddSol(     9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0)
    AddSol(    -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1)
    AddSol(    -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2)
    AddSol(     0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4)
    AddSol(    -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0)
    AddSol(    -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2)
    AddSol(    -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4)
    AddSol(    +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2)
    AddSol(    +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0)
    AddSol(    +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2)
    AddSol(    -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2)
    AddSol(    -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2)
    AddSol(   -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0)
    AddSol(    -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2)
    AddSol(    -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4)
    AddSol(    -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2)
    AddSol(    39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0)
    AddSol(     9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2)
    AddSol(     0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4)

    AddSol(     0.415,    0.10, 0.0  ,  0.0013,0, 1, 2, 0)
    AddSol(    -2.152,   -2.26, 0.0  , -0.0066,0, 1, 2,-2)
    AddSol(    -1.440,   -1.30, 0.0  , +0.0014,0, 1,-2, 2)
    AddSol(     0.384,   -0.04, 0.0  ,  0.0   ,0, 1,-2,-2)
    AddSol(    +1.938,   +3.60,-0.145, +0.0401,4, 0, 0, 0)
    AddSol(    -0.952,   -1.58,+0.052, -0.0130,4, 0, 0,-2)
    AddSol(    -0.551,   -0.94,+0.032, -0.0097,3, 1, 0, 0)
    AddSol(    -0.482,   -0.57,+0.005, -0.0045,3, 1, 0,-2)
    AddSol(     0.681,    0.96,-0.026,  0.0115,3,-1, 0, 0)
    AddSol(    -0.297,   -0.27, 0.002, -0.0009,2, 2, 0,-2)
    AddSol(     0.254,   +0.21,-0.003,  0.0   ,2,-2, 0,-2)
    AddSol(    -0.250,   -0.22, 0.004,  0.0014,1, 3, 0,-2)
    AddSol(    -3.996,    0.0 , 0.0  , +0.0004,2, 0, 2, 0)
    AddSol(     0.557,   -0.75, 0.0  , -0.0090,2, 0, 2,-2)
    AddSol(    -0.459,   -0.38, 0.0  , -0.0053,2, 0,-2, 2)
    AddSol(    -1.298,    0.74, 0.0  , +0.0004,2, 0,-2, 0)
    AddSol(     0.538,    1.14, 0.0  , -0.0141,2, 0,-2,-2)
    AddSol(     0.263,    0.02, 0.0  ,  0.0   ,1, 1, 2, 0)
    AddSol(     0.426,   +0.07, 0.0  , -0.0006,1, 1,-2,-2)
    AddSol(    -0.304,   +0.03, 0.0  , +0.0003,1,-1, 2, 0)
    AddSol(    -0.372,   -0.19, 0.0  , -0.0027,1,-1,-2, 2)
    AddSol(    +0.418,    0.0 , 0.0  ,  0.0   ,0, 0, 4, 0)
    AddSol(    -0.330,   -0.04, 0.0  ,  0.0   ,3, 0, 2, 0)

    def ADDN(coeffn, p, q, r, s):
        return coeffn * Term(p, q, r, s)[1]

    N = 0
    N += ADDN(-526.069, 0, 0,1,-2)
    N += ADDN(  -3.352, 0, 0,1,-4)
    N += ADDN( +44.297,+1, 0,1,-2)
    N += ADDN(  -6.000,+1, 0,1,-4)
    N += ADDN( +20.599,-1, 0,1, 0)
    N += ADDN( -30.598,-1, 0,1,-2)
    N += ADDN( -24.649,-2, 0,1, 0)
    N += ADDN(  -2.000,-2, 0,1,-2)
    N += ADDN( -22.571, 0,+1,1,-2)
    N += ADDN( +10.985, 0,-1,1,-2)

    DLAM += (
        +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
        +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
        +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
        +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
        +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
        +0.33*Sine(0.3132   +6.3368*T)
    )
    S = F + DS/_ARC
    lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*math.sin(S) - 6.24*math.sin(3*S) + N
    return _moonpos(
        _PI2 * Frac((L0+DLAM/_ARC) / _PI2),
        (math.pi / (180 * 3600)) * lat_seconds,
        (_ARC * (_ERAD / _AU)) / (0.999953253 * SINPI)
    )

def GeoMoon(time):
    m = _CalcMoon(time)

    # Convert geocentric ecliptic spherical coordinates to Cartesian coordinates.
    dist_cos_lat = m.distance_au * math.cos(m.geo_eclip_lat)
    gepos = [
        dist_cos_lat * math.cos(m.geo_eclip_lon),
        dist_cos_lat * math.sin(m.geo_eclip_lon),
        m.distance_au * math.sin(m.geo_eclip_lat)
    ]

    # Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
    mpos1 = _ecl2equ_vec(time, gepos)

    # Convert from mean equinox of date to J2000.
    mpos2 = _precession(time.tt, mpos1, 0)
    return Vector(mpos2[0], mpos2[1], mpos2[2], time)

# END CalcMoon
#----------------------------------------------------------------------------
# BEGIN VSOP

_vsop = [
    # Mercury
    $ASTRO_LIST_VSOP(Mercury),

    # Venus
    $ASTRO_LIST_VSOP(Venus),

    # Earth
    $ASTRO_LIST_VSOP(Earth),

    # Mars
    $ASTRO_LIST_VSOP(Mars),

    # Jupiter
    $ASTRO_LIST_VSOP(Jupiter),

    # Saturn
    $ASTRO_LIST_VSOP(Saturn),

    # Uranus
    $ASTRO_LIST_VSOP(Uranus),

    # Neptune
    $ASTRO_LIST_VSOP(Neptune),
]

# END VSOP
#----------------------------------------------------------------------------
# BEGIN CHEBYSHEV

_pluto = $ASTRO_LIST_CHEBYSHEV(8)

# END CHEBYSHEV
#----------------------------------------------------------------------------
