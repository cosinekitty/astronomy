#!/usr/bin/env python3

import math
import datetime

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

def VectorLength(vector):
    return math.sqrt(vector.x**2 + vector.y**2 + vector.z**2)

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
    r = VectorLength(a) * VectorLength(b)
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

_DT = [
_delta_t_entry_t(-72638.0, 38),
_delta_t_entry_t(-65333.0, 26),
_delta_t_entry_t(-58028.0, 21),
_delta_t_entry_t(-50724.0, 21.1),
_delta_t_entry_t(-43419.0, 13.5),
_delta_t_entry_t(-39766.0, 13.7),
_delta_t_entry_t(-36114.0, 14.8),
_delta_t_entry_t(-32461.0, 15.7),
_delta_t_entry_t(-28809.0, 15.6),
_delta_t_entry_t(-25156.0, 13.3),
_delta_t_entry_t(-21504.0, 12.6),
_delta_t_entry_t(-17852.0, 11.2),
_delta_t_entry_t(-14200.0, 11.13),
_delta_t_entry_t(-10547.0, 7.95),
_delta_t_entry_t(-6895.0, 6.22),
_delta_t_entry_t(-3242.0, 6.55),
_delta_t_entry_t(-1416.0, 7.26),
_delta_t_entry_t(410.0, 7.35),
_delta_t_entry_t(2237.0, 5.92),
_delta_t_entry_t(4063.0, 1.04),
_delta_t_entry_t(5889.0, -3.19),
_delta_t_entry_t(7715.0, -5.36),
_delta_t_entry_t(9542.0, -5.74),
_delta_t_entry_t(11368.0, -5.86),
_delta_t_entry_t(13194.0, -6.41),
_delta_t_entry_t(15020.0, -2.70),
_delta_t_entry_t(16846.0, 3.92),
_delta_t_entry_t(18672.0, 10.38),
_delta_t_entry_t(20498.0, 17.19),
_delta_t_entry_t(22324.0, 21.41),
_delta_t_entry_t(24151.0, 23.63),
_delta_t_entry_t(25977.0, 24.02),
_delta_t_entry_t(27803.0, 23.91),
_delta_t_entry_t(29629.0, 24.35),
_delta_t_entry_t(31456.0, 26.76),
_delta_t_entry_t(33282.0, 29.15),
_delta_t_entry_t(35108.0, 31.07),
_delta_t_entry_t(36934.0, 33.150),
_delta_t_entry_t(38761.0, 35.738),
_delta_t_entry_t(40587.0, 40.182),
_delta_t_entry_t(42413.0, 45.477),
_delta_t_entry_t(44239.0, 50.540),
_delta_t_entry_t(44605.0, 51.3808),
_delta_t_entry_t(44970.0, 52.1668),
_delta_t_entry_t(45335.0, 52.9565),
_delta_t_entry_t(45700.0, 53.7882),
_delta_t_entry_t(46066.0, 54.3427),
_delta_t_entry_t(46431.0, 54.8712),
_delta_t_entry_t(46796.0, 55.3222),
_delta_t_entry_t(47161.0, 55.8197),
_delta_t_entry_t(47527.0, 56.3000),
_delta_t_entry_t(47892.0, 56.8553),
_delta_t_entry_t(48257.0, 57.5653),
_delta_t_entry_t(48622.0, 58.3092),
_delta_t_entry_t(48988.0, 59.1218),
_delta_t_entry_t(49353.0, 59.9845),
_delta_t_entry_t(49718.0, 60.7853),
_delta_t_entry_t(50083.0, 61.6287),
_delta_t_entry_t(50449.0, 62.2950),
_delta_t_entry_t(50814.0, 62.9659),
_delta_t_entry_t(51179.0, 63.4673),
_delta_t_entry_t(51544.0, 63.8285),
_delta_t_entry_t(51910.0, 64.0908),
_delta_t_entry_t(52275.0, 64.2998),
_delta_t_entry_t(52640.0, 64.4734),
_delta_t_entry_t(53005.0, 64.5736),
_delta_t_entry_t(53371.0, 64.6876),
_delta_t_entry_t(53736.0, 64.8452),
_delta_t_entry_t(54101.0, 65.1464),
_delta_t_entry_t(54466.0, 65.4573),
_delta_t_entry_t(54832.0, 65.7768),
_delta_t_entry_t(55197.0, 66.0699),
_delta_t_entry_t(55562.0, 66.3246),
_delta_t_entry_t(55927.0, 66.6030),
_delta_t_entry_t(56293.0, 66.9069),
_delta_t_entry_t(56658.0, 67.2810),
_delta_t_entry_t(57023.0, 67.6439),
_delta_t_entry_t(57388.0, 68.1024),
_delta_t_entry_t(57754.0, 68.5927),
_delta_t_entry_t(58119.0, 68.9676),
_delta_t_entry_t(58484.0, 69.2201),
_delta_t_entry_t(58849.0, 69.87),
_delta_t_entry_t(59214.0, 70.39),
_delta_t_entry_t(59580.0, 70.91),
_delta_t_entry_t(59945.0, 71.40),
_delta_t_entry_t(60310.0, 71.88),
_delta_t_entry_t(60675.0, 72.36),
_delta_t_entry_t(61041.0, 72.83),
_delta_t_entry_t(61406.0, 73.32),
_delta_t_entry_t(61680.0, 73.66)
]

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

class astro_time_t:
    def __init__(self, ut):
        self.ut = ut
        self.tt = _TerrestrialTime(ut)

    def AddDays(self, days):
        return astro_time_t(self.ut + days)

    def __str__(self):
        millis = round(self.ut * 86400000.0)
        n = _EPOCH + datetime.timedelta(milliseconds=millis)
        return '{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}Z'.format(n.year, n.month, n.day, n.hour, n.minute, n.second, math.floor(n.microsecond / 1000))

    def Utc(self):
        return _EPOCH + datetime.timedelta(days=self.ut)

_EPOCH = datetime.datetime(2000, 1, 1, 12)

def CurrentTime():
    ut = (datetime.datetime.utcnow() - _EPOCH).total_seconds() / 86400.0
    return astro_time_t(ut)

def MakeTime(year, month, day, hour, minute, second):
    micro = round((second % 1) * 1000000)
    second = math.floor(second - micro/1000000)
    d = datetime.datetime(year, month, day, hour, minute, second, micro)
    ut = (d - _EPOCH).total_seconds() / 86400
    return astro_time_t(ut)
