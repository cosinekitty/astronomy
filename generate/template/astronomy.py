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

class astro_time_t:
    def __init__(self, ut):
        self.ut = ut
        self.tt = _TerrestrialTime(ut)

    def AddDays(self, days):
        return astro_time_t(self.ut + days)

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
