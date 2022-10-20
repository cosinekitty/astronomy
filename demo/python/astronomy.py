#!/usr/bin/env python3
#
#    MIT License
#
#    Copyright (c) 2019-2022 Don Cross <cosinekitty@gmail.com>
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#
"""Astronomy Engine by Don Cross

See the GitHub project page for full documentation, examples,
and other information:

https://github.com/cosinekitty/astronomy

"""

import math
import datetime
import enum
import re
import abc

def _cbrt(x):
    if x < 0.0:
        return -((-x) ** (1.0 / 3.0))
    return x ** (1.0 / 3.0)

KM_PER_AU = 1.4959787069098932e+8   #<const> The number of kilometers per astronomical unit.
C_AUDAY   = 173.1446326846693       #<const> The speed of light expressed in astronomical units per day.

# Jupiter radius data are nominal values obtained from:
# https://www.iau.org/static/resolutions/IAU2015_English.pdf
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

JUPITER_EQUATORIAL_RADIUS_KM = 71492.0  #<const> The equatorial radius of Jupiter, expressed in kilometers.
JUPITER_POLAR_RADIUS_KM = 66854.0       #<const> The polar radius of Jupiter, expressed in kilometers.
JUPITER_MEAN_RADIUS_KM = 69911.0        #<const> The volumetric mean radius of Jupiter, expressed in kilometers.

# The radii of Jupiter's 4 largest moons were obtained from:
# https://ssd.jpl.nasa.gov/?sat_phys_par
IO_RADIUS_KM       = 1821.6     #<const> The mean radius of Jupiter's moon Io, expressed in kilometers.
EUROPA_RADIUS_KM   = 1560.8     #<const> The mean radius of Jupiter's moon Europa, expressed in kilometers.
GANYMEDE_RADIUS_KM = 2631.2     #<const> The mean radius of Jupiter's moon Ganymede, expressed in kilometers.
CALLISTO_RADIUS_KM = 2410.3     #<const> The mean radius of Jupiter's moon Callisto, expressed in kilometers.

_CalcMoonCount = 0

_RAD2HOUR  =  3.819718634205488         # 12/pi = factor to convert radians to sidereal hours
_HOUR2RAD  =  0.2617993877991494365     # pi/12 = factor to convert sidereal hours to radians
_DAYS_PER_TROPICAL_YEAR = 365.24217
_PI2 = 2.0 * math.pi
_EPOCH = datetime.datetime(2000, 1, 1, 12)
_ASEC360 = 1296000.0
_ASEC2RAD = 4.848136811095359935899141e-6
_ARC = 3600.0 * 180.0 / math.pi     # arcseconds per radian
_ANGVEL = 7.2921150e-5
_SECONDS_PER_DAY = 24.0 * 3600.0
_SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592
_MEAN_SYNODIC_MONTH = 29.530588
_EARTH_ORBITAL_PERIOD = 365.256
_NEPTUNE_ORBITAL_PERIOD = 60189.0
_REFRACTION_NEAR_HORIZON = 34.0 / 60.0

_SUN_RADIUS_KM = 695700.0
_SUN_RADIUS_AU  = _SUN_RADIUS_KM / KM_PER_AU

_EARTH_FLATTENING = 0.996647180302104
_EARTH_FLATTENING_SQUARED = _EARTH_FLATTENING ** 2
_EARTH_EQUATORIAL_RADIUS_KM = 6378.1366
_EARTH_POLAR_RADIUS_KM = _EARTH_EQUATORIAL_RADIUS_KM * _EARTH_FLATTENING
_EARTH_EQUATORIAL_RADIUS_AU = _EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU
_EARTH_MEAN_RADIUS_KM = 6371.0      # mean radius of the Earth's geoid, without atmosphere
_EARTH_ATMOSPHERE_KM = 88.0         # effective atmosphere thickness for lunar eclipses
_EARTH_ECLIPSE_RADIUS_KM = _EARTH_MEAN_RADIUS_KM + _EARTH_ATMOSPHERE_KM

_MOON_EQUATORIAL_RADIUS_KM = 1738.1
_MOON_EQUATORIAL_RADIUS_AU = (_MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU)
_MOON_MEAN_RADIUS_KM       = 1737.4
_MOON_POLAR_RADIUS_KM      = 1736.0
_MOON_POLAR_RADIUS_AU      = (_MOON_POLAR_RADIUS_KM / KM_PER_AU)

_ASEC180 = 180.0 * 60.0 * 60.0
_AU_PER_PARSEC = _ASEC180 / math.pi
_EARTH_MOON_MASS_RATIO = 81.30056

#
#    Masses of the Sun and outer planets, used for:
#    (1) Calculating the Solar System Barycenter
#    (2) Integrating the movement of Pluto
#
#    https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
#
#    Page 10 in the above document describes the constants used in the DE405 ephemeris.
#    The following are G*M values (gravity constant * mass) in [au^3 / day^2].
#    This side-steps issues of not knowing the exact values of G and masses M[i];
#    the products GM[i] are known extremely accurately.
#
_SUN_GM     = 0.2959122082855911e-03
_MERCURY_GM = 0.4912547451450812e-10
_VENUS_GM   = 0.7243452486162703e-09
_EARTH_GM   = 0.8887692390113509e-09
_MARS_GM    = 0.9549535105779258e-10
_JUPITER_GM = 0.2825345909524226e-06
_SATURN_GM  = 0.8459715185680659e-07
_URANUS_GM  = 0.1292024916781969e-07
_NEPTUNE_GM = 0.1524358900784276e-07
_PLUTO_GM   = 0.2188699765425970e-11

_MOON_GM = _EARTH_GM / _EARTH_MOON_MASS_RATIO


def MassProduct(body):
    """Returns the product of mass and universal gravitational constant of a Solar System body.

    For problems involving the gravitational interactions of Solar System bodies,
    it is helpful to know the product GM, where G = the universal gravitational constant
    and M = the mass of the body. In practice, GM is known to a higher precision than
    either G or M alone, and thus using the product results in the most accurate results.
    This function returns the product GM in the units au^3/day^2.
    The values come from page 10 of a
    [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).

    Parameters
    ----------
    body : Body
        The body for which to find the GM product.
        Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet.
        Any other value will cause an exception to be thrown.

    Returns
    -------
    float
        The mass product of the given body in au^3/day^2.
    """
    if body == Body.Sun:      return _SUN_GM
    if body == Body.Mercury:  return _MERCURY_GM
    if body == Body.Venus:    return _VENUS_GM
    if body == Body.Earth:    return _EARTH_GM
    if body == Body.Moon:     return _MOON_GM
    if body == Body.EMB:      return _EARTH_GM + _MOON_GM
    if body == Body.Mars:     return _MARS_GM
    if body == Body.Jupiter:  return _JUPITER_GM
    if body == Body.Saturn:   return _SATURN_GM
    if body == Body.Uranus:   return _URANUS_GM
    if body == Body.Neptune:  return _NEPTUNE_GM
    if body == Body.Pluto:    return _PLUTO_GM
    raise InvalidBodyError()

@enum.unique
class _PrecessDir(enum.Enum):
    From2000 = 0
    Into2000 = 1

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
    """A Cartesian vector with 3 space coordinates and 1 time coordinate.

    The vector's space coordinates are measured in astronomical units (AU).
    The coordinate system varies and depends on context.
    The vector also includes a time stamp.

    Attributes
    ----------
    x : float
        The x-coordinate of the vector, measured in AU.
    y : float
        The y-coordinate of the vector, measured in AU.
    z : float
        The z-coordinate of the vector, measured in AU.
    t : Time
        The date and time at which the coordinate is valid.
    """
    def __init__(self, x, y, z, t):
        self.x = x
        self.y = y
        self.z = z
        self.t = t

    def __repr__(self):
        return 'Vector({}, {}, {}, {})'.format(self.x, self.y, self.z, repr(self.t))

    def Length(self):
        """Returns the length of the vector in AU."""
        # It would be nice to use math.hypot() here,
        # but before Python 3.8, it only accepts 2 arguments.
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z, self.t)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z, self.t)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z, self.t)

    def format(self, coord_format):
        """Returns a custom format string representation of the vector."""
        layout = '({:' + coord_format + '}, {:' + coord_format + '}, {:' + coord_format + '}, {})'
        return layout.format(self.x, self.y, self.z, str(self.t))

class StateVector:
    """A combination of a position vector, a velocity vector, and a time.

    The position (x, y, z) is measured in astronomical units (AU).
    The velocity (vx, vy, vz) is measured in AU/day.
    The coordinate system varies and depends on context.
    The state vector also includes a time stamp.

    Attributes
    ----------
    x : float
        The x-coordinate of the position, measured in AU.
    y : float
        The y-coordinate of the position, measured in AU.
    z : float
        The z-coordinate of the position, measured in AU.
    vx : float
        The x-component of the velocity, measured in AU/day.
    vy : float
        The y-component of the velocity, measured in AU/day.
    vz : float
        The z-component of the velocity, measured in AU/day.
    t : Time
        The date and time at which the position and velocity vectors are valid.
    """
    def __init__(self, x, y, z, vx, vy, vz, t):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.t = t

    def __repr__(self):
        return 'StateVector(x={}, y={}, z={}, vx={}, vy={}, vz={}, t={})'.format(
            self.x, self.y, self.z,
            self.vx, self.vy, self.vz,
            repr(self.t))

    def __add__(self, other):
        return StateVector(
            self.x  + other.x,
            self.y  + other.y,
            self.z  + other.z,
            self.vx + other.vx,
            self.vy + other.vy,
            self.vz + other.vz,
            self.t
        )

    def __sub__(self, other):
        return StateVector(
            self.x  - other.x,
            self.y  - other.y,
            self.z  - other.z,
            self.vx - other.vx,
            self.vy - other.vy,
            self.vz - other.vz,
            self.t
        )

@enum.unique
class Body(enum.Enum):
    """The celestial bodies supported by Astronomy Engine calculations.

    Values
    ------
    Invalid: An unknown, invalid, or undefined celestial body.
    Mercury: The planet Mercury.
    Venus: The planet Venus.
    Earth: The planet Earth.
    Mars: The planet Mars.
    Jupiter: The planet Jupiter.
    Saturn: The planet Saturn.
    Uranus: The planet Uranus.
    Neptune: The planet Neptune.
    Pluto: The planet Pluto.
    Sun: The Sun.
    Moon: The Earth's moon.
    EMB: The Earth/Moon Barycenter.
    SSB: The Solar System Barycenter.
    """
    Invalid = -1
    Mercury = 0
    Venus = 1
    Earth = 2
    Mars = 3
    Jupiter = 4
    Saturn = 5
    Uranus = 6
    Neptune = 7
    Pluto = 8
    Sun = 9
    Moon = 10
    EMB = 11
    SSB = 12

def BodyCode(name):
    """Finds the Body enumeration value, given the name of a body.

    Parameters
    ----------
    name: str
        The common English name of a supported celestial body.

    Returns
    -------
    Body
        If `name` is a valid body name, returns the enumeration
        value associated with that body.
        Otherwise, returns `Body.Invalid`.

    Example
    -------

    >>> astronomy.BodyCode('Mars')
    <Body.Mars: 3>

    """
    if name not in Body.__members__:
        return Body.Invalid
    return Body[name]

def _IsSuperiorPlanet(body):
    return body in [Body.Mars, Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto]

_PlanetOrbitalPeriod = [
    87.969,
    224.701,
    _EARTH_ORBITAL_PERIOD,
    686.980,
    4332.589,
    10759.22,
    30685.4,
    _NEPTUNE_ORBITAL_PERIOD,
    90560.0
]

class Error(Exception):
    """Indicates an error in an astronomical calculation."""
    def __init__(self, message):
        Exception.__init__(self, message)

class DateTimeFormatError(Error):
    """The syntax of a UTC date/time string was not valid, or it contains invalid values."""
    def __init__(self, text):
        Error.__init__(self, 'The date/time string is not valid: "{}"'.format(text))

class EarthNotAllowedError(Error):
    """The Earth is not allowed as the celestial body in this calculation."""
    def __init__(self):
        Error.__init__(self, 'The Earth is not allowed as the body.')

class InvalidBodyError(Error):
    """The celestial body is not allowed for this calculation."""
    def __init__(self):
        Error.__init__(self, 'Invalid astronomical body.')

class BadVectorError(Error):
    """A vector magnitude is too small to have a direction in space."""
    def __init__(self):
        Error.__init__(self, 'Vector is too small to have a direction.')

class InternalError(Error):
    """An internal error occured that should be reported as a bug.

    Indicates an unexpected and unrecoverable condition occurred.
    If you encounter this error using Astronomy Engine, it would be very
    helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
    page on GitHub. Please include a copy of the stack trace, along with a description
    of how to reproduce the error. This will help improve the quality of
    Astronomy Engine for everyone! (Thank you in advance from the author.)
    """
    def __init__(self):
        Error.__init__(self, 'Internal error - please report issue, including stack trace, at https://github.com/cosinekitty/astronomy/issues')

class NoConvergeError(Error):
    """A numeric solver did not converge.

    Indicates that there was a failure of a numeric solver to converge.
    If you encounter this error using Astronomy Engine, it would be very
    helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
    page on GitHub. Please include a copy of the stack trace, along with a description
    of how to reproduce the error. This will help improve the quality of
    Astronomy Engine for everyone! (Thank you in advance from the author.)
    """
    def __init__(self):
        Error.__init__(self, 'Numeric solver did not converge - please report issue at https://github.com/cosinekitty/astronomy/issues')

def _SynodicPeriod(body):
    if body == Body.Earth:
        raise EarthNotAllowedError()
    if body.value < 0 or body.value >= len(_PlanetOrbitalPeriod):
        raise InvalidBodyError()
    if body == Body.Moon:
        return _MEAN_SYNODIC_MONTH
    return abs(_EARTH_ORBITAL_PERIOD / (_EARTH_ORBITAL_PERIOD/_PlanetOrbitalPeriod[body.value] - 1.0))

def AngleBetween(a, b):
    """Calculates the angle in degrees between two vectors.

    Given a pair of vectors, this function returns the angle in degrees
    between the two vectors in 3D space.
    The angle is measured in the plane that contains both vectors.

    Parameters
    ----------
    a : Vector
        The first of a pair of vectors between which to measure an angle.
    b : Vector
        The second of a pair of vectors between which to measure an angle.

    Returns
    -------
    float
        The angle between the two vectors expressed in degrees.
        The value is in the range [0, 180].
    """
    r = a.Length() * b.Length()
    if r < 1.0e-8:
        return BadVectorError()
    dot = (a.x*b.x + a.y*b.y + a.z*b.z) / r
    if dot <= -1.0:
        return 180.0
    if dot >= +1.0:
        return 0.0
    return math.degrees(math.acos(dot))


def DeltaT_EspenakMeeus(ut):
    """The default Delta T function used by Astronomy Engine.

    Espenak and Meeus use a series of piecewise polynomials to
    approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
    See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    This is the default Delta T function used by Astronomy Engine.

    Parameters
    ----------
    ut: float
        The floating point number of days since noon UTC on January 1, 2000.

    Returns
    -------
    float
        The estimated difference TT-UT on the given date, expressed in seconds.
    """
    # Fred Espenak writes about Delta-T generically here:
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html
    # He provides polynomial approximations for distant years here:
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    # They start with a year value 'y' such that y=2000 corresponds
    # to the UTC Date 15-January-2000. Convert difference in days
    # to mean tropical years.

    y = 2000 + ((ut - 14) / _DAYS_PER_TROPICAL_YEAR)

    if y < -500:
        u = (y - 1820) / 100
        return -20 + (32 * u*u)

    if y < 500:
        u = y / 100
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3
        return 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3 - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6

    if y < 1600:
        u = (y - 1000) / 100
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3
        return 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3 - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6

    if y < 1700:
        u = y - 1600
        u2 = u*u; u3 = u*u2
        return 120 - 0.9808*u - 0.01532*u2 + u3/7129.0

    if y < 1800:
        u = y - 1700
        u2 = u*u; u3 = u*u2; u4 = u2*u2
        return 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000

    if y < 1860:
        u = y - 1800
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3; u7 = u3*u4
        return 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4 + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7

    if y < 1900:
        u = y - 1860
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3
        return 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3 - 0.0004473624*u4 + u5/233174

    if y < 1920:
        u = y - 1900
        u2 = u*u; u3 = u*u2; u4 = u2*u2
        return -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4

    if y < 1941:
        u = y - 1920
        u2 = u*u; u3 = u*u2
        return 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3

    if y < 1961:
        u = y - 1950
        u2 = u*u; u3 = u*u2
        return 29.07 + 0.407*u - u2/233 + u3/2547

    if y < 1986:
        u = y - 1975
        u2 = u*u; u3 = u*u2
        return 45.45 + 1.067*u - u2/260 - u3/718

    if y < 2005:
        u = y - 2000
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3
        return 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4 + 0.00002373599*u5

    if y < 2050:
        u = y - 2000
        return 62.92 + 0.32217*u + 0.005589*u*u

    if y < 2150:
        u = (y-1820)/100
        return -20 + 32*u*u - 0.5628*(2150 - y)

    # all years after 2150
    u = (y - 1820) / 100
    return -20 + (32 * u*u)


_DeltaT = DeltaT_EspenakMeeus


def _TerrestrialTime(ut):
    return ut + _DeltaT(ut) / 86400.0

def _UniversalTime(tt):
    # This is the inverse function of _TerrestrialTime.
    # This is an iterative numerical solver, but because
    # the relationship between UT and TT is almost perfectly linear,
    # it converges extremely fast (never more than 3 iterations).
    dt = _TerrestrialTime(tt) - tt      # first approximation of dt = tt - ut
    while True:
        ut = tt - dt
        tt_check = _TerrestrialTime(ut)
        err = tt_check - tt
        if abs(err) < 1.0e-12:
            return ut
        dt += err

_TimeRegex = re.compile(r'^([\+\-]?[0-9]+)-([0-9]{2})-([0-9]{2})(T([0-9]{2}):([0-9]{2})(:([0-9]{2}(\.[0-9]+)?))?Z)?$')

class Time:
    """Represents a date and time used for performing astronomy calculations.

    All calculations performed by Astronomy Engine are based on
    dates and times represented by `Time` objects.

    Parameters
    ----------
    ut : float
        UT1/UTC number of days since noon on January 1, 2000.
        See the `ut` attribute of this class for more details.

    Attributes
    ----------
    ut : float
        The floating point number of days of Universal Time since noon UTC January 1, 2000.
        Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
        not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds.
        This approximation is sufficient for the accuracy requirements of Astronomy Engine.
        Universal Time Coordinate (UTC) is the international standard for legal and civil
        timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
        UTC is kept in sync with unpredictable observed changes in the Earth's rotation
        by occasionally adding leap seconds as needed.
        UT1 is an idealized time scale based on observed rotation of the Earth, which
        gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
        large scale weather events like hurricanes, and internal seismic and convection effects.
        Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
        is adjusted by a scheduled whole number of leap seconds as needed.
        The value in `ut` is appropriate for any calculation involving the Earth's rotation,
        such as calculating rise/set times, culumination, and anything involving apparent
        sidereal time.
        Before the era of atomic timekeeping, days based on the Earth's rotation
        were often known as *mean solar days*.
    tt : float
        Terrestrial Time days since noon on January 1, 2000.
        Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
        In this system, days are not based on Earth rotations, but instead by
        the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
        divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
        for changes in the Earth's rotation.
        The value in `tt` is used for calculations of movements not involving the Earth's rotation,
        such as the orbits of planets around the Sun, or the Moon around the Earth.
        Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
    """
    def __init__(self, ut, tt = None):
        if isinstance(ut, str):
            # Undocumented hack, to make repr(time) reversible.
            other = Time.Parse(ut)
            self.ut = other.ut
            self.tt = other.tt
        else:
            self.ut = ut
            if tt is None:
                self.tt = _TerrestrialTime(ut)
            else:
                self.tt = tt
        self._et = None     # lazy-cache for earth tilt
        self._st = None     # lazy-cache for sidereal time

    @staticmethod
    def FromTerrestrialTime(tt):
        """Creates a #Time object from a Terrestrial Time day value.

        Parameters
        ----------
        tt : float
            The number of days after the J2000 epoch.

        Returns
        -------
        Time
        """
        return Time(_UniversalTime(tt), tt)

    @staticmethod
    def Parse(text):
        """Creates a #Time object from a string of the form 'yyyy-mm-ddThh:mm:ss.sssZ'

        Parses a UTC date and time from a string and returns a #Time object.
        Permits a subset of ISO 8601 format.
        The year, month, and day are required.
        Hours, minutes, seconds, and fractions of a second are optional.
        If time is specified, there must be a 'T' between the date and the time
        and a 'Z' at the end of the time.

        Parameters
        ----------
        text : string
            A string of the following formats:
            `yyyy-mm-dd`
            `yyyy-mm-ddThh:mmZ`
            `yyyy-mm-ddThh:mm:ssZ`
            `yyyy-mm-ddThh:mm:ss.sssZ`

        Returns
        -------
        Time
        """
        m = _TimeRegex.match(text)
        if m is None:
            raise DateTimeFormatError(text)
        year = int(m.group(1))
        month = int(m.group(2))
        if not (1 <= month <= 12):
            raise DateTimeFormatError(text)
        day = int(m.group(3))
        if not (1 <= day <= 31):
            raise DateTimeFormatError(text)
        hour = int(m.group(5) or '0')
        if not (0 <= hour <= 23):
            raise DateTimeFormatError(text)
        minute = int(m.group(6) or '0')
        if not (0 <= minute <= 59):
            raise DateTimeFormatError(text)
        second = float(m.group(8) or '0')
        if not (0.0 <= second < 60.0):
            raise DateTimeFormatError(text)
        return Time.Make(year, month, day, hour, minute, second)

    @staticmethod
    def Make(year, month, day, hour, minute, second):
        """Creates a #Time object from a UTC calendar date and time.

        Parameters
        ----------
        year : int
            The UTC year value, e.g. 2019.
        month : int
            The UTC month in the range 1..12.
        day : int
            The UTC day of the month, in the range 1..31.
        hour : int
            The UTC hour, in the range 0..23.
        minute : int
            The UTC minute, in the range 0..59.
        second : float
            The real-valued UTC second, in the range [0, 60).

        Returns
        -------
        Time
        """
        # This formula is adapted from NOVAS C 3.1 function julian_date().
        y = int(year)
        m = int(month)
        d = int(day)
        y2000 = float(
            (d - 2483620)
            + 1461 * (y + 4800 - (14 - m) // 12) // 4
            + 367 * (m - 2 + (14 - m) // 12 * 12) // 12
            - 3 * ((y + 4900 - (14 - m) // 12) // 100) // 4
        )
        ut = (y2000 - 0.5) + (hour / 24.0) + (minute / 1440.0) + (second / 86400.0)
        return Time(ut)

    @staticmethod
    def Now():
        """Returns the computer's current date and time in the form of a #Time object.

        Uses the computer's system clock to find the current UTC date and time.
        Converts that date and time to a #Time value and returns the result.
        Callers can pass this value to other Astronomy Engine functions to
        calculate current observational conditions.

        Returns
        -------
        Time
        """
        ut = (datetime.datetime.utcnow() - _EPOCH).total_seconds() / 86400.0
        return Time(ut)

    def AddDays(self, days):
        """Calculates the sum or difference of a #Time with a specified real-valued number of days.

        Sometimes we need to adjust a given #Time value by a certain amount of time.
        This function adds the given real number of days in `days` to the date and time
        in the calling object.

        More precisely, the result's Universal Time field `ut` is exactly adjusted by `days`
        and the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
        using a best-fit piecewise polynomial model devised by
        [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).

        The value of the calling object is not modified. This function creates a brand new
        #Time object and returns it.

        Parameters
        ----------
        days : float
            A floating point number of days by which to adjust `time`.
            May be negative, 0, or positive.

        Returns
        -------
        Time
        """
        return Time(self.ut + days)

    def __repr__(self):
        return 'Time(\'' + str(self) + '\')'

    def __str__(self):
        # Adapted from the NOVAS C 3.1 function cal_date().
        djd = self.ut + 2451545.5
        jd = int(djd)
        x = 24.0 * math.fmod(djd, 1.0)
        hour = int(x)
        x = 60.0 * math.fmod(x, 1.0)
        minute = int(x)
        second = round(60.0 * math.fmod(x, 1.0), 3)
        if second >= 60.0:
            second -= 60.0
            minute += 1
            if minute >= 60:
                minute -= 60
                hour += 1
                if hour >= 24:
                    hour -= 24
                    jd += 1
        k = jd + 68569
        n = 4 * k // 146097
        k = k - (146097 * n + 3) // 4
        m = 4000 * (k + 1) // 1461001
        k = k - 1461 * m // 4 + 31
        month = (80 * k) // 2447
        day = k - 2447 * month // 80
        k = month // 11
        month = month + 2 - 12 * k
        year = 100 * (n - 49) + m + k
        millis = max(0, min(59999, round(1000.0 * second)))
        if year < 0:
            text = '-{:06d}'.format(-year)
        elif year <= 9999:
            text = '{:04d}'.format(year)
        else:
            text = '+{:06d}'.format(year)
        text += '-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}Z'.format(month, day, hour, minute, millis // 1000, millis % 1000)
        return text

    def Utc(self):
        """Returns the UTC date and time as a `datetime` object.

        Uses the standard [`datetime`](https://docs.python.org/3/library/datetime.html) class
        to represent the date and time in this Time object.

        Returns
        -------
        datetime
        """
        return _EPOCH + datetime.timedelta(days=self.ut)

    def _etilt(self):
        # Calculates precession and nutation of the Earth's axis.
        # The calculations are very expensive, so lazy-evaluate and cache
        # the result inside this Time object.
        if self._et is None:
            self._et = _e_tilt(self)
        return self._et

    def __lt__(self, other):
        return self.tt < other.tt

    def __eq__(self, other):
        return self.tt == other.tt

    def __le__(self, other):
        return self.tt <= other.tt

    def __ne__(self, other):
        return self.tt != other.tt

    def __gt__(self, other):
        return self.tt > other.tt

    def __ge__(self, other):
        return self.tt >= other.tt


class Observer:
    """Represents the geographic location of an observer on the surface of the Earth.

    Parameters
    ----------
    latitude : float
        Geographic latitude in degrees north of the equator.
    longitude : float
        Geographic longitude in degrees east of the prime meridian at Greenwich, England.
    height : float
        Elevation above sea level in meters.
    """
    def __init__(self, latitude, longitude, height=0.0):
        self.latitude = latitude
        self.longitude = longitude
        self.height = height

    def __repr__(self):
        return 'Observer(latitude={}, longitude={}, height={})'.format(self.latitude, self.longitude, self.height)

    def __str__(self):
        text = '('
        text += 'S' if (self.latitude < 0) else 'N'
        text += '{:0.8f}, '.format(abs(self.latitude))
        text += 'W' if (self.longitude < 0) else 'E'
        text += '{:0.8f}, '.format(abs(self.longitude))
        text += '{:0.3f}m'.format(self.height)
        text += ')'
        return text

class RotationMatrix:
    """Contains a rotation matrix that can be used to transform one
    coordinate system into another.

    Parameters
    ----------
    rot : float[3][3]
        A normalized 3x3 rotation matrix.
    """
    def __init__(self, rot):
        self.rot = rot

    def __repr__(self):
        return 'RotationMatrix({})'.format(self.rot)

class Spherical:
    """Holds spherical coordinates: latitude, longitude, distance.

    Parameters
    ----------
    lat : float
        The latitude angle: -90..+90 degrees.
    lon : float
        The longitude angle: 0..360 degrees.
    dist : float
        Distance in AU.
    """
    def __init__(self, lat, lon, dist):
        self.lat = lat
        self.lon = lon
        self.dist = dist

    def __repr__(self):
        return 'Spherical(lat={}, lon={}, dist={})'.format(self.lat, self.lon, self.dist)

class _iau2000b:
    def __init__(self, time):
        t = time.tt / 36525.0
        el  = math.fmod((485868.249036 + t*1717915923.2178), _ASEC360) * _ASEC2RAD
        elp = math.fmod((1287104.79305 + t*129596581.0481),  _ASEC360) * _ASEC2RAD
        f   = math.fmod((335779.526232 + t*1739527262.8478), _ASEC360) * _ASEC2RAD
        d   = math.fmod((1072260.70369 + t*1602961601.2090), _ASEC360) * _ASEC2RAD
        om  = math.fmod((450160.398036 - t*6962890.5431),    _ASEC360) * _ASEC2RAD
        dp = 0
        de = 0


        sarg = math.sin(om)
        carg = math.cos(om)
        dp += (-172064161.0 - 174666.0*t)*sarg + 33386.0*carg
        de += (92052331.0 + 9086.0*t)*carg + 15377.0*sarg


        arg = 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-13170906.0 - 1675.0*t)*sarg - 13696.0*carg
        de += (5730336.0 - 3015.0*t)*carg - 4587.0*sarg


        arg = 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-2276413.0 - 234.0*t)*sarg + 2796.0*carg
        de += (978459.0 - 485.0*t)*carg + 1374.0*sarg


        arg = 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (2074554.0 + 207.0*t)*sarg - 698.0*carg
        de += (-897492.0 + 470.0*t)*carg - 291.0*sarg


        sarg = math.sin(elp)
        carg = math.cos(elp)
        dp += (1475877.0 - 3633.0*t)*sarg + 11817.0*carg
        de += (73871.0 - 184.0*t)*carg - 1924.0*sarg


        arg = elp + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-516821.0 + 1226.0*t)*sarg - 524.0*carg
        de += (224386.0 - 677.0*t)*carg - 174.0*sarg


        sarg = math.sin(el)
        carg = math.cos(el)
        dp += (711159.0 + 73.0*t)*sarg - 872.0*carg
        de += (-6750.0)*carg + 358.0*sarg


        arg = 2.0*f + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-387298.0 - 367.0*t)*sarg + 380.0*carg
        de += (200728.0 + 18.0*t)*carg + 318.0*sarg


        arg = el + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-301461.0 - 36.0*t)*sarg + 816.0*carg
        de += (129025.0 - 63.0*t)*carg + 367.0*sarg


        arg = -elp + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (215829.0 - 494.0*t)*sarg + 111.0*carg
        de += (-95929.0 + 299.0*t)*carg + 132.0*sarg


        arg = 2.0*f - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (128227.0 + 137.0*t)*sarg + 181.0*carg
        de += (-68982.0 - 9.0*t)*carg + 39.0*sarg


        arg = -el + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (123457.0 + 11.0*t)*sarg + 19.0*carg
        de += (-53311.0 + 32.0*t)*carg - 4.0*sarg


        arg = -el + 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (156994.0 + 10.0*t)*sarg - 168.0*carg
        de += (-1235.0)*carg + 82.0*sarg


        arg = el + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (63110.0 + 63.0*t)*sarg + 27.0*carg
        de += (-33228.0)*carg - 9.0*sarg


        arg = -el + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-57976.0 - 63.0*t)*sarg - 189.0*carg
        de += (31429.0)*carg - 75.0*sarg


        arg = -el + 2.0*f + 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-59641.0 - 11.0*t)*sarg + 149.0*carg
        de += (25543.0 - 11.0*t)*carg + 66.0*sarg


        arg = el + 2.0*f + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-51613.0 - 42.0*t)*sarg + 129.0*carg
        de += (26366.0)*carg + 78.0*sarg


        arg = -2.0*el + 2.0*f + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (45893.0 + 50.0*t)*sarg + 31.0*carg
        de += (-24236.0 - 10.0*t)*carg + 20.0*sarg


        arg = 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (63384.0 + 11.0*t)*sarg - 150.0*carg
        de += (-1220.0)*carg + 29.0*sarg


        arg = 2.0*f + 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-38571.0 - 1.0*t)*sarg + 158.0*carg
        de += (16452.0 - 11.0*t)*carg + 68.0*sarg


        arg = -2.0*elp + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (32481.0)*sarg
        de += (-13870.0)*carg


        arg = -2.0*el + 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-47722.0)*sarg - 18.0*carg
        de += (477.0)*carg - 25.0*sarg


        arg = 2.0*el + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-31046.0 - 1.0*t)*sarg + 131.0*carg
        de += (13238.0 - 11.0*t)*carg + 59.0*sarg


        arg = el + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (28593.0)*sarg - carg
        de += (-12338.0 + 10.0*t)*carg - 3.0*sarg


        arg = -el + 2.0*f + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (20441.0 + 21.0*t)*sarg + 10.0*carg
        de += (-10758.0)*carg - 3.0*sarg


        arg = 2.0*el
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (29243.0)*sarg - 74.0*carg
        de += (-609.0)*carg + 13.0*sarg


        arg = 2.0*f
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (25887.0)*sarg - 66.0*carg
        de += (-550.0)*carg + 11.0*sarg


        arg = elp + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-14053.0 - 25.0*t)*sarg + 79.0*carg
        de += (8551.0 - 2.0*t)*carg - 45.0*sarg


        arg = -el + 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (15164.0 + 10.0*t)*sarg + 11.0*carg
        de += (-8001.0)*carg - sarg


        arg = 2.0*elp + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-15794.0 + 72.0*t)*sarg - 16.0*carg
        de += (6850.0 - 42.0*t)*carg - 5.0*sarg


        arg = -2.0*f + 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (21783.0)*sarg + 13.0*carg
        de += (-167.0)*carg + 13.0*sarg


        arg = el - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-12873.0 - 10.0*t)*sarg - 37.0*carg
        de += (6953.0)*carg - 14.0*sarg


        arg = -elp + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-12654.0 + 11.0*t)*sarg + 63.0*carg
        de += (6415.0)*carg + 26.0*sarg


        arg = -el + 2.0*f + 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-10204.0)*sarg + 25.0*carg
        de += (5222.0)*carg + 15.0*sarg


        arg = 2.0*elp
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (16707.0 - 85.0*t)*sarg - 10.0*carg
        de += (168.0 - 1.0*t)*carg + 10.0*sarg


        arg = el + 2.0*f + 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-7691.0)*sarg + 44.0*carg
        de += (3268.0)*carg + 19.0*sarg


        arg = -2.0*el + 2.0*f
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-11024.0)*sarg - 14.0*carg
        de += (104.0)*carg + 2.0*sarg


        arg = elp + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (7566.0 - 21.0*t)*sarg - 11.0*carg
        de += (-3250.0)*carg - 5.0*sarg


        arg = 2.0*f + 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-6637.0 - 11.0*t)*sarg + 25.0*carg
        de += (3353.0)*carg + 14.0*sarg


        arg = -elp + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-7141.0 + 21.0*t)*sarg + 8.0*carg
        de += (3070.0)*carg + 4.0*sarg


        arg = 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-6302.0 - 11.0*t)*sarg + 2.0*carg
        de += (3272.0)*carg + 4.0*sarg


        arg = el + 2.0*f - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (5800.0 + 10.0*t)*sarg + 2.0*carg
        de += (-3045.0)*carg - sarg


        arg = 2.0*el + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (6443.0)*sarg - 7.0*carg
        de += (-2768.0)*carg - 4.0*sarg


        arg = -2.0*el + 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-5774.0 - 11.0*t)*sarg - 15.0*carg
        de += (3041.0)*carg - 5.0*sarg


        arg = 2.0*el + 2.0*f + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-5350.0)*sarg + 21.0*carg
        de += (2695.0)*carg + 12.0*sarg


        arg = -elp + 2.0*f - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-4752.0 - 11.0*t)*sarg - 3.0*carg
        de += (2719.0)*carg - 3.0*sarg


        arg = -2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-4940.0 - 11.0*t)*sarg - 21.0*carg
        de += (2720.0)*carg - 9.0*sarg


        arg = -el - elp + 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (7350.0)*sarg - 8.0*carg
        de += (-51.0)*carg + 4.0*sarg


        arg = 2.0*el - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (4065.0)*sarg + 6.0*carg
        de += (-2206.0)*carg + sarg


        arg = el + 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (6579.0)*sarg - 24.0*carg
        de += (-199.0)*carg + 2.0*sarg


        arg = elp + 2.0*f - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (3579.0)*sarg + 5.0*carg
        de += (-1900.0)*carg + sarg


        arg = el - elp
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (4725.0)*sarg - 6.0*carg
        de += (-41.0)*carg + 3.0*sarg


        arg = -2.0*el + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-3075.0)*sarg - 2.0*carg
        de += (1313.0)*carg - sarg


        arg = 3.0*el + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-2904.0)*sarg + 15.0*carg
        de += (1233.0)*carg + 7.0*sarg


        arg = -elp + 2.0*d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (4348.0)*sarg - 10.0*carg
        de += (-81.0)*carg + 2.0*sarg


        arg = el - elp + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-2878.0)*sarg + 8.0*carg
        de += (1232.0)*carg + 4.0*sarg


        sarg = math.sin(d)
        carg = math.cos(d)
        dp += (-4230.0)*sarg + 5.0*carg
        de += (-20.0)*carg - 2.0*sarg


        arg = -el - elp + 2.0*f + 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-2819.0)*sarg + 7.0*carg
        de += (1207.0)*carg + 3.0*sarg


        arg = -el + 2.0*f
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-4056.0)*sarg + 5.0*carg
        de += (40.0)*carg - 2.0*sarg


        arg = -elp + 2.0*f + 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-2647.0)*sarg + 11.0*carg
        de += (1129.0)*carg + 5.0*sarg


        arg = -2.0*el + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-2294.0)*sarg - 10.0*carg
        de += (1266.0)*carg - 4.0*sarg


        arg = el + elp + 2.0*f + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (2481.0)*sarg - 7.0*carg
        de += (-1062.0)*carg - 3.0*sarg


        arg = 2.0*el + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (2179.0)*sarg - 2.0*carg
        de += (-1129.0)*carg - 2.0*sarg


        arg = -el + elp + d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (3276.0)*sarg + carg
        de += (-9.0)*carg


        arg = el + elp
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-3389.0)*sarg + 5.0*carg
        de += (35.0)*carg - 2.0*sarg


        arg = el + 2.0*f
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (3339.0)*sarg - 13.0*carg
        de += (-107.0)*carg + sarg


        arg = -el + 2.0*f - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-1987.0)*sarg - 6.0*carg
        de += (1073.0)*carg - 2.0*sarg


        arg = el + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-1981.0)*sarg
        de += (854.0)*carg


        arg = -el + d
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (4026.0)*sarg - 353.0*carg
        de += (-553.0)*carg - 139.0*sarg


        arg = 2.0*f + d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (1660.0)*sarg - 5.0*carg
        de += (-710.0)*carg - 2.0*sarg


        arg = -el + 2.0*f + 4.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-1521.0)*sarg + 9.0*carg
        de += (647.0)*carg + 4.0*sarg


        arg = -el + elp + d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (1314.0)*sarg
        de += (-700.0)*carg


        arg = -2.0*elp + 2.0*f - 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-1283.0)*sarg
        de += (672.0)*carg


        arg = el + 2.0*f + 2.0*d + om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (-1331.0)*sarg + 8.0*carg
        de += (663.0)*carg + 4.0*sarg


        arg = -2.0*el + 2.0*f + 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (1383.0)*sarg - 2.0*carg
        de += (-594.0)*carg - 2.0*sarg


        arg = -el + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (1405.0)*sarg + 4.0*carg
        de += (-610.0)*carg + 2.0*sarg


        arg = el + elp + 2.0*f - 2.0*d + 2.0*om
        sarg = math.sin(arg)
        carg = math.cos(arg)
        dp += (1290.0)*sarg
        de += (-556.0)*carg


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
        self.ee = e.dpsi * math.cos(math.radians(self.mobl)) / 15.0

def _ecl2equ_vec(time, ecl):
    obl = math.radians(_mean_obliq(time.tt))
    cos_obl = math.cos(obl)
    sin_obl = math.sin(obl)
    return [
        ecl[0],
        ecl[1]*cos_obl - ecl[2]*sin_obl,
        ecl[1]*sin_obl + ecl[2]*cos_obl
    ]

def _precession_rot(time, direction):
    eps0 = 84381.406
    t = time.tt / 36525

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
    if direction == _PrecessDir.Into2000:
        # Perform rotation from other epoch to J2000.0.
        return RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ])

    if direction == _PrecessDir.From2000:
        # Perform rotation from J2000.0 to other epoch.
        return RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ])

    raise Error('Inalid precession direction')

def _rotate(rot, vec):
    return [
        rot.rot[0][0]*vec[0] + rot.rot[1][0]*vec[1] + rot.rot[2][0]*vec[2],
        rot.rot[0][1]*vec[0] + rot.rot[1][1]*vec[1] + rot.rot[2][1]*vec[2],
        rot.rot[0][2]*vec[0] + rot.rot[1][2]*vec[1] + rot.rot[2][2]*vec[2]
    ]

def _precession(pos, time, direction):
    r = _precession_rot(time, direction)
    return _rotate(r, pos)

def _precession_posvel(state, time, direction):
    r = _precession_rot(time, direction)
    return RotateState(r, state)

class Equatorial:
    """Equatorial angular coordinates

    Coordinates of a celestial body as seen from the Earth.
    Can be geocentric or topocentric, depending on context.
    The coordinates are oriented with respect to the Earth's
    equator projected onto the sky.

    Attributes
    ----------
    ra : float
        Right ascension in sidereal hours.
    dec : float
        Declination in degrees.
    dist : float
        Distance to the celestial body in AU.
    vec : Vector
        The equatorial coordinates in cartesian form, using AU distance units.
        x = direction of the March equinox,
        y = direction of the June solstice,
        z = north.
    """
    def __init__(self, ra, dec, dist, vec):
        self.ra = ra
        self.dec = dec
        self.dist = dist
        self.vec = vec

    def __repr__(self):
        return 'Equatorial(ra={}, dec={}, dist={}, vec={})'.format(self.ra, self.dec, self.dist, repr(self.vec))


def _vector2radec(pos, time):
    xyproj = pos[0]*pos[0] + pos[1]*pos[1]
    dist = math.sqrt(xyproj + pos[2]*pos[2])
    if xyproj == 0.0:
        if pos[2] == 0.0:
            # Indeterminate coordinates: pos vector has zero length.
            raise BadVectorError()
        ra = 0.0
        if pos[2] < 0.0:
            dec = -90.0
        else:
            dec = +90.0
    else:
        ra = _RAD2HOUR * math.atan2(pos[1], pos[0])
        if ra < 0:
            ra += 24
        dec = math.degrees(math.atan2(pos[2], math.sqrt(xyproj)))
    vec = Vector(pos[0], pos[1], pos[2], time)
    return Equatorial(ra, dec, dist, vec)


def _nutation_rot(time, direction):
    tilt = time._etilt()
    oblm = math.radians(tilt.mobl)
    oblt = math.radians(tilt.tobl)
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

    if direction == _PrecessDir.From2000:
        # convert J2000 to of-date
        return RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ])

    if direction == _PrecessDir.Into2000:
        # convert of-date to J2000
        return RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ])

    raise Error('Invalid nutation direction')

def _nutation(pos, time, direction):
    r = _nutation_rot(time, direction)
    return _rotate(r, pos)

def _nutation_posvel(state, time, direction):
    r = _nutation_rot(time, direction)
    return RotateState(r, state)

def _era(time):        # Earth Rotation Angle
    thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut
    thet3 = math.fmod(time.ut, 1.0)
    theta = 360.0 * math.fmod((thet1 + thet3), 1.0)
    if theta < 0.0:
        theta += 360.0
    return theta

def SiderealTime(time):
    """Calculates Greenwich Apparent Sidereal Time (GAST).

    Given a date and time, this function calculates the rotation of the
    Earth, represented by the equatorial angle of the Greenwich prime meridian
    with respect to distant stars (not the Sun, which moves relative to background
    stars by almost one degree per day).
    This angle is called Greenwich Apparent Sidereal Time (GAST).
    GAST is measured in sidereal hours in the half-open range [0, 24).
    When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
    corrected at that time for precession and nutation of the Earth's axis.
    In this context, the "equinox" is the direction in space where the Earth's
    orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
    at the location on the Earth's orbit of the (seasonal) March equinox.
    As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
    then starts over at 0.
    To convert to degrees, multiply the return value by 15.

    Parameters
    ----------
    time : Time
        The date and time for which to find GAST.
        As an optimization, this function caches the sideral time value in `time`,
        unless it has already been cached, in which case the cached value is reused.

    Returns
    -------
    float
        GAST expressed in sidereal hours.
    """
    if time._st is None:
        t = time.tt / 36525.0
        eqeq = 15.0 * time._etilt().ee    # Replace with eqeq=0 to get GMST instead of GAST (if we ever need it)
        theta = _era(time)
        st = (eqeq + 0.014506 +
            (((( -    0.0000000368   * t
                -    0.000029956  ) * t
                -    0.00000044   ) * t
                +    1.3915817    ) * t
                + 4612.156534     ) * t)
        gst = math.fmod((st/3600.0 + theta), 360.0) / 15.0
        if gst < 0.0:
            gst += 24.0
        time._st = gst
    # return sidereal hours in the half-open range [0, 24).
    return time._st

def _inverse_terra(ovec, st):
    # Convert from AU to kilometers
    x = ovec[0] * KM_PER_AU
    y = ovec[1] * KM_PER_AU
    z = ovec[2] * KM_PER_AU
    p = math.hypot(x, y)
    if p < 1.0e-6:
        # Special case: within 1 millimeter of a pole!
        # Use arbitrary longitude, and latitude determined by polarity of z.
        lon_deg = 0.0
        if z > 0.0:
            lat_deg = +90.0
        else:
            lat_deg = -90.0
        # Elevation is calculated directly from z
        height_km = abs(z) - _EARTH_POLAR_RADIUS_KM
    else:
        stlocl = math.atan2(y, x)
        # Calculate exact longitude.
        lon_deg = math.degrees(stlocl) - (15.0 * st)
        # Normalize longitude to the range (-180, +180].
        while lon_deg <= -180.0:
            lon_deg += 360.0
        while lon_deg > +180.0:
            lon_deg -= 360.0
        # Numerically solve for exact latitude, using Newton's Method.
        # Start with initial latitude estimate, based on a spherical Earth.
        lat = math.atan2(z, p)
        while True:
            # Calculate the error function W(lat).
            # We try to find the root of W, meaning where the error is 0.
            cos = math.cos(lat)
            sin = math.sin(lat)
            factor = (_EARTH_FLATTENING_SQUARED - 1)*_EARTH_EQUATORIAL_RADIUS_KM
            cos2 = cos*cos
            sin2 = sin*sin
            radicand = cos2 + _EARTH_FLATTENING_SQUARED*sin2
            denom = math.sqrt(radicand)
            W = (factor*sin*cos)/denom - z*cos + p*sin
            if abs(W) < 1.0e-12:
                # The error is now negligible
                break
            # Error is still too large. Find the next estimate.
            # Calculate D = the derivative of W with respect to lat.
            D = factor*((cos2 - sin2)/denom - sin2*cos2*(_EARTH_FLATTENING_SQUARED - 1)/(factor*radicand)) + z*sin + p*cos
            lat -= W/D
        # We now have a solution for the latitude in radians.
        lat_deg = math.degrees(lat)
        # Solve for exact height in meters.
        # There are two formulas I can use. Use whichever has the less risky denominator.
        adjust = _EARTH_EQUATORIAL_RADIUS_KM / denom
        if abs(sin) > abs(cos):
            height_km = z/sin - _EARTH_FLATTENING_SQUARED*adjust
        else:
            height_km = p/cos - adjust
    return Observer(lat_deg, lon_deg, 1000*height_km)

def _terra_posvel(observer, st):
    phi = math.radians(observer.latitude)
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    c = 1.0 / math.hypot(cosphi, sinphi*_EARTH_FLATTENING)
    s = _EARTH_FLATTENING_SQUARED * c
    ht_km = observer.height / 1000.0
    ach = _EARTH_EQUATORIAL_RADIUS_KM*c + ht_km
    ash = _EARTH_EQUATORIAL_RADIUS_KM*s + ht_km
    stlocl = math.radians(15.0*st + observer.longitude)
    sinst = math.sin(stlocl)
    cosst = math.cos(stlocl)
    return [
        ach * cosphi * cosst / KM_PER_AU,
        ach * cosphi * sinst / KM_PER_AU,
        ash * sinphi / KM_PER_AU,
        -_ANGVEL * ach * cosphi * sinst * 86400 / KM_PER_AU,
        +_ANGVEL * ach * cosphi * cosst * 86400 / KM_PER_AU,
        0.0
    ]

def _terra(observer, st):
    return _terra_posvel(observer, st)[0:3]

def _geo_pos(time, observer):
    gast = SiderealTime(time)
    pos1 = _terra(observer, gast)
    pos2 = _nutation(pos1, time, _PrecessDir.Into2000)
    outpos = _precession(pos2, time, _PrecessDir.Into2000)
    return outpos

def _spin(angle, pos1):
    angr = math.radians(angle)
    cosang = math.cos(angr)
    sinang = math.sin(angr)
    return [
        +cosang*pos1[0] + sinang*pos1[1],
        -sinang*pos1[0] + cosang*pos1[1],
        pos1[2]
    ]

#----------------------------------------------------------------------------
# BEGIN CalcMoon

def _Array1(xmin, xmax):
    return dict((key, 0j) for key in range(xmin, 1+xmax))

def _Array2(xmin, xmax, ymin, ymax):
    return dict((key, _Array1(ymin, ymax)) for key in range(xmin, 1+xmax))

class _moonpos:
    def __init__(self, lon, lat, dist):
        self.geo_eclip_lon = lon
        self.geo_eclip_lat = lat
        self.distance_au = dist

def _CalcMoon(time):
    global _CalcMoonCount
    _CalcMoonCount += 1

    T = time.tt / 36525
    ex = _Array2(-6, 6, 1, 4)

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

        ex[0][I] = complex(1, 0)
        ex[1][I] = complex(FAC * math.cos(ARG), FAC * math.sin(ARG))

        J = 2
        while J <= MAX:
            ex[J][I] = ex[J-1][I] * ex[1][I]
            J += 1

        J = 1
        while J <= MAX:
            ex[-J][I] = ex[J][I].conjugate()
            J += 1

        I += 1



    # AddSol(13.902000, 14.060000, -0.001000, 0.260700, 0.000000, 0.000000, 0.000000, 4.000000)
    z = ex[4][4]
    DLAM  += 13.902 * z.imag
    DS    += 14.06 * z.imag
    GAM1C += -0.001 * z.real
    SINPI += 0.2607 * z.real

    # AddSol(0.403000, -4.010000, 0.394000, 0.002300, 0.000000, 0.000000, 0.000000, 3.000000)
    z = ex[3][4]
    DLAM  += 0.403 * z.imag
    DS    += -4.01 * z.imag
    GAM1C += 0.394 * z.real
    SINPI += 0.0023 * z.real

    # AddSol(2369.912000, 2373.360000, 0.601000, 28.233300, 0.000000, 0.000000, 0.000000, 2.000000)
    z = ex[2][4]
    DLAM  += 2369.912 * z.imag
    DS    += 2373.36 * z.imag
    GAM1C += 0.601 * z.real
    SINPI += 28.2333 * z.real

    # AddSol(-125.154000, -112.790000, -0.725000, -0.978100, 0.000000, 0.000000, 0.000000, 1.000000)
    z = ex[1][4]
    DLAM  += -125.154 * z.imag
    DS    += -112.79 * z.imag
    GAM1C += -0.725 * z.real
    SINPI += -0.9781 * z.real

    # AddSol(1.979000, 6.980000, -0.445000, 0.043300, 1.000000, 0.000000, 0.000000, 4.000000)
    z = ex[1][1] * ex[4][4]
    DLAM  += 1.979 * z.imag
    DS    += 6.98 * z.imag
    GAM1C += -0.445 * z.real
    SINPI += 0.0433 * z.real

    # AddSol(191.953000, 192.720000, 0.029000, 3.086100, 1.000000, 0.000000, 0.000000, 2.000000)
    z = ex[1][1] * ex[2][4]
    DLAM  += 191.953 * z.imag
    DS    += 192.72 * z.imag
    GAM1C += 0.029 * z.real
    SINPI += 3.0861 * z.real

    # AddSol(-8.466000, -13.510000, 0.455000, -0.109300, 1.000000, 0.000000, 0.000000, 1.000000)
    z = ex[1][1] * ex[1][4]
    DLAM  += -8.466 * z.imag
    DS    += -13.51 * z.imag
    GAM1C += 0.455 * z.real
    SINPI += -0.1093 * z.real

    # AddSol(22639.500000, 22609.070000, 0.079000, 186.539800, 1.000000, 0.000000, 0.000000, 0.000000)
    z = ex[1][1]
    DLAM  += 22639.500 * z.imag
    DS    += 22609.07 * z.imag
    GAM1C += 0.079 * z.real
    SINPI += 186.5398 * z.real

    # AddSol(18.609000, 3.590000, -0.094000, 0.011800, 1.000000, 0.000000, 0.000000, -1.000000)
    z = ex[1][1] * ex[-1][4]
    DLAM  += 18.609 * z.imag
    DS    += 3.59 * z.imag
    GAM1C += -0.094 * z.real
    SINPI += 0.0118 * z.real

    # AddSol(-4586.465000, -4578.130000, -0.077000, 34.311700, 1.000000, 0.000000, 0.000000, -2.000000)
    z = ex[1][1] * ex[-2][4]
    DLAM  += -4586.465 * z.imag
    DS    += -4578.13 * z.imag
    GAM1C += -0.077 * z.real
    SINPI += 34.3117 * z.real

    # AddSol(3.215000, 5.440000, 0.192000, -0.038600, 1.000000, 0.000000, 0.000000, -3.000000)
    z = ex[1][1] * ex[-3][4]
    DLAM  += 3.215 * z.imag
    DS    += 5.44 * z.imag
    GAM1C += 0.192 * z.real
    SINPI += -0.0386 * z.real

    # AddSol(-38.428000, -38.640000, 0.001000, 0.600800, 1.000000, 0.000000, 0.000000, -4.000000)
    z = ex[1][1] * ex[-4][4]
    DLAM  += -38.428 * z.imag
    DS    += -38.64 * z.imag
    GAM1C += 0.001 * z.real
    SINPI += 0.6008 * z.real

    # AddSol(-0.393000, -1.430000, -0.092000, 0.008600, 1.000000, 0.000000, 0.000000, -6.000000)
    z = ex[1][1] * ex[-6][4]
    DLAM  += -0.393 * z.imag
    DS    += -1.43 * z.imag
    GAM1C += -0.092 * z.real
    SINPI += 0.0086 * z.real

    # AddSol(-0.289000, -1.590000, 0.123000, -0.005300, 0.000000, 1.000000, 0.000000, 4.000000)
    z = ex[1][2] * ex[4][4]
    DLAM  += -0.289 * z.imag
    DS    += -1.59 * z.imag
    GAM1C += 0.123 * z.real
    SINPI += -0.0053 * z.real

    # AddSol(-24.420000, -25.100000, 0.040000, -0.300000, 0.000000, 1.000000, 0.000000, 2.000000)
    z = ex[1][2] * ex[2][4]
    DLAM  += -24.420 * z.imag
    DS    += -25.10 * z.imag
    GAM1C += 0.040 * z.real
    SINPI += -0.3000 * z.real

    # AddSol(18.023000, 17.930000, 0.007000, 0.149400, 0.000000, 1.000000, 0.000000, 1.000000)
    z = ex[1][2] * ex[1][4]
    DLAM  += 18.023 * z.imag
    DS    += 17.93 * z.imag
    GAM1C += 0.007 * z.real
    SINPI += 0.1494 * z.real

    # AddSol(-668.146000, -126.980000, -1.302000, -0.399700, 0.000000, 1.000000, 0.000000, 0.000000)
    z = ex[1][2]
    DLAM  += -668.146 * z.imag
    DS    += -126.98 * z.imag
    GAM1C += -1.302 * z.real
    SINPI += -0.3997 * z.real

    # AddSol(0.560000, 0.320000, -0.001000, -0.003700, 0.000000, 1.000000, 0.000000, -1.000000)
    z = ex[1][2] * ex[-1][4]
    DLAM  += 0.560 * z.imag
    DS    += 0.32 * z.imag
    GAM1C += -0.001 * z.real
    SINPI += -0.0037 * z.real

    # AddSol(-165.145000, -165.060000, 0.054000, 1.917800, 0.000000, 1.000000, 0.000000, -2.000000)
    z = ex[1][2] * ex[-2][4]
    DLAM  += -165.145 * z.imag
    DS    += -165.06 * z.imag
    GAM1C += 0.054 * z.real
    SINPI += 1.9178 * z.real

    # AddSol(-1.877000, -6.460000, -0.416000, 0.033900, 0.000000, 1.000000, 0.000000, -4.000000)
    z = ex[1][2] * ex[-4][4]
    DLAM  += -1.877 * z.imag
    DS    += -6.46 * z.imag
    GAM1C += -0.416 * z.real
    SINPI += 0.0339 * z.real

    # AddSol(0.213000, 1.020000, -0.074000, 0.005400, 2.000000, 0.000000, 0.000000, 4.000000)
    z = ex[2][1] * ex[4][4]
    DLAM  += 0.213 * z.imag
    DS    += 1.02 * z.imag
    GAM1C += -0.074 * z.real
    SINPI += 0.0054 * z.real

    # AddSol(14.387000, 14.780000, -0.017000, 0.283300, 2.000000, 0.000000, 0.000000, 2.000000)
    z = ex[2][1] * ex[2][4]
    DLAM  += 14.387 * z.imag
    DS    += 14.78 * z.imag
    GAM1C += -0.017 * z.real
    SINPI += 0.2833 * z.real

    # AddSol(-0.586000, -1.200000, 0.054000, -0.010000, 2.000000, 0.000000, 0.000000, 1.000000)
    z = ex[2][1] * ex[1][4]
    DLAM  += -0.586 * z.imag
    DS    += -1.20 * z.imag
    GAM1C += 0.054 * z.real
    SINPI += -0.0100 * z.real

    # AddSol(769.016000, 767.960000, 0.107000, 10.165700, 2.000000, 0.000000, 0.000000, 0.000000)
    z = ex[2][1]
    DLAM  += 769.016 * z.imag
    DS    += 767.96 * z.imag
    GAM1C += 0.107 * z.real
    SINPI += 10.1657 * z.real

    # AddSol(1.750000, 2.010000, -0.018000, 0.015500, 2.000000, 0.000000, 0.000000, -1.000000)
    z = ex[2][1] * ex[-1][4]
    DLAM  += 1.750 * z.imag
    DS    += 2.01 * z.imag
    GAM1C += -0.018 * z.real
    SINPI += 0.0155 * z.real

    # AddSol(-211.656000, -152.530000, 5.679000, -0.303900, 2.000000, 0.000000, 0.000000, -2.000000)
    z = ex[2][1] * ex[-2][4]
    DLAM  += -211.656 * z.imag
    DS    += -152.53 * z.imag
    GAM1C += 5.679 * z.real
    SINPI += -0.3039 * z.real

    # AddSol(1.225000, 0.910000, -0.030000, -0.008800, 2.000000, 0.000000, 0.000000, -3.000000)
    z = ex[2][1] * ex[-3][4]
    DLAM  += 1.225 * z.imag
    DS    += 0.91 * z.imag
    GAM1C += -0.030 * z.real
    SINPI += -0.0088 * z.real

    # AddSol(-30.773000, -34.070000, -0.308000, 0.372200, 2.000000, 0.000000, 0.000000, -4.000000)
    z = ex[2][1] * ex[-4][4]
    DLAM  += -30.773 * z.imag
    DS    += -34.07 * z.imag
    GAM1C += -0.308 * z.real
    SINPI += 0.3722 * z.real

    # AddSol(-0.570000, -1.400000, -0.074000, 0.010900, 2.000000, 0.000000, 0.000000, -6.000000)
    z = ex[2][1] * ex[-6][4]
    DLAM  += -0.570 * z.imag
    DS    += -1.40 * z.imag
    GAM1C += -0.074 * z.real
    SINPI += 0.0109 * z.real

    # AddSol(-2.921000, -11.750000, 0.787000, -0.048400, 1.000000, 1.000000, 0.000000, 2.000000)
    z = ex[1][1] * ex[1][2] * ex[2][4]
    DLAM  += -2.921 * z.imag
    DS    += -11.75 * z.imag
    GAM1C += 0.787 * z.real
    SINPI += -0.0484 * z.real

    # AddSol(1.267000, 1.520000, -0.022000, 0.016400, 1.000000, 1.000000, 0.000000, 1.000000)
    z = ex[1][1] * ex[1][2] * ex[1][4]
    DLAM  += 1.267 * z.imag
    DS    += 1.52 * z.imag
    GAM1C += -0.022 * z.real
    SINPI += 0.0164 * z.real

    # AddSol(-109.673000, -115.180000, 0.461000, -0.949000, 1.000000, 1.000000, 0.000000, 0.000000)
    z = ex[1][1] * ex[1][2]
    DLAM  += -109.673 * z.imag
    DS    += -115.18 * z.imag
    GAM1C += 0.461 * z.real
    SINPI += -0.9490 * z.real

    # AddSol(-205.962000, -182.360000, 2.056000, 1.443700, 1.000000, 1.000000, 0.000000, -2.000000)
    z = ex[1][1] * ex[1][2] * ex[-2][4]
    DLAM  += -205.962 * z.imag
    DS    += -182.36 * z.imag
    GAM1C += 2.056 * z.real
    SINPI += 1.4437 * z.real

    # AddSol(0.233000, 0.360000, 0.012000, -0.002500, 1.000000, 1.000000, 0.000000, -3.000000)
    z = ex[1][1] * ex[1][2] * ex[-3][4]
    DLAM  += 0.233 * z.imag
    DS    += 0.36 * z.imag
    GAM1C += 0.012 * z.real
    SINPI += -0.0025 * z.real

    # AddSol(-4.391000, -9.660000, -0.471000, 0.067300, 1.000000, 1.000000, 0.000000, -4.000000)
    z = ex[1][1] * ex[1][2] * ex[-4][4]
    DLAM  += -4.391 * z.imag
    DS    += -9.66 * z.imag
    GAM1C += -0.471 * z.real
    SINPI += 0.0673 * z.real

    # AddSol(0.283000, 1.530000, -0.111000, 0.006000, 1.000000, -1.000000, 0.000000, 4.000000)
    z = ex[1][1] * ex[-1][2] * ex[4][4]
    DLAM  += 0.283 * z.imag
    DS    += 1.53 * z.imag
    GAM1C += -0.111 * z.real
    SINPI += 0.0060 * z.real

    # AddSol(14.577000, 31.700000, -1.540000, 0.230200, 1.000000, -1.000000, 0.000000, 2.000000)
    z = ex[1][1] * ex[-1][2] * ex[2][4]
    DLAM  += 14.577 * z.imag
    DS    += 31.70 * z.imag
    GAM1C += -1.540 * z.real
    SINPI += 0.2302 * z.real

    # AddSol(147.687000, 138.760000, 0.679000, 1.152800, 1.000000, -1.000000, 0.000000, 0.000000)
    z = ex[1][1] * ex[-1][2]
    DLAM  += 147.687 * z.imag
    DS    += 138.76 * z.imag
    GAM1C += 0.679 * z.real
    SINPI += 1.1528 * z.real

    # AddSol(-1.089000, 0.550000, 0.021000, 0.000000, 1.000000, -1.000000, 0.000000, -1.000000)
    z = ex[1][1] * ex[-1][2] * ex[-1][4]
    DLAM  += -1.089 * z.imag
    DS    += 0.55 * z.imag
    GAM1C += 0.021 * z.real

    # AddSol(28.475000, 23.590000, -0.443000, -0.225700, 1.000000, -1.000000, 0.000000, -2.000000)
    z = ex[1][1] * ex[-1][2] * ex[-2][4]
    DLAM  += 28.475 * z.imag
    DS    += 23.59 * z.imag
    GAM1C += -0.443 * z.real
    SINPI += -0.2257 * z.real

    # AddSol(-0.276000, -0.380000, -0.006000, -0.003600, 1.000000, -1.000000, 0.000000, -3.000000)
    z = ex[1][1] * ex[-1][2] * ex[-3][4]
    DLAM  += -0.276 * z.imag
    DS    += -0.38 * z.imag
    GAM1C += -0.006 * z.real
    SINPI += -0.0036 * z.real

    # AddSol(0.636000, 2.270000, 0.146000, -0.010200, 1.000000, -1.000000, 0.000000, -4.000000)
    z = ex[1][1] * ex[-1][2] * ex[-4][4]
    DLAM  += 0.636 * z.imag
    DS    += 2.27 * z.imag
    GAM1C += 0.146 * z.real
    SINPI += -0.0102 * z.real

    # AddSol(-0.189000, -1.680000, 0.131000, -0.002800, 0.000000, 2.000000, 0.000000, 2.000000)
    z = ex[2][2] * ex[2][4]
    DLAM  += -0.189 * z.imag
    DS    += -1.68 * z.imag
    GAM1C += 0.131 * z.real
    SINPI += -0.0028 * z.real

    # AddSol(-7.486000, -0.660000, -0.037000, -0.008600, 0.000000, 2.000000, 0.000000, 0.000000)
    z = ex[2][2]
    DLAM  += -7.486 * z.imag
    DS    += -0.66 * z.imag
    GAM1C += -0.037 * z.real
    SINPI += -0.0086 * z.real

    # AddSol(-8.096000, -16.350000, -0.740000, 0.091800, 0.000000, 2.000000, 0.000000, -2.000000)
    z = ex[2][2] * ex[-2][4]
    DLAM  += -8.096 * z.imag
    DS    += -16.35 * z.imag
    GAM1C += -0.740 * z.real
    SINPI += 0.0918 * z.real

    # AddSol(-5.741000, -0.040000, 0.000000, -0.000900, 0.000000, 0.000000, 2.000000, 2.000000)
    z = ex[2][3] * ex[2][4]
    DLAM  += -5.741 * z.imag
    DS    += -0.04 * z.imag
    SINPI += -0.0009 * z.real

    # AddSol(0.255000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 1.000000)
    z = ex[2][3] * ex[1][4]
    DLAM  += 0.255 * z.imag

    # AddSol(-411.608000, -0.200000, 0.000000, -0.012400, 0.000000, 0.000000, 2.000000, 0.000000)
    z = ex[2][3]
    DLAM  += -411.608 * z.imag
    DS    += -0.20 * z.imag
    SINPI += -0.0124 * z.real

    # AddSol(0.584000, 0.840000, 0.000000, 0.007100, 0.000000, 0.000000, 2.000000, -1.000000)
    z = ex[2][3] * ex[-1][4]
    DLAM  += 0.584 * z.imag
    DS    += 0.84 * z.imag
    SINPI += 0.0071 * z.real

    # AddSol(-55.173000, -52.140000, 0.000000, -0.105200, 0.000000, 0.000000, 2.000000, -2.000000)
    z = ex[2][3] * ex[-2][4]
    DLAM  += -55.173 * z.imag
    DS    += -52.14 * z.imag
    SINPI += -0.1052 * z.real

    # AddSol(0.254000, 0.250000, 0.000000, -0.001700, 0.000000, 0.000000, 2.000000, -3.000000)
    z = ex[2][3] * ex[-3][4]
    DLAM  += 0.254 * z.imag
    DS    += 0.25 * z.imag
    SINPI += -0.0017 * z.real

    # AddSol(0.025000, -1.670000, 0.000000, 0.003100, 0.000000, 0.000000, 2.000000, -4.000000)
    z = ex[2][3] * ex[-4][4]
    DLAM  += 0.025 * z.imag
    DS    += -1.67 * z.imag
    SINPI += 0.0031 * z.real

    # AddSol(1.060000, 2.960000, -0.166000, 0.024300, 3.000000, 0.000000, 0.000000, 2.000000)
    z = ex[3][1] * ex[2][4]
    DLAM  += 1.060 * z.imag
    DS    += 2.96 * z.imag
    GAM1C += -0.166 * z.real
    SINPI += 0.0243 * z.real

    # AddSol(36.124000, 50.640000, -1.300000, 0.621500, 3.000000, 0.000000, 0.000000, 0.000000)
    z = ex[3][1]
    DLAM  += 36.124 * z.imag
    DS    += 50.64 * z.imag
    GAM1C += -1.300 * z.real
    SINPI += 0.6215 * z.real

    # AddSol(-13.193000, -16.400000, 0.258000, -0.118700, 3.000000, 0.000000, 0.000000, -2.000000)
    z = ex[3][1] * ex[-2][4]
    DLAM  += -13.193 * z.imag
    DS    += -16.40 * z.imag
    GAM1C += 0.258 * z.real
    SINPI += -0.1187 * z.real

    # AddSol(-1.187000, -0.740000, 0.042000, 0.007400, 3.000000, 0.000000, 0.000000, -4.000000)
    z = ex[3][1] * ex[-4][4]
    DLAM  += -1.187 * z.imag
    DS    += -0.74 * z.imag
    GAM1C += 0.042 * z.real
    SINPI += 0.0074 * z.real

    # AddSol(-0.293000, -0.310000, -0.002000, 0.004600, 3.000000, 0.000000, 0.000000, -6.000000)
    z = ex[3][1] * ex[-6][4]
    DLAM  += -0.293 * z.imag
    DS    += -0.31 * z.imag
    GAM1C += -0.002 * z.real
    SINPI += 0.0046 * z.real

    # AddSol(-0.290000, -1.450000, 0.116000, -0.005100, 2.000000, 1.000000, 0.000000, 2.000000)
    z = ex[2][1] * ex[1][2] * ex[2][4]
    DLAM  += -0.290 * z.imag
    DS    += -1.45 * z.imag
    GAM1C += 0.116 * z.real
    SINPI += -0.0051 * z.real

    # AddSol(-7.649000, -10.560000, 0.259000, -0.103800, 2.000000, 1.000000, 0.000000, 0.000000)
    z = ex[2][1] * ex[1][2]
    DLAM  += -7.649 * z.imag
    DS    += -10.56 * z.imag
    GAM1C += 0.259 * z.real
    SINPI += -0.1038 * z.real

    # AddSol(-8.627000, -7.590000, 0.078000, -0.019200, 2.000000, 1.000000, 0.000000, -2.000000)
    z = ex[2][1] * ex[1][2] * ex[-2][4]
    DLAM  += -8.627 * z.imag
    DS    += -7.59 * z.imag
    GAM1C += 0.078 * z.real
    SINPI += -0.0192 * z.real

    # AddSol(-2.740000, -2.540000, 0.022000, 0.032400, 2.000000, 1.000000, 0.000000, -4.000000)
    z = ex[2][1] * ex[1][2] * ex[-4][4]
    DLAM  += -2.740 * z.imag
    DS    += -2.54 * z.imag
    GAM1C += 0.022 * z.real
    SINPI += 0.0324 * z.real

    # AddSol(1.181000, 3.320000, -0.212000, 0.021300, 2.000000, -1.000000, 0.000000, 2.000000)
    z = ex[2][1] * ex[-1][2] * ex[2][4]
    DLAM  += 1.181 * z.imag
    DS    += 3.32 * z.imag
    GAM1C += -0.212 * z.real
    SINPI += 0.0213 * z.real

    # AddSol(9.703000, 11.670000, -0.151000, 0.126800, 2.000000, -1.000000, 0.000000, 0.000000)
    z = ex[2][1] * ex[-1][2]
    DLAM  += 9.703 * z.imag
    DS    += 11.67 * z.imag
    GAM1C += -0.151 * z.real
    SINPI += 0.1268 * z.real

    # AddSol(-0.352000, -0.370000, 0.001000, -0.002800, 2.000000, -1.000000, 0.000000, -1.000000)
    z = ex[2][1] * ex[-1][2] * ex[-1][4]
    DLAM  += -0.352 * z.imag
    DS    += -0.37 * z.imag
    GAM1C += 0.001 * z.real
    SINPI += -0.0028 * z.real

    # AddSol(-2.494000, -1.170000, -0.003000, -0.001700, 2.000000, -1.000000, 0.000000, -2.000000)
    z = ex[2][1] * ex[-1][2] * ex[-2][4]
    DLAM  += -2.494 * z.imag
    DS    += -1.17 * z.imag
    GAM1C += -0.003 * z.real
    SINPI += -0.0017 * z.real

    # AddSol(0.360000, 0.200000, -0.012000, -0.004300, 2.000000, -1.000000, 0.000000, -4.000000)
    z = ex[2][1] * ex[-1][2] * ex[-4][4]
    DLAM  += 0.360 * z.imag
    DS    += 0.20 * z.imag
    GAM1C += -0.012 * z.real
    SINPI += -0.0043 * z.real

    # AddSol(-1.167000, -1.250000, 0.008000, -0.010600, 1.000000, 2.000000, 0.000000, 0.000000)
    z = ex[1][1] * ex[2][2]
    DLAM  += -1.167 * z.imag
    DS    += -1.25 * z.imag
    GAM1C += 0.008 * z.real
    SINPI += -0.0106 * z.real

    # AddSol(-7.412000, -6.120000, 0.117000, 0.048400, 1.000000, 2.000000, 0.000000, -2.000000)
    z = ex[1][1] * ex[2][2] * ex[-2][4]
    DLAM  += -7.412 * z.imag
    DS    += -6.12 * z.imag
    GAM1C += 0.117 * z.real
    SINPI += 0.0484 * z.real

    # AddSol(-0.311000, -0.650000, -0.032000, 0.004400, 1.000000, 2.000000, 0.000000, -4.000000)
    z = ex[1][1] * ex[2][2] * ex[-4][4]
    DLAM  += -0.311 * z.imag
    DS    += -0.65 * z.imag
    GAM1C += -0.032 * z.real
    SINPI += 0.0044 * z.real

    # AddSol(0.757000, 1.820000, -0.105000, 0.011200, 1.000000, -2.000000, 0.000000, 2.000000)
    z = ex[1][1] * ex[-2][2] * ex[2][4]
    DLAM  += 0.757 * z.imag
    DS    += 1.82 * z.imag
    GAM1C += -0.105 * z.real
    SINPI += 0.0112 * z.real

    # AddSol(2.580000, 2.320000, 0.027000, 0.019600, 1.000000, -2.000000, 0.000000, 0.000000)
    z = ex[1][1] * ex[-2][2]
    DLAM  += 2.580 * z.imag
    DS    += 2.32 * z.imag
    GAM1C += 0.027 * z.real
    SINPI += 0.0196 * z.real

    # AddSol(2.533000, 2.400000, -0.014000, -0.021200, 1.000000, -2.000000, 0.000000, -2.000000)
    z = ex[1][1] * ex[-2][2] * ex[-2][4]
    DLAM  += 2.533 * z.imag
    DS    += 2.40 * z.imag
    GAM1C += -0.014 * z.real
    SINPI += -0.0212 * z.real

    # AddSol(-0.344000, -0.570000, -0.025000, 0.003600, 0.000000, 3.000000, 0.000000, -2.000000)
    z = ex[3][2] * ex[-2][4]
    DLAM  += -0.344 * z.imag
    DS    += -0.57 * z.imag
    GAM1C += -0.025 * z.real
    SINPI += 0.0036 * z.real

    # AddSol(-0.992000, -0.020000, 0.000000, 0.000000, 1.000000, 0.000000, 2.000000, 2.000000)
    z = ex[1][1] * ex[2][3] * ex[2][4]
    DLAM  += -0.992 * z.imag
    DS    += -0.02 * z.imag

    # AddSol(-45.099000, -0.020000, 0.000000, -0.001000, 1.000000, 0.000000, 2.000000, 0.000000)
    z = ex[1][1] * ex[2][3]
    DLAM  += -45.099 * z.imag
    DS    += -0.02 * z.imag
    SINPI += -0.0010 * z.real

    # AddSol(-0.179000, -9.520000, 0.000000, -0.083300, 1.000000, 0.000000, 2.000000, -2.000000)
    z = ex[1][1] * ex[2][3] * ex[-2][4]
    DLAM  += -0.179 * z.imag
    DS    += -9.52 * z.imag
    SINPI += -0.0833 * z.real

    # AddSol(-0.301000, -0.330000, 0.000000, 0.001400, 1.000000, 0.000000, 2.000000, -4.000000)
    z = ex[1][1] * ex[2][3] * ex[-4][4]
    DLAM  += -0.301 * z.imag
    DS    += -0.33 * z.imag
    SINPI += 0.0014 * z.real

    # AddSol(-6.382000, -3.370000, 0.000000, -0.048100, 1.000000, 0.000000, -2.000000, 2.000000)
    z = ex[1][1] * ex[-2][3] * ex[2][4]
    DLAM  += -6.382 * z.imag
    DS    += -3.37 * z.imag
    SINPI += -0.0481 * z.real

    # AddSol(39.528000, 85.130000, 0.000000, -0.713600, 1.000000, 0.000000, -2.000000, 0.000000)
    z = ex[1][1] * ex[-2][3]
    DLAM  += 39.528 * z.imag
    DS    += 85.13 * z.imag
    SINPI += -0.7136 * z.real

    # AddSol(9.366000, 0.710000, 0.000000, -0.011200, 1.000000, 0.000000, -2.000000, -2.000000)
    z = ex[1][1] * ex[-2][3] * ex[-2][4]
    DLAM  += 9.366 * z.imag
    DS    += 0.71 * z.imag
    SINPI += -0.0112 * z.real

    # AddSol(0.202000, 0.020000, 0.000000, 0.000000, 1.000000, 0.000000, -2.000000, -4.000000)
    z = ex[1][1] * ex[-2][3] * ex[-4][4]
    DLAM  += 0.202 * z.imag
    DS    += 0.02 * z.imag

    # AddSol(0.415000, 0.100000, 0.000000, 0.001300, 0.000000, 1.000000, 2.000000, 0.000000)
    z = ex[1][2] * ex[2][3]
    DLAM  += 0.415 * z.imag
    DS    += 0.10 * z.imag
    SINPI += 0.0013 * z.real

    # AddSol(-2.152000, -2.260000, 0.000000, -0.006600, 0.000000, 1.000000, 2.000000, -2.000000)
    z = ex[1][2] * ex[2][3] * ex[-2][4]
    DLAM  += -2.152 * z.imag
    DS    += -2.26 * z.imag
    SINPI += -0.0066 * z.real

    # AddSol(-1.440000, -1.300000, 0.000000, 0.001400, 0.000000, 1.000000, -2.000000, 2.000000)
    z = ex[1][2] * ex[-2][3] * ex[2][4]
    DLAM  += -1.440 * z.imag
    DS    += -1.30 * z.imag
    SINPI += 0.0014 * z.real

    # AddSol(0.384000, -0.040000, 0.000000, 0.000000, 0.000000, 1.000000, -2.000000, -2.000000)
    z = ex[1][2] * ex[-2][3] * ex[-2][4]
    DLAM  += 0.384 * z.imag
    DS    += -0.04 * z.imag

    # AddSol(1.938000, 3.600000, -0.145000, 0.040100, 4.000000, 0.000000, 0.000000, 0.000000)
    z = ex[4][1]
    DLAM  += 1.938 * z.imag
    DS    += 3.60 * z.imag
    GAM1C += -0.145 * z.real
    SINPI += 0.0401 * z.real

    # AddSol(-0.952000, -1.580000, 0.052000, -0.013000, 4.000000, 0.000000, 0.000000, -2.000000)
    z = ex[4][1] * ex[-2][4]
    DLAM  += -0.952 * z.imag
    DS    += -1.58 * z.imag
    GAM1C += 0.052 * z.real
    SINPI += -0.0130 * z.real

    # AddSol(-0.551000, -0.940000, 0.032000, -0.009700, 3.000000, 1.000000, 0.000000, 0.000000)
    z = ex[3][1] * ex[1][2]
    DLAM  += -0.551 * z.imag
    DS    += -0.94 * z.imag
    GAM1C += 0.032 * z.real
    SINPI += -0.0097 * z.real

    # AddSol(-0.482000, -0.570000, 0.005000, -0.004500, 3.000000, 1.000000, 0.000000, -2.000000)
    z = ex[3][1] * ex[1][2] * ex[-2][4]
    DLAM  += -0.482 * z.imag
    DS    += -0.57 * z.imag
    GAM1C += 0.005 * z.real
    SINPI += -0.0045 * z.real

    # AddSol(0.681000, 0.960000, -0.026000, 0.011500, 3.000000, -1.000000, 0.000000, 0.000000)
    z = ex[3][1] * ex[-1][2]
    DLAM  += 0.681 * z.imag
    DS    += 0.96 * z.imag
    GAM1C += -0.026 * z.real
    SINPI += 0.0115 * z.real

    # AddSol(-0.297000, -0.270000, 0.002000, -0.000900, 2.000000, 2.000000, 0.000000, -2.000000)
    z = ex[2][1] * ex[2][2] * ex[-2][4]
    DLAM  += -0.297 * z.imag
    DS    += -0.27 * z.imag
    GAM1C += 0.002 * z.real
    SINPI += -0.0009 * z.real

    # AddSol(0.254000, 0.210000, -0.003000, 0.000000, 2.000000, -2.000000, 0.000000, -2.000000)
    z = ex[2][1] * ex[-2][2] * ex[-2][4]
    DLAM  += 0.254 * z.imag
    DS    += 0.21 * z.imag
    GAM1C += -0.003 * z.real

    # AddSol(-0.250000, -0.220000, 0.004000, 0.001400, 1.000000, 3.000000, 0.000000, -2.000000)
    z = ex[1][1] * ex[3][2] * ex[-2][4]
    DLAM  += -0.250 * z.imag
    DS    += -0.22 * z.imag
    GAM1C += 0.004 * z.real
    SINPI += 0.0014 * z.real

    # AddSol(-3.996000, 0.000000, 0.000000, 0.000400, 2.000000, 0.000000, 2.000000, 0.000000)
    z = ex[2][1] * ex[2][3]
    DLAM  += -3.996 * z.imag
    SINPI += 0.0004 * z.real

    # AddSol(0.557000, -0.750000, 0.000000, -0.009000, 2.000000, 0.000000, 2.000000, -2.000000)
    z = ex[2][1] * ex[2][3] * ex[-2][4]
    DLAM  += 0.557 * z.imag
    DS    += -0.75 * z.imag
    SINPI += -0.0090 * z.real

    # AddSol(-0.459000, -0.380000, 0.000000, -0.005300, 2.000000, 0.000000, -2.000000, 2.000000)
    z = ex[2][1] * ex[-2][3] * ex[2][4]
    DLAM  += -0.459 * z.imag
    DS    += -0.38 * z.imag
    SINPI += -0.0053 * z.real

    # AddSol(-1.298000, 0.740000, 0.000000, 0.000400, 2.000000, 0.000000, -2.000000, 0.000000)
    z = ex[2][1] * ex[-2][3]
    DLAM  += -1.298 * z.imag
    DS    += 0.74 * z.imag
    SINPI += 0.0004 * z.real

    # AddSol(0.538000, 1.140000, 0.000000, -0.014100, 2.000000, 0.000000, -2.000000, -2.000000)
    z = ex[2][1] * ex[-2][3] * ex[-2][4]
    DLAM  += 0.538 * z.imag
    DS    += 1.14 * z.imag
    SINPI += -0.0141 * z.real

    # AddSol(0.263000, 0.020000, 0.000000, 0.000000, 1.000000, 1.000000, 2.000000, 0.000000)
    z = ex[1][1] * ex[1][2] * ex[2][3]
    DLAM  += 0.263 * z.imag
    DS    += 0.02 * z.imag

    # AddSol(0.426000, 0.070000, 0.000000, -0.000600, 1.000000, 1.000000, -2.000000, -2.000000)
    z = ex[1][1] * ex[1][2] * ex[-2][3] * ex[-2][4]
    DLAM  += 0.426 * z.imag
    DS    += 0.07 * z.imag
    SINPI += -0.0006 * z.real

    # AddSol(-0.304000, 0.030000, 0.000000, 0.000300, 1.000000, -1.000000, 2.000000, 0.000000)
    z = ex[1][1] * ex[-1][2] * ex[2][3]
    DLAM  += -0.304 * z.imag
    DS    += 0.03 * z.imag
    SINPI += 0.0003 * z.real

    # AddSol(-0.372000, -0.190000, 0.000000, -0.002700, 1.000000, -1.000000, -2.000000, 2.000000)
    z = ex[1][1] * ex[-1][2] * ex[-2][3] * ex[2][4]
    DLAM  += -0.372 * z.imag
    DS    += -0.19 * z.imag
    SINPI += -0.0027 * z.real

    # AddSol(0.418000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 4.000000, 0.000000)
    z = ex[4][3]
    DLAM  += 0.418 * z.imag

    # AddSol(-0.330000, -0.040000, 0.000000, 0.000000, 3.000000, 0.000000, 2.000000, 0.000000)
    z = ex[3][1] * ex[2][3]
    DLAM  += -0.330 * z.imag
    DS    += -0.04 * z.imag


    def ADDN(coeffn, p, q, r, s):
        return coeffn * (ex[p][1] * ex[q][2] * ex[r][3] * ex[s][4]).imag

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
        (_ARC * _EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI)
    )

def GeoMoon(time):
    """Calculates equatorial geocentric position of the Moon at a given time.

    Given a time of observation, calculates the Moon's position as a vector.
    The vector gives the location of the Moon's center relative to the Earth's center
    with x-, y-, and z-components measured in astronomical units.
    The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
    In Astronomy Engine, this orientation is called EQJ.

    This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
    which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
    It is adapted from Turbo Pascal code from the book
    [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
    by Montenbruck and Pfleger.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Moon's position.

    Returns
    -------
    Vector
        The Moon's position as a vector in J2000 Cartesian equatorial coordinates (EQJ).
    """
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
    mpos2 = _precession(mpos1, time, _PrecessDir.Into2000)
    return Vector(mpos2[0], mpos2[1], mpos2[2], time)


def EclipticGeoMoon(time):
    """Calculates spherical ecliptic geocentric position of the Moon.

    Given a time of observation, calculates the Moon's geocentric position
    in ecliptic spherical coordinates. Provides the ecliptic latitude and
    longitude in degrees, and the geocentric distance in astronomical units (AU).
    The ecliptic longitude is measured relative to the equinox of date.

    This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
    which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
    It is adapted from Turbo Pascal code from the book
    [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
    by Montenbruck and Pfleger.

    To calculate an equatorial J2000 vector instead, use #GeoMoon.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Moon's position.

    Returns
    -------
    Spherical
        The Moon's position as a distance, ecliptic latitude, and ecliptic longitude.
    """
    moon = _CalcMoon(time)
    return Spherical(
        math.degrees(moon.geo_eclip_lat),
        math.degrees(moon.geo_eclip_lon),
        moon.distance_au
    )


def GeoMoonState(time):
    """Calculates equatorial geocentric position and velocity of the Moon at a given time.

    Given a time of observation, calculates the Moon's position and velocity vectors.
    The position and velocity are of the Moon's center relative to the Earth's center.
    The position (x, y, z) components are expressed in AU (astronomical units).
    The velocity (vx, vy, vz) components are expressed in AU/day.
    The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
    In Astronomy Engine, this orientation is called EQJ.
    If you need the Moon's position only, and not its velocity,
    it is much more efficient to use #GeoMoon instead.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Moon's position and velocity.

    Returns
    -------
    StateVector
        The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).
    """
    # This is a hack, because trying to figure out how to derive a time
    # derivative for CalcMoon() would be extremely painful!
    # Calculate just before and just after the given time.
    # Average to find position, subtract to find velocity.
    dt = 1.0e-5   # 0.864 seconds
    t1 = time.AddDays(-dt)
    t2 = time.AddDays(+dt)
    r1 = GeoMoon(t1)
    r2 = GeoMoon(t2)
    return StateVector(
        (r1.x + r2.x) / 2,
        (r1.y + r2.y) / 2,
        (r1.z + r2.z) / 2,
        (r2.x - r1.x) / (2 * dt),
        (r2.y - r1.y) / (2 * dt),
        (r2.z - r1.z) / (2 * dt),
        time
    )


def GeoEmbState(time):
    """Calculates the geocentric position and velocity of the Earth/Moon barycenter.

    Given a time of observation, calculates the geocentric position and velocity vectors
    of the Earth/Moon barycenter (EMB).
    The position (x, y, z) components are expressed in AU (astronomical units).
    The velocity (vx, vy, vz) components are expressed in AU/day.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the EMB's geocentric state.

    Returns
    -------
    StateVector
        The EMB's position and velocity vectors in J2000 equatorial coordinates.
    """
    s = GeoMoonState(time)
    d = 1.0 + _EARTH_MOON_MASS_RATIO
    return StateVector(s.x/d, s.y/d, s.z/d, s.vx/d, s.vy/d, s.vz/d, time)

# END CalcMoon
#----------------------------------------------------------------------------
# BEGIN VSOP

_vsop = [
    # Mercury
    [
  [
    [
      [4.40250710144, 0.00000000000, 0.00000000000],
      [0.40989414977, 1.48302034195, 26087.90314157420],
      [0.05046294200, 4.47785489551, 52175.80628314840],
      [0.00855346844, 1.16520322459, 78263.70942472259],
      [0.00165590362, 4.11969163423, 104351.61256629678],
      [0.00034561897, 0.77930768443, 130439.51570787099],
      [0.00007583476, 3.71348404924, 156527.41884944518]
    ],
    [
      [26087.90313685529, 0.00000000000, 0.00000000000],
      [0.01131199811, 6.21874197797, 26087.90314157420],
      [0.00292242298, 3.04449355541, 52175.80628314840],
      [0.00075775081, 6.08568821653, 78263.70942472259],
      [0.00019676525, 2.80965111777, 104351.61256629678]
    ]
  ],
  [
    [
      [0.11737528961, 1.98357498767, 26087.90314157420],
      [0.02388076996, 5.03738959686, 52175.80628314840],
      [0.01222839532, 3.14159265359, 0.00000000000],
      [0.00543251810, 1.79644363964, 78263.70942472259],
      [0.00129778770, 4.83232503958, 104351.61256629678],
      [0.00031866927, 1.58088495658, 130439.51570787099],
      [0.00007963301, 4.60972126127, 156527.41884944518]
    ],
    [
      [0.00274646065, 3.95008450011, 26087.90314157420],
      [0.00099737713, 3.14159265359, 0.00000000000]
    ]
  ],
  [
    [
      [0.39528271651, 0.00000000000, 0.00000000000],
      [0.07834131818, 6.19233722598, 26087.90314157420],
      [0.00795525558, 2.95989690104, 52175.80628314840],
      [0.00121281764, 6.01064153797, 78263.70942472259],
      [0.00021921969, 2.77820093972, 104351.61256629678],
      [0.00004354065, 5.82894543774, 130439.51570787099]
    ],
    [
      [0.00217347740, 4.65617158665, 26087.90314157420],
      [0.00044141826, 1.42385544001, 52175.80628314840]
    ]
  ]
],

    # Venus
    [
  [
    [
      [3.17614666774, 0.00000000000, 0.00000000000],
      [0.01353968419, 5.59313319619, 10213.28554621100],
      [0.00089891645, 5.30650047764, 20426.57109242200],
      [0.00005477194, 4.41630661466, 7860.41939243920],
      [0.00003455741, 2.69964447820, 11790.62908865880],
      [0.00002372061, 2.99377542079, 3930.20969621960],
      [0.00001317168, 5.18668228402, 26.29831979980],
      [0.00001664146, 4.25018630147, 1577.34354244780],
      [0.00001438387, 4.15745084182, 9683.59458111640],
      [0.00001200521, 6.15357116043, 30639.85663863300]
    ],
    [
      [10213.28554621638, 0.00000000000, 0.00000000000],
      [0.00095617813, 2.46406511110, 10213.28554621100],
      [0.00007787201, 0.62478482220, 20426.57109242200]
    ]
  ],
  [
    [
      [0.05923638472, 0.26702775812, 10213.28554621100],
      [0.00040107978, 1.14737178112, 20426.57109242200],
      [0.00032814918, 3.14159265359, 0.00000000000]
    ],
    [
      [0.00287821243, 1.88964962838, 10213.28554621100]
    ]
  ],
  [
    [
      [0.72334820891, 0.00000000000, 0.00000000000],
      [0.00489824182, 4.02151831717, 10213.28554621100],
      [0.00001658058, 4.90206728031, 20426.57109242200],
      [0.00001378043, 1.12846591367, 11790.62908865880],
      [0.00001632096, 2.84548795207, 7860.41939243920],
      [0.00000498395, 2.58682193892, 9683.59458111640],
      [0.00000221985, 2.01346696541, 19367.18916223280],
      [0.00000237454, 2.55136053886, 15720.83878487840]
    ],
    [
      [0.00034551041, 0.89198706276, 10213.28554621100]
    ]
  ]
],

    # Earth
    [
  [
    [
      [1.75347045673, 0.00000000000, 0.00000000000],
      [0.03341656453, 4.66925680415, 6283.07584999140],
      [0.00034894275, 4.62610242189, 12566.15169998280],
      [0.00003417572, 2.82886579754, 3.52311834900],
      [0.00003497056, 2.74411783405, 5753.38488489680],
      [0.00003135899, 3.62767041756, 77713.77146812050],
      [0.00002676218, 4.41808345438, 7860.41939243920],
      [0.00002342691, 6.13516214446, 3930.20969621960],
      [0.00001273165, 2.03709657878, 529.69096509460],
      [0.00001324294, 0.74246341673, 11506.76976979360],
      [0.00000901854, 2.04505446477, 26.29831979980],
      [0.00001199167, 1.10962946234, 1577.34354244780],
      [0.00000857223, 3.50849152283, 398.14900340820],
      [0.00000779786, 1.17882681962, 5223.69391980220],
      [0.00000990250, 5.23268072088, 5884.92684658320],
      [0.00000753141, 2.53339052847, 5507.55323866740],
      [0.00000505267, 4.58292599973, 18849.22754997420],
      [0.00000492392, 4.20505711826, 775.52261132400],
      [0.00000356672, 2.91954114478, 0.06731030280],
      [0.00000284125, 1.89869240932, 796.29800681640],
      [0.00000242879, 0.34481445893, 5486.77784317500],
      [0.00000317087, 5.84901948512, 11790.62908865880],
      [0.00000271112, 0.31486255375, 10977.07880469900],
      [0.00000206217, 4.80646631478, 2544.31441988340],
      [0.00000205478, 1.86953770281, 5573.14280143310],
      [0.00000202318, 2.45767790232, 6069.77675455340],
      [0.00000126225, 1.08295459501, 20.77539549240],
      [0.00000155516, 0.83306084617, 213.29909543800]
    ],
    [
      [6283.07584999140, 0.00000000000, 0.00000000000],
      [0.00206058863, 2.67823455808, 6283.07584999140],
      [0.00004303419, 2.63512233481, 12566.15169998280]
    ],
    [
      [0.00008721859, 1.07253635559, 6283.07584999140]
    ]
  ],
  [
    [
    ],
    [
      [0.00227777722, 3.41376620530, 6283.07584999140],
      [0.00003805678, 3.37063423795, 12566.15169998280]
    ]
  ],
  [
    [
      [1.00013988784, 0.00000000000, 0.00000000000],
      [0.01670699632, 3.09846350258, 6283.07584999140],
      [0.00013956024, 3.05524609456, 12566.15169998280],
      [0.00003083720, 5.19846674381, 77713.77146812050],
      [0.00001628463, 1.17387558054, 5753.38488489680],
      [0.00001575572, 2.84685214877, 7860.41939243920],
      [0.00000924799, 5.45292236722, 11506.76976979360],
      [0.00000542439, 4.56409151453, 3930.20969621960],
      [0.00000472110, 3.66100022149, 5884.92684658320],
      [0.00000085831, 1.27079125277, 161000.68573767410],
      [0.00000057056, 2.01374292245, 83996.84731811189],
      [0.00000055736, 5.24159799170, 71430.69561812909],
      [0.00000174844, 3.01193636733, 18849.22754997420],
      [0.00000243181, 4.27349530790, 11790.62908865880]
    ],
    [
      [0.00103018607, 1.10748968172, 6283.07584999140],
      [0.00001721238, 1.06442300386, 12566.15169998280]
    ],
    [
      [0.00004359385, 5.78455133808, 6283.07584999140]
    ]
  ]
],

    # Mars
    [
  [
    [
      [6.20347711581, 0.00000000000, 0.00000000000],
      [0.18656368093, 5.05037100270, 3340.61242669980],
      [0.01108216816, 5.40099836344, 6681.22485339960],
      [0.00091798406, 5.75478744667, 10021.83728009940],
      [0.00027744987, 5.97049513147, 3.52311834900],
      [0.00010610235, 2.93958560338, 2281.23049651060],
      [0.00012315897, 0.84956094002, 2810.92146160520],
      [0.00008926784, 4.15697846427, 0.01725365220],
      [0.00008715691, 6.11005153139, 13362.44970679920],
      [0.00006797556, 0.36462229657, 398.14900340820],
      [0.00007774872, 3.33968761376, 5621.84292321040],
      [0.00003575078, 1.66186505710, 2544.31441988340],
      [0.00004161108, 0.22814971327, 2942.46342329160],
      [0.00003075252, 0.85696614132, 191.44826611160],
      [0.00002628117, 0.64806124465, 3337.08930835080],
      [0.00002937546, 6.07893711402, 0.06731030280],
      [0.00002389414, 5.03896442664, 796.29800681640],
      [0.00002579844, 0.02996736156, 3344.13554504880],
      [0.00001528141, 1.14979301996, 6151.53388830500],
      [0.00001798806, 0.65634057445, 529.69096509460],
      [0.00001264357, 3.62275122593, 5092.15195811580],
      [0.00001286228, 3.06796065034, 2146.16541647520],
      [0.00001546404, 2.91579701718, 1751.53953141600],
      [0.00001024902, 3.69334099279, 8962.45534991020],
      [0.00000891566, 0.18293837498, 16703.06213349900],
      [0.00000858759, 2.40093811940, 2914.01423582380],
      [0.00000832715, 2.46418619474, 3340.59517304760],
      [0.00000832720, 4.49495782139, 3340.62968035200],
      [0.00000712902, 3.66335473479, 1059.38193018920],
      [0.00000748723, 3.82248614017, 155.42039943420],
      [0.00000723861, 0.67497311481, 3738.76143010800],
      [0.00000635548, 2.92182225127, 8432.76438481560],
      [0.00000655162, 0.48864064125, 3127.31333126180],
      [0.00000550474, 3.81001042328, 0.98032106820],
      [0.00000552750, 4.47479317037, 1748.01641306700],
      [0.00000425966, 0.55364317304, 6283.07584999140],
      [0.00000415131, 0.49662285038, 213.29909543800],
      [0.00000472167, 3.62547124025, 1194.44701022460],
      [0.00000306551, 0.38052848348, 6684.74797174860],
      [0.00000312141, 0.99853944405, 6677.70173505060],
      [0.00000293198, 4.22131299634, 20.77539549240],
      [0.00000302375, 4.48618007156, 3532.06069281140],
      [0.00000274027, 0.54222167059, 3340.54511639700],
      [0.00000281079, 5.88163521788, 1349.86740965880],
      [0.00000231183, 1.28242156993, 3870.30339179440],
      [0.00000283602, 5.76885434940, 3149.16416058820],
      [0.00000236117, 5.75503217933, 3333.49887969900],
      [0.00000274033, 0.13372524985, 3340.67973700260],
      [0.00000299395, 2.78323740866, 6254.62666252360]
    ],
    [
      [3340.61242700512, 0.00000000000, 0.00000000000],
      [0.01457554523, 3.60433733236, 3340.61242669980],
      [0.00168414711, 3.92318567804, 6681.22485339960],
      [0.00020622975, 4.26108844583, 10021.83728009940],
      [0.00003452392, 4.73210393190, 3.52311834900],
      [0.00002586332, 4.60670058555, 13362.44970679920],
      [0.00000841535, 4.45864030426, 2281.23049651060]
    ],
    [
      [0.00058152577, 2.04961712429, 3340.61242669980],
      [0.00013459579, 2.45738706163, 6681.22485339960]
    ]
  ],
  [
    [
      [0.03197134986, 3.76832042431, 3340.61242669980],
      [0.00298033234, 4.10616996305, 6681.22485339960],
      [0.00289104742, 0.00000000000, 0.00000000000],
      [0.00031365539, 4.44651053090, 10021.83728009940],
      [0.00003484100, 4.78812549260, 13362.44970679920]
    ],
    [
      [0.00217310991, 6.04472194776, 3340.61242669980],
      [0.00020976948, 3.14159265359, 0.00000000000],
      [0.00012834709, 1.60810667915, 6681.22485339960]
    ]
  ],
  [
    [
      [1.53033488271, 0.00000000000, 0.00000000000],
      [0.14184953160, 3.47971283528, 3340.61242669980],
      [0.00660776362, 3.81783443019, 6681.22485339960],
      [0.00046179117, 4.15595316782, 10021.83728009940],
      [0.00008109733, 5.55958416318, 2810.92146160520],
      [0.00007485318, 1.77239078402, 5621.84292321040],
      [0.00005523191, 1.36436303770, 2281.23049651060],
      [0.00003825160, 4.49407183687, 13362.44970679920],
      [0.00002306537, 0.09081579001, 2544.31441988340],
      [0.00001999396, 5.36059617709, 3337.08930835080],
      [0.00002484394, 4.92545639920, 2942.46342329160],
      [0.00001960195, 4.74249437639, 3344.13554504880],
      [0.00001167119, 2.11260868341, 5092.15195811580],
      [0.00001102816, 5.00908403998, 398.14900340820],
      [0.00000899066, 4.40791133207, 529.69096509460],
      [0.00000992252, 5.83861961952, 6151.53388830500],
      [0.00000807354, 2.10217065501, 1059.38193018920],
      [0.00000797915, 3.44839203899, 796.29800681640],
      [0.00000740975, 1.49906336885, 2146.16541647520]
    ],
    [
      [0.01107433345, 2.03250524857, 3340.61242669980],
      [0.00103175887, 2.37071847807, 6681.22485339960],
      [0.00012877200, 0.00000000000, 0.00000000000],
      [0.00010815880, 2.70888095665, 10021.83728009940]
    ],
    [
      [0.00044242249, 0.47930604954, 3340.61242669980],
      [0.00008138042, 0.86998389204, 6681.22485339960]
    ]
  ]
],

    # Jupiter
    [
  [
    [
      [0.59954691494, 0.00000000000, 0.00000000000],
      [0.09695898719, 5.06191793158, 529.69096509460],
      [0.00573610142, 1.44406205629, 7.11354700080],
      [0.00306389205, 5.41734730184, 1059.38193018920],
      [0.00097178296, 4.14264726552, 632.78373931320],
      [0.00072903078, 3.64042916389, 522.57741809380],
      [0.00064263975, 3.41145165351, 103.09277421860],
      [0.00039806064, 2.29376740788, 419.48464387520],
      [0.00038857767, 1.27231755835, 316.39186965660],
      [0.00027964629, 1.78454591820, 536.80451209540],
      [0.00013589730, 5.77481040790, 1589.07289528380],
      [0.00008246349, 3.58227925840, 206.18554843720],
      [0.00008768704, 3.63000308199, 949.17560896980],
      [0.00007368042, 5.08101194270, 735.87651353180],
      [0.00006263150, 0.02497628807, 213.29909543800],
      [0.00006114062, 4.51319998626, 1162.47470440780],
      [0.00004905396, 1.32084470588, 110.20632121940],
      [0.00005305285, 1.30671216791, 14.22709400160],
      [0.00005305441, 4.18625634012, 1052.26838318840],
      [0.00004647248, 4.69958103684, 3.93215326310],
      [0.00003045023, 4.31676431084, 426.59819087600],
      [0.00002609999, 1.56667394063, 846.08283475120],
      [0.00002028191, 1.06376530715, 3.18139373770],
      [0.00001764763, 2.14148655117, 1066.49547719000],
      [0.00001722972, 3.88036268267, 1265.56747862640],
      [0.00001920945, 0.97168196472, 639.89728631400],
      [0.00001633223, 3.58201833555, 515.46387109300],
      [0.00001431999, 4.29685556046, 625.67019231240],
      [0.00000973272, 4.09764549134, 95.97922721780]
    ],
    [
      [529.69096508814, 0.00000000000, 0.00000000000],
      [0.00489503243, 4.22082939470, 529.69096509460],
      [0.00228917222, 6.02646855621, 7.11354700080],
      [0.00030099479, 4.54540782858, 1059.38193018920],
      [0.00020720920, 5.45943156902, 522.57741809380],
      [0.00012103653, 0.16994816098, 536.80451209540],
      [0.00006067987, 4.42422292017, 103.09277421860],
      [0.00005433968, 3.98480737746, 419.48464387520],
      [0.00004237744, 5.89008707199, 14.22709400160]
    ],
    [
      [0.00047233601, 4.32148536482, 7.11354700080],
      [0.00030649436, 2.92977788700, 529.69096509460],
      [0.00014837605, 3.14159265359, 0.00000000000]
    ]
  ],
  [
    [
      [0.02268615702, 3.55852606721, 529.69096509460],
      [0.00109971634, 3.90809347197, 1059.38193018920],
      [0.00110090358, 0.00000000000, 0.00000000000],
      [0.00008101428, 3.60509572885, 522.57741809380],
      [0.00006043996, 4.25883108339, 1589.07289528380],
      [0.00006437782, 0.30627119215, 536.80451209540]
    ],
    [
      [0.00078203446, 1.52377859742, 529.69096509460]
    ]
  ],
  [
    [
      [5.20887429326, 0.00000000000, 0.00000000000],
      [0.25209327119, 3.49108639871, 529.69096509460],
      [0.00610599976, 3.84115365948, 1059.38193018920],
      [0.00282029458, 2.57419881293, 632.78373931320],
      [0.00187647346, 2.07590383214, 522.57741809380],
      [0.00086792905, 0.71001145545, 419.48464387520],
      [0.00072062974, 0.21465724607, 536.80451209540],
      [0.00065517248, 5.97995884790, 316.39186965660],
      [0.00029134542, 1.67759379655, 103.09277421860],
      [0.00030135335, 2.16132003734, 949.17560896980],
      [0.00023453271, 3.54023522184, 735.87651353180],
      [0.00022283743, 4.19362594399, 1589.07289528380],
      [0.00023947298, 0.27458037480, 7.11354700080],
      [0.00013032614, 2.96042965363, 1162.47470440780],
      [0.00009703360, 1.90669633585, 206.18554843720],
      [0.00012749023, 2.71550286592, 1052.26838318840],
      [0.00007057931, 2.18184839926, 1265.56747862640],
      [0.00006137703, 6.26418240033, 846.08283475120],
      [0.00002616976, 2.00994012876, 1581.95934828300]
    ],
    [
      [0.01271801520, 2.64937512894, 529.69096509460],
      [0.00061661816, 3.00076460387, 1059.38193018920],
      [0.00053443713, 3.89717383175, 522.57741809380],
      [0.00031185171, 4.88276958012, 536.80451209540],
      [0.00041390269, 0.00000000000, 0.00000000000]
    ]
  ]
],

    # Saturn
    [
  [
    [
      [0.87401354025, 0.00000000000, 0.00000000000],
      [0.11107659762, 3.96205090159, 213.29909543800],
      [0.01414150957, 4.58581516874, 7.11354700080],
      [0.00398379389, 0.52112032699, 206.18554843720],
      [0.00350769243, 3.30329907896, 426.59819087600],
      [0.00206816305, 0.24658372002, 103.09277421860],
      [0.00079271300, 3.84007056878, 220.41264243880],
      [0.00023990355, 4.66976924553, 110.20632121940],
      [0.00016573588, 0.43719228296, 419.48464387520],
      [0.00014906995, 5.76903183869, 316.39186965660],
      [0.00015820290, 0.93809155235, 632.78373931320],
      [0.00014609559, 1.56518472000, 3.93215326310],
      [0.00013160301, 4.44891291899, 14.22709400160],
      [0.00015053543, 2.71669915667, 639.89728631400],
      [0.00013005299, 5.98119023644, 11.04570026390],
      [0.00010725067, 3.12939523827, 202.25339517410],
      [0.00005863206, 0.23656938524, 529.69096509460],
      [0.00005227757, 4.20783365759, 3.18139373770],
      [0.00006126317, 1.76328667907, 277.03499374140],
      [0.00005019687, 3.17787728405, 433.71173787680],
      [0.00004592550, 0.61977744975, 199.07200143640],
      [0.00004005867, 2.24479718502, 63.73589830340],
      [0.00002953796, 0.98280366998, 95.97922721780],
      [0.00003873670, 3.22283226966, 138.51749687070],
      [0.00002461186, 2.03163875071, 735.87651353180],
      [0.00003269484, 0.77492638211, 949.17560896980],
      [0.00001758145, 3.26580109940, 522.57741809380],
      [0.00001640172, 5.50504453050, 846.08283475120],
      [0.00001391327, 4.02333150505, 323.50541665740],
      [0.00001580648, 4.37265307169, 309.27832265580],
      [0.00001123498, 2.83726798446, 415.55249061210],
      [0.00001017275, 3.71700135395, 227.52618943960],
      [0.00000848642, 3.19150170830, 209.36694217490]
    ],
    [
      [213.29909521690, 0.00000000000, 0.00000000000],
      [0.01297370862, 1.82834923978, 213.29909543800],
      [0.00564345393, 2.88499717272, 7.11354700080],
      [0.00093734369, 1.06311793502, 426.59819087600],
      [0.00107674962, 2.27769131009, 206.18554843720],
      [0.00040244455, 2.04108104671, 220.41264243880],
      [0.00019941774, 1.27954390470, 103.09277421860],
      [0.00010511678, 2.74880342130, 14.22709400160],
      [0.00006416106, 0.38238295041, 639.89728631400],
      [0.00004848994, 2.43037610229, 419.48464387520],
      [0.00004056892, 2.92133209468, 110.20632121940],
      [0.00003768635, 3.64965330780, 3.93215326310]
    ],
    [
      [0.00116441330, 1.17988132879, 7.11354700080],
      [0.00091841837, 0.07325195840, 213.29909543800],
      [0.00036661728, 0.00000000000, 0.00000000000],
      [0.00015274496, 4.06493179167, 206.18554843720]
    ]
  ],
  [
    [
      [0.04330678039, 3.60284428399, 213.29909543800],
      [0.00240348302, 2.85238489373, 426.59819087600],
      [0.00084745939, 0.00000000000, 0.00000000000],
      [0.00030863357, 3.48441504555, 220.41264243880],
      [0.00034116062, 0.57297307557, 206.18554843720],
      [0.00014734070, 2.11846596715, 639.89728631400],
      [0.00009916667, 5.79003188904, 419.48464387520],
      [0.00006993564, 4.73604689720, 7.11354700080],
      [0.00004807588, 5.43305312061, 316.39186965660]
    ],
    [
      [0.00198927992, 4.93901017903, 213.29909543800],
      [0.00036947916, 3.14159265359, 0.00000000000],
      [0.00017966989, 0.51979431110, 426.59819087600]
    ]
  ],
  [
    [
      [9.55758135486, 0.00000000000, 0.00000000000],
      [0.52921382865, 2.39226219573, 213.29909543800],
      [0.01873679867, 5.23549604660, 206.18554843720],
      [0.01464663929, 1.64763042902, 426.59819087600],
      [0.00821891141, 5.93520042303, 316.39186965660],
      [0.00547506923, 5.01532618980, 103.09277421860],
      [0.00371684650, 2.27114821115, 220.41264243880],
      [0.00361778765, 3.13904301847, 7.11354700080],
      [0.00140617506, 5.70406606781, 632.78373931320],
      [0.00108974848, 3.29313390175, 110.20632121940],
      [0.00069006962, 5.94099540992, 419.48464387520],
      [0.00061053367, 0.94037691801, 639.89728631400],
      [0.00048913294, 1.55733638681, 202.25339517410],
      [0.00034143772, 0.19519102597, 277.03499374140],
      [0.00032401773, 5.47084567016, 949.17560896980],
      [0.00020936596, 0.46349251129, 735.87651353180],
      [0.00009796004, 5.20477537945, 1265.56747862640],
      [0.00011993338, 5.98050967385, 846.08283475120],
      [0.00020839300, 1.52102476129, 433.71173787680],
      [0.00015298404, 3.05943814940, 529.69096509460],
      [0.00006465823, 0.17732249942, 1052.26838318840],
      [0.00011380257, 1.73105427040, 522.57741809380],
      [0.00003419618, 4.94550542171, 1581.95934828300]
    ],
    [
      [0.06182981340, 0.25843511480, 213.29909543800],
      [0.00506577242, 0.71114625261, 206.18554843720],
      [0.00341394029, 5.79635741658, 426.59819087600],
      [0.00188491195, 0.47215589652, 220.41264243880],
      [0.00186261486, 3.14159265359, 0.00000000000],
      [0.00143891146, 1.40744822888, 7.11354700080]
    ],
    [
      [0.00436902572, 4.78671677509, 213.29909543800]
    ]
  ]
],

    # Uranus
    [
  [
    [
      [5.48129294297, 0.00000000000, 0.00000000000],
      [0.09260408234, 0.89106421507, 74.78159856730],
      [0.01504247898, 3.62719260920, 1.48447270830],
      [0.00365981674, 1.89962179044, 73.29712585900],
      [0.00272328168, 3.35823706307, 149.56319713460],
      [0.00070328461, 5.39254450063, 63.73589830340],
      [0.00068892678, 6.09292483287, 76.26607127560],
      [0.00061998615, 2.26952066061, 2.96894541660],
      [0.00061950719, 2.85098872691, 11.04570026390],
      [0.00026468770, 3.14152083966, 71.81265315070],
      [0.00025710476, 6.11379840493, 454.90936652730],
      [0.00021078850, 4.36059339067, 148.07872442630],
      [0.00017818647, 1.74436930289, 36.64856292950],
      [0.00014613507, 4.73732166022, 3.93215326310],
      [0.00011162509, 5.82681796350, 224.34479570190],
      [0.00010997910, 0.48865004018, 138.51749687070],
      [0.00009527478, 2.95516862826, 35.16409022120],
      [0.00007545601, 5.23626582400, 109.94568878850],
      [0.00004220241, 3.23328220918, 70.84944530420],
      [0.00004051900, 2.27755017300, 151.04766984290],
      [0.00003354596, 1.06549007380, 4.45341812490],
      [0.00002926718, 4.62903718891, 9.56122755560],
      [0.00003490340, 5.48306144511, 146.59425171800],
      [0.00003144069, 4.75199570434, 77.75054398390],
      [0.00002922333, 5.35235361027, 85.82729883120],
      [0.00002272788, 4.36600400036, 70.32818044240],
      [0.00002051219, 1.51773566586, 0.11187458460],
      [0.00002148602, 0.60745949945, 38.13303563780],
      [0.00001991643, 4.92437588682, 277.03499374140],
      [0.00001376226, 2.04283539351, 65.22037101170],
      [0.00001666902, 3.62744066769, 380.12776796000],
      [0.00001284107, 3.11347961505, 202.25339517410],
      [0.00001150429, 0.93343589092, 3.18139373770],
      [0.00001533221, 2.58594681212, 52.69019803950],
      [0.00001281604, 0.54271272721, 222.86032299360],
      [0.00001372139, 4.19641530878, 111.43016149680],
      [0.00001221029, 0.19900650030, 108.46121608020],
      [0.00000946181, 1.19253165736, 127.47179660680],
      [0.00001150989, 4.17898916639, 33.67961751290]
    ],
    [
      [74.78159860910, 0.00000000000, 0.00000000000],
      [0.00154332863, 5.24158770553, 74.78159856730],
      [0.00024456474, 1.71260334156, 1.48447270830],
      [0.00009258442, 0.42829732350, 11.04570026390],
      [0.00008265977, 1.50218091379, 63.73589830340],
      [0.00009150160, 1.41213765216, 149.56319713460]
    ]
  ],
  [
    [
      [0.01346277648, 2.61877810547, 74.78159856730],
      [0.00062341400, 5.08111189648, 149.56319713460],
      [0.00061601196, 3.14159265359, 0.00000000000],
      [0.00009963722, 1.61603805646, 76.26607127560],
      [0.00009926160, 0.57630380333, 73.29712585900]
    ],
    [
      [0.00034101978, 0.01321929936, 74.78159856730]
    ]
  ],
  [
    [
      [19.21264847206, 0.00000000000, 0.00000000000],
      [0.88784984413, 5.60377527014, 74.78159856730],
      [0.03440836062, 0.32836099706, 73.29712585900],
      [0.02055653860, 1.78295159330, 149.56319713460],
      [0.00649322410, 4.52247285911, 76.26607127560],
      [0.00602247865, 3.86003823674, 63.73589830340],
      [0.00496404167, 1.40139935333, 454.90936652730],
      [0.00338525369, 1.58002770318, 138.51749687070],
      [0.00243509114, 1.57086606044, 71.81265315070],
      [0.00190522303, 1.99809394714, 1.48447270830],
      [0.00161858838, 2.79137786799, 148.07872442630],
      [0.00143706183, 1.38368544947, 11.04570026390],
      [0.00093192405, 0.17437220467, 36.64856292950],
      [0.00071424548, 4.24509236074, 224.34479570190],
      [0.00089806014, 3.66105364565, 109.94568878850],
      [0.00039009723, 1.66971401684, 70.84944530420],
      [0.00046677296, 1.39976401694, 35.16409022120],
      [0.00039025624, 3.36234773834, 277.03499374140],
      [0.00036755274, 3.88649278513, 146.59425171800],
      [0.00030348723, 0.70100838798, 151.04766984290],
      [0.00029156413, 3.18056336700, 77.75054398390],
      [0.00022637073, 0.72518687029, 529.69096509460],
      [0.00011959076, 1.75043392140, 984.60033162190],
      [0.00025620756, 5.25656086672, 380.12776796000]
    ],
    [
      [0.01479896629, 3.67205697578, 74.78159856730]
    ]
  ]
],

    # Neptune
    [
  [
    [
      [5.31188633046, 0.00000000000, 0.00000000000],
      [0.01798475530, 2.90101273890, 38.13303563780],
      [0.01019727652, 0.48580922867, 1.48447270830],
      [0.00124531845, 4.83008090676, 36.64856292950],
      [0.00042064466, 5.41054993053, 2.96894541660],
      [0.00037714584, 6.09221808686, 35.16409022120],
      [0.00033784738, 1.24488874087, 76.26607127560],
      [0.00016482741, 0.00007727998, 491.55792945680],
      [0.00009198584, 4.93747051954, 39.61750834610],
      [0.00008994250, 0.27462171806, 175.16605980020]
    ],
    [
      [38.13303563957, 0.00000000000, 0.00000000000],
      [0.00016604172, 4.86323329249, 1.48447270830],
      [0.00015744045, 2.27887427527, 38.13303563780]
    ]
  ],
  [
    [
      [0.03088622933, 1.44104372644, 38.13303563780],
      [0.00027780087, 5.91271884599, 76.26607127560],
      [0.00027623609, 0.00000000000, 0.00000000000],
      [0.00015355489, 2.52123799551, 36.64856292950],
      [0.00015448133, 3.50877079215, 39.61750834610]
    ]
  ],
  [
    [
      [30.07013205828, 0.00000000000, 0.00000000000],
      [0.27062259632, 1.32999459377, 38.13303563780],
      [0.01691764014, 3.25186135653, 36.64856292950],
      [0.00807830553, 5.18592878704, 1.48447270830],
      [0.00537760510, 4.52113935896, 35.16409022120],
      [0.00495725141, 1.57105641650, 491.55792945680],
      [0.00274571975, 1.84552258866, 175.16605980020],
      [0.00012012320, 1.92059384991, 1021.24889455140],
      [0.00121801746, 5.79754470298, 76.26607127560],
      [0.00100896068, 0.37702724930, 73.29712585900],
      [0.00135134092, 3.37220609835, 39.61750834610],
      [0.00007571796, 1.07149207335, 388.46515523820]
    ]
  ]
],
]

def _VsopFormula(formula, t, clamp_angle):
    tpower = 1.0
    coord = 0.0
    for series in formula:
        incr = tpower * sum(A * math.cos(B + C*t) for (A, B, C) in series)
        if clamp_angle:
            # Longitude angles can be hundreds of radians.
            # Improve precision by keeping each increment within [-2*pi, +2*pi].
            incr = math.fmod(incr, _PI2)
        coord += incr
        tpower *= t
    return coord

def _VsopDeriv(formula, t):
    tpower = 1      # t**s
    dpower = 0      # t**(s-1)
    deriv = 0
    s = 0
    for series in formula:
        sin_sum = 0
        cos_sum = 0
        for (ampl, phas, freq) in series:
            angle = phas + (t * freq)
            sin_sum += ampl * freq * math.sin(angle)
            if s > 0:
                cos_sum += ampl * math.cos(angle)
        deriv += (s * dpower * cos_sum) - (tpower * sin_sum)
        dpower = tpower
        tpower *= t
        s += 1
    return deriv

_DAYS_PER_MILLENNIUM = 365250.0
_LON_INDEX = 0
_LAT_INDEX = 1
_RAD_INDEX = 2


def _VsopRotate(eclip):
    # Convert ecliptic cartesian coordinates to equatorial cartesian coordinates.
    x = eclip.x + 0.000000440360*eclip.y - 0.000000190919*eclip.z
    y = -0.000000479966*eclip.x + 0.917482137087*eclip.y - 0.397776982902*eclip.z
    z = 0.397776982902*eclip.y + 0.917482137087*eclip.z
    return _TerseVector(x, y, z)

def _VsopSphereToRect(lon, lat, rad):
    # Convert spherical coordinates to cartesian coordinates.
    r_coslat = rad * math.cos(lat)
    return _TerseVector(
        r_coslat * math.cos(lon),
        r_coslat * math.sin(lon),
        rad * math.sin(lat)
    )

def _CalcVsop(model, time):
    t = time.tt / _DAYS_PER_MILLENNIUM
    lon = _VsopFormula(model[0], t, True)
    lat = _VsopFormula(model[1], t, False)
    rad = _VsopFormula(model[2], t, False)
    eclip = _VsopSphereToRect(lon, lat, rad)
    return _VsopRotate(eclip).ToAstroVector(time)

class _body_state_t:
    def __init__(self, tt, r, v):
        self.tt  = tt
        self.r = r
        self.v = v

    def clone(self):
        '''Make a copy of this body state.'''
        return _body_state_t(self.tt, self.r.clone(), self.v.clone())

    def __sub__(self, other):
        return _body_state_t(self.tt, self.r - other.r, self.v - other.v)

def _CalcVsopPosVel(model, tt):
    t = tt / _DAYS_PER_MILLENNIUM

    lon = _VsopFormula(model[0], t, True)
    lat = _VsopFormula(model[1], t, False)
    rad = _VsopFormula(model[2], t, False)

    (dlon_dt, dlat_dt, drad_dt) = [_VsopDeriv(formula, t) for formula in model]

    # Use spherical coords and spherical derivatives to calculate
    # the velocity vector in rectangular coordinates.

    coslon = math.cos(lon)
    sinlon = math.sin(lon)
    coslat = math.cos(lat)
    sinlat = math.sin(lat)

    vx = (
        + (drad_dt * coslat * coslon)
        - (rad * sinlat * coslon * dlat_dt)
        - (rad * coslat * sinlon * dlon_dt)
    )

    vy = (
        + (drad_dt * coslat * sinlon)
        - (rad * sinlat * sinlon * dlat_dt)
        + (rad * coslat * coslon * dlon_dt)
    )

    vz = (
        + (drad_dt * sinlat)
        + (rad * coslat * dlat_dt)
    )

    eclip_pos = _VsopSphereToRect(lon, lat, rad)

    # Convert speed units from [AU/millennium] to [AU/day].
    eclip_vel = _TerseVector(vx, vy, vz) / _DAYS_PER_MILLENNIUM

    # Rotate the vectors from ecliptic to equatorial coordinates.
    equ_pos = _VsopRotate(eclip_pos)
    equ_vel = _VsopRotate(eclip_vel)
    return _body_state_t(tt, equ_pos, equ_vel)


def _AdjustBarycenter(ssb, time, body, pmass):
    shift = pmass / (pmass + _SUN_GM)
    planet = _CalcVsop(_vsop[body.value], time)
    ssb.x += shift * planet.x
    ssb.y += shift * planet.y
    ssb.z += shift * planet.z


def _CalcSolarSystemBarycenter(time):
    ssb = Vector(0.0, 0.0, 0.0, time)
    _AdjustBarycenter(ssb, time, Body.Jupiter, _JUPITER_GM)
    _AdjustBarycenter(ssb, time, Body.Saturn,  _SATURN_GM)
    _AdjustBarycenter(ssb, time, Body.Uranus,  _URANUS_GM)
    _AdjustBarycenter(ssb, time, Body.Neptune, _NEPTUNE_GM)
    return ssb


def _VsopHelioDistance(model, time):
    # The caller only wants to know the distance between the planet and the Sun.
    # So we only need to calculate the radial component of the spherical coordinates.
    # There is no need to translate coordinates.
    return _VsopFormula(model[2], time.tt / _DAYS_PER_MILLENNIUM, False)

def _CalcEarth(time):
    return _CalcVsop(_vsop[Body.Earth.value], time)

# END VSOP
#----------------------------------------------------------------------------
# BEGIN Pluto Integrator

_PLUTO_NUM_STATES = 51
_PLUTO_TIME_STEP  = 29200
_PLUTO_DT         = 146
_PLUTO_NSTEPS     = 201

_PlutoStateTable = [
    [ -730000.0, [-26.118207232108, -14.376168177825,   3.384402515299], [ 1.6339372163656e-03, -2.7861699588508e-03, -1.3585880229445e-03]]
,   [ -700800.0, [ 41.974905202127,  -0.448502952929, -12.770351505989], [ 7.3458569351457e-04,  2.2785014891658e-03,  4.8619778602049e-04]]
,   [ -671600.0, [ 14.706930780744,  44.269110540027,   9.353698474772], [-2.1000147999800e-03,  2.2295915939915e-04,  7.0143443551414e-04]]
,   [ -642400.0, [-29.441003929957,  -6.430161530570,   6.858481011305], [ 8.4495803960544e-04, -3.0783914758711e-03, -1.2106305981192e-03]]
,   [ -613200.0, [ 39.444396946234,  -6.557989760571, -13.913760296463], [ 1.1480029005873e-03,  2.2400006880665e-03,  3.5168075922288e-04]]
,   [ -584000.0, [ 20.230380950700,  43.266966657189,   7.382966091923], [-1.9754081700585e-03,  5.3457141292226e-04,  7.5929169129793e-04]]
,   [ -554800.0, [-30.658325364620,   2.093818874552,   9.880531138071], [ 6.1010603013347e-05, -3.1326500935382e-03, -9.9346125151067e-04]]
,   [ -525600.0, [ 35.737703251673, -12.587706024764, -14.677847247563], [ 1.5802939375649e-03,  2.1347678412429e-03,  1.9074436384343e-04]]
,   [ -496400.0, [ 25.466295188546,  41.367478338417,   5.216476873382], [-1.8054401046468e-03,  8.3283083599510e-04,  8.0260156912107e-04]]
,   [ -467200.0, [-29.847174904071,  10.636426313081,  12.297904180106], [-6.3257063052907e-04, -2.9969577578221e-03, -7.4476074151596e-04]]
,   [ -438000.0, [ 30.774692107687, -18.236637015304, -14.945535879896], [ 2.0113162005465e-03,  1.9353827024189e-03, -2.0937793168297e-06]]
,   [ -408800.0, [ 30.243153324028,  38.656267888503,   2.938501750218], [-1.6052508674468e-03,  1.1183495337525e-03,  8.3333973416824e-04]]
,   [ -379600.0, [-27.288984772533,  18.643162147874,  14.023633623329], [-1.1856388898191e-03, -2.7170609282181e-03, -4.9015526126399e-04]]
,   [ -350400.0, [ 24.519605196774, -23.245756064727, -14.626862367368], [ 2.4322321483154e-03,  1.6062008146048e-03, -2.3369181613312e-04]]
,   [ -321200.0, [ 34.505274805875,  35.125338586954,   0.557361475637], [-1.3824391637782e-03,  1.3833397561817e-03,  8.4823598806262e-04]]
,   [ -292000.0, [-23.275363915119,  25.818514298769,  15.055381588598], [-1.6062295460975e-03, -2.3395961498533e-03, -2.4377362639479e-04]]
,   [ -262800.0, [ 17.050384798092, -27.180376290126, -13.608963321694], [ 2.8175521080578e-03,  1.1358749093955e-03, -4.9548725258825e-04]]
,   [ -233600.0, [ 38.093671910285,  30.880588383337,  -1.843688067413], [-1.1317697153459e-03,  1.6128814698472e-03,  8.4177586176055e-04]]
,   [ -204400.0, [-18.197852930878,  31.932869934309,  15.438294826279], [-1.9117272501813e-03, -1.9146495909842e-03, -1.9657304369835e-05]]
,   [ -175200.0, [  8.528924039997, -29.618422200048, -11.805400994258], [ 3.1034370787005e-03,  5.1393633292430e-04, -7.7293066202546e-04]]
,   [ -146000.0, [ 40.946857258640,  25.904973592021,  -4.256336240499], [-8.3652705194051e-04,  1.8129497136404e-03,  8.1564228273060e-04]]
,   [ -116800.0, [-12.326958895325,  36.881883446292,  15.217158258711], [-2.1166103705038e-03, -1.4814420035990e-03,  1.7401209844705e-04]]
,   [  -87600.0, [ -0.633258375909, -30.018759794709,  -9.171932874950], [ 3.2016994581737e-03, -2.5279858672148e-04, -1.0411088271861e-03]]
,   [  -58400.0, [ 42.936048423883,  20.344685584452,  -6.588027007912], [-5.0525450073192e-04,  1.9910074335507e-03,  7.7440196540269e-04]]
,   [  -29200.0, [ -5.975910552974,  40.611809958460,  14.470131723673], [-2.2184202156107e-03, -1.0562361130164e-03,  3.3652250216211e-04]]
,   [       0.0, [ -9.875369580774, -27.978926224737,  -5.753711824704], [ 3.0287533248818e-03, -1.1276087003636e-03, -1.2651326732361e-03]]
,   [   29200.0, [ 43.958831986165,  14.214147973292,  -8.808306227163], [-1.4717608981871e-04,  2.1404187242141e-03,  7.1486567806614e-04]]
,   [   58400.0, [  0.678136763520,  43.094461639362,  13.243238780721], [-2.2358226110718e-03, -6.3233636090933e-04,  4.7664798895648e-04]]
,   [   87600.0, [-18.282602096834, -23.305039586660,  -1.766620508028], [ 2.5567245263557e-03, -1.9902940754171e-03, -1.3943491701082e-03]]
,   [  116800.0, [ 43.873338744526,   7.700705617215, -10.814273666425], [ 2.3174803055677e-04,  2.2402163127924e-03,  6.2988756452032e-04]]
,   [  146000.0, [  7.392949027906,  44.382678951534,  11.629500214854], [-2.1932815453830e-03, -2.1751799585364e-04,  5.9556516201114e-04]]
,   [  175200.0, [-24.981690229261, -16.204012851426,   2.466457544298], [ 1.8193989149580e-03, -2.6765419531201e-03, -1.3848283502247e-03]]
,   [  204400.0, [ 42.530187039511,   0.845935508021, -12.554907527683], [ 6.5059779150669e-04,  2.2725657282262e-03,  5.1133743202822e-04]]
,   [  233600.0, [ 13.999526486822,  44.462363044894,   9.669418486465], [-2.1079296569252e-03,  1.7533423831993e-04,  6.9128485798076e-04]]
,   [  262800.0, [-29.184024803031,  -7.371243995762,   6.493275957928], [ 9.3581363109681e-04, -3.0610357109184e-03, -1.2364201089345e-03]]
,   [  292000.0, [ 39.831980671753,  -6.078405766765, -13.909815358656], [ 1.1117769689167e-03,  2.2362097830152e-03,  3.6230548231153e-04]]
,   [  321200.0, [ 20.294955108476,  43.417190420251,   7.450091985932], [-1.9742157451535e-03,  5.3102050468554e-04,  7.5938408813008e-04]]
,   [  350400.0, [-30.669992302160,   2.318743558955,   9.973480913858], [ 4.5605107450676e-05, -3.1308219926928e-03, -9.9066533301924e-04]]
,   [  379600.0, [ 35.626122155983, -12.897647509224, -14.777586508444], [ 1.6015684949743e-03,  2.1171931182284e-03,  1.8002516202204e-04]]
,   [  408800.0, [ 26.133186148561,  41.232139187599,   5.006401326220], [-1.7857704419579e-03,  8.6046232702817e-04,  8.0614690298954e-04]]
,   [  438000.0, [-29.576740229230,  11.863535943587,  12.631323039872], [-7.2292830060955e-04, -2.9587820140709e-03, -7.0824296450300e-04]]
,   [  467200.0, [ 29.910805787391, -19.159019294000, -15.013363865194], [ 2.0871080437997e-03,  1.8848372554514e-03, -3.8528655083926e-05]]
,   [  496400.0, [ 31.375957451819,  38.050372720763,   2.433138343754], [-1.5546055556611e-03,  1.1699815465629e-03,  8.3565439266001e-04]]
,   [  525600.0, [-26.360071336928,  20.662505904952,  14.414696258958], [-1.3142373118349e-03, -2.6236647854842e-03, -4.2542017598193e-04]]
,   [  554800.0, [ 22.599441488648, -24.508879898306, -14.484045731468], [ 2.5454108304806e-03,  1.4917058755191e-03, -3.0243665086079e-04]]
,   [  584000.0, [ 35.877864013014,  33.894226366071,  -0.224524636277], [-1.2941245730845e-03,  1.4560427668319e-03,  8.4762160640137e-04]]
,   [  613200.0, [-21.538149762417,  28.204068269761,  15.321973799534], [-1.7312117409010e-03, -2.1939631314577e-03, -1.6316913275180e-04]]
,   [  642400.0, [ 13.971521374415, -28.339941764789, -13.083792871886], [ 2.9334630526035e-03,  9.1860931752944e-04, -5.9939422488627e-04]]
,   [  671600.0, [ 39.526942044143,  28.939897360110,  -2.872799527539], [-1.0068481658095e-03,  1.7021132888090e-03,  8.3578230511981e-04]]
,   [  700800.0, [-15.576200701394,  34.399412961275,  15.466033737854], [-2.0098814612884e-03, -1.7191109825989e-03,  7.0414782780416e-05]]
,   [  730000.0, [  4.243252837090, -30.118201690825, -10.707441231349], [ 3.1725847067411e-03,  1.6098461202270e-04, -9.0672150593868e-04]]
]


class _TerseVector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def clone(self):
        '''Create a copy of this vector.'''
        return _TerseVector(self.x, self.y, self.z)

    @staticmethod
    def zero():
        '''Return a zero vector.'''
        return _TerseVector(0.0, 0.0, 0.0)

    def ToAstroVector(self, time):
        '''Convert _TerseVector object to Vector object.'''
        return Vector(self.x, self.y, self.z, time)

    def quadrature(self):
        '''Return magnitude squared of this vector.'''
        return self.x**2 + self.y**2 + self.z**2

    def mean(self, other):
        '''Return the average of this vector and another vector.'''
        return _TerseVector((self.x + other.x)/2.0, (self.y + other.y)/2.0, (self.z + other.z)/2.0)

    def __add__(self, other):
        return _TerseVector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return _TerseVector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar):
        return _TerseVector(scalar * self.x, scalar * self.y, scalar * self.z)

    def __rmul__(self, scalar):
        return _TerseVector(scalar * self.x, scalar * self.y, scalar * self.z)

    def __truediv__(self, scalar):
        return _TerseVector(self.x / scalar, self.y / scalar, self.z / scalar)


def _BodyStateFromTable(entry):
    [ tt, [rx, ry, rz], [vx, vy, vz] ] = entry
    return _body_state_t(tt, _TerseVector(rx, ry, rz), _TerseVector(vx, vy, vz))


def _AdjustBarycenterPosVel(ssb, tt, body, planet_gm):
    shift = planet_gm / (planet_gm + _SUN_GM)
    planet = _CalcVsopPosVel(_vsop[body.value], tt)
    ssb.r += shift * planet.r
    ssb.v += shift * planet.v
    return planet

def _AccelerationIncrement(small_pos, gm, major_pos):
    delta = major_pos - small_pos
    r2 = delta.quadrature()
    return (gm / (r2 * math.sqrt(r2))) * delta


class _major_bodies_t:
    def __init__(self, tt):
        # Accumulate the Solar System Barycenter position.
        ssb = _body_state_t(tt, _TerseVector(0,0,0), _TerseVector(0,0,0))
        # Calculate the position and velocity vectors of the 4 major planets.
        self.Jupiter = _AdjustBarycenterPosVel(ssb, tt, Body.Jupiter, _JUPITER_GM)
        self.Saturn  = _AdjustBarycenterPosVel(ssb, tt, Body.Saturn,  _SATURN_GM)
        self.Uranus  = _AdjustBarycenterPosVel(ssb, tt, Body.Uranus,  _URANUS_GM)
        self.Neptune = _AdjustBarycenterPosVel(ssb, tt, Body.Neptune, _NEPTUNE_GM)
        # Convert the planets' vectors from heliocentric to barycentric.
        self.Jupiter.r -= ssb.r
        self.Jupiter.v -= ssb.v
        self.Saturn.r  -= ssb.r
        self.Saturn.v  -= ssb.v
        self.Uranus.r  -= ssb.r
        self.Uranus.v  -= ssb.v
        self.Neptune.r -= ssb.r
        self.Neptune.v -= ssb.v
        # Convert heliocentric SSB to barycentric Sun.
        self.Sun = _body_state_t(tt, -1*ssb.r, -1*ssb.v)

    def Acceleration(self, pos):
        '''Use barycentric coordinates of the Sun and major planets to calculate
        the gravitational acceleration vector experienced at location 'pos'.'''
        acc  = _AccelerationIncrement(pos, _SUN_GM,     self.Sun.r)
        acc += _AccelerationIncrement(pos, _JUPITER_GM, self.Jupiter.r)
        acc += _AccelerationIncrement(pos, _SATURN_GM,  self.Saturn.r)
        acc += _AccelerationIncrement(pos, _URANUS_GM,  self.Uranus.r)
        acc += _AccelerationIncrement(pos, _NEPTUNE_GM, self.Neptune.r)
        return acc


class _body_grav_calc_t:
    def __init__(self, tt, r, v, a):
        self.tt = tt    # J2000 terrestrial time [days]
        self.r = r      # position [au]
        self.v = v      # velocity [au/day]
        self.a = a      # acceleration [au/day^2]

    def clone(self):
        '''Creates a copy of this gravity simulation state.'''
        return _body_grav_calc_t(self.tt, self.r.clone(), self.v.clone(), self.a.clone())


class _grav_sim_t:
    def __init__(self, bary, grav):
        self.bary = bary
        self.grav = grav


def _UpdatePosition(dt, r, v, a):
    return _TerseVector(
        r.x + dt*(v.x + dt*a.x/2.0),
        r.y + dt*(v.y + dt*a.y/2.0),
        r.z + dt*(v.z + dt*a.z/2.0)
    )

def _UpdateVelocity(dt, v, a):
    return _TerseVector(
        v.x + dt*a.x,
        v.y + dt*a.y,
        v.z + dt*a.z
    )


def _GravSim(tt2, calc1):
    dt = tt2 - calc1.tt

    # Calculate where the major bodies (Sun, Jupiter...Neptune) will be at tt2.
    bary2 = _major_bodies_t(tt2)

    # Estimate position of small body as if current acceleration applies across the whole time interval.
    approx_pos = _UpdatePosition(dt, calc1.r, calc1.v, calc1.a)

    # Calculate the average acceleration of the endpoints.
    # This becomes our estimate of the mean effective acceleration over the whole interval.
    mean_acc = bary2.Acceleration(approx_pos).mean(calc1.a)

    # Refine the estimates of [pos, vel, acc] at tt2 using the mean acceleration.
    pos = _UpdatePosition(dt, calc1.r, calc1.v, mean_acc)
    vel = calc1.v + dt*mean_acc
    acc = bary2.Acceleration(pos)
    grav = _body_grav_calc_t(tt2, pos, vel, acc)
    return _grav_sim_t(bary2, grav)


_pluto_cache = [None] * (_PLUTO_NUM_STATES - 1)


def _ClampIndex(frac, nsteps):
    index = math.floor(frac)
    if index < 0:
        return 0
    if index >= nsteps:
        return nsteps-1
    return index


def _GravFromState(entry):
    state = _BodyStateFromTable(entry)
    bary = _major_bodies_t(state.tt)
    r = state.r + bary.Sun.r
    v = state.v + bary.Sun.v
    a = bary.Acceleration(r)
    grav = _body_grav_calc_t(state.tt, r, v, a)
    return _grav_sim_t(bary, grav)


def _GetSegment(cache, tt):
    if (tt < _PlutoStateTable[0][0]) or (tt > _PlutoStateTable[_PLUTO_NUM_STATES-1][0]):
        # Don't bother calculating a segment. Let the caller crawl backward/forward to this time.
        return None

    seg_index = _ClampIndex((tt - _PlutoStateTable[0][0]) / _PLUTO_TIME_STEP, _PLUTO_NUM_STATES-1)
    if cache[seg_index] is None:
        seg = cache[seg_index] = [None] * _PLUTO_NSTEPS

        # Each endpoint is exact.
        seg[0] = _GravFromState(_PlutoStateTable[seg_index]).grav
        seg[_PLUTO_NSTEPS-1] = _GravFromState(_PlutoStateTable[seg_index + 1]).grav

        # Simulate forwards from the lower time bound.
        step_tt = seg[0].tt
        i = 1
        while i < _PLUTO_NSTEPS-1:
            step_tt += _PLUTO_DT
            seg[i] = _GravSim(step_tt, seg[i-1]).grav
            i += 1

        # Simulate backwards from the upper time bound.
        step_tt = seg[_PLUTO_NSTEPS-1].tt
        reverse = [None] * _PLUTO_NSTEPS
        reverse[_PLUTO_NSTEPS-1] = seg[_PLUTO_NSTEPS-1]
        i = _PLUTO_NSTEPS - 2
        while i > 0:
            step_tt -= _PLUTO_DT
            reverse[i] = _GravSim(step_tt, reverse[i+1]).grav
            i -= 1

        # Fade-mix the two series so that there are no discontinuities.
        i = _PLUTO_NSTEPS - 2
        while i > 0:
            ramp = i / (_PLUTO_NSTEPS-1)
            seg[i].r = seg[i].r*(1 - ramp) + reverse[i].r*ramp
            seg[i].v = seg[i].v*(1 - ramp) + reverse[i].v*ramp
            seg[i].a = seg[i].a*(1 - ramp) + reverse[i].a*ramp
            i -= 1

    return cache[seg_index]


def _CalcPlutoOneWay(entry, target_tt, dt):
    sim = _GravFromState(entry)
    n = math.ceil((target_tt - sim.grav.tt) / dt)
    for i in range(n):
        sim = _GravSim(target_tt if (i+1 == n) else (sim.grav.tt + dt), sim.grav)
    return sim


def _CalcPluto(time, helio):
    bary = None
    seg = _GetSegment(_pluto_cache, time.tt)
    if seg is None:
        # The target time is outside the year range 0000..4000.
        # Calculate it by crawling backward from 0000 or forward from 4000.
        # FIXFIXFIX - This is super slow. Could optimize this with extra caching if needed.
        if time.tt < _PlutoStateTable[0][0]:
            sim = _CalcPlutoOneWay(_PlutoStateTable[0], time.tt, -_PLUTO_DT)
        else:
            sim = _CalcPlutoOneWay(_PlutoStateTable[_PLUTO_NUM_STATES-1], time.tt, +_PLUTO_DT)
        r = sim.grav.r
        v = sim.grav.v
        bary = sim.bary
    else:
        left = _ClampIndex((time.tt - seg[0].tt) / _PLUTO_DT, _PLUTO_NSTEPS-1)
        s1 = seg[left]
        s2 = seg[left+1]

        # Find mean acceleration vector over the interval.
        acc = s1.a.mean(s2.a)

        # Use Newtonian mechanics to extrapolate away from t1 in the positive time direction.
        ra = _UpdatePosition(time.tt - s1.tt, s1.r, s1.v, acc)
        va = _UpdateVelocity(time.tt - s1.tt, s1.v, acc)

        # Use Newtonian mechanics to extrapolate away from t2 in the negative time direction.
        rb = _UpdatePosition(time.tt - s2.tt, s2.r, s2.v, acc)
        vb = _UpdateVelocity(time.tt - s2.tt, s2.v, acc)

        # Use fade in/out idea to blend the two position estimates.
        ramp = (time.tt - s1.tt)/_PLUTO_DT
        r = ra*(1 - ramp) + rb*ramp
        v = va*(1 - ramp) + vb*ramp

    if helio:
        # Convert barycentric vectors to heliocentric vectors.
        if bary is None:
            bary = _major_bodies_t(time.tt)
        r -= bary.Sun.r
        v -= bary.Sun.v

    return StateVector(r.x, r.y, r.z, v.x, v.y, v.z, time)


# END Pluto Integrator
#----------------------------------------------------------------------------
# BEGIN Jupiter Moons

_Rotation_JUP_EQJ = RotationMatrix([
    [  9.99432765338654e-01, -3.36771074697641e-02,  0.00000000000000e+00 ],
    [  3.03959428906285e-02,  9.02057912352809e-01,  4.30543388542295e-01 ],
    [ -1.44994559663353e-02, -4.30299169409101e-01,  9.02569881273754e-01 ]
])

_JupiterMoonModel = [
    # [0] Io
    [
         2.8248942843381399e-07,  1.4462132960212239e+00,  3.5515522861824000e+00, # mu, al0, al1
        [   # a
            [  0.0028210960212903,  0.0000000000000000e+00,  0.0000000000000000e+00 ]
        ],
        [   # l
            [ -0.0001925258348666,  4.9369589722644998e+00,  1.3584836583050000e-02 ],
            [ -0.0000970803596076,  4.3188796477322002e+00,  1.3034138432430000e-02 ],
            [ -0.0000898817416500,  1.9080016428616999e+00,  3.0506486715799999e-03 ],
            [ -0.0000553101050262,  1.4936156681568999e+00,  1.2938928911549999e-02 ]
        ],
        [   # z
            [  0.0041510849668155,  4.0899396355450000e+00, -1.2906864146660001e-02 ],
            [  0.0006260521444113,  1.4461888986270000e+00,  3.5515522949801999e+00 ],
            [  0.0000352747346169,  2.1256287034577999e+00,  1.2727416566999999e-04 ]
        ],
        [   # zeta
            [  0.0003142172466014,  2.7964219722923001e+00, -2.3150960980000000e-03 ],
            [  0.0000904169207946,  1.0477061879627001e+00, -5.6920638196000003e-04 ]
        ]
    ],

    # [1] Europa
    [
         2.8248327439289299e-07, -3.7352634374713622e-01,  1.7693227111234699e+00, # mu, al0, al1
        [   # a
            [  0.0044871037804314,  0.0000000000000000e+00,  0.0000000000000000e+00 ],
            [  0.0000004324367498,  1.8196456062910000e+00,  1.7822295777568000e+00 ]
        ],
        [   # l
            [  0.0008576433172936,  4.3188693178264002e+00,  1.3034138308049999e-02 ],
            [  0.0004549582875086,  1.4936531751079001e+00,  1.2938928819619999e-02 ],
            [  0.0003248939825174,  1.8196494533458001e+00,  1.7822295777568000e+00 ],
            [ -0.0003074250079334,  4.9377037005910998e+00,  1.3584832867240000e-02 ],
            [  0.0001982386144784,  1.9079869054759999e+00,  3.0510121286900001e-03 ],
            [  0.0001834063551804,  2.1402853388529000e+00,  1.4500978933800000e-03 ],
            [ -0.0001434383188452,  5.6222140366630002e+00,  8.9111478887838003e-01 ],
            [ -0.0000771939140944,  4.3002724372349999e+00,  2.6733443704265998e+00 ]
        ],
        [   # z
            [ -0.0093589104136341,  4.0899396509038999e+00, -1.2906864146660001e-02 ],
            [  0.0002988994545555,  5.9097265185595003e+00,  1.7693227079461999e+00 ],
            [  0.0002139036390350,  2.1256289300016000e+00,  1.2727418406999999e-04 ],
            [  0.0001980963564781,  2.7435168292649998e+00,  6.7797343008999997e-04 ],
            [  0.0001210388158965,  5.5839943711203004e+00,  3.2056614899999997e-05 ],
            [  0.0000837042048393,  1.6094538368039000e+00, -9.0402165808846002e-01 ],
            [  0.0000823525166369,  1.4461887708689001e+00,  3.5515522949801999e+00 ]
        ],
        [   # zeta
            [  0.0040404917832303,  1.0477063169425000e+00, -5.6920640539999997e-04 ],
            [  0.0002200421034564,  3.3368857864364001e+00, -1.2491307306999999e-04 ],
            [  0.0001662544744719,  2.4134862374710999e+00,  0.0000000000000000e+00 ],
            [  0.0000590282470983,  5.9719930968366004e+00, -3.0561602250000000e-05 ]
        ]
    ],

    # [2] Ganymede
    [
         2.8249818418472298e-07,  2.8740893911433479e-01,  8.7820792358932798e-01, # mu, al0, al1
        [   # a
            [  0.0071566594572575,  0.0000000000000000e+00,  0.0000000000000000e+00 ],
            [  0.0000013930299110,  1.1586745884981000e+00,  2.6733443704265998e+00 ]
        ],
        [   # l
            [  0.0002310797886226,  2.1402987195941998e+00,  1.4500978438400001e-03 ],
            [ -0.0001828635964118,  4.3188672736968003e+00,  1.3034138282630000e-02 ],
            [  0.0001512378778204,  4.9373102372298003e+00,  1.3584834812520000e-02 ],
            [ -0.0001163720969778,  4.3002659861490002e+00,  2.6733443704265998e+00 ],
            [ -0.0000955478069846,  1.4936612842567001e+00,  1.2938928798570001e-02 ],
            [  0.0000815246854464,  5.6222137132535002e+00,  8.9111478887838003e-01 ],
            [ -0.0000801219679602,  1.2995922951532000e+00,  1.0034433456728999e+00 ],
            [ -0.0000607017260182,  6.4978769669238001e-01,  5.0172167043264004e-01 ]
        ],
        [   # z
            [  0.0014289811307319,  2.1256295942738999e+00,  1.2727413029000001e-04 ],
            [  0.0007710931226760,  5.5836330003496002e+00,  3.2064341100000001e-05 ],
            [  0.0005925911780766,  4.0899396636447998e+00, -1.2906864146660001e-02 ],
            [  0.0002045597496146,  5.2713683670371996e+00, -1.2523544076106000e-01 ],
            [  0.0001785118648258,  2.8743156721063001e-01,  8.7820792442520001e-01 ],
            [  0.0001131999784893,  1.4462127277818000e+00,  3.5515522949801999e+00 ],
            [ -0.0000658778169210,  2.2702423990985001e+00, -1.7951364394536999e+00 ],
            [  0.0000497058888328,  5.9096792204858000e+00,  1.7693227129285001e+00 ]
        ],
        [   # zeta
            [  0.0015932721570848,  3.3368862796665000e+00, -1.2491307058000000e-04 ],
            [  0.0008533093128905,  2.4133881688166001e+00,  0.0000000000000000e+00 ],
            [  0.0003513347911037,  5.9720789850126996e+00, -3.0561017709999999e-05 ],
            [ -0.0001441929255483,  1.0477061764435001e+00, -5.6920632124000004e-04 ]
        ]
    ],

    # [3] Callisto
    [
         2.8249214488990899e-07, -3.6203412913757038e-01,  3.7648623343382798e-01, # mu, al0, al1
        [   # a
            [  0.0125879701715314,  0.0000000000000000e+00,  0.0000000000000000e+00 ],
            [  0.0000035952049470,  6.4965776007116005e-01,  5.0172168165034003e-01 ],
            [  0.0000027580210652,  1.8084235781510001e+00,  3.1750660413359002e+00 ]
        ],
        [   # l
            [  0.0005586040123824,  2.1404207189814999e+00,  1.4500979323100001e-03 ],
            [ -0.0003805813868176,  2.7358844897852999e+00,  2.9729650620000000e-05 ],
            [  0.0002205152863262,  6.4979652596399995e-01,  5.0172167243580001e-01 ],
            [  0.0001877895151158,  1.8084787604004999e+00,  3.1750660413359002e+00 ],
            [  0.0000766916975242,  6.2720114319754998e+00,  1.3928364636651001e+00 ],
            [  0.0000747056855106,  1.2995916202344000e+00,  1.0034433456728999e+00 ]
        ],
        [   # z
            [  0.0073755808467977,  5.5836071576083999e+00,  3.2065099140000001e-05 ],
            [  0.0002065924169942,  5.9209831565786004e+00,  3.7648624194703001e-01 ],
            [  0.0001589869764021,  2.8744006242622999e-01,  8.7820792442520001e-01 ],
            [ -0.0001561131605348,  2.1257397865089001e+00,  1.2727441285000001e-04 ],
            [  0.0001486043380971,  1.4462134301023000e+00,  3.5515522949801999e+00 ],
            [  0.0000635073108731,  5.9096803285953996e+00,  1.7693227129285001e+00 ],
            [  0.0000599351698525,  4.1125517584797997e+00, -2.7985797954588998e+00 ],
            [  0.0000540660842731,  5.5390350845569003e+00,  2.8683408228299999e-03 ],
            [ -0.0000489596900866,  4.6218149483337996e+00, -6.2695712529518999e-01 ]
        ],
        [   # zeta
            [  0.0038422977898495,  2.4133922085556998e+00,  0.0000000000000000e+00 ],
            [  0.0022453891791894,  5.9721736773277003e+00, -3.0561255249999997e-05 ],
            [ -0.0002604479450559,  3.3368746306408998e+00, -1.2491309972000001e-04 ],
            [  0.0000332112143230,  5.5604137742336999e+00,  2.9003768850700000e-03 ]
        ]
    ]
]

class JupiterMoonsInfo:
    """Holds the positions and velocities of Jupiter's major 4 moons.

    The #JupiterMoons function returns an object of this type
    to report position and velocity vectors for Jupiter's largest 4 moons
    Io, Europa, Ganymede, and Callisto. Each position vector is relative
    to the center of Jupiter. Both position and velocity are oriented in
    the EQJ system (that is, using Earth's equator at the J2000 epoch).
    The positions are expressed in astronomical units (AU),
    and the velocities in AU/day.

    Attributes
    ----------
    io : StateVector
        The position and velocity of Jupiter's moon Io.
    europa : StateVector
        The position and velocity of Jupiter's moon Europa.
    ganymede : StateVector
        The position and velocity of Jupiter's moon Ganymede.
    callisto : StateVector
        The position and velocity of Jupiter's moon Callisto.
    """
    def __init__(self, moon):
        self.io = moon[0]
        self.europa = moon[1]
        self.ganymede = moon[2]
        self.callisto = moon[3]

    def __repr__(self):
        return 'JupiterMoonsInfo(io={}, europa={}, ganymede={}, callisto={})'.format(
            repr(self.io),
            repr(self.europa),
            repr(self.ganymede),
            repr(self.callisto)
        )


def _JupiterMoon_elem2pv(time, mu, A, AL, K, H, Q, P):
    # Translation of FORTRAN subroutine ELEM2PV from:
    # https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/
    AN = math.sqrt(mu / (A*A*A))
    EE = AL + K*math.sin(AL) - H*math.cos(AL)
    DE = 1
    while abs(DE) >= 1.0e-12:
        CE = math.cos(EE)
        SE = math.sin(EE)
        DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE)
        EE += DE
    CE = math.cos(EE)
    SE = math.sin(EE)
    DLE = H*CE - K*SE
    RSAM1 = -K*CE - H*SE
    ASR = 1.0/(1.0 + RSAM1)
    PHI = math.sqrt(1.0 - K*K - H*H)
    PSI = 1.0/(1.0 + PHI)
    X1 = A*(CE - K - PSI*H*DLE)
    Y1 = A*(SE - H + PSI*K*DLE)
    VX1 = AN*ASR*A*(-SE - PSI*H*RSAM1)
    VY1 = AN*ASR*A*(+CE + PSI*K*RSAM1)
    F2 = 2.0*math.sqrt(1.0 - Q*Q - P*P)
    P2 = 1.0 - 2.0*P*P
    Q2 = 1.0 - 2.0*Q*Q
    PQ = 2.0*P*Q
    return StateVector(
        X1*P2 + Y1*PQ,
        X1*PQ + Y1*Q2,
        (Q*Y1 - X1*P)*F2,
        VX1*P2 + VY1*PQ,
        VX1*PQ + VY1*Q2,
        (Q*VY1 - VX1*P)*F2,
        time
    )


def _CalcJupiterMoon(time, mu, al0, al1, a, l, z, zeta):
    # This is a translation of FORTRAN code by Duriez, Lainey, and Vienne:
    # https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

    t = time.tt + 18262.5    # number of days since 1950-01-01T00:00:00Z

    # Calculate 6 orbital elements at the given time t.
    elem0 = 0.0
    for (amplitude, phase, frequency) in a:
        elem0 += amplitude * math.cos(phase + (t * frequency))

    elem1 = al0 + (t * al1)
    for (amplitude, phase, frequency) in l:
        elem1 += amplitude * math.sin(phase + (t * frequency))

    elem1 = math.fmod(elem1, _PI2)
    if elem1 < 0:
        elem1 += _PI2

    elem2 = 0.0
    elem3 = 0.0
    for (amplitude, phase, frequency) in z:
        arg = phase + (t * frequency)
        elem2 += amplitude * math.cos(arg)
        elem3 += amplitude * math.sin(arg)

    elem4 = 0.0
    elem5 = 0.0
    for (amplitude, phase, frequency) in zeta:
        arg = phase + (t * frequency)
        elem4 += amplitude * math.cos(arg)
        elem5 += amplitude * math.sin(arg)

    # Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP).
    state = _JupiterMoon_elem2pv(time, mu, elem0, elem1, elem2, elem3, elem4, elem5)

    # Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ).
    return RotateState(_Rotation_JUP_EQJ, state)


def JupiterMoons(time):
    """Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.

    Calculates position and velocity vectors for Jupiter's moons
    Io, Europa, Ganymede, and Callisto, at the given date and time.
    The vectors are jovicentric (relative to the center of Jupiter).
    Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
    The position components are expressed in astronomical units (AU), and the
    velocity components are in AU/day.

    To convert to heliocentric vectors, call #HelioVector
    with `Body.Jupiter` to get Jupiter's heliocentric position, then
    add the jovicentric vectors. Likewise, you can call #GeoVector
    to convert to geocentric vectors.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate Jupiter's moons.

    Returns
    -------
    JupiterMoonsInfo
        The positions and velocities of Jupiter's 4 largest moons.
    """
    infolist = []
    for (mu, al0, al1, a, l, z, zeta) in _JupiterMoonModel:
        infolist.append(_CalcJupiterMoon(time, mu, al0, al1, a, l, z, zeta))
    return JupiterMoonsInfo(infolist)

# END Jupiter Moons
#----------------------------------------------------------------------------
# BEGIN Search

def _QuadInterp(tm, dt, fa, fm, fb):
    Q = (fb + fa)/2 - fm
    R = (fb - fa)/2
    S = fm

    if Q == 0:
        # This is a line, not a parabola.
        if R == 0:
            # This is a HORIZONTAL line... can't make progress!
            return None
        x = -S / R
        if not (-1 <= x <= +1):
            return None  # out of bounds
    else:
        # It really is a parabola. Find roots x1, x2.
        u = R*R - 4*Q*S
        if u <= 0:
            return None
        ru = math.sqrt(u)
        x1 = (-R + ru) / (2 * Q)
        x2 = (-R - ru) / (2 * Q)

        if -1 <= x1 <= +1:
            if -1 <= x2 <= +1:
                # Two solutions... so parabola intersects twice.
                return None
            x = x1
        elif -1 <= x2 <= +1:
            x = x2
        else:
            return None

    t = tm + x*dt
    df_dt = (2*Q*x + R) / dt
    return (t, df_dt)

def Search(func, context, t1, t2, dt_tolerance_seconds):
    """Searches for a time at which a function's value increases through zero.

    Certain astronomy calculations involve finding a time when an event occurs.
    Often such events can be defined as the root of a function:
    the time at which the function's value becomes zero.

    `Search` finds the *ascending root* of a function: the time at which
    the function's value becomes zero while having a positive slope. That is, as time increases,
    the function transitions from a negative value, through zero at a specific moment,
    to a positive value later. The goal of the search is to find that specific moment.

    The search function is specified by two parameters: `func` and `context`.
    The `func` parameter is a function itself that accepts a time
    and a context containing any other arguments needed to evaluate the function.
    The `context` parameter supplies that context for the given search.
    As an example, a caller may wish to find the moment a celestial body reaches a certain
    ecliptic longitude. In that case, the caller might create a type (class, tuple, whatever)
    that contains a #Body member to specify the body and a numeric value to hold the target longitude.
    A different function might use a completely different context type.

    Every time it is called, `func` returns a `float` value or it raises an exception.
    If `func` raises an exception, the search immediately fails and the exception is
    propagated back to the caller.
    Otherwise, the search proceeds until it either finds the ascending root or fails for some reason.

    The search calls `func` repeatedly to rapidly narrow in on any ascending
    root within the time window specified by `t1` and `t2`. The search never
    reports a solution outside this time window.

    `Search` uses a combination of bisection and quadratic interpolation
    to minimize the number of function calls. However, it is critical that the
    supplied time window be small enough that there cannot be more than one root
    (ascedning or descending) within it; otherwise the search can fail.
    Beyond that, it helps to make the time window as small as possible, ideally
    such that the function itself resembles a smooth parabolic curve within that window.

    If an ascending root is not found, or more than one root
    (ascending and/or descending) exists within the window `t1`..`t2`,
    `Search` will return `None` to indicate a normal search failure.

    If the search does not converge within 20 iterations, it will raise
    an #Error exception.

    Parameters
    ----------
    func : function(context, Time)
        A function that takes an arbitrary context parameter and a #Time parameter.
        Returns a float value.  See remarks above for more details.
    context : object
        An arbitrary data structure needed to be passed to the function `func`
        every time it is called.
    t1 : float
        The lower time bound of the search window.
        See remarks above for more details.
    t2 : float
        The upper time bound of the search window.
        See remarks above for more details.
    dt_tolerance_seconds : float
        Specifies an amount of time in seconds within which a bounded ascending root
        is considered accurate enough to stop. A typical value is 1 second.

    Returns
    -------
    Time or `None`
        If the search is successful, returns a #Time object that is within
        `dt_tolerance_seconds` of an ascending root.
        In this case, the returned time value will always be within the
        inclusive range [`t1`, `t2`].
        If there is no ascending root, or there is more than one ascending root,
        the function returns `None`.

    """
    dt_days = abs(dt_tolerance_seconds / _SECONDS_PER_DAY)
    f1 = func(context, t1)
    f2 = func(context, t2)
    iter_count = 0
    iter_limit = 20
    calc_fmid = True
    while True:
        iter_count += 1
        if iter_count > iter_limit:
            raise Error('Excessive iteration in Search')

        dt = (t2.tt - t1.tt) / 2.0
        tmid = t1.AddDays(dt)
        if abs(dt) < dt_days:
            # We are close enough to the event to stop the search.
            return tmid

        if calc_fmid:
            fmid = func(context, tmid)
        else:
            # We already have the correct value of fmid from the previous loop.
            calc_fmid = True

        # Quadratic interpolation:
        # Try to find a parabola that passes through the 3 points we have sampled:
        # (t1,f1), (tmid,fmid), (t2,f2).
        q = _QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2)
        if q:
            (q_ut, q_df_dt) = q
            tq = Time(q_ut)
            fq = func(context, tq)
            if q_df_dt != 0.0:
                dt_guess = abs(fq / q_df_dt)
                if dt_guess < dt_days:
                    # The estimated time error is small enough that we can quit now.
                    return tq

                # Try guessing a tighter boundary with the interpolated root at the center.
                dt_guess *= 1.2
                if dt_guess < dt / 10.0:
                    tleft = tq.AddDays(-dt_guess)
                    tright = tq.AddDays(+dt_guess)
                    if (tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0.0:
                        if (tright.ut - t1.ut)*(tright.ut - t2.ut) < 0.0:
                            fleft = func(context, tleft)
                            fright = func(context, tright)
                            if fleft < 0.0 and fright >= 0.0:
                                f1 = fleft
                                f2 = fright
                                t1 = tleft
                                t2 = tright
                                fmid = fq
                                calc_fmid = False
                                continue

        # Quadratic interpolation attempt did not work out.
        # Just divide the region in two parts and pick whichever one appears to contain a root.
        if f1 < 0.0 and fmid >= 0.0:
            t2 = tmid
            f2 = fmid
            continue

        if fmid < 0.0 and f2 >= 0.0:
            t1 = tmid
            f1 = fmid
            continue

        # Either there is no ascending zero-crossing in this range
        # or the search window is too wide (more than one zero-crossing).
        return None

# END Search
#----------------------------------------------------------------------------


def HelioVector(body, time):
    """Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.

    This function calculates the position of the given celestial body as a vector,
    using the center of the Sun as the origin.  The result is expressed as a Cartesian
    vector in the J2000 equatorial system: the coordinates are based on the mean equator
    of the Earth at noon UTC on 1 January 2000.

    The position is not corrected for light travel time or aberration.
    This is different from the behavior of #GeoVector.

    If given an invalid value for `body`, this function raises an exception.

    Parameters
    ----------
    body : Body
        The celestial body whose heliocentric position is to be calculated:
        The Sun, Moon, EMB, SSB, or any of the planets.
    time : Time
        The time at which to calculate the heliocentric position.

    Returns
    -------
    Vector
        A heliocentric position vector of the center of the given body
        at the given time.
    """
    if body == Body.Pluto:
        planet = _CalcPluto(time, True)
        return Vector(planet.x, planet.y, planet.z, time)

    if 0 <= body.value < len(_vsop):
        return _CalcVsop(_vsop[body.value], time)

    if body == Body.Sun:
        return Vector(0.0, 0.0, 0.0, time)

    if body == Body.Moon:
        e = _CalcEarth(time)
        m = GeoMoon(time)
        return Vector(e.x+m.x, e.y+m.y, e.z+m.z, time)

    if body == Body.EMB:
        e = _CalcEarth(time)
        m = GeoMoon(time)
        d = 1.0 + _EARTH_MOON_MASS_RATIO
        return Vector(e.x+(m.x/d), e.y+(m.y/d), e.z+(m.z/d), time)

    if body == Body.SSB:
        return _CalcSolarSystemBarycenter(time)

    raise InvalidBodyError()


def HelioDistance(body, time):
    """Calculates the distance between a body and the Sun at a given time.

    Given a date and time, this function calculates the distance between
    the center of `body` and the center of the Sun.
    For the planets Mercury through Neptune, this function is significantly
    more efficient than calling #HelioVector followed by taking the length
    of the resulting vector.

    Parameters
    ----------
    body : Body
        A body for which to calculate a heliocentric distance:
        the Sun, Moon, or any of the planets.
    time : Time
        The date and time for which to calculate the heliocentric distance.

    Returns
    -------
    float
        The heliocentric distance in AU.
    """
    if body == Body.Sun:
        return 0.0

    if 0 <= body.value < len(_vsop):
        return _VsopHelioDistance(_vsop[body.value], time)

    return HelioVector(body, time).Length()


class PositionFunction(abc.ABC):
    """A function for which to solve a light-travel time problem.

    This abstract class defines the contract for wrapping a
    position vector as a function of time. A class derived from
    `PositionFunction` must define a `Position` method that
    returns a position vector for a given time.

    The function #CorrectLightTravel solves a generalized
    problem of deducing how far in the past light must have
    left a target object to be seen by an observer at a
    specified time. It is passed an instance of `PositionFunction`
    that expresses a relative position vector function.
    """
    def __init__(self):
        pass

    @abc.abstractmethod
    def Position(self, time):
        """Returns a relative position vector for a given time.

        Parameters
        ----------
        time : Time
            The time at which to evaluate a relative position vector.

        Returns
        -------
        Vector
        """

def CorrectLightTravel(func, time):
    """Solve for light travel time of a vector function.

    When observing a distant object, for example Jupiter as seen from Earth,
    the amount of time it takes for light to travel from the object to the
    observer can significantly affect the object's apparent position.
    This function is a generic solver that figures out how long in the
    past light must have left the observed object to reach the observer
    at the specified observation time. It uses #PositionFunction
    to express an arbitrary position vector as a function of time.

    This function repeatedly calls `func.Position`, passing a series of time
    estimates in the past. Then `func.Position` must return a relative state vector between
    the observer and the target. `CorrectLightTravel` keeps calling
    `func.Position` with more and more refined estimates of the time light must have
    left the target to arrive at the observer.

    For common use cases, it is simpler to use #BackdatePosition
    for calculating the light travel time correction of one body observing another body.

    For geocentric calculations, #GeoVector also backdates the returned
    position vector for light travel time, only it returns the observation time in
    the returned vector's `t` field rather than the backdated time.

    Parameters
    ----------
    func : PositionFunction
         An arbitrary position vector as a function of time.

    time : Time
        The observation time for which to solve for light travel delay.

    Returns
    -------
    Vector
        The position vector at the solved backdated time.
        The `t` field holds the time that light left the observed
        body to arrive at the observer at the observation time.
    """
    ltime = time
    for _ in range(10):
        pos = func.Position(ltime)
        ltime2 = time.AddDays(-pos.Length() / C_AUDAY)
        dt = abs(ltime2.tt - ltime.tt)
        if dt < 1.0e-9:     # 86.4 microseconds
            return pos
        ltime = ltime2
    raise NoConvergeError()


class _BodyPosition(PositionFunction):
    def __init__(self, observerBody, targetBody, aberration, observerPos):
        super().__init__()
        self.observerBody = observerBody
        self.targetBody = targetBody
        self.aberration = aberration
        self.observerPos = observerPos

    def Position(self, time):
        if self.aberration:
            # The following discussion is worded with the observer body being the Earth,
            # which is often the case. However, the same reasoning applies to any observer body
            # without loss of generality.
            #
            # To include aberration, make a good first-order approximation
            # by backdating the Earth's position also.
            # This is confusing, but it works for objects within the Solar System
            # because the distance the Earth moves in that small amount of light
            # travel time (a few minutes to a few hours) is well approximated
            # by a line segment that substends the angle seen from the remote
            # body viewing Earth. That angle is pretty close to the aberration
            # angle of the moving Earth viewing the remote body.
            # In other words, both of the following approximate the aberration angle:
            #     (transverse distance Earth moves) / (distance to body)
            #     (transverse speed of Earth) / (speed of light).
            observerPos = HelioVector(self.observerBody, time)
        else:
            # No aberration, so use the pre-calculated initial position of
            # the observer body that is already stored in this object.
            observerPos = self.observerPos
        # Subtract the bodies' heliocentric positions to obtain a relative position vector.
        return HelioVector(self.targetBody, time) - observerPos


def BackdatePosition(time, observerBody, targetBody, aberration):
    """Solve for light travel time correction of apparent position.

    When observing a distant object, for example Jupiter as seen from Earth,
    the amount of time it takes for light to travel from the object to the
    observer can significantly affect the object's apparent position.

    This function solves the light travel time correction for the apparent
    relative position vector of a target body as seen by an observer body
    at a given observation time.

    For geocentric calculations, #GeoVector also includes light
    travel time correction, but the time `t` embedded in its returned vector
    refers to the observation time, not the backdated time that light left
    the observed body. Thus `BackdatePosition` provides direct
    access to the light departure time for callers that need it.

    For a more generalized light travel correction solver, see #CorrectLightTravel.

    Parameters
    ----------
    time : Time
        The time of observation.
    observerBody : Body
        The body to be used as the observation location.
    targetBody : Body
        The body to be observed.
    aberration : bool
        `True` to correct for aberration, or `False` to leave uncorrected.

    Returns
    -------
    Vector
        The position vector at the solved backdated time.
        Its `t` field holds the time that light left the observed
        body to arrive at the observer at the observation time.
    """
    if aberration:
        # With aberration, `BackdatePosition` will calculate `observerPos` at different times.
        # Therefore, do not waste time calculating it now.
        # Provide a placeholder value.
        observerPos = None
    else:
        # Without aberration, we need the observer body position at the observation time only.
        # For efficiency, calculate it once and hold onto it, so `BodyPosition` can keep using it.
        observerPos = HelioVector(observerBody, time)
    func = _BodyPosition(observerBody, targetBody, aberration, observerPos)
    return CorrectLightTravel(func, time)


def GeoVector(body, time, aberration):
    """Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.

    This function calculates the position of the given celestial body as a vector,
    using the center of the Earth as the origin.  The result is expressed as a Cartesian
    vector in the J2000 equatorial system: the coordinates are based on the mean equator
    of the Earth at noon UTC on 1 January 2000.

    If given an invalid value for `body`, this function will raise an exception.

    Unlike #HelioVector, this function corrects for light travel time.
    This means the position of the body is "back-dated" by the amount of time it takes
    light to travel from that body to an observer on the Earth.

    Also, the position can optionally be corrected for
    [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
    causing the apparent direction of the body to be shifted due to transverse
    movement of the Earth with respect to the rays of light coming from that body.

    Parameters
    ----------
    body : Body
        A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.
    time : Time
        The date and time for which to calculate the position.
    aberration : bool
        A boolean value indicating whether to correct for aberration.

    Returns
    -------
    Vector
        A geocentric position vector of the center of the given body.
    """
    if body == Body.Moon:
        return GeoMoon(time)

    if body == Body.Earth:
        return Vector(0.0, 0.0, 0.0, time)

    # Correct for light-travel time, to get position of body as seen from Earth's center.
    vec = BackdatePosition(time, Body.Earth, body, aberration)

    # Tricky: return the observation time, not the backdated time.
    vec.t = time
    return vec


def _ExportState(terse, time):
    return StateVector(
        terse.r.x, terse.r.y, terse.r.z,
        terse.v.x, terse.v.y, terse.v.z,
        time
    )


def BaryState(body, time):
    """Calculates barycentric position and velocity vectors for the given body.

    Given a body and a time, calculates the barycentric position and velocity
    vectors for the center of that body at that time.
    The vectors are expressed in equatorial J2000 coordinates (EQJ).

    Parameters
    ----------
    body : Body
        The celestial body whose barycentric state vector is to be calculated.
        Supported values are `Body.Sun`, `Body.SSB`, `Body.Moon`, `Body.EMB`, and all planets:
        `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
        `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
    time : Time
        The date and time for which to calculate position and velocity.

    Returns
    -------
    StateVector
        An object that contains barycentric position and velocity vectors.
    """
    # Trivial case: the solar sytem barycenter itself.
    if body == Body.SSB:
        return StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time)

    if body == Body.Pluto:
        return _CalcPluto(time, False)

    # Find the barycentric positions and velocities for the 5 major bodies.
    bary = _major_bodies_t(time.tt)

    # If the caller is asking for one of the major bodies,
    # we can immediately return the answer.
    if body == Body.Sun:
        return _ExportState(bary.Sun, time)

    if body == Body.Jupiter:
        return _ExportState(bary.Jupiter, time)

    if body == Body.Saturn:
        return _ExportState(bary.Saturn, time)

    if body == Body.Uranus:
        return _ExportState(bary.Uranus, time)

    if body == Body.Neptune:
        return _ExportState(bary.Neptune, time)

    if body in [Body.Moon, Body.EMB]:
        earth = _CalcVsopPosVel(_vsop[Body.Earth.value], time.tt)
        state = GeoMoonState(time) if body == Body.Moon else GeoEmbState(time)
        return StateVector(
            state.x  + bary.Sun.r.x + earth.r.x,
            state.y  + bary.Sun.r.y + earth.r.y,
            state.z  + bary.Sun.r.z + earth.r.z,
            state.vx + bary.Sun.v.x + earth.v.x,
            state.vy + bary.Sun.v.y + earth.v.y,
            state.vz + bary.Sun.v.z + earth.v.z,
            time
        )

    if 0 <= body.value < len(_vsop):
        # Handle the remaining VSOP bodies: Mercury, Venus, Earth, Mars.
        planet = _CalcVsopPosVel(_vsop[body.value], time.tt)
        return StateVector(
            bary.Sun.r.x + planet.r.x,
            bary.Sun.r.y + planet.r.y,
            bary.Sun.r.z + planet.r.z,
            bary.Sun.v.x + planet.v.x,
            bary.Sun.v.y + planet.v.y,
            bary.Sun.v.z + planet.v.z,
            time
        )

    raise InvalidBodyError()


def HelioState(body, time):
    """Calculates heliocentric position and velocity vectors for the given body.

    Given a body and a time, calculates the position and velocity
    vectors for the center of that body at that time, relative to the center of the Sun.
    The vectors are expressed in equatorial J2000 coordinates (EQJ).
    If you need the position vector only, it is more efficient to call #HelioVector.
    The Sun's center is a non-inertial frame of reference. In other words, the Sun
    experiences acceleration due to gravitational forces, mostly from the larger
    planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
    kinetic energy, or other quantities that require a non-accelerating frame
    of reference, consider using #BaryState instead.

    Parameters
    ----------
    body : Body
        The celestial body whose heliocentric state vector is to be calculated.
        Supported values are `Body.Sun`, `Body.SSB`, `Body.Moon`, `Body.EMB`, and all planets:
        `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
        `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
    time : Time
        The date and time for which to calculate position and velocity.

    Returns
    -------
    StateVector
        An object that contains heliocentric position and velocity vectors.
    """
    if body == Body.Sun:
        # Trivial case: the Sun is the origin of the heliocentric frame.
        return StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time)

    if body == Body.SSB:
        # Calculate the barycentric Sun. Then the negative of that is the heliocentric SSB.
        bary = _major_bodies_t(time.tt)
        return StateVector(
            -bary.Sun.r.x,
            -bary.Sun.r.y,
            -bary.Sun.r.z,
            -bary.Sun.v.x,
            -bary.Sun.v.y,
            -bary.Sun.v.z,
            time
        )

    if 0 <= body.value < len(_vsop):
        # Planets included in the VSOP87 model.
        planet = _CalcVsopPosVel(_vsop[body.value], time.tt)
        return _ExportState(planet, time)

    if body == Body.Pluto:
        return _CalcPluto(time, True)

    if body in [Body.Moon, Body.EMB]:
        earth = _CalcVsopPosVel(_vsop[Body.Earth.value], time.tt)
        state = GeoMoonState(time) if body == Body.Moon else GeoEmbState(time)
        return StateVector(
            state.x  + earth.r.x,
            state.y  + earth.r.y,
            state.z  + earth.r.z,
            state.vx + earth.v.x,
            state.vy + earth.v.y,
            state.vz + earth.v.z,
            time
        )

    raise InvalidBodyError()


def Equator(body, time, observer, ofdate, aberration):
    """Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.

    Calculates topocentric equatorial coordinates in one of two different systems:
    J2000 or true-equator-of-date, depending on the value of the `ofdate` parameter.
    Equatorial coordinates include right ascension, declination, and distance in astronomical units.

    This function corrects for light travel time: it adjusts the apparent location
    of the observed body based on how long it takes for light to travel from the body to the Earth.

    This function corrects for *topocentric parallax*, meaning that it adjusts for the
    angular shift depending on where the observer is located on the Earth. This is most
    significant for the Moon, because it is so close to the Earth. However, parallax corection
    has a small effect on the apparent positions of other bodies.

    Correction for aberration is optional, using the `aberration` parameter.

    Parameters
    ----------
    body : Body
        The celestial body to be observed. Not allowed to be `Body.Earth`.
    time : Time
        The date and time at which the observation takes place.
    observer : Observer
        A location on or near the surface of the Earth.
    ofdate : bool
        Selects the date of the Earth's equator in which to express the equatorial coordinates.
        If `True`, returns coordinates using the equator and equinox of date.
        If `False`, returns coordinates converted to the J2000 system.
    aberration : bool
        If `True`, corrects for aberration of light based on the motion of the Earth
        with respect to the heliocentric origin.
        If `False`, does not correct for aberration.

    Returns
    -------
    Equatorial
        Equatorial coordinates in the specified frame of reference.
    """
    gc_observer = _geo_pos(time, observer)
    gc = GeoVector(body, time, aberration)
    j2000 = [
        gc.x - gc_observer[0],
        gc.y - gc_observer[1],
        gc.z - gc_observer[2]
    ]
    if not ofdate:
        return _vector2radec(j2000, time)
    temp = _precession(j2000, time, _PrecessDir.From2000)
    datevect = _nutation(temp, time, _PrecessDir.From2000)
    return _vector2radec(datevect, time)


def ObserverVector(time, observer, ofdate):
    """Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.

    This function calculates a vector from the center of the Earth to
    a point on or near the surface of the Earth, expressed in equatorial
    coordinates. It takes into account the rotation of the Earth at the given
    time, along with the given latitude, longitude, and elevation of the observer.

    The caller may pass `ofdate` as `True` to return coordinates relative to the Earth's
    equator at the specified time, or `False` to use the J2000 equator.

    The returned vector has components expressed in astronomical units (AU).
    To convert to kilometers, multiply the `x`, `y`, and `z` values by
    the constant value #KM_PER_AU.

    The inverse of this function is also available: #VectorObserver.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the observer's position vector.
    observer : Observer
        The geographic location of a point on or near the surface of the Earth.
    ofdate : bool
        Selects the date of the Earth's equator in which to express the equatorial coordinates.
        The caller may pass `False` to use the orientation of the Earth's equator
        at noon UTC on January 1, 2000, in which case this function corrects for precession
        and nutation of the Earth as it was at the moment specified by the `time` parameter.
        Or the caller may pass `True` to use the Earth's equator at `time`
        as the orientation.

    Returns
    -------
    Vector
        An equatorial vector from the center of the Earth to the specified location
        on (or near) the Earth's surface.
    """
    gast = SiderealTime(time)
    ovec = _terra(observer, gast)
    if not ofdate:
        ovec = _nutation(ovec, time, _PrecessDir.Into2000)
        ovec = _precession(ovec, time, _PrecessDir.Into2000)
    return Vector(ovec[0], ovec[1], ovec[2], time)

def ObserverState(time, observer, ofdate):
    """Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.

    This function calculates position and velocity vectors of an observer
    on or near the surface of the Earth, expressed in equatorial
    coordinates. It takes into account the rotation of the Earth at the given
    time, along with the given latitude, longitude, and elevation of the observer.

    The caller may pass `ofdate` as `True` to return coordinates relative to the Earth's
    equator at the specified time, or `False` to use the J2000 equator.

    The returned position vector has components expressed in astronomical units (AU).
    To convert to kilometers, multiply the `x`, `y`, and `z` values by
    the constant value #KM_PER_AU.
    The returned velocity vector has components expressed in AU/day.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the observer's position and velocity vectors.
    observer : Observer
        The geographic location of a point on or near the surface of the Earth.
    ofdate : bool
        Selects the date of the Earth's equator in which to express the equatorial coordinates.
        The caller may pass `False` to use the orientation of the Earth's equator
        at noon UTC on January 1, 2000, in which case this function corrects for precession
        and nutation of the Earth as it was at the moment specified by the `time` parameter.
        Or the caller may pass `True` to use the Earth's equator at `time`
        as the orientation.

    Returns
    -------
    StateVector
        An equatorial position vector and velocity vector relative to the center of the Earth.
    """
    gast = SiderealTime(time)
    ovec = _terra_posvel(observer, gast)
    state = StateVector(
        ovec[0], ovec[1], ovec[2],
        ovec[3], ovec[4], ovec[5],
        time
    )
    if not ofdate:
        state = _nutation_posvel(state, time, _PrecessDir.Into2000)
        state = _precession_posvel(state, time, _PrecessDir.Into2000)
    return state

def VectorObserver(vector, ofdate):
    """Calculates the geographic location corresponding to an equatorial vector.

    This is the inverse function of #ObserverVector.
    Given a geocentric equatorial vector, it returns the geographic
    latitude, longitude, and elevation for that vector.

    Parameters
    ----------
    vector : Vector
        The geocentric equatorial position vector for which to find geographic coordinates.
        The components are expressed in astronomical units (AU).
        The time `vector.t` determines the Earth's rotation.
    ofdate : bool
        Selects the date of the Earth's equator in which `vector` is expressed.
        The caller may pass `False` to use the orientation of the Earth's equator
        at noon UTC on January 1, 2000, in which case this function corrects for precession
        and nutation of the Earth as it was at the moment specified by the the time `vector.t`.
        Or the caller may pass `True` to use the Earth's equator at `vector.t` as the orientation.

    Returns
    -------
    Observer
        The geographic latitude, longitude, and elevation above sea level
        that corresponds to the given equatorial vector.
    """
    gast = SiderealTime(vector.t)
    ovec = [vector.x, vector.y, vector.z]
    if not ofdate:
        ovec = _precession(ovec, vector.t, _PrecessDir.From2000)
        ovec = _nutation(ovec, vector.t, _PrecessDir.From2000)
    return _inverse_terra(ovec, gast)

def ObserverGravity(latitude, height):
    """Calculates the gravitational acceleration experienced by an observer on the Earth.

    This function implements the WGS 84 Ellipsoidal Gravity Formula.
    The result is a combination of inward gravitational acceleration
    with outward centrifugal acceleration, as experienced by an observer
    in the Earth's rotating frame of reference.
    The resulting value increases toward the Earth's poles and decreases
    toward the equator, consistent with changes of the weight measured
    by a spring scale of a fixed mass moved to different latitudes and heights
    on the Earth.

    Parameters
    ----------
    latitude : float
        The latitude of the observer in degrees north or south of the equator.
        By formula symmetry, positive latitudes give the same answer as negative
        latitudes, so the sign does not matter.
    height : float
        The height above the sea level geoid in meters.
        No range checking is done; however, accuracy is only valid in the
        range 0 to 100000 meters.

    Returns
    -------
    float
        The effective gravitational acceleration expressed in meters per second squared [m/s^2].
    """
    s2 = math.sin(math.radians(latitude)) ** 2
    g0 = 9.7803253359 * (1.0 + 0.00193185265241*s2) / math.sqrt(1.0 - 0.00669437999013*s2)
    return g0 * (1.0 - (3.15704e-07 - 2.10269e-09*s2)*height + 7.37452e-14*height*height)

@enum.unique
class Refraction(enum.Enum):
    """Selects if/how to correct for atmospheric refraction.

    Some functions allow enabling or disabling atmospheric refraction
    for the calculated apparent position of a celestial body
    as seen by an observer on the surface of the Earth.

    Values
    ------
    Airless:      No atmospheric refraction correction.
    Normal:       Recommended correction for standard atmospheric refraction.
    JplHorizons:  Used only for compatibility testing with JPL Horizons online tool.
    """
    Airless = 0
    Normal = 1
    JplHorizons = 2

class HorizontalCoordinates:
    """Coordinates of a celestial body as seen by a topocentric observer.

    Contains horizontal and equatorial coordinates as seen by an observer
    on or near the surface of the Earth (a topocentric observer).
    All coordinates are optionally corrected for atmospheric refraction.

    Attributes
    ----------
    azimuth : float
        The compass direction laterally around the observer's horizon,
        measured in degrees.
        North is 0 degrees, east is 90 degrees, south is 180 degrees, etc.
    altitude : float
        The angle in degrees above (positive) or below (negative) the observer's horizon.
    ra : float
        The right ascension in sidereal hours.
    dec : float
        The declination in degrees.
    """
    def __init__(self, azimuth, altitude, ra, dec):
        self.azimuth = azimuth
        self.altitude = altitude
        self.ra = ra
        self.dec = dec

    def __repr__(self):
        return 'HorizontalCoordinates(azimuth={}, altitude={}, ra={}, dec={})'.format(
            self.azimuth,
            self.altitude,
            self.ra,
            self.dec
        )

def Horizon(time, observer, ra, dec, refraction):
    """Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.

    Given a date and time, the geographic location of an observer on the Earth, and
    equatorial coordinates (right ascension and declination) of a celestial body,
    this function returns horizontal coordinates (azimuth and altitude angles) for the body
    relative to the horizon at the geographic location.

    The right ascension `ra` and declination `dec` passed in must be *equator of date*
    coordinates, based on the Earth's true equator at the date and time of the observation.
    Otherwise the resulting horizontal coordinates will be inaccurate.
    Equator of date coordinates can be obtained by calling #Equator, passing in
    `True` as its `ofdate` parameter. It is also recommended to enable
    aberration correction by passing in `True` for the `aberration` parameter.

    This function optionally corrects for atmospheric refraction.
    For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
    correct for optical lensing of the Earth's atmosphere that causes objects
    to appear somewhat higher above the horizon than they actually are.
    However, callers may choose to avoid this correction by passing in `Refraction.Airless`.
    If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
    in the #HorizontalCoordinates object returned by this function will all be corrected for refraction.
    If refraction is disabled, none of these four coordinates will be corrected; in that case,
    the right ascension and declination in the returned object will be numerically identical
    to the respective `ra` and `dec` values passed in.

    Parameters
    ----------
    time : Time
        The date and time for which to find horizontal coordinates.
    observer : Observer
        The location of the observer for which to find horizontal coordinates.
    ra : float
        Right ascension in sidereal hours of the celestial object,
        referred to the mean equinox of date for the J2000 epoch.
    dec : float
        Declination in degrees of the celestial object,
        referred to the mean equator of date for the J2000 epoch.
        Positive values are north of the celestial equator
        and negative values are south of it.
    refraction : Refraction
        The option for selecting whether to correct for atmospheric lensing.
        If `Refraction.Normal`, a well-behaved refraction model is used.
        If `Refraction.Airless`, no refraction correct is performed.
        `Refraction.JplHorizons` is used only for compatibility testing
        with the JPL Horizons online tool.

    Returns
    -------
    HorizontalCoordinates
        The horizontal coordinates (altitude and azimuth), along with
        equatorial coordinates (right ascension and declination), all
        optionally corrected for atmospheric refraction. See remarks above
        for more details.
    """
    if not (Refraction.Airless.value <= refraction.value <= Refraction.JplHorizons.value):
        raise Error('Invalid refraction type')

    latrad = math.radians(observer.latitude)
    lonrad = math.radians(observer.longitude)
    decrad = math.radians(dec)
    rarad = ra * _HOUR2RAD

    sinlat = math.sin(latrad)
    coslat = math.cos(latrad)
    sinlon = math.sin(lonrad)
    coslon = math.cos(lonrad)
    sindc = math.sin(decrad)
    cosdc = math.cos(decrad)
    sinra = math.sin(rarad)
    cosra = math.cos(rarad)

    # Calculate three mutually perpendicular unit vectors
    # in equatorial coordinates: uze, une, uwe.
    #
    # uze = The direction of the observer's local zenith (straight up).
    # une = The direction toward due north on the observer's horizon.
    # uwe = The direction toward due west on the observer's horizon.
    #
    # HOWEVER, these are uncorrected for the Earth's rotation due to the time of day.
    #
    # The components of these 3 vectors are as follows:
    # [0] = x = direction from center of Earth toward 0 degrees longitude (the prime meridian) on equator.
    # [1] = y = direction from center of Earth toward 90 degrees west longitude on equator.
    # [2] = z = direction from center of Earth toward the north pole.

    uze = [coslat*coslon, coslat*sinlon, sinlat]
    une = [-sinlat*coslon, -sinlat*sinlon, coslat]
    uwe = [sinlon, -coslon, 0.0]

    # Correct the vectors uze, une, uwe for the Earth's rotation by calculating
    # sideral time. Call spin() for each uncorrected vector to rotate about
    # the Earth's axis to yield corrected unit vectors uz, un, uw.
    # Multiply sidereal hours by -15 to convert to degrees and flip eastward
    # rotation of the Earth to westward apparent movement of objects with time.

    angle = -15.0 * SiderealTime(time)
    uz = _spin(angle, uze)
    un = _spin(angle, une)
    uw = _spin(angle, uwe)

    # Convert angular equatorial coordinates (RA, DEC) to
    # cartesian equatorial coordinates in 'p', using the
    # same orientation system as uze, une, uwe.

    p = [cosdc*cosra, cosdc*sinra, sindc]

    # Use dot products of p with the zenith, north, and west
    # vectors to obtain the cartesian coordinates of the body in
    # the observer's horizontal orientation system.
    #
    # pz = zenith component [-1, +1]
    # pn = north  component [-1, +1]
    # pw = west   component [-1, +1]

    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2]
    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2]
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2]

    # proj is the "shadow" of the body vector along the observer's flat ground.
    proj = math.hypot(pn, pw)

    # Calculate az = azimuth (compass direction clockwise from East.)
    if proj > 0.0:
        # If the body is not exactly straight up/down, it has an azimuth.
        # Invert the angle to produce degrees eastward from north.
        az = math.degrees(-math.atan2(pw, pn))
        if az < 0:
            az += 360
    else:
        # The body is straight up/down, so it does not have an azimuth.
        # Report an arbitrary but reasonable value.
        az = 0.0

    # zd = the angle of the body away from the observer's zenith.
    zd = math.degrees(math.atan2(proj, pz))
    hor_ra = ra
    hor_dec = dec

    if refraction != Refraction.Airless:
        zd0 = zd
        refr = RefractionAngle(refraction, 90.0 - zd)
        zd -= refr
        if refr > 0.0 and zd > 3.0e-4:
            zdrad = math.radians(zd)
            sinzd = math.sin(zdrad)
            coszd = math.cos(zdrad)
            zd0rad = math.radians(zd0)
            sinzd0 = math.sin(zd0rad)
            coszd0 = math.cos(zd0rad)

            pr = [(((p[j] - coszd0 * uz[j]) / sinzd0)*sinzd + uz[j]*coszd) for j in range(3)]
            proj = math.hypot(pr[0], pr[1])
            if proj > 0:
                hor_ra = _RAD2HOUR * math.atan2(pr[1], pr[0])
                if hor_ra < 0:
                    hor_ra += 24
            else:
                hor_ra = 0
            hor_dec = math.degrees(math.atan2(pr[2], proj))

    return HorizontalCoordinates(az, 90.0 - zd, hor_ra, hor_dec)

def RefractionAngle(refraction, altitude):
    """Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.

    Given an altitude angle and a refraction option, calculates
    the amount of "lift" caused by atmospheric refraction.
    This is the number of degrees higher in the sky an object appears
    due to lensing of the Earth's atmosphere.

    Parameters
    ----------
    refraction : Refraction
        The option for selecting whether to correct for atmospheric lensing.
        If `Refraction.Normal`, a well-behaved refraction model is used.
        If `Refraction.Airless`, no refraction correct is performed.
        `Refraction.JplHorizons` is used only for compatibility testing
        with the JPL Horizons online tool.
        Any other value raises an exception.
    altitude : float
        The number of degrees above (positive) or below (negative) the
        horizon an object is, before being corrected for refraction.

    Returns
    -------
    float
        The number of additional degrees of altitude an object appears
        to have, due to atmospheric refraction, depending on the
        option selected by the `refraction` parameter.
    """
    if altitude < -90.0 or altitude > +90.0:
        return 0.0      # No attempt to correct an invalid altitude

    if refraction in (Refraction.Normal, Refraction.JplHorizons):
        # http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        # JPL Horizons says it uses refraction algorithm from
        # Meeus "Astronomical Algorithms", 1991, p. 101-102.
        # I found the following Go implementation:
        # https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        # This is a translation from the function "Saemundsson" there.
        # I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        # This is important because the 'refr' formula below goes crazy near hd = -5.11.
        hd = max(altitude, -1.0)
        refr = (1.02 / math.tan(math.radians((hd+10.3/(hd+5.11))))) / 60.0

        if refraction == Refraction.Normal and altitude < -1.0:
            # In "normal" mode we gradually reduce refraction toward the nadir
            # so that we never get an altitude angle less than -90 degrees.
            # When horizon angle is -1 degrees, the factor is exactly 1.
            # As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
            refr *= (altitude + 90.0) / 89.0
    elif refraction == Refraction.Airless:
        # The caller does not want refraction correction.
        refr = 0.0
    else:
        raise Error('Inalid refraction option')
    return refr

def InverseRefractionAngle(refraction, bent_altitude):
    """Calculates the inverse of an atmospheric refraction angle.

    Given an observed altitude angle that includes atmospheric refraction,
    calculates the negative angular correction to obtain the unrefracted
    altitude. This is useful for cases where observed horizontal
    coordinates are to be converted to another orientation system,
    but refraction first must be removed from the observed position.

    Parameters
    ----------
    refraction : Refraction
        `Refraction.Normal` - corrects for atmospheric refraction (recommended).
        `Refraction.Airless` - no correction is performed.
        `Refraction.JplHorizons` - For JPL Horizons compatibility testing only.
    bent_altitude : float
        The apparent altitude that includes atmospheric refraction.

    Returns
    -------
    float
        The angular adjustment in degrees, to be added to the
        altitude angle to correct for atmospheric lensing.
        This will be less than or equal to zero.
    """
    if bent_altitude < -90.0 or bent_altitude > +90.0:
        return 0.0      # No attempt to correct an invalid altitude
    # Find the pre-adjusted altitude whose refraction correction leads to 'altitude'.
    altitude = bent_altitude - RefractionAngle(refraction, bent_altitude)
    while True:
        # See how close we got. Keep iterating until the solution converges.
        diff = (altitude + RefractionAngle(refraction, altitude)) - bent_altitude
        if abs(diff) < 1.0e-14:
            return altitude - bent_altitude
        altitude -= diff

class EclipticCoordinates:
    """Ecliptic angular and Cartesian coordinates.

    Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
    oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).

    Attributes
    ----------
    vec : Vector
        Ecliptic cartesian vector with the following components:
        x: in the direction of the equinox along the ecliptic plane.
        y: Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox.
        z: Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north.
    elat : float
        Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
    elon : float
        Longitude in degrees around the ecliptic plane prograde from the equinox.
    """
    def __init__(self, vec, elat, elon):
        self.vec = vec
        self.elat = elat
        self.elon = elon

    def __repr__(self):
        return 'EclipticCoordinates({}, elat={}, elon={})'.format(repr(self.vec), self.elat, self.elon)

def _RotateEquatorialToEcliptic(pos, obliq_radians, time):
    cos_ob = math.cos(obliq_radians)
    sin_ob = math.sin(obliq_radians)
    ex = +pos[0]
    ey = +pos[1]*cos_ob + pos[2]*sin_ob
    ez = -pos[1]*sin_ob + pos[2]*cos_ob
    xyproj = math.hypot(ex, ey)
    if xyproj > 0.0:
        elon = math.degrees(math.atan2(ey, ex))
        if elon < 0.0:
            elon += 360.0
    else:
        elon = 0.0
    elat = math.degrees(math.atan2(ez, xyproj))
    vec = Vector(ex, ey, ez, time)
    return EclipticCoordinates(vec, elat, elon)

def SunPosition(time):
    """Calculates geocentric ecliptic coordinates for the Sun.

    This function calculates the position of the Sun as seen from the Earth.
    The returned value includes both Cartesian and spherical coordinates.
    The x-coordinate and longitude values in the returned object are based
    on the *true equinox of date*: one of two points in the sky where the instantaneous
    plane of the Earth's equator at the given date and time (the *equatorial plane*)
    intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
    By convention, the apparent location of the Sun at the March equinox is chosen
    as the longitude origin and x-axis direction, instead of the one for September.

    `SunPosition` corrects for precession and nutation of the Earth's axis
    in order to obtain the exact equatorial plane at the given time.

    This function can be used for calculating changes of seasons: equinoxes and solstices.
    In fact, the function #Seasons does use this function for that purpose.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Sun's position.

    Returns
    -------
    EclipticCoordinates
        The ecliptic coordinates of the Sun using the Earth's true equator of date.
    """
    # Correct for light travel time from the Sun.
    # Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes!
    adjusted_time = time.AddDays(-1.0 / C_AUDAY)
    earth2000 = _CalcEarth(adjusted_time)
    sun2000 = [-earth2000.x, -earth2000.y, -earth2000.z]

    # Convert to equatorial Cartesian coordinates of date.
    stemp = _precession(sun2000, adjusted_time, _PrecessDir.From2000)
    sun_ofdate = _nutation(stemp, adjusted_time, _PrecessDir.From2000)

    # Convert equatorial coordinates to ecliptic coordinates.
    true_obliq = math.radians(adjusted_time._etilt().tobl)
    return _RotateEquatorialToEcliptic(sun_ofdate, true_obliq, time)

def Ecliptic(equ):
    """Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.

    Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
    on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
    which are relative to the plane of the Earth's orbit around the Sun.

    Parameters
    ----------
    equ : Equatorial
        Equatorial coordinates in the J2000 frame of reference.

    Returns
    -------
    EclipticCoordinates
        Ecliptic coordinates in the J2000 frame of reference.
    """
    # Based on NOVAS functions equ2ecl() and equ2ecl_vec().
    ob2000 = 0.40909260059599012   # mean obliquity of the J2000 ecliptic in radians
    return _RotateEquatorialToEcliptic([equ.x, equ.y, equ.z], ob2000, equ.t)

def EclipticLongitude(body, time):
    """Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.

    This function calculates the angle around the plane of the Earth's orbit
    of a celestial body, as seen from the center of the Sun.
    The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
    in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).

    Parameters
    ----------
    body : Body
        A body other than the Sun.
    time : Time
        The date and time at which the body's ecliptic longitude is to be calculated.

    Returns
    -------
    float
        An angular value in degrees indicating the ecliptic longitude of the body.
    """
    if body == Body.Sun:
        raise InvalidBodyError()
    hv = HelioVector(body, time)
    eclip = Ecliptic(hv)
    return eclip.elon

def AngleFromSun(body, time):
    """Returns the angle between the given body and the Sun, as seen from the Earth.

    This function calculates the angular separation between the given body and the Sun,
    as seen from the center of the Earth. This angle is helpful for determining how
    easy it is to see the body away from the glare of the Sun.

    Parameters
    ----------
    body : Body
        The celestial body whose angle from the Sun is to be measured.
        Not allowed to be `Body.Earth`.
    time : Time
        The time at which the observation is made.

    Returns
    -------
    float
        A numeric value indicating the angle in degrees between the Sun
        and the specified body as seen from the center of the Earth.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    sv = GeoVector(Body.Sun, time, True)
    bv = GeoVector(body, time, True)
    return AngleBetween(sv, bv)

def PairLongitude(body1, body2, time):
    """Returns one body's ecliptic longitude with respect to another, as seen from the Earth.

    This function determines where one body appears around the ecliptic plane
    (the plane of the Earth's orbit around the Sun) as seen from the Earth,
    relative to the another body's apparent position.
    The function returns an angle in the half-open range [0, 360) degrees.
    The value is the ecliptic longitude of `body1` relative to the ecliptic
    longitude of `body2`.

    The angle is 0 when the two bodies are at the same ecliptic longitude
    as seen from the Earth. The angle increases in the prograde direction
    (the direction that the planets orbit the Sun and the Moon orbits the Earth).

    When the angle is 180 degrees, it means the two bodies appear on opposite sides
    of the sky for an Earthly observer.

    Neither `body1` nor `body2` is allowed to be `Body.Earth`.
    If this happens, the function throws an exception.

    Parameters
    ----------
    body1 : Body
        The first body, whose longitude is to be found relative to the second body.
    body2 : Body
        The second body, relative to which the longitude of the first body is to be found.
    time : Time
        The date and time of the observation.

    Returns
    -------
    float
        An angle in degrees in the range [0, 360).
    """
    if body1 == Body.Earth or body2 == Body.Earth:
        raise EarthNotAllowedError()
    vector1 = GeoVector(body1, time, False)
    eclip1 = Ecliptic(vector1)
    vector2 = GeoVector(body2, time, False)
    eclip2 = Ecliptic(vector2)
    return _NormalizeLongitude(eclip1.elon - eclip2.elon)

class ElongationEvent:
    """Contains information about the visibility of a celestial body at a given date and time.

    See the #Elongation function for more detailed information about the members of this class.
    See also #SearchMaxElongation for how to search for maximum elongation events.

    Attributes
    ----------
    time : Time
        The date and time of the observation.
    visibility : Visibility
        Whether the body is best seen in the morning or the evening.
    elongation : float
        The angle in degrees between the body and the Sun, as seen from the Earth.
    ecliptic_separation : float
        The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.
    """
    def __init__(self, time, visibility, elongation, ecliptic_separation):
        self.time = time
        self.visibility = visibility
        self.elongation = elongation
        self.ecliptic_separation = ecliptic_separation

    def __repr__(self):
        return 'ElongationEvent({}, {}, elongation={}, ecliptic_separation={})'.format(
            repr(self.time),
            self.visibility,
            self.elongation,
            self.ecliptic_separation
        )

@enum.unique
class Visibility(enum.Enum):
    """Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.

    Values
    ------
    Morning : The body is best visible in the morning, before sunrise.
    Evening : The body is best visible in the evening, after sunset.
    """
    Morning = 0
    Evening = 1

def Elongation(body, time):
    """Determines visibility of a celestial body relative to the Sun, as seen from the Earth.

    This function returns an #ElongationEvent object, which provides the following
    information about the given celestial body at the given time:

    - `visibility` is an enumerated type that specifies whether the body is more
      easily seen in the morning before sunrise, or in the evening after sunset.

    - `elongation` is the angle in degrees between two vectors: one from the center
      of the Earth to the center of the Sun, the other from the center of the Earth
      to the center of the specified body. This angle indicates how far away the body
      is from the glare of the Sun. The elongation angle is always in the range [0, 180].

    - `ecliptic_separation` is the absolute value of the difference between the body's
      ecliptic longitude and the Sun's ecliptic longitude, both as seen from the center
      of the Earth. This angle measures around the plane of the Earth's orbit, and ignores
      how far above or below that plane the body is.
      The ecliptic separation is measured in degrees and is always in the range [0, 180].

    Parameters
    ----------
    body : Body
        The celestial body whose visibility is to be calculated.
    time : Time
        The date and time of the observation.

    Returns
    -------
    ElongationEvent
    """
    angle = PairLongitude(body, Body.Sun, time)
    if angle > 180.0:
        visibility = Visibility.Morning
        esep = 360.0 - angle
    else:
        visibility = Visibility.Evening
        esep = angle
    angle = AngleFromSun(body, time)
    return ElongationEvent(time, visibility, angle, esep)

def _rlon_offset(body, time, direction, targetRelLon):
    plon = EclipticLongitude(body, time)
    elon = EclipticLongitude(Body.Earth, time)
    diff = direction * (elon - plon)
    return _LongitudeOffset(diff - targetRelLon)

def SearchRelativeLongitude(body, targetRelLon, startTime):
    """Searches for when the Earth and another planet are separated by a certain ecliptic longitude.

    Searches for the time when the Earth and another planet are separated by a specified angle
    in ecliptic longitude, as seen from the Sun.

    A relative longitude is the angle between two bodies measured in the plane of the
    Earth's orbit (the ecliptic plane). The distance of the bodies above or below the ecliptic
    plane is ignored. If you imagine the shadow of the body cast onto the ecliptic plane,
    and the angle measured around that plane from one body to the other in the direction
    the planets orbit the Sun, you will get an angle somewhere between 0 and 360 degrees.
    This is the relative longitude.

    Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
    this function searches for the next time that the relative longitude measured from the
    planet to the Earth is `targetRelLon`.

    Certain astronomical events are defined in terms of relative longitude between
    the Earth and another planet:

    - When the relative longitude is 0 degrees, it means both planets are in the same
      direction from the Sun. For planets that orbit closer to the Sun (Mercury and Venus),
      this is known as *inferior conjunction*, a time when the other planet becomes very
      difficult to see because of being lost in the Sun's glare.
      (The only exception is in the rare event of a transit, when we see the silhouette
      of the planet passing between the Earth and the Sun.)

    - When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
      this is known as *opposition*. Opposition is when the planet is closest to the Earth,
      and also when it is visible for most of the night, so it is considered the best time
      to observe the planet.

    - When the relative longitude is 180 degrees, it means the other planet is on the opposite
      side of the Sun from the Earth.  This is called *superior conjunction*.  Like inferior
      conjunction, the planet is very difficult to see from the Earth.
      Superior conjunction is possible for any planet other than the Earth.

    Parameters
    ----------
    body : Body
        A planet other than the Earth. If `body` is not a planet, or if it is `Body.Earth`, an error occurs.
    targetRelLon : float
        The desired relative longitude, expressed in degrees. Must be in the range [0, 360).
    startTime : Time
        The date and time at which to begin the search.

    Returns
    -------
    Time
        The date and time of the relative longitude event.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    if body in (Body.Moon, Body.Sun):
        raise InvalidBodyError()
    syn = _SynodicPeriod(body)
    direction = +1 if _IsSuperiorPlanet(body) else -1
    # Iterate until we converge on the desired event.
    # Calculate the error angle, which will be a negative number of degrees,
    # meaning we are "behind" the target relative longitude.
    error_angle = _rlon_offset(body, startTime, direction, targetRelLon)
    if error_angle > 0.0:
        error_angle -= 360.0    # force searching forward in time
    time = startTime
    iter_count = 0
    while iter_count < 100:
        # Estimate how many days in the future (positive) or past (negative)
        # we have to go to get closer to the target relative longitude.
        day_adjust = (-error_angle/360.0) * syn
        time = time.AddDays(day_adjust)
        if abs(day_adjust) * _SECONDS_PER_DAY < 1.0:
            return time
        prev_angle = error_angle
        error_angle = _rlon_offset(body, time, direction, targetRelLon)
        if abs(prev_angle) < 30.0 and prev_angle != error_angle:
            # Improve convergence for Mercury/Mars (eccentric orbits)
            # by adjusting the synodic period to more closely match the
            # variable speed of both planets in this part of their respective orbits.
            ratio = prev_angle / (prev_angle - error_angle)
            if 0.5 < ratio < 2.0:
                syn *= ratio
        iter_count += 1
    raise NoConvergeError()

def _neg_elong_slope(body, time):
    dt = 0.1
    t1 = time.AddDays(-dt/2.0)
    t2 = time.AddDays(+dt/2.0)
    e1 = AngleFromSun(body, t1)
    e2 = AngleFromSun(body, t2)
    return (e1 - e2)/dt

def SearchMaxElongation(body, startTime):
    """Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.

    Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
    Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
    The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
    a telescope without atmospheric interference, are when these planets reach maximum elongation.
    These are events where the planets reach the maximum angle from the Sun as seen from the Earth.

    This function solves for those times, reporting the next maximum elongation event's date and time,
    the elongation value itself, the relative longitude with the Sun, and whether the planet is best
    observed in the morning or evening. See #ElongationEvent for more details about the returned object.

    Parameters
    ----------
    body : Body
        Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception.
        To find the best viewing opportunities for planets farther from the Sun than the
        Earth is (Mars through Pluto), use #SearchRelativeLongitude to find the next opposition event.
    startTime : Time
        The date and time at which to begin the search. The maximum elongation event
        found will always be the first one that occurs after this date and time.

    Returns
    -------
    ElongationEvent
    """
    if body == Body.Mercury:
        s1 = 50.0
        s2 = 85.0
    elif body == Body.Venus:
        s1 = 40.0
        s2 = 50.0
    else:
        raise InvalidBodyError()
    syn = _SynodicPeriod(body)
    iter_count = 1
    while iter_count <= 2:
        plon = EclipticLongitude(body, startTime)
        elon = EclipticLongitude(Body.Earth, startTime)
        rlon = _LongitudeOffset(plon - elon)    # clamp to (-180, +180]

        # The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        # because there is a cusp there that causes a discontinuity in the derivative.
        # So we need to guard against searching near such times.
        if -s1 <= rlon < +s1:
            # Seek to the window [+s1, +s2].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = +s1.
            rlon_lo = +s1
            # Search forward for the time t2 when rel lon = +s2.
            rlon_hi = +s2
        elif rlon > +s2 or rlon < -s2:
            # Seek to the next search window at [-s2, -s1].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = -s2.
            rlon_lo = -s2
            # Search forward for the time t2 when rel lon = -s1.
            rlon_hi = -s1
        elif rlon >= 0.0:
            # rlon must be in the middle of the window [+s1, +s2].
            # Search BACKWARD for the time t1 when rel lon = +s1.
            adjust_days = -syn / 4.0
            rlon_lo = +s1
            rlon_hi = +s2
            # Search forward from t1 to find t2 such that rel lon = +s2.
        else:
            # rlon must be in the middle of the window [-s2, -s1].
            # Search BACKWARD for the time t1 when rel lon = -s2.
            adjust_days = -syn / 4.0
            rlon_lo = -s2
            # Search forward from t1 to find t2 such that rel lon = -s1.
            rlon_hi = -s1

        t_start = startTime.AddDays(adjust_days)
        t1 = SearchRelativeLongitude(body, rlon_lo, t_start)
        if t1 is None:
            return None

        t2 = SearchRelativeLongitude(body, rlon_hi, t1)
        if t2 is None:
            return None

        # Now we have a time range [t1,t2] that brackets a maximum elongation event.
        # Confirm the bracketing.
        m1 = _neg_elong_slope(body, t1)
        if m1 >= 0.0:
            raise InternalError()   # there is a bug in the bracketing algorithm!

        m2 = _neg_elong_slope(body, t2)
        if m2 <= 0.0:
            raise InternalError()   # there is a bug in the bracketing algorithm!

        # Use the generic search algorithm to home in on where the slope crosses from negative to positive.
        tx = Search(_neg_elong_slope, body, t1, t2, 10.0)
        if tx is None:
            return None

        if tx.tt >= startTime.tt:
            return Elongation(body, tx)

        # This event is in the past (earlier than startTime).
        # We need to search forward from t2 to find the next possible window.
        # We never need to search more than twice.
        startTime = t2.AddDays(1.0)
        iter_count += 1


def _sun_offset(targetLon, time):
    ecl = SunPosition(time)
    return _LongitudeOffset(ecl.elon - targetLon)

def SearchSunLongitude(targetLon, startTime, limitDays):
    """Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.

    This function finds the moment in time, if any exists in the given time window,
    that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.

    This function can be used to determine equinoxes and solstices.
    However, it is usually more convenient and efficient to call #Seasons
    to calculate all equinoxes and solstices for a given calendar year.

    The function searches the window of time specified by `startTime` and `startTime+limitDays`.
    The search will return `None` if the Sun never reaches the longitude `targetLon` or
    if the window is so large that the longitude ranges more than 180 degrees within it.
    It is recommended to keep the window smaller than 10 days when possible.

    Parameters
    ----------
    targetLon : float
         The desired ecliptic longitude in degrees, relative to the true equinox of date.
         This may be any value in the range [0, 360), although certain values have
         conventional meanings:
         0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice.
    startTime : Time
         The date and time for starting the search for the desired longitude event.
    limitDays : float
         The real-valued number of days, which when added to `startTime`, limits the
         range of time over which the search looks.
         It is recommended to keep this value between 1 and 10 days.
         See remarks above for more details.

    Returns
    -------
    Time or `None`
    """
    t2 = startTime.AddDays(limitDays)
    return Search(_sun_offset, targetLon, startTime, t2, 0.01)

def MoonPhase(time):
    """Returns the Moon's phase as an angle from 0 to 360 degrees.

    This function determines the phase of the Moon using its apparent
    ecliptic longitude relative to the Sun, as seen from the center of the Earth.
    Certain values of the angle have conventional definitions:

    - 0 = new moon
    - 90 = first quarter
    - 180 = full moon
    - 270 = third quarter

    Parameters
    ----------
    time : Time
         The date and time of the observation.

    Returns
    -------
    float
    """
    return PairLongitude(Body.Moon, Body.Sun, time)

def _moon_offset(targetLon, time):
    angle = MoonPhase(time)
    return _LongitudeOffset(angle - targetLon)

def SearchMoonPhase(targetLon, startTime, limitDays):
    """Searches for the time that the Moon reaches a specified phase.

    Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
    longitude with respect to the Sun's geocentric ecliptic longitude.
    When the Moon and the Sun have the same longitude, that is defined as a new moon.
    When their longitudes are 180 degrees apart, that is defined as a full moon.

    This function searches for any value of the lunar phase expressed as an
    angle in degrees in the range [0, 360).

    If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
    it is much easier to call the functions #SearchMoonQuarter and #NextMoonQuarter.
    This function is useful for finding general phase angles outside those four quarters.

    Parameters
    ----------
    targetLon : float
         The difference in geocentric longitude between the Sun and Moon
         that specifies the lunar phase being sought. This can be any value
         in the range [0, 360).  Certain values have conventional names:
         0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter.
    startTime : Time
         The beginning of the time window in which to search for the Moon reaching the specified phase.
    limitDays : float
         The number of days away from `startTime` that limits the time window for the search.
         If the value is negative, the search is performed into the past from `startTime`.
         Otherwise, the search is performed into the future from `startTime`.

    Returns
    -------
    Time or `None`
    """
    # To avoid discontinuities in the _moon_offset function causing problems,
    # we need to approximate when that function will next return 0.
    # We probe it with the start time and take advantage of the fact
    # that every lunar phase repeats roughly every 29.5 days.
    # There is a surprising uncertainty in the quarter timing,
    # due to the eccentricity of the moon's orbit.
    # I have seen more than 0.9 days away from the simple prediction.
    # To be safe, we take the predicted time of the event and search
    # +/-1.5 days around it (a 3-day wide window).
    # But we must return None if the final result goes beyond limitDays after startTime.
    uncertainty = 1.5
    ya = _moon_offset(targetLon, startTime)
    if limitDays < 0.0:
        # Search backward in time.
        if ya < 0.0:
            ya += 360.0
        est_dt = -(_MEAN_SYNODIC_MONTH * ya) / 360.0
        dt2 = est_dt + uncertainty
        if dt2 < limitDays:
            return None     # not possible for moon phase to occur within the specified window
        dt1 = max(limitDays, est_dt - uncertainty)
    else:
        # Search forward in time.
        if ya > 0.0:
            ya -= 360.0
        est_dt = -(_MEAN_SYNODIC_MONTH * ya) / 360.0
        dt1 = est_dt - uncertainty
        if dt1 > limitDays:
            return None     # not possible for moon phase to occur within the specified window
        dt2 = min(limitDays, est_dt + uncertainty)
    t1 = startTime.AddDays(dt1)
    t2 = startTime.AddDays(dt2)
    return Search(_moon_offset, targetLon, t1, t2, 0.1)

class MoonQuarter:
    """A lunar quarter event along with its date and time.

    An object of this type represents one of the four major
    lunar phases that appear on calendars:
    new moon, first quarter, full moon, or third quarter.
    Along with the `quarter` attribute that specifies the
    type of quarter, it contains a `time` field that indicates
    when the lunar quarter event happens.

    Attributes
    ----------
    quarter : int
        0=new moon, 1=first quarter, 2=full moon, 3=third quarter.
    time : Time
        The date and time of the lunar quarter.
    """
    def __init__(self, quarter, time):
        self.quarter = quarter
        self.time = time

    def __repr__(self):
        return 'MoonQuarter({}, {})'.format(self.quarter, repr(self.time))

def SearchMoonQuarter(startTime):
    """Finds the first lunar quarter after the specified date and time.

    A lunar quarter is one of the following four lunar phase events:
    new moon, first quarter, full moon, third quarter.
    This function finds the lunar quarter that happens soonest
    after the specified date and time.

    To continue iterating through consecutive lunar quarters, call this function once,
    followed by calls to #NextMoonQuarter as many times as desired.

    Parameters
    ----------
    startTime : Time
        The date and time at which to start the search.

    Returns
    -------
    MoonQuarter
    """
    angle = MoonPhase(startTime)
    quarter = int(1 + math.floor(angle / 90.0)) % 4
    time = SearchMoonPhase(90.0 * quarter, startTime, 10.0)
    if time is None:
        # The search should never fail. We should always find another lunar quarter.
        raise InternalError()
    return MoonQuarter(quarter, time)

def NextMoonQuarter(mq):
    """Continues searching for lunar quarters from a previous search.

    After calling #SearchMoonQuarter, this function can be called
    one or more times to continue finding consecutive lunar quarters.
    This function finds the next consecutive moon quarter event after
    the one passed in as the parameter `mq`.

    Parameters
    ----------
    mq : MoonQuarter
        A value returned by a prior call to #SearchMoonQuarter or #NextMoonQuarter.

    Returns
    -------
    MoonQuarter
    """
    # Skip 6 days past the previous found moon quarter to find the next one.
    # This is less than the minimum possible increment.
    # So far I have seen the interval well contained by the range (6.5, 8.3) days.
    time = mq.time.AddDays(6.0)
    next_mq = SearchMoonQuarter(time)
    # Verify that we found the expected moon quarter.
    if next_mq.quarter != (1 + mq.quarter) % 4:
        raise InternalError()
    return next_mq


class IlluminationInfo:
    """Information about the brightness and illuminated shape of a celestial body.

    Returned by functions #Illumination and #SearchPeakMagnitude
    to report the visual magnitude and illuminated fraction of a celestial
    body at a given date and time.

    Attributes
    ----------
    time : Time
        The date and time of the observation.
    mag : float
        The visual magnitude of the body. Smaller values are brighter.
    phase_angle : float
        The angle in degrees between the Sun and the Earth, as seen from the body.
        Indicates the body's phase as seen from the Earth.
    phase_fraction : float
        A value in the range [0.0, 1.0] indicating what fraction of the
        body's apparent disc is illuminated, as seen from the Earth.
    helio_dist : float
        The distance between the Sun and the body at the observation time, in AU.
    geo_dist : dist
        The distance between the Earth and the both at the observation time, in AU.
    hc : Vector
        The body's heliocentric vector.
    gc : Vector
        The body's geocentric vector.
    ring_tilt : float
        For Saturn, the tilt angle in degrees of its rings as seen from Earth.
        When the `ring_tilt` is very close to 0, it means the rings are edge-on
        as seen from observers on the Earth, and are thus very difficult to see.
        For bodies other than Saturn, `ring_tilt` is `None`.
    """
    def __init__(self, time, mag, phase, helio_dist, geo_dist, hc, gc, ring_tilt):
        self.time = time
        self.mag = mag
        self.phase_angle = phase
        self.phase_fraction = (1.0 + math.cos(math.radians(phase))) / 2.0
        self.helio_dist = helio_dist
        self.geo_dist = geo_dist
        self.hc = hc
        self.gc = gc
        self.ring_tilt = ring_tilt

    def __repr__(self):
        return 'IlluminationInfo({}, mag={}, phase_angle={}, helio_dist={}, geo_dist={}, hc={}, gc={}, ring_tilt={})'.format(
            repr(self.time),
            self.mag,
            self.phase_angle,
            self.helio_dist,
            self.geo_dist,
            repr(self.hc),
            repr(self.gc),
            repr(self.ring_tilt)
        )

def _MoonMagnitude(phase, helio_dist, geo_dist):
    # https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
    rad = math.radians(phase)
    mag = -12.717 + 1.49*abs(rad) + 0.0431*(rad**4)
    moon_mean_distance_au = 385000.6 / KM_PER_AU
    geo_au = geo_dist / moon_mean_distance_au
    mag += 5.0 * math.log10(helio_dist * geo_au)
    return mag

def _SaturnMagnitude(phase, helio_dist, geo_dist, gc, time):
    # Based on formulas by Paul Schlyter found here:
    # http://www.stjarnhimlen.se/comp/ppcomp.html#15

    # We must handle Saturn's rings as a major component of its visual magnitude.
    # Find geocentric ecliptic coordinates of Saturn.
    eclip = Ecliptic(gc)

    ir = math.radians(28.06)   # tilt of Saturn's rings to the ecliptic, in radians
    Nr = math.radians(169.51 + (3.82e-5 * time.tt))    # ascending node of Saturn's rings, in radians

    # Find tilt of Saturn's rings, as seen from Earth.
    lat = math.radians(eclip.elat)
    lon = math.radians(eclip.elon)
    tilt = math.asin(math.sin(lat)*math.cos(ir) - math.cos(lat)*math.sin(ir)*math.sin(lon-Nr))
    sin_tilt = math.sin(abs(tilt))

    mag = -9.0 + 0.044*phase
    mag += sin_tilt*(-2.6 + 1.2*sin_tilt)
    mag += 5.0 * math.log10(helio_dist * geo_dist)
    ring_tilt = math.degrees(tilt)
    return (mag, ring_tilt)

def _VisualMagnitude(body, phase, helio_dist, geo_dist):
    # For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212
    c0 = c1 = c2 = c3 = 0
    if body == Body.Mercury:
        c0 = -0.60; c1 = +4.98; c2 = -4.88; c3 = +3.02
    elif body == Body.Venus:
        if phase < 163.6:
            c0 = -4.47; c1 = +1.03; c2 = +0.57; c3 = +0.13
        else:
            c0 = +0.98; c1 = -1.02
    elif body == Body.Mars:
        c0 = -1.52; c1 = +1.60
    elif body == Body.Jupiter:
        c0 = -9.40; c1 = +0.50
    elif body == Body.Uranus:
        c0 = -7.19; c1 = +0.25
    elif body == Body.Neptune:
        c0 = -6.87
    elif body == Body.Pluto:
        c0 = -1.00; c1 = +4.00
    else:
        raise InvalidBodyError()

    x = phase / 100.0
    mag = c0 + x*(c1 + x*(c2 + x*c3))
    mag += 5.0 * math.log10(helio_dist * geo_dist)
    return mag

def Illumination(body, time):
    """Finds visual magnitude, phase angle, and other illumination information about a celestial body.

    This function calculates information about how bright a celestial body appears from the Earth,
    reported as visual magnitude, which is a smaller (or even negative) number for brighter objects,
    and a larger number for dimmer objects.

    For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
    the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
    of the body appears illuminated as seen from the Earth. For example, when the phase angle is
    near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
    180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
    of 90 degrees means the body appears "half full".
    For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
    so it doesn't have a phase angle.

    When the body is Saturn, the returned object contains a field `ring_tilt` that holds
    the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
    the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
    0 for all bodies other than Saturn.

    Parameters
    ----------
    body : Body
        The Sun, Moon, or any planet other than the Earth.
    time : Time
        The date and time of the observation.

    Returns
    -------
    IlluminationInfo
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    earth = _CalcEarth(time)
    if body == Body.Sun:
        gc = Vector(-earth.x, -earth.y, -earth.z, time)
        hc = Vector(0.0, 0.0, 0.0, time)
        phase = 0.0     # placeholder value; the Sun does not have a phase angle.
    else:
        if body == Body.Moon:
            # For extra numeric precision, use geocentric moon formula directly.
            gc = GeoMoon(time)
            hc = Vector(earth.x + gc.x, earth.y + gc.y, earth.z + gc.z, time)
        else:
            # For planets, heliocentric vector is most direct to calculate.
            hc = HelioVector(body, time)
            gc = Vector(hc.x - earth.x, hc.y - earth.y, hc.z - earth.z, time)
        phase = AngleBetween(gc, hc)

    geo_dist = gc.Length()      # distance from body to center of Earth
    helio_dist = hc.Length()    # distance from body to center of Sun
    ring_tilt = None            # only reported for Saturn
    if body == Body.Sun:
        mag = -0.17 + 5.0*math.log10(geo_dist / _AU_PER_PARSEC)
    elif body == Body.Moon:
        mag = _MoonMagnitude(phase, helio_dist, geo_dist)
    elif body == Body.Saturn:
        mag, ring_tilt = _SaturnMagnitude(phase, helio_dist, geo_dist, gc, time)
    else:
        mag = _VisualMagnitude(body, phase, helio_dist, geo_dist)
    return IlluminationInfo(time, mag, phase, helio_dist, geo_dist, hc, gc, ring_tilt)

def _mag_slope(body, time):
    # The Search() function finds a transition from negative to positive values.
    # The derivative of magnitude y with respect to time t (dy/dt)
    # is negative as an object gets brighter, because the magnitude numbers
    # get smaller. At peak magnitude dy/dt = 0, then as the object gets dimmer,
    # dy/dt > 0.
    dt = 0.01
    t1 = time.AddDays(-dt/2)
    t2 = time.AddDays(+dt/2)
    y1 = Illumination(body, t1)
    y2 = Illumination(body, t2)
    return (y2.mag - y1.mag) / dt

def SearchPeakMagnitude(body, startTime):
    """Searches for the date and time Venus will next appear brightest as seen from the Earth.

    This function searches for the date and time Venus appears brightest as seen from the Earth.
    Currently only Venus is supported for the `body` parameter, though this could change in the future.
    Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see
    from the Earth, so peak magnitude events have little practical value for that planet.
    Planets other than Venus and Mercury reach peak magnitude at opposition, which can
    be found using #SearchRelativeLongitude.
    The Moon reaches peak magnitude at full moon, which can be found using
    #SearchMoonQuarter or #SearchMoonPhase.
    The Sun reaches peak magnitude at perihelion, which occurs each year in January.
    However, the difference is minor and has little practical value.

    Parameters
    ----------
    body : Body
        Currently only `Body.Venus` is allowed. Any other value results in an exception.
        See remarks above for more details.
    startTime : Time
        The date and time to start searching for the next peak magnitude event.

    Returns
    -------
    IlluminationInfo
    """
    # s1 and s2 are relative longitudes within which peak magnitude of Venus can occur.
    s1 = 10.0
    s2 = 30.0
    if body != Body.Venus:
        raise InvalidBodyError()

    iter_count = 1
    while iter_count <= 2:
        # Find current heliocentric relative longitude between the
        # inferior planet and the Earth.
        plon = EclipticLongitude(body, startTime)
        elon = EclipticLongitude(Body.Earth, startTime)
        rlon = _LongitudeOffset(plon - elon)
        # The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        # because there is a cusp there that causes a discontinuity in the derivative.
        # So we need to guard against searching near such times.
        if -s1 <= rlon < +s1:
            # Seek to the window [+s1, +s2].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = +s1.
            rlon_lo = +s1
            # Search forward for the time t2 when rel lon = +s2.
            rlon_hi = +s2
        elif rlon >= +s2 or rlon < -s2:
            # Seek to the next search window at [-s2, -s1].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = -s2.
            rlon_lo = -s2
            # Search forward for the time t2 when rel lon = -s1.
            rlon_hi = -s1
        elif rlon >= 0:
            # rlon must be in the middle of the window [+s1, +s2].
            # Search BACKWARD for the time t1 when rel lon = +s1.
            syn = _SynodicPeriod(body)
            adjust_days = -syn / 4
            rlon_lo = +s1
            # Search forward from t1 to find t2 such that rel lon = +s2.
            rlon_hi = +s2
        else:
            # rlon must be in the middle of the window [-s2, -s1].
            # Search BACKWARD for the time t1 when rel lon = -s2.
            syn = _SynodicPeriod(body)
            adjust_days = -syn / 4
            rlon_lo = -s2
            # Search forward from t1 to find t2 such that rel lon = -s1.
            rlon_hi = -s1

        t_start = startTime.AddDays(adjust_days)
        t1 = SearchRelativeLongitude(body, rlon_lo, t_start)
        t2 = SearchRelativeLongitude(body, rlon_hi, t1)

        # Now we have a time range [t1,t2] that brackets a maximum magnitude event.
        # Confirm the bracketing.
        m1 = _mag_slope(body, t1)
        if m1 >= 0.0:
            raise InternalError()

        m2 = _mag_slope(body, t2)
        if m2 <= 0.0:
            raise InternalError()

        # Use the generic search algorithm to home in on where the slope crosses from negative to positive.
        tx = Search(_mag_slope, body, t1, t2, 10.0)
        if tx is None:
            # The search should have found the ascending root in the interval [t1, t2].
            raise InternalError()

        if tx.tt >= startTime.tt:
            return Illumination(body, tx)

        # This event is in the past (earlier than startTime).
        # We need to search forward from t2 to find the next possible window.
        # We never need to search more than twice.
        startTime = t2.AddDays(1.0)
        iter_count += 1

    # We should have found the peak magnitude in at most 2 iterations.
    raise InternalError()

class HourAngleEvent:
    """Information about a celestial body crossing a specific hour angle.

    Returned by the function #SearchHourAngle to report information about
    a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

    Attributes
    ----------
    time : Time
        The date and time when the body crosses the specified hour angle.
    hor : HorizontalCoordinates
        Apparent coordinates of the body at the time it crosses the specified hour angle.
    """
    def __init__(self, time, hor):
        self.time = time
        self.hor = hor

    def __repr__(self):
        return 'HourAngleEvent({}, {})'.format(repr(self.time), repr(self.hor))

def SearchHourAngle(body, observer, hourAngle, startTime, direction = +1):
    """Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.

    The *hour angle* of a celestial body indicates its position in the sky with respect
    to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
    The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
    The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
    to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
    the number of hours that have passed since the most recent time that the body has culminated,
    or reached its highest point.

    This function searches for the next time a celestial body reaches the given hour angle
    after the date and time specified by `startTime`.
    To find when a body culminates, pass 0 for `hourAngle`.
    To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.

    Note that, especially close to the Earth's poles, a body as seen on a given day
    may always be above the horizon or always below the horizon, so the caller cannot
    assume that a culminating object is visible nor that an object is below the horizon
    at its minimum altitude.

    On success, the function reports the date and time, along with the horizontal coordinates
    of the body at that time, as seen by the given observer.

    Parameters
    ----------
    body : Body
         The celestial body, which can the Sun, the Moon, or any planet other than the Earth.
    observer : Observer
         Indicates a location on or near the surface of the Earth where the observer is located.
    hourAngle : float
         An hour angle value in the range [0.0, 24.0) indicating the number of sidereal hours after the
         body's most recent culmination.
    startTime : Time
         The date and time at which to start the search.
    direction : int
        The direction in time to perform the search: a positive value
        searches forward in time, a negative value searches backward in time.
        The function throws an exception if `direction` is zero.

    Returns
    -------
    HourAngleEvent
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if hourAngle < 0.0 or hourAngle >= 24.0:
        raise Error('Invalid hour angle.')

    if direction == 0:
        raise Error('Direction must be positive or negative.')

    iter_count = 0
    time = startTime
    while True:
        iter_count += 1
        # Calculate Greenwich Apparent Sidereal Time (GAST) at the given time.
        gast = SiderealTime(time)
        ofdate = Equator(body, time, observer, True, True)

        # Calculate the adjustment needed in sidereal time to bring
        # the hour angle to the desired value.
        delta_sidereal_hours = math.fmod(((hourAngle + ofdate.ra - observer.longitude/15) - gast), 24.0)
        if iter_count == 1:
            # On the first iteration, always search in the requested time direction.
            if direction > 0:
                # Search forward in time.
                if delta_sidereal_hours < 0.0:
                    delta_sidereal_hours += 24.0
            else:
                # Search backward in time.
                if delta_sidereal_hours > 0.0:
                    delta_sidereal_hours -= 24.0
        else:
            # On subsequent iterations, we make the smallest possible adjustment,
            # either forward or backward in time.
            if delta_sidereal_hours < -12.0:
                delta_sidereal_hours += 24.0
            elif delta_sidereal_hours > +12.0:
                delta_sidereal_hours -= 24.0

        # If the error is tolerable (less than 0.1 seconds), stop searching.
        if abs(delta_sidereal_hours) * 3600.0 < 0.1:
            hor = Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.Normal)
            return HourAngleEvent(time, hor)

        # We need to loop another time to get more accuracy.
        # Update the terrestrial time (in solar days) adjusting by sidereal time.
        delta_days = (delta_sidereal_hours / 24.0) * _SOLAR_DAYS_PER_SIDEREAL_DAY
        time = time.AddDays(delta_days)

@enum.unique
class Direction(enum.Enum):
    """Indicates whether a body is rising above or setting below the horizon.

    Specifies the direction of a rising or setting event for a body.
    For example, `Direction.Rise` is used to find sunrise times,
    and `Direction.Set` is used to find sunset times.

    Values
    ------
    Rise:   First appearance of a body as it rises above the horizon.
    Set:    Last appearance of a body as it sinks below the horizon.
    """
    Rise = +1
    Set  = -1

def _ForwardSearchAltitude(body, observer, direction, startTime, limitDays, altitude_error_func, altitude_error_context):
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if direction == Direction.Rise:
        ha_before = 12.0    # minimum altitude (bottom) happens BEFORE the body rises.
        ha_after  =  0.0    # maximum altitude (culmination) happens AFTER the body rises.
    elif direction == Direction.Set:
        ha_before =  0.0    # culmination happens BEFORE the body sets.
        ha_after  = 12.0    # bottom happens AFTER the body sets.
    else:
        raise Error('Invalid value for direction parameter')

    # We cannot possibly satisfy a forward search without a positive time limit.
    if limitDays <= 0.0:
        return None

    # See if the body is currently above/below the horizon.
    # If we are looking for next rise time and the body is below the horizon,
    # we use the current time as the lower time bound and the next culmination
    # as the upper bound.
    # If the body is above the horizon, we search for the next bottom and use it
    # as the lower bound and the next culmination after that bottom as the upper bound.
    # The same logic applies for finding set times, only we swap the hour angles.
    time_start = startTime
    alt_before = altitude_error_func(altitude_error_context, time_start)
    if alt_before > 0.0:
        # We are past the sought event, so we have to wait for the next "before" event (culm/bottom).
        evt_before = SearchHourAngle(body, observer, ha_before, time_start, +1)
        time_before = evt_before.time
        alt_before = altitude_error_func(altitude_error_context, time_before)
    else:
        # We are before or at the sought ebvent, so we find the next "after" event (bottom/culm),
        # and use the current time as the "before" event.
        time_before = time_start

    evt_after = SearchHourAngle(body, observer, ha_after, time_before, +1)
    alt_after = altitude_error_func(altitude_error_context, evt_after.time)

    while True:
        if alt_before <= 0.0 and alt_after > 0.0:
            # Search between the "before time" and the "after time" for the desired event.
            event_time = Search(altitude_error_func, altitude_error_context, time_before, evt_after.time, 1.0)
            if event_time is not None:
                # If we found the rise/set time, but it falls outside limitDays, fail the search.
                if event_time.ut > startTime.ut + limitDays:
                    return None
                # The search succeeded.
                return event_time
        # We didn't find the desired event, so use the "after" time to find the next "before" event.
        evt_before = SearchHourAngle(body, observer, ha_before, evt_after.time, +1)
        if evt_before.time.ut >= time_start.ut + limitDays:
            return None
        evt_after = SearchHourAngle(body, observer, ha_after, evt_before.time, +1)
        time_before = evt_before.time
        alt_before = altitude_error_func(altitude_error_context, evt_before.time)
        alt_after = altitude_error_func(altitude_error_context, evt_after.time)

def _BackwardSearchAltitude(body, observer, direction, startTime, limitDays, altitude_error_func, altitude_error_context):
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if direction == Direction.Rise:
        ha_before = 12.0    # minimum altitude (bottom) happens BEFORE the body rises.
        ha_after  =  0.0    # maximum altitude (culmination) happens AFTER the body rises.
    elif direction == Direction.Set:
        ha_before =  0.0    # culmination happens BEFORE the body sets.
        ha_after  = 12.0    # bottom happens AFTER the body sets.
    else:
        raise Error('Invalid value for direction parameter')

    # We cannot possibly satisfy a backward search without a negative time limit.
    if limitDays >= 0.0:
        return None

    # See if the body is currently above/below the horizon.
    # If we are looking for previous rise time and the body is above the horizon,
    # we use the current time as the upper time bound and the previous bottom as the lower time bound.
    # If the body is below the horizon, we search for the previous culmination and use it
    # as the upper time bound. Then we search for the bottom before that culmination and
    # use it as the lower time bound.
    # The same logic applies for finding set times; altitude_error_func and
    # altitude_error_context ensure that the desired event is represented
    # by ascending through zero, so the Search function works correctly.
    time_start = startTime
    alt_after = altitude_error_func(altitude_error_context, time_start)
    if alt_after < 0.0:
        evt_after = SearchHourAngle(body, observer, ha_after, time_start, -1)
        time_after = evt_after.time
        alt_after = altitude_error_func(altitude_error_context, time_after)
    else:
        time_after = time_start

    evt_before = SearchHourAngle(body, observer, ha_before, time_after, -1)
    alt_before = altitude_error_func(altitude_error_context, evt_before.time)

    while True:
        if alt_before <= 0.0 and alt_after > 0.0:
            # Search between the "before time" and the "after time" for the desired event.
            event_time = Search(altitude_error_func, altitude_error_context, evt_before.time, time_after, 1.0)
            if event_time is not None:
                # If we found the rise/set time, but it falls outside limitDays, fail the search.
                if event_time.ut < startTime.ut + limitDays:
                    return None
                # The search succeeded.
                return event_time
        evt_after = SearchHourAngle(body, observer, ha_after, evt_before.time, -1)
        if evt_after.time.ut <= time_start.ut + limitDays:
            return None
        evt_before = SearchHourAngle(body, observer, ha_before, evt_after.time, -1)
        time_after = evt_after.time
        alt_before = altitude_error_func(altitude_error_context, evt_before.time)
        alt_after = altitude_error_func(altitude_error_context, evt_after.time)


class _peak_altitude_context:
    def __init__(self, body, direction, observer, body_radius_au):
        self.body = body
        self.direction = direction
        self.observer = observer
        self.body_radius_au = body_radius_au


def _peak_altitude(context, time):
    # Return the angular altitude above or below the horizon
    # of the highest part (the peak) of the given object.
    # This is defined as the apparent altitude of the center of the body plus
    # the body's angular radius.
    # The 'direction' parameter controls whether the angle is measured
    # positive above the horizon or positive below the horizon,
    # depending on whether the caller wants rise times or set times, respectively.

    ofdate = Equator(context.body, time, context.observer, True, True)

    # We calculate altitude without refraction, then add fixed refraction near the horizon.
    # This gives us the time of rise/set without the extra work.
    hor = Horizon(time, context.observer, ofdate.ra, ofdate.dec, Refraction.Airless)
    alt = hor.altitude + math.degrees(context.body_radius_au / ofdate.dist)
    return context.direction.value * (alt + _REFRACTION_NEAR_HORIZON)


def SearchRiseSet(body, observer, direction, startTime, limitDays):
    """Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.

    This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
    Rise time is when the body first starts to be visible above the horizon.
    For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
    Set time is the moment when the body appears to vanish below the horizon.

    This function corrects for typical atmospheric refraction, which causes celestial
    bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
    It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).

    Note that rise or set may not occur in every 24 hour period.
    For example, near the Earth's poles, there are long periods of time where
    the Sun stays below the horizon, never rising.
    Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
    This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
    significant amount during each rotation of the Earth.
    Therefore callers must not assume that the function will always succeed.

    Parameters
    ----------
    body : Body
        The Sun, Moon, or any planet other than the Earth.
    observer : Observer
        The location where observation takes place.
    direction : Direction
        Either `Direction.Rise` to find a rise time or `Direction.Set` to find a set time.
    startTime : Time
        The date and time at which to start the search.
    limitDays : float
        Limits how many days to search for a rise or set time, and defines
        the direction in time to search. When `limitDays` is positive, the
        search is performed into the future, after `startTime`.
        When negative, the search is performed into the past, before `startTime`.
        To limit a rise or set time to the same day, you can use a value of 1 day.
        In cases where you want to find the next rise or set time no matter how far
        in the future (for example, for an observer near the south pole), you can
        pass in a larger value like 365.

    Returns
    -------
    Time or `None`
        If the rise or set time is found within the specified time window,
        this function returns that time. Otherwise, it returns `None`.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if body == Body.Sun:
        body_radius = _SUN_RADIUS_AU
    elif body == Body.Moon:
        body_radius = _MOON_EQUATORIAL_RADIUS_AU
    else:
        body_radius = 0.0

    context = _peak_altitude_context(body, direction, observer, body_radius)
    if limitDays < 0.0:
        return _BackwardSearchAltitude(body, observer, direction, startTime, limitDays, _peak_altitude, context)
    return _ForwardSearchAltitude(body, observer, direction, startTime, limitDays, _peak_altitude, context)


class _altitude_error_context:
    def __init__(self, body, direction, observer, altitude):
        self.body = body
        self.direction = direction
        self.observer = observer
        self.altitude = altitude

def _altitude_error_func(context, time):
    ofdate = Equator(context.body, time, context.observer, True, True)
    hor = Horizon(time, context.observer, ofdate.ra, ofdate.dec, Refraction.Airless)
    return context.direction.value * (hor.altitude - context.altitude)


def SearchAltitude(body, observer, direction, dateStart, limitDays, altitude):
    """Finds the next time a body reaches a given altitude.

    Finds when the given body ascends or descends through a given
    altitude angle, as seen by an observer at the specified location on the Earth.
    By using the appropriate combination of `direction` and `altitude` parameters,
    this function can be used to find when civil, nautical, or astronomical twilight
    begins (dawn) or ends (dusk).

    Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
    the horizon. To find civil dawn, pass `Direction.Rise` for `direction` and -6 for `altitude`.

    Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
    To find civil dusk, pass `Direction.Set` for `direction` and -6 for `altitude`.

    Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.

    Astronomical twilight uses -18 degrees as the `altitude` value.

    Parameters
    ----------
    body : Body
        The Sun, Moon, or any planet other than the Earth.
    observer : Observer
        The location where observation takes place.
    direction : Direction
        Either `Direction.Rise` to find an ascending altitude event
        or `Direction.Set` to find a descending altitude event.
    startTime : Time
        The date and time at which to start the search.
    limitDays : float
        Limits how many days to search for the body reaching the altitude angle,
        and defines the direction in time to search. When `limitDays` is positive, the
        search is performed into the future, after `startTime`.
        When negative, the search is performed into the past, before `startTime`.
        To limit the search to the same day, you can use a value of 1 day.
        In cases where you want to find the altitude event no matter how far
        in the future (for example, for an observer near the south pole), you can
        pass in a larger value like 365.
    altitude : float
        The desired altitude angle of the body's center above (positive)
        or below (negative) the observer's local horizon, expressed in degrees.
        Must be in the range [-90, +90].

    Returns
    -------
    Time or `None`
        If the altitude event time is found within the specified time window,
        this function returns that time. Otherwise, it returns `None`.
    """
    if not (-90.0 <= altitude <= +90.0):
        raise Error('Invalid altitude: {}'.format(altitude))
    context = _altitude_error_context(body, direction, observer, altitude)
    if limitDays < 0.0:
        return _BackwardSearchAltitude(body, observer, direction, dateStart, limitDays, _altitude_error_func, context)
    return _ForwardSearchAltitude(body, observer, direction, dateStart, limitDays, _altitude_error_func, context)

class SeasonInfo:
    """The dates and times of changes of season for a given calendar year.

    Call #Seasons to calculate this data structure for a given year.

    Attributes
    ----------
    mar_equinox : Time
        The date and time of the March equinox for the specified year.
    jun_solstice : Time
        The date and time of the June solstice for the specified year.
    sep_equinox : Time
        The date and time of the September equinox for the specified year.
    dec_solstice : Time
        The date and time of the December solstice for the specified year.
    """
    def __init__(self, mar_equinox, jun_solstice, sep_equinox, dec_solstice):
        self.mar_equinox = mar_equinox
        self.jun_solstice = jun_solstice
        self.sep_equinox = sep_equinox
        self.dec_solstice = dec_solstice

    def __repr__(self):
        return 'SeasonInfo(mar_equinox={}, jun_solstice={}, sep_equinox={}, dec_solstice={})'.format(
            repr(self.mar_equinox),
            repr(self.jun_solstice),
            repr(self.sep_equinox),
            repr(self.dec_solstice)
        )

def _FindSeasonChange(targetLon, year, month, day):
    startTime = Time.Make(year, month, day, 0, 0, 0)
    time = SearchSunLongitude(targetLon, startTime, 20.0)
    if time is None:
        # We should always be able to find a season change.
        raise InternalError()
    return time

def Seasons(year):
    """Finds both equinoxes and both solstices for a given calendar year.

    The changes of seasons are defined by solstices and equinoxes.
    Given a calendar year number, this function calculates the
    March and September equinoxes and the June and December solstices.

    The equinoxes are the moments twice each year when the plane of the
    Earth's equator passes through the center of the Sun. In other words,
    the Sun's declination is zero at both equinoxes.
    The March equinox defines the beginning of spring in the northern hemisphere
    and the beginning of autumn in the southern hemisphere.
    The September equinox defines the beginning of autumn in the northern hemisphere
    and the beginning of spring in the southern hemisphere.

    The solstices are the moments twice each year when one of the Earth's poles
    is most tilted toward the Sun. More precisely, the Sun's declination reaches
    its minimum value at the December solstice, which defines the beginning of
    winter in the northern hemisphere and the beginning of summer in the southern
    hemisphere. The Sun's declination reaches its maximum value at the June solstice,
    which defines the beginning of summer in the northern hemisphere and the beginning
    of winter in the southern hemisphere.

    Parameters
    ----------
    year : int
        The calendar year number for which to calculate equinoxes and solstices.
        The value may be any integer, but only the years 1800 through 2100 have
        been validated for accuracy: unit testing against data from the
        United States Naval Observatory confirms that all equinoxes and solstices
        for that range of years are within 2 minutes of the correct time.

    Returns
    -------
    SeasonInfo
    """
    # https://github.com/cosinekitty/astronomy/issues/187
    # Solstices and equinoxes drift over long spans of time,
    # due to precession of the Earth's axis.
    # Therefore, we have to search a wider range of time than
    # one might expect. It turns out this has very little
    # effect on efficiency, thanks to the quick convergence
    # of quadratic interpolation inside the `Search` function.
    mar_equinox  = _FindSeasonChange(  0, year,  3, 10)
    jun_solstice = _FindSeasonChange( 90, year,  6, 10)
    sep_equinox  = _FindSeasonChange(180, year,  9, 10)
    dec_solstice = _FindSeasonChange(270, year, 12, 10)
    return SeasonInfo(mar_equinox, jun_solstice, sep_equinox, dec_solstice)

def _MoonDistance(time):
    return _CalcMoon(time).distance_au

def _moon_distance_slope(direction, time):
    dt = 0.001
    t1 = time.AddDays(-dt/2.0)
    t2 = time.AddDays(+dt/2.0)
    dist1 = _MoonDistance(t1)
    dist2 = _MoonDistance(t2)
    return direction * (dist2 - dist1) / dt

@enum.unique
class ApsisKind(enum.Enum):
    """Represents whether a satellite is at a closest or farthest point in its orbit.

    An apsis is a point in a satellite's orbit that is closest to,
    or farthest from, the body it orbits (its primary).
    `ApsisKind` is an enumerated type that indicates which of these
    two cases applies to a particular apsis event.

    Values
    ------
    Pericenter: The satellite is at its closest point to its primary.
    Apocenter: The satellite is at its farthest point from its primary.
    Invalid: A placeholder for an undefined, unknown, or invalid apsis.
    """
    Pericenter = 0
    Apocenter  = 1
    Invalid    = 2

class Apsis:
    """An event where a satellite is closest to or farthest from the body it orbits.

    For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
    event where the orbiting body reaches its closest or farthest point from the primary body.
    The closest approach is called *pericenter* and the farthest point is *apocenter*.

    More specific terminology is common for particular orbiting bodies.
    The Moon's closest approach to the Earth is called *perigee* and its furthest
    point is called *apogee*. The closest approach of a planet to the Sun is called
    *perihelion* and the furthest point is called *aphelion*.

    This data structure is returned by #SearchLunarApsis and #NextLunarApsis
    to iterate through consecutive alternating perigees and apogees.

    Attributes
    ----------
    time : Time
        The date and time of the apsis.
    kind : ApsisKind
        Whether this is a pericenter or apocenter event.
    dist_au : float
        The distance between the centers of the bodies in astronomical units.
    dist_km : float
        The distance between the centers of the bodies in kilometers.
    """
    def __init__(self, time, kind, dist_au):
        self.time = time
        self.kind = kind
        self.dist_au = dist_au
        self.dist_km = dist_au * KM_PER_AU

    def __repr__(self):
        return 'Apsis({}, {}, dist_au={})'.format(
            repr(self.time),
            self.kind,
            self.dist_au
        )

def SearchLunarApsis(startTime):
    """Finds the time of the first lunar apogee or perigee after the given time.

    Given a date and time to start the search in `startTime`, this function finds
    the next date and time that the center of the Moon reaches the closest or
    farthest point in its orbit with respect to the center of the Earth, whichever
    comes first after `startTime`.  The return value (of type #Apsis) also
    contains an indicator of whether the event is apogee or perigee.

    The closest point is called *perigee* and the farthest point is called *apogee*.
    The word *apsis* refers to either event.

    To iterate through consecutive alternating perigee and apogee events,
    call #SearchLunarApsis once, then use the return value to call #NextLunarApsis.
    After that, keep feeding the previous return value from `NextLunarApsis` into
    another call of `NextLunarApsis` as many times as desired.

    Parameters
    ----------
    startTime : Time
        The date and time at which to start searching for the next perigee or apogee.

    Returns
    -------
    Apsis
    """
    increment = 5.0     # number of days to skip on each iteration
    t1 = startTime
    m1 = _moon_distance_slope(+1, t1)
    iter_count = 0
    while iter_count * increment < 2.0 * _MEAN_SYNODIC_MONTH:
        t2 = t1.AddDays(increment)
        m2 = _moon_distance_slope(+1, t2)
        if m1 * m2 <= 0.0:
            # There is a change of slope polarity within the time range [t1, t2].
            # Therefore this time range contains an apsis.
            # Figure out whether it is perigee or apogee.
            if m1 < 0.0 or m2 > 0.0:
                # We found a minimum-distance event: perigee.
                # Search the time range for the time when the slope goes from negative to positive.
                apsis_time = Search(_moon_distance_slope, +1, t1, t2, 1.0)
                kind = ApsisKind.Pericenter
            elif m1 > 0.0 or m2 < 0.0:
                # We found a maximum-distance event: apogee.
                # Search the time range for the time when the slope goes from positive to negative.
                apsis_time = Search(_moon_distance_slope, -1, t1, t2, 1.0)
                kind = ApsisKind.Apocenter
            else:
                # This should never happen. It should not be possible for both slopes to be zero.
                raise InternalError()

            if apsis_time is None:
                # We should always be able to find a lunar apsis!
                raise InternalError()

            dist = _MoonDistance(apsis_time)
            return Apsis(apsis_time, kind, dist)

        # We have not yet found a slope polarity change. Keep searching.
        t1 = t2
        m1 = m2
        iter_count += 1

    # It should not be possible to fail to find an apsis within 2 synodic months.
    raise InternalError()


def NextLunarApsis(apsis):
    """Finds the next lunar perigee or apogee in a series.

    This function requires an #Apsis value obtained from a call to
    #SearchLunarApsis or `NextLunarApsis`.
    Given an apogee event, this function finds the next perigee event,
    and vice versa.

    See #SearchLunarApsis for more details.

    Parameters
    ----------
    apsis : Apsis

    Returns
    -------
    Apsis
    """
    skip = 11.0     # number of days to skip to start looking for next apsis event
    time = apsis.time.AddDays(skip)
    next_apsis = SearchLunarApsis(time)
    # Verify that we found the opposite apsis from the previous one.
    if apsis.kind not in [ApsisKind.Apocenter, ApsisKind.Pericenter]:
        raise Error('Parameter "apsis" contains an invalid "kind" value.')
    if next_apsis.kind.value + apsis.kind.value != 1:
        raise InternalError()   # should have found opposite apsis kind
    return next_apsis


def _planet_distance_slope(context, time):
    (direction, body) = context
    dt = 0.001
    t1 = time.AddDays(-dt/2.0)
    t2 = time.AddDays(+dt/2.0)
    dist1 = HelioDistance(body, t1)
    dist2 = HelioDistance(body, t2)
    return direction * (dist2 - dist1) / dt


def SearchPlanetApsis(body, startTime):
    """Finds the next planet perihelion or aphelion, after a given time.

    Given a date and time to start the search in `startTime`, this function finds the
    next date and time that the center of the specified planet reaches the closest or farthest point
    in its orbit with respect to the center of the Sun, whichever comes first after `startTime`.

    The closest point is called *perihelion* and the farthest point is called *aphelion*.
    The word *apsis* refers to either event.

    To iterate through consecutive alternating perihelion and aphelion events,
    call `SearchPlanetApsis` once, then use the return value to call #NextPlanetApsis.
    After that, keep feeding the previous return value from `NextPlanetApsis`
    into another call of `NextPlanetApsis` as many times as desired.

    Parameters
    ----------
    body : Body
        The planet for which to find the next perihelion/aphelion event.
        Not allowed to be `Body.Sun` or `Body.Moon`.
    startTime : Time
        The date and time at which to start searching for the next perihelion or aphelion.

    Returns
    -------
    Apsis
    """
    if body == Body.Neptune or body == Body.Pluto:
        return _BruteSearchPlanetApsis(body, startTime)
    positive_slope = (+1.0, body)
    negative_slope = (-1.0, body)
    orbit_period_days = _PlanetOrbitalPeriod[body.value]
    increment = orbit_period_days / 6.0
    t1 = startTime
    m1 = _planet_distance_slope(positive_slope, t1)
    iter_count = 0
    while iter_count * increment < 2 * orbit_period_days:
        t2 = t1.AddDays(increment)
        m2 = _planet_distance_slope(positive_slope, t2)
        if m1 * m2 <= 0.0:
            # There is a change of slope polarity within the time range [t1, t2].
            # Therefore this time range contains an apsis.
            # Figure out whether it is perihelion or aphelion.
            if m1 < 0.0 or m2 > 0.0:
                # We found a minimum-distance event: perihelion.
                # Search the time range for the time when the slope goes from negative to positive.
                context = positive_slope
                kind = ApsisKind.Pericenter
            elif m1 > 0.0 or m2 < 0.0:
                # We found a maximum-distance event: aphelion.
                # Search the time range for the time when the slope goes from positive to negative.
                context = negative_slope
                kind = ApsisKind.Apocenter
            else:
                raise InternalError()   # at least one of the planet distance slopes should have been nonzero
            search = Search(_planet_distance_slope, context, t1, t2, 1.0)
            if search is None:
                raise InternalError()   # failed to find where planet distance slope passed through zero

            dist = HelioDistance(body, search)
            return Apsis(search, kind, dist)
        t1 = t2
        m1 = m2
        iter_count += 1
    raise InternalError()   # should have found planet apsis within 2 planet orbits


def NextPlanetApsis(body, apsis):
    """Finds the next planetary perihelion or aphelion event in a series.

    This function requires an #Apsis value obtained from a call
    to #SearchPlanetApsis or `NextPlanetApsis`.
    Given an aphelion event, this function finds the next perihelion event, and vice versa.
    See #SearchPlanetApsis for more details.

    Parameters
    ----------
    body : Body
        The planet for which to find the next perihelion/aphelion event.
        Not allowed to be `Body.Sun` or `Body.Moon`.
        Must match the body passed into the call that produced the `apsis` parameter.
    apsis : Apsis
        An apsis event obtained from a call to #SearchPlanetApsis or `NextPlanetApsis`.

    Returns
    -------
    Apsis
    """
    if apsis.kind not in [ApsisKind.Apocenter, ApsisKind.Pericenter]:
        raise Error('Parameter "apsis" contains an invalid "kind" value.')
    # Skip 1/4 of an orbit before starting search again.
    skip = 0.25 * _PlanetOrbitalPeriod[body.value]
    time = apsis.time.AddDays(skip)
    next_apsis = SearchPlanetApsis(body, time)
    # Verify that we found the opposite apsis from the previous one.
    if next_apsis.kind.value + apsis.kind.value != 1:
        raise InternalError()   # should have found opposite planetary apsis type
    return next_apsis


def _PlanetExtreme(body, kind, start_time, dayspan):
    direction = +1.0 if (kind == ApsisKind.Apocenter) else -1.0
    npoints = 10
    while True:
        interval = dayspan / (npoints - 1)
        # Iterate until uncertainty is less than one minute.
        if interval < 1/1440:
            apsis_time = start_time.AddDays(interval/2)
            dist_au = HelioDistance(body, apsis_time)
            return Apsis(apsis_time, kind, dist_au)
        for i in range(npoints):
            time = start_time.AddDays(i * interval)
            dist = direction * HelioDistance(body, time)
            if i==0 or dist > best_dist:
                best_i = i
                best_dist = dist
        # Narrow in on the extreme point.
        start_time = start_time.AddDays((best_i - 1) * interval)
        dayspan = 2 * interval


def _BruteSearchPlanetApsis(body, startTime):
    # Neptune is a special case for two reasons:
    # 1. Its orbit is nearly circular (low orbital eccentricity).
    # 2. It is so distant from the Sun that the orbital period is very long.
    # Put together, this causes wobbling of the Sun around the Solar System Barycenter (SSB)
    # to be so significant that there are 3 local minima in the distance-vs-time curve
    # near each apsis. Therefore, unlike for other planets, we can't use an optimized
    # algorithm for finding dr/dt = 0.
    # Instead, we use a dumb, brute-force algorithm of sampling and finding min/max
    # heliocentric distance.
    #
    # There is a similar problem in the TOP2013 model for Pluto:
    # Its position vector has high-frequency oscillations that confuse the
    # slope-based determination of apsides.
    #
    # Rewind approximately 30 degrees in the orbit,
    # then search forward for 270 degrees.
    # This is a very cautious way to prevent missing an apsis.
    # Typically we will find two apsides, and we pick whichever
    # apsis is ealier, but after startTime.
    # Sample points around this orbital arc and find when the distance
    # is greatest and smallest.
    period = _PlanetOrbitalPeriod[body.value]
    t1 = startTime.AddDays(period * ( -30 / 360))
    t2 = startTime.AddDays(period * (+270 / 360))
    t_min = t_max = t1
    npoints = 100
    interval = (t2.ut - t1.ut) / (npoints - 1)

    for i in range(npoints):
        time = t1.AddDays(i * interval)
        dist = HelioDistance(body, time)
        if i == 0:
            min_dist = max_dist = dist
        else:
            if dist > max_dist:
                max_dist = dist
                t_max = time
            if dist < min_dist:
                min_dist = dist
                t_min = time

    perihelion = _PlanetExtreme(body, ApsisKind.Pericenter, t_min.AddDays(-2*interval), 4*interval)
    aphelion = _PlanetExtreme(body, ApsisKind.Apocenter, t_max.AddDays(-2*interval), 4*interval)
    if perihelion.time.tt >= startTime.tt:
        if startTime.tt <= aphelion.time.tt < perihelion.time.tt:
            return aphelion
        return perihelion
    if aphelion.time.tt >= startTime.tt:
        return aphelion
    raise InternalError()   # failed to find Neptune apsis


def VectorFromSphere(sphere, time):
    """Converts spherical coordinates to Cartesian coordinates.

    Given spherical coordinates and a time at which they are valid,
    returns a vector of Cartesian coordinates. The returned value
    includes the time, as required by all `Time` objects.

    Parameters
    ----------
    sphere : Spherical
        Spherical coordinates to be converted.
    time : Time
        The time that should be included in the returned vector.

    Returns
    -------
    Vector
        The vector form of the supplied spherical coordinates.
    """
    radlat = math.radians(sphere.lat)
    radlon = math.radians(sphere.lon)
    rcoslat = sphere.dist * math.cos(radlat)
    return Vector(
        rcoslat * math.cos(radlon),
        rcoslat * math.sin(radlon),
        sphere.dist * math.sin(radlat),
        time
    )


def EquatorFromVector(vec):
    """Given an equatorial vector, calculates equatorial angular coordinates.

    Parameters
    ----------
    vec : Vector
        A vector in an equatorial coordinate system.

    Returns
    -------
    Equatorial
        Angular coordinates expressed in the same equatorial system as `vec`.
    """
    sphere = SphereFromVector(vec)
    return Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist, vec)


def SphereFromVector(vector):
    """Converts Cartesian coordinates to spherical coordinates.

    Given a Cartesian vector, returns latitude, longitude, and distance.

    Parameters
    ----------
    vector : Vector
        Cartesian vector to be converted to spherical coordinates.

    Returns
    -------
    Spherical
        Spherical coordinates that are equivalent to the given vector.
    """
    xyproj = vector.x*vector.x + vector.y*vector.y
    dist = math.sqrt(xyproj + vector.z*vector.z)
    if xyproj == 0.0:
        if vector.z == 0.0:
            raise Exception('Zero-length vector not allowed.')
        lon = 0.0
        if vector.z < 0.0:
            lat = -90.0
        else:
            lat = +90.0
    else:
        lon = math.degrees(math.atan2(vector.y, vector.x))
        if lon < 0.0:
            lon += 360.0
        lat = math.degrees(math.atan2(vector.z, math.sqrt(xyproj)))
    return Spherical(lat, lon, dist)


def _ToggleAzimuthDirection(az):
    az = 360.0 - az
    if az >= 360.0:
        az -= 360.0
    elif az < 0.0:
        az += 360.0
    return az


def VectorFromHorizon(sphere, time, refraction):
    """Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.

    Parameters
    ----------
    sphere : Spherical
        A structure that contains apparent horizontal coordinates:
        `lat` holds the refracted azimuth angle,
        `lon` holds the azimuth in degrees clockwise from north,
        and `dist` holds the distance from the observer to the object in AU.
    time : Time
        The date and time of the observation. This is needed because the returned
        vector object requires a valid time value when passed to certain other functions.
    refraction : Refraction
        See remarks in function #RefractionAngle.

    Returns
    -------
    Vector
        A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
    """
    lon = _ToggleAzimuthDirection(sphere.lon)
    lat = sphere.lat + InverseRefractionAngle(refraction, sphere.lat)
    xsphere = Spherical(lat, lon, sphere.dist)
    return VectorFromSphere(xsphere, time)


def HorizonFromVector(vector, refraction):
    """Converts Cartesian coordinates to horizontal coordinates.

    Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.

    *IMPORTANT:* This function differs from `SphereFromVector` in two ways:
    - `SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
      from north (e.g., west = +90), but this function represents a clockwise rotation
      (e.g., east = +90). The difference is because `SphereFromVector` is intended
      to preserve the vector "right-hand rule", while this function defines azimuth in a more
      traditional way as used in navigation and cartography.
    - This function optionally corrects for atmospheric refraction, while `SphereFromVector` does not.

    The returned object contains the azimuth in `lon`.
    It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.

    The altitude is stored in `lat`.

    The distance to the observed object is stored in `dist`,
    and is expressed in astronomical units (AU).

    Parameters
    ----------
    vector : Vector
        Cartesian vector to be converted to horizontal angular coordinates.
    refraction : Refraction
        See comments in the #RefractionAngle function.
    """
    sphere = SphereFromVector(vector)
    sphere.lon = _ToggleAzimuthDirection(sphere.lon)
    sphere.lat += RefractionAngle(refraction, sphere.lat)
    return sphere


def InverseRotation(rotation):
    """Calculates the inverse of a rotation matrix.

    Given a rotation matrix that performs some coordinate transform,
    this function returns the matrix that reverses that transform.

    Parameters
    ----------
    rotation : RotationMatrix
        The rotation matrix to be inverted.

    Returns
    -------
    RotationMatrix
        The inverse rotation matrix.
    """
    return RotationMatrix([
        [rotation.rot[0][0], rotation.rot[1][0], rotation.rot[2][0]],
        [rotation.rot[0][1], rotation.rot[1][1], rotation.rot[2][1]],
        [rotation.rot[0][2], rotation.rot[1][2], rotation.rot[2][2]]
    ])


def CombineRotation(a, b):
    """Creates a rotation based on applying one rotation followed by another.

    Given two rotation matrices, returns a combined rotation matrix that is
    equivalent to rotating based on the first matrix, followed by the second.

    Parameters
    ----------
    a : RotationMatrix
        The first rotation to apply.
    b : RotationMatrix
        The second rotation to apply.

    Returns
    -------
    RotationMatrix
        The combined rotation matrix.
    """
    # Use matrix multiplication: c = b*a.
    # We put 'b' on the left and 'a' on the right because,
    # just like when you use a matrix M to rotate a vector V,
    # you put the M on the left in the product M*V.
    # We can think of this as 'b' rotating all the 3 column vectors in 'a'.

    return RotationMatrix([
        [
            b.rot[0][0]*a.rot[0][0] + b.rot[1][0]*a.rot[0][1] + b.rot[2][0]*a.rot[0][2],
            b.rot[0][1]*a.rot[0][0] + b.rot[1][1]*a.rot[0][1] + b.rot[2][1]*a.rot[0][2],
            b.rot[0][2]*a.rot[0][0] + b.rot[1][2]*a.rot[0][1] + b.rot[2][2]*a.rot[0][2]
        ],
        [
            b.rot[0][0]*a.rot[1][0] + b.rot[1][0]*a.rot[1][1] + b.rot[2][0]*a.rot[1][2],
            b.rot[0][1]*a.rot[1][0] + b.rot[1][1]*a.rot[1][1] + b.rot[2][1]*a.rot[1][2],
            b.rot[0][2]*a.rot[1][0] + b.rot[1][2]*a.rot[1][1] + b.rot[2][2]*a.rot[1][2]
        ],
        [
            b.rot[0][0]*a.rot[2][0] + b.rot[1][0]*a.rot[2][1] + b.rot[2][0]*a.rot[2][2],
            b.rot[0][1]*a.rot[2][0] + b.rot[1][1]*a.rot[2][1] + b.rot[2][1]*a.rot[2][2],
            b.rot[0][2]*a.rot[2][0] + b.rot[1][2]*a.rot[2][1] + b.rot[2][2]*a.rot[2][2]
        ]
    ])


def IdentityMatrix():
    """Creates an identity rotation matrix.

    Returns a rotation matrix that has no effect on orientation.
    This matrix can be the starting point for other operations,
    such as using a series of calls to #Pivot to
    create a custom rotation matrix.

    Returns
    -------
    RotationMatrix
        The identity rotation matrix.
    """
    return RotationMatrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])


def Pivot(rotation, axis, angle):
    """Re-orients a rotation matrix by pivoting it by an angle around one of its axes.

    Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
    this function pivots the rotation matrix by that angle around that coordinate axis.

    For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
    to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
    of a telescope camera pointed at a given body, you can use `Pivot` twice:
    (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
    western axis by the body's altitude angle. The resulting rotation matrix will then
    reorient ECL coordinates to the orientation of your telescope camera.

    Parameters
    ----------
    rotation : RotationMatrix
        The input rotation matrix.
    axis : int
        An integer that selects which coordinate axis to rotate around:
        0 = x, 1 = y, 2 = z. Any other value will cause an exception.
    angle : float
        An angle in degrees indicating the amount of rotation around the specified axis.
        Positive angles indicate rotation counterclockwise as seen from the positive
        direction along that axis, looking towards the origin point of the orientation system.
        Any finite number of degrees is allowed, but best precision will result from
        keeping `angle` in the range [-360, +360].

    Returns
    -------
    RotationMatrix
        A pivoted matrix object.
    """
    # Check for an invalid coordinate axis.
    if axis not in [0, 1, 2]:
        raise Error('Invalid axis {}. Must be [0, 1, 2].'.format(axis))

    radians = math.radians(angle)
    c = math.cos(radians)
    s = math.sin(radians)

    # We need to maintain the "right-hand" rule, no matter which
    # axis was selected. That means we pick (i, j, k) axis order
    # such that the following vector cross product is satisfied:
    # i x j = k
    i = (axis + 1) % 3
    j = (axis + 2) % 3
    k = axis

    rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    rot[i][i] = c*rotation.rot[i][i] - s*rotation.rot[i][j]
    rot[i][j] = s*rotation.rot[i][i] + c*rotation.rot[i][j]
    rot[i][k] = rotation.rot[i][k]

    rot[j][i] = c*rotation.rot[j][i] - s*rotation.rot[j][j]
    rot[j][j] = s*rotation.rot[j][i] + c*rotation.rot[j][j]
    rot[j][k] = rotation.rot[j][k]

    rot[k][i] = c*rotation.rot[k][i] - s*rotation.rot[k][j]
    rot[k][j] = s*rotation.rot[k][i] + c*rotation.rot[k][j]
    rot[k][k] = rotation.rot[k][k]

    return RotationMatrix(rot)


def RotateVector(rotation, vector):
    """Applies a rotation to a vector, yielding a rotated vector.

    This function transforms a vector in one orientation to a vector
    in another orientation.

    Parameters
    ----------
    rotation : RotationMatrix
        A rotation matrix that specifies how the orientation of the vector is to be changed.
    vector : Vector
        The vector whose orientation is to be changed.

    Returns
    -------
    Vector
        A vector in the orientation specified by `rotation`.
    """
    return Vector(
        rotation.rot[0][0]*vector.x + rotation.rot[1][0]*vector.y + rotation.rot[2][0]*vector.z,
        rotation.rot[0][1]*vector.x + rotation.rot[1][1]*vector.y + rotation.rot[2][1]*vector.z,
        rotation.rot[0][2]*vector.x + rotation.rot[1][2]*vector.y + rotation.rot[2][2]*vector.z,
        vector.t
    )


def RotateState(rotation, state):
    """Applies a rotation to a state vector, yielding a rotated state vector.

    This function transforms a state vector in one orientation to a
    state vector in another orientation. Both the position and velocity
    vectors are rotated the same way.

    Parameters
    ----------
    rotation : RotationMatrix
        A rotation matrix that specifies how the orientation of the vector is to be changed.
    state : StateVector
        The state vector whose orientation is to be changed.

    Returns
    -------
    StateVector
        A state vector in the orientation specified by `rotation`.
    """
    return StateVector(
        rotation.rot[0][0]*state.x + rotation.rot[1][0]*state.y + rotation.rot[2][0]*state.z,
        rotation.rot[0][1]*state.x + rotation.rot[1][1]*state.y + rotation.rot[2][1]*state.z,
        rotation.rot[0][2]*state.x + rotation.rot[1][2]*state.y + rotation.rot[2][2]*state.z,
        rotation.rot[0][0]*state.vx + rotation.rot[1][0]*state.vy + rotation.rot[2][0]*state.vz,
        rotation.rot[0][1]*state.vx + rotation.rot[1][1]*state.vy + rotation.rot[2][1]*state.vz,
        rotation.rot[0][2]*state.vx + rotation.rot[1][2]*state.vy + rotation.rot[2][2]*state.vz,
        state.t
    )


def Rotation_EQJ_ECL():
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using equator at J2000 epoch.
    Target: ECL = ecliptic system, using equator at J2000 epoch.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to ECL.
    """
    # ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
    c = 0.9174821430670688    # cos(ob)
    s = 0.3977769691083922    # sin(ob)
    return RotationMatrix([
        [ 1,  0,  0],
        [ 0, +c, -s],
        [ 0, +s, +c]
    ])


def Rotation_ECL_EQJ():
    """Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: ECL = ecliptic system, using equator at J2000 epoch.
    Target: EQJ = equatorial system, using equator at J2000 epoch.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts ECL to EQJ.
    """
    # ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
    c = 0.9174821430670688    # cos(ob)
    s = 0.3977769691083922    # sin(ob)
    return RotationMatrix([
        [ 1,  0,  0],
        [ 0, +c, +s],
        [ 0, -s, +c]
    ])

def Rotation_EQJ_EQD(time):
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using equator at J2000 epoch.
    Target: EQD = equatorial system, using equator of the specified date/time.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator defines the target orientation.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to EQD at `time`.
    """
    prec = _precession_rot(time, _PrecessDir.From2000)
    nut = _nutation_rot(time, _PrecessDir.From2000)
    return CombineRotation(prec, nut)


def Rotation_EQD_EQJ(time):
    """Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQD = equatorial system, using equator of the specified date/time.
    Target: EQJ = equatorial system, using equator at J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator defines the source orientation.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQD at `time` to EQJ.
    """
    nut = _nutation_rot(time, _PrecessDir.Into2000)
    prec = _precession_rot(time, _PrecessDir.Into2000)
    return CombineRotation(nut, prec)


def Rotation_EQD_HOR(time, observer):
    """Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQD = equatorial system, using equator of the specified date/time.
    Target: HOR = horizontal system.

    Use #HorizonFromVector to convert the return value
    to a traditional altitude/azimuth pair.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator applies.
    observer: Observer
        A location near the Earth's mean sea level that defines the observer's location.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQD to HOR at `time` and for `observer`.
        The components of the horizontal vector are:
        x = north, y = west, z = zenith (straight up from the observer).
        These components are chosen so that the "right-hand rule" works for the vector
        and so that north represents the direction where azimuth = 0.
    """
    sinlat = math.sin(math.radians(observer.latitude))
    coslat = math.cos(math.radians(observer.latitude))
    sinlon = math.sin(math.radians(observer.longitude))
    coslon = math.cos(math.radians(observer.longitude))
    uze = [coslat * coslon, coslat * sinlon, sinlat]
    une = [-sinlat * coslon, -sinlat * sinlon, coslat]
    uwe = [sinlon, -coslon, 0.0]
    spin_angle = -15.0 * SiderealTime(time)
    uz = _spin(spin_angle, uze)
    un = _spin(spin_angle, une)
    uw = _spin(spin_angle, uwe)
    return RotationMatrix([
        [un[0], uw[0], uz[0]],
        [un[1], uw[1], uz[1]],
        [un[2], uw[2], uz[2]],
    ])


def Rotation_HOR_EQD(time, observer):
    """Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: HOR = horizontal system (x=North, y=West, z=Zenith).
    Target: EQD = equatorial system, using equator of the specified date/time.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator applies.
    observer : Observer
        A location near the Earth's mean sea level that defines the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts HOR to EQD at `time` and for `observer`.
    """
    rot = Rotation_EQD_HOR(time, observer)
    return InverseRotation(rot)


def Rotation_HOR_EQJ(time, observer):
    """Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: HOR = horizontal system (x=North, y=West, z=Zenith).
    Target: EQJ = equatorial system, using equator at the J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time of the observation.
    observer : Observer
        A location near the Earth's mean sea level that define's the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts HOR to EQJ at `time` and for `observer`.
    """
    hor_eqd = Rotation_HOR_EQD(time, observer)
    eqd_eqj = Rotation_EQD_EQJ(time)
    return CombineRotation(hor_eqd, eqd_eqj)


def Rotation_EQJ_HOR(time, observer):
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using the equator at the J2000 epoch.
    Target: HOR = horizontal system.

    Use #HorizonFromVector to convert the return value to
    a traditional altitude/azimuth pair.

    Parameters
    ----------
    time : Time
        The date and time of the desired horizontal orientation.
    observer : Observer
        A location near the Earth's mean sea level that defines the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
        The components of the horizontal vector are:
        x = north, y = west, z = zenith (straight up from the observer).
        These components are chosen so that the "right-hand rule" works for the vector
        and so that north represents the direction where azimuth = 0.
    """
    rot = Rotation_HOR_EQJ(time, observer)
    return InverseRotation(rot)


def Rotation_EQD_ECL(time):
    """Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQD = equatorial system, using equator of date.
    Target: ECL = ecliptic system, using equator at J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time of the source equator.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQD to ECL.
    """
    eqd_eqj = Rotation_EQD_EQJ(time)
    eqj_ecl = Rotation_EQJ_ECL()
    return CombineRotation(eqd_eqj, eqj_ecl)


def Rotation_ECL_EQD(time):
    """Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: ECL = ecliptic system, using equator at J2000 epoch.
    Target: EQD = equatorial system, using equator of date.

    Parameters
    ----------
    time : Time
        The date and time of the desired equator.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts ECL to EQD.
    """
    rot = Rotation_EQD_ECL(time)
    return InverseRotation(rot)


def Rotation_ECL_HOR(time, observer):
    """Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: ECL = ecliptic system, using equator at J2000 epoch.
    Target: HOR = horizontal system.

    Use #HorizonFromVector to convert the return value
    to a traditional altitude/azimuth pair.

    Parameters
    ----------
    time : Time
        The date and time of the desired horizontal orientation.
    observer : Observer
        A location near the Earth's mean sea level that defines the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts ECL to HOR at `time` and for `observer`.
        The components of the horizontal vector are:
        x = north, y = west, z = zenith (straight up from the observer).
        These components are chosen so that the "right-hand rule" works for the vector
        and so that north represents the direction where azimuth = 0.
    """
    ecl_eqd = Rotation_ECL_EQD(time)
    eqd_hor = Rotation_EQD_HOR(time, observer)
    return CombineRotation(ecl_eqd, eqd_hor)


def Rotation_HOR_ECL(time, observer):
    """Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: HOR = horizontal system.
    Target: ECL = ecliptic system, using equator at J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time of the horizontal observation.
    observer : Observer
        The location of the horizontal observer.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts HOR to ECL.
    """
    rot = Rotation_ECL_HOR(time, observer)
    return InverseRotation(rot)

def Rotation_EQJ_GAL():
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to galactic (GAL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using the equator at the J2000 epoch.
    Target: GAL = galactic system (IAU 1958 definition).

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to GAL.
    """
    # This rotation matrix was calculated by the following script
    # in this same source code repository:
    # demo/python/galeqj_matrix.py
    return RotationMatrix([
        [-0.0548624779711344, +0.4941095946388765, -0.8676668813529025],
        [-0.8734572784246782, -0.4447938112296831, -0.1980677870294097],
        [-0.4838000529948520, +0.7470034631630423, +0.4559861124470794]
    ])

def Rotation_GAL_EQJ():
    """Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: GAL = galactic system (IAU 1958 definition).
    Target: EQJ = equatorial system, using the equator at the J2000 epoch.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts GAL to EQJ.
    """
    # This rotation matrix was calculated by the following script
    # in this same source code repository:
    # demo/python/galeqj_matrix.py
    return RotationMatrix([
        [-0.0548624779711344, -0.8734572784246782, -0.4838000529948520],
        [+0.4941095946388765, -0.4447938112296831, +0.7470034631630423],
        [-0.8676668813529025, -0.1980677870294097, +0.4559861124470794]
    ])

class ConstellationInfo:
    """Reports the constellation that a given celestial point lies within.

    The #Constellation function returns a `ConstellationInfo` object
    to report which constellation corresponds with a given point in the sky.
    Constellations are defined with respect to the B1875 equatorial system
    per IAU standard. Although the `Constellation` function requires J2000 equatorial
    coordinates as input, the returned object contains converted B1875 coordinates for reference.

    Attributes
    ----------
    symbol : string
        3-character mnemonic symbol for the constellation, e.g. "Ori".
    name : string
        Full name of constellation, e.g. "Orion".
    ra1875 : float
        Right ascension expressed in B1875 coordinates.
    dec1875 : float
        Declination expressed in B1875 coordinates.
    """
    def __init__(self, symbol, name, ra1875, dec1875):
        self.symbol = symbol
        self.name = name
        self.ra1875 = ra1875
        self.dec1875 = dec1875

    def __repr__(self):
        return 'ConstellationInfo(symbol={}, name={}, ra1875={}, dec1875={})'.format(
            repr(self.symbol),
            repr(self.name),
            self.ra1875,
            self.dec1875
        )


_ConstelRot = None
_Epoch2000 = None
_ConstelNames = (
    ('And', 'Andromeda'           )  #  0
,   ('Ant', 'Antila'              )  #  1
,   ('Aps', 'Apus'                )  #  2
,   ('Aql', 'Aquila'              )  #  3
,   ('Aqr', 'Aquarius'            )  #  4
,   ('Ara', 'Ara'                 )  #  5
,   ('Ari', 'Aries'               )  #  6
,   ('Aur', 'Auriga'              )  #  7
,   ('Boo', 'Bootes'              )  #  8
,   ('Cae', 'Caelum'              )  #  9
,   ('Cam', 'Camelopardis'        )  # 10
,   ('Cap', 'Capricornus'         )  # 11
,   ('Car', 'Carina'              )  # 12
,   ('Cas', 'Cassiopeia'          )  # 13
,   ('Cen', 'Centaurus'           )  # 14
,   ('Cep', 'Cepheus'             )  # 15
,   ('Cet', 'Cetus'               )  # 16
,   ('Cha', 'Chamaeleon'          )  # 17
,   ('Cir', 'Circinus'            )  # 18
,   ('CMa', 'Canis Major'         )  # 19
,   ('CMi', 'Canis Minor'         )  # 20
,   ('Cnc', 'Cancer'              )  # 21
,   ('Col', 'Columba'             )  # 22
,   ('Com', 'Coma Berenices'      )  # 23
,   ('CrA', 'Corona Australis'    )  # 24
,   ('CrB', 'Corona Borealis'     )  # 25
,   ('Crt', 'Crater'              )  # 26
,   ('Cru', 'Crux'                )  # 27
,   ('Crv', 'Corvus'              )  # 28
,   ('CVn', 'Canes Venatici'      )  # 29
,   ('Cyg', 'Cygnus'              )  # 30
,   ('Del', 'Delphinus'           )  # 31
,   ('Dor', 'Dorado'              )  # 32
,   ('Dra', 'Draco'               )  # 33
,   ('Equ', 'Equuleus'            )  # 34
,   ('Eri', 'Eridanus'            )  # 35
,   ('For', 'Fornax'              )  # 36
,   ('Gem', 'Gemini'              )  # 37
,   ('Gru', 'Grus'                )  # 38
,   ('Her', 'Hercules'            )  # 39
,   ('Hor', 'Horologium'          )  # 40
,   ('Hya', 'Hydra'               )  # 41
,   ('Hyi', 'Hydrus'              )  # 42
,   ('Ind', 'Indus'               )  # 43
,   ('Lac', 'Lacerta'             )  # 44
,   ('Leo', 'Leo'                 )  # 45
,   ('Lep', 'Lepus'               )  # 46
,   ('Lib', 'Libra'               )  # 47
,   ('LMi', 'Leo Minor'           )  # 48
,   ('Lup', 'Lupus'               )  # 49
,   ('Lyn', 'Lynx'                )  # 50
,   ('Lyr', 'Lyra'                )  # 51
,   ('Men', 'Mensa'               )  # 52
,   ('Mic', 'Microscopium'        )  # 53
,   ('Mon', 'Monoceros'           )  # 54
,   ('Mus', 'Musca'               )  # 55
,   ('Nor', 'Norma'               )  # 56
,   ('Oct', 'Octans'              )  # 57
,   ('Oph', 'Ophiuchus'           )  # 58
,   ('Ori', 'Orion'               )  # 59
,   ('Pav', 'Pavo'                )  # 60
,   ('Peg', 'Pegasus'             )  # 61
,   ('Per', 'Perseus'             )  # 62
,   ('Phe', 'Phoenix'             )  # 63
,   ('Pic', 'Pictor'              )  # 64
,   ('PsA', 'Pisces Austrinus'    )  # 65
,   ('Psc', 'Pisces'              )  # 66
,   ('Pup', 'Puppis'              )  # 67
,   ('Pyx', 'Pyxis'               )  # 68
,   ('Ret', 'Reticulum'           )  # 69
,   ('Scl', 'Sculptor'            )  # 70
,   ('Sco', 'Scorpius'            )  # 71
,   ('Sct', 'Scutum'              )  # 72
,   ('Ser', 'Serpens'             )  # 73
,   ('Sex', 'Sextans'             )  # 74
,   ('Sge', 'Sagitta'             )  # 75
,   ('Sgr', 'Sagittarius'         )  # 76
,   ('Tau', 'Taurus'              )  # 77
,   ('Tel', 'Telescopium'         )  # 78
,   ('TrA', 'Triangulum Australe' )  # 79
,   ('Tri', 'Triangulum'          )  # 80
,   ('Tuc', 'Tucana'              )  # 81
,   ('UMa', 'Ursa Major'          )  # 82
,   ('UMi', 'Ursa Minor'          )  # 83
,   ('Vel', 'Vela'                )  # 84
,   ('Vir', 'Virgo'               )  # 85
,   ('Vol', 'Volans'              )  # 86
,   ('Vul', 'Vulpecula'           )  # 87
)

_ConstelBounds = (
    ( 83,      0,   8640,   2112 )    # UMi
,   ( 83,   2880,   5220,   2076 )    # UMi
,   ( 83,   7560,   8280,   2068 )    # UMi
,   ( 83,   6480,   7560,   2064 )    # UMi
,   ( 15,      0,   2880,   2040 )    # Cep
,   ( 10,   3300,   3840,   1968 )    # Cam
,   ( 15,      0,   1800,   1920 )    # Cep
,   ( 10,   3840,   5220,   1920 )    # Cam
,   ( 83,   6300,   6480,   1920 )    # UMi
,   ( 33,   7260,   7560,   1920 )    # Dra
,   ( 15,      0,   1263,   1848 )    # Cep
,   ( 10,   4140,   4890,   1848 )    # Cam
,   ( 83,   5952,   6300,   1800 )    # UMi
,   ( 15,   7260,   7440,   1800 )    # Cep
,   ( 10,   2868,   3300,   1764 )    # Cam
,   ( 33,   3300,   4080,   1764 )    # Dra
,   ( 83,   4680,   5952,   1680 )    # UMi
,   ( 13,   1116,   1230,   1632 )    # Cas
,   ( 33,   7350,   7440,   1608 )    # Dra
,   ( 33,   4080,   4320,   1596 )    # Dra
,   ( 15,      0,    120,   1584 )    # Cep
,   ( 83,   5040,   5640,   1584 )    # UMi
,   ( 15,   8490,   8640,   1584 )    # Cep
,   ( 33,   4320,   4860,   1536 )    # Dra
,   ( 33,   4860,   5190,   1512 )    # Dra
,   ( 15,   8340,   8490,   1512 )    # Cep
,   ( 10,   2196,   2520,   1488 )    # Cam
,   ( 33,   7200,   7350,   1476 )    # Dra
,   ( 15, 7393.2,   7416,   1462 )    # Cep
,   ( 10,   2520,   2868,   1440 )    # Cam
,   ( 82,   2868,   3030,   1440 )    # UMa
,   ( 33,   7116,   7200,   1428 )    # Dra
,   ( 15,   7200, 7393.2,   1428 )    # Cep
,   ( 15,   8232,   8340,   1418 )    # Cep
,   ( 13,      0,    876,   1404 )    # Cas
,   ( 33,   6990,   7116,   1392 )    # Dra
,   ( 13,    612,    687,   1380 )    # Cas
,   ( 13,    876,   1116,   1368 )    # Cas
,   ( 10,   1116,   1140,   1368 )    # Cam
,   ( 15,   8034,   8232,   1350 )    # Cep
,   ( 10,   1800,   2196,   1344 )    # Cam
,   ( 82,   5052,   5190,   1332 )    # UMa
,   ( 33,   5190,   6990,   1332 )    # Dra
,   ( 10,   1140,   1200,   1320 )    # Cam
,   ( 15,   7968,   8034,   1320 )    # Cep
,   ( 15,   7416,   7908,   1316 )    # Cep
,   ( 13,      0,    612,   1296 )    # Cas
,   ( 50,   2196,   2340,   1296 )    # Lyn
,   ( 82,   4350,   4860,   1272 )    # UMa
,   ( 33,   5490,   5670,   1272 )    # Dra
,   ( 15,   7908,   7968,   1266 )    # Cep
,   ( 10,   1200,   1800,   1260 )    # Cam
,   ( 13,   8232,   8400,   1260 )    # Cas
,   ( 33,   5670,   6120,   1236 )    # Dra
,   ( 62,    735,    906,   1212 )    # Per
,   ( 33,   6120,   6564,   1212 )    # Dra
,   ( 13,      0,    492,   1200 )    # Cas
,   ( 62,    492,    600,   1200 )    # Per
,   ( 50,   2340,   2448,   1200 )    # Lyn
,   ( 13,   8400,   8640,   1200 )    # Cas
,   ( 82,   4860,   5052,   1164 )    # UMa
,   ( 13,      0,    402,   1152 )    # Cas
,   ( 13,   8490,   8640,   1152 )    # Cas
,   ( 39,   6543,   6564,   1140 )    # Her
,   ( 33,   6564,   6870,   1140 )    # Dra
,   ( 30,   6870,   6900,   1140 )    # Cyg
,   ( 62,    600,    735,   1128 )    # Per
,   ( 82,   3030,   3300,   1128 )    # UMa
,   ( 13,     60,    312,   1104 )    # Cas
,   ( 82,   4320,   4350,   1080 )    # UMa
,   ( 50,   2448,   2652,   1068 )    # Lyn
,   ( 30,   7887,   7908,   1056 )    # Cyg
,   ( 30,   7875,   7887,   1050 )    # Cyg
,   ( 30,   6900,   6984,   1044 )    # Cyg
,   ( 82,   3300,   3660,   1008 )    # UMa
,   ( 82,   3660,   3882,    960 )    # UMa
,   (  8,   5556,   5670,    960 )    # Boo
,   ( 39,   5670,   5880,    960 )    # Her
,   ( 50,   3330,   3450,    954 )    # Lyn
,   (  0,      0,    906,    882 )    # And
,   ( 62,    906,    924,    882 )    # Per
,   ( 51,   6969,   6984,    876 )    # Lyr
,   ( 62,   1620,   1689,    864 )    # Per
,   ( 30,   7824,   7875,    864 )    # Cyg
,   ( 44,   7875,   7920,    864 )    # Lac
,   (  7,   2352,   2652,    852 )    # Aur
,   ( 50,   2652,   2790,    852 )    # Lyn
,   (  0,      0,    720,    840 )    # And
,   ( 44,   7920,   8214,    840 )    # Lac
,   ( 44,   8214,   8232,    828 )    # Lac
,   (  0,   8232,   8460,    828 )    # And
,   ( 62,    924,    978,    816 )    # Per
,   ( 82,   3882,   3960,    816 )    # UMa
,   ( 29,   4320,   4440,    816 )    # CVn
,   ( 50,   2790,   3330,    804 )    # Lyn
,   ( 48,   3330,   3558,    804 )    # LMi
,   (  0,    258,    507,    792 )    # And
,   (  8,   5466,   5556,    792 )    # Boo
,   (  0,   8460,   8550,    770 )    # And
,   ( 29,   4440,   4770,    768 )    # CVn
,   (  0,   8550,   8640,    752 )    # And
,   ( 29,   5025,   5052,    738 )    # CVn
,   ( 80,    870,    978,    736 )    # Tri
,   ( 62,    978,   1620,    736 )    # Per
,   (  7,   1620,   1710,    720 )    # Aur
,   ( 51,   6543,   6969,    720 )    # Lyr
,   ( 82,   3960,   4320,    696 )    # UMa
,   ( 30,   7080,   7530,    696 )    # Cyg
,   (  7,   1710,   2118,    684 )    # Aur
,   ( 48,   3558,   3780,    684 )    # LMi
,   ( 29,   4770,   5025,    684 )    # CVn
,   (  0,      0,     24,    672 )    # And
,   ( 80,    507,    600,    672 )    # Tri
,   (  7,   2118,   2352,    672 )    # Aur
,   ( 37,   2838,   2880,    672 )    # Gem
,   ( 30,   7530,   7824,    672 )    # Cyg
,   ( 30,   6933,   7080,    660 )    # Cyg
,   ( 80,    690,    870,    654 )    # Tri
,   ( 25,   5820,   5880,    648 )    # CrB
,   (  8,   5430,   5466,    624 )    # Boo
,   ( 25,   5466,   5820,    624 )    # CrB
,   ( 51,   6612,   6792,    624 )    # Lyr
,   ( 48,   3870,   3960,    612 )    # LMi
,   ( 51,   6792,   6933,    612 )    # Lyr
,   ( 80,    600,    690,    600 )    # Tri
,   ( 66,    258,    306,    570 )    # Psc
,   ( 48,   3780,   3870,    564 )    # LMi
,   ( 87,   7650,   7710,    564 )    # Vul
,   ( 77,   2052,   2118,    548 )    # Tau
,   (  0,     24,     51,    528 )    # And
,   ( 73,   5730,   5772,    528 )    # Ser
,   ( 37,   2118,   2238,    516 )    # Gem
,   ( 87,   7140,   7290,    510 )    # Vul
,   ( 87,   6792,   6930,    506 )    # Vul
,   (  0,     51,    306,    504 )    # And
,   ( 87,   7290,   7404,    492 )    # Vul
,   ( 37,   2811,   2838,    480 )    # Gem
,   ( 87,   7404,   7650,    468 )    # Vul
,   ( 87,   6930,   7140,    460 )    # Vul
,   (  6,   1182,   1212,    456 )    # Ari
,   ( 75,   6792,   6840,    444 )    # Sge
,   ( 59,   2052,   2076,    432 )    # Ori
,   ( 37,   2238,   2271,    420 )    # Gem
,   ( 75,   6840,   7140,    388 )    # Sge
,   ( 77,   1788,   1920,    384 )    # Tau
,   ( 39,   5730,   5790,    384 )    # Her
,   ( 75,   7140,   7290,    378 )    # Sge
,   ( 77,   1662,   1788,    372 )    # Tau
,   ( 77,   1920,   2016,    372 )    # Tau
,   ( 23,   4620,   4860,    360 )    # Com
,   ( 39,   6210,   6570,    344 )    # Her
,   ( 23,   4272,   4620,    336 )    # Com
,   ( 37,   2700,   2811,    324 )    # Gem
,   ( 39,   6030,   6210,    308 )    # Her
,   ( 61,      0,     51,    300 )    # Peg
,   ( 77,   2016,   2076,    300 )    # Tau
,   ( 37,   2520,   2700,    300 )    # Gem
,   ( 61,   7602,   7680,    300 )    # Peg
,   ( 37,   2271,   2496,    288 )    # Gem
,   ( 39,   6570,   6792,    288 )    # Her
,   ( 31,   7515,   7578,    284 )    # Del
,   ( 61,   7578,   7602,    284 )    # Peg
,   ( 45,   4146,   4272,    264 )    # Leo
,   ( 59,   2247,   2271,    240 )    # Ori
,   ( 37,   2496,   2520,    240 )    # Gem
,   ( 21,   2811,   2853,    240 )    # Cnc
,   ( 61,   8580,   8640,    240 )    # Peg
,   (  6,    600,   1182,    238 )    # Ari
,   ( 31,   7251,   7308,    204 )    # Del
,   (  8,   4860,   5430,    192 )    # Boo
,   ( 61,   8190,   8580,    180 )    # Peg
,   ( 21,   2853,   3330,    168 )    # Cnc
,   ( 45,   3330,   3870,    168 )    # Leo
,   ( 58,   6570, 6718.4,    150 )    # Oph
,   (  3, 6718.4,   6792,    150 )    # Aql
,   ( 31,   7500,   7515,    144 )    # Del
,   ( 20,   2520,   2526,    132 )    # CMi
,   ( 73,   6570,   6633,    108 )    # Ser
,   ( 39,   5790,   6030,     96 )    # Her
,   ( 58,   6570,   6633,     72 )    # Oph
,   ( 61,   7728,   7800,     66 )    # Peg
,   ( 66,      0,    720,     48 )    # Psc
,   ( 73,   6690,   6792,     48 )    # Ser
,   ( 31,   7308,   7500,     48 )    # Del
,   ( 34,   7500,   7680,     48 )    # Equ
,   ( 61,   7680,   7728,     48 )    # Peg
,   ( 61,   7920,   8190,     48 )    # Peg
,   ( 61,   7800,   7920,     42 )    # Peg
,   ( 20,   2526,   2592,     36 )    # CMi
,   ( 77,   1290,   1662,      0 )    # Tau
,   ( 59,   1662,   1680,      0 )    # Ori
,   ( 20,   2592,   2910,      0 )    # CMi
,   ( 85,   5280,   5430,      0 )    # Vir
,   ( 58,   6420,   6570,      0 )    # Oph
,   ( 16,    954,   1182,    -42 )    # Cet
,   ( 77,   1182,   1290,    -42 )    # Tau
,   ( 73,   5430,   5856,    -78 )    # Ser
,   ( 59,   1680,   1830,    -96 )    # Ori
,   ( 59,   2100,   2247,    -96 )    # Ori
,   ( 73,   6420,   6468,    -96 )    # Ser
,   ( 73,   6570,   6690,    -96 )    # Ser
,   (  3,   6690,   6792,    -96 )    # Aql
,   ( 66,   8190,   8580,    -96 )    # Psc
,   ( 45,   3870,   4146,   -144 )    # Leo
,   ( 85,   4146,   4260,   -144 )    # Vir
,   ( 66,      0,    120,   -168 )    # Psc
,   ( 66,   8580,   8640,   -168 )    # Psc
,   ( 85,   5130,   5280,   -192 )    # Vir
,   ( 58,   5730,   5856,   -192 )    # Oph
,   (  3,   7200,   7392,   -216 )    # Aql
,   (  4,   7680,   7872,   -216 )    # Aqr
,   ( 58,   6180,   6468,   -240 )    # Oph
,   ( 54,   2100,   2910,   -264 )    # Mon
,   ( 35,   1770,   1830,   -264 )    # Eri
,   ( 59,   1830,   2100,   -264 )    # Ori
,   ( 41,   2910,   3012,   -264 )    # Hya
,   ( 74,   3450,   3870,   -264 )    # Sex
,   ( 85,   4260,   4620,   -264 )    # Vir
,   ( 58,   6330,   6360,   -280 )    # Oph
,   (  3,   6792,   7200, -288.8 )    # Aql
,   ( 35,   1740,   1770,   -348 )    # Eri
,   (  4,   7392,   7680,   -360 )    # Aqr
,   ( 73,   6180,   6570,   -384 )    # Ser
,   ( 72,   6570,   6792,   -384 )    # Sct
,   ( 41,   3012,   3090,   -408 )    # Hya
,   ( 58,   5856,   5895,   -438 )    # Oph
,   ( 41,   3090,   3270,   -456 )    # Hya
,   ( 26,   3870,   3900,   -456 )    # Crt
,   ( 71,   5856,   5895,   -462 )    # Sco
,   ( 47,   5640,   5730,   -480 )    # Lib
,   ( 28,   4530,   4620,   -528 )    # Crv
,   ( 85,   4620,   5130,   -528 )    # Vir
,   ( 41,   3270,   3510,   -576 )    # Hya
,   ( 16,    600,    954, -585.2 )    # Cet
,   ( 35,    954,   1350, -585.2 )    # Eri
,   ( 26,   3900,   4260,   -588 )    # Crt
,   ( 28,   4260,   4530,   -588 )    # Crv
,   ( 47,   5130,   5370,   -588 )    # Lib
,   ( 58,   5856,   6030,   -590 )    # Oph
,   ( 16,      0,    600,   -612 )    # Cet
,   ( 11,   7680,   7872,   -612 )    # Cap
,   (  4,   7872,   8580,   -612 )    # Aqr
,   ( 16,   8580,   8640,   -612 )    # Cet
,   ( 41,   3510,   3690,   -636 )    # Hya
,   ( 35,   1692,   1740,   -654 )    # Eri
,   ( 46,   1740,   2202,   -654 )    # Lep
,   ( 11,   7200,   7680,   -672 )    # Cap
,   ( 41,   3690,   3810,   -700 )    # Hya
,   ( 41,   4530,   5370,   -708 )    # Hya
,   ( 47,   5370,   5640,   -708 )    # Lib
,   ( 71,   5640,   5760,   -708 )    # Sco
,   ( 35,   1650,   1692,   -720 )    # Eri
,   ( 58,   6030,   6336,   -720 )    # Oph
,   ( 76,   6336,   6420,   -720 )    # Sgr
,   ( 41,   3810,   3900,   -748 )    # Hya
,   ( 19,   2202,   2652,   -792 )    # CMa
,   ( 41,   4410,   4530,   -792 )    # Hya
,   ( 41,   3900,   4410,   -840 )    # Hya
,   ( 36,   1260,   1350,   -864 )    # For
,   ( 68,   3012,   3372,   -882 )    # Pyx
,   ( 35,   1536,   1650,   -888 )    # Eri
,   ( 76,   6420,   6900,   -888 )    # Sgr
,   ( 65,   7680,   8280,   -888 )    # PsA
,   ( 70,   8280,   8400,   -888 )    # Scl
,   ( 36,   1080,   1260,   -950 )    # For
,   (  1,   3372,   3960,   -954 )    # Ant
,   ( 70,      0,    600,   -960 )    # Scl
,   ( 36,    600,   1080,   -960 )    # For
,   ( 35,   1392,   1536,   -960 )    # Eri
,   ( 70,   8400,   8640,   -960 )    # Scl
,   ( 14,   5100,   5370,  -1008 )    # Cen
,   ( 49,   5640,   5760,  -1008 )    # Lup
,   ( 71,   5760, 5911.5,  -1008 )    # Sco
,   (  9,   1740,   1800,  -1032 )    # Cae
,   ( 22,   1800,   2370,  -1032 )    # Col
,   ( 67,   2880,   3012,  -1032 )    # Pup
,   ( 35,   1230,   1392,  -1056 )    # Eri
,   ( 71, 5911.5,   6420,  -1092 )    # Sco
,   ( 24,   6420,   6900,  -1092 )    # CrA
,   ( 76,   6900,   7320,  -1092 )    # Sgr
,   ( 53,   7320,   7680,  -1092 )    # Mic
,   ( 35,   1080,   1230,  -1104 )    # Eri
,   (  9,   1620,   1740,  -1116 )    # Cae
,   ( 49,   5520,   5640,  -1152 )    # Lup
,   ( 63,      0,    840,  -1156 )    # Phe
,   ( 35,    960,   1080,  -1176 )    # Eri
,   ( 40,   1470,   1536,  -1176 )    # Hor
,   (  9,   1536,   1620,  -1176 )    # Cae
,   ( 38,   7680,   7920,  -1200 )    # Gru
,   ( 67,   2160,   2880,  -1218 )    # Pup
,   ( 84,   2880,   2940,  -1218 )    # Vel
,   ( 35,    870,    960,  -1224 )    # Eri
,   ( 40,   1380,   1470,  -1224 )    # Hor
,   ( 63,      0,    660,  -1236 )    # Phe
,   ( 12,   2160,   2220,  -1260 )    # Car
,   ( 84,   2940,   3042,  -1272 )    # Vel
,   ( 40,   1260,   1380,  -1276 )    # Hor
,   ( 32,   1380,   1440,  -1276 )    # Dor
,   ( 63,      0,    570,  -1284 )    # Phe
,   ( 35,    780,    870,  -1296 )    # Eri
,   ( 64,   1620,   1800,  -1296 )    # Pic
,   ( 49,   5418,   5520,  -1296 )    # Lup
,   ( 84,   3042,   3180,  -1308 )    # Vel
,   ( 12,   2220,   2340,  -1320 )    # Car
,   ( 14,   4260,   4620,  -1320 )    # Cen
,   ( 49,   5100,   5418,  -1320 )    # Lup
,   ( 56,   5418,   5520,  -1320 )    # Nor
,   ( 32,   1440,   1560,  -1356 )    # Dor
,   ( 84,   3180,   3960,  -1356 )    # Vel
,   ( 14,   3960,   4050,  -1356 )    # Cen
,   (  5,   6300,   6480,  -1368 )    # Ara
,   ( 78,   6480,   7320,  -1368 )    # Tel
,   ( 38,   7920,   8400,  -1368 )    # Gru
,   ( 40,   1152,   1260,  -1380 )    # Hor
,   ( 64,   1800,   1980,  -1380 )    # Pic
,   ( 12,   2340,   2460,  -1392 )    # Car
,   ( 63,      0,    480,  -1404 )    # Phe
,   ( 35,    480,    780,  -1404 )    # Eri
,   ( 63,   8400,   8640,  -1404 )    # Phe
,   ( 32,   1560,   1650,  -1416 )    # Dor
,   ( 56,   5520, 5911.5,  -1440 )    # Nor
,   ( 43,   7320,   7680,  -1440 )    # Ind
,   ( 64,   1980,   2160,  -1464 )    # Pic
,   ( 18,   5460,   5520,  -1464 )    # Cir
,   (  5, 5911.5,   5970,  -1464 )    # Ara
,   ( 18,   5370,   5460,  -1526 )    # Cir
,   (  5,   5970,   6030,  -1526 )    # Ara
,   ( 64,   2160,   2460,  -1536 )    # Pic
,   ( 12,   2460,   3252,  -1536 )    # Car
,   ( 14,   4050,   4260,  -1536 )    # Cen
,   ( 27,   4260,   4620,  -1536 )    # Cru
,   ( 14,   4620,   5232,  -1536 )    # Cen
,   ( 18,   4860,   4920,  -1560 )    # Cir
,   (  5,   6030,   6060,  -1560 )    # Ara
,   ( 40,    780,   1152,  -1620 )    # Hor
,   ( 69,   1152,   1650,  -1620 )    # Ret
,   ( 18,   5310,   5370,  -1620 )    # Cir
,   (  5,   6060,   6300,  -1620 )    # Ara
,   ( 60,   6300,   6480,  -1620 )    # Pav
,   ( 81,   7920,   8400,  -1620 )    # Tuc
,   ( 32,   1650,   2370,  -1680 )    # Dor
,   ( 18,   4920,   5310,  -1680 )    # Cir
,   ( 79,   5310,   6120,  -1680 )    # TrA
,   ( 81,      0,    480,  -1800 )    # Tuc
,   ( 42,   1260,   1650,  -1800 )    # Hyi
,   ( 86,   2370,   3252,  -1800 )    # Vol
,   ( 12,   3252,   4050,  -1800 )    # Car
,   ( 55,   4050,   4920,  -1800 )    # Mus
,   ( 60,   6480,   7680,  -1800 )    # Pav
,   ( 43,   7680,   8400,  -1800 )    # Ind
,   ( 81,   8400,   8640,  -1800 )    # Tuc
,   ( 81,    270,    480,  -1824 )    # Tuc
,   ( 42,      0,   1260,  -1980 )    # Hyi
,   ( 17,   2760,   4920,  -1980 )    # Cha
,   (  2,   4920,   6480,  -1980 )    # Aps
,   ( 52,   1260,   2760,  -2040 )    # Men
,   ( 57,      0,   8640,  -2160 )    # Oct
)




def Constellation(ra, dec):
    """Determines the constellation that contains the given point in the sky.

    Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
    constellation that contains that point.

    Parameters
    ----------
    ra : float
        The right ascension (RA) of a point in the sky, using the J2000 equatorial system.
    dec : float
        The declination (DEC) of a point in the sky, using the J2000 equatorial system.

    Returns
    -------
    ConstellationInfo
        A structure that contains the 3-letter abbreviation and full name
        of the constellation that contains the given (ra,dec), along with
        the converted B1875 (ra,dec) for that point.
    """
    global _ConstelRot, _Epoch2000

    if dec < -90.0 or dec > +90.0:
        raise Error('Invalid declination angle. Must be -90..+90.')

    # Clamp right ascension to [0, 24) sidereal hours.
    ra = math.fmod(ra, 24.0)
    if ra < 0.0:
        ra += 24.0

    # Lazy-initialize rotation matrix.
    if _ConstelRot is None:
        # Need to calculate the B1875 epoch. Based on this:
        # https://en.wikipedia.org/wiki/Epoch_(astronomy)#Besselian_years
        # B = 1900 + (JD - 2415020.31352) / 365.242198781
        # I'm interested in using TT instead of JD, giving:
        # B = 1900 + ((TT+2451545) - 2415020.31352) / 365.242198781
        # B = 1900 + (TT + 36524.68648) / 365.242198781
        # TT = 365.242198781*(B - 1900) - 36524.68648 = -45655.741449525
        # But the Time constructor wants UT, not TT.
        # Near that date, I get a historical correction of ut-tt = 3.2 seconds.
        # That gives UT = -45655.74141261017 for the B1875 epoch,
        # or 1874-12-31T18:12:21.950Z.
        _ConstelRot = Rotation_EQJ_EQD(Time(-45655.74141261017))
        _Epoch2000 = Time(0.0)

    # Convert coordinates from J2000 to B1875.
    sph2000 = Spherical(dec, 15.0 * ra, 1.0)
    vec2000 = VectorFromSphere(sph2000, _Epoch2000)
    vec1875 = RotateVector(_ConstelRot, vec2000)
    equ1875 = EquatorFromVector(vec1875)

    # Convert DEC from degrees and RA from hours, into compact angle units used in the _ConstelBounds table.
    x_dec = 24.0 * equ1875.dec
    x_ra = (24.0 * 15.0) * equ1875.ra

    # Search for the constellation using the B1875 coordinates.
    for b in _ConstelBounds:
        index, ra_lo, ra_hi, b_dec = b
        if (b_dec <= x_dec) and (ra_lo <= x_ra < ra_hi):
            symbol, name = _ConstelNames[index]
            return ConstellationInfo(symbol, name, equ1875.ra, equ1875.dec)

    # This should never happen!
    raise Error('Unable to find constellation for given coordinates.')


@enum.unique
class EclipseKind(enum.Enum):
    """The different kinds of lunar/solar eclipses.

    Values
    ------
    Invalid: No eclipse found.
    Penumbral: A penumbral lunar eclipse. (Never used for a solar eclipse.)
    Partial: A partial lunar/solar eclipse.
    Annular: An annular solar eclipse. (Never used for a lunar eclipse.)
    Total: A total lunar/solar eclipse.
    """
    Invalid = 0
    Penumbral = 1
    Partial = 2
    Annular = 3
    Total = 4


class _ShadowInfo:
    def __init__(self, time, u, r, k, p, target, direction):
        self.time = time
        self.u = u   # dot product of (heliocentric earth) and (geocentric moon): defines the shadow plane where the Moon is
        self.r = r   # km distance between center of Moon and the line passing through the centers of the Sun and Earth.
        self.k = k   # umbra radius in km, at the shadow plane
        self.p = p   # penumbra radius in km, at the shadow plane
        self.target = target        # vector from center of shadow-casting body to target location that might receive the shadow
        self.dir = direction        # vector from center of Sun to center of shadow-casting body


def _CalcShadow(body_radius_km, time, target, sdir):
    u = (sdir.x*target.x + sdir.y*target.y + sdir.z*target.z) / (sdir.x*sdir.x + sdir.y*sdir.y + sdir.z*sdir.z)
    dx = (u * sdir.x) - target.x
    dy = (u * sdir.y) - target.y
    dz = (u * sdir.z) - target.z
    r = KM_PER_AU * math.sqrt(dx*dx + dy*dy + dz*dz)
    k = +_SUN_RADIUS_KM - (1.0 + u)*(_SUN_RADIUS_KM - body_radius_km)
    p = -_SUN_RADIUS_KM + (1.0 + u)*(_SUN_RADIUS_KM + body_radius_km)
    return _ShadowInfo(time, u, r, k, p, target, sdir)


def _EarthShadow(time):
    # This function helps find when the Earth's shadow falls upon the Moon.
    # Light-travel and aberration corrected vector from the Earth to the Sun.
    # The negative vector -s is thus the path of sunlight through the center of the Earth.
    s = GeoVector(Body.Sun, time, True)
    m = GeoMoon(time)
    return _CalcShadow(_EARTH_ECLIPSE_RADIUS_KM, time, m, -s)


def _MoonShadow(time):
    s = GeoVector(Body.Sun, time, True)
    m = GeoMoon(time)       # geocentric Moon
    # -m  = lunacentric Earth
    # m-s = heliocentric Moon
    return _CalcShadow(_MOON_MEAN_RADIUS_KM, time, -m, m-s)


def _LocalMoonShadow(time, observer):
    # Calculate observer's geocentric position.
    pos = _geo_pos(time, observer)
    s = GeoVector(Body.Sun, time, True)
    m = GeoMoon(time)        # geocentric Moon
    # Calculate lunacentric location of an observer on the Earth's surface.
    lo = Vector(pos[0] - m.x, pos[1] - m.y, pos[2] - m.z, time)
    # m-s = heliocentric Moon
    return _CalcShadow(_MOON_MEAN_RADIUS_KM, time, lo, m-s)


def _PlanetShadow(body, planet_radius_km, time):
    # Calculate light-travel-corrected vector from Earth to planet.
    p = GeoVector(body, time, True)
    # Calculate light-travel-corrected vector from Earth to Sun.
    s = GeoVector(Body.Sun, time, True)
    # -p  = planetocentric Earth
    # p-s = heliocentric planet
    return _CalcShadow(planet_radius_km, time, -p, p-s)


def _ShadowDistanceSlope(shadowfunc, time):
    dt = 1.0 / 86400.0
    t1 = time.AddDays(-dt)
    t2 = time.AddDays(+dt)
    shadow1 = shadowfunc(t1)
    shadow2 = shadowfunc(t2)
    return (shadow2.r - shadow1.r) / dt


def _PlanetShadowSlope(context, time):
    (body, planet_radius_km) = context
    dt = 1.0 / 86400.0
    shadow1 = _PlanetShadow(body, planet_radius_km, time.AddDays(-dt))
    shadow2 = _PlanetShadow(body, planet_radius_km, time.AddDays(+dt))
    return (shadow2.r - shadow1.r) / dt


def _PeakEarthShadow(search_center_time):
    window = 0.03        # initial search window, in days, before/after given time
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_ShadowDistanceSlope, _EarthShadow, t1, t2, 1.0)
    return _EarthShadow(tx)


def _PeakMoonShadow(search_center_time):
    window = 0.03        # initial search window, in days, before/after given time
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_ShadowDistanceSlope, _MoonShadow, t1, t2, 1.0)
    return _MoonShadow(tx)

def _PeakLocalMoonShadow(search_center_time, observer):
    # Search for the time near search_center_time that the Moon's shadow comes
    # closest to the given observer.
    window = 0.2
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_ShadowDistanceSlope, lambda time: _LocalMoonShadow(time, observer), t1, t2, 1.0)
    return _LocalMoonShadow(tx, observer)

def _PeakPlanetShadow(body, planet_radius_km, search_center_time):
    # Search for when the body's shadow is closest to the center of the Earth.
    window = 1.0     # days before/after inferior conjunction to search for minimum shadow distance.
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_PlanetShadowSlope, (body, planet_radius_km), t1, t2, 1.0)
    return _PlanetShadow(body, planet_radius_km, tx)

class _ShadowDiffContext:
    def __init__(self, radius_limit, direction):
        self.radius_limit = radius_limit
        self.direction = direction

def _ShadowDiff(context, time):
    return context.direction * (_EarthShadow(time).r - context.radius_limit)


class LunarEclipseInfo:
    """Returns information about a lunar eclipse.

    Returned by #SearchLunarEclipse or #NextLunarEclipse
    to report information about a lunar eclipse event.
    When a lunar eclipse is found, it is classified as penumbral, partial, or total.
    Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
    by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
    Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
    Total eclipses occur when the entire Moon passes into the Earth's umbra.

    The `kind` field thus holds one of the values `EclipseKind.Penumbral`, `EclipseKind.Partial`,
    or `EclipseKind.Total`, depending on the kind of lunar eclipse found.

    The `obscuration` field holds a value in the range [0, 1] that indicates what fraction
    of the Moon's apparent disc area is covered by the Earth's umbra at the eclipse's peak.
    This indicates how dark the peak eclipse appears. For penumbral eclipses, the obscuration
    is 0, because the Moon does not pass through the Earth's umbra. For partial eclipses,
    the obscuration is somewhere between 0 and 1. For total lunar eclipses, the obscuration is 1.

    Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.

    Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
    of the eclipse, which is half of the amount of time the eclipse spends in each
    phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
    By converting from minutes to days, and subtracting/adding with `peak`, the caller
    may determine the date and time of the beginning/end of each eclipse phase.

    Attributes
    ----------
    kind : EclipseKind
         The type of lunar eclipse found.
    obscuration : float
        The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra.
    peak : Time
         The time of the eclipse at its peak.
    sd_penum : float
         The semi-duration of the penumbral phase in minutes.
    sd_partial : float
         The semi-duration of the penumbral phase in minutes, or 0.0 if none.
    sd_total : float
         The semi-duration of the penumbral phase in minutes, or 0.0 if none.
    """
    def __init__(self, kind, obscuration, peak, sd_penum, sd_partial, sd_total):
        self.kind = kind
        self.obscuration = obscuration
        self.peak = peak
        self.sd_penum = sd_penum
        self.sd_partial = sd_partial
        self.sd_total = sd_total

    def __repr__(self):
        return 'LunarEclipseInfo({}, obscuration={}, peak={}, sd_penum={}, sd_partial={}, sd_total={})'.format(
            self.kind,
            self.obscuration,
            repr(self.peak),
            self.sd_penum,
            self.sd_partial,
            self.sd_total
        )


class GlobalSolarEclipseInfo:
    """Reports the time and geographic location of the peak of a solar eclipse.

    Returned by #SearchGlobalSolarEclipse or #NextGlobalSolarEclipse
    to report information about a solar eclipse event.

    The eclipse is classified as partial, annular, or total, depending on the
    maximum amount of the Sun's disc obscured, as seen at the peak location
    on the surface of the Earth.

    The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
    An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
    to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
    A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
    observes either a total or annular eclipse.

    If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
    fields give the geographic coordinates of the center of the Moon's shadow projected
    onto the daytime side of the Earth at the instant of the eclipse's peak.
    If `kind` has any other value, `latitude` and `longitude` are undefined and should
    not be used.

    For total or annular eclipses, the `obscuration` field holds the fraction (0, 1]
    of the Sun's apparent disc area that is blocked from view by the Moon's silhouette,
    as seen by an observer located at the geographic coordinates `latitude`, `longitude`
    at the darkest time `peak`. The value will always be 1 for total eclipses, and less than
    1 for annular eclipses.
    For partial eclipses, `obscuration` holds the value `None`.
    This is because there is little practical use for an obscuration value of
    a partial eclipse without supplying a particular observation location.
    Developers who wish to find an obscuration value for partial solar eclipses should therefore use
    #SearchLocalSolarEclipse and provide the geographic coordinates of an observer.

    Attributes
    ----------
    kind : EclipseKind
        The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    obscuration : float
        The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only).
    peak : Time
        The date and time when the solar eclipse is darkest.
        This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.
    distance : float
        The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.
    latitude : float
        The geographic latitude at the center of the peak eclipse shadow.
    longitude : float
        The geographic longitude at the center of the peak eclipse shadow.
    """
    def __init__(self, kind, obscuration, peak, distance, latitude, longitude):
        self.kind = kind
        self.obscuration = obscuration
        self.peak = peak
        self.distance = distance
        self.latitude = latitude
        self.longitude = longitude

    def __repr__(self):
        return 'GlobalSolarEclipseInfo({}, obscuration={}, peak={}, distance={}, latitude={}, longitude={})'.format(
            self.kind,
            self.obscuration,
            repr(self.peak),
            self.distance,
            self.latitude,
            self.longitude
        )


class EclipseEvent:
    """Holds a time and the observed altitude of the Sun at that time.

    When reporting a solar eclipse observed at a specific location on the Earth
    (a "local" solar eclipse), a series of events occur. In addition
    to the time of each event, it is important to know the altitude of the Sun,
    because each event may be invisible to the observer if the Sun is below
    the horizon.

    If `altitude` is negative, the event is theoretical only; it would be
    visible if the Earth were transparent, but the observer cannot actually see it.
    If `altitude` is positive but less than a few degrees, visibility will be impaired by
    atmospheric interference (sunrise or sunset conditions).

    Attributes
    ----------
    time : Time
        The date and time of the event.
    altitude : float
        The angular altitude of the center of the Sun above/below the horizon, at `time`,
        corrected for atmospheric refraction and expressed in degrees.
    """
    def __init__(self, time, altitude):
        self.time = time
        self.altitude = altitude

    def __repr__(self):
        return 'EclipseEvent({}, altitude={})'.format(
            repr(self.time),
            self.altitude
        )


class LocalSolarEclipseInfo:
    """Information about a solar eclipse as seen by an observer at a given time and geographic location.

    Returned by #SearchLocalSolarEclipse or #NextLocalSolarEclipse
    to report information about a solar eclipse as seen at a given geographic location.

    When a solar eclipse is found, it is classified as partial, annular, or total.
    The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    A partial solar eclipse is when the Moon does not line up directly enough with the Sun
    to completely block the Sun's light from reaching the observer.
    An annular eclipse occurs when the Moon's disc is completely visible against the Sun
    but the Moon is too far away to completely block the Sun's light; this leaves the
    Sun with a ring-like appearance.
    A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
    Sun just right to completely block all sunlight from reaching the observer.

    The `obscuration` field reports what fraction of the Sun's disc appears blocked
    by the Moon when viewed by the observer at the peak eclipse time.
    This is a value that ranges from 0 (no blockage) to 1 (total eclipse).
    The obscuration value will be between 0 and 1 for partial eclipses and annular eclipses.
    The value will be exactly 1 for total eclipses. Obscuration gives an indication
    of how dark the eclipse appears.

    There are 5 "event" fields, each of which contains a time and a solar altitude.
    Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
    The fields `partial_begin` and `partial_end` are always set, and indicate when
    the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
    `total_begin` and `total_end` indicate when the total/annular phase begins/ends.
    When an event field is valid, the caller must also check its `altitude` field to
    see whether the Sun is above the horizon at the time indicated by the `time` field.
    See #EclipseEvent for more information.

    Attributes
    ----------
    kind : EclipseKind
        The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    obscuration : float
        The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak.
    partial_begin : EclipseEvent
        The time and Sun altitude at the beginning of the eclipse.
    total_begin : EclipseEvent
        If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise `None`.
    peak : EclipseEvent
        The time and Sun altitude when the eclipse reaches its peak.
    total_end : EclipseEvent
        If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise `None`.
    partial_end : EclipseEvent
        The time and Sun altitude at the end of the eclipse.
    """
    def __init__(self, kind, obscuration, partial_begin, total_begin, peak, total_end, partial_end):
        self.kind = kind
        self.obscuration = obscuration
        self.partial_begin = partial_begin
        self.total_begin = total_begin
        self.peak = peak
        self.total_end = total_end
        self.partial_end = partial_end

    def __repr__(self):
        return 'LocalSolarEclipseInfo({}, obscuration={}, partial_begin={}, total_begin={}, peak={}, total_end={}, partial_end={})'.format(
            self.kind,
            self.obscuration,
            repr(self.partial_begin),
            repr(self.total_begin),
            repr(self.peak),
            repr(self.total_end),
            repr(self.partial_end)
        )


def _EclipseKindFromUmbra(k):
    # The umbra radius tells us what kind of eclipse the observer sees.
    # If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular.
    # HACK: I added a tiny bias (14 meters) to match Espenak test data.
    if k > 0.014:
        return EclipseKind.Total
    return EclipseKind.Annular

class _LocalTransitionContext:
    def __init__(self, observer, direction, func):
        self.observer = observer
        self.direction = direction
        self.func = func


def _LocalTransitionFunc(context, time):
    shadow = _LocalMoonShadow(time, context.observer)
    return context.direction * context.func(shadow)


def _LocalEclipseTransition(observer, direction, func, t1, t2):
    context = _LocalTransitionContext(observer, direction, func)
    search = Search(_LocalTransitionFunc, context, t1, t2, 1.0)
    if search is None:
        raise Error('Local eclipse transition search failed')
    return _CalcEvent(observer, search)


def _CalcEvent(observer, time):
    altitude = _SunAltitude(time, observer)
    return EclipseEvent(time, altitude)


def _SunAltitude(time, observer):
    equ = Equator(Body.Sun, time, observer, True, True)
    hor = Horizon(time, observer, equ.ra, equ.dec, Refraction.Normal)
    return hor.altitude


def _local_partial_distance(shadow):
    # Must take the absolute value of the umbra radius 'k'
    # because it can be negative for an annular eclipse.
    return shadow.p - shadow.r


def _local_total_distance(shadow):
    return abs(shadow.k) - shadow.r


def _LocalEclipse(shadow, observer):
    PARTIAL_WINDOW = 0.2
    TOTAL_WINDOW = 0.01
    peak = _CalcEvent(observer, shadow.time)
    t1 = shadow.time.AddDays(-PARTIAL_WINDOW)
    t2 = shadow.time.AddDays(+PARTIAL_WINDOW)
    partial_begin = _LocalEclipseTransition(observer, +1.0, _local_partial_distance, t1, shadow.time)
    partial_end   = _LocalEclipseTransition(observer, -1.0, _local_partial_distance, shadow.time, t2)
    if shadow.r < abs(shadow.k):    # take absolute value of 'k' to handle annular eclipses too.
        t1 = shadow.time.AddDays(-TOTAL_WINDOW)
        t2 = shadow.time.AddDays(+TOTAL_WINDOW)
        total_begin = _LocalEclipseTransition(observer, +1.0, _local_total_distance, t1, shadow.time)
        total_end   = _LocalEclipseTransition(observer, -1.0, _local_total_distance, shadow.time, t2)
        kind = _EclipseKindFromUmbra(shadow.k)
    else:
        total_begin = None
        total_end = None
        kind = EclipseKind.Partial
    obscuration = 1.0 if (kind == EclipseKind.Total) else _SolarEclipseObscuration(shadow.dir, shadow.target)
    return LocalSolarEclipseInfo(kind, obscuration, partial_begin, total_begin, peak, total_end, partial_end)


def _GeoidIntersect(shadow):
    kind = EclipseKind.Partial
    peak = shadow.time
    distance = shadow.r
    latitude = longitude = math.nan

    # We want to calculate the intersection of the shadow axis with the Earth's geoid.
    # First we must convert EQJ (equator of J2000) coordinates to EQD (equator of date)
    # coordinates that are perfectly aligned with the Earth's equator at this
    # moment in time.
    rot = Rotation_EQJ_EQD(shadow.time)
    v = RotateVector(rot, shadow.dir)        # shadow-axis vector in equator-of-date coordinates
    e = RotateVector(rot, shadow.target)     # lunacentric Earth in equator-of-date coordinates

    # Convert all distances from AU to km.
    # But dilate the z-coordinates so that the Earth becomes a perfect sphere.
    # Then find the intersection of the vector with the sphere.
    # See p 184 in Montenbruck & Pfleger's "Astronomy on the Personal Computer", second edition.
    v.x *= KM_PER_AU
    v.y *= KM_PER_AU
    v.z *= KM_PER_AU / _EARTH_FLATTENING
    e.x *= KM_PER_AU
    e.y *= KM_PER_AU
    e.z *= KM_PER_AU / _EARTH_FLATTENING

    # Solve the quadratic equation that finds whether and where
    # the shadow axis intersects with the Earth in the dilated coordinate system.
    R = _EARTH_EQUATORIAL_RADIUS_KM
    A = v.x*v.x + v.y*v.y + v.z*v.z
    B = -2.0 * (v.x*e.x + v.y*e.y + v.z*e.z)
    C = (e.x*e.x + e.y*e.y + e.z*e.z) - R*R
    radic = B*B - 4*A*C

    if radic > 0.0:
        # Calculate the closer of the two intersection points.
        # This will be on the day side of the Earth.
        u = (-B - math.sqrt(radic)) / (2 * A)

        # Convert lunacentric dilated coordinates to geocentric coordinates.
        px = u*v.x - e.x
        py = u*v.y - e.y
        pz = (u*v.z - e.z) * _EARTH_FLATTENING

        # Convert cartesian coordinates into geodetic latitude/longitude.
        proj = math.hypot(px, py) * _EARTH_FLATTENING_SQUARED
        if proj == 0.0:
            if pz > 0.0:
                latitude = +90.0
            else:
                latitude = -90.0
        else:
            latitude = math.degrees(math.atan(pz / proj))

        # Adjust longitude for Earth's rotation at the given UT.
        gast = SiderealTime(peak)
        longitude = math.fmod(math.degrees(math.atan2(py, px)) - (15*gast), 360.0)
        if longitude <= -180.0:
            longitude += 360.0
        elif longitude > +180.0:
            longitude -= 360.0

        # We want to determine whether the observer sees a total eclipse or an annular eclipse.
        # We need to perform a series of vector calculations...
        # Calculate the inverse rotation matrix, so we can convert EQD to EQJ.
        inv = InverseRotation(rot)

        # Put the EQD geocentric coordinates of the observer into the vector 'o'.
        # Also convert back from kilometers to astronomical units.
        o = Vector(px / KM_PER_AU, py / KM_PER_AU, pz / KM_PER_AU, shadow.time)

        # Rotate the observer's geocentric EQD back to the EQJ system.
        o = RotateVector(inv, o)

        # Convert geocentric vector to lunacentric vector.
        o.x += shadow.target.x
        o.y += shadow.target.y
        o.z += shadow.target.z

        # Recalculate the shadow using a vector from the Moon's center toward the observer.
        surface = _CalcShadow(_MOON_POLAR_RADIUS_KM, shadow.time, o, shadow.dir)

        # If we did everything right, the shadow distance should be very close to zero.
        # That's because we already determined the observer 'o' is on the shadow axis!
        if surface.r > 1.0e-9 or surface.r < 0.0:
            raise Error('Unexpected shadow distance from geoid intersection = {}'.format(surface.r))

        kind = _EclipseKindFromUmbra(surface.k)
        obscuration = 1.0 if (kind == EclipseKind.Total) else _SolarEclipseObscuration(shadow.dir, o)
    else:
        # This is a partial solar eclipse. It does not make practical sense to calculate obscuration.
        # Anyone who wants obscuration should use SearchLocalSolarEclipse for a specific location on the Earth.
        obscuration = None

    return GlobalSolarEclipseInfo(kind, obscuration, peak, distance, latitude, longitude)


def _ShadowSemiDurationMinutes(center_time, radius_limit, window_minutes):
    # Search backwards and forwards from the center time until shadow axis distance crosses radius limit.
    window = window_minutes / (24.0 * 60.0)
    before = center_time.AddDays(-window)
    after  = center_time.AddDays(+window)
    t1 = Search(_ShadowDiff, _ShadowDiffContext(radius_limit, -1.0), before, center_time, 1.0)
    t2 = Search(_ShadowDiff, _ShadowDiffContext(radius_limit, +1.0), center_time, after, 1.0)
    if (t1 is None) or (t2 is None):
        raise Error('Failed to find shadow semiduration')
    return (t2.ut - t1.ut) * ((24.0 * 60.0) / 2.0)   # convert days to minutes and average the semi-durations.


def _MoonEclipticLatitudeDegrees(time):
    moon = _CalcMoon(time)
    return math.degrees(moon.geo_eclip_lat)

def _Obscuration(a, b, c):
    # a = radius of first disc
    # b = radius of second disc
    # c = distance between the centers of the discs
    if a <= 0.0:
        raise Error('Radius of first disc must be positive.')

    if b <= 0.0:
        raise Error('Radius of second disc must be positive.')

    if c < 0.0:
        raise Error('Distance between discs is not allowed to be negative.')

    if c >= a + b:
        # The discs are too far apart to have any overlapping area
        return 0.0

    if c == 0.0:
        # The discs have a common center. Therefore, one disc is inside the other.
        return 1.0 if (a <= b) else (b*b)/(a*a)

    x = (a*a - b*b + c*c) / (2.0*c)
    radicand = a*a - x*x
    if radicand <= 0.0:
        # The circumferences do not intersect, or are tangent.
        # We already ruled out the case of non-overlapping discs.
        # Therefore, one disc is inside the other.
        return 1.0 if (a <= b) else (b*b)/(a*a)

    # The discs overlap fractionally in a pair of lens-shaped areas.
    y = math.sqrt(radicand)

    # Return the overlapping fractional area.
    # There are two lens-shaped areas, one to the left of x, the other to the right of x.
    # Each part is calculated by subtracting a triangular area from a sector's area.
    lens1 = a*a*math.acos(x/a) - x*y
    lens2 = b*b*math.acos((c-x)/b) - (c-x)*y

    # Find the fractional area with respect to the first disc.
    return (lens1 + lens2) / (math.pi*a*a)


def _SolarEclipseObscuration(hm, lo):
    # hm = heliocentric Moon
    # lo = lunacentric observer
    # Find heliocentric observer
    ho = hm + lo
    # Calculate the apparent angular radius of the Sun for the observer.
    sun_radius = math.asin(_SUN_RADIUS_AU / ho.Length())

    # Calculate the apparent angular radius of the Moon for the observer.
    moon_radius = math.asin(_MOON_POLAR_RADIUS_AU / lo.Length())

    # Calculate the apparent angular separation between the Sun's center and the Moon's center.
    sun_moon_separation = math.radians(AngleBetween(lo, ho))

    # Find the fraction of the Sun's apparent disc area that is covered by the Moon.
    obscuration = _Obscuration(sun_radius, moon_radius, sun_moon_separation)

    # HACK: In marginal cases, we need to clamp obscuration to less than 1.0.
    # This function is never called for total eclipses, so it should never return 1.0.
    return min(0.9999, obscuration)


def SearchLunarEclipse(startTime):
    """Searches for a lunar eclipse.

    This function finds the first lunar eclipse that occurs after `startTime`.
    A lunar eclipse may be penumbral, partial, or total.
    See #LunarEclipseInfo for more information.
    To find a series of lunar eclipses, call this function once,
    then keep calling #NextLunarEclipse as many times as desired,
    passing in the `peak` value returned from the previous call.

    Parameters
    ----------
    startTime : Time
         The date and time for starting the search for a lunar eclipse.

    Returns
    -------
    LunarEclipseInfo
    """
    PruneLatitude = 1.8   # full Moon's ecliptic latitude above which eclipse is impossible
    fmtime = startTime
    for _ in range(12):
        # Search for the next full moon. Any eclipse will be near it.
        fullmoon = SearchMoonPhase(180, fmtime, 40)
        if fullmoon is None:
            raise InternalError()   # should have always found the next full moon

        # Pruning: if the full Moon's ecliptic latitude is too large,
        # a lunar eclipse is not possible. Avoid needless work searching for
        # the minimum moon distance.
        eclip_lat = _MoonEclipticLatitudeDegrees(fullmoon)
        if abs(eclip_lat) < PruneLatitude:
            # Search near the full moon for the time when the center of the Moon
            # is closest to the line passing through the centers of the Sun and Earth.
            shadow = _PeakEarthShadow(fullmoon)
            if shadow.r < shadow.p + _MOON_MEAN_RADIUS_KM:
                # This is at least a penumbral eclipse. We will return a result.
                kind = EclipseKind.Penumbral
                obscuration = 0.0
                sd_total = 0.0
                sd_partial = 0.0
                sd_penum = _ShadowSemiDurationMinutes(shadow.time, shadow.p + _MOON_MEAN_RADIUS_KM, 200.0)

                if shadow.r < shadow.k + _MOON_MEAN_RADIUS_KM:
                    # This is at least a partial eclipse.
                    kind = EclipseKind.Partial
                    sd_partial = _ShadowSemiDurationMinutes(shadow.time, shadow.k + _MOON_MEAN_RADIUS_KM, sd_penum)

                    if shadow.r + _MOON_MEAN_RADIUS_KM < shadow.k:
                        # This is a total eclipse.
                        kind = EclipseKind.Total
                        obscuration = 1.0
                        sd_total = _ShadowSemiDurationMinutes(shadow.time, shadow.k - _MOON_MEAN_RADIUS_KM, sd_partial)
                    else:
                        obscuration = _Obscuration(_MOON_MEAN_RADIUS_KM, shadow.k, shadow.r)

                return LunarEclipseInfo(kind, obscuration, shadow.time, sd_penum, sd_partial, sd_total)

        # We didn't find an eclipse on this full moon, so search for the next one.
        fmtime = fullmoon.AddDays(10)

    # This should never happen because there are always at least 2 full moons per year.
    raise Error('Failed to find lunar eclipse within 12 full moons.')


def NextLunarEclipse(prevEclipseTime):
    """Searches for the next lunar eclipse in a series.

     After using #SearchLunarEclipse to find the first lunar eclipse
     in a series, you can call this function to find the next consecutive lunar eclipse.
     Pass in the `peak` value from the #LunarEclipseInfo returned by the
     previous call to `SearchLunarEclipse` or `NextLunarEclipse`
     to find the next lunar eclipse.

    Parameters
    ----------
    prevEclipseTime : Time
        A date and time near a full moon. Lunar eclipse search will start at the next full moon.

    Returns
    -------
    LunarEclipseInfo
    """
    startTime = prevEclipseTime.AddDays(10.0)
    return SearchLunarEclipse(startTime)


def SearchGlobalSolarEclipse(startTime):
    """Searches for a solar eclipse visible anywhere on the Earth's surface.

    This function finds the first solar eclipse that occurs after `startTime`.
    A solar eclipse may be partial, annular, or total.
    See #GlobalSolarEclipseInfo for more information.
    To find a series of solar eclipses, call this function once,
    then keep calling #NextGlobalSolarEclipse as many times as desired,
    passing in the `peak` value returned from the previous call.

    Parameters
    ----------
    startTime : Time
        The date and time for starting the search for a solar eclipse.

    Returns
    -------
    GlobalSolarEclipseInfo
    """
    PruneLatitude = 1.8     # Moon's ecliptic latitude beyond which eclipse is impossible
    # Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
    nmtime = startTime
    for _ in range(12):
        # Search for the next new moon. Any eclipse will be near it.
        newmoon = SearchMoonPhase(0.0, nmtime, 40.0)
        if newmoon is None:
            raise InternalError()   # should always find the next new moon

        # Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
        eclip_lat = _MoonEclipticLatitudeDegrees(newmoon)
        if abs(eclip_lat) < PruneLatitude:
            # Search near the new moon for the time when the center of the Earth
            # is closest to the line passing through the centers of the Sun and Moon.
            shadow = _PeakMoonShadow(newmoon)
            if shadow.r < shadow.p + _EARTH_MEAN_RADIUS_KM:
                # This is at least a partial solar eclipse visible somewhere on Earth.
                # Try to find an intersection between the shadow axis and the Earth's oblate geoid.
                return _GeoidIntersect(shadow)

        # We didn't find an eclipse on this new moon, so search for the next one.
        nmtime = newmoon.AddDays(10.0)

    # Safety valve to prevent infinite loop.
    # This should never happen, because at least 2 solar eclipses happen per year.
    raise Error('Failed to find solar eclipse within 12 full moons.')


def NextGlobalSolarEclipse(prevEclipseTime):
    """Searches for the next global solar eclipse in a series.

    After using #SearchGlobalSolarEclipse to find the first solar eclipse
    in a series, you can call this function to find the next consecutive solar eclipse.
    Pass in the `peak` value from the #GlobalSolarEclipseInfo returned by the
    previous call to `SearchGlobalSolarEclipse` or `NextGlobalSolarEclipse`
    to find the next solar eclipse.

    Parameters
    ----------
    prevEclipseTime : Time
        A date and time near a new moon. Solar eclipse search will start at the next new moon.

    Returns
    -------
    GlobalSolarEclipseInfo
    """
    startTime = prevEclipseTime.AddDays(10.0)
    return SearchGlobalSolarEclipse(startTime)


def SearchLocalSolarEclipse(startTime, observer):
    """Searches for a solar eclipse visible at a specific location on the Earth's surface.
    This function finds the first solar eclipse that occurs after `startTime`.
    A solar eclipse may be partial, annular, or total.
    See #LocalSolarEclipseInfo for more information.

    To find a series of solar eclipses, call this function once,
    then keep calling #NextLocalSolarEclipse as many times as desired,
    passing in the `peak` value returned from the previous call.

    IMPORTANT: An eclipse reported by this function might be partly or
    completely invisible to the observer due to the time of day.
    See #LocalSolarEclipseInfo for more information about this topic.

    Parameters
    ----------
    startTime : Time
        The date and time for starting the search for a solar eclipse.
    observer : Observer
        The geographic location of the observer.

    Returns
    -------
    LocalSolarEclipseInfo
    """
    PruneLatitude = 1.8   # Moon's ecliptic latitude beyond which eclipse is impossible

    # Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
    nmtime = startTime
    while True:
        # Search for the next new moon. Any eclipse will be near it.
        newmoon = SearchMoonPhase(0.0, nmtime, 40.0)
        if newmoon is None:
            raise InternalError()   # should always find the next new moon

        # Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
        eclip_lat = _MoonEclipticLatitudeDegrees(newmoon)
        if abs(eclip_lat) < PruneLatitude:
            # Search near the new moon for the time when the observer
            # is closest to the line passing through the centers of the Sun and Moon.
            shadow = _PeakLocalMoonShadow(newmoon, observer)
            if shadow.r < shadow.p:
                # This is at least a partial solar eclipse for the observer.
                eclipse = _LocalEclipse(shadow, observer)

                # Ignore any eclipse that happens completely at night.
                # More precisely, the center of the Sun must be above the horizon
                # at the beginning or the end of the eclipse, or we skip the event.
                if eclipse.partial_begin.altitude > 0.0 or eclipse.partial_end.altitude > 0.0:
                    return eclipse

        # We didn't find an eclipse on this new moon, so search for the next one.
        nmtime = newmoon.AddDays(10.0)


def NextLocalSolarEclipse(prevEclipseTime, observer):
    """Searches for the next local solar eclipse in a series.

    After using #SearchLocalSolarEclipse to find the first solar eclipse
    in a series, you can call this function to find the next consecutive solar eclipse.
    Pass in the `peak` value from the #LocalSolarEclipseInfo returned by the
    previous call to `SearchLocalSolarEclipse` or `NextLocalSolarEclipse`
    to find the next solar eclipse.

    Parameters
    ----------
    prevEclipseTime : Time
        A date and time near a new moon. Solar eclipse search will start at the next new moon.
    observer : Observer
        The geographic location of the observer.

    Returns
    -------
    LocalSolarEclipseInfo
    """
    startTime = prevEclipseTime.AddDays(10.0)
    return SearchLocalSolarEclipse(startTime, observer)


class TransitInfo:
    """Information about a transit of Mercury or Venus, as seen from the Earth.

    Returned by #SearchTransit or #NextTransit to report
    information about a transit of Mercury or Venus.
    A transit is when Mercury or Venus passes between the Sun and Earth so that
    the other planet is seen in silhouette against the Sun.

    The calculations are performed from the point of view of a geocentric observer.

    Attributes
    ----------
    start : Time
        The date and time at the beginning of the transit.
        This is the moment the planet first becomes visible against the Sun in its background.
    peak : Time
        When the planet is most aligned with the Sun, as seen from the Earth.
    finish : Time
        The date and time at the end of the transit.
        This is the moment the planet is last seen against the Sun in its background.
    separation : float
        The minimum angular separation, in arcminutes, between the centers of the Sun and the planet.
        This angle pertains to the time stored in `peak`.
    """
    def __init__(self, start, peak, finish, separation):
        self.start = start
        self.peak = peak
        self.finish = finish
        self.separation = separation

    def __repr__(self):
        return 'TransitInfo(start={}, peak={}, finish={}, separation={})'.format(
            repr(self.start),
            repr(self.peak),
            repr(self.finish),
            self.separation
        )


def _PlanetShadowBoundary(context, time):
    (body, planet_radius_km, direction) = context
    shadow = _PlanetShadow(body, planet_radius_km, time)
    return direction * (shadow.r - shadow.p)


def _PlanetTransitBoundary(body, planet_radius_km, t1, t2, direction):
    # Search for the time the planet's penumbra begins/ends making contact with the center of the Earth.
    # context = new SearchContext_PlanetShadowBoundary(body, planet_radius_km, direction);
    tx = Search(_PlanetShadowBoundary, (body, planet_radius_km, direction), t1, t2, 1.0)
    if tx is None:
        raise Error('Planet transit boundary search failed')
    return tx


def SearchTransit(body, startTime):
    """Searches for the first transit of Mercury or Venus after a given date.

    Finds the first transit of Mercury or Venus after a specified date.
    A transit is when an inferior planet passes between the Sun and the Earth
    so that the silhouette of the planet is visible against the Sun in the background.
    To continue the search, pass the `finish` time in the returned structure to
    #NextTransit.

    Parameters
    ----------
    body : Body
        The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
    startTime : Time
        The date and time for starting the search for a transit.

    Returns
    -------
    TransitInfo
    """
    threshold_angle = 0.4     # maximum angular separation to attempt transit calculation
    dt_days = 1.0

    # Validate the planet and find its mean radius.
    if body == Body.Mercury:
        planet_radius_km = 2439.7
    elif body == Body.Venus:
        planet_radius_km = 6051.8
    else:
        raise InvalidBodyError()

    search_time = startTime
    while True:
        # Search for the next inferior conjunction of the given planet.
        # This is the next time the Earth and the other planet have the same
        # ecliptic longitude as seen from the Sun.
        conj = SearchRelativeLongitude(body, 0.0, search_time)

        # Calculate the angular separation between the body and the Sun at this time.
        conj_separation = AngleFromSun(body, conj)

        if conj_separation < threshold_angle:
            # The planet's angular separation from the Sun is small enough
            # to consider it a transit candidate.
            # Search for the moment when the line passing through the Sun
            # and planet are closest to the Earth's center.
            shadow = _PeakPlanetShadow(body, planet_radius_km, conj)
            if shadow.r < shadow.p:      # does the planet's penumbra touch the Earth's center?
                # Find the beginning and end of the penumbral contact.
                time_before = shadow.time.AddDays(-dt_days)
                start = _PlanetTransitBoundary(body, planet_radius_km, time_before, shadow.time, -1.0)
                time_after = shadow.time.AddDays(+dt_days)
                finish = _PlanetTransitBoundary(body, planet_radius_km, shadow.time, time_after, +1.0)
                min_separation = 60.0 * AngleFromSun(body, shadow.time)
                return TransitInfo(start, shadow.time, finish, min_separation)

        # This inferior conjunction was not a transit. Try the next inferior conjunction.
        search_time = conj.AddDays(10.0)


def NextTransit(body, prevTransitTime):
    """Searches for another transit of Mercury or Venus.

    After calling #SearchTransit to find a transit of Mercury or Venus,
    this function finds the next transit after that.
    Keep calling this function as many times as you want to keep finding more transits.

    Parameters
    ----------
    body : Body
        The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
    prevTransitTime : Time
        A date and time near the previous transit.

    Returns
    -------
    TransitInfo
    """
    startTime = prevTransitTime.AddDays(100.0)
    return SearchTransit(body, startTime)


@enum.unique
class NodeEventKind(enum.Enum):
    """Indicates whether a crossing through the ecliptic plane is ascending or descending.

    Values
    ------
    Invalid: A placeholder for an invalid or undefined node.
    Ascending: indicates a body passing through the ecliptic plane from south to north.
    Descending: indicates a body passing through the ecliptic plane from north to south.
    """
    Invalid = 0
    Ascending = +1
    Descending = -1

class NodeEventInfo:
    """Information about an ascending or descending node of a body.

    This object is returned by #SearchMoonNode and #NextMoonNode
    to report information about the center of the Moon passing through the ecliptic plane.

    Attributes
    ----------
    kind : NodeEventKind
        Whether the node is ascending (south to north) or descending (north to south).
    time : Time
        The time when the body passes through the ecliptic plane.
    """
    def __init__(self, kind, time):
        self.kind = kind
        self.time = time

    def __repr__(self):
        return 'NodeEventInfo({}, {})'.format(self.kind, repr(self.time))

_MoonNodeStepDays = +10.0   # a safe number of days to step without missing a Moon node

def _MoonNodeSearchFunc(direction, time):
    return direction * EclipticGeoMoon(time).lat

def SearchMoonNode(startTime):
    """Searches for a time when the Moon's center crosses through the ecliptic plane.

    Searches for the first ascending or descending node of the Moon after `startTime`.
    An ascending node is when the Moon's center passes through the ecliptic plane
    (the plane of the Earth's orbit around the Sun) from south to north.
    A descending node is when the Moon's center passes through the ecliptic plane
    from north to south. Nodes indicate possible times of solar or lunar eclipses,
    if the Moon also happens to be in the correct phase (new or full, respectively).
    Call `SearchMoonNode` to find the first of a series of nodes.
    Then call #NextMoonNode to find as many more consecutive nodes as desired.

    Parameters
    ----------
    startTime : Time
        The date and time for starting the search for an ascending or descending node of the Moon.

    Returns
    -------
    NodeEventInfo
    """
    # Start at the given moment in time and sample the Moon's ecliptic latitude.
    # Step 10 days at a time, searching for an interval where that latitude crosses zero.
    time1 = startTime
    eclip1 = EclipticGeoMoon(time1)
    while True:
        time2 = time1.AddDays(_MoonNodeStepDays)
        eclip2 = EclipticGeoMoon(time2)
        if eclip1.lat * eclip2.lat <= 0.0:
            # There is a node somewhere inside this closed time interval.
            # Figure out whether it is an ascending node or a descending node.
            kind = NodeEventKind.Ascending if (eclip2.lat > eclip1.lat) else NodeEventKind.Descending
            result = Search(_MoonNodeSearchFunc, kind.value, time1, time2, 1.0)
            if result is None:
                raise InternalError()   # should always find the next lunar node
            return NodeEventInfo(kind, result)
        time1 = time2
        eclip1 = eclip2


def NextMoonNode(prevNode):
    """Searches for the next time when the Moon's center crosses through the ecliptic plane.

    Call #SearchMoonNode to find the first of a series of nodes.
    Then call `NextMoonNode` to find as many more consecutive nodes as desired.

    Parameters
    ----------
    prevNode : NodeEventInfo
        The previous node find from calling #SearchMoonNode or `NextMoonNode`.

    Returns
    -------
    NodeEventInfo
    """
    time = prevNode.time.AddDays(_MoonNodeStepDays)
    node = SearchMoonNode(time)
    if prevNode.kind == NodeEventKind.Ascending:
        if node.kind != NodeEventKind.Descending:
            raise InternalError()
    elif prevNode.kind == NodeEventKind.Descending:
        if node.kind != NodeEventKind.Ascending:
            raise InternalError()
    else:
        raise Error('prevNode contains an invalid node kind: {}'.format(prevNode.kind))
    return node


class LibrationInfo:
    """Lunar libration angles, returned by #Libration.

    Contains lunar libration angles and lunar position information
    for a given moment in time. See #Libration for more details.

    Attributes
    ----------
    elat : float
        Sub-Earth libration ecliptic latitude angle, in degrees.
    elon : float
        Sub-Earth libration ecliptic longitude angle, in degrees.
    mlat : float
        Moon's geocentric ecliptic latitude, in degrees.
    mlon : float
        Moon's geocentric ecliptic longitude, in degrees.
    dist_km : float
        Distance between the centers of the Earth and Moon in kilometers.
    diam_deg : float
        The apparent angular diameter of the Moon as seen from the center of the Earth.
    """
    def __init__(self, elat, elon, mlat, mlon, dist_km, diam_deg):
        self.elat = elat
        self.elon = elon
        self.mlat = mlat
        self.mlon = mlon
        self.dist_km = dist_km
        self.diam_deg = diam_deg

    def __repr__(self):
        return 'LibrationInfo(elat={}, elon={}, mlat={}, mlon={}, dist_km={}, diam_deg={})'.format(
            self.elat,
            self.elon,
            self.mlat,
            self.mlon,
            self.dist_km,
            self.diam_deg
        )


def Libration(time):
    """Calculates the Moon's libration angles at a given moment in time.

    Libration is an observed back-and-forth wobble of the portion of the
    Moon visible from the Earth. It is caused by the imperfect tidal locking
    of the Moon's fixed rotation rate, compared to its variable angular speed
    of orbit around the Earth.

    This function calculates a pair of perpendicular libration angles,
    one representing rotation of the Moon in eclitpic longitude `elon`, the other
    in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.

    This function also returns the geocentric position of the Moon
    expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
    distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
    and the apparent angular diameter of the Moon `diam_deg`.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Moon's libration angles.

    Returns
    -------
    LibrationInfo
    """
    t = time.tt / 36525.0
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    moon = _CalcMoon(time)
    mlon = moon.geo_eclip_lon
    mlat = moon.geo_eclip_lat
    dist_km = moon.distance_au * KM_PER_AU
    diam_deg = 2.0 * math.degrees(math.atan(_MOON_MEAN_RADIUS_KM / math.sqrt(dist_km*dist_km - _MOON_MEAN_RADIUS_KM*_MOON_MEAN_RADIUS_KM)))

    # Inclination angle
    I = math.radians(1.543)

    # Moon's argument of latitude in radians.
    f = math.radians(_NormalizeLongitude(93.2720950 + 483202.0175233*t - 0.0036539*t2 - t3/3526000 + t4/863310000))

    # Moon's ascending node's mean longitude in radians.
    omega = math.radians(_NormalizeLongitude(125.0445479 - 1934.1362891*t + 0.0020754*t2 + t3/467441 - t4/60616000))

    # Sun's mean anomaly.
    m = math.radians(_NormalizeLongitude(357.5291092 + 35999.0502909*t - 0.0001536*t2 + t3/24490000))

    # Moon's mean anomaly.
    mdash = math.radians(_NormalizeLongitude(134.9633964 + 477198.8675055*t + 0.0087414*t2 + t3/69699 - t4/14712000))

    # Moon's mean elongation.
    d = math.radians(_NormalizeLongitude(297.8501921 + 445267.1114034*t - 0.0018819*t2 + t3/545868 - t4/113065000))

    # Eccentricity of the Earth's orbit.
    e = 1.0 - 0.002516*t - 0.0000074*t2

    # Optical librations
    w = mlon - omega
    a = math.atan2(math.sin(w)*math.cos(mlat)*math.cos(I) - math.sin(mlat)*math.sin(I), math.cos(w)*math.cos(mlat))
    ldash = _LongitudeOffset(math.degrees(a - f))
    bdash = math.asin(-math.sin(w)*math.cos(mlat)*math.sin(I) - math.sin(mlat)*math.cos(I))

    # Physical librations
    k1 = math.radians(119.75 + 131.849*t)
    k2 = math.radians(72.56 + 20.186*t)

    rho = (
        -0.02752*math.cos(mdash) +
        -0.02245*math.sin(f) +
        +0.00684*math.cos(mdash - 2*f) +
        -0.00293*math.cos(2*f) +
        -0.00085*math.cos(2*f - 2*d) +
        -0.00054*math.cos(mdash - 2*d) +
        -0.00020*math.sin(mdash + f) +
        -0.00020*math.cos(mdash + 2*f) +
        -0.00020*math.cos(mdash - f) +
        +0.00014*math.cos(mdash + 2*f - 2*d)
    )

    sigma = (
        -0.02816*math.sin(mdash) +
        +0.02244*math.cos(f) +
        -0.00682*math.sin(mdash - 2*f) +
        -0.00279*math.sin(2*f) +
        -0.00083*math.sin(2*f - 2*d) +
        +0.00069*math.sin(mdash - 2*d) +
        +0.00040*math.cos(mdash + f) +
        -0.00025*math.sin(2*mdash) +
        -0.00023*math.sin(mdash + 2*f) +
        +0.00020*math.cos(mdash - f) +
        +0.00019*math.sin(mdash - f) +
        +0.00013*math.sin(mdash + 2*f - 2*d) +
        -0.00010*math.cos(mdash - 3*f)
    )

    tau = (
        +0.02520*e*math.sin(m) +
        +0.00473*math.sin(2*mdash - 2*f) +
        -0.00467*math.sin(mdash) +
        +0.00396*math.sin(k1) +
        +0.00276*math.sin(2*mdash - 2*d) +
        +0.00196*math.sin(omega) +
        -0.00183*math.cos(mdash - f) +
        +0.00115*math.sin(mdash - 2*d) +
        -0.00096*math.sin(mdash - d) +
        +0.00046*math.sin(2*f - 2*d) +
        -0.00039*math.sin(mdash - f) +
        -0.00032*math.sin(mdash - m - d) +
        +0.00027*math.sin(2*mdash - m - 2*d) +
        +0.00023*math.sin(k2) +
        -0.00014*math.sin(2*d) +
        +0.00014*math.cos(2*mdash - 2*f) +
        -0.00012*math.sin(mdash - 2*f) +
        -0.00012*math.sin(2*mdash) +
        +0.00011*math.sin(2*mdash - 2*m - 2*d)
    )

    ldash2 = -tau + (rho*math.cos(a) + sigma*math.sin(a))*math.tan(bdash)
    bdash = math.degrees(bdash)
    bdash2 = sigma*math.cos(a) - rho*math.sin(a)
    return LibrationInfo(bdash + bdash2, ldash + ldash2, math.degrees(mlat), math.degrees(mlon), dist_km, diam_deg)


class AxisInfo:
    """Information about a body's rotation axis at a given time.

    This structure is returned by #RotationAxis to report
    the orientation of a body's rotation axis at a given moment in time.
    The axis is specified by the direction in space that the body's north pole
    points, using angular equatorial coordinates in the J2000 system (EQJ).

    Thus `ra` is the right ascension, and `dec` is the declination, of the
    body's north pole vector at the given moment in time. The north pole
    of a body is defined as the pole that lies on the north side of the
    [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
    regardless of the body's direction of rotation.

    The `spin` field indicates the angular position of a prime meridian
    arbitrarily recommended for the body by the International Astronomical
    Union (IAU).

    The fields `ra`, `dec`, and `spin` correspond to the variables
    α0, δ0, and W, respectively, from
    [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).

    The field `north` is a unit vector pointing in the direction of the body's north pole.
    It is expressed in the equatorial J2000 system (EQJ).

    Attributes
    ----------
    ra : float
        The J2000 right ascension of the body's north pole direction, in sidereal hours.
    dec : float
        The J2000 declination of the body's north pole direction, in degrees.
    spin : float
        Rotation angle of the body's prime meridian, in degrees.
    north : Vector
        A J2000 dimensionless unit vector pointing in the direction of the body's north pole.
    """
    def __init__(self, ra, dec, spin, north):
        self.ra = ra
        self.dec = dec
        self.spin = spin
        self.north = north

    def __repr__(self):
        return 'AxisInfo(ra={}, dec={}, spin={}, north={})'.format(
            self.ra,
            self.dec,
            self.spin,
            repr(self.north)
        )


def _EarthRotationAxis(time):
    # Unlike the other planets, we have a model of precession and nutation
    # for the Earth's axis that provides a north pole vector.
    # So calculate the vector first, then derive the (RA,DEC) angles from the vector.

    # Start with a north pole vector in equator-of-date coordinates: (0,0,1).
    # Convert the vector into J2000 coordinates.
    pos2 = _nutation([0, 0, 1], time, _PrecessDir.Into2000)
    nvec = _precession(pos2, time, _PrecessDir.Into2000)
    north = Vector(nvec[0], nvec[1], nvec[2], time)

    # Derive angular values: right ascension and declination.
    equ = EquatorFromVector(north)

    # Use a modified version of the era() function that does not trim to 0..360 degrees.
    # This expression is also corrected to give the correct angle at the J2000 epoch.
    spin = 190.41375788700253 + (360.9856122880876 * time.ut)

    return AxisInfo(equ.ra, equ.dec, spin, north)


def RotationAxis(body, time):
    """Calculates information about a body's rotation axis at a given time.

    Calculates the orientation of a body's rotation axis, along with
    the rotation angle of its prime meridian, at a given moment in time.

    This function uses formulas standardized by the IAU Working Group
    on Cartographics and Rotational Elements 2015 report, as described
    in the following document:

    https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf

    See #AxisInfo for more detailed information.

    Parameters
    ----------
    body : Body
        One of the following values:
        `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`,
        `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
    time : Time
        The time at which to calculate the body's rotation axis.

    Returns
    -------
    AxisInfo
        The body's north pole direction and angle of its prime meridian.
    """
    d = time.tt
    T = d / 36525.0
    if body == Body.Sun:
        ra = 286.13
        dec = 63.87
        w = 84.176 + (14.1844 * d)
    elif body == Body.Mercury:
        ra = 281.0103 - (0.0328 * T)
        dec = 61.4155 - (0.0049 * T)
        w = (
            329.5988
            + (6.1385108 * d)
            + (0.01067257 * math.sin(math.radians(174.7910857 + 4.092335*d)))
            - (0.00112309 * math.sin(math.radians(349.5821714 + 8.184670*d)))
            - (0.00011040 * math.sin(math.radians(164.3732571 + 12.277005*d)))
            - (0.00002539 * math.sin(math.radians(339.1643429 + 16.369340*d)))
            - (0.00000571 * math.sin(math.radians(153.9554286 + 20.461675*d)))
        )
    elif body == Body.Venus:
        ra = 272.76
        dec = 67.16
        w = 160.20 - (1.4813688 * d)
    elif body == Body.Earth:
        return _EarthRotationAxis(time)
    elif body == Body.Moon:
        # See page 8, Table 2 in:
        # https://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf
        E1  = math.radians(125.045 -  0.0529921*d)
        E2  = math.radians(250.089 -  0.1059842*d)
        E3  = math.radians(260.008 + 13.0120009*d)
        E4  = math.radians(176.625 + 13.3407154*d)
        E5  = math.radians(357.529 +  0.9856003*d)
        E6  = math.radians(311.589 + 26.4057084*d)
        E7  = math.radians(134.963 + 13.0649930*d)
        E8  = math.radians(276.617 +  0.3287146*d)
        E9  = math.radians(34.226  +  1.7484877*d)
        E10 = math.radians(15.134  -  0.1589763*d)
        E11 = math.radians(119.743 +  0.0036096*d)
        E12 = math.radians(239.961 +  0.1643573*d)
        E13 = math.radians(25.053  + 12.9590088*d)

        ra = (
            269.9949 + 0.0031*T
            - 3.8787*math.sin(E1)
            - 0.1204*math.sin(E2)
            + 0.0700*math.sin(E3)
            - 0.0172*math.sin(E4)
            + 0.0072*math.sin(E6)
            - 0.0052*math.sin(E10)
            + 0.0043*math.sin(E13)
        )

        dec = (
            66.5392 + 0.0130*T
            + 1.5419*math.cos(E1)
            + 0.0239*math.cos(E2)
            - 0.0278*math.cos(E3)
            + 0.0068*math.cos(E4)
            - 0.0029*math.cos(E6)
            + 0.0009*math.cos(E7)
            + 0.0008*math.cos(E10)
            - 0.0009*math.cos(E13)
        )

        w = (
            38.3213 + (13.17635815 - 1.4e-12*d)*d
            + 3.5610*math.sin(E1)
            + 0.1208*math.sin(E2)
            - 0.0642*math.sin(E3)
            + 0.0158*math.sin(E4)
            + 0.0252*math.sin(E5)
            - 0.0066*math.sin(E6)
            - 0.0047*math.sin(E7)
            - 0.0046*math.sin(E8)
            + 0.0028*math.sin(E9)
            + 0.0052*math.sin(E10)
            + 0.0040*math.sin(E11)
            + 0.0019*math.sin(E12)
            - 0.0044*math.sin(E13)
        )
    elif body == Body.Mars:
        ra = (
            317.269202 - 0.10927547*T
            + 0.000068 * math.sin(math.radians(198.991226 + 19139.4819985*T))
            + 0.000238 * math.sin(math.radians(226.292679 + 38280.8511281*T))
            + 0.000052 * math.sin(math.radians(249.663391 + 57420.7251593*T))
            + 0.000009 * math.sin(math.radians(266.183510 + 76560.6367950*T))
            + 0.419057 * math.sin(math.radians(79.398797 + 0.5042615*T))
        )

        dec = (
            54.432516 - 0.05827105*T
            + 0.000051*math.cos(math.radians(122.433576 + 19139.9407476*T))
            + 0.000141*math.cos(math.radians(43.058401 + 38280.8753272*T))
            + 0.000031*math.cos(math.radians(57.663379 + 57420.7517205*T))
            + 0.000005*math.cos(math.radians(79.476401 + 76560.6495004*T))
            + 1.591274*math.cos(math.radians(166.325722 + 0.5042615*T))
        )

        w = (
            176.049863 + 350.891982443297*d
            + 0.000145*math.sin(math.radians(129.071773 + 19140.0328244*T))
            + 0.000157*math.sin(math.radians(36.352167 + 38281.0473591*T))
            + 0.000040*math.sin(math.radians(56.668646 + 57420.9295360*T))
            + 0.000001*math.sin(math.radians(67.364003 + 76560.2552215*T))
            + 0.000001*math.sin(math.radians(104.792680 + 95700.4387578*T))
            + 0.584542*math.sin(math.radians(95.391654 + 0.5042615*T))
        )

    elif body == Body.Jupiter:
        Ja = math.radians(99.360714 + 4850.4046*T)
        Jb = math.radians(175.895369 + 1191.9605*T)
        Jc = math.radians(300.323162 + 262.5475*T)
        Jd = math.radians(114.012305 + 6070.2476*T)
        Je = math.radians(49.511251 + 64.3000*T)

        ra = (
            268.056595 - 0.006499*T
            + 0.000117*math.sin(Ja)
            + 0.000938*math.sin(Jb)
            + 0.001432*math.sin(Jc)
            + 0.000030*math.sin(Jd)
            + 0.002150*math.sin(Je)
        )

        dec = (
            64.495303 + 0.002413*T
            + 0.000050*math.cos(Ja)
            + 0.000404*math.cos(Jb)
            + 0.000617*math.cos(Jc)
            - 0.000013*math.cos(Jd)
            + 0.000926*math.cos(Je)
        )

        w = 284.95 + 870.536*d

    elif body == Body.Saturn:
        ra = 40.589 - 0.036*T
        dec = 83.537 - 0.004*T
        w = 38.90 + 810.7939024*d

    elif body == Body.Uranus:
        ra = 257.311
        dec = -15.175
        w = 203.81 - 501.1600928*d

    elif body == Body.Neptune:
        N = math.radians(357.85 + 52.316*T)
        ra = 299.36 + 0.70*math.sin(N)
        dec = 43.46 - 0.51*math.cos(N)
        w = 249.978 + 541.1397757*d - 0.48*math.sin(N)

    elif body == Body.Pluto:
        ra = 132.993
        dec = -6.163
        w = 302.695 + 56.3625225*d

    else:
        raise InvalidBodyError()

    # Calculate the north pole vector using the given angles.
    radlat = math.radians(dec)
    radlon = math.radians(ra)
    rcoslat = math.cos(radlat)
    north = Vector(
        rcoslat * math.cos(radlon),
        rcoslat * math.sin(radlon),
        math.sin(radlat),
        time
    )
    return AxisInfo(ra/15.0, dec, w, north)


def LagrangePoint(point, time, major_body, minor_body):
    """Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.

    Given a more massive "major" body and a much less massive "minor" body,
    calculates one of the five Lagrange points in relation to the minor body's
    orbit around the major body. The parameter `point` is an integer that
    selects the Lagrange point as follows:

    1 = the Lagrange point between the major body and minor body.
    2 = the Lagrange point on the far side of the minor body.
    3 = the Lagrange point on the far side of the major body.
    4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
    5 = the Lagrange point 60 degrees behind the minor body's orbital position.

    The function returns the state vector for the selected Lagrange point
    in equatorial J2000 coordinates (EQJ), with respect to the center of the
    major body.

    To calculate Sun/Earth Lagrange points, pass in `Body.Sun` for `major_body`
    and `Body.EMB` (Earth/Moon barycenter) for `minor_body`.
    For Lagrange points of the Sun and any other planet, pass in just that planet
    (e.g. `Body.Jupiter`) for `minor_body`.
    To calculate Earth/Moon Lagrange points, pass in `Body.Earth` and `Body.Moon`
    for the major and minor bodies respectively.

    In some cases, it may be more efficient to call #LagrangePointFast,
    especially when the state vectors have already been calculated, or are needed
    for some other purpose.

    Parameters
    ----------
    point : int
        An integer 1..5 that selects which of the Lagrange points to calculate.
    time : Time
        The time for which the Lagrange point is to be calculated.
    major_body : Body
        The more massive of the co-orbiting bodies: `Body.Sun` or `Body.Earth`.
    minor_body : Body
        The less massive of the co-orbiting bodies. See main remarks.

    Returns
    -------
    StateVector
        The position and velocity of the selected Lagrange point with respect to the major body's center.
    """
    major_mass = MassProduct(major_body)
    minor_mass = MassProduct(minor_body)

    # Calculate the state vectors for the major and minor bodies.
    if major_body == Body.Earth and minor_body == Body.Moon:
        # Use geocentric calculations for more precision.
        major_state = StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time)
        minor_state = GeoMoonState(time)
    else:
        major_state = HelioState(major_body, time)
        minor_state = HelioState(minor_body, time)

    return LagrangePointFast(
        point,
        major_state,
        major_mass,
        minor_state,
        minor_mass
    )


def LagrangePointFast(point, major_state, major_mass, minor_state, minor_mass):
    """Calculates one of the 5 Lagrange points from body masses and state vectors.

    Given a more massive "major" body and a much less massive "minor" body,
    calculates one of the five Lagrange points in relation to the minor body's
    orbit around the major body. The parameter `point` is an integer that
    selects the Lagrange point as follows:

    1 = the Lagrange point between the major body and minor body.
    2 = the Lagrange point on the far side of the minor body.
    3 = the Lagrange point on the far side of the major body.
    4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
    5 = the Lagrange point 60 degrees behind the minor body's orbital position.

    The caller passes in the state vector and mass for both bodies.
    The state vectors can be in any orientation and frame of reference.
    The body masses are expressed as GM products, where G = the universal
    gravitation constant and M = the body's mass. Thus the units for
    `major_mass` and `minor_mass` must be au^3/day^2.
    Use #MassProduct to obtain GM values for various solar system bodies.

    The function returns the state vector for the selected Lagrange point
    using the same orientation as the state vector parameters `major_state` and `minor_state`,
    and the position and velocity components are with respect to the major body's center.

    Consider calling #LagrangePoint, instead of this function, for simpler usage in most cases.

    Parameters
    ----------
    point : int
        An integer 1..5 that selects which of the Lagrange points to calculate.
    major_state : StateVector
        The state vector of the major (more massive) of the pair of bodies.
    major_mass : float
        The mass product GM of the major body.
    minor_state : StateVector
        The state vector of the minor (less massive) of the pair of bodies.
    minor_mass : float
        The mass product GM of the minor body.

    Returns
    -------
    StateVector
        The position and velocity of the selected Lagrange point with respect to the major body's center.
    """
    cos_60 = 0.5
    sin_60 = 0.8660254037844386   # sqrt(3) / 2

    if point not in [1, 2, 3, 4, 5]:
        raise Error('Invalid Lagrange point [{}]. Must be an integer 1..5.'.format(point))

    if major_mass <= 0.0:
        raise Error('Major mass must be a positive number.')

    if minor_mass <= 0.0:
        raise Error('Minor mass must be a positive number.')

    # Find the relative position vector <dx, dy, dz>.
    dx = minor_state.x - major_state.x
    dy = minor_state.y - major_state.y
    dz = minor_state.z - major_state.z
    R2 = (dx*dx + dy*dy + dz*dz)

    # R = Total distance between the bodies.
    R = math.sqrt(R2)

    # Find the velocity vector <vx, vy, vz>.
    vx = minor_state.vx - major_state.vx
    vy = minor_state.vy - major_state.vy
    vz = minor_state.vz - major_state.vz

    if point in [4, 5]:
        # For L4 and L5, we need to find points 60 degrees away from the
        # line connecting the two bodies and in the instantaneous orbital plane.
        # Define the instantaneous orbital plane as the unique plane that contains
        # both the relative position vector and the relative velocity vector.

        # Take the cross product of position and velocity to find a normal vector <nx, ny, nz>.
        nx = dy*vz - dz*vy
        ny = dz*vx - dx*vz
        nz = dx*vy - dy*vx

        # Take the cross product normal*position to get a tangential vector <ux, uy, uz>.
        ux = ny*dz - nz*dy
        uy = nz*dx - nx*dz
        uz = nx*dy - ny*dx

        # Convert the tangential direction vector to a unit vector.
        U = math.sqrt(ux*ux + uy*uy + uz*uz)
        ux /= U
        uy /= U
        uz /= U

        # Convert the relative position vector into a unit vector.
        dx /= R
        dy /= R
        dz /= R

        # Now we have two perpendicular unit vectors in the orbital plane: 'd' and 'u'.

        # Create new unit vectors rotated (+/-)60 degrees from the radius/tangent directions.
        if point == 4:
            vert = +sin_60
        else:
            vert = -sin_60

        # Rotated radial vector
        Dx = cos_60*dx + vert*ux
        Dy = cos_60*dy + vert*uy
        Dz = cos_60*dz + vert*uz

        # Rotated tangent vector
        Ux = cos_60*ux - vert*dx
        Uy = cos_60*uy - vert*dy
        Uz = cos_60*uz - vert*dz

        # Calculate L4/L5 positions relative to the major body.
        px = R * Dx
        py = R * Dy
        pz = R * Dz

        # Use dot products to find radial and tangential components of the relative velocity.
        vrad = vx*dx + vy*dy + vz*dz
        vtan = vx*ux + vy*uy + vz*uz

        # Calculate L4/L5 velocities.
        pvx = vrad*Dx + vtan*Ux
        pvy = vrad*Dy + vtan*Uy
        pvz = vrad*Dz + vtan*Uz

        return StateVector(px, py, pz, pvx, pvy, pvz, major_state.t)

    # Handle Langrange points 1, 2, 3.
    # Calculate the distances of each body from their mutual barycenter.
    # r1 = negative distance of major mass from barycenter (e.g. Sun to the left of barycenter)
    # r2 = positive distance of minor mass from barycenter (e.g. Earth to the right of barycenter)
    r1 = -R * (minor_mass / (major_mass + minor_mass))
    r2 = +R * (major_mass / (major_mass + minor_mass))

    # Calculate the square of the angular orbital speed in [rad^2 / day^2].
    omega2 = (major_mass + minor_mass) / (R2*R)

    # Use Newton's Method to numerically solve for the location where
    # outward centrifugal acceleration in the rotating frame of reference
    # is equal to net inward gravitational acceleration.
    # First derive a good initial guess based on approximate analysis.
    if point in [1, 2]:
        scale = (major_mass / (major_mass + minor_mass)) * _cbrt(minor_mass / (3.0 * major_mass))
        numer1 = -major_mass    # The major mass is to the left of L1 and L2
        if point == 1:
            scale = 1.0 - scale
            numer2 = +minor_mass    # The minor mass is to the right of L1.
        else:
            scale = 1.0 + scale
            numer2 = -minor_mass    # The minor mass is to the left of L2.
    else:  # point == 3
        scale = ((7.0/12.0)*minor_mass - major_mass) / (minor_mass + major_mass)
        numer1 = +major_mass    # major mass is to the right of L3.
        numer2 = +minor_mass    # minor mass is to the right of L3.

    # Iterate Newton's Method until it converges.
    x = R*scale - r1
    while True:
        dr1 = x - r1
        dr2 = x - r2
        accel = omega2*x + numer1/(dr1*dr1) + numer2/(dr2*dr2)
        deriv = omega2 - 2*numer1/(dr1*dr1*dr1) - 2*numer2/(dr2*dr2*dr2)
        deltax = accel/deriv
        x -= deltax
        if abs(deltax/R) <= 1.0e-14:
            break
    scale = (x - r1) / R
    return StateVector(scale*dx, scale*dy, scale*dz, scale*vx, scale*vy, scale*vz, major_state.t)

#--------------------------------------------------------------------------------------------------

class GravitySimulator:
    """A simulation of zero or more small bodies moving through the Solar System.

    This class calculates the movement of arbitrary small bodies,
    such as asteroids or comets, that move through the Solar System.
    It does so by calculating the gravitational forces on the bodies
    from the Sun and planets. The user of this class supplies a
    list of initial positions and velocities for the small bodies.
    Then the class can update the positions and velocities over small
    time steps.
    """

    def __init__(self, originBody, time, bodyStates):
        """Creates a gravity simulation object.

        Parameters
        ----------
        originBody: Body
            Specifies the origin of the reference frame.
            All position vectors and velocity vectors will use `originBody`
            as the origin of the coordinate system.
            This origin applies to all the input vectors provided in the
            `bodyStates` parameter of this function, along with all
            output vectors returned by #GravitySimulator.Update.
            Most callers will want to provide one of the following:
            `Body.Sun` for heliocentric coordinates,
            `Body.SSB` for solar system barycentric coordinates,
            or `Body.Earth` for geocentric coordinates. Note that the
            gravity simulator does not correct for light travel time;
            all state vectors are tied to a Newtonian "instantaneous" time.
        time: Time
            The initial time at which to start the simulation.
        bodyStates: StateVector[]
            An array of zero or more initial state vectors (positions and velocities)
            of the small bodies to be simulated.
            The caller must know the positions and velocities of the small bodies at an initial moment in time.
            Their positions and velocities are expressed with respect to `originBody`,
            using equatorial J2000 orientation (EQJ).
            Positions are expressed in astronomical units (AU).
            Velocities are expressed in AU/day.
            All the times embedded within the state vectors must exactly match `time`,
            or this constructor will throw an exception.
        """
        self._originBody = originBody
        # Verify that the state vectors have matching times.
        for b in bodyStates:
            if b.t.tt != time.tt:
                raise Error('Inconsistent times in bodyStates')

        # Create a stub list of small body states that we will append to later.
        # We just need the stub to put into `self.curr`
        smallBodyList = []

        # Calculate the states of the Sun and planets at the initial time.
        largeBodyDict = _CalcSolarSystem(time)

        # Create a simulation endpoint for the initial time.
        self.curr = _GravSimEndpoint(time, largeBodyDict, smallBodyList)

        # Convert origin-centric bodyStates vectors into a barycentric `_body_grav_calc_t` array.
        o = self._InternalBodyState(originBody)
        for b in bodyStates:
            r = _TerseVector(b.x + o.r.x, b.y + o.r.y, b.z + o.r.z)
            v = _TerseVector(b.vx + o.v.x, b.vy + o.v.y, b.vz + o.v.z)
            a = _TerseVector.zero()
            smallBodyList.append(_body_grav_calc_t(time.tt, r, v, a))

        # Calculate the net acceleration experienced by the small bodies.
        self._CalcBodyAccelerations()

        # To prepare for a possible swap operation, duplicate the current state into the previous state.
        self.prev = self._Duplicate()

    def Time(self):
        """The time represented by the current step of the gravity simulation.

        Returns
        -------
        Time
        """
        return self.curr.time

    def OriginBody(self):
        """The origin of the reference frame. See constructor for more info.

        Returns
        -------
        Body
        """
        return self._originBody

    def Update(self, time):
        """Advances the gravity simulation by a small time step.

        Updates the simulation of the user-supplied small bodies
        to the time indicated by the `time` parameter.
        Returns an array of state vectors for the simulated bodies.
        The array is in the same order as the original array that
        was used to construct this simulator object.
        The positions and velocities in the returned array are
        referenced to the `originBody` that was used to construct
        this simulator.

        Parameters
        ----------
        time : Time
            A time that is a small increment away from the current simulation time.
            It is up to the developer to figure out an appropriate time increment.
            Depending on the trajectories, a smaller or larger increment
            may be needed for the desired accuracy. Some experimentation may be needed.
            Generally, bodies that stay in the outer Solar System and move slowly can
            use larger time steps. Bodies that pass into the inner Solar System and
            move faster will need a smaller time step to maintain accuracy.
            The `time` value may be after or before the current simulation time
            to move forward or backward in time.

        Returns
        -------
        StateVector[]
            An array of state vectors, one for each small body.
        """
        dt = time.tt - self.curr.time.tt
        if dt == 0.0:
            # Special case: the time has not changed, so skip the usual physics calculations.
            # This allows another way for the caller to query the current body states.
            # It is also necessary to avoid dividing by `dt` if `dt` is zero.
            # To prepare for a possible swap operation, duplicate the current state into the previous state.
            self.prev = self._Duplicate()
        else:
            # Exchange the current state with the previous state. Then calculate the new current state.
            self.Swap()

            # Update the current time
            self.curr.time = time

            # Calculate the positions and velocities of the Sun and planets at the given time.
            self.curr.gravitators = _CalcSolarSystem(time)

            # Estimate the positions of the small bodies as if their existing
            # existing accelerations apply across the whole time interval.
            nbodies = len(self.curr.bodies)
            for i in range(nbodies):
                p = self.prev.bodies[i]
                c = self.curr.bodies[i]
                c.r = _UpdatePosition(dt, p.r, p.v, p.a)

            # Calculate the acceleration experienced by the small bodies at
            # their respective approximate next locations.
            self._CalcBodyAccelerations()

            for i in range(nbodies):
                # Calculate the average of the acceleration vectors
                # experienced by the previous body positions and
                # their estimated next positions.
                # These become estimates of the mean effective accelerations
                # over the whole interval.
                p = self.prev.bodies[i]
                c = self.curr.bodies[i]
                acc = p.a.mean(c.a)
                # Refine the estimates of position and velocity at the next time step,
                # using the mean acceleration as a better approximation of the
                # continuously changing acceleration acting on each body.
                c.tt = time.tt
                c.r = _UpdatePosition(dt, p.r, p.v, acc)
                c.v = _UpdateVelocity(dt, p.v, acc)

            # Re-calculate accelerations experienced by each body.
            # These will be needed for the next simulation step (if any).
            # Also, they will be potentially useful if some day we add
            # a function to query the acceleration vectors for the bodies.
            self._CalcBodyAccelerations()

        # Translate our internal calculations of body positions and velocities
        # into state vectors that the caller can understand.
        # We have to convert the internal type _body_grav_calc_t to the public type StateVector.
        # Also convert from barycentric coordinates to coordinates based on the selected origin body.
        bodyStates = []
        ostate = self._InternalBodyState(self._originBody)
        for bcalc in self.curr.bodies:
            bodyStates.append(StateVector(
                bcalc.r.x - ostate.r.x,
                bcalc.r.y - ostate.r.y,
                bcalc.r.z - ostate.r.z,
                bcalc.v.x - ostate.v.x,
                bcalc.v.y - ostate.v.y,
                bcalc.v.z - ostate.v.z,
                time
            ))
        return bodyStates


    def Swap(self):
        """Exchange the current time step with the previous time step.

        Sometimes it is helpful to "explore" various times near a given
        simulation time step, while repeatedly returning to the original
        time step. For example, when backdating a position for light travel
        time, the caller may wish to repeatedly try different amounts of
        backdating. When the backdating solver has converged, the caller
        wants to leave the simulation in its original state.

        This function allows a single "undo" of a simulation, and does so
        very efficiently.

        Usually this function will be called immediately after a matching
        call to #GravitySimulator.Update. It has the effect of rolling
        back the most recent update. If called twice in a row, it reverts
        the swap and thus has no net effect.

        The constructor initializes the current state and previous
        state to be identical. Both states represent the `time` parameter that was
        passed into the constructor. Therefore, `Swap` will
        have no effect from the caller's point of view when passed a simulator
        that has not yet been updated by a call to #GravitySimulator.Update.
        """
        (self.curr, self.prev) = (self.prev, self.curr)

    def SolarSystemBodyState(self, body):
        """Get the position and velocity of a Solar System body included in the simulation.

        In order to simulate the movement of small bodies through the Solar System,
        the simulator needs to calculate the state vectors for the Sun and planets.

        If an application wants to know the positions of one or more of the planets
        in addition to the small bodies, this function provides a way to obtain
        their state vectors. This is provided for the sake of efficiency, to avoid
        redundant calculations.

        The state vector is returned relative to the position and velocity
        of the `originBody` parameter that was passed to this object's constructor.

        Parameters
        ----------
        body : Body
            The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune.

        Returns
        -------
        StateVector
            The state vector of the requested Solar System body.
        """
        bstate = self._InternalBodyState(body)
        ostate = self._InternalBodyState(self._originBody)
        return _ExportState(bstate - ostate, self.curr.time)

    def _InternalBodyState(self, body):
        if body == Body.SSB:
            return _body_state_t(self.curr.time.tt, _TerseVector.zero(), _TerseVector.zero())
        if body in self.curr.gravitators:
            return self.curr.gravitators[body]
        raise Error('Invalid body: {}'.format(body))

    def _CalcBodyAccelerations(self):
        for b in self.curr.bodies:
            b.a = _TerseVector.zero()
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Sun    ].r, _SUN_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Mercury].r, _MERCURY_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Venus  ].r, _VENUS_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Earth  ].r, _EARTH_GM + _MOON_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Mars   ].r, _MARS_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Jupiter].r, _JUPITER_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Saturn ].r, _SATURN_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Uranus ].r, _URANUS_GM)
            _AddAcceleration(b.a, b.r, self.curr.gravitators[Body.Neptune].r, _NEPTUNE_GM)

    def _Duplicate(self):
        # Copy the current stateinto the previous state, so that both become the same moment in time.
        gravitators = {}
        for body, grav in self.curr.gravitators.items():
            gravitators[body] = grav.clone()

        bodies = []
        for b in self.curr.bodies:
            bodies.append(b.clone())

        return _GravSimEndpoint(self.curr.time, gravitators, bodies)

class _GravSimEndpoint:
    def __init__(self, time, gravitators, bodies):
        self.time = time
        self.gravitators = gravitators
        self.bodies = bodies

def _CalcSolarSystem(time):
    d = {}
    # Start with the SSB at zero position and velocity.
    ssb = _body_state_t(time.tt, _TerseVector.zero(), _TerseVector.zero())

    # Calculate the heliocentric position of each planet, and adjust the SSB
    # based on each planet's pull on the Sun.
    d[Body.Mercury] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Mercury, _MERCURY_GM)
    d[Body.Venus  ] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Venus,   _VENUS_GM)
    d[Body.Earth  ] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Earth,   _EARTH_GM + _MOON_GM)
    d[Body.Mars   ] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Mars,    _MARS_GM)
    d[Body.Jupiter] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Jupiter, _JUPITER_GM)
    d[Body.Saturn ] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Saturn,  _SATURN_GM)
    d[Body.Uranus ] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Uranus,  _URANUS_GM)
    d[Body.Neptune] = _AdjustBarycenterPosVel(ssb, time.tt, Body.Neptune, _NEPTUNE_GM)

    # Convert planet states from heliocentric to barycentric.
    for b in d.values():
        b.r -= ssb.r
        b.v -= ssb.v

    # Convert heliocentric SSB to barycentric Sun.
    d[Body.Sun] = _body_state_t(time.tt, -1.0 * ssb.r, -1.0 * ssb.v)
    return d

def _AddAcceleration(acc, smallPos, majorPos, gm):
    dx = majorPos.x - smallPos.x
    dy = majorPos.y - smallPos.y
    dz = majorPos.z - smallPos.z
    r2 = dx*dx + dy*dy + dz*dz
    pull = gm / (r2 * math.sqrt(r2))
    acc.x += dx * pull
    acc.y += dy * pull
    acc.z += dz * pull

#--------------------------------------------------------------------------------------------------
