#!/usr/bin/env python3
#
#    MIT License
#
#    Copyright (c) 2019-2020 Don Cross <cosinekitty@gmail.com>
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

_CalcMoonCount = 0

_DAYS_PER_TROPICAL_YEAR = 365.24217
_PI2 = 2.0 * math.pi
_EPOCH = datetime.datetime(2000, 1, 1, 12)
_ASEC360 = 1296000.0
_ASEC2RAD = 4.848136811095359935899141e-6
_ARC = 3600.0 * 180.0 / math.pi     # arcseconds per radian
_C_AUDAY = 173.1446326846693        # speed of light in AU/day
_KM_PER_AU = 1.4959787069098932e+8
_METERS_PER_AU = _KM_PER_AU * 1000.0
_ANGVEL = 7.2921150e-5
_SECONDS_PER_DAY = 24.0 * 3600.0
_SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592
_MEAN_SYNODIC_MONTH = 29.530588
_EARTH_ORBITAL_PERIOD = 365.256
_NEPTUNE_ORBITAL_PERIOD = 60189.0
_REFRACTION_NEAR_HORIZON = 34.0 / 60.0

_SUN_RADIUS_KM = 695700.0
_SUN_RADIUS_AU  = _SUN_RADIUS_KM / _KM_PER_AU

_EARTH_FLATTENING = 0.996647180302104
_EARTH_EQUATORIAL_RADIUS_KM = 6378.1366
_EARTH_EQUATORIAL_RADIUS_AU = _EARTH_EQUATORIAL_RADIUS_KM / _KM_PER_AU
_EARTH_MEAN_RADIUS_KM = 6371.0      # mean radius of the Earth's geoid, without atmosphere
_EARTH_ATMOSPHERE_KM = 88.0         # effective atmosphere thickness for lunar eclipses
_EARTH_ECLIPSE_RADIUS_KM = _EARTH_MEAN_RADIUS_KM + _EARTH_ATMOSPHERE_KM

_MOON_EQUATORIAL_RADIUS_KM = 1738.1
_MOON_MEAN_RADIUS_KM       = 1737.4
_MOON_POLAR_RADIUS_KM      = 1736.0
_MOON_EQUATORIAL_RADIUS_AU = (_MOON_EQUATORIAL_RADIUS_KM / _KM_PER_AU)

_ASEC180 = 180.0 * 60.0 * 60.0
_AU_PER_PARSEC = _ASEC180 / math.pi
_EARTH_MOON_MASS_RATIO = 81.30056
_SUN_MASS     = 333054.25318      # Sun's mass relative to Earth.
_JUPITER_MASS =    317.84997      # Jupiter's mass relative to Earth.
_SATURN_MASS  =     95.16745      # Saturn's mass relative to Earth.
_URANUS_MASS  =     14.53617      # Uranus's mass relative to Earth.
_NEPTUNE_MASS =     17.14886      # Neptune's mass relative to Earth.


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
        return 'Vector({}, {}, {}, {})'.format(self.x, self.y, self.z, str(self.t))

    def __str__(self):
        return '({}, {}, {}, {})'.format(self.x, self.y, self.z, str(self.t))

    def Length(self):
        """Returns the length of the vector in AU."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

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

class BadTimeError(Error):
    """Cannot calculate Pluto position for this date/time."""
    def __init__(self, time):
        Error.__init__(self, 'Cannot calculate Pluto position for time: {}'.format(time))

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
        Error.__init__(self, 'Internal error - please report issue at https://github.com/cosinekitty/astronomy/issues')

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

def _AngleBetween(a, b):
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

_TimeRegex = re.compile(r'^([0-9]{1,4})-([0-9]{2})-([0-9]{2})(T([0-9]{2}):([0-9]{2})(:([0-9]{2}(\.[0-9]+)?))?Z)?$')

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
    def __init__(self, ut):
        self.ut = ut
        self.tt = _TerrestrialTime(ut)
        self.etilt = None

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
            The UTC 4-digit year value, e.g. 2019.
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
        micro = round(math.fmod(second, 1.0) * 1000000)
        second = math.floor(second - micro/1000000)
        d = datetime.datetime(year, month, day, hour, minute, second, micro)
        ut = (d - _EPOCH).total_seconds() / 86400
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
        and the Terrestrial Time field `tt` is adjusted correctly for the resulting UTC date and time,
        according to the historical and predictive Delta-T model provided by the
        [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).

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
        return 'Time(' + str(self) + ')'

    def __str__(self):
        millis = round(self.ut * 86400000.0)
        n = _EPOCH + datetime.timedelta(milliseconds=millis)
        return '{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}Z'.format(n.year, n.month, n.day, n.hour, n.minute, n.second, math.floor(n.microsecond / 1000))

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
        if self.etilt is None:
            self.etilt = _e_tilt(self)
        return self.etilt


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
        return 'Observer({}, {}, {})'.format(self.latitude, self.longitude, self.height)

    def __str__(self):
        text = '('
        text += 'S' if (self.latitude < 0) else 'N'
        text += '{}, '.format(abs(self.latitude))
        text += 'W' if (self.longitude < 0) else 'E'
        text += '{}, '.format(abs(self.longitude))
        text += '{}m'.format(self.height)
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

def _precession_rot(tt1, tt2):
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
        return RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ])

    # Perform rotation from J2000.0 to other epoch.
    return RotationMatrix([
        [xx, xy, xz],
        [yx, yy, yz],
        [zx, zy, zz]
    ])

def _precession(tt1, pos1, tt2):
    r = _precession_rot(tt1, tt2)
    return [
        r.rot[0][0]*pos1[0] + r.rot[1][0]*pos1[1] + r.rot[2][0]*pos1[2],
        r.rot[0][1]*pos1[0] + r.rot[1][1]*pos1[1] + r.rot[2][1]*pos1[2],
        r.rot[0][2]*pos1[0] + r.rot[1][2]*pos1[1] + r.rot[2][2]*pos1[2]
    ]

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
    """
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
            raise BadVectorError()
        ra = 0.0
        if pos[2] < 0.0:
            dec = -90.0
        else:
            dec = +90.0
    else:
        ra = math.degrees(math.atan2(pos[1], pos[0])) / 15.0
        if ra < 0:
            ra += 24
        dec = math.degrees(math.atan2(pos[2], math.sqrt(xyproj)))
    return Equatorial(ra, dec, dist)


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

    if direction == 0:
        # forward rotation
        return RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ])

    # inverse rotation
    return RotationMatrix([
        [xx, yx, zx],
        [xy, yy, zy],
        [xz, yz, zz]
    ])

def _nutation(time, direction, pos):
    r = _nutation_rot(time, direction)
    return [
        r.rot[0][0]*pos[0] + r.rot[1][0]*pos[1] + r.rot[2][0]*pos[2],
        r.rot[0][1]*pos[0] + r.rot[1][1]*pos[1] + r.rot[2][1]*pos[2],
        r.rot[0][2]*pos[0] + r.rot[1][2]*pos[1] + r.rot[2][2]*pos[2]
    ]

def _era(time):        # Earth Rotation Angle
    thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut
    thet3 = math.fmod(time.ut, 1.0)
    theta = 360.0 * math.fmod((thet1 + thet3), 1.0)
    if theta < 0.0:
        theta += 360.0
    return theta


def _sidereal_time(time):
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
    return gst


def _terra(observer, st):
    df = 1.0 - 0.003352819697896    # flattening of the Earth
    df2 = df * df
    phi = math.radians(observer.latitude)
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    c = 1.0 / math.sqrt(cosphi*cosphi + df2*sinphi*sinphi)
    s = df2 * c
    ht_km = observer.height / 1000.0
    ach = _EARTH_EQUATORIAL_RADIUS_KM*c + ht_km
    ash = _EARTH_EQUATORIAL_RADIUS_KM*s + ht_km
    stlocl = math.radians(15.0*st + observer.longitude)
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
    """Calculates the geocentric position of the Moon at a given time.

    Given a time of observation, calculates the Moon's position as a vector.
    The vector gives the location of the Moon's center relative to the Earth's center
    with x-, y-, and z-components measured in astronomical units.

    This algorithm is based on Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
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
        The Moon's position as a vector in J2000 Cartesian equatorial coordinates.

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
    mpos2 = _precession(time.tt, mpos1, 0)
    return Vector(mpos2[0], mpos2[1], mpos2[2], time)

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

def _VsopFormula(formula, t):
    tpower = 1.0
    coord = 0.0
    for series in formula:
        coord += tpower * sum(A * math.cos(B + C*t) for (A, B, C) in series)
        tpower *= t
    return coord

def _CalcVsop(model, time):
    t = time.tt / 365250.0
    spher = [_VsopFormula(formula, t) for formula in model]

    # Convert spherical coordinates to ecliptic cartesian coordinates.
    r_coslat = spher[2] * math.cos(spher[1])
    ex = r_coslat * math.cos(spher[0])
    ey = r_coslat * math.sin(spher[0])
    ez = spher[2] * math.sin(spher[1])

    # Convert ecliptic cartesian coordinates to equatorial cartesian coordinates.
    vx = ex + 0.000000440360*ey - 0.000000190919*ez
    vy = -0.000000479966*ex + 0.917482137087*ey - 0.397776982902*ez
    vz = 0.397776982902*ey + 0.917482137087*ez
    return Vector(vx, vy, vz, time)

def _VsopHelioDistance(model, time):
    # The caller only wants to know the distance between the planet and the Sun.
    # So we only need to calculate the radial component of the spherical coordinates.
    # There is no need to translate coordinates.
    return _VsopFormula(model[2], time.tt / 365250.0)

def _CalcEarth(time):
    return _CalcVsop(_vsop[Body.Earth.value], time)

# END VSOP
#----------------------------------------------------------------------------
# BEGIN TOP2013

_pluto = [
    [  # f=0
        [  # f=0, s= 0
[        0,  3.9544617144029999e+01,  0.0000000000000000e+00],   # f=0, s= 0, t=   0
[     1402, -1.8891373533434089e-01, -8.5258197635470073e-02],   # f=0, s= 0, t=   1
[     1331, -4.1495877833812339e-02, -3.3387415274886263e-02],   # f=0, s= 0, t=   2
[      522, -4.8502474919249819e-02, -7.3455412272547278e-03],   # f=0, s= 0, t=   3
[       71,  2.8102132918948221e-02, -5.3468437660152152e-03],   # f=0, s= 0, t=   4
[     1261, -9.1600251304608978e-03, -1.2309204554431390e-02],   # f=0, s= 0, t=   5
[      452, -1.2202161344831991e-02, -5.2519856329982890e-03],   # f=0, s= 0, t=   6
[     2875, -9.7845229475467185e-03, -7.2968560644718955e-04],   # f=0, s= 0, t=   7
[       35,  4.8494518209585983e-03, -6.8918970374425084e-03],   # f=0, s= 0, t=   8
[      141,  5.8271375488234932e-03, -3.1946778653436391e-03],   # f=0, s= 0, t=   9
[      137,  1.5300059509576150e-03, -6.0327729791525954e-03],   # f=0, s= 0, t=  10
[        4, -1.9494897408412360e-03,  4.5717130739994704e-03],   # f=0, s= 0, t=  11
[     1190, -1.7524220672664240e-03, -4.2980683502911454e-03],   # f=0, s= 0, t=  12
[      381, -3.1062775803702681e-03, -2.4728667551542258e-03],   # f=0, s= 0, t=  13
[        8, -1.7188663433411050e-03,  2.9756270077158122e-03],   # f=0, s= 0, t=  14
[     1543, -7.7472653128184826e-04,  2.6514626782777680e-03],   # f=0, s= 0, t=  15
[     1115, -1.5405722125111840e-03, -2.0778390548994150e-03],   # f=0, s= 0, t=  16
[     2804, -2.0048397209869230e-03, -4.1957951179189120e-04],   # f=0, s= 0, t=  17
[       67,  9.6850762192148931e-04, -1.5811913714969829e-03],   # f=0, s= 0, t=  18
[      212,  1.5466715821083480e-03, -9.6836654994834209e-04],   # f=0, s= 0, t=  19
[     1119, -1.8820367463891121e-04, -1.4293834479379090e-03],   # f=0, s= 0, t=  20
[    17405, -1.0738845199599739e-03,  9.5985010997943349e-04],   # f=0, s= 0, t=  21
[    28337,  7.5821211083786067e-04,  1.1416213940389449e-03],   # f=0, s= 0, t=  22
[      310, -6.9650370158153983e-04, -1.0024667762364200e-03],   # f=0, s= 0, t=  23
[     1044, -3.2801771454650589e-04, -6.4947140155397116e-04],   # f=0, s= 0, t=  24
[       63,  4.0424767075158291e-04,  5.7315886355325109e-04],   # f=0, s= 0, t=  25
[     1614, -1.3238448498587051e-04,  6.7949492229369074e-04],   # f=0, s= 0, t=  26
[       12, -6.5048480404874143e-04, -1.0430021074697129e-04],   # f=0, s= 0, t=  27
[      283,  3.3929768009878552e-04, -5.1573678840263412e-04],   # f=0, s= 0, t=  28
[      133,  4.0138197254613279e-04,  4.5284627712770291e-04],   # f=0, s= 0, t=  29
[     1421,  5.0823717785117468e-04,  3.1010622577270548e-04],   # f=0, s= 0, t=  30
[     1383, -5.6813906164891585e-04, -1.7225147327178090e-04],   # f=0, s= 0, t=  31
[     4348, -5.1094469101100064e-04,  1.4416513132178369e-04],   # f=0, s= 0, t=  32
[     2733, -4.8192307867155672e-04, -1.8631709444892481e-04],   # f=0, s= 0, t=  33
[      664,  6.3695948832436563e-05,  5.0429756322537705e-04],   # f=0, s= 0, t=  34
[      177, -3.0532744995821309e-04,  4.0226535349675440e-04],   # f=0, s= 0, t=  35
[      204,  3.4677834246970749e-04,  3.2588064363496143e-04],   # f=0, s= 0, t=  36
[     1048,  5.6219036569383140e-05, -4.5145044715373130e-04],   # f=0, s= 0, t=  37
[     1406,  3.6764565360298558e-04, -2.4793326161876619e-04],   # f=0, s= 0, t=  38
[     1398, -6.1399335668996843e-05, -4.0368084856225791e-04],   # f=0, s= 0, t=  39
[      880,  3.8206123841964649e-04,  1.0727867737651879e-04],   # f=0, s= 0, t=  40
[      239, -9.2482984334176681e-05, -3.6357760642322443e-04],   # f=0, s= 0, t=  41
[      541, -3.5492935880988093e-04, -7.2696664262252120e-05],   # f=0, s= 0, t=  42
[    17334, -3.1608420086033548e-04,  1.6282567766535830e-04],   # f=0, s= 0, t=  43
[      503,  3.4204624270103872e-04,  6.9800261361389600e-06],   # f=0, s= 0, t=  44
[    28266,  1.1026801272880990e-04,  3.1928911807394130e-04],   # f=0, s= 0, t=  45
[      974, -1.5403436032839399e-04, -2.6908883645786191e-04],   # f=0, s= 0, t=  46
[      345, -2.1738414747494940e-04,  1.5563797862269319e-04],   # f=0, s= 0, t=  47
[     1924,  2.1398279518720060e-04,  1.4058024849360011e-04],   # f=0, s= 0, t=  48
[      275,  1.7682640345886169e-04,  1.7533174872185701e-04],   # f=0, s= 0, t=  49
[      106,  2.0541555697206479e-04, -1.2965626091678131e-04],   # f=0, s= 0, t=  50
[     1335,  2.1493444414672229e-04, -1.0499061765897920e-04],   # f=0, s= 0, t=  51
[      271,  2.1059070752496800e-04, -9.6647526271973576e-05]    # f=0, s= 0, t=  52
        ],
        [  # f=0, s= 1
[     1402, -2.4541007710854022e-02,  5.3822529675651168e-02],   # f=0, s= 1, t=   0
[        0,  3.7890000000000000e-02,  0.0000000000000000e+00],   # f=0, s= 1, t=   1
[     1331, -1.5981733929364941e-02,  1.9623571765255570e-02],   # f=0, s= 1, t=   2
[      522, -2.1667846285401051e-03,  1.3851523320136261e-02],   # f=0, s= 1, t=   3
[       71,  8.6723020043621556e-04,  4.7241967094509424e-03],   # f=0, s= 1, t=   4
[     1261, -3.8008804453790821e-03,  2.7392780267271222e-03],   # f=0, s= 1, t=   5
[     2875, -5.9015441259769620e-04,  3.3190792202711962e-03],   # f=0, s= 1, t=   6
[       35, -2.0646844712970550e-03, -2.6221746776036530e-03],   # f=0, s= 1, t=   7
[        4, -2.7341177070724500e-03,  8.4711655196688536e-04],   # f=0, s= 1, t=   8
[     1190, -2.1427361371102339e-03,  8.3829844178717706e-04],   # f=0, s= 1, t=   9
[      452, -6.5191457442425422e-04,  1.3962477909825410e-03],   # f=0, s= 1, t=  10
[        8, -1.4174353266390050e-03, -4.8800609526062797e-04],   # f=0, s= 1, t=  11
[      381, -7.7874017402358557e-04,  9.3105764052437356e-04],   # f=0, s= 1, t=  12
[      137, -1.1495119655719700e-03, -2.9568842990989599e-04],   # f=0, s= 1, t=  13
[     2804, -2.9320339800219448e-04,  1.0245691529492771e-03],   # f=0, s= 1, t=  14
[     1119, -9.8344008047702861e-04,  1.2122440012927610e-04],   # f=0, s= 1, t=  15
[     1543,  7.2662531437559094e-04,  2.2380632879329121e-04],   # f=0, s= 1, t=  16
[     1115, -4.9059020449034973e-04,  5.4125269506251448e-04],   # f=0, s= 1, t=  17
[      310, -5.0858314998190228e-04,  3.3071892812221339e-04]    # f=0, s= 1, t=  18
        ],
        [  # f=0, s= 2
[     1402,  6.1564734183986881e-03,  6.9727544398706350e-03],   # f=0, s= 2, t=   0
[     1331,  3.4473718686991932e-03,  5.3547458162748890e-03],   # f=0, s= 2, t=   1
[        0, -6.0199059744408881e-03,  0.0000000000000000e+00],   # f=0, s= 2, t=   2
[      522,  1.8506017571807240e-03,  1.2255527089933200e-03],   # f=0, s= 2, t=   3
[        4, -9.0791196955079877e-04, -1.0597464135497440e-03],   # f=0, s= 2, t=   4
[     1261, -2.5360705302150178e-04,  1.1006253876392330e-03]    # f=0, s= 2, t=   5
        ]
    ],
    [  # f=1
        [  # f=1, s= 0
[        0,  4.1654711248260003e+00,  0.0000000000000000e+00],   # f=1, s= 0, t=   0
[     1402,  2.0610516934972331e-03, -4.5672163937433728e-03],   # f=1, s= 0, t=   1
[        4,  3.4901370659179681e-03,  2.2214931280208532e-03],   # f=1, s= 0, t=   2
[     1473, -2.5274121559055551e-04,  1.4438005411011950e-03],   # f=1, s= 0, t=   3
[        8,  1.1811406701290150e-03,  5.6294240646172070e-04],   # f=1, s= 0, t=   4
[      522,  1.6107717821433110e-04, -1.0623010008063969e-03],   # f=1, s= 0, t=   5
[     1331,  2.7651217171242138e-04, -3.4287860661004818e-04],   # f=1, s= 0, t=   6
[      593,  3.3064342491574279e-05,  3.2281574175741223e-04],   # f=1, s= 0, t=   7
[     2875,  1.7391859631880900e-05, -2.4371046197993790e-04],   # f=1, s= 0, t=   8
[       71,  4.3540142673103847e-05, -2.3324350638665871e-04],   # f=1, s= 0, t=   9
[       12, -1.4733740390542089e-05,  1.6458293257830351e-04],   # f=1, s= 0, t=  10
[       35, -7.7366158840780485e-05, -1.2373453663453069e-04],   # f=1, s= 0, t=  11
[      137,  1.0874981420582371e-04,  2.9642106239502481e-05],   # f=1, s= 0, t=  12
[      452,  3.2837879995135691e-05, -8.0083302669255936e-05],   # f=1, s= 0, t=  13
[      106, -7.1058821662524569e-05, -2.9131217624735831e-05],   # f=1, s= 0, t=  14
[     2945,  1.2419927037655510e-05,  6.9699905790296052e-05],   # f=1, s= 0, t=  15
[     1115,  4.9015184839992363e-05, -3.7626776449287022e-05],   # f=1, s= 0, t=  16
[     1261,  3.7725355301707661e-05, -2.7730601457362580e-05],   # f=1, s= 0, t=  17
[       63,  3.4620953450636643e-05, -2.5727054077819990e-05],   # f=1, s= 0, t=  18
[    17405, -2.4737834776937812e-05, -2.7668964560844859e-05],   # f=1, s= 0, t=  19
[      141, -1.9152627579774910e-05,  3.0608840803230902e-05],   # f=1, s= 0, t=  20
[     1543,  3.4403440209432860e-05,  8.2302490865683481e-06],   # f=1, s= 0, t=  21
[    28337, -2.9450824049753070e-05,  1.9509947544013829e-05],   # f=1, s= 0, t=  22
[      208, -3.3402595486149828e-05,  6.9721024372198099e-07],   # f=1, s= 0, t=  23
[       16, -2.5708226560219951e-05,  5.4690438392817274e-06],   # f=1, s= 0, t=  24
[     2804,  8.2434428446966487e-06, -2.2056820101697141e-05],   # f=1, s= 0, t=  25
[      133,  1.7535496693018641e-05, -1.5317433153649361e-05],   # f=1, s= 0, t=  26
[      177,  9.4209450622295435e-06,  1.8609791865429379e-05],   # f=1, s= 0, t=  27
[      204,  1.3418024861975440e-05, -1.5519809537438971e-05],   # f=1, s= 0, t=  28
[       67,  1.9542040889474939e-05, -2.0438794816549699e-07],   # f=1, s= 0, t=  29
[     1186, -1.1497747150912189e-05,  1.4072246679624251e-05],   # f=1, s= 0, t=  30
[     1421, -7.4843280005381670e-06,  1.2303197921178691e-05],   # f=1, s= 0, t=  31
[     1383,  4.1628765439336240e-06, -1.3719681259667719e-05],   # f=1, s= 0, t=  32
[      275,  7.1803653660818622e-06, -1.2018565522581090e-05],   # f=1, s= 0, t=  33
[      212,  1.3642921909731060e-05,  2.1147480912091179e-06],   # f=1, s= 0, t=  34
[       59, -2.8996742569855378e-06, -1.3093850574478501e-05],   # f=1, s= 0, t=  35
[     4348, -3.6912767653086940e-06, -1.2768209379424250e-05],   # f=1, s= 0, t=  36
[    17476,  8.8159150487736961e-06,  5.9364512059773267e-06],   # f=1, s= 0, t=  37
[       27,  6.2980935945243974e-06, -8.4472797839010716e-06],   # f=1, s= 0, t=  38
[     1406,  5.7734139509693619e-06,  8.7248712048903149e-06],   # f=1, s= 0, t=  39
[     1398,  1.0115989678361140e-05, -1.3220473394572090e-06],   # f=1, s= 0, t=  40
[    28407,  6.7642022362205018e-06, -7.5210649765782616e-06],   # f=1, s= 0, t=  41
[      664,  7.9736063591689162e-06, -2.2024063975673750e-06],   # f=1, s= 0, t=  42
[      381,  3.1761072779613261e-06, -7.3448935507901516e-06],   # f=1, s= 0, t=  43
[      541,  2.5898481389558310e-06, -7.0323742443685127e-06],   # f=1, s= 0, t=  44
[      503, -1.5979599372675240e-07,  7.4586091247535634e-06],   # f=1, s= 0, t=  45
[      271, -3.1728511373578081e-06, -6.6603677156048786e-06],   # f=1, s= 0, t=  46
[      129, -1.8397756711585990e-06, -7.0451043684758007e-06],   # f=1, s= 0, t=  47
[      247,  3.1533300480278400e-07, -7.1604266714610849e-06],   # f=1, s= 0, t=  48
[      200, -2.2707299296152400e-06, -6.6783936815674601e-06],   # f=1, s= 0, t=  49
[       31, -4.2146281840702766e-06, -5.6388188974303981e-06],   # f=1, s= 0, t=  50
[      341, -3.9640903174182786e-06, -5.5968952746387096e-06],   # f=1, s= 0, t=  51
[       20, -2.1086487296606532e-06, -4.2993743513743323e-06],   # f=1, s= 0, t=  52
[      974, -1.9337296464942458e-06,  4.1280092704873532e-06],   # f=1, s= 0, t=  53
[     1492,  1.2386884071711840e-06, -4.0279260453190260e-06],   # f=1, s= 0, t=  54
[     1454, -1.8674475925452479e-07,  4.1994903608977473e-06],   # f=1, s= 0, t=  55
[       55, -4.1054351042334211e-06, -5.3592659852261218e-07],   # f=1, s= 0, t=  56
[     1044,  3.5262751998119361e-06,  1.9367410649084749e-06],   # f=1, s= 0, t=  57
[       39, -7.0347175355875052e-07,  3.9228440423731852e-06],   # f=1, s= 0, t=  58
[      345, -2.7714162211172452e-06, -2.8012028299863999e-06],   # f=1, s= 0, t=  59
[      283, -3.8855399506750791e-06,  4.2353837211089759e-07],   # f=1, s= 0, t=  60
[     4418,  1.9380717555480920e-06,  3.3362906758656841e-06],   # f=1, s= 0, t=  61
[     1708,  3.8363066095252423e-06,  4.0071146882428590e-07],   # f=1, s= 0, t=  62
[      903, -2.4535489999476891e-06,  2.7713257447896069e-06],   # f=1, s= 0, t=  63
[      412, -2.4194398659193562e-06, -2.6818837171571071e-06],   # f=1, s= 0, t=  64
[    17334, -1.5878732753673419e-06, -3.0847459710695832e-06],   # f=1, s= 0, t=  65
[     1410,  2.1528516528967772e-06,  2.5955495586096718e-06],   # f=1, s= 0, t=  66
[      318, -1.6007757166478231e-06,  2.9567402302706461e-06],   # f=1, s= 0, t=  67
[    28266, -3.1282886569147322e-06,  1.0733981547714400e-06],   # f=1, s= 0, t=  68
[     1394,  3.2570224915190650e-06, -5.9986669416035846e-08],   # f=1, s= 0, t=  69
[    72490,  7.8958524750142168e-07,  3.1029527324934821e-06],   # f=1, s= 0, t=  70
[     9220,  2.8325155820290858e-06, -1.4397983817771050e-06],   # f=1, s= 0, t=  71
[     1614,  3.0748937595146630e-06,  5.0326493238056022e-07],   # f=1, s= 0, t=  72
[      408, -3.0722440820616450e-06,  4.9345659005564487e-07],   # f=1, s= 0, t=  73
[      337, -2.9632105019568382e-06,  9.3548023528643972e-08],   # f=1, s= 0, t=  74
[       79,  2.7845485358443480e-06,  7.4025133333557920e-07],   # f=1, s= 0, t=  75
[      354,  1.9360680354843408e-06,  2.0552415375649160e-06],   # f=1, s= 0, t=  76
[      612,  5.3703132984668158e-07,  2.6669817089723940e-06],   # f=1, s= 0, t=  77
[      279, -2.5413275396432281e-07,  2.6502544396045040e-06],   # f=1, s= 0, t=  78
[       75, -2.8830829576295098e-07,  2.6106707165907741e-06],   # f=1, s= 0, t=  79
[      479, -2.4319941706034390e-06,  7.6942886776602959e-07],   # f=1, s= 0, t=  80
[      416,  2.0096007406999570e-06,  1.5520759224768389e-06],   # f=1, s= 0, t=  81
[     1096, -2.3257110689311771e-06,  9.7831153778294669e-07],   # f=1, s= 0, t=  82
[      267, -2.5074175872274202e-06, -1.7701423128546709e-07],   # f=1, s= 0, t=  83
[     3162,  2.0698330022272070e-06, -1.3093627746428469e-06],   # f=1, s= 0, t=  84
[      832, -1.9285373422971391e-06,  1.4160461431375550e-06],   # f=1, s= 0, t=  85
[     1190,  2.2341390259691221e-06, -7.2482790590962755e-07],   # f=1, s= 0, t=  86
[      526,  1.7171371907181510e-06,  1.4529625626913769e-06],   # f=1, s= 0, t=  87
[     2733,  1.0269773239622380e-06, -2.0007533912028411e-06],   # f=1, s= 0, t=  88
[      574, -5.1232176836951890e-07, -2.1792410718247190e-06]    # f=1, s= 0, t=  89
        ],
        [  # f=1, s= 1
[        0,  2.5335660204370001e+01,  0.0000000000000000e+00],   # f=1, s= 1, t=   0
[        4, -2.1897916824442529e-04,  1.7741955864290151e-03],   # f=1, s= 1, t=   1
[     1402, -1.3029527067277161e-03, -5.8987171410210541e-04],   # f=1, s= 1, t=   2
[        8, -5.7857466108113462e-05,  6.2766572590063064e-04],   # f=1, s= 1, t=   3
[      522, -3.0338218770459020e-04, -4.6787236333551149e-05],   # f=1, s= 1, t=   4
[     1331, -1.6234312756262061e-04, -1.3226650376835280e-04],   # f=1, s= 1, t=   5
[       35, -1.6846578574212679e-04,  5.1474337733700072e-05],   # f=1, s= 1, t=   6
[     1473,  1.3716672333829441e-04,  2.7605165443775319e-05],   # f=1, s= 1, t=   7
[       12, -1.1426985570984979e-04,  1.7062621376287921e-05],   # f=1, s= 1, t=   8
[     2875, -8.2654240464229338e-05, -1.4230432684507880e-05],   # f=1, s= 1, t=   9
[     2945,  3.6084370654255798e-05, -3.8286603564628976e-06],   # f=1, s= 1, t=  10
[      593,  3.1050555447062818e-05, -2.3120488185252789e-06],   # f=1, s= 1, t=  11
[       16, -8.7245446550032187e-06, -2.2703299980076952e-05],   # f=1, s= 1, t=  12
[      137,  5.5783261240786444e-06, -2.0642718216052691e-05],   # f=1, s= 1, t=  13
[     1115, -1.3513503169043030e-05, -1.1380785261757990e-05],   # f=1, s= 1, t=  14
[       71,  1.0574660674287450e-05,  1.3281184472910570e-05],   # f=1, s= 1, t=  15
[     1261, -8.2709801392281214e-06, -1.1688074277507450e-05],   # f=1, s= 1, t=  16
[     2804, -1.1505534938221521e-05, -5.1883886955282330e-06],   # f=1, s= 1, t=  17
[       27, -9.0174735475683120e-06, -5.2285448475511598e-06],   # f=1, s= 1, t=  18
[      452, -9.3099342968225354e-06, -3.8372164516008766e-06],   # f=1, s= 1, t=  19
[     1543,  2.4828575108302229e-06, -9.6169483680059550e-06],   # f=1, s= 1, t=  20
[       63, -6.5055787623425069e-06, -7.4385376578198243e-06],   # f=1, s= 1, t=  21
[      106, -4.5512156458035857e-06,  8.6658259433823040e-06],   # f=1, s= 1, t=  22
[      133, -6.3072609842624613e-06, -6.9148116074105151e-06],   # f=1, s= 1, t=  23
[    28337, -4.4530735415054223e-06, -6.7124377273408592e-06],   # f=1, s= 1, t=  24
[      141,  6.9466051844107804e-06,  4.0449257224296616e-06],   # f=1, s= 1, t=  25
[     1398,  1.0522658062293921e-06, -7.0563446922276261e-06],   # f=1, s= 1, t=  26
[       59, -5.7183266993529056e-06,  1.6273101086020201e-06],   # f=1, s= 1, t=  27
[     1383, -5.6807273131313079e-06, -9.5798703870467336e-07],   # f=1, s= 1, t=  28
[       20,  4.8130148431204634e-06, -2.5299767388566989e-06],   # f=1, s= 1, t=  29
[     4348, -5.3269871823774314e-06,  6.0216192993201565e-07],   # f=1, s= 1, t=  30
[      129, -4.2816786763917080e-06,  1.2524207693334740e-06]    # f=1, s= 1, t=  31
        ],
        [  # f=1, s= 2
[        0, -1.8272218839163919e-02,  0.0000000000000000e+00],   # f=1, s= 2, t=   0
[        4, -4.2382205514535369e-04,  6.0951225139627217e-05],   # f=1, s= 2, t=   1
[        8, -2.4214950161520300e-04,  1.3555498866999700e-04],   # f=1, s= 2, t=   2
[     1402, -1.6731953542354990e-04,  1.4877590360992571e-04],   # f=1, s= 2, t=   3
[       12, -4.0514580902466651e-05, -4.6859429953625132e-05],   # f=1, s= 2, t=   4
[     1331, -4.4255479363996823e-05,  2.8523966290637259e-05],   # f=1, s= 2, t=   5
[      522, -2.6626163791203940e-05,  4.0414506964273448e-05],   # f=1, s= 2, t=   6
[       35, -2.3510812229573040e-06,  1.7218835328103482e-05],   # f=1, s= 2, t=   7
[     2875, -7.9375576768093987e-06,  1.3976740463566499e-05],   # f=1, s= 2, t=   8
[       16,  9.7305978044426803e-06, -1.0841117038814680e-05],   # f=1, s= 2, t=   9
[     2945, -4.5741854272088092e-07, -9.4120776625816716e-06]    # f=1, s= 2, t=  10
        ],
        [  # f=1, s= 3
[        0,  1.9409931667071581e-03,  0.0000000000000000e+00],   # f=1, s= 3, t=   0
[        8, -5.2528177259165953e-05, -6.1055439661395280e-05],   # f=1, s= 3, t=   1
[        4, -6.0738316603187738e-05,  3.5387155556321078e-05],   # f=1, s= 3, t=   2
[     1402,  1.5967276563962451e-05,  3.5538891371840628e-05],   # f=1, s= 3, t=   3
[       12,  1.5332065032938759e-05, -2.1295985440213221e-05]    # f=1, s= 3, t=   4
        ],
        [  # f=1, s= 4
[        0,  8.6099959150566781e-05,  0.0000000000000000e+00],   # f=1, s= 4, t=   0
[        4, -5.3544872037241911e-05, -3.2079251781364850e-05]    # f=1, s= 4, t=   1
        ]
    ],
    [  # f=2
        [  # f=2, s= 0
[        0, -1.7873895940349999e-01,  0.0000000000000000e+00],   # f=2, s= 0, t=   0
[     1473,  3.1629832749992988e-03, -2.0985870294942081e-03],   # f=2, s= 0, t=   1
[       71, -6.8180357663860398e-04,  1.1113519163940930e-03],   # f=2, s= 0, t=   2
[     1331,  1.6911781266743710e-04,  1.1885125258969610e-03],   # f=2, s= 0, t=   3
[      593,  5.4774166542017916e-04, -6.3844158559171749e-04],   # f=2, s= 0, t=   4
[     1402, -2.1438473609329679e-05, -5.2026929830140383e-04],   # f=2, s= 0, t=   5
[     1261, -5.7368849879358177e-05,  4.6888445010058260e-04],   # f=2, s= 0, t=   6
[      141, -9.6040157452198532e-05,  3.0415286563684890e-04],   # f=2, s= 0, t=   7
[      452,  1.2124121714193030e-04,  2.7470351829330812e-04],   # f=2, s= 0, t=   8
[     2945,  1.1035693239374779e-04, -1.4867545478165761e-04],   # f=2, s= 0, t=   9
[     1190, -5.9767419503565077e-05,  1.5052450817429230e-04],   # f=2, s= 0, t=  10
[     1543, -1.0081841648710751e-04,  1.1874253775582299e-04],   # f=2, s= 0, t=  11
[      522, -3.7998691236472437e-05, -1.1710923038633279e-04],   # f=2, s= 0, t=  12
[      381,  1.7939534333569889e-05,  1.2165766508553840e-04],   # f=2, s= 0, t=  13
[        8, -5.1065826531599063e-05, -1.0380178618772450e-04],   # f=2, s= 0, t=  14
[        4, -9.6101012457693182e-05, -3.7409967137145340e-05],   # f=2, s= 0, t=  15
[      212, -1.0011374095546651e-05,  9.1861092027193070e-05],   # f=2, s= 0, t=  16
[      208,  6.3118700287276614e-05,  6.6551813665925998e-05],   # f=2, s= 0, t=  17
[      106,  4.6362839297978442e-05,  7.1895265631125496e-05],   # f=2, s= 0, t=  18
[     2804,  2.5951137645879960e-05,  4.8810749391069221e-05],   # f=2, s= 0, t=  19
[     1119, -3.1798344550251363e-05,  4.3442443246298743e-05],   # f=2, s= 0, t=  20
[       35, -2.6365309213148771e-05, -3.7318748312868848e-05],   # f=2, s= 0, t=  21
[     1186,  4.5373467873211630e-05, -3.9315678706168707e-06],   # f=2, s= 0, t=  22
[      310, -5.7440916735318752e-06,  4.2608008327900133e-05],   # f=2, s= 0, t=  23
[       67, -3.2628206231552743e-05,  1.4188485967101389e-05],   # f=2, s= 0, t=  24
[      664, -1.3138986439777871e-05,  2.8745825031870391e-05],   # f=2, s= 0, t=  25
[    17476, -4.5717070878285386e-06, -2.7229104695341480e-05],   # f=2, s= 0, t=  26
[      283,  4.1649905603620409e-06,  2.5778764278582031e-05],   # f=2, s= 0, t=  27
[    28407, -2.6086891392275169e-05,  6.5019643042548361e-07],   # f=2, s= 0, t=  28
[     2875, -1.2604181155312930e-05, -2.2155626126791911e-05],   # f=2, s= 0, t=  29
[     2733,  4.6350251894966598e-06,  2.1642396781861640e-05],   # f=2, s= 0, t=  30
[       12,  1.3787761564064851e-05, -1.6882625493196851e-05],   # f=2, s= 0, t=  31
[       63,  7.3243482301750183e-06, -2.0430415756143061e-05],   # f=2, s= 0, t=  32
[      137,  1.7688914790024399e-05, -2.6466558968384412e-06],   # f=2, s= 0, t=  33
[     1048, -1.3743478634746570e-05,  1.1178051438197709e-05],   # f=2, s= 0, t=  34
[     1614, -1.4527685432017941e-05, -2.5021418914365550e-06],   # f=2, s= 0, t=  35
[     1044, -5.5188522930647246e-06,  1.3429997819426629e-05],   # f=2, s= 0, t=  36
[      239, -6.3224380303283023e-06,  1.2677314304305631e-05],   # f=2, s= 0, t=  37
[      133,  2.9376054111895790e-06, -1.3077911011252820e-05],   # f=2, s= 0, t=  38
[     1492, -9.7624713292315362e-06,  4.8639192065138997e-06],   # f=2, s= 0, t=  39
[     1454,  8.1915641060518482e-06, -7.1484181790528786e-06],   # f=2, s= 0, t=  40
[     4418,  2.9126133839920300e-06, -9.7148097661500440e-06],   # f=2, s= 0, t=  41
[      177, -8.1108828396444682e-06, -5.8479936691072763e-06],   # f=2, s= 0, t=  42
[     2663, -7.4109426380852051e-07,  8.1320485449777430e-06],   # f=2, s= 0, t=  43
[    17334,  7.7062782105104396e-06,  2.1832991677727231e-06],   # f=2, s= 0, t=  44
[     3016, -2.7097054305923721e-06,  7.2985152622098958e-06],   # f=2, s= 0, t=  45
[    28266,  3.0978305562474201e-06, -6.9523478329398256e-06],   # f=2, s= 0, t=  46
[       75, -7.5005963705916821e-06,  5.8552558907898478e-07],   # f=2, s= 0, t=  47
[      354,  2.2767150385710288e-06,  7.1310390742790356e-06],   # f=2, s= 0, t=  48
[      974, -2.8814074492595782e-06,  6.7739175425797296e-06],   # f=2, s= 0, t=  49
[       59, -4.8255490447201409e-06, -5.3633858025303309e-06],   # f=2, s= 0, t=  50
[     1115,  2.4165363833927252e-06, -5.9760260615249186e-06],   # f=2, s= 0, t=  51
[      204, -1.5420460797033710e-07, -6.3645330720578637e-06],   # f=2, s= 0, t=  52
[      612,  4.7791204064686441e-06, -4.1989180399775663e-06],   # f=2, s= 0, t=  53
[      574, -3.2533619740383079e-06,  4.8942036839749464e-06],   # f=2, s= 0, t=  54
[      978, -5.2816161970547299e-06,  2.4194318535884879e-06],   # f=2, s= 0, t=  55
[      129, -4.1704589160563732e-06, -3.8074741128128591e-06],   # f=2, s= 0, t=  56
[      200, -4.3860533690857434e-06, -2.8426817721421870e-06],   # f=2, s= 0, t=  57
[     1685, -4.4937857977759297e-06, -2.6664132359831861e-06],   # f=2, s= 0, t=  58
[     1327, -4.0394363815867892e-06,  3.0479150446235980e-06],   # f=2, s= 0, t=  59
[     1257, -4.4570496522011644e-06,  1.8417244601604851e-06]    # f=2, s= 0, t=  60
        ],
        [  # f=2, s= 1
[        0, -6.1339663802007860e-04,  0.0000000000000000e+00],   # f=2, s= 1, t=   0
[     1331,  5.6664110172313892e-04, -8.1133128701713292e-05],   # f=2, s= 1, t=   1
[     1473, -1.9632805654171590e-04, -2.9843120844069362e-04],   # f=2, s= 1, t=   2
[       71, -1.9721831884265259e-04, -1.2980246221969450e-04],   # f=2, s= 1, t=   3
[     1402, -1.4909525632022481e-04,  5.5470899732068926e-06],   # f=2, s= 1, t=   4
[     1261,  1.4379398981882231e-04,  1.8301247697520619e-05],   # f=2, s= 1, t=   5
[     2945, -7.2143717725975376e-05, -6.1374574108146777e-05],   # f=2, s= 1, t=   6
[     1190,  7.4611359353024315e-05,  3.0289539777815409e-05],   # f=2, s= 1, t=   7
[      593, -5.9939495914062063e-05, -5.1836225465555612e-05],   # f=2, s= 1, t=   8
[        8,  3.9984073707170713e-05, -3.0799116570219068e-05],   # f=2, s= 1, t=   9
[     1543,  3.1094179263405468e-05,  2.6861784078602320e-05],   # f=2, s= 1, t=  10
[      381,  3.7436713283265682e-05, -5.1576484771147526e-06],   # f=2, s= 1, t=  11
[     1119,  2.9759191459225299e-05,  2.2175266912892650e-05],   # f=2, s= 1, t=  12
[      452,  3.2286698454013622e-05, -1.4330028054859710e-05],   # f=2, s= 1, t=  13
[      522, -3.3610639556594409e-05,  1.0732560764842151e-05],   # f=2, s= 1, t=  14
[     2804,  2.6892727962532569e-05, -1.2202088159044841e-05],   # f=2, s= 1, t=  15
[        4,  1.6145362682689730e-05, -2.1535180644982459e-05],   # f=2, s= 1, t=  16
[       35,  1.1619842301292181e-05,  2.1259793782462840e-05],   # f=2, s= 1, t=  17
[      310,  2.1159419772930708e-05,  3.2365387780938209e-06],   # f=2, s= 1, t=  18
[      212, -1.7383881837993001e-05, -1.2570139978872199e-07],   # f=2, s= 1, t=  19
[     2733,  1.5857546421435610e-05, -2.6232455173985150e-06],   # f=2, s= 1, t=  20
[     1048,  9.8372722970562282e-06,  1.2182701032807050e-05],   # f=2, s= 1, t=  21
[       12,  1.1952587407675961e-05,  7.5377841861106607e-06],   # f=2, s= 1, t=  22
[      283, -1.0224581925761621e-05,  1.2303488561812810e-06],   # f=2, s= 1, t=  23
[      239,  8.6293198115122691e-06,  4.6174857696641109e-06]    # f=2, s= 1, t=  24
        ],
        [  # f=2, s= 2
[     1331,  2.3820657956875651e-05, -1.4117217046583029e-04],   # f=2, s= 2, t=   0
[        0,  5.7169784317623151e-05,  0.0000000000000000e+00],   # f=2, s= 2, t=   1
[     1261,  2.8584658843704731e-05, -1.9153125455819128e-05],   # f=2, s= 2, t=   2
[       71, -6.3863669075967202e-06, -3.0182516098252429e-05],   # f=2, s= 2, t=   3
[     2945, -1.6831473476249729e-05,  1.7778057659448761e-05],   # f=2, s= 2, t=   4
[     1402, -8.7282659443619542e-06,  2.1969730180690399e-05],   # f=2, s= 2, t=   5
[     1190,  1.8655664239774939e-05, -1.4275349333837280e-05],   # f=2, s= 2, t=   6
[        8,  2.0572423386833478e-05,  2.9970650695923720e-06]    # f=2, s= 2, t=   7
        ]
    ],
    [  # f=3
        [  # f=3, s= 0
[        0, -1.7340471864230000e-01,  0.0000000000000000e+00],   # f=3, s= 0, t=   0
[     1473,  2.1234813115076170e-03,  3.0130058621335031e-03],   # f=3, s= 0, t=   1
[       71, -1.0505499200359431e-03, -7.0849649152681434e-04],   # f=3, s= 0, t=   2
[     1331,  1.1926139685020671e-03, -1.2467514555563171e-04],   # f=3, s= 0, t=   3
[      593,  6.3328818822552336e-04,  5.1838245948530805e-04],   # f=3, s= 0, t=   4
[     1402, -4.0198902705438299e-04,  3.4879411820834029e-04],   # f=3, s= 0, t=   5
[     1261,  4.6754614871725769e-04,  6.5799333440031267e-05],   # f=3, s= 0, t=   6
[      141, -3.1574545029112658e-04, -1.1057917696717649e-04],   # f=3, s= 0, t=   7
[      452,  2.7842572840314938e-04, -1.1076531434037220e-04],   # f=3, s= 0, t=   8
[     2945,  1.4727652777389929e-04,  1.0316600256389621e-04],   # f=3, s= 0, t=   9
[     1190,  1.4974726446281261e-04,  6.1544751468406387e-05],   # f=3, s= 0, t=  10
[     1543, -1.1418008898386521e-04, -9.8159556344747838e-05],   # f=3, s= 0, t=  11
[      522, -6.9432011255429296e-05,  1.0500648997119170e-04],   # f=3, s= 0, t=  12
[      381,  1.2168346615540400e-04, -1.5724958971846889e-05],   # f=3, s= 0, t=  13
[        4,  3.5694868279606838e-05, -1.1563400718698000e-04],   # f=3, s= 0, t=  14
[        8,  1.0601495232676501e-04, -5.3792770194159239e-05],   # f=3, s= 0, t=  15
[      212, -8.8206829660081605e-05, -5.6186044809908770e-06],   # f=3, s= 0, t=  16
[      208, -6.2116275994912254e-05,  6.2681468828265711e-05],   # f=3, s= 0, t=  17
[      106, -5.9881700950469839e-05,  1.8173312312951261e-05],   # f=3, s= 0, t=  18
[     2804,  4.9734881336053070e-05, -2.4663551556821301e-05],   # f=3, s= 0, t=  19
[     1119,  4.3142176920002318e-05,  3.2166118196149167e-05],   # f=3, s= 0, t=  20
[     1186,  5.3160131153765831e-06,  4.6681999520445390e-05],   # f=3, s= 0, t=  21
[      310,  4.2381917128956993e-05,  6.2114355060802278e-06],   # f=3, s= 0, t=  22
[       67,  9.7801450369382322e-06,  3.8008739374555199e-05],   # f=3, s= 0, t=  23
[       35,  3.4777523130477253e-05,  8.0807437621346330e-06],   # f=3, s= 0, t=  24
[      664, -2.7772391338347379e-05, -1.2933056019872690e-05],   # f=3, s= 0, t=  25
[    17476,  2.6116611382398249e-05, -5.2512707918252350e-06],   # f=3, s= 0, t=  26
[      283, -2.5901894192097461e-05,  7.6613799535167943e-07],   # f=3, s= 0, t=  27
[    28407, -1.3727149833568929e-06, -2.5571899427495732e-05],   # f=3, s= 0, t=  28
[       12,  1.7840323650284740e-05,  1.4017491541647661e-05],   # f=3, s= 0, t=  29
[     2875, -1.2752702589739361e-05,  1.8674367648851481e-05],   # f=3, s= 0, t=  30
[       63, -2.1141833621778831e-05, -6.5305966080253417e-06],   # f=3, s= 0, t=  31
[     2733,  2.1699001266656308e-05, -4.2325482084098446e-06],   # f=3, s= 0, t=  32
[     1048,  1.1071845430997370e-05,  1.3823004095404599e-05],   # f=3, s= 0, t=  33
[      137, -1.3154331920868999e-05, -9.3564722901572934e-06],   # f=3, s= 0, t=  34
[     1044,  1.3602483536936790e-05,  6.6866532885834768e-06],   # f=3, s= 0, t=  35
[     1614,  3.0521160394305390e-06, -1.4295846223173739e-05],   # f=3, s= 0, t=  36
[      239,  1.2491201410693211e-05,  6.3469005652658392e-06],   # f=3, s= 0, t=  37
[      133, -1.3409179117794041e-05, -2.6113349355778870e-06],   # f=3, s= 0, t=  38
[     1492, -4.9926812258289137e-06, -9.3442347365693383e-06],   # f=3, s= 0, t=  39
[     1454,  7.1182846121023137e-06,  7.7835699509488784e-06],   # f=3, s= 0, t=  40
[     4418,  9.5085270250543767e-06,  2.5720788995216492e-06],   # f=3, s= 0, t=  41
[      354, -7.7105865781446408e-06,  3.5195735852526461e-06],   # f=3, s= 0, t=  42
[     2663,  8.0823599962064017e-06,  8.8289336639692517e-07],   # f=3, s= 0, t=  43
[       75, -6.1089854320482681e-07, -8.0284473434343767e-06],   # f=3, s= 0, t=  44
[    17334,  2.4666413292489791e-06, -7.6105265820617346e-06],   # f=3, s= 0, t=  45
[      974,  6.9130746927758178e-06,  3.5048585089334658e-06],   # f=3, s= 0, t=  46
[    28266, -6.8235793370187026e-06, -3.3487790585021339e-06],   # f=3, s= 0, t=  47
[     3016, -7.0427593206860361e-06, -2.6580682576526168e-06],   # f=3, s= 0, t=  48
[       59, -5.5456494164143571e-06,  4.9667283828556126e-06],   # f=3, s= 0, t=  49
[      204, -7.0017722378178446e-06, -1.6817751033546979e-06],   # f=3, s= 0, t=  50
[     1115, -6.4024937184914472e-06,  2.1192315231268140e-06],   # f=3, s= 0, t=  51
[      978,  2.3623585675585379e-06,  5.3260537432133881e-06],   # f=3, s= 0, t=  52
[      129, -3.9867297079366098e-06,  4.2341276287979132e-06],   # f=3, s= 0, t=  53
[      574, -4.8314043929159330e-06, -2.9544972252560498e-06],   # f=3, s= 0, t=  54
[      612,  4.0663883671886490e-06,  3.8452439907704313e-06],   # f=3, s= 0, t=  55
[      200, -3.2782990554685469e-06,  4.3746397303989681e-06],   # f=3, s= 0, t=  56
[      247,  4.4086088444075100e-06, -3.0262449322858620e-06],   # f=3, s= 0, t=  57
[     1685,  2.7665561507113492e-06, -4.4616466628081997e-06],   # f=3, s= 0, t=  58
[     1335, -2.3134382826680641e-06,  4.5013914052228881e-06],   # f=3, s= 0, t=  59
[     1327,  3.1175433103076271e-06,  3.7924860549293981e-06],   # f=3, s= 0, t=  60
[      903,  4.4741994730303496e-06,  1.7822080022821150e-06],   # f=3, s= 0, t=  61
[       16, -1.4927342901310379e-06,  4.4884086825602869e-06],   # f=3, s= 0, t=  62
[      169,  2.5354018835911621e-06,  3.5730000800789471e-06],   # f=3, s= 0, t=  63
[      271, -2.1688629242033941e-06,  3.7820832945120280e-06],   # f=3, s= 0, t=  64
[       79,  3.8979055894631164e-06, -1.8923407502835230e-06],   # f=3, s= 0, t=  65
[      416,  3.9376438761867171e-06, -5.2538703904847715e-07],   # f=3, s= 0, t=  66
[    17405,  1.0904344065136349e-06,  3.5409195078133210e-06],   # f=3, s= 0, t=  67
[    28337,  3.4736718287656798e-06, -5.2407682497120616e-07],   # f=3, s= 0, t=  68
[     1351, -3.4436468577006500e-06, -7.4140986208801818e-08],   # f=3, s= 0, t=  69
[     1312,  3.3433376977429141e-06, -7.9547674738811048e-07],   # f=3, s= 0, t=  70
[      832,  3.0647904819735712e-06,  1.2969722855848811e-06],   # f=3, s= 0, t=  71
[      275,  1.5394331893610239e-06, -2.7751033895625021e-06],   # f=3, s= 0, t=  72
[     2592,  2.8273974520595848e-06,  1.3333677990888271e-06],   # f=3, s= 0, t=  73
[    17264,  1.6047868553221060e-06, -2.5731990436490451e-06]    # f=3, s= 0, t=  74
        ],
        [  # f=3, s= 1
[     1331, -6.0044174932694893e-05, -5.6841814115066798e-04],   # f=3, s= 1, t=   0
[     1473,  2.8299001111328808e-04, -1.9899683117934169e-04],   # f=3, s= 1, t=   1
[       71,  1.3109847892346901e-04, -1.9645547029607899e-04],   # f=3, s= 1, t=   2
[     1402,  9.8753503746441094e-05,  1.1624141830784459e-04],   # f=3, s= 1, t=   3
[     1261,  2.0857393493693379e-05, -1.4333407898431410e-04],   # f=3, s= 1, t=   4
[     2945,  5.7600265965337831e-05, -7.1668848804472786e-05],   # f=3, s= 1, t=   5
[     1190,  3.1158128562030883e-05, -7.4221508541133025e-05],   # f=3, s= 1, t=   6
[      593,  4.9568824746840162e-05, -5.9443745710499677e-05],   # f=3, s= 1, t=   7
[        0, -6.3916210233836622e-05,  0.0000000000000000e+00],   # f=3, s= 1, t=   8
[        8,  3.2326459763508557e-05,  4.1097067153751401e-05],   # f=3, s= 1, t=   9
[        4,  3.8505374914800768e-05, -2.4690750688454580e-05],   # f=3, s= 1, t=  10
[     1543, -2.6127774810849349e-05,  2.9873942333388641e-05],   # f=3, s= 1, t=  11
[      381, -4.4772545431527864e-06, -3.7396119251182150e-05],   # f=3, s= 1, t=  12
[     1119,  2.2430310753365549e-05, -2.9557761386193919e-05],   # f=3, s= 1, t=  13
[      522,  2.9886224978272480e-05,  2.0259791365064171e-05],   # f=3, s= 1, t=  14
[      452, -1.3106498259692980e-05, -3.2637289766253343e-05],   # f=3, s= 1, t=  15
[     2804, -1.1508149132666279e-05, -2.7316790371599908e-05],   # f=3, s= 1, t=  16
[      310,  3.4738350608151140e-06, -2.1017437982180379e-05],   # f=3, s= 1, t=  17
[      212,  2.1402716283134239e-06, -1.9178262093378069e-05],   # f=3, s= 1, t=  18
[     2733, -2.3283762042050202e-06, -1.5882216399278110e-05],   # f=3, s= 1, t=  19
[     1048,  1.2246829111767810e-05, -9.7512980427413643e-06],   # f=3, s= 1, t=  20
[       12, -7.6135340345472092e-06,  1.2831304341856990e-05],   # f=3, s= 1, t=  21
[       35, -5.7732549173295510e-06,  1.0016273100035810e-05],   # f=3, s= 1, t=  22
[      283, -2.2363814820487809e-06, -9.6218108123675468e-06],   # f=3, s= 1, t=  23
[      239,  4.6765179736328208e-06, -8.5111533123445213e-06],   # f=3, s= 1, t=  24
[      106, -1.9836445213635969e-06,  7.8016311549488939e-06],   # f=3, s= 1, t=  25
[     2875,  5.8334802150818700e-06,  5.1997278670085084e-06],   # f=3, s= 1, t=  26
[     1044,  2.4978011849958161e-06, -6.6276491176842909e-06]    # f=3, s= 1, t=  27
        ],
        [  # f=3, s= 2
[     1331, -1.3993900805601181e-04, -2.9006879879181548e-05],   # f=3, s= 2, t=   0
[        0,  6.7829783611580237e-05,  0.0000000000000000e+00],   # f=3, s= 2, t=   1
[     1261, -1.8612445466711130e-05, -2.8908877766737650e-05],   # f=3, s= 2, t=   2
[       71,  2.9853088727300880e-05, -7.1012915175069879e-06],   # f=3, s= 2, t=   3
[        4,  2.6321188906063949e-05,  1.3464234343799510e-05],   # f=3, s= 2, t=   4
[     1402,  2.3284325058403039e-05, -6.8144038337101704e-06],   # f=3, s= 2, t=   5
[     2945, -1.7707706192590389e-05, -1.5833741011362669e-05],   # f=3, s= 2, t=   6
[     1190, -1.4044910483808361e-05, -1.8816964596903641e-05],   # f=3, s= 2, t=   7
[        8, -3.2757789019448882e-06,  2.1268129579065079e-05],   # f=3, s= 2, t=   8
[     1473, -9.8933125150553904e-06, -1.3083286541079480e-05],   # f=3, s= 2, t=   9
[     1119, -7.2140288335202669e-06, -1.1729524860644351e-05]    # f=3, s= 2, t=  10
        ],
        [  # f=3, s= 3
[     1331, -1.7737816039767501e-05,  2.8055420335966871e-05]    # f=3, s= 3, t=   0
        ]
    ],
    [  # f=4
        [  # f=4, s= 0
[        0, -5.1702307822780000e-02,  0.0000000000000000e+00],   # f=4, s= 0, t=   0
[     1402, -1.3296787157385281e-04,  1.3523784099823399e-04],   # f=4, s= 0, t=   1
[     1543,  1.6300899640352639e-04,  5.5927755421839437e-05],   # f=4, s= 0, t=   2
[     1473, -2.3179975765177140e-05, -9.8184622098344091e-05],   # f=4, s= 0, t=   3
[      522, -1.9097390976957099e-05,  3.6607666661423369e-05],   # f=4, s= 0, t=   4
[      664,  3.2543564829276270e-05,  1.0672339397019160e-06],   # f=4, s= 0, t=   5
[     1331, -2.0841818514170151e-05,  1.2410223144740100e-05],   # f=4, s= 0, t=   6
[        4, -8.0356073429830634e-06,  1.9857752253856971e-05],   # f=4, s= 0, t=   7
[      593, -1.0459691293338759e-05, -1.7857031017673810e-05],   # f=4, s= 0, t=   8
[     1614,  2.0002336741252269e-05,  1.4700260033157791e-06],   # f=4, s= 0, t=   9
[     2875, -3.5880881863000781e-06,  8.0256130260543578e-06],   # f=4, s= 0, t=  10
[     3016,  8.5254320868063325e-06, -1.8479247218425171e-07],   # f=4, s= 0, t=  11
[        8, -5.4541260810660720e-06,  3.9805913226076404e-06],   # f=4, s= 0, t=  12
[      452, -3.7355080445397231e-06,  4.0025013123176398e-06],   # f=4, s= 0, t=  13
[       35, -2.1277230177098351e-06, -4.9607027782605714e-06],   # f=4, s= 0, t=  14
[      137, -3.8233051639183096e-06, -3.1849844866680141e-06]    # f=4, s= 0, t=  15
        ],
        [  # f=4, s= 1
[        0,  1.9166844047183211e-04,  0.0000000000000000e+00],   # f=4, s= 1, t=   0
[     1402,  3.9019701387765488e-05,  3.8671378812923788e-05],   # f=4, s= 1, t=   1
[     1543,  1.5119207079965939e-05, -4.3332491992581899e-05],   # f=4, s= 1, t=   2
[     1331,  5.8922122899993431e-06,  1.0029385340576861e-05],   # f=4, s= 1, t=   3
[      522,  1.0280806202111021e-05,  5.2916455231996232e-06]    # f=4, s= 1, t=   4
        ]
    ],
    [  # f=5
        [  # f=5, s= 0
[        0,  1.3977992515640000e-01,  0.0000000000000000e+00],   # f=5, s= 0, t=   0
[     1402,  1.2883200572372361e-04,  1.3471156433694759e-04],   # f=5, s= 0, t=   1
[     1543, -4.9979310098879121e-05,  1.6180896042292819e-04],   # f=5, s= 0, t=   2
[     1473, -2.2060691195269469e-05, -9.3668360843677731e-05],   # f=5, s= 0, t=   3
[      522,  3.5537595493916502e-05,  1.9797390613745271e-05],   # f=5, s= 0, t=   4
[      664, -7.4709042862279015e-08,  3.1983280221550278e-05],   # f=5, s= 0, t=   5
[     1331,  1.1695376261855580e-05,  2.0798562232570440e-05],   # f=5, s= 0, t=   6
[      593, -9.9523808317349038e-06, -1.7074224315439731e-05],   # f=5, s= 0, t=   7
[     1614, -9.4118113483089596e-07,  1.9684674172624949e-05],   # f=5, s= 0, t=   8
[        4, -1.5932093431444279e-05,  9.4977011164352140e-06],   # f=5, s= 0, t=   9
[     2875,  7.8417800319575014e-06,  3.7228551547920001e-06],   # f=5, s= 0, t=  10
[     3016,  4.3444197029476981e-07,  8.3694483911743903e-06],   # f=5, s= 0, t=  11
[        8, -5.9902893836260643e-06, -4.3498982251067908e-06],   # f=5, s= 0, t=  12
[      137, -3.0046395968887028e-06,  4.7462972937701063e-06],   # f=5, s= 0, t=  13
[      452,  3.9819715036126689e-06,  3.7540995120279268e-06],   # f=5, s= 0, t=  14
[       35, -4.0365965056395762e-08,  5.3159364396026402e-06],   # f=5, s= 0, t=  15
[      141,  1.2649492629014499e-07,  4.8536627302669076e-06],   # f=5, s= 0, t=  16
[     1261,  1.1901325673451860e-06,  4.3994722705730970e-06],   # f=5, s= 0, t=  17
[     2945, -2.7264749792923272e-06, -3.6454310686215340e-06],   # f=5, s= 0, t=  18
[     1685,  6.9410338646847810e-07,  3.3566050467992739e-06],   # f=5, s= 0, t=  19
[      734,  8.2952840171332113e-07,  3.2356209378486682e-06],   # f=5, s= 0, t=  20
[      208,  2.4256846662067760e-06, -9.5181484992359002e-07],   # f=5, s= 0, t=  21
[      279, -2.5271926303381422e-06, -3.7922031983877281e-07],   # f=5, s= 0, t=  22
[     1115,  5.6515925473070844e-07,  2.1078436053637770e-06],   # f=5, s= 0, t=  23
[     1257, -1.3558265037946960e-06,  1.4215822110290420e-06],   # f=5, s= 0, t=  24
[       12,  1.3204389683748120e-07, -1.8835278008457151e-06],   # f=5, s= 0, t=  25
[       71, -1.9042145379339269e-07, -1.8630215116559800e-06],   # f=5, s= 0, t=  26
[      212,  1.6053042001407380e-06,  8.4671630560881122e-07],   # f=5, s= 0, t=  27
[    17405,  1.3113174080732559e-06, -4.9311123653462177e-07],   # f=5, s= 0, t=  28
[    17547,  1.0602949344425471e-06,  8.9833952832208403e-07],   # f=5, s= 0, t=  29
[     2804,  9.2561678972169579e-07,  9.7236819416524451e-07],   # f=5, s= 0, t=  30
[    28337, -9.4657252626238933e-08, -1.1160255453211769e-06],   # f=5, s= 0, t=  31
[    28478,  9.2430781502463316e-07, -6.2265716348343233e-07],   # f=5, s= 0, t=  32
[       75, -9.2130387409465058e-07,  5.9403283989994324e-07],   # f=5, s= 0, t=  33
[      106, -2.0304770919600219e-07, -1.0388792453459521e-06],   # f=5, s= 0, t=  34
[      381,  5.9915368756519242e-07,  8.3180916441107633e-07],   # f=5, s= 0, t=  35
[     1190,  2.7201162062442490e-08,  1.0190741100380801e-06],   # f=5, s= 0, t=  36
[     3087,  2.9873098731161120e-07,  9.6483981362737291e-07],   # f=5, s= 0, t=  37
[       63,  5.8956418822020129e-07, -7.7745585971202028e-07],   # f=5, s= 0, t=  38
[     1186,  2.4117256641652101e-07, -9.2801537744071498e-07],   # f=5, s= 0, t=  39
[       67, -3.2941946279216570e-07,  7.9963658829793748e-07],   # f=5, s= 0, t=  40
[    17476, -7.4361213537121858e-07, -1.2891468267067801e-07],   # f=5, s= 0, t=  41
[      283,  2.8076344536095831e-08,  7.3853891759583075e-07],   # f=5, s= 0, t=  42
[     1756,  3.1161152193772379e-07,  6.3104215205486464e-07],   # f=5, s= 0, t=  43
[      247, -3.0658465478480348e-07, -6.2524043128321752e-07],   # f=5, s= 0, t=  44
[    28407, -2.6107491970896922e-07,  5.4499817140536954e-07],   # f=5, s= 0, t=  45
[      354,  4.9460726403829815e-07,  2.9447867230664199e-07],   # f=5, s= 0, t=  46
[     1421, -3.1957781759109811e-07, -4.3033184225872391e-07],   # f=5, s= 0, t=  47
[     1383,  4.0796989384481828e-07,  3.2647659230997358e-07],   # f=5, s= 0, t=  48
[      133,  1.1815490710962760e-07, -5.0524398677283852e-07],   # f=5, s= 0, t=  49
[      177, -4.6755777925983522e-07,  2.0641392675487111e-07],   # f=5, s= 0, t=  50
[      805,  2.4510909183821792e-07,  4.3308098936875489e-07],   # f=5, s= 0, t=  51
[     1563,  2.0041064699998999e-07, -4.4445832345851582e-07],   # f=5, s= 0, t=  52
[       59, -1.2658294965005670e-07, -4.6671870382330612e-07],   # f=5, s= 0, t=  53
[     1524, -7.8614887212231827e-08,  4.6760439016103602e-07],   # f=5, s= 0, t=  54
[     4348,  4.6767825285878958e-07,  3.5290700804681861e-08],   # f=5, s= 0, t=  55
[       16,  4.2127049178471667e-07, -1.8684252468076170e-07],   # f=5, s= 0, t=  56
[     4489,  1.8055735116211371e-07,  4.2107641051871891e-07],   # f=5, s= 0, t=  57
[       79, -4.0548063730912672e-07, -2.1065225409216059e-07],   # f=5, s= 0, t=  58
[     1398, -1.2435454796683390e-07,  3.9919922981210330e-07],   # f=5, s= 0, t=  59
[      145, -3.0487589061275071e-07,  2.8316787420643652e-07],   # f=5, s= 0, t=  60
[      345,  3.4373900983948560e-07, -2.2193498580565619e-07],   # f=5, s= 0, t=  61
[      275, -7.6720036516022037e-09, -4.0650745658300768e-07],   # f=5, s= 0, t=  62
[     1406, -3.7816142635595027e-07,  9.1996049592574050e-08],   # f=5, s= 0, t=  63
[     1044,  7.4671401050977241e-08,  3.7193907363226299e-07],   # f=5, s= 0, t=  64
[      318,  7.1927751245935557e-08,  3.6928457890682560e-07],   # f=5, s= 0, t=  65
[      416,  4.1798174498739052e-08, -3.6796074151811958e-07],   # f=5, s= 0, t=  66
[      204, -1.5624747021437499e-08, -3.6071195814331651e-07],   # f=5, s= 0, t=  67
[     1539,  3.3199296973509958e-07, -1.0438655407426570e-07],   # f=5, s= 0, t=  68
[     1547,  2.1465474690281490e-07,  2.6396972342688539e-07],   # f=5, s= 0, t=  69
[     1327, -1.5719732197123969e-07,  2.8693854391332578e-07],   # f=5, s= 0, t=  70
[      541,  2.6198971229951878e-07,  1.7511268865594321e-07]    # f=5, s= 0, t=  71
        ],
        [  # f=5, s= 1
[        0,  7.7993297015907634e-05,  0.0000000000000000e+00],   # f=5, s= 1, t=   0
[     1402,  3.9152687150310170e-05, -3.7172486736475712e-05],   # f=5, s= 1, t=   1
[     1543,  4.3031136375758243e-05,  1.3517906141260201e-05],   # f=5, s= 1, t=   2
[     1331,  1.0005173914454110e-05, -5.5510992766576859e-06],   # f=5, s= 1, t=   3
[      522,  5.4799934835759310e-06, -9.9803123122233115e-06],   # f=5, s= 1, t=   4
[     1473, -9.2379228749259814e-06,  2.0049455605464851e-06],   # f=5, s= 1, t=   5
[      664, -3.2895419619690672e-06, -4.5974891215164072e-08],   # f=5, s= 1, t=   6
[        8,  1.3908106265782320e-06, -2.9292856582721550e-06],   # f=5, s= 1, t=   7
[     2875,  1.5562630279974180e-06, -2.5756858183805818e-06],   # f=5, s= 1, t=   8
[        4, -2.4245388783587939e-06, -1.5068482062721909e-06],   # f=5, s= 1, t=   9
[     3016,  2.7144173403334680e-06,  1.5168858882407361e-07],   # f=5, s= 1, t=  10
[     2945, -1.9556025171954529e-06,  1.2707620314125440e-06],   # f=5, s= 1, t=  11
[       35,  1.3344840689621260e-06,  1.4001396639490341e-06],   # f=5, s= 1, t=  12
[      593, -1.5045227426363480e-06,  9.0293863845346919e-07],   # f=5, s= 1, t=  13
[     1614,  1.4931074066475380e-06,  1.3901882598918031e-07],   # f=5, s= 1, t=  14
[     1261,  1.3653087138211129e-06, -3.4669288959066708e-07],   # f=5, s= 1, t=  15
[       12,  1.2315835838874919e-06, -5.9561839689487617e-08],   # f=5, s= 1, t=  16
[      137,  9.0428610682743405e-07,  5.7865365381345505e-07],   # f=5, s= 1, t=  17
[     2804,  5.6172637460901586e-07, -4.7331838000806132e-07],   # f=5, s= 1, t=  18
[      212, -6.0813020695312143e-07,  1.1354550736721950e-07]    # f=5, s= 1, t=  19
        ],
        [  # f=5, s= 2
[     1402, -2.9412500484509639e-06, -8.0598049136548207e-06],   # f=5, s= 2, t=   0
[     1543, -1.1369975198335059e-06, -6.6800431913703534e-06],   # f=5, s= 2, t=   1
[     1331, -5.6564967655275266e-07, -2.8456118566588939e-06],   # f=5, s= 2, t=   2
[        0,  2.4849999999999999e-06,  0.0000000000000000e+00],   # f=5, s= 2, t=   3
[      522, -1.0188268325841671e-06, -1.4049711115470319e-06],   # f=5, s= 2, t=   4
[        4, -1.2789041308060199e-06, -7.7265243668823065e-07],   # f=5, s= 2, t=   5
[        8,  1.2656119872224440e-06, -4.3448350220923409e-07]    # f=5, s= 2, t=   6
        ]
    ]
]

_top_freq = [
    0.5296909622785881e+03,     # Jupiter
    0.2132990811942489e+03,     # Saturn
    0.7478166163181234e+02,     # Uranus
    0.3813297236217556e+02,     # Neptune
    0.2533566020437000e+02      # Pluto
]


def _calc_elliptical_coord(formula, dmu, f, t1):
    el = 0.0
    tpower = 1.0
    s = 0
    for series in formula:
        for term in series:
            term_k, term_c, term_s = term
            if f==1 and s==1 and term_k==0:
                continue
            arg = term_k * dmu * t1
            el += tpower * (term_c*math.cos(arg) + term_s*math.sin(arg))
        tpower *= t1
        s += 1
    return el


def _TopCalcElliptical(body, model, tt):
    # Translated from: TOP2013.f
    # See: https://github.com/cosinekitty/ephemeris/tree/master/top2013
    # Copied from: ftp://ftp.imcce.fr/pub/ephem/planets/top2013
    t1 = tt / 365250.0
    dmu = (_top_freq[0] - _top_freq[1]) / 880.0

    ellip = [
        _calc_elliptical_coord(model[0], dmu, 0, t1),    # a
        _calc_elliptical_coord(model[1], dmu, 1, t1),    # lambda
        _calc_elliptical_coord(model[2], dmu, 2, t1),    # k
        _calc_elliptical_coord(model[3], dmu, 3, t1),    # h
        _calc_elliptical_coord(model[4], dmu, 4, t1),    # q
        _calc_elliptical_coord(model[5], dmu, 5, t1)     # p
    ]

    xl = (ellip[1] + _top_freq[body.value - 4] * t1) % _PI2
    if xl < 0.0:
        xl += _PI2
    ellip[1] = xl
    return ellip


def _TopEcliptic(ellip, time):
    xa, xl, xk, xh, xq, xp = ellip
    xfi = math.sqrt(1.0 - xk*xk - xh*xh)
    xki = math.sqrt(1.0 - xq*xq - xp*xp)
    zr = xk
    zi = xh
    u = 1.0 / (1.0 + xfi)
    ex2 = zr*zr + zi*zi
    ex = math.sqrt(ex2)
    ex3 = ex * ex2
    z1r = zr
    z1i = -zi
    gl = xl % _PI2
    gm = gl - math.atan2(xh, xk)
    e = gl + (ex - 0.125*ex3)*math.sin(gm) + 0.5*ex2*math.sin(2.0*gm) + 0.375*ex3*math.sin(3.0*gm)

    while True:
        z2r = 0.0
        z2i = e
        zteta_r = math.cos(z2i)
        zteta_i = math.sin(z2i)
        z3r = z1r*zteta_r - z1i*zteta_i
        z3i = z1r*zteta_i + z1i*zteta_r
        dl = gl - e + z3i
        rsa = 1.0 - z3r
        e += dl/rsa
        if abs(dl) < 1.0e-15:
            break

    z1r = z3i * u * zr
    z1i = z3i * u * zi
    z2r = +z1i
    z2i = -z1r
    zto_r = (-zr + zteta_r + z2r) / rsa
    zto_i = (-zi + zteta_i + z2i) / rsa
    xcw = zto_r
    xsw = zto_i
    xm = xp*xcw - xq*xsw
    xr = xa*rsa

    return Vector(
        xr*(xcw - 2.0*xp*xm),
        xr*(xsw + 2.0*xq*xm),
        -2.0*xr*xki*xm,
        time
    )


_top_eq_rot = None


def _TopEquatorial(ecl):
    global _top_eq_rot

    if _top_eq_rot is None:
        sdrad = math.radians(1.0 / 3600.0)
        eps = math.radians(23.0 + 26.0/60.0 + 21.41136/3600.0)
        phi = -0.05188 * sdrad
        ceps = math.cos(eps)
        seps = math.sin(eps)
        cphi = math.cos(phi)
        sphi = math.sin(phi)

        _top_eq_rot = [
            [cphi, -sphi*ceps, sphi*seps],
            [sphi, cphi*ceps, -cphi*seps],
            [0.0, seps, ceps]
        ]

    return Vector(
        (_top_eq_rot[0][0] * ecl.x) + (_top_eq_rot[0][1] * ecl.y) + (_top_eq_rot[0][2] * ecl.z),
        (_top_eq_rot[1][0] * ecl.x) + (_top_eq_rot[1][1] * ecl.y) + (_top_eq_rot[1][2] * ecl.z),
        (_top_eq_rot[2][0] * ecl.x) + (_top_eq_rot[2][1] * ecl.y) + (_top_eq_rot[2][2] * ecl.z),
        ecl.t
    )


def _TopPosition(model, body, time):
    ellip = _TopCalcElliptical(body, model, time.tt)
    ecl = _TopEcliptic(ellip, time)
    return _TopEquatorial(ecl)


# END TOP2013
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
    return (x, t, df_dt)

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
    iter = 0
    iter_limit = 20
    calc_fmid = True
    while True:
        iter += 1
        if iter > iter_limit:
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
            (q_x, q_ut, q_df_dt) = q
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

def _AdjustBarycenter(ssb, time, body, pmass):
    shift = pmass / (pmass + _SUN_MASS)
    planet = _CalcVsop(_vsop[body.value], time)
    ssb.x += shift * planet.x
    ssb.y += shift * planet.y
    ssb.z += shift * planet.z

def _CalcSolarSystemBarycenter(time):
    ssb = Vector(0.0, 0.0, 0.0, time)
    _AdjustBarycenter(ssb, time, Body.Jupiter, _JUPITER_MASS)
    _AdjustBarycenter(ssb, time, Body.Saturn,  _SATURN_MASS)
    _AdjustBarycenter(ssb, time, Body.Uranus,  _URANUS_MASS)
    _AdjustBarycenter(ssb, time, Body.Neptune, _NEPTUNE_MASS)
    return ssb

def HelioVector(body, time):
    """Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.

    This function calculates the position of the given celestial body as a vector,
    using the center of the Sun as the origin.  The result is expressed as a Cartesian
    vector in the J2000 equatorial system: the coordinates are based on the mean equator
    of the Earth at noon UTC on 1 January 2000.

    The position is not corrected for light travel time or aberration.
    This is different from the behavior of #GeoVector.

    If given an invalid value for `body`, or the body is `Body.Pluto` and `time` is outside
    the year range 1700..2200, this function raise an exception.

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
        return _TopPosition(_pluto, body, time)

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


def GeoVector(body, time, aberration):
    """Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.

    This function calculates the position of the given celestial body as a vector,
    using the center of the Earth as the origin.  The result is expressed as a Cartesian
    vector in the J2000 equatorial system: the coordinates are based on the mean equator
    of the Earth at noon UTC on 1 January 2000.

    If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
    the year range 1700..2200, this function will raise an exception.

    Unlike #HelioVector, this function always corrects for light travel time.
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

    if not aberration:
        # No aberration, so calculate Earth's position once, at the time of observation.
        earth = _CalcEarth(time)

    # Correct for light-travel time, to get position of body as seen from Earth's center.
    ltime = time
    for iter in range(10):
        h = HelioVector(body, ltime)
        if aberration:
            # Include aberration, so make a good first-order approximation
            # by backdating the Earth's position also.
            # This is confusing, but it works for objects within the Solar System
            # because the distance the Earth moves in that small amount of light
            # travel time (a few minutes to a few hours) is well approximated
            # by a line segment that substends the angle seen from the remote
            # body viewing Earth. That angle is pretty close to the aberration
            # angle of the moving Earth viewing the remote body.
            # In other words, both of the following approximate the aberration angle:
            #    (transverse distance Earth moves) / (distance to body)
            #    (transverse speed of Earth) / (speed of light).
            earth = _CalcEarth(ltime)

        geo = Vector(h.x-earth.x, h.y-earth.y, h.z-earth.z, time)
        ltime2 = time.AddDays(-geo.Length() / _C_AUDAY)
        dt = abs(ltime2.tt - ltime.tt)
        if dt < 1.0e-9:
            return geo

        ltime = ltime2

    raise Error('Light-travel time solver did not converge: dt={}'.format(dt))


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
        return _vector2radec(j2000)
    temp = _precession(0, j2000, time.tt)
    datevect = _nutation(time, 0, temp)
    return _vector2radec(datevect)

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
        If `Refraction.None`, no refraction correct is performed.
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
    rarad = math.radians(ra * 15.0)

    sinlat = math.sin(latrad)
    coslat = math.cos(latrad)
    sinlon = math.sin(lonrad)
    coslon = math.cos(lonrad)
    sindc = math.sin(decrad)
    cosdc = math.cos(decrad)
    sinra = math.sin(rarad)
    cosra = math.cos(rarad)

    uze = [coslat*coslon, coslat*sinlon, sinlat]
    une = [-sinlat*coslon, -sinlat*sinlon, coslat]
    uwe = [sinlon, -coslon, 0.0]

    angle = -15.0 * _sidereal_time(time)
    uz = _spin(angle, uze)
    un = _spin(angle, une)
    uw = _spin(angle, uwe)

    p = [cosdc*cosra, cosdc*sinra, sindc]

    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2]
    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2]
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2]

    proj = math.sqrt(pn*pn + pw*pw)
    az = 0.0
    if proj > 0.0:
        az = math.degrees(-math.atan2(pw, pn))
        if az < 0:
            az += 360
        if az >= 360:
            az -= 360
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
            proj = math.sqrt(pr[0]*pr[0] + pr[1]*pr[1])
            if proj > 0:
                hor_ra = math.degrees(math.atan2(pr[1], pr[0])) / 15
                if hor_ra < 0:
                    hor_ra += 24
                if hor_ra >= 24:
                    hor_ra -= 24
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

    if refraction == Refraction.Normal or refraction == Refraction.JplHorizons:
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
    else:
        # No refraction, or the refraction option is invalid.
        refr = 0.0
    return refr

def InverseRefractionAngle(refraction, bent_altitude):
    """Calculates the inverse of an atmospheric refraction angle.

    Given an observed altitude angle that includes atmospheric refraction,
    calculate the negative angular correction to obtain the unrefracted
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
    ex : float
        Cartesian x-coordinate: in the direction of the equinox along the ecliptic plane.
    ey : float
        Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox.
    ez : float
        Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north.
    elat : float
        Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
    elon : float
        Longitude in degrees around the ecliptic plane prograde from the equinox.
    """
    def __init__(self, ex, ey, ez, elat, elon):
        self.ex = ex
        self.ey = ey
        self.ez = ez
        self.elat = elat
        self.elon = elon

def _RotateEquatorialToEcliptic(pos, obliq_radians):
    cos_ob = math.cos(obliq_radians)
    sin_ob = math.sin(obliq_radians)
    ex = +pos[0]
    ey = +pos[1]*cos_ob + pos[2]*sin_ob
    ez = -pos[1]*sin_ob + pos[2]*cos_ob
    xyproj = math.sqrt(ex*ex + ey*ey)
    if xyproj > 0.0:
        elon = math.degrees(math.atan2(ey, ex))
        if elon < 0.0:
            elon += 360.0
    else:
        elon = 0.0
    elat = math.degrees(math.atan2(ez, xyproj))
    return EclipticCoordinates(ex, ey, ez, elat, elon)

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
    adjusted_time = time.AddDays(-1.0 / _C_AUDAY)
    earth2000 = _CalcEarth(adjusted_time)
    sun2000 = [-earth2000.x, -earth2000.y, -earth2000.z]

    # Convert to equatorial Cartesian coordinates of date.
    stemp = _precession(0.0, sun2000, adjusted_time.tt)
    sun_ofdate = _nutation(adjusted_time, 0, stemp)

    # Convert equatorial coordinates to ecliptic coordinates.
    true_obliq = math.radians(adjusted_time._etilt().tobl)
    return _RotateEquatorialToEcliptic(sun_ofdate, true_obliq)

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
    return _RotateEquatorialToEcliptic([equ.x, equ.y, equ.z], ob2000)

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
    return _AngleBetween(sv, bv)

def LongitudeFromSun(body, time):
    """Returns a body's ecliptic longitude with respect to the Sun, as seen from the Earth.

    This function can be used to determine where a planet appears around the ecliptic plane
    (the plane of the Earth's orbit around the Sun) as seen from the Earth,
    relative to the Sun's apparent position.

    The angle starts at 0 when the body and the Sun are at the same ecliptic longitude
    as seen from the Earth. The angle increases in the prograde direction
    (the direction that the planets orbit the Sun and the Moon orbits the Earth).

    When the angle is 180 degrees, it means the Sun and the body appear on opposite sides
    of the sky for an Earthly observer. When `body` is a planet whose orbit around the
    Sun is farther than the Earth's, 180 degrees indicates opposition. For the Moon,
    it indicates a full moon.

    The angle keeps increasing up to 360 degrees as the body's apparent prograde
    motion continues relative to the Sun. When the angle reaches 360 degrees, it starts
    over at 0 degrees.

    Values between 0 and 180 degrees indicate that the body is visible in the evening sky
    after sunset.  Values between 180 degrees and 360 degrees indicate that the body
    is visible in the morning sky before sunrise.

    Parameters
    ----------
    body : Body
        The celestial body for which to find longitude from the Sun.

    time : Time
        The date and time of the observation.

    Returns
    -------
    float
        An angle in degrees in the range [0, 360).
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    sv = GeoVector(Body.Sun, time, False)
    se = Ecliptic(sv)
    bv = GeoVector(body, time, False)
    be = Ecliptic(bv)
    return _NormalizeLongitude(be.elon - se.elon)

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
    angle = LongitudeFromSun(body, time)
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
    if body == Body.Moon or body == Body.Sun:
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
    iter = 0
    while iter < 100:
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
        iter += 1
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
    iter = 1
    while iter <= 2:
        plon = EclipticLongitude(body, startTime)
        elon = EclipticLongitude(Body.Earth, startTime)
        rlon = _LongitudeOffset(plon - elon)    # clamp to (-180, +180]

        # The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        # because there is a cusp there that causes a discontinuity in the derivative.
        # So we need to guard against searching near such times.
        if rlon >= -s1 and rlon < +s1:
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
        iter += 1


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
    return Search(_sun_offset, targetLon, startTime, t2, 1.0)

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
    return LongitudeFromSun(Body.Moon, time)

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
         The number of days after `startTime` that limits the time window for the search.

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
    # I have seen up to 0.826 days away from the simple prediction.
    # To be safe, we take the predicted time of the event and search
    # +/-0.9 days around it (a 1.8-day wide window).
    # But we must return None if the final result goes beyond limitDays after startTime.
    uncertainty = 0.9
    ya = _moon_offset(targetLon, startTime)
    if ya > 0.0:
        ya -= 360.0     # force searching forward in time, not backward
    est_dt = -(_MEAN_SYNODIC_MONTH * ya) / 360.0
    dt1 = est_dt - uncertainty
    if dt1 > limitDays:
        return None     # not possible for moon phase to occur within the specified window
    dt2 = min(limitDays, est_dt + uncertainty)
    t1 = startTime.AddDays(dt1)
    t2 = startTime.AddDays(dt2)
    return Search(_moon_offset, targetLon, t1, t2, 1.0)

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
    ring_tilt : float
        For Saturn, the tilt angle in degrees of its rings as seen from Earth.
        When the `ring_tilt` is very close to 0, it means the rings are edge-on
        as seen from observers on the Earth, and are thus very difficult to see.
        For bodies other than Saturn, `ring_tilt` is `None`.
    """
    def __init__(self, time, mag, phase, helio_dist, geo_dist, gc, hc, ring_tilt):
        self.time = time
        self.mag = mag
        self.phase_angle = phase
        self.phase_fraction = (1.0 + math.cos(math.radians(phase))) / 2.0
        self.helio_dist = helio_dist
        self.geo_dist = geo_dist
        self.gc = gc
        self.hc = hc
        self.ring_tilt = ring_tilt

def _MoonMagnitude(phase, helio_dist, geo_dist):
    # https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
    rad = math.radians(phase)
    mag = -12.717 + 1.49*abs(rad) + 0.0431*(rad**4)
    moon_mean_distance_au = 385000.6 / _KM_PER_AU
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
        phase = _AngleBetween(gc, hc)

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
    return IlluminationInfo(time, mag, phase, helio_dist, geo_dist, gc, hc, ring_tilt)

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

    iter = 1
    while iter <= 2:
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
        iter += 1

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

def SearchHourAngle(body, observer, hourAngle, startTime):
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

    Returns
    -------
    HourAngleEvent
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if hourAngle < 0.0 or hourAngle >= 24.0:
        raise Error('Invalid hour angle.')

    iter = 0
    time = startTime
    while True:
        iter += 1
        # Calculate Greenwich Apparent Sidereal Time (GAST) at the given time.
        gast = _sidereal_time(time)
        ofdate = Equator(body, time, observer, True, True)

        # Calculate the adjustment needed in sidereal time to bring
        # the hour angle to the desired value.
        delta_sidereal_hours = math.fmod(((hourAngle + ofdate.ra - observer.longitude/15) - gast), 24.0)
        if iter == 1:
            # On the first iteration, always search forward in time.
            if delta_sidereal_hours < 0.0:
                delta_sidereal_hours += 24.0
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
        Limits how many days to search for a rise or set time.
        To limit a rise or set time to the same day, you can use a value of 1 day.
        In cases where you want to find the next rise or set time no matter how far
        in the future (for example, for an observer near the south pole), you can pass
        in a larger value like 365.

    Returns
    -------
    Time or `None`
        If the rise or set time is found within the specified time window,
        this function returns that time. Otherwise, it returns `None`.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    elif body == Body.Sun:
        body_radius = _SUN_RADIUS_AU
    elif body == Body.Moon:
        body_radius = _MOON_EQUATORIAL_RADIUS_AU
    else:
        body_radius = 0.0

    if direction == Direction.Rise:
        ha_before = 12.0    # minimum altitude (bottom) happens BEFORE the body rises.
        ha_after  =  0.0    # maximum altitude (culmination) happens AFTER the body rises.
    elif direction == Direction.Set:
        ha_before =  0.0    # culmination happens BEFORE the body sets.
        ha_after  = 12.0    # bottom happens AFTER the body sets.
    else:
        raise Error('Invalid value for direction parameter')

    context = _peak_altitude_context(body, direction, observer, body_radius)

    # See if the body is currently above/below the horizon.
    # If we are looking for next rise time and the body is below the horizon,
    # we use the current time as the lower time bound and the next culmination
    # as the upper bound.
    # If the body is above the horizon, we search for the next bottom and use it
    # as the lower bound and the next culmination after that bottom as the upper bound.
    # The same logic applies for finding set times, only we swap the hour angles.
    time_start = startTime
    alt_before = _peak_altitude(context, time_start)
    if alt_before > 0.0:
        # We are past the sought event, so we have to wait for the next "before" event (culm/bottom).
        evt_before = SearchHourAngle(body, observer, ha_before, time_start)
        time_before = evt_before.time
        alt_before = _peak_altitude(context, time_before)
    else:
        # We are before or at the sought ebvent, so we find the next "after" event (bottom/culm),
        # and use the current time as the "before" event.
        time_before = time_start

    evt_after = SearchHourAngle(body, observer, ha_after, time_before)
    alt_after = _peak_altitude(context, evt_after.time)

    while True:
        if alt_before <= 0.0 and alt_after > 0.0:
            # Search between the "before time" and the "after time" for the desired event.
            event_time = Search(_peak_altitude, context, time_before, evt_after.time, 1.0)
            if event_time is not None:
                return event_time
        # We didn't find the desired event, so use the "after" time to find the next "before" event.
        evt_before = SearchHourAngle(body, observer, ha_before, evt_after.time)
        evt_after = SearchHourAngle(body, observer, ha_after, evt_before.time)
        if evt_before.time.ut >= time_start.ut + limitDays:
            return None
        time_before = evt_before.time
        alt_before = _peak_altitude(context, evt_before.time)
        alt_after = _peak_altitude(context, evt_after.time)

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

def _FindSeasonChange(targetLon, year, month, day):
    startTime = Time.Make(year, month, day, 0, 0, 0)
    time = SearchSunLongitude(targetLon, startTime, 4.0)
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
    mar_equinox = _FindSeasonChange(0, year, 3, 19)
    jun_solstice = _FindSeasonChange(90, year, 6, 19)
    sep_equinox = _FindSeasonChange(180, year, 9, 21)
    dec_solstice = _FindSeasonChange(270, year, 12, 20)
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
        self.dist_km = dist_au * _KM_PER_AU

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
    iter = 0
    while iter * increment < 2.0 * _MEAN_SYNODIC_MONTH:
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
        iter += 1

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
    next = SearchLunarApsis(time)
    # Verify that we found the opposite apsis from the previous one.
    if apsis.kind not in [ApsisKind.Apocenter, ApsisKind.Pericenter]:
        raise Error('Parameter "apsis" contains an invalid "kind" value.')
    if next.kind.value + apsis.kind.value != 1:
        raise InternalError()   # should have found opposite apsis kind
    return next


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
    iter = 0
    while iter * increment < 2 * orbit_period_days:
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
        iter += 1
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
    next = SearchPlanetApsis(body, time)
    # Verify that we found the opposite apsis from the previous one.
    if next.kind.value + apsis.kind.value != 1:
        raise InternalError()   # should have found opposite planetary apsis type
    return next


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


def VectorFromEquator(equ, time):
    """Given angular equatorial coordinates in `equ`, calculates equatorial vector.

    Parameters
    ----------
    equ : Equatorial
        Angular equatorial coordinates to be converted to a vector.
    time : Time
        The date and time of the observation. This is needed because the returned
        vector object requires a valid time value when passed to certain other functions.

    Returns
    -------
    Vector
        A vector in the equatorial system.
    """
    return VectorFromSphere(Spherical(equ.dec, 15.0 * equ.ra, equ.dist), time)


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
    return Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist)


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
    this function returns the matrix that reverses that trasnform.

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
    prec = _precession_rot(0.0, time.tt)
    nut = _nutation_rot(time, 0)
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
    nut = _nutation_rot(time, 1)
    prec = _precession_rot(time.tt, 0.0)
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
    spin_angle = -15.0 * _sidereal_time(time)
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
        A rotation matrix that converts HOR to EQD at `time` and for `observer`.
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
    RotationMatrix
        A rotation matrix that converts HOR to ECL.
    -------
    """
    rot = Rotation_ECL_HOR(time, observer)
    return InverseRotation(rot)

class ConstellationInfo:
    """Reports the constellation that a given celestial point lies within.

    The #Constellation function returns this struct
    to report which constellation corresponds with a given point in the sky.
    Constellations are defined with respect to the B1875 equatorial system
    per IAU standard. Although `Constellation` requires J2000 equatorial
    coordinates, the struct contains converted B1875 coordinates for reference.

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
    ( 83,  0.00000000000000, 24.00000000000000, 88.00000000000000 )    # UMi
,   ( 83,  8.00000000000000, 14.50000000000000, 86.50000000000000 )    # UMi
,   ( 83, 21.00000000000000, 23.00000000000000, 86.16666666666667 )    # UMi
,   ( 83, 18.00000000000000, 21.00000000000000, 86.00000000000000 )    # UMi
,   ( 15,  0.00000000000000,  8.00000000000000, 85.00000000000000 )    # Cep
,   ( 10,  9.16666666666667, 10.66666666666667, 82.00000000000000 )    # Cam
,   ( 15,  0.00000000000000,  5.00000000000000, 80.00000000000000 )    # Cep
,   ( 10, 10.66666666666667, 14.50000000000000, 80.00000000000000 )    # Cam
,   ( 83, 17.50000000000000, 18.00000000000000, 80.00000000000000 )    # UMi
,   ( 33, 20.16666666666667, 21.00000000000000, 80.00000000000000 )    # Dra
,   ( 15,  0.00000000000000,  3.50833333333333, 77.00000000000000 )    # Cep
,   ( 10, 11.50000000000000, 13.58333333333333, 77.00000000000000 )    # Cam
,   ( 83, 16.53333333333333, 17.50000000000000, 75.00000000000000 )    # UMi
,   ( 15, 20.16666666666667, 20.66666666666667, 75.00000000000000 )    # Cep
,   ( 10,  7.96666666666667,  9.16666666666667, 73.50000000000000 )    # Cam
,   ( 33,  9.16666666666667, 11.33333333333333, 73.50000000000000 )    # Dra
,   ( 83, 13.00000000000000, 16.53333333333333, 70.00000000000000 )    # UMi
,   ( 13,  3.10000000000000,  3.41666666666667, 68.00000000000000 )    # Cas
,   ( 33, 20.41666666666667, 20.66666666666667, 67.00000000000000 )    # Dra
,   ( 33, 11.33333333333333, 12.00000000000000, 66.50000000000000 )    # Dra
,   ( 15,  0.00000000000000,  0.33333333333333, 66.00000000000000 )    # Cep
,   ( 83, 14.00000000000000, 15.66666666666667, 66.00000000000000 )    # UMi
,   ( 15, 23.58333333333333, 24.00000000000000, 66.00000000000000 )    # Cep
,   ( 33, 12.00000000000000, 13.50000000000000, 64.00000000000000 )    # Dra
,   ( 33, 13.50000000000000, 14.41666666666667, 63.00000000000000 )    # Dra
,   ( 15, 23.16666666666667, 23.58333333333333, 63.00000000000000 )    # Cep
,   ( 10,  6.10000000000000,  7.00000000000000, 62.00000000000000 )    # Cam
,   ( 33, 20.00000000000000, 20.41666666666667, 61.50000000000000 )    # Dra
,   ( 15, 20.53666666666667, 20.60000000000000, 60.91666666666666 )    # Cep
,   ( 10,  7.00000000000000,  7.96666666666667, 60.00000000000000 )    # Cam
,   ( 82,  7.96666666666667,  8.41666666666667, 60.00000000000000 )    # UMa
,   ( 33, 19.76666666666667, 20.00000000000000, 59.50000000000000 )    # Dra
,   ( 15, 20.00000000000000, 20.53666666666667, 59.50000000000000 )    # Cep
,   ( 15, 22.86666666666667, 23.16666666666667, 59.08333333333334 )    # Cep
,   ( 13,  0.00000000000000,  2.43333333333333, 58.50000000000000 )    # Cas
,   ( 33, 19.41666666666667, 19.76666666666667, 58.00000000000000 )    # Dra
,   ( 13,  1.70000000000000,  1.90833333333333, 57.50000000000000 )    # Cas
,   ( 13,  2.43333333333333,  3.10000000000000, 57.00000000000000 )    # Cas
,   ( 10,  3.10000000000000,  3.16666666666667, 57.00000000000000 )    # Cam
,   ( 15, 22.31666666666667, 22.86666666666667, 56.25000000000000 )    # Cep
,   ( 10,  5.00000000000000,  6.10000000000000, 56.00000000000000 )    # Cam
,   ( 82, 14.03333333333333, 14.41666666666667, 55.50000000000000 )    # UMa
,   ( 33, 14.41666666666667, 19.41666666666667, 55.50000000000000 )    # Dra
,   ( 10,  3.16666666666667,  3.33333333333333, 55.00000000000000 )    # Cam
,   ( 15, 22.13333333333333, 22.31666666666667, 55.00000000000000 )    # Cep
,   ( 15, 20.60000000000000, 21.96666666666667, 54.83333333333334 )    # Cep
,   ( 13,  0.00000000000000,  1.70000000000000, 54.00000000000000 )    # Cas
,   ( 50,  6.10000000000000,  6.50000000000000, 54.00000000000000 )    # Lyn
,   ( 82, 12.08333333333333, 13.50000000000000, 53.00000000000000 )    # UMa
,   ( 33, 15.25000000000000, 15.75000000000000, 53.00000000000000 )    # Dra
,   ( 15, 21.96666666666667, 22.13333333333333, 52.75000000000000 )    # Cep
,   ( 10,  3.33333333333333,  5.00000000000000, 52.50000000000000 )    # Cam
,   ( 13, 22.86666666666667, 23.33333333333333, 52.50000000000000 )    # Cas
,   ( 33, 15.75000000000000, 17.00000000000000, 51.50000000000000 )    # Dra
,   ( 62,  2.04166666666667,  2.51666666666667, 50.50000000000000 )    # Per
,   ( 33, 17.00000000000000, 18.23333333333333, 50.50000000000000 )    # Dra
,   ( 13,  0.00000000000000,  1.36666666666667, 50.00000000000000 )    # Cas
,   ( 62,  1.36666666666667,  1.66666666666667, 50.00000000000000 )    # Per
,   ( 50,  6.50000000000000,  6.80000000000000, 50.00000000000000 )    # Lyn
,   ( 13, 23.33333333333333, 24.00000000000000, 50.00000000000000 )    # Cas
,   ( 82, 13.50000000000000, 14.03333333333333, 48.50000000000000 )    # UMa
,   ( 13,  0.00000000000000,  1.11666666666667, 48.00000000000000 )    # Cas
,   ( 13, 23.58333333333333, 24.00000000000000, 48.00000000000000 )    # Cas
,   ( 39, 18.17500000000000, 18.23333333333333, 47.50000000000000 )    # Her
,   ( 33, 18.23333333333333, 19.08333333333333, 47.50000000000000 )    # Dra
,   ( 30, 19.08333333333333, 19.16666666666667, 47.50000000000000 )    # Cyg
,   ( 62,  1.66666666666667,  2.04166666666667, 47.00000000000000 )    # Per
,   ( 82,  8.41666666666667,  9.16666666666667, 47.00000000000000 )    # UMa
,   ( 13,  0.16666666666667,  0.86666666666667, 46.00000000000000 )    # Cas
,   ( 82, 12.00000000000000, 12.08333333333333, 45.00000000000000 )    # UMa
,   ( 50,  6.80000000000000,  7.36666666666667, 44.50000000000000 )    # Lyn
,   ( 30, 21.90833333333333, 21.96666666666667, 44.00000000000000 )    # Cyg
,   ( 30, 21.87500000000000, 21.90833333333333, 43.75000000000000 )    # Cyg
,   ( 30, 19.16666666666667, 19.40000000000000, 43.50000000000000 )    # Cyg
,   ( 82,  9.16666666666667, 10.16666666666667, 42.00000000000000 )    # UMa
,   ( 82, 10.16666666666667, 10.78333333333333, 40.00000000000000 )    # UMa
,   (  8, 15.43333333333333, 15.75000000000000, 40.00000000000000 )    # Boo
,   ( 39, 15.75000000000000, 16.33333333333333, 40.00000000000000 )    # Her
,   ( 50,  9.25000000000000,  9.58333333333333, 39.75000000000000 )    # Lyn
,   (  0,  0.00000000000000,  2.51666666666667, 36.75000000000000 )    # And
,   ( 62,  2.51666666666667,  2.56666666666667, 36.75000000000000 )    # Per
,   ( 51, 19.35833333333333, 19.40000000000000, 36.50000000000000 )    # Lyr
,   ( 62,  4.50000000000000,  4.69166666666667, 36.00000000000000 )    # Per
,   ( 30, 21.73333333333333, 21.87500000000000, 36.00000000000000 )    # Cyg
,   ( 44, 21.87500000000000, 22.00000000000000, 36.00000000000000 )    # Lac
,   (  7,  6.53333333333333,  7.36666666666667, 35.50000000000000 )    # Aur
,   ( 50,  7.36666666666667,  7.75000000000000, 35.50000000000000 )    # Lyn
,   (  0,  0.00000000000000,  2.00000000000000, 35.00000000000000 )    # And
,   ( 44, 22.00000000000000, 22.81666666666667, 35.00000000000000 )    # Lac
,   ( 44, 22.81666666666667, 22.86666666666667, 34.50000000000000 )    # Lac
,   (  0, 22.86666666666667, 23.50000000000000, 34.50000000000000 )    # And
,   ( 62,  2.56666666666667,  2.71666666666667, 34.00000000000000 )    # Per
,   ( 82, 10.78333333333333, 11.00000000000000, 34.00000000000000 )    # UMa
,   ( 29, 12.00000000000000, 12.33333333333333, 34.00000000000000 )    # CVn
,   ( 50,  7.75000000000000,  9.25000000000000, 33.50000000000000 )    # Lyn
,   ( 48,  9.25000000000000,  9.88333333333333, 33.50000000000000 )    # LMi
,   (  0,  0.71666666666667,  1.40833333333333, 33.00000000000000 )    # And
,   (  8, 15.18333333333333, 15.43333333333333, 33.00000000000000 )    # Boo
,   (  0, 23.50000000000000, 23.75000000000000, 32.08333333333334 )    # And
,   ( 29, 12.33333333333333, 13.25000000000000, 32.00000000000000 )    # CVn
,   (  0, 23.75000000000000, 24.00000000000000, 31.33333333333333 )    # And
,   ( 29, 13.95833333333333, 14.03333333333333, 30.75000000000000 )    # CVn
,   ( 80,  2.41666666666667,  2.71666666666667, 30.66666666666667 )    # Tri
,   ( 62,  2.71666666666667,  4.50000000000000, 30.66666666666667 )    # Per
,   (  7,  4.50000000000000,  4.75000000000000, 30.00000000000000 )    # Aur
,   ( 51, 18.17500000000000, 19.35833333333333, 30.00000000000000 )    # Lyr
,   ( 82, 11.00000000000000, 12.00000000000000, 29.00000000000000 )    # UMa
,   ( 30, 19.66666666666667, 20.91666666666667, 29.00000000000000 )    # Cyg
,   (  7,  4.75000000000000,  5.88333333333333, 28.50000000000000 )    # Aur
,   ( 48,  9.88333333333333, 10.50000000000000, 28.50000000000000 )    # LMi
,   ( 29, 13.25000000000000, 13.95833333333333, 28.50000000000000 )    # CVn
,   (  0,  0.00000000000000,  0.06666666666667, 28.00000000000000 )    # And
,   ( 80,  1.40833333333333,  1.66666666666667, 28.00000000000000 )    # Tri
,   (  7,  5.88333333333333,  6.53333333333333, 28.00000000000000 )    # Aur
,   ( 37,  7.88333333333333,  8.00000000000000, 28.00000000000000 )    # Gem
,   ( 30, 20.91666666666667, 21.73333333333333, 28.00000000000000 )    # Cyg
,   ( 30, 19.25833333333333, 19.66666666666667, 27.50000000000000 )    # Cyg
,   ( 80,  1.91666666666667,  2.41666666666667, 27.25000000000000 )    # Tri
,   ( 25, 16.16666666666667, 16.33333333333333, 27.00000000000000 )    # CrB
,   (  8, 15.08333333333333, 15.18333333333333, 26.00000000000000 )    # Boo
,   ( 25, 15.18333333333333, 16.16666666666667, 26.00000000000000 )    # CrB
,   ( 51, 18.36666666666667, 18.86666666666667, 26.00000000000000 )    # Lyr
,   ( 48, 10.75000000000000, 11.00000000000000, 25.50000000000000 )    # LMi
,   ( 51, 18.86666666666667, 19.25833333333333, 25.50000000000000 )    # Lyr
,   ( 80,  1.66666666666667,  1.91666666666667, 25.00000000000000 )    # Tri
,   ( 66,  0.71666666666667,  0.85000000000000, 23.75000000000000 )    # Psc
,   ( 48, 10.50000000000000, 10.75000000000000, 23.50000000000000 )    # LMi
,   ( 87, 21.25000000000000, 21.41666666666667, 23.50000000000000 )    # Vul
,   ( 77,  5.70000000000000,  5.88333333333333, 22.83333333333333 )    # Tau
,   (  0,  0.06666666666667,  0.14166666666667, 22.00000000000000 )    # And
,   ( 73, 15.91666666666667, 16.03333333333333, 22.00000000000000 )    # Ser
,   ( 37,  5.88333333333333,  6.21666666666667, 21.50000000000000 )    # Gem
,   ( 87, 19.83333333333333, 20.25000000000000, 21.25000000000000 )    # Vul
,   ( 87, 18.86666666666667, 19.25000000000000, 21.08333333333333 )    # Vul
,   (  0,  0.14166666666667,  0.85000000000000, 21.00000000000000 )    # And
,   ( 87, 20.25000000000000, 20.56666666666667, 20.50000000000000 )    # Vul
,   ( 37,  7.80833333333333,  7.88333333333333, 20.00000000000000 )    # Gem
,   ( 87, 20.56666666666667, 21.25000000000000, 19.50000000000000 )    # Vul
,   ( 87, 19.25000000000000, 19.83333333333333, 19.16666666666667 )    # Vul
,   (  6,  3.28333333333333,  3.36666666666667, 19.00000000000000 )    # Ari
,   ( 75, 18.86666666666667, 19.00000000000000, 18.50000000000000 )    # Sge
,   ( 59,  5.70000000000000,  5.76666666666667, 18.00000000000000 )    # Ori
,   ( 37,  6.21666666666667,  6.30833333333333, 17.50000000000000 )    # Gem
,   ( 75, 19.00000000000000, 19.83333333333333, 16.16666666666667 )    # Sge
,   ( 77,  4.96666666666667,  5.33333333333333, 16.00000000000000 )    # Tau
,   ( 39, 15.91666666666667, 16.08333333333333, 16.00000000000000 )    # Her
,   ( 75, 19.83333333333333, 20.25000000000000, 15.75000000000000 )    # Sge
,   ( 77,  4.61666666666667,  4.96666666666667, 15.50000000000000 )    # Tau
,   ( 77,  5.33333333333333,  5.60000000000000, 15.50000000000000 )    # Tau
,   ( 23, 12.83333333333333, 13.50000000000000, 15.00000000000000 )    # Com
,   ( 39, 17.25000000000000, 18.25000000000000, 14.33333333333333 )    # Her
,   ( 23, 11.86666666666667, 12.83333333333333, 14.00000000000000 )    # Com
,   ( 37,  7.50000000000000,  7.80833333333333, 13.50000000000000 )    # Gem
,   ( 39, 16.75000000000000, 17.25000000000000, 12.83333333333333 )    # Her
,   ( 61,  0.00000000000000,  0.14166666666667, 12.50000000000000 )    # Peg
,   ( 77,  5.60000000000000,  5.76666666666667, 12.50000000000000 )    # Tau
,   ( 37,  7.00000000000000,  7.50000000000000, 12.50000000000000 )    # Gem
,   ( 61, 21.11666666666667, 21.33333333333333, 12.50000000000000 )    # Peg
,   ( 37,  6.30833333333333,  6.93333333333333, 12.00000000000000 )    # Gem
,   ( 39, 18.25000000000000, 18.86666666666667, 12.00000000000000 )    # Her
,   ( 31, 20.87500000000000, 21.05000000000000, 11.83333333333333 )    # Del
,   ( 61, 21.05000000000000, 21.11666666666667, 11.83333333333333 )    # Peg
,   ( 45, 11.51666666666667, 11.86666666666667, 11.00000000000000 )    # Leo
,   ( 59,  6.24166666666667,  6.30833333333333, 10.00000000000000 )    # Ori
,   ( 37,  6.93333333333333,  7.00000000000000, 10.00000000000000 )    # Gem
,   ( 21,  7.80833333333333,  7.92500000000000, 10.00000000000000 )    # Cnc
,   ( 61, 23.83333333333333, 24.00000000000000, 10.00000000000000 )    # Peg
,   (  6,  1.66666666666667,  3.28333333333333,  9.91666666666667 )    # Ari
,   ( 31, 20.14166666666667, 20.30000000000000,  8.50000000000000 )    # Del
,   (  8, 13.50000000000000, 15.08333333333333,  8.00000000000000 )    # Boo
,   ( 61, 22.75000000000000, 23.83333333333333,  7.50000000000000 )    # Peg
,   ( 21,  7.92500000000000,  9.25000000000000,  7.00000000000000 )    # Cnc
,   ( 45,  9.25000000000000, 10.75000000000000,  7.00000000000000 )    # Leo
,   ( 58, 18.25000000000000, 18.66222222222222,  6.25000000000000 )    # Oph
,   (  3, 18.66222222222222, 18.86666666666667,  6.25000000000000 )    # Aql
,   ( 31, 20.83333333333333, 20.87500000000000,  6.00000000000000 )    # Del
,   ( 20,  7.00000000000000,  7.01666666666667,  5.50000000000000 )    # CMi
,   ( 73, 18.25000000000000, 18.42500000000000,  4.50000000000000 )    # Ser
,   ( 39, 16.08333333333333, 16.75000000000000,  4.00000000000000 )    # Her
,   ( 58, 18.25000000000000, 18.42500000000000,  3.00000000000000 )    # Oph
,   ( 61, 21.46666666666667, 21.66666666666667,  2.75000000000000 )    # Peg
,   ( 66,  0.00000000000000,  2.00000000000000,  2.00000000000000 )    # Psc
,   ( 73, 18.58333333333333, 18.86666666666667,  2.00000000000000 )    # Ser
,   ( 31, 20.30000000000000, 20.83333333333333,  2.00000000000000 )    # Del
,   ( 34, 20.83333333333333, 21.33333333333333,  2.00000000000000 )    # Equ
,   ( 61, 21.33333333333333, 21.46666666666667,  2.00000000000000 )    # Peg
,   ( 61, 22.00000000000000, 22.75000000000000,  2.00000000000000 )    # Peg
,   ( 61, 21.66666666666667, 22.00000000000000,  1.75000000000000 )    # Peg
,   ( 20,  7.01666666666667,  7.20000000000000,  1.50000000000000 )    # CMi
,   ( 77,  3.58333333333333,  4.61666666666667,  0.00000000000000 )    # Tau
,   ( 59,  4.61666666666667,  4.66666666666667,  0.00000000000000 )    # Ori
,   ( 20,  7.20000000000000,  8.08333333333333,  0.00000000000000 )    # CMi
,   ( 85, 14.66666666666667, 15.08333333333333,  0.00000000000000 )    # Vir
,   ( 58, 17.83333333333333, 18.25000000000000,  0.00000000000000 )    # Oph
,   ( 16,  2.65000000000000,  3.28333333333333, -1.75000000000000 )    # Cet
,   ( 77,  3.28333333333333,  3.58333333333333, -1.75000000000000 )    # Tau
,   ( 73, 15.08333333333333, 16.26666666666667, -3.25000000000000 )    # Ser
,   ( 59,  4.66666666666667,  5.08333333333333, -4.00000000000000 )    # Ori
,   ( 59,  5.83333333333333,  6.24166666666667, -4.00000000000000 )    # Ori
,   ( 73, 17.83333333333333, 17.96666666666667, -4.00000000000000 )    # Ser
,   ( 73, 18.25000000000000, 18.58333333333333, -4.00000000000000 )    # Ser
,   (  3, 18.58333333333333, 18.86666666666667, -4.00000000000000 )    # Aql
,   ( 66, 22.75000000000000, 23.83333333333333, -4.00000000000000 )    # Psc
,   ( 45, 10.75000000000000, 11.51666666666667, -6.00000000000000 )    # Leo
,   ( 85, 11.51666666666667, 11.83333333333333, -6.00000000000000 )    # Vir
,   ( 66,  0.00000000000000,  0.33333333333333, -7.00000000000000 )    # Psc
,   ( 66, 23.83333333333333, 24.00000000000000, -7.00000000000000 )    # Psc
,   ( 85, 14.25000000000000, 14.66666666666667, -8.00000000000000 )    # Vir
,   ( 58, 15.91666666666667, 16.26666666666667, -8.00000000000000 )    # Oph
,   (  3, 20.00000000000000, 20.53333333333333, -9.00000000000000 )    # Aql
,   (  4, 21.33333333333333, 21.86666666666667, -9.00000000000000 )    # Aqr
,   ( 58, 17.16666666666667, 17.96666666666667, -10.00000000000000 )    # Oph
,   ( 54,  5.83333333333333,  8.08333333333333, -11.00000000000000 )    # Mon
,   ( 35,  4.91666666666667,  5.08333333333333, -11.00000000000000 )    # Eri
,   ( 59,  5.08333333333333,  5.83333333333333, -11.00000000000000 )    # Ori
,   ( 41,  8.08333333333333,  8.36666666666667, -11.00000000000000 )    # Hya
,   ( 74,  9.58333333333333, 10.75000000000000, -11.00000000000000 )    # Sex
,   ( 85, 11.83333333333333, 12.83333333333333, -11.00000000000000 )    # Vir
,   ( 58, 17.58333333333333, 17.66666666666667, -11.66666666666667 )    # Oph
,   (  3, 18.86666666666667, 20.00000000000000, -12.03333333333333 )    # Aql
,   ( 35,  4.83333333333333,  4.91666666666667, -14.50000000000000 )    # Eri
,   (  4, 20.53333333333333, 21.33333333333333, -15.00000000000000 )    # Aqr
,   ( 73, 17.16666666666667, 18.25000000000000, -16.00000000000000 )    # Ser
,   ( 72, 18.25000000000000, 18.86666666666667, -16.00000000000000 )    # Sct
,   ( 41,  8.36666666666667,  8.58333333333333, -17.00000000000000 )    # Hya
,   ( 58, 16.26666666666667, 16.37500000000000, -18.25000000000000 )    # Oph
,   ( 41,  8.58333333333333,  9.08333333333333, -19.00000000000000 )    # Hya
,   ( 26, 10.75000000000000, 10.83333333333333, -19.00000000000000 )    # Crt
,   ( 71, 16.26666666666667, 16.37500000000000, -19.25000000000000 )    # Sco
,   ( 47, 15.66666666666667, 15.91666666666667, -20.00000000000000 )    # Lib
,   ( 28, 12.58333333333333, 12.83333333333333, -22.00000000000000 )    # Crv
,   ( 85, 12.83333333333333, 14.25000000000000, -22.00000000000000 )    # Vir
,   ( 41,  9.08333333333333,  9.75000000000000, -24.00000000000000 )    # Hya
,   ( 16,  1.66666666666667,  2.65000000000000, -24.38333333333333 )    # Cet
,   ( 35,  2.65000000000000,  3.75000000000000, -24.38333333333333 )    # Eri
,   ( 26, 10.83333333333333, 11.83333333333333, -24.50000000000000 )    # Crt
,   ( 28, 11.83333333333333, 12.58333333333333, -24.50000000000000 )    # Crv
,   ( 47, 14.25000000000000, 14.91666666666667, -24.50000000000000 )    # Lib
,   ( 58, 16.26666666666667, 16.75000000000000, -24.58333333333333 )    # Oph
,   ( 16,  0.00000000000000,  1.66666666666667, -25.50000000000000 )    # Cet
,   ( 11, 21.33333333333333, 21.86666666666667, -25.50000000000000 )    # Cap
,   (  4, 21.86666666666667, 23.83333333333333, -25.50000000000000 )    # Aqr
,   ( 16, 23.83333333333333, 24.00000000000000, -25.50000000000000 )    # Cet
,   ( 41,  9.75000000000000, 10.25000000000000, -26.50000000000000 )    # Hya
,   ( 35,  4.70000000000000,  4.83333333333333, -27.25000000000000 )    # Eri
,   ( 46,  4.83333333333333,  6.11666666666667, -27.25000000000000 )    # Lep
,   ( 11, 20.00000000000000, 21.33333333333333, -28.00000000000000 )    # Cap
,   ( 41, 10.25000000000000, 10.58333333333333, -29.16666666666667 )    # Hya
,   ( 41, 12.58333333333333, 14.91666666666667, -29.50000000000000 )    # Hya
,   ( 47, 14.91666666666667, 15.66666666666667, -29.50000000000000 )    # Lib
,   ( 71, 15.66666666666667, 16.00000000000000, -29.50000000000000 )    # Sco
,   ( 35,  4.58333333333333,  4.70000000000000, -30.00000000000000 )    # Eri
,   ( 58, 16.75000000000000, 17.60000000000000, -30.00000000000000 )    # Oph
,   ( 76, 17.60000000000000, 17.83333333333333, -30.00000000000000 )    # Sgr
,   ( 41, 10.58333333333333, 10.83333333333333, -31.16666666666667 )    # Hya
,   ( 19,  6.11666666666667,  7.36666666666667, -33.00000000000000 )    # CMa
,   ( 41, 12.25000000000000, 12.58333333333333, -33.00000000000000 )    # Hya
,   ( 41, 10.83333333333333, 12.25000000000000, -35.00000000000000 )    # Hya
,   ( 36,  3.50000000000000,  3.75000000000000, -36.00000000000000 )    # For
,   ( 68,  8.36666666666667,  9.36666666666667, -36.75000000000000 )    # Pyx
,   ( 35,  4.26666666666667,  4.58333333333333, -37.00000000000000 )    # Eri
,   ( 76, 17.83333333333333, 19.16666666666667, -37.00000000000000 )    # Sgr
,   ( 65, 21.33333333333333, 23.00000000000000, -37.00000000000000 )    # PsA
,   ( 70, 23.00000000000000, 23.33333333333333, -37.00000000000000 )    # Scl
,   ( 36,  3.00000000000000,  3.50000000000000, -39.58333333333334 )    # For
,   (  1,  9.36666666666667, 11.00000000000000, -39.75000000000000 )    # Ant
,   ( 70,  0.00000000000000,  1.66666666666667, -40.00000000000000 )    # Scl
,   ( 36,  1.66666666666667,  3.00000000000000, -40.00000000000000 )    # For
,   ( 35,  3.86666666666667,  4.26666666666667, -40.00000000000000 )    # Eri
,   ( 70, 23.33333333333333, 24.00000000000000, -40.00000000000000 )    # Scl
,   ( 14, 14.16666666666667, 14.91666666666667, -42.00000000000000 )    # Cen
,   ( 49, 15.66666666666667, 16.00000000000000, -42.00000000000000 )    # Lup
,   ( 71, 16.00000000000000, 16.42083333333333, -42.00000000000000 )    # Sco
,   (  9,  4.83333333333333,  5.00000000000000, -43.00000000000000 )    # Cae
,   ( 22,  5.00000000000000,  6.58333333333333, -43.00000000000000 )    # Col
,   ( 67,  8.00000000000000,  8.36666666666667, -43.00000000000000 )    # Pup
,   ( 35,  3.41666666666667,  3.86666666666667, -44.00000000000000 )    # Eri
,   ( 71, 16.42083333333333, 17.83333333333333, -45.50000000000000 )    # Sco
,   ( 24, 17.83333333333333, 19.16666666666667, -45.50000000000000 )    # CrA
,   ( 76, 19.16666666666667, 20.33333333333333, -45.50000000000000 )    # Sgr
,   ( 53, 20.33333333333333, 21.33333333333333, -45.50000000000000 )    # Mic
,   ( 35,  3.00000000000000,  3.41666666666667, -46.00000000000000 )    # Eri
,   (  9,  4.50000000000000,  4.83333333333333, -46.50000000000000 )    # Cae
,   ( 49, 15.33333333333333, 15.66666666666667, -48.00000000000000 )    # Lup
,   ( 63,  0.00000000000000,  2.33333333333333, -48.16666666666666 )    # Phe
,   ( 35,  2.66666666666667,  3.00000000000000, -49.00000000000000 )    # Eri
,   ( 40,  4.08333333333333,  4.26666666666667, -49.00000000000000 )    # Hor
,   (  9,  4.26666666666667,  4.50000000000000, -49.00000000000000 )    # Cae
,   ( 38, 21.33333333333333, 22.00000000000000, -50.00000000000000 )    # Gru
,   ( 67,  6.00000000000000,  8.00000000000000, -50.75000000000000 )    # Pup
,   ( 84,  8.00000000000000,  8.16666666666667, -50.75000000000000 )    # Vel
,   ( 35,  2.41666666666667,  2.66666666666667, -51.00000000000000 )    # Eri
,   ( 40,  3.83333333333333,  4.08333333333333, -51.00000000000000 )    # Hor
,   ( 63,  0.00000000000000,  1.83333333333333, -51.50000000000000 )    # Phe
,   ( 12,  6.00000000000000,  6.16666666666667, -52.50000000000000 )    # Car
,   ( 84,  8.16666666666667,  8.45000000000000, -53.00000000000000 )    # Vel
,   ( 40,  3.50000000000000,  3.83333333333333, -53.16666666666666 )    # Hor
,   ( 32,  3.83333333333333,  4.00000000000000, -53.16666666666666 )    # Dor
,   ( 63,  0.00000000000000,  1.58333333333333, -53.50000000000000 )    # Phe
,   ( 35,  2.16666666666667,  2.41666666666667, -54.00000000000000 )    # Eri
,   ( 64,  4.50000000000000,  5.00000000000000, -54.00000000000000 )    # Pic
,   ( 49, 15.05000000000000, 15.33333333333333, -54.00000000000000 )    # Lup
,   ( 84,  8.45000000000000,  8.83333333333333, -54.50000000000000 )    # Vel
,   ( 12,  6.16666666666667,  6.50000000000000, -55.00000000000000 )    # Car
,   ( 14, 11.83333333333333, 12.83333333333333, -55.00000000000000 )    # Cen
,   ( 49, 14.16666666666667, 15.05000000000000, -55.00000000000000 )    # Lup
,   ( 56, 15.05000000000000, 15.33333333333333, -55.00000000000000 )    # Nor
,   ( 32,  4.00000000000000,  4.33333333333333, -56.50000000000000 )    # Dor
,   ( 84,  8.83333333333333, 11.00000000000000, -56.50000000000000 )    # Vel
,   ( 14, 11.00000000000000, 11.25000000000000, -56.50000000000000 )    # Cen
,   (  5, 17.50000000000000, 18.00000000000000, -57.00000000000000 )    # Ara
,   ( 78, 18.00000000000000, 20.33333333333333, -57.00000000000000 )    # Tel
,   ( 38, 22.00000000000000, 23.33333333333333, -57.00000000000000 )    # Gru
,   ( 40,  3.20000000000000,  3.50000000000000, -57.50000000000000 )    # Hor
,   ( 64,  5.00000000000000,  5.50000000000000, -57.50000000000000 )    # Pic
,   ( 12,  6.50000000000000,  6.83333333333333, -58.00000000000000 )    # Car
,   ( 63,  0.00000000000000,  1.33333333333333, -58.50000000000000 )    # Phe
,   ( 35,  1.33333333333333,  2.16666666666667, -58.50000000000000 )    # Eri
,   ( 63, 23.33333333333333, 24.00000000000000, -58.50000000000000 )    # Phe
,   ( 32,  4.33333333333333,  4.58333333333333, -59.00000000000000 )    # Dor
,   ( 56, 15.33333333333333, 16.42083333333333, -60.00000000000000 )    # Nor
,   ( 43, 20.33333333333333, 21.33333333333333, -60.00000000000000 )    # Ind
,   ( 64,  5.50000000000000,  6.00000000000000, -61.00000000000000 )    # Pic
,   ( 18, 15.16666666666667, 15.33333333333333, -61.00000000000000 )    # Cir
,   (  5, 16.42083333333333, 16.58333333333333, -61.00000000000000 )    # Ara
,   ( 18, 14.91666666666667, 15.16666666666667, -63.58333333333334 )    # Cir
,   (  5, 16.58333333333333, 16.75000000000000, -63.58333333333334 )    # Ara
,   ( 64,  6.00000000000000,  6.83333333333333, -64.00000000000000 )    # Pic
,   ( 12,  6.83333333333333,  9.03333333333333, -64.00000000000000 )    # Car
,   ( 14, 11.25000000000000, 11.83333333333333, -64.00000000000000 )    # Cen
,   ( 27, 11.83333333333333, 12.83333333333333, -64.00000000000000 )    # Cru
,   ( 14, 12.83333333333333, 14.53333333333333, -64.00000000000000 )    # Cen
,   ( 18, 13.50000000000000, 13.66666666666667, -65.00000000000000 )    # Cir
,   (  5, 16.75000000000000, 16.83333333333333, -65.00000000000000 )    # Ara
,   ( 40,  2.16666666666667,  3.20000000000000, -67.50000000000000 )    # Hor
,   ( 69,  3.20000000000000,  4.58333333333333, -67.50000000000000 )    # Ret
,   ( 18, 14.75000000000000, 14.91666666666667, -67.50000000000000 )    # Cir
,   (  5, 16.83333333333333, 17.50000000000000, -67.50000000000000 )    # Ara
,   ( 60, 17.50000000000000, 18.00000000000000, -67.50000000000000 )    # Pav
,   ( 81, 22.00000000000000, 23.33333333333333, -67.50000000000000 )    # Tuc
,   ( 32,  4.58333333333333,  6.58333333333333, -70.00000000000000 )    # Dor
,   ( 18, 13.66666666666667, 14.75000000000000, -70.00000000000000 )    # Cir
,   ( 79, 14.75000000000000, 17.00000000000000, -70.00000000000000 )    # TrA
,   ( 81,  0.00000000000000,  1.33333333333333, -75.00000000000000 )    # Tuc
,   ( 42,  3.50000000000000,  4.58333333333333, -75.00000000000000 )    # Hyi
,   ( 86,  6.58333333333333,  9.03333333333333, -75.00000000000000 )    # Vol
,   ( 12,  9.03333333333333, 11.25000000000000, -75.00000000000000 )    # Car
,   ( 55, 11.25000000000000, 13.66666666666667, -75.00000000000000 )    # Mus
,   ( 60, 18.00000000000000, 21.33333333333333, -75.00000000000000 )    # Pav
,   ( 43, 21.33333333333333, 23.33333333333333, -75.00000000000000 )    # Ind
,   ( 81, 23.33333333333333, 24.00000000000000, -75.00000000000000 )    # Tuc
,   ( 81,  0.75000000000000,  1.33333333333333, -76.00000000000000 )    # Tuc
,   ( 42,  0.00000000000000,  3.50000000000000, -82.50000000000000 )    # Hyi
,   ( 17,  7.66666666666667, 13.66666666666667, -82.50000000000000 )    # Cha
,   (  2, 13.66666666666667, 18.00000000000000, -82.50000000000000 )    # Aps
,   ( 52,  3.50000000000000,  7.66666666666667, -85.00000000000000 )    # Men
,   ( 57,  0.00000000000000, 24.00000000000000, -90.00000000000000 )    # Oct
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
    equ2000 = Equatorial(ra, dec, 1.0)
    vec2000 = VectorFromEquator(equ2000, _Epoch2000)
    vec1875 = RotateVector(_ConstelRot, vec2000)
    equ1875 = EquatorFromVector(vec1875)

    # Search for the constellation using the B1875 coordinates.
    for b in _ConstelBounds:
        index, ra_lo, ra_hi, dec = b
        if (dec <= equ1875.dec) and (ra_lo <= equ1875.ra < ra_hi):
            symbol, name = _ConstelNames[index]
            return ConstellationInfo(symbol, name, equ1875.ra, equ1875.dec)

    # This should never happen!
    raise Error('Unable to find constellation for given coordinates.')


@enum.unique
class EclipseKind(enum.Enum):
    """Selects if/how to correct for atmospheric refraction.

    Some functions allow enabling or disabling atmospheric refraction
    for the calculated apparent position of a celestial body
    as seen by an observer on the surface of the Earth.

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
    def __init__(self, time, u, r, k, p, target, dir):
        self.time = time
        self.u = u   # dot product of (heliocentric earth) and (geocentric moon): defines the shadow plane where the Moon is
        self.r = r   # km distance between center of Moon and the line passing through the centers of the Sun and Earth.
        self.k = k   # umbra radius in km, at the shadow plane
        self.p = p   # penumbra radius in km, at the shadow plane
        self.target = target        # vector from center of shadow-casting body to target location that might receive the shadow
        self.dir = dir              # vector from center of Sun to center of shadow-casting body


def _CalcShadow(body_radius_km, time, target, dir):
    u = (dir.x*target.x + dir.y*target.y + dir.z*target.z) / (dir.x*dir.x + dir.y*dir.y + dir.z*dir.z)
    dx = (u * dir.x) - target.x
    dy = (u * dir.y) - target.y
    dz = (u * dir.z) - target.z
    r = _KM_PER_AU * math.sqrt(dx*dx + dy*dy + dz*dz)
    k = +_SUN_RADIUS_KM - (1.0 + u)*(_SUN_RADIUS_KM - body_radius_km)
    p = -_SUN_RADIUS_KM + (1.0 + u)*(_SUN_RADIUS_KM + body_radius_km)
    return _ShadowInfo(time, u, r, k, p, target, dir)


def _EarthShadow(time):
    e = _CalcEarth(time)
    m = GeoMoon(time)
    return _CalcShadow(_EARTH_ECLIPSE_RADIUS_KM, time, m, e)


def _MoonShadow(time):
    # This is a variation on the logic in _EarthShadow().
    # Instead of a heliocentric Earth and a geocentric Moon,
    # we want a heliocentric Moon and a lunacentric Earth.
    h = _CalcEarth(time)    # heliocentric Earth
    m = GeoMoon(time)       # geocentric Moon
    # Calculate lunacentric Earth.
    e = Vector(-m.x, -m.y, -m.z, m.t)
    # Convert geocentric moon to heliocentric Moon.
    m.x += h.x
    m.y += h.y
    m.z += h.z
    return _CalcShadow(_MOON_MEAN_RADIUS_KM, time, e, m)


def _LocalMoonShadow(time, observer):
    # Calculate observer's geocentric position.
    # For efficiency, do this first, to populate the earth rotation parameters in 'time'.
    # That way they can be recycled instead of recalculated.
    pos = _geo_pos(time, observer)
    h = _CalcEarth(time)     # heliocentric Earth
    m = GeoMoon(time)        # geocentric Moon
    # Calculate lunacentric location of an observer on the Earth's surface.
    o = Vector(pos[0] - m.x, pos[1] - m.y, pos[2] - m.z, time)
    # Convert geocentric moon to heliocentric Moon.
    m.x += h.x
    m.y += h.y
    m.z += h.z
    return _CalcShadow(_MOON_MEAN_RADIUS_KM, time, o, m)


def _PlanetShadow(body, planet_radius_km, time):
    # Calculate light-travel-corrected vector from Earth to planet.
    g = GeoVector(body, time, False)
    # Calculate light-travel-corrected vector from Earth to Sun.
    e = GeoVector(Body.Sun, time, False)
    # Deduce light-travel-corrected vector from Sun to planet.
    p = Vector(g.x - e.x, g.y - e.y, g.z - e.z, time)
    # Calcluate Earth's position from the planet's point of view.
    e.x = -g.x
    e.y = -g.y
    e.z = -g.z
    return _CalcShadow(planet_radius_km, time, e, p)


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
    Penumbral eclipses are difficult to observe, because the moon is only slightly dimmed
    by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
    Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
    Total eclipses occur when the entire Moon passes into the Earth's umbra.

    The `kind` field thus holds one of the values `EclipseKind.Penumbral`, `EclipseKind.Partial`,
    or `EclipseKind.Total`, depending on the kind of lunar eclipse found.

    Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.

    Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
    of the eclipse, which is half of the amount of time the eclipse spends in each
    phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
    By converting from minutes to days, and subtracting/adding with `peak`, the caller
    may determine the date and time of the beginning/end of each eclipse phase.

    Attributes
    ----------
    kind : string
         The type of lunar eclipse found.
    peak : Time
         The time of the eclipse at its peak.
    sd_penum : float
         The semi-duration of the penumbral phase in minutes.
    sd_partial : float
         The semi-duration of the penumbral phase in minutes, or 0.0 if none.
    sd_total : float
         The semi-duration of the penumbral phase in minutes, or 0.0 if none.
    """
    def __init__(self, kind, peak, sd_penum, sd_partial, sd_total):
        self.kind = kind
        self.peak = peak
        self.sd_penum = sd_penum
        self.sd_partial = sd_partial
        self.sd_total = sd_total


class GlobalSolarEclipseInfo:
    """Reports the time and geographic location of the peak of a solar eclipse.

    Returned by #SearchGlobalSolarEclipse or #NextGlobalSolarEclipse
    to report information about a solar eclipse event.

    Field `peak` holds the date and time of the peak of the eclipse, defined as
    the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.

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

    Attributes
    ----------
    kind : EclipseKind
        The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    peak : Time
        The date and time of the eclipse at its peak.
    distance : float
        The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.
    latitude : float
        The geographic latitude at the center of the peak eclipse shadow.
    longitude : float
        The geographic longitude at the center of the peak eclipse shadow.
    """
    def __init__(self, kind, peak, distance, latitude, longitude):
        self.kind = kind
        self.peak = peak
        self.distance = distance
        self.latitude = latitude
        self.longitude = longitude


class EclipseEvent:
    """Holds a time and the observed altitude of the Sun at that time.

    When reporting a solar eclipse observed at a specific location on the Earth
    (a "local" solar eclipse), a series of events occur. In addition
    to the time of each event, it is important to know the altitude of the Sun,
    because each event may be invisible to the observer if the Sun is below
    the horizon (i.e. it at night).

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
    def __init__(self, kind, partial_begin, total_begin, peak, total_end, partial_end):
        self.kind = kind
        self.partial_begin = partial_begin
        self.total_begin = total_begin
        self.peak = peak
        self.total_end = total_end
        self.partial_end = partial_end


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
    return LocalSolarEclipseInfo(kind, partial_begin, total_begin, peak, total_end, partial_end)


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
    v.x *= _KM_PER_AU
    v.y *= _KM_PER_AU
    v.z *= _KM_PER_AU / _EARTH_FLATTENING
    e.x *= _KM_PER_AU
    e.y *= _KM_PER_AU
    e.z *= _KM_PER_AU / _EARTH_FLATTENING

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
        proj = math.sqrt(px*px + py*py) * (_EARTH_FLATTENING * _EARTH_FLATTENING)
        if proj == 0.0:
            if pz > 0.0:
                latitude = +90.0
            else:
                latitude = -90.0
        else:
            latitude = math.degrees(math.atan(pz / proj))

        # Adjust longitude for Earth's rotation at the given UT.
        gast = _sidereal_time(peak)
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
        o = Vector(px / _KM_PER_AU, py / _KM_PER_AU, pz / _KM_PER_AU, shadow.time)

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

    return GlobalSolarEclipseInfo(kind, peak, distance, latitude, longitude)


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
    for fmcount in range(12):
        # Search for the next full moon. Any eclipse will be near it.
        fullmoon = SearchMoonPhase(180, fmtime, 40)
        if fullmoon is None:
            raise Error('Cannot find full moon.')

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
                        sd_total = _ShadowSemiDurationMinutes(shadow.time, shadow.k - _MOON_MEAN_RADIUS_KM, sd_partial)

                return LunarEclipseInfo(kind, shadow.time, sd_penum, sd_partial, sd_total)

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
    for nmcount in range(12):
        # Search for the next new moon. Any eclipse will be near it.
        newmoon = SearchMoonPhase(0.0, nmtime, 40.0)
        if newmoon is None:
            raise Error('Cannot find new moon')

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
