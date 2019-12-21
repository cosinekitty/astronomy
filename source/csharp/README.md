# Astronomy Engine (C#)

This is the complete programming reference for the C# version of Astronomy Engine.
Other programming languages are supported.
See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start
To get started quickly, here are some [examples](../../demo/csharp/).

---

## Contents

- [Topic Index](#topics)
- [Functions](#functions)
- [Enumerated Types](#enums)
- [Structures](#structs)
- [Type Definitions](#typedefs)

---

<a name="topics"></a>
## Topic Index

### Dates and times

| Function | Description |
| -------- | ----------- |
| [CurrentTime](#Astronomy.CurrentTime) | Obtains the current date and time of the computer's clock in the form of an [`astro_time_t`](#astro_time_t) that can be used for astronomy calculations. |
| [MakeTime](#Astronomy.MakeTime) | Converts a UTC calendar date and time given as separate numeric parameters into an [`astro_time_t`](#astro_time_t) that can be used for astronomy calculations. |
| [AddDays](#Astronomy.AddDays) | Adds or subtracts an amount of time to an [`astro_time_t`](#astro_time_t) to get another [`astro_time_t`](#astro_time_t). |
| [TimeFromUtc](#Astronomy.TimeFromUtc) | Converts UTC calendar date and time from an [`astro_utc_t`](#astro_utc_t) structure to an [`astro_time_t`](#astro_time_t) structure that can be used for astronomy calculations. |
| [UtcFromTime](#Astronomy.UtcFromTime) | Converts an astronomical [`astro_time_t`](#astro_time_t) time value to an [`astro_utc_t`](#astro_utc_t) structure that can be used for displaying a UTC calendar date and time. |

### Celestial bodies

| Function | Description |
| -------- | ----------- |
| [BodyCode](#Astronomy.BodyCode) | Converts the English name of a celestial body to its equivalent [`astro_body_t`](#astro_body_t) enumeration value. |
| [BodyName](#Astronomy.BodyName) | Converts an [`astro_body_t`](#astro_body_t) enumeration value to its equivalent English name as a string. |

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#Astronomy.HelioVector) | Calculates vector with respect to the center of the Sun. |
| [GeoVector](#Astronomy.GeoVector)     | Calculates vector with respect to the center of the Earth. |
| [Equator](#Astronomy.Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Astronomy.Ecliptic)       | Converts J2000 equatorial coordinates to J2000 ecliptic coordinates. |
| [EclipticLongitude](#Astronomy.EclipticLongitude) | Calculates ecliptic longitude of a body in the J2000 system. |
| [Horizon](#Astronomy.Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |
| [LongitudeFromSun](#Astronomy.LongitudeFromSun) | Calculates a body's apparent ecliptic longitude difference from the Sun, as seen by an observer on the Earth. |

### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSet](#Astronomy.SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchHourAngle](#Astronomy.SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#Astronomy.MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonPhase](#Astronomy.SearchMoonPhase) | Finds the next instance of the Moon reaching a specific ecliptic longitude separation from the Sun. |
| [SearchMoonQuarter](#Astronomy.SearchMoonQuarter) | Finds the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#Astronomy.NextMoonQuarter) | Finds the next quarter moon phase after a previous one that has been found. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#Astronomy.SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#Astronomy.NextLunarApsis) | Given an already-found apsis, finds the next perigee or apogee of the Moon. |

### Visual magnitude and elongation

| Function | Description |
| -------- | ----------- |
| [Illumination](#Astronomy.Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |
| [SearchPeakMagnitude](#Astronomy.SearchPeakMagnitude) | Searches for the date and time Venus will next appear brightest as seen from the Earth. |
| [AngleFromSun](#Astronomy.AngleFromSun) | Returns full angle seen from Earth between body and Sun. |
| [Elongation](#Astronomy.Elongation) | Calculates ecliptic longitude angle between a body and the Sun, as seen from the Earth. |
| [SearchMaxElongation](#Astronomy.SearchMaxElongation) | Searches for the next maximum elongation event for Mercury or Venus that occurs after the given date. |

### Oppositions and conjunctions

| Function | Description |
| -------- | ----------- |
| [SearchRelativeLongitude](#Astronomy.SearchRelativeLongitude) | Finds oppositions and conjunctions of planets. |

### Equinoxes, solstices, and apparent solar motion

| Function | Description |
| -------- | ----------- |
| [SearchSunLongitude](#Astronomy.SearchSunLongitude) | Finds the next time the Sun reaches a specified apparent ecliptic longitude in the *true equator of date* system. |
| [Seasons](#Astronomy.Seasons) | Finds the equinoxes and solstices for a given calendar year. |
| [SunPosition](#Astronomy.SunPosition) | Calculates the Sun's apparent ecliptic coordinates as seen from the Earth. |

### Coordinate transforms

The following four orientation systems are supported.
Astronomy Engine can convert a vector from any of these orientations to any of the others.
It also allows converting from a vector to spherical (angular) coordinates and back,
within a given orientation. Note the 3-letter codes for each of the orientation systems;
these are used in function and type names.

- **EQJ = Equatorial J2000**: Uses the Earth's equator on January 1, 2000, at noon UTC.
- **EQD = Equator of-date**: Uses the Earth's equator on a given date and time, adjusted for precession and nutation.
- **ECL = Ecliptic**: Uses the mean plane of the Earth's orbit around the Sun. The x-axis is referenced against the J2000 equinox.
- **HOR = Horizontal**: Uses the viewpoint of an observer at a specific location on the Earth at a given date and time.

| Function | Description |
| -------- | ----------- |
| [RotateVector](#Astronomy.RotateVector) | Applies a rotation matrix to a vector, yielding a vector in another orientation system. |
| [InverseRotation](#Astronomy.InverseRotation) | Given a rotation matrix, finds the inverse rotation matrix that does the opposite transformation. |
| [CombineRotation](#Astronomy.CombineRotation) | Given two rotation matrices, returns a rotation matrix that combines them into a net transformation. |
| [VectorFromSphere](#Astronomy.VectorFromSphere) | Converts spherical coordinates to Cartesian coordinates. |
| [SphereFromVector](#Astronomy.SphereFromVector) | Converts Cartesian coordinates to spherical coordinates. |
| [VectorFromEquator](#Astronomy.VectorFromEquator) | Given angular equatorial coordinates, calculates equatorial vector. |
| [EquatorFromVector](#Astronomy.EquatorFromVector) | Given an equatorial vector, calculates equatorial angular coordinates. |
| [VectorFromHorizon](#Astronomy.VectorFromHorizon) | Given apparent angular horizontal coordinates, calculates horizontal vector. |
| [HorizonFromVector](#Astronomy.HorizonFromVector) | Given a vector in horizontal orientation, calculates horizontal angular coordinates. |
| [Rotation_EQD_EQJ](#Astronomy.Rotation_EQD_EQJ) | Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ). |
| [Rotation_EQD_ECL](#Astronomy.Rotation_EQD_ECL) | Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL). |
| [Rotation_EQD_HOR](#Astronomy.Rotation_EQD_HOR) | Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR). |
| [Rotation_EQJ_EQD](#Astronomy.Rotation_EQJ_EQD) | Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD). |
| [Rotation_EQJ_ECL](#Astronomy.Rotation_EQJ_ECL) | Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL). |
| [Rotation_EQJ_HOR](#Astronomy.Rotation_EQJ_HOR) | Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR). |
| [Rotation_ECL_EQD](#Astronomy.Rotation_ECL_EQD) | Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD). |
| [Rotation_ECL_EQJ](#Astronomy.Rotation_ECL_EQJ) | Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ). |
| [Rotation_ECL_HOR](#Astronomy.Rotation_ECL_HOR) | Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR). |
| [Rotation_HOR_EQD](#Astronomy.Rotation_HOR_EQD) | Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD). |
| [Rotation_HOR_EQJ](#Astronomy.Rotation_HOR_EQJ) | Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ). |
| [Rotation_HOR_ECL](#Astronomy.Rotation_HOR_ECL) | Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL). |

---

<a name="functions"></a>
## Functions

---

<a name="Astronomy.AngleFromSun"></a>
### Astronomy.AngleFromSun(body, time) &#8658; double

**Returns the angle between the given body and the Sun, as seen from the Earth.**

This function calculates the angular separation between the given body and the Sun,
as seen from the center of the Earth. This angle is helpful for determining how
easy it is to see the body away from the glare of the Sun.

<a name="Astronomy.EclipticLongitude"></a>
### Astronomy.EclipticLongitude(body, time) &#8658; double

**Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.**

This function calculates the angle around the plane of the Earth's orbit
of a celestial body, as seen from the center of the Sun.
The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).

<a name="Astronomy.Elongation"></a>
### Astronomy.Elongation(body, time) &#8658; [`ElongationInfo`](#ElongationInfo)

**Determines visibility of a celestial body relative to the Sun, as seen from the Earth.**

This function returns an #ElongationInfo structure, which provides the following
information about the given celestial body at the given time:

- `visibility` is an enumerated type that specifies whether the body is more easily seen
   in the morning before sunrise, or in the evening after sunset.

- `elongation` is the angle in degrees between two vectors: one from the center of the Earth to the
   center of the Sun, the other from the center of the Earth to the center of the specified body.
   This angle indicates how far away the body is from the glare of the Sun.
   The elongation angle is always in the range [0, 180].

- `ecliptic_separation` is the absolute value of the difference between the body's ecliptic longitude
  and the Sun's ecliptic longitude, both as seen from the center of the Earth. This angle measures
  around the plane of the Earth's orbit, and ignores how far above or below that plane the body is.
  The ecliptic separation is measured in degrees and is always in the range [0, 180].

<a name="Astronomy.Equator"></a>
### Astronomy.Equator(body, time, observer, equdate, aberration) &#8658; [`Equatorial`](#Equatorial)

**Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.**

Calculates topocentric equatorial coordinates in one of two different systems:
J2000 or true-equator-of-date, depending on the value of the `equdate` parameter.
Equatorial coordinates include right ascension, declination, and distance in astronomical units.

This function corrects for light travel time: it adjusts the apparent location
of the observed body based on how long it takes for light to travel from the body to the Earth.

This function corrects for *topocentric parallax*, meaning that it adjusts for the
angular shift depending on where the observer is located on the Earth. This is most
significant for the Moon, because it is so close to the Earth. However, parallax corection
has a small effect on the apparent positions of other bodies.

Correction for aberration is optional, using the `aberration` parameter.

<a name="Astronomy.EquatorialToEcliptic"></a>
### Astronomy.EquatorialToEcliptic(equ) &#8658; [`Ecliptic`](#Ecliptic)

**Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.**

Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
which are relative to the plane of the Earth's orbit around the Sun.

<a name="Astronomy.GeoVector"></a>
### Astronomy.GeoVector(body, time, aberration) &#8658; [`AstroVector`](#AstroVector)

**Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Earth as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.

If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
the year range 1700..2200, this function will throw an exception.

Unlike #HelioVector, this function always corrects for light travel time.
This means the position of the body is "back-dated" by the amount of time it takes
light to travel from that body to an observer on the Earth.

Also, the position can optionally be corrected for
[aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
causing the apparent direction of the body to be shifted due to transverse
movement of the Earth with respect to the rays of light coming from that body.

<a name="Astronomy.HelioVector"></a>
### Astronomy.HelioVector(body, time) &#8658; [`AstroVector`](#AstroVector)

**Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Sun as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.

The position is not corrected for light travel time or aberration.
This is different from the behavior of #GeoVector.

If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
the year range 1700..2200, this function will throw an `ArgumentException`.

<a name="Astronomy.Horizon"></a>
### Astronomy.Horizon(time, observer, ra, dec, refraction) &#8658; [`Topocentric`](#Topocentric)

**Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.**

Given a date and time, the geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial body,
this function returns horizontal coordinates (azimuth and altitude angles) for the body
relative to the horizon at the geographic location.

The right ascension `ra` and declination `dec` passed in must be *equator of date*
coordinates, based on the Earth's true equator at the date and time of the observation.
Otherwise the resulting horizontal coordinates will be inaccurate.
Equator of date coordinates can be obtained by calling #Equator, passing in
`EquatorEpoch.OfDate` as its `equdate` parameter. It is also recommended to enable
aberration correction by passing in `Aberration.Corrected` as the `aberration` parameter.

This function optionally corrects for atmospheric refraction.
For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
correct for optical lensing of the Earth's atmosphere that causes objects
to appear somewhat higher above the horizon than they actually are.
However, callers may choose to avoid this correction by passing in `Refraction.None`.
If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
in the #Topocentric structure returned by this function will all be corrected for refraction.
If refraction is disabled, none of these four coordinates will be corrected; in that case,
the right ascension and declination in the returned structure will be numerically identical
to the respective `ra` and `dec` values passed in.

<a name="Astronomy.Illumination"></a>
### Astronomy.Illumination(body, time) &#8658; [`IllumInfo`](#IllumInfo)

**Finds visual magnitude, phase angle, and other illumination information about a celestial body.**

This function calculates information about how bright a celestial body appears from the Earth,
reported as visual magnitude, which is a smaller (or even negative) number for brighter objects
and a larger number for dimmer objects.

For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
of the body appears illuminated as seen from the Earth. For example, when the phase angle is
near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
of 90 degrees means the body appears "half full".
For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
so it doesn't have a phase angle.

When the body is Saturn, the returned structure contains a field `ring_tilt` that holds
the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
0 for all bodies other than Saturn.

<a name="Astronomy.LongitudeFromSun"></a>
### Astronomy.LongitudeFromSun(body, time) &#8658; double

**Returns a body's ecliptic longitude with respect to the Sun, as seen from the Earth.**

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

<a name="Astronomy.MoonPhase"></a>
### Astronomy.MoonPhase(time) &#8658; double

**Returns the Moon's phase as an angle from 0 to 360 degrees.**

This function determines the phase of the Moon using its apparent
ecliptic longitude relative to the Sun, as seen from the center of the Earth.
Certain values of the angle have conventional definitions:

- 0 = new moon
- 90 = first quarter
- 180 = full moon
- 270 = third quarter

<a name="Astronomy.NextLunarApsis"></a>
### Astronomy.NextLunarApsis(apsis) &#8658; [`ApsisInfo`](#ApsisInfo)

**Finds the next lunar perigee or apogee event in a series.**

This function requires an #ApsisInfo value obtained from a call
to #Astronomy.SearchLunarApsis or `Astronomy.NextLunarApsis`. Given
an apogee event, this function finds the next perigee event, and vice versa.

See #Astronomy_SearchLunarApsis for more details.

<a name="Astronomy.NextMoonQuarter"></a>
### Astronomy.NextMoonQuarter(mq) &#8658; [`MoonQuarterInfo`](#MoonQuarterInfo)

**Continues searching for lunar quarters from a previous search.**

After calling #SearchMoonQuarter, this function can be called
one or more times to continue finding consecutive lunar quarters.
This function finds the next consecutive moon quarter event after
the one passed in as the parameter `mq`.

<a name="Astronomy.Search"></a>
### Astronomy.Search(func, t1, t2, dt_tolerance_seconds) &#8658; [`AstroTime`](#AstroTime)

**Searches for a time at which a function's value increases through zero.**

Certain astronomy calculations involve finding a time when an event occurs.
Often such events can be defined as the root of a function:
the time at which the function's value becomes zero.

`Search` finds the *ascending root* of a function: the time at which
the function's value becomes zero while having a positive slope. That is, as time increases,
the function transitions from a negative value, through zero at a specific moment,
to a positive value later. The goal of the search is to find that specific moment.

The `func` parameter is an instance of the abstract class #SearchContext.
As an example, a caller may wish to find the moment a celestial body reaches a certain
ecliptic longitude. In that case, the caller might derive a class that contains
a #Body member to specify the body and a `double` to hold the target longitude.
It could subtract the target longitude from the actual longitude at a given time;
thus the difference would equal zero at the moment in time the planet reaches the
desired longitude.

Every call to `func.Eval` must either return a valid #AstroTime or throw an exception.

The search calls `func.Eval` repeatedly to rapidly narrow in on any ascending
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
the search will return `null`.

If the search does not converge within 20 iterations, it will throw an exception.

<a name="Astronomy.SearchHourAngle"></a>
### Astronomy.SearchHourAngle(body, observer, hourAngle, startTime) &#8658; [`HourAngleInfo`](#HourAngleInfo)

**Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.**

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

<a name="Astronomy.SearchLunarApsis"></a>
### Astronomy.SearchLunarApsis(startTime) &#8658; [`ApsisInfo`](#ApsisInfo)

**Finds the date and time of the Moon's closest distance (perigee)
or farthest distance (apogee) with respect to the Earth.**

Given a date and time to start the search in `startTime`, this function finds the
next date and time that the center of the Moon reaches the closest or farthest point
in its orbit with respect to the center of the Earth, whichever comes first
after `startTime`.

The closest point is called *perigee* and the farthest point is called *apogee*.
The word *apsis* refers to either event.

To iterate through consecutive alternating perigee and apogee events, call `Astronomy.SearchLunarApsis`
once, then use the return value to call #Astronomy.NextLunarApsis. After that,
keep feeding the previous return value from `Astronomy_NextLunarApsis` into another
call of `Astronomy_NextLunarApsis` as many times as desired.

<a name="Astronomy.SearchMaxElongation"></a>
### Astronomy.SearchMaxElongation(body, startTime) &#8658; [`ElongationInfo`](#ElongationInfo)

**Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.**

Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
a telescope without atmospheric interference, are when these planets reach maximum elongation.
These are events where the planets reach the maximum angle from the Sun as seen from the Earth.

This function solves for those times, reporting the next maximum elongation event's date and time,
the elongation value itself, the relative longitude with the Sun, and whether the planet is best
observed in the morning or evening. See #Astronomy_Elongation for more details about the returned structure.

<a name="Astronomy.SearchMoonPhase"></a>
### Astronomy.SearchMoonPhase(targetLon, startTime, limitDays) &#8658; [`AstroTime`](#AstroTime)

**Searches for the time that the Moon reaches a specified phase.**

Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
longitude with respect to the Sun's geocentric ecliptic longitude.
When the Moon and the Sun have the same longitude, that is defined as a new moon.
When their longitudes are 180 degrees apart, that is defined as a full moon.

This function searches for any value of the lunar phase expressed as an
angle in degrees in the range [0, 360).

If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
it is much easier to call the functions #SearchMoonQuarter and #NextMoonQuarter.
This function is useful for finding general phase angles outside those four quarters.

<a name="Astronomy.SearchMoonQuarter"></a>
### Astronomy.SearchMoonQuarter(startTime) &#8658; [`MoonQuarterInfo`](#MoonQuarterInfo)

**Finds the first lunar quarter after the specified date and time.**

A lunar quarter is one of the following four lunar phase events:
new moon, first quarter, full moon, third quarter.
This function finds the lunar quarter that happens soonest
after the specified date and time.

To continue iterating through consecutive lunar quarters, call this function once,
followed by calls to #NextMoonQuarter as many times as desired.

<a name="Astronomy.SearchPeakMagnitude"></a>
### Astronomy.SearchPeakMagnitude(body, startTime) &#8658; [`IllumInfo`](#IllumInfo)

**Searches for the date and time Venus will next appear brightest as seen from the Earth.**

This function searches for the date and time Venus appears brightest as seen from the Earth.
Currently only Venus is supported for the `body` parameter, though this could change in the future.
Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from the Earth,
so peak magnitude events have little practical value for that planet.
Planets other than Venus and Mercury reach peak magnitude at opposition, which can
be found using #Astronomy.SearchRelativeLongitude.
The Moon reaches peak magnitude at full moon, which can be found using
#Astronomy.SearchMoonQuarter or #Astronomy.SearchMoonPhase.
The Sun reaches peak magnitude at perihelion, which occurs each year in January.
However, the difference is minor and has little practical value.

<a name="Astronomy.SearchRelativeLongitude"></a>
### Astronomy.SearchRelativeLongitude(body, targetRelLon, startTime) &#8658; [`AstroTime`](#AstroTime)

**Searches for the time when the Earth and another planet are separated by a specified angle
in ecliptic longitude, as seen from the Sun.**

A relative longitude is the angle between two bodies measured in the plane of the Earth's orbit
(the ecliptic plane). The distance of the bodies above or below the ecliptic plane is ignored.
If you imagine the shadow of the body cast onto the ecliptic plane, and the angle measured around
that plane from one body to the other in the direction the planets orbit the Sun, you will get an
angle somewhere between 0 and 360 degrees. This is the relative longitude.

Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
this function searches for the next time that the relative longitude measured from the planet
to the Earth is `targetRelLon`.

Certain astronomical events are defined in terms of relative longitude between the Earth and another planet:

- When the relative longitude is 0 degrees, it means both planets are in the same direction from the Sun.
  For planets that orbit closer to the Sun (Mercury and Venus), this is known as *inferior conjunction*,
  a time when the other planet becomes very difficult to see because of being lost in the Sun's glare.
  (The only exception is in the rare event of a transit, when we see the silhouette of the planet passing
  between the Earth and the Sun.)

- When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
  this is known as *opposition*.  Opposition is when the planet is closest to the Earth, and
  also when it is visible for most of the night, so it is considered the best time to observe the planet.

- When the relative longitude is 180 degrees, it means the other planet is on the opposite side of the Sun
  from the Earth. This is called *superior conjunction*. Like inferior conjunction, the planet is
  very difficult to see from the Earth. Superior conjunction is possible for any planet other than the Earth.

<a name="Astronomy.SearchRiseSet"></a>
### Astronomy.SearchRiseSet(body, observer, direction, startTime, limitDays) &#8658; [`AstroTime`](#AstroTime)

**Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.**

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

<a name="Astronomy.SearchSunLongitude"></a>
### Astronomy.SearchSunLongitude(targetLon, startTime, limitDays) &#8658; [`AstroTime`](#AstroTime)

**Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.**

This function finds the moment in time, if any exists in the given time window,
that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.

This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call #Seasons
to calculate all equinoxes and solstices for a given calendar year.

The function searches the window of time specified by `startTime` and `startTime+limitDays`.
The search will return an error if the Sun never reaches the longitude `targetLon` or
if the window is so large that the longitude ranges more than 180 degrees within it.
It is recommended to keep the window smaller than 10 days when possible.

<a name="Astronomy.Seasons"></a>
### Astronomy.Seasons(year) &#8658; [`SeasonsInfo`](#SeasonsInfo)

**Finds both equinoxes and both solstices for a given calendar year.**

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

<a name="Astronomy.SunPosition"></a>
### Astronomy.SunPosition(time) &#8658; [`Ecliptic`](#Ecliptic)

**Calculates geocentric ecliptic coordinates for the Sun.**

This function calculates the position of the Sun as seen from the Earth.
The returned value includes both Cartesian and spherical coordinates.
The x-coordinate and longitude values in the returned structure are based
on the *true equinox of date*: one of two points in the sky where the instantaneous
plane of the Earth's equator at the given date and time (the *equatorial plane*)
intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
By convention, the apparent location of the Sun at the March equinox is chosen
as the longitude origin and x-axis direction, instead of the one for September.

`SunPosition` corrects for precession and nutation of the Earth's axis
in order to obtain the exact equatorial plane at the given time.

This function can be used for calculating changes of seasons: equinoxes and solstices.
In fact, the function #Seasons does use this function for that purpose.

