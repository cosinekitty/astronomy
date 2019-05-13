# Astronomy Engine (JavaScript)

This is the complete programming reference for the JavaScript version
of the [Astronomy Engine](../../). It supports client side programming
in the browser and backend use of [Node.js](https://nodejs.org).
Other programming languages are supported also. See the [home page](../../) for more info.

See here for [examples](../../demo/js) to get started quickly.

# Usage Guide

Position vectors of Sun, Moon, and planets:

| [HelioVector](#Astronomy.HelioVector) | Calculates vector with respect to the center of the Sun.   |
| [GeoVector](#Astronomy.GeoVector)     | Calculates vector with respect to the center of the Earth. |

Visual magnitude:

| [Illumination](#Astronomy.Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |

# API Reference
<a name="Astronomy"></a>

## Astronomy : <code>object</code>
**Kind**: global namespace  

* * *

<a name="Astronomy.PerformanceInfo"></a>

### Astronomy.PerformanceInfo
Holds performance metrics for developers to optimize execution speed.
Most users can safely ignore this class.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| search_func | <code>number</code> | Number of times [Search](#Astronomy.Search) called a <code>func</code> passed to it. |
| search | <code>number</code> | Number of times [Search](#Astronomy.Search) was called. |
| longitude_search | <code>number</code> | Number of times [SearchRelativeLongitude](#Astronomy.SearchRelativeLongitude) was called. |
| longitude_iter | <code>number</code> | The total number of iterations executed inside [SearchRelativeLongitude](#Astronomy.SearchRelativeLongitude). |
| lunar_apsis_calls | <code>number</code> | The number of times [SearchLunarApsis](#Astronomy.SearchLunarApsis) was called. |
| lunar_apsis_iter | <code>number</code> | The number of search iterations inside [SearchLunarApsis](#Astronomy.SearchLunarApsis). |
| calcmoon | <code>number</code> | The number of times the Moon's position was calculated. (This is an expensive operation.) |


* * *

<a name="Astronomy.PerformanceInfo+Clone"></a>

#### performanceInfo.Clone() ⇒ [<code>PerformanceInfo</code>](#Astronomy.PerformanceInfo)
Creates a copy of a <code>PerformanceInfo</code> object.
This allows us to create a snapshot of the performance metrics
that can be handed back to outside code that will not change
as the Astronomy code continues to execute and change the metrics.

**Kind**: instance method of [<code>PerformanceInfo</code>](#Astronomy.PerformanceInfo)  

* * *

<a name="Astronomy.Time"></a>

### Astronomy.Time
The date and time of an astronomical observation.
Objects of this type are used throughout the internals
of the Astronomy library, and are included in certain return objects.
The constructor is not accessible outside the Astronomy library;
outside users should call the [MakeTime](#Astronomy.MakeTime) function
to create a <code>Time</code> object.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> | The JavaScript Date object for the given date and time.      This Date corresponds to the numeric day value stored in the <code>ut</code> property. |
| ut | <code>number</code> | Universal Time (UT1/UTC) in fractional days since the J2000 epoch.      Universal Time represents time measured with respect to the Earth's rotation,      tracking mean solar days.      The Astronomy library approximates UT1 and UTC as being the same thing.      This gives sufficient accuracy for the precision requirements of this project. |
| tt | <code>number</code> | Terrestrial Time in fractional days since the J2000 epoch.      TT represents a continuously flowing ephemeris timescale independent of      any variations of the Earth's rotation, and is adjusted from UT      using historical and predictive models of those variations. |


* [.Time](#Astronomy.Time)
    * [new Time(date)](#new_Astronomy.Time_new)
    * [.toString()](#Astronomy.Time+toString) ⇒ <code>string</code>
    * [.AddDays(days)](#Astronomy.Time+AddDays) ⇒ [<code>Time</code>](#Astronomy.Time)


* * *

<a name="new_Astronomy.Time_new"></a>

#### new Time(date)

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> | A JavaScript Date object or a numeric UTC value expressed in J2000 days. |


* * *

<a name="Astronomy.Time+toString"></a>

#### time.toString() ⇒ <code>string</code>
Formats a <code>Time</code> object as an 
<a href="https://en.wikipedia.org/wiki/ISO_8601">ISO 8601</a>
date/time string in UTC, to millisecond resolution.
Example: 
<pre>
<code>2018-08-17T17:22:04.050Z</code>
</pre>

**Kind**: instance method of [<code>Time</code>](#Astronomy.Time)  

* * *

<a name="Astronomy.Time+AddDays"></a>

#### time.AddDays(days) ⇒ [<code>Time</code>](#Astronomy.Time)
Returns a new <code>Time</code> object adjusted by the floating point number of days.
Does NOT modify the original <code>Time</code> object.

**Kind**: instance method of [<code>Time</code>](#Astronomy.Time)  

| Param | Type | Description |
| --- | --- | --- |
| days | <code>number</code> | The floating point number of days by which to adjust the given date and time.      Positive values adjust the date toward the future, and      negative values adjust the date toward the past. |


* * *

<a name="Astronomy.Vector"></a>

### Astronomy.Vector
Holds the Cartesian coordinates of a vector in 3D space,
along with the time at which the vector is valid.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| x | <code>number</code> | The x-coordinate expressed in astronomical units (AU). |
| y | <code>number</code> | The y-coordinate expressed in astronomical units (AU). |
| z | <code>number</code> | The z-coordinate expressed in astronomical units (AU). |
| t | [<code>Time</code>](#Astronomy.Time) | The time at which the vector is valid. |


* * *

<a name="Astronomy.Vector+Length"></a>

#### vector.Length() ⇒ <code>number</code>
Returns the length of the vector in astronomical units (AU).

**Kind**: instance method of [<code>Vector</code>](#Astronomy.Vector)  

* * *

<a name="Astronomy.EquatorialCoordinates"></a>

### Astronomy.EquatorialCoordinates
Holds right ascension, declination, and distance of a celestial object.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| ra | <code>number</code> | Right ascension in sidereal hours: [0, 24). |
| dec | <code>number</code> | Declination in degrees: [-90, +90]. |
| dist | <code>number</code> | Distance to the celestial object expressed in       <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU). |


* * *

<a name="Astronomy.SkyCoordinates"></a>

### Astronomy.SkyCoordinates
Holds topocentric equatorial coordinates (right ascension and declination)
simultaneously in two different systems: J2000 and true-equator-of-date.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| t | [<code>Time</code>](#Astronomy.Time) | The date and time at which the coordinates are valid. |
| j2000 | [<code>EquatorialCoordinates</code>](#Astronomy.EquatorialCoordinates) | Equatorial coordinates referenced to the J2000 coordinate system. |
| ofdate | [<code>EquatorialCoordinates</code>](#Astronomy.EquatorialCoordinates) | Equatorial coordinates referenced to the true equator and equinox      at the specified date and time stored in <code>t</code>.      These coordinates are corrected for precession and nutation of the      Earth's axis of rotation at time <code>t</code>. |


* * *

<a name="Astronomy.HorizontalCoordinates"></a>

### Astronomy.HorizontalCoordinates
Holds azimuth (compass direction) and altitude (angle above/below the horizon)
of a celestial object as seen by an observer at a particular location on the Earth's surface.
Also holds right ascension and declination of the same object.
All of these coordinates are optionally adjusted for atmospheric refraction;
therefore the right ascension and declination values may not exactly match
those found inside a corresponding [EquatorialCoordinates](#Astronomy.EquatorialCoordinates) object.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| azimuth | <code>number</code> | A horizontal compass direction angle in degrees measured starting at north       and increasing positively toward the east.      The value is in the range [0, 360).      North = 0, east = 90, south = 180, west = 270. |
| altitude | <code>number</code> | A vertical angle in degrees above (positive) or below (negative) the horizon.      The value is in the range [-90, +90].      The altitude angle is optionally adjusted upward due to atmospheric refraction. |
| ra | <code>number</code> | The right ascension of the celestial body in sidereal hours.      The value is in the reange [0, 24).      If <code>altitude</code> was adjusted for atmospheric reaction, <code>ra</code>      is likewise adjusted. |
| dec | <code>number</code> | The declination of of the celestial body in degrees.      The value in the range [-90, +90].      If <code>altitude</code> was adjusted for atmospheric reaction, <code>dec</code>      is likewise adjusted. |


* * *

<a name="Astronomy.EclipticCoordinates"></a>

### Astronomy.EclipticCoordinates
Holds ecliptic coordinates of a celestial body.
The origin and date of the coordinate system may vary depending on the caller's usage.
In general, ecliptic coordinates are measured with respect to the mean plane of the Earth's 
orbit around the Sun.
Includes Cartesian coordinates <code>(ex, ey, ez)</code> measured in 
<a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU)
and spherical coordinates <code>(elon, elat)</code> measured in degrees.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| ex | <code>number</code> | The Cartesian x-coordinate of the body in astronomical units (AU).      The x-axis is within the ecliptic plane and is oriented in the direction of the       <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>. |
| ey | <code>number</code> | The Cartesian y-coordinate of the body in astronomical units (AU).      The y-axis is within the ecliptic plane and is oriented 90 degrees       counterclockwise from the equinox, as seen from above the Sun's north pole. |
| ez | <code>number</code> | The Cartesian z-coordinate of the body in astronomical units (AU).      The z-axis is oriented perpendicular to the ecliptic plane,      along the direction of the Sun's north pole. |
| elat | <code>number</code> | The ecliptic latitude of the body in degrees.      This is the angle north or south of the ecliptic plane.      The value is in the range [-90, +90].      Positive values are north and negative values are south. |
| elon | <code>number</code> | The ecliptic longitude of the body in degrees.      This is the angle measured counterclockwise around the ecliptic plane,      as seen from above the Sun's north pole.      This is the same direction that the Earth orbits around the Sun.      The angle is measured starting at 0 from the equinox and increases      up to 360 degrees. |


* * *

<a name="Astronomy.Observer"></a>

### Astronomy.Observer
Represents the geographic location of an observer on the surface of the Earth.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| latitude_degrees | <code>number</code> | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude_degrees | <code>number</code> | The observer's geographic longitude in degrees east of the prime meridian       passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height_in_meters | <code>number</code> | The observer's elevation above mean sea level, expressed in meters. |


* * *

<a name="Astronomy.IlluminationInfo"></a>

### Astronomy.IlluminationInfo
Contains information about the apparent brightness and sunlit phase of a celestial object.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | The date and time pertaining to the other calculated values in this object. |
| mag | <code>number</code> | The <a href="https://en.wikipedia.org/wiki/Apparent_magnitude">apparent visual magnitude</a> of the celestial body. |
| phase_angle | <code>number</code> | The angle in degrees as seen from the center of the celestial body between the Sun and the Earth.      The value is always in the range 0 to 180.      The phase angle provides a measure of what fraction of the body's face appears       illuminated by the Sun as seen from the Earth.      When the observed body is the Sun, the <code>phase</code> property is set to 0,      although this has no physical meaning because the Sun emits, rather than reflects, light.      When the phase is near 0 degrees, the body appears "full".      When it is 90 degrees, the body appears "half full".       And when it is 180 degrees, the body appears "new" and is very difficult to see      because it is both dim and lost in the Sun's glare as seen from the Earth. |
| phase_fraction | <code>number</code> | The fraction of the body's face that is illuminated by the Sun, as seen from the Earth.      Calculated from <code>phase_angle</code> for convenience.      This value ranges from 0 to 1. |
| helio_dist | <code>number</code> | The distance between the center of the Sun and the center of the body in       <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU). |
| geo_dist | <code>number</code> | The distance between the center of the Earth and the center of the body in AU. |
| gc | [<code>Vector</code>](#Astronomy.Vector) | Geocentric coordinates: the 3D vector from the center of the Earth to the center of the body.      The components are in expressed in AU and are oriented with respect to the J2000 equatorial plane. |
| hc | [<code>Vector</code>](#Astronomy.Vector) | Heliocentric coordinates: The 3D vector from the center of the Sun to the center of the body.      Like <code>gc</code>, <code>hc</code> is expressed in AU and oriented with respect      to the J2000 equatorial plane. |
| ring_tilt | <code>number</code> \| <code>null</code> | For Saturn, this is the angular tilt of the planet's rings in degrees away      from the line of sight from the Earth. When the value is near 0, the rings      appear edge-on from the Earth and are therefore difficult to see.      When <code>ring_tilt</code> approaches its maximum value (about 27 degrees),      the rings appear widest and brightest from the Earth.      Unlike the <a href="https://ssd.jpl.nasa.gov/horizons.cgi">JPL Horizons</a> online tool,       this library includes the effect of the ring tilt angle in the calculated value       for Saturn's visual magnitude.      For all bodies other than Saturn, the value of <code>ring_tilt</code> is <code>null</code>. |


* * *

<a name="Astronomy.MoonQuarter"></a>

### Astronomy.MoonQuarter
Represents a quarter lunar phase, along with when it occurs.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| quarter | <code>number</code> | An integer as follows:      0 = new moon,       1 = first quarter,      2 = full moon,      3 = third quarter. |
| time | [<code>Time</code>](#Astronomy.Time) | The date and time of the quarter lunar phase. |


* * *

<a name="Astronomy.HourAngleEvent"></a>

### Astronomy.HourAngleEvent
Returns information about an occurrence of a celestial body
reaching a given hour angle as seen by an observer at a given
location on the surface of the Earth.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | The date and time of the celestial body reaching the hour angle. |
| pos | [<code>Vector</code>](#Astronomy.Vector) | Geocentric Cartesian coordinates for the body in the J2000 equatorial system      at the time indicated by the <code>time</code> property. |
| sky | [<code>SkyCoordinates</code>](#Astronomy.SkyCoordinates) | Topocentric equatorial coordinates for the body      at the time indicated by the <code>time</code> property. |
| hor | [<code>HorizontalCoordinates</code>](#Astronomy.HorizontalCoordinates) | Topocentric horizontal coordinates for the body      at the time indicated by the <code>time</code> property. |
| iter | <code>number</code> | The positive integer number of iterations required by      <code>SearchHourAngle</code> to converge on the hour angle      solution. |


* * *

<a name="Astronomy.SeasonInfo"></a>

### Astronomy.SeasonInfo
Represents the dates and times of the two solstices
and the two equinoxes in a given calendar year.
These four events define the changing of the seasons on the Earth.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  

* * *

<a name="new_Astronomy.SeasonInfo_new"></a>

#### new SeasonInfo(mar_equinox, jun_solstice, sep_equinox, dec_solstice)

| Param | Type | Description |
| --- | --- | --- |
| mar_equinox | [<code>Time</code>](#Astronomy.Time) | The date and time of the March equinox in the given calendar year.      This is the moment in March that the plane of the Earth's equator passes      through the center of the Sun; thus the Sun's declination      changes from a negative number to a positive number.      The March equinox defines      the beginning of spring in the northern hemisphere and      the beginning of autumn in the southern hemisphere. |
| jun_solstice | [<code>Time</code>](#Astronomy.Time) | The date and time of the June solstice in the given calendar year.      This is the moment in June that the Sun reaches its most positive      declination value.      At this moment the Earth's north pole is most tilted most toward the Sun.      The June solstice defines      the beginning of summer in the northern hemisphere and      the beginning of winter in the southern hemisphere. |
| sep_equinox | [<code>Time</code>](#Astronomy.Time) | The date and time of the September equinox in the given calendar year.      This is the moment in September that the plane of the Earth's equator passes      through the center of the Sun; thus the Sun's declination      changes from a positive number to a negative number.      The September equinox defines      the beginning of autumn in the northern hemisphere and      the beginning of spring in the southern hemisphere. |
| dec_solstice | [<code>Time</code>](#Astronomy.Time) | The date and time of the December solstice in the given calendar year.      This is the moment in December that the Sun reaches its most negative      declination value.      At this moment the Earth's south pole is tilted most toward the Sun.      The December solstice defines      the beginning of winter in the northern hemisphere and      the beginning of summer in the southern hemisphere. |


* * *

<a name="Astronomy.ElongationEvent"></a>

### Astronomy.ElongationEvent
Represents the visibility of a planet or the Moon relative to the Sun.
Includes angular separation from the Sun and whether visibility is
best in the morning or the evening.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**See**: [Elongation](#Astronomy.Elongation)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | When the event occurs. |
| visibility | <code>string</code> | Either <code>"morning"</code> or <code>"evening"</code>,      indicating when the body is most easily seen. |
| elongation | <code>number</code> | The angle in degrees, as seen from the center of the Earth,      of the apparent separation between the body and the Sun.      This angle is measured in 3D space and is not projected onto the ecliptic plane. |
| relative_longitude | <code>number</code> | The angle in degrees, as seen from the Sun, between the      observed body and the Earth. This value is always between      0 and 180. More precisely, <code>relative_longitude</code> is the absolute      value of the difference between the heliocentric ecliptic longitudes of      the centers of the observed body and the Earth. |


* * *

<a name="Astronomy.Apsis"></a>

### Astronomy.Apsis
Represents a closest or farthest point in a body's orbit around its primary.
For a planet orbiting the Sun, this is a perihelion or aphelion, respectively.
For the Moon orbiting the Earth, this is a perigee or apogee, respectively.

**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**See**

- [SearchLunarApsis](#Astronomy.SearchLunarApsis)
- [NextLunarApsis](#Astronomy.NextLunarApsis)

**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | The date and time of the apsis. |
| apsisType | <code>number</code> | For a closest approach (perigee or perihelion), <code>apsisType</code> is 0.      For a farthest distance event (apogee or aphelion), <code>apsisType</code> is 1. |
| dist_au | <code>number</code> | The distance between the centers of the two bodies in astronomical units (AU). |
| dist_km | <code>number</code> | The distance between the centers of the two bodies in kilometers. |


* * *

<a name="Astronomy.Bodies"></a>

### Astronomy.Bodies : <code>Array.&lt;string&gt;</code>
An array of strings, each a name of a supported astronomical body.
     Not all bodies are valid for all functions, but any string not in this
     list is not supported at all.

**Kind**: static constant of [<code>Astronomy</code>](#Astronomy)  

* * *

<a name="Astronomy.GetPerformanceMetrics"></a>

### Astronomy.GetPerformanceMetrics() ⇒ [<code>PerformanceInfo</code>](#Astronomy.PerformanceInfo)
Takes a snapshot of the current state of the performance metrics.
The metrics inside the returned object will not change and can be retained by calling code
to be compared with later snapshots.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

* * *

<a name="Astronomy.ResetPerformanceMetrics"></a>

### Astronomy.ResetPerformanceMetrics()
Resets the internal performance metrics back to their initial states.
You can call this before starting a new series of performance tests.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

* * *

<a name="Astronomy.MakeTime"></a>

### Astronomy.MakeTime(date) ⇒ [<code>Time</code>](#Astronomy.Time)
Given a Date object or a number days since noon (12:00) on January 1, 2000 (UTC),
this function creates an [Time](#Astronomy.Time) object.
Given an [Time](#Astronomy.Time) object, returns the same object unmodified.
Use of this function is not required for any of the other exposed functions in this library,
because they all guarantee converting date/time parameters to Astronomy.Time
as needed. However, it may be convenient for callers who need to understand
the difference between UTC and TT (Terrestrial Time). In some use cases,
converting once to Astronomy.Time format and passing the result into multiple
function calls may be more efficient than passing in native JavaScript Date objects.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | A Date object, a number of UTC days since the J2000 epoch (noon on January 1, 2000),      or an Astronomy.Time object. See remarks above. |


* * *

<a name="Astronomy.Horizon"></a>

### Astronomy.Horizon(date, location, ra, dec, refraction) ⇒ [<code>HorizontalCoordinates</code>](#Astronomy.HorizontalCoordinates)
Given a date and time, a geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial object,
returns horizontal coordinates (azimuth and altitude angles) for that object
as seen by that observer. Allows optional correction for atmospheric refraction.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to find horizontal coordinates. |
| location | [<code>Observer</code>](#Astronomy.Observer) | The location of the observer for which to find horizontal coordinates. |
| ra | <code>number</code> | Right ascension in sidereal hours of the celestial object,       referred to the mean equinox of date for the J2000 epoch. |
| dec | <code>number</code> | Declination in degrees of the celestial object,       referred to the mean equator of date for the J2000 epoch.      Positive values are north of the celestial equator and negative values are south. |
| refraction | <code>string</code> | If omitted or has a false-like value (false, null, undefined, etc.)      the calculations are performed without any correction for atmospheric      refraction. If the value is the string <code>"normal"</code>,      uses the recommended refraction correction based on Meeus "Astronomical Algorithms"      with a linear taper more than 1 degree below the horizon. The linear      taper causes the refraction to linearly approach 0 as the altitude of the      body approaches the nadir (-90 degrees).      If the value is the string <code>"jplhor"</code>, uses a JPL Horizons      compatible formula. This is the same algorithm as <code>"normal"</code>,       only without linear tapering; this can result in physically impossible      altitudes of less than -90 degrees, which may cause problems for some applications.      (The <code>"jplhor"</code> option was created for unit testing against data      generated by JPL Horizons, and is otherwise not recommended for use.) |


* * *

<a name="Astronomy.MakeObserver"></a>

### Astronomy.MakeObserver(latitude_degrees, longitude_degrees, height_in_meters)
Creates an [Observer](#Astronomy.Observer) object that represents a location
on the surface of the Earth from which observations are made.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| latitude_degrees | <code>number</code> | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude_degrees | <code>number</code> | The observer's geographic longitude in degrees east of the prime meridian       passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height_in_meters | <code>number</code> | The observer's elevation above mean sea level, expressed in meters.      If omitted, the elevation is assumed to be 0 meters. |


* * *

<a name="Astronomy.SunPosition"></a>

### Astronomy.SunPosition(date) ⇒ [<code>EclipticCoordinates</code>](#Astronomy.EclipticCoordinates)
Returns apparent geocentric true ecliptic coordinates of date for the Sun.
<i>Geocentric</i> means coordinates as the Sun would appear to a hypothetical observer
at the center of the Earth.
<i>Ecliptic coordinates of date</i> are measured along the plane of the Earth's mean
orbit around the Sun, using the 
<a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>
of the Earth as adjusted for precession and nutation of the Earth's
axis of rotation on the given date.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time at which to calculate the Sun's apparent location as seen from      the center of the Earth. |


* * *

<a name="Astronomy.SkyPos"></a>

### Astronomy.SkyPos(gc_vector, observer) ⇒ [<code>SkyCoordinates</code>](#Astronomy.SkyCoordinates)
Returns topocentric equatorial coordinates (right ascension and declination)
simultaneously in two different systems: J2000 and true-equator-of-date.
<i>Topocentric</i> refers to a position as seen by an observer on the surface of the Earth.
This function corrects for
<a href="https://en.wikipedia.org/wiki/Parallax">parallax</a> 
of the object between a geocentric observer and a topocentric observer.
This is most significant for the Moon, because it is so close to the Earth.
However, it can have a small effect on the apparent positions of other bodies.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: [<code>SkyCoordinates</code>](#Astronomy.SkyCoordinates) - The topocentric coordinates of the body as adjusted for the given observer.  

| Param | Type | Description |
| --- | --- | --- |
| gc_vector | [<code>Vector</code>](#Astronomy.Vector) | A geocentric vector in the J2000 equatorial system.      <i>Geocentric</i> refers to a position seen by a hypothetical observer at the center of the Earth. |
| observer | [<code>Observer</code>](#Astronomy.Observer) | The location on the Earth of the observer. |


* * *

<a name="Astronomy.Ecliptic"></a>

### Astronomy.Ecliptic(gx, gy, gz) ⇒ [<code>EclipticCoordinates</code>](#Astronomy.EclipticCoordinates)
Given J2000 equatorial Cartesian coordinates, 
returns J2000 ecliptic latitude, longitude, and cartesian coordinates.
You can call [GeoVector](#Astronomy.GeoVector) and use its (x, y, z) return values
to pass into this function.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| gx | <code>number</code> | The x-coordinate of a 3D vector in the J2000 equatorial coordinate system. |
| gy | <code>number</code> | The y-coordinate of a 3D vector in the J2000 equatorial coordinate system. |
| gz | <code>number</code> | The z-coordinate of a 3D vector in the J2000 equatorial coordinate system. |


* * *

<a name="Astronomy.GeoMoon"></a>

### Astronomy.GeoMoon(date) ⇒ [<code>Vector</code>](#Astronomy.Vector)
Calculates the geocentric Cartesian coordinates for the Moon in the J2000 equatorial system.
Based on the Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
which in turn derives from E. W. Brown's lunar theories.
Adapted from Turbo Pascal code from the book 
<a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a> 
by Montenbruck and Pfleger.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to calculate the Moon's geocentric position. |


* * *

<a name="Astronomy.HelioVector"></a>

### Astronomy.HelioVector(body, date) ⇒ [<code>Vector</code>](#Astronomy.Vector)
Calculates heliocentric (i.e., with respect to the center of the Sun)
Cartesian coordinates in the J2000 equatorial system of a celestial
body at a specified time.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | One of the strings       <code>"Sun"</code>, <code>"Moon"</code>, <code>"Mercury"</code>, <code>"Venus"</code>,       <code>"Earth"</code>, <code>"Mars"</code>, <code>"Jupiter"</code>, <code>"Saturn"</code>,       <code>"Uranus"</code>, <code>"Neptune"</code>, or <code>"Pluto"</code>. |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which the body's position is to be calculated. |


* * *

<a name="Astronomy.GeoVector"></a>

### Astronomy.GeoVector(body, date) ⇒ [<code>Vector</code>](#Astronomy.Vector)
Calculates geocentric (i.e., with respect to the center of the Earth)
Cartesian coordinates in the J2000 equatorial system of a celestial
body at a specified time.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | One of the strings       <code>"Sun"</code>, <code>"Moon"</code>, <code>"Mercury"</code>, <code>"Venus"</code>,       <code>"Earth"</code>, <code>"Mars"</code>, <code>"Jupiter"</code>, <code>"Saturn"</code>,       <code>"Uranus"</code>, <code>"Neptune"</code>, or <code>"Pluto"</code>. |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which the body's position is to be calculated. |


* * *

<a name="Astronomy.Search"></a>

### Astronomy.Search(func, t1, t2, options) ⇒ <code>null</code> \| [<code>Time</code>](#Astronomy.Time)
Search for next time <i>t</i> (such that <i>t</i> is between <code>t1</code> and <code>t2</code>)
that <code>func(t)</code> crosses from a negative value to a non-negative value.
The given function must have "smooth" behavior over the entire inclusive range [<code>t1</code>, <code>t2</code>],
meaning that it behaves like a continuous differentiable function.
It is not required that <code>t1</code> &lt; <code>t2</code>; <code>t1</code> &gt; <code>t2</code> 
allows searching backward in time.
Note: <code>t1</code> and <code>t2</code> must be chosen such that there is no possibility
of more than one zero-crossing (ascending or descending), or it is possible
that the "wrong" event will be found (i.e. not the first event after t1)
or even that the function will return null, indicating that no event was found.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: <code>null</code> \| [<code>Time</code>](#Astronomy.Time) - If the search is successful, returns the date and time of the solution.
     If the search fails, returns null.  

| Param | Type | Description |
| --- | --- | --- |
| func | [<code>ContinuousFunction</code>](#Astronomy.ContinuousFunction) | The function to find an ascending zero crossing for.      The function must accept a single parameter of type [Time](#Astronomy.Time)      and return a numeric value. |
| t1 | [<code>Time</code>](#Astronomy.Time) | The lower time bound of a search window. |
| t2 | [<code>Time</code>](#Astronomy.Time) | The upper time bound of a search window. |
| options | <code>null</code> \| [<code>SearchOptions</code>](#Astronomy.SearchOptions) | Options that can tune the behavior of the search.      Most callers can omit this argument or pass in <code>null</code>. |


* * *

<a name="Astronomy.SearchSunLongitude"></a>

### Astronomy.SearchSunLongitude(targetLon, dateStart, limitDays) ⇒ [<code>Time</code>](#Astronomy.Time) \| <code>null</code>
Searches for the moment in time when the center of the Sun reaches a given apparent
ecliptic longitude, as seen from the center of the Earth, within a given range of dates.
This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call [Seasons](#Astronomy.Seasons) 
to calculate equinoxes and solstices for a given calendar year.
<code>SearchSunLongitude</code> is more general in that it allows searching for arbitrary longitude values.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: [<code>Time</code>](#Astronomy.Time) \| <code>null</code> - The date and time when the Sun reaches the apparent ecliptic longitude <code>targetLon</code>
     within the range of times specified by <code>dateStart</code> and <code>limitDays</code>.
     If the Sun does not reach the target longitude within the specified time range, or the
     time range is excessively wide, the return value is <code>null</code>.
     To avoid a <code>null</code> return value, the caller must pick a time window around
     the event that is within a few days but not so small that the event might fall outside the window.  

| Param | Type | Description |
| --- | --- | --- |
| targetLon | <code>number</code> | The desired ecliptic longitude of date in degrees.      This may be any value in the range [0, 360), although certain      values have conventional meanings:      When <code>targetLon</code> is 0, finds the March equinox,      which is the moment spring begins in the northern hemisphere      and the beginning of autumn in the southern hemisphere.      When <code>targetLon</code> is 180, finds the September equinox,      which is the moment autumn begins in the northern hemisphere and      spring begins in the southern hemisphere.      When <code>targetLon</code> is 90, finds the northern solstice, which is the      moment summer begins in the northern hemisphere and winter      begins in the southern hemisphere.      When <code>targetLon</code> is 270, finds the southern solstice, which is the      moment winter begins in the northern hemisphere and summer      begins in the southern hemisphere. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | A date and time known to be earlier than the desired longitude event. |
| limitDays | <code>number</code> | A floating point number of days, which when added to <code>dateStart</code>,      yields a date and time known to be after the desired longitude event. |


* * *

<a name="Astronomy.LongitudeFromSun"></a>

### Astronomy.LongitudeFromSun(body, date) ⇒ <code>number</code>
Calculates the ecliptic longitude difference 
between the given body and the Sun as seen from 
the Earth at a given moment in time.
The returned value ranges [0, 360) degrees.
By definition, the Earth and the Sun are both in the plane of the ecliptic.
Ignores the height of the <code>body</code> above or below the ecliptic plane;
the resulting angle is measured around the ecliptic plane for the "shadow"
of the body onto that plane.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: <code>number</code> - An angle in degrees in the range [0, 360).
     Values less than 180 indicate that the body is to the east
     of the Sun as seen from the Earth; that is, the body sets after
     the Sun does and is visible in the evening sky.
     Values greater than 180 indicate that the body is to the west of
     the Sun and is visible in the morning sky.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a supported celestial body other than the Earth. |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The time at which the relative longitude is to be found. |


* * *

<a name="Astronomy.AngleFromSun"></a>

### Astronomy.AngleFromSun(body, date) ⇒ <code>number</code>
Returns the full angle seen from
the Earth, between the given body and the Sun.
Unlike [LongitudeFromSun](#Astronomy.LongitudeFromSun), this function does not
project the body's "shadow" onto the ecliptic; 
the angle is measured in 3D space around the plane that 
contains the centers of the Earth, the Sun, and <code>body</code>.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: <code>number</code> - An angle in degrees in the range [0, 180].  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a supported celestial body other than the Earth. |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The time at which the angle from the Sun is to be found. |


* * *

<a name="Astronomy.EclipticLongitude"></a>

### Astronomy.EclipticLongitude(body, date) ⇒ <code>number</code>
Calculates heliocentric ecliptic longitude based on the J2000 equinox.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: <code>number</code> - The ecliptic longitude angle of the body in degrees measured counterclockwise around the mean
     plane of the Earth's orbit, as seen from above the Sun's north pole.
     Ecliptic longitude starts at 0 at the J2000
     <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a> and
     increases in the same direction the Earth orbits the Sun.
     The returned value is always in the range [0, 360).  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a celestial body other than the Sun. |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to calculate the ecliptic longitude. |


* * *

<a name="Astronomy.Illumination"></a>

### Astronomy.Illumination(body, date) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
Calculates the phase angle, visual maginitude, 
and other values relating to the body's illumination
at the given date and time, as seen from the Earth.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the celestial body being observed.      Not allowed to be <code>"Earth"</code>. |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to calculate the illumination data for the given body. |


* * *

<a name="Astronomy.SearchRelativeLongitude"></a>

### Astronomy.SearchRelativeLongitude(body, targetRelLon, startDate) ⇒ [<code>Time</code>](#Astronomy.Time)
Searches for the date and time the relative ecliptic longitudes of
the specified body and the Earth, as seen from the Sun, reach a certain
difference. This function is useful for finding conjunctions and oppositions
of the planets. For the opposition of a superior planet (Mars, Jupiter, ..., Pluto),
or the inferior conjunction of an inferior planet (Mercury, Venus),
call with <code>targetRelLon</code> = 0. The 0 value indicates that both
planets are on the same ecliptic longitude line, ignoring the other planet's
distance above or below the plane of the Earth's orbit.
For superior conjunctions, call with <code>targetRelLon</code> = 180.
This means the Earth and the other planet are on opposite sides of the Sun.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: [<code>Time</code>](#Astronomy.Time) - The time when the Earth and the body next reach the specified relative longitudes.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a planet other than the Earth. |
| targetRelLon | <code>number</code> | The desired angular difference in degrees between the ecliptic longitudes      of <code>body</code> and the Earth. Must be in the range (-180, +180]. |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which to find the next occurrence of the      body and the Earth reaching the desired relative longitude. |


* * *

<a name="Astronomy.MoonPhase"></a>

### Astronomy.MoonPhase(date) ⇒ <code>number</code>
Determines the moon's phase expressed as an ecliptic longitude.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: <code>number</code> - A value in the range [0, 360) indicating the difference
     in ecliptic longitude between the center of the Sun and the
     center of the Moon, as seen from the center of the Earth.
     Certain longitude values have conventional meanings:

* 0 = new moon
* 90 = first quarter
* 180 = full moon
* 270 = third quarter  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to calculate the moon's phase. |


* * *

<a name="Astronomy.SearchMoonPhase"></a>

### Astronomy.SearchMoonPhase(targetLon, dateStart, limitDays) ⇒ [<code>Time</code>](#Astronomy.Time) \| <code>null</code>
Searches for the date and time that the Moon reaches a specified phase.
Lunar phases are defined in terms of geocentric ecliptic longitudes
with respect to the Sun.  When the Moon and the Sun have the same ecliptic
longitude, that is defined as a new moon. When the two ecliptic longitudes
are 180 degrees apart, that is defined as a full moon.
To enumerate quarter lunar phases, it is simpler to call
[SearchMoonQuarter](#Astronomy.SearchMoonQuarter) once, followed by repeatedly calling
[NextMoonQuarter](#Astronomy.NextMoonQuarter). <code>SearchMoonPhase</code> is only
necessary for finding other lunar phases than the usual quarter phases.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: [<code>Time</code>](#Astronomy.Time) \| <code>null</code> - If the specified lunar phase occurs after <code>dateStart</code>
     and before <code>limitDays</code> days after <code>dateStart</code>,
     this function returns the date and time of the first such occurrence.
     Otherwise, it returns <code>null</code>.  

| Param | Type | Description |
| --- | --- | --- |
| targetLon | <code>number</code> | The difference in geocentric ecliptic longitude between the Sun and Moon      that specifies the lunar phase being sought. This can be any value      in the range [0, 360). Here are some helpful examples:      0 = new moon,      90 = first quarter,      180 = full moon,       270 = third quarter. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The beginning of the window of time in which to search. |
| limitDays | <code>number</code> | The floating point number of days after <code>dateStart</code>      that limits the window of time in which to search. |


* * *

<a name="Astronomy.SearchMoonQuarter"></a>

### Astronomy.SearchMoonQuarter(dateStart) ⇒ [<code>MoonQuarter</code>](#Astronomy.MoonQuarter)
Finds the first quarter lunar phase after the specified date and time.
The quarter lunar phases are: new moon, first quarter, full moon, and third quarter.
To enumerate quarter lunar phases, call <code>SearchMoonQuarter</code> once,
then pass its return value to [NextMoonQuarter](#Astronomy.NextMoonQuarter) to find the next
<code>MoonQuarter</code>. Keep calling <code>NextMoonQuarter</code> in a loop,
passing the previous return value as the argument to the next call.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which to find the first quarter lunar phase. |


* * *

<a name="Astronomy.NextMoonQuarter"></a>

### Astronomy.NextMoonQuarter(mq)
Given a [MoonQuarter](#Astronomy.MoonQuarter) object, finds the next consecutive
quarter lunar phase. See remarks in [SearchMoonQuarter](#Astronomy.SearchMoonQuarter)
for explanation of usage.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| mq | [<code>MoonQuarter</code>](#Astronomy.MoonQuarter) | The return value of a prior call to [MoonQuarter](#Astronomy.MoonQuarter) or <code>NextMoonQuarter</code>. |


* * *

<a name="Astronomy.SearchRiseSet"></a>

### Astronomy.SearchRiseSet(body, observer, direction, dateStart, limitDays) ⇒ [<code>Time</code>](#Astronomy.Time) \| <code>null</code>
Finds a rise or set time for the given body as 
seen by an observer at the specified location on the Earth.
Rise time is defined as the moment when the top of the body
is observed to first appear above the horizon in the east.
Set time is defined as the moment the top of the body
is observed to sink below the horizon in the west.
The times are adjusted for typical atmospheric refraction conditions.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: [<code>Time</code>](#Astronomy.Time) \| <code>null</code> - The date and time of the rise or set event, or null if no such event
     occurs within the specified time window.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the body to find the rise or set time for. |
| observer | [<code>Observer</code>](#Astronomy.Observer) | Specifies the geographic coordinates and elevation above sea level of the observer.      Call [MakeObserver](#Astronomy.MakeObserver) to create an observer object. |
| direction | <code>number</code> | Either +1 to find rise time or -1 to find set time.      Any other value will cause an exception to be thrown. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which the specified rise or set time is to be found. |
| limitDays | <code>number</code> | The fractional number of days after <code>dateStart</code> that limits      when the rise or set time is to be found. |


* * *

<a name="Astronomy.SearchHourAngle"></a>

### Astronomy.SearchHourAngle(body, observer, hourAngle, dateStart) ⇒ [<code>HourAngleEvent</code>](#Astronomy.HourAngleEvent)
Finds the next time the given body is seen to reach the specified 
<a href="https://en.wikipedia.org/wiki/Hour_angle">hour angle</a>
by the given observer.
Providing <code>hourAngle</code> = 0 finds the next maximum altitude event (culmination).
Providing <code>hourAngle</code> = 12 finds the next minimum altitude event.
Note that, especially close to the Earth's poles, a body as seen on a given day
may always be above the horizon or always below the horizon, so the caller cannot
assume that a culminating object is visible nor that an object is below the horizon
at its minimum altitude.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a celestial body other than the Earth. |
| observer | [<code>Observer</code>](#Astronomy.Observer) | Specifies the geographic coordinates and elevation above sea level of the observer.      Call [MakeObserver](#Astronomy.MakeObserver) to create an observer object. |
| hourAngle | <code>number</code> | The hour angle expressed in       <a href="https://en.wikipedia.org/wiki/Sidereal_time">sidereal</a>       hours for which the caller seeks to find the body attain.       The value must be in the range [0, 24).      The hour angle represents the number of sidereal hours that have       elapsed since the most recent time the body crossed the observer's local      <a href="https://en.wikipedia.org/wiki/Meridian_(astronomy)">meridian</a>.      This specifying <code>hourAngle</code> = 0 finds the moment in time      the body reaches the highest angular altitude in a given sidereal day. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which the desired hour angle crossing event      is to be found. |


* * *

<a name="Astronomy.Seasons"></a>

### Astronomy.Seasons(year) ⇒ [<code>SeasonInfo</code>](#Astronomy.SeasonInfo)
Find the equinoxes and solstices for a given calendar year.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| year | <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The integer value or <code>Time</code> object that specifies      the UTC calendar year for which to find equinoxes and solstices. |


* * *

<a name="Astronomy.Elongation"></a>

### Astronomy.Elongation(body) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
Calculates the absolute value of the angle between the centers 
of the given body and the Sun as seen from the center of the Earth at the given date.
The angle is measured along the plane of the Earth's orbit (i.e. the ecliptic) 
and ranges [0, 180] degrees. This function is helpful for determining how easy 
it is to view Mercury or Venus away from the Sun's glare on a given date.
The function also determines whether the object is visible in the morning or evening; 
this is more important the smaller the elongation is.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the observed body. Not allowed to be <code>"Earth"</code>. |


* * *

<a name="Astronomy.SearchMaxElongation"></a>

### Astronomy.SearchMaxElongation(body, startDate) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
Searches for the next maximum elongation event for Mercury or Venus 
that occurs after the given start date. Calling with other values 
of <code>body</code> will result in an exception.
Maximum elongation occurs when the body has the greatest
angular separation from the Sun, as seen from the Earth.
Returns an <code>ElongationEvent</code> object containing the date and time of the next
maximum elongation, the elongation in degrees, and whether
the body is visible in the morning or evening.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | Either <code>"Mercury"</code> or <code>"Venus"</code>. |
| startDate | <code>Date</code> | The date and time after which to search for the next maximum elongation event. |


* * *

<a name="Astronomy.SearchPeakMagnitude"></a>

### Astronomy.SearchPeakMagnitude(body, startDate) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
Searches for the date and time Venus will next appear brightest as seen from the Earth.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | Currently only <code>"Venus"</code> is supported.      Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from Earth,      so peak magnitude events have little practical value for that planet.      The Moon reaches peak magnitude very close to full moon, which can be found using       [SearchMoonQuarter](#Astronomy.SearchMoonQuarter) or [SearchMoonPhase](#Astronomy.SearchMoonPhase).      The other planets reach peak magnitude very close to opposition,       which can be found using [SearchRelativeLongitude](#Astronomy.SearchRelativeLongitude). |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which to find the next peak magnitude event. |


* * *

<a name="Astronomy.SearchLunarApsis"></a>

### Astronomy.SearchLunarApsis(startDate) ⇒ [<code>Apsis</code>](#Astronomy.Apsis)
Finds the next perigee (closest approach) or apogee (farthest remove) of the Moon
that occurs after the specified date and time.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which to find the next perigee or apogee. |


* * *

<a name="Astronomy.NextLunarApsis"></a>

### Astronomy.NextLunarApsis(apsis) ⇒ [<code>Apsis</code>](#Astronomy.Apsis)
Given a lunar apsis returned by an initial call to [SearchLunarApsis](SearchLunarApsis), 
or a previous call to <code>NextLunarApsis</code>, finds the next lunar apsis.
If the given apsis is a perigee, this function finds the next apogee, and vice versa.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: [<code>Apsis</code>](#Astronomy.Apsis) - The successor apogee for the given perigee, or the successor perigee for the given apogee.  

| Param | Type | Description |
| --- | --- | --- |
| apsis | [<code>Apsis</code>](#Astronomy.Apsis) | A lunar perigee or apogee event. |


* * *

<a name="Astronomy.ContinuousFunction"></a>

### Astronomy.ContinuousFunction ⇒ <code>number</code>
A continuous function of time used in a call to the <code>Search</code> function.

**Kind**: static typedef of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| t | [<code>Time</code>](#Astronomy.Time) | The time at which to evaluate the function. |


* * *

<a name="Astronomy.SearchOptions"></a>

### Astronomy.SearchOptions : <code>Object</code>
Options for the [Search](#Astronomy.Search) function.

**Kind**: static typedef of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| dt_tolerance_seconds | <code>number</code> \| <code>null</code> | The number of seconds for a time window smaller than which the search      is considered successful.  Using too large a tolerance can result in      an inaccurate time estimate.  Using too small a tolerance can cause      excessive computation, or can even cause the search to fail because of      limited floating-point resolution.  Defaults to 1 second. |
| init_f1 | <code>number</code> \| <code>null</code> | As an optimization, if the caller of [Search](#Astronomy.Search)       has already calculated the value of the function being searched (the parameter <code>func</code>)       at the time coordinate <code>t1</code>, it can pass in that value as <code>init_f1</code>.      For very expensive calculations, this can measurably improve performance. |
| init_f2 | <code>number</code> \| <code>null</code> | The same as <code>init_f1</code>, except this is the optional initial value of <code>func(t2)</code>      instead of <code>func(t1)</code>. |


* * *

