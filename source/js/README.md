# API Reference
<a name="Astronomy"></a>

## Astronomy : <code>object</code>
**Kind**: global namespace  

* [Astronomy](#Astronomy) : <code>object</code>
    * [.PerformanceInfo](#Astronomy.PerformanceInfo)
        * [new PerformanceInfo()](#new_Astronomy.PerformanceInfo_new)
        * [.Clone()](#Astronomy.PerformanceInfo+Clone) ⇒ [<code>PerformanceInfo</code>](#Astronomy.PerformanceInfo)
    * [.Time](#Astronomy.Time)
        * [new Time(date)](#new_Astronomy.Time_new)
        * [.toString()](#Astronomy.Time+toString) ⇒ <code>string</code>
        * [.AddDays(days)](#Astronomy.Time+AddDays) ⇒ [<code>Time</code>](#Astronomy.Time)
    * [.Vector](#Astronomy.Vector)
        * [new Vector()](#new_Astronomy.Vector_new)
    * [.EquatorialCoordinates](#Astronomy.EquatorialCoordinates)
        * [new EquatorialCoordinates()](#new_Astronomy.EquatorialCoordinates_new)
    * [.HorizontalCoordinates](#Astronomy.HorizontalCoordinates)
        * [new HorizontalCoordinates()](#new_Astronomy.HorizontalCoordinates_new)
    * [.EclipticCoordinates](#Astronomy.EclipticCoordinates)
        * [new EclipticCoordinates()](#new_Astronomy.EclipticCoordinates_new)
    * [.Observer](#Astronomy.Observer)
        * [new Observer()](#new_Astronomy.Observer_new)
    * [.IlluminationInfo](#Astronomy.IlluminationInfo)
        * [new IlluminationInfo()](#new_Astronomy.IlluminationInfo_new)
    * [.ElongationEvent](#Astronomy.ElongationEvent)
        * [new ElongationEvent()](#new_Astronomy.ElongationEvent_new)
    * [.Bodies](#Astronomy.Bodies) : <code>Array.&lt;string&gt;</code>
    * [.GetPerformanceMetrics()](#Astronomy.GetPerformanceMetrics) ⇒ [<code>PerformanceInfo</code>](#Astronomy.PerformanceInfo)
    * [.ResetPerformanceMetrics()](#Astronomy.ResetPerformanceMetrics)
    * [.MakeTime(date)](#Astronomy.MakeTime) ⇒ [<code>Time</code>](#Astronomy.Time)
    * [.Horizon(date, location, ra, dec)](#Astronomy.Horizon) ⇒ [<code>HorizontalCoordinates</code>](#Astronomy.HorizontalCoordinates)
    * [.MakeObserver(latitude_degrees, longitude_degrees, height_in_meters)](#Astronomy.MakeObserver)
    * [.SunPosition(date)](#Astronomy.SunPosition) ⇒ [<code>EclipticCoordinates</code>](#Astronomy.EclipticCoordinates)
    * [.SkyPos(gc_vector, observer)](#Astronomy.SkyPos) ⇒ <code>Astronomy.SkyCoordinates</code>
    * [.Ecliptic(gx, gy, gz)](#Astronomy.Ecliptic) ⇒ [<code>EclipticCoordinates</code>](#Astronomy.EclipticCoordinates)
    * [.GeoMoon(date)](#Astronomy.GeoMoon) ⇒ [<code>Vector</code>](#Astronomy.Vector)
    * [.HelioVector(body, date)](#Astronomy.HelioVector) ⇒ [<code>Vector</code>](#Astronomy.Vector)
    * [.GeoVector(body, date)](#Astronomy.GeoVector) ⇒ [<code>Vector</code>](#Astronomy.Vector)
    * [.Search(func, t1, t2, options)](#Astronomy.Search)
    * [.LongitudeFromSun(body, date)](#Astronomy.LongitudeFromSun) ⇒ <code>number</code>
    * [.AngleFromSun(body, date)](#Astronomy.AngleFromSun) ⇒ <code>number</code>
    * [.EclipticLongitude(body, date)](#Astronomy.EclipticLongitude) ⇒ <code>number</code>
    * [.Illumination(body, date)](#Astronomy.Illumination) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
    * [.Elongation(body)](#Astronomy.Elongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
    * [.SearchMaxElongation(body, startDate)](#Astronomy.SearchMaxElongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
    * [.SearchPeakMagnitude(body, startDate)](#Astronomy.SearchPeakMagnitude) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
    * [.ContinuousFunction](#Astronomy.ContinuousFunction) ⇒ <code>number</code>
    * [.SearchOptions](#Astronomy.SearchOptions) : <code>Object</code>


* * *

<a name="Astronomy.PerformanceInfo"></a>

### Astronomy.PerformanceInfo
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| search_func | <code>number</code> | Number of times [Search](#Astronomy.Search) called a <code>func</code> passed to it. |
| search | <code>number</code> | Number of times [Search](#Astronomy.Search) was called. |
| longitude_search | <code>number</code> | Number of times [Astronomy.SearchRelativeLongitude](Astronomy.SearchRelativeLongitude) was called. |
| longitude_iter | <code>number</code> | The total number of iterations executed inside [Astronomy.SearchRelativeLongitude](Astronomy.SearchRelativeLongitude). |


* [.PerformanceInfo](#Astronomy.PerformanceInfo)
    * [new PerformanceInfo()](#new_Astronomy.PerformanceInfo_new)
    * [.Clone()](#Astronomy.PerformanceInfo+Clone) ⇒ [<code>PerformanceInfo</code>](#Astronomy.PerformanceInfo)


* * *

<a name="new_Astronomy.PerformanceInfo_new"></a>

#### new PerformanceInfo()
Holds performance metrics for developers to optimize execution speed.
Most users can safely ignore this class.


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
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> | The JavaScript Date object for the given date and time.      This Date corresponds to the numeric day value stored in the <code>ut</code> property. |
| ut | <code>number</code> | Universal Time (UT1/UTC) in fractional days since the J2000 epoch.      Universal Time represents time measured with respect to the Earth's rotation,      tracking mean solar days.      The Astronomy library approximates UT1 and UTC as being the same thing.      This gives sufficient accuracy for the 1-arcminute angular resolution requirement      of this project. |
| tt | <code>number</code> | Terrestrial Time in fractional days since the J2000 epoch.      TT represents a continuously flowing ephemeris timescale independent of      any variations of the Earth's rotation, and is adjusted from UT      using historical and predictive models of those variations. |


* [.Time](#Astronomy.Time)
    * [new Time(date)](#new_Astronomy.Time_new)
    * [.toString()](#Astronomy.Time+toString) ⇒ <code>string</code>
    * [.AddDays(days)](#Astronomy.Time+AddDays) ⇒ [<code>Time</code>](#Astronomy.Time)


* * *

<a name="new_Astronomy.Time_new"></a>

#### new Time(date)
The date and time of an astronomical observation.
Objects of this type are used throughout the internals
of the Astronomy library, and are included in certain return objects.
The constructor is not accessible outside the Astronomy library;
outside users should call the [MakeTime](#Astronomy.MakeTime) function
to create a <code>Time</code> object.


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
| days | <code>number</code> | The floating point numbers of days by which to adjust the given date and time.      Positive values adjust the date toward the future, and      negative values adjust the date toward the past. |


* * *

<a name="Astronomy.Vector"></a>

### Astronomy.Vector
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| x | <code>number</code> | The x-coordinate expressed in astronomical units (AU). |
| y | <code>number</code> | The y-coordinate expressed in astronomical units (AU). |
| z | <code>number</code> | The z-coordinate expressed in astronomical units (AU). |
| t | [<code>Time</code>](#Astronomy.Time) | The time at which the vector is valid. |


* * *

<a name="new_Astronomy.Vector_new"></a>

#### new Vector()
Holds the Cartesian coordinates of a vector in 3D space,
along with the time at which the vector is valid.


* * *

<a name="Astronomy.EquatorialCoordinates"></a>

### Astronomy.EquatorialCoordinates
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| ra | <code>number</code> | Right ascension in sidereal hours: [0, 24). |
| dec | <code>number</code> | Declination in degrees: [-90, +90]. |
| dist | <code>number</code> | Distance to the celestial object expressed in       <a href="https://en.wikipedia.org/wiki/Astronomical_unit">Astronomical Units</a> (AU). |


* * *

<a name="new_Astronomy.EquatorialCoordinates_new"></a>

#### new EquatorialCoordinates()
Holds right ascension, declination, and distance of a celestial object.


* * *

<a name="Astronomy.HorizontalCoordinates"></a>

### Astronomy.HorizontalCoordinates
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| azimuth | <code>number</code> | An angle in degrees measured starting at north and increasing positively toward the east.      The value is in the range [0, 360).      North = 0, east = 90, south = 180, west = 270. |
| altitude | <code>number</code> | An angle in degrees above (positive) or below (negative) the horizon.      The value is in the range [-90, +90].      The altitude angle is optionally adjusted upward due to atmospheric refraction. |
| ra | <code>number</code> | The right ascension of the celestial body in sidereal hours.      The value is in the reange [0, 24).      If <code>altitude</code> was adjusted for atmospheric reaction, <code>ra</code>      is likewise adjusted. |
| dec | <code>number</code> | The declination of of the celestial body in degrees.      The value in the range [-90, +90].      If <code>altitude</code> was adjusted for atmospheric reaction, <code>dec</code>      is likewise adjusted. |


* * *

<a name="new_Astronomy.HorizontalCoordinates_new"></a>

#### new HorizontalCoordinates()
Holds azimuth (compass direction) and altitude (angle above/below the horizon)
of a celestial object as seen by an observer at a particular location on the Earth's surface.
Also holds right ascension and declination of the same object.
All of these coordinates are optionally adjusted for atmospheric refraction;
therefore the right ascension and declination values may not exactly match
those found inside a corresponding [EquatorialCoordinates](#Astronomy.EquatorialCoordinates) object.


* * *

<a name="Astronomy.EclipticCoordinates"></a>

### Astronomy.EclipticCoordinates
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| ex | <code>number</code> | The Cartesian x-coordinate of the body in astronomical units (AU).      The x-axis is oriented in the direction of the       <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>. |
| ey | <code>number</code> | The Cartesian y-coordinate of the body in astronomical units (AU).      The y-axis is oriented 90 degrees counterclockwise from the equinox,      as seen from above the Sun's north pole. |
| ez | <code>number</code> | The Cartesian z-coordinate of the body in astronomical units (AU).      The z-axis is oriented perpendicular to the mean plane of the Earth's orbit,      along the direction of the Sun's north pole. |
| elat | <code>number</code> | The ecliptic latitude of the body in degrees.      This is the angle north or south of the mean plane of the Earth's orbit.      The value is in the range [-90, +90].      Positive values are north and negative values are south. |
| elon | <code>number</code> | The ecliptic longitude of the body in degrees.      This is the angle measured counterclockwise around the mean plane      of the Earth's orbit, as seen from above the Sun's north pole.      This is the same direction that the Earth orbits around the Sun.      The angle is measured starting at 0 from the equinox and increases      up to 360 degrees. |


* * *

<a name="new_Astronomy.EclipticCoordinates_new"></a>

#### new EclipticCoordinates()
Holds ecliptic coordinates of a celestial body.
The origin and date of the coordinate system may vary depending on the caller's usage.
In general, ecliptic coordinates are measured with respect to the mean plane of the Earth's 
orbit around the Sun.
Includes Cartesian coordinates <code>(ex, ey, ez)</code> measured in 
<a href="https://en.wikipedia.org/wiki/Astronomical_unit">Astronomical Units</a> (AU)
and spherical coordinates <code>(elon, elat)</code> measured in degrees.


* * *

<a name="Astronomy.Observer"></a>

### Astronomy.Observer
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| latitude_degrees | <code>number</code> | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude_degrees | <code>number</code> | The observer's geographic longitude in degrees east of the prime meridian       passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height_in_meters | <code>number</code> | The observer's elevation above mean sea level, expressed in meters. |


* * *

<a name="new_Astronomy.Observer_new"></a>

#### new Observer()
Represents the geographic location of an observer on the surface of the Earth.


* * *

<a name="Astronomy.IlluminationInfo"></a>

### Astronomy.IlluminationInfo
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | The date and time pertaining to the other calculated values in this object. |
| mag | <code>number</code> | The <a href="https://en.wikipedia.org/wiki/Apparent_magnitude">apparent visual magnitude</a> of the celestial body. |
| phase_angle | <code>number</code> | The angle in degrees as seen from the center of the celestial body between the Sun and the Earth.      The value is always in the range 0 to 180.      The phase angle provides a measure of what fraction of the body's face appears       illuminated by the Sun as seen from the Earth.      When the observed body is the Sun, the <code>phase</code> property is set to 0,      although this has no physical meaning because the Sun emits, rather than reflects, light.      When the phase is near 0 degrees, the body appears "full".      When it is 90 degrees, the body appears "half full".       And when it is 180 degrees, the body appears "new" and is very difficult to see      because it is both dim and lost in the Sun's glare as seen from the Earth. |
| phase_fraction | <code>number</code> | The fraction of the body's face that is illuminated by the Sun, as seen from the Earth.      Calculated from <code>phase_angle</code> for convenience. |
| helio_dist | <code>number</code> | The distance between the center of the Sun and the center of the body in       <a href="https://en.wikipedia.org/wiki/Astronomical_unit">Astronomical Units</a> (AU). |
| geo_dist | <code>number</code> | The distance between the center of the Earth and the center of the body in AU. |
| gc | [<code>Vector</code>](#Astronomy.Vector) | Geocentric coordinates: the 3D vector from the center of the Earth to the center of the body.      The components are in expressed in AU and are oriented with respect to the J2000 equatorial plane. |
| hc | [<code>Vector</code>](#Astronomy.Vector) | Heliocentric coordinates: The 3D vector from the center of the Sun to the center of the body.      Like <code>gc</code>, <code>hc</code> is expressed in AU and oriented with respect      to the J2000 equatorial plane. |
| ring_tilt | <code>number</code> \| <code>null</code> | For Saturn, this is the angular tilt of the planet's rings in degrees away      from the line of sight from the Earth. When the value is near 0, the rings      appear edge-on from the Earth and are therefore difficult to see.      When <code>ring_tilt</code> approaches its maximum value (about 27 degrees),      the rings appear widest and brightest from the Earth.      Unlike the <a href="https://ssd.jpl.nasa.gov/horizons.cgi">JPL Horizons</a> online tool,       this library includes the effect of the ring tilt angle in the calculated value       for Saturn's visual magnitude.      For all bodies other than Saturn, the value of <code>ring_tilt</code> is <code>null</code>. |


* * *

<a name="new_Astronomy.IlluminationInfo_new"></a>

#### new IlluminationInfo()
Contains information about the apparent brightness and sunlit phase of a celestial object.


* * *

<a name="Astronomy.ElongationEvent"></a>

### Astronomy.ElongationEvent
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**See**: Astronomy.Elongation  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | When the event occurs. |
| visibility | <code>string</code> | Either <code>"morning"</code> or <code>"evening"</code>,       indicating when the body is most easily seen. |
| elongation | <code>number</code> | The angle in degrees, as seen from the center of the Earth,       of the apparent separation between the body and the Sun.      This angle is measured in 3D space and is not projected onto the ecliptic plane. |
| relative_longitude | <code>number</code> | The angle in degrees, as seen from the Sun, between the      observed body and the Earth. This value is always between      0 and 180. More precisely, <code>relative_longitude</code> is the absolute      value of the difference between the heliocentric ecliptic longitudes of      the centers of the observed body and the Earth. |


* * *

<a name="new_Astronomy.ElongationEvent_new"></a>

#### new ElongationEvent()
Represents the visibility of a planet or the Moon relative to the Sun.
Includes angular separation from the Sun and whether visibility is
best in the morning or the evening.


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
Resets the internal performance metrics back to zero values.
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

### Astronomy.Horizon(date, location, ra, dec) ⇒ [<code>HorizontalCoordinates</code>](#Astronomy.HorizontalCoordinates)
Given a date and time, a geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial object,
returns horizontal coordinates (azimuth and altitude angles) for that object
as seen by that observer.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to find horizontal coordinates. |
| location | [<code>Observer</code>](#Astronomy.Observer) | The location of the observer for which to find horizontal coordinates. |
| ra | <code>number</code> | Right ascension in sidereal hours of the celestial object,       referred to the mean equinox of date for the J2000 epoch. |
| dec | <code>number</code> | Declination in degrees of the celestial object,       referred to the mean equator of date for the J2000 epoch.      Positive values are north of the celestial equator and negative values are south. |


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

### Astronomy.SkyPos(gc_vector, observer) ⇒ <code>Astronomy.SkyCoordinates</code>
Returns topocentric equatorial coordinates (right ascension and declination)
simultaneously in two different systems: J2000 and true-equator-of-date.
<i>Topocentric</i> refers to a position as seen by an observer on the surface of the Earth.
This function corrects for
<a href="https://en.wikipedia.org/wiki/Parallax">parallax</a> 
of the object between a geocentric observer and a topocentric observer.
This is most significant for the Moon, because it is so close to the Earth.
However, it can have a small effect on the apparent positions of other bodies.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
**Returns**: <code>Astronomy.SkyCoordinates</code> - The topocentric coordinates of the body as adjusted for the given observer.  

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

### Astronomy.Search(func, t1, t2, options)
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

| Param | Type |
| --- | --- |
| func | [<code>ContinuousFunction</code>](#Astronomy.ContinuousFunction) | 
| t1 | [<code>Time</code>](#Astronomy.Time) | 
| t2 | [<code>Time</code>](#Astronomy.Time) | 
| options | <code>null</code> \| [<code>SearchOptions</code>](#Astronomy.SearchOptions) | 


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
| body | <code>string</code> | Currently only <code>"Venus"</code> is supported.      Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from Earth,      so peak magnitude events have little practical value for that planet.      The Moon reaches peak magnitude very close to full moon, which can be found using       [Astronomy.SearchMoonQuarter](Astronomy.SearchMoonQuarter) or [Astronomy.SearchMoonPhase](Astronomy.SearchMoonPhase).      The other planets reach peak magnitude very close to opposition,       which can be found using [Astronomy.SearchRelativeLongitude](Astronomy.SearchRelativeLongitude). |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which to find the next peak magnitude event. |


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

