# API Reference
<a name="Astronomy"></a>

## Astronomy : <code>object</code>
**Kind**: global namespace  

* [Astronomy](#Astronomy) : <code>object</code>
    * [.Time](#Astronomy.Time)
        * [new Time(date)](#new_Astronomy.Time_new)
        * [.toString()](#Astronomy.Time+toString) ⇒ <code>string</code>
        * [.AddDays()](#Astronomy.Time+AddDays) ⇒ [<code>Time</code>](#Astronomy.Time)
    * [.Observer](#Astronomy.Observer)
        * [new Observer()](#new_Astronomy.Observer_new)
    * [.IlluminationInfo](#Astronomy.IlluminationInfo)
        * [new IlluminationInfo()](#new_Astronomy.IlluminationInfo_new)
    * [.ElongationEvent](#Astronomy.ElongationEvent)
        * [new ElongationEvent()](#new_Astronomy.ElongationEvent_new)
    * [.Bodies](#Astronomy.Bodies) : <code>Array.&lt;string&gt;</code>
    * [.MakeTime(date)](#Astronomy.MakeTime) ⇒ [<code>Time</code>](#Astronomy.Time)
    * [.Horizon(date, location, ra, dec)](#Astronomy.Horizon) ⇒ <code>Astronomy.HorizontalCoordinates</code>
    * [.MakeObserver(latitude_degrees, longitude_degrees, height_in_meters)](#Astronomy.MakeObserver)
    * [.SunPosition()](#Astronomy.SunPosition)
    * [.SkyPos()](#Astronomy.SkyPos)
    * [.Ecliptic()](#Astronomy.Ecliptic)
    * [.GeoMoon(The)](#Astronomy.GeoMoon) ⇒ <code>Astronomy.Vector</code>
    * [.Illumination(body, date)](#Astronomy.Illumination) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
    * [.Elongation(body)](#Astronomy.Elongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
    * [.SearchMaxElongation(body, startDate)](#Astronomy.SearchMaxElongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
    * [.SearchPeakMagnitude(body, startDate)](#Astronomy.SearchPeakMagnitude) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)

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
    * [.AddDays()](#Astronomy.Time+AddDays) ⇒ [<code>Time</code>](#Astronomy.Time)

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

<a name="Astronomy.Time+toString"></a>

#### time.toString() ⇒ <code>string</code>
Formats a <code>Time</code> object as an 
<a href="https://en.wikipedia.org/wiki/ISO_8601">ISO 8601</a>
date/time string in UTC, to millisecond resolution.
Example: <code>2018-08-17T17:22:04.050Z</code>

**Kind**: instance method of [<code>Time</code>](#Astronomy.Time)  
<a name="Astronomy.Time+AddDays"></a>

#### time.AddDays() ⇒ [<code>Time</code>](#Astronomy.Time)
Returns a new <code>Time</code> object adjusted by the floating point number of days.
Does NOT modify the original <code>Time</code> object.

**Kind**: instance method of [<code>Time</code>](#Astronomy.Time)  
<a name="Astronomy.Observer"></a>

### Astronomy.Observer
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| latitude_degrees | <code>Number</code> | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude_degrees | <code>Number</code> | The observer's geographic longitude in degrees east of the prime meridian       passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height_in_meters | <code>Number</code> | The observer's elevation above mean sea level, expressed in meters. |

<a name="new_Astronomy.Observer_new"></a>

#### new Observer()
Represents the geographic location of an observer on the surface of the Earth.

<a name="Astronomy.IlluminationInfo"></a>

### Astronomy.IlluminationInfo
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | The date and time pertaining to the other calculated values in this object. |
| mag | <code>Number</code> | The <a href="https://en.wikipedia.org/wiki/Apparent_magnitude">apparent visual magnitude</a> of the celestial body. |
| phase | <code>Number</code> | The angle in degrees as seen from the center of the celestial body between the Sun and the Earth.      The value is always in the range 0 to 180.      The phase angle provides a measure of what fraction of the body's face appears       illuminated by the Sun as seen from the Earth.      When the observed body is the Sun, the <code>phase</code> property is set to 0,      although this has no physical meaning because the Sun emits, rather than reflects, light.      To calculate the illuminated fraction, use the formula $f = \frac{1}{2} \left( 1 + cos(\phi) \right)$      where $f$ is the illuminated fraction and $\phi$ is the phase angle.      When the phase is near 0 degrees, the body appears "full".      When it is 90 degrees, the body appears "half full".       And when it is 180 degrees, the body appears "new" and is very difficult to see      because it is both dim and lost in the Sun's glare as seen from the Earth. |
| helio_dist | <code>Number</code> | The distance between the center of the Sun and the center of the body in       <a href="https://en.wikipedia.org/wiki/Astronomical_unit">Astronomical Units</a> (AU). |
| geo_dist | <code>Number</code> | The distance between the center of the Earth and the center of the body in AU. |
| gc | <code>Astronomy.Vector</code> | Geocentric coordinates: the 3D vector from the center of the Earth to the center of the body.      The components are in expressed in AU and the oriented with respect to the J2000 equatorial plane. |
| hc | <code>Astronomy.Vector</code> | Heliocentric coordinates: The 3D vector from the center of the Sun to the center of the body.      Like <code>gc</code>, <code>hc</code> is expressed in AU and oriented with respect      to the J2000 equatorial plane. |
| ring_tilt | <code>Number</code> \| <code>null</code> | For Saturn, this is the angular tilt of the planet's rings in degrees away      from the line of sight from the Earth. When the value is near 0, the rings      appear edge-on from the Earth and are therefore difficult to see.      When <code>ring_tilt</code> approaches its maximum value (about 27 degrees),      the rings appear widest and brightest from the Earth.      Unlike the <a href="https://ssd.jpl.nasa.gov/horizons.cgi">JPL Horizons</a> online tool,       this library includes the effect of the ring tilt angle in the calculated value       for Saturn's visual magnitude.      For all bodies other than Saturn, the value of <code>ring_tilt</code> is <code>null</code>. |

<a name="new_Astronomy.IlluminationInfo_new"></a>

#### new IlluminationInfo()
Contains information about the apparent brightness and sunlit phase of a celestial object.

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

<a name="new_Astronomy.ElongationEvent_new"></a>

#### new ElongationEvent()
Represents the visibility of a planet or the Moon relative to the Sun.
Includes angular separation from the Sun and whether visibility is
best in the morning or the evening.

<a name="Astronomy.Bodies"></a>

### Astronomy.Bodies : <code>Array.&lt;string&gt;</code>
An array of strings, each a name of a supported astronomical body.
     Not all bodies are valid for all functions, but any string not in this
     list is not supported at all.

**Kind**: static constant of [<code>Astronomy</code>](#Astronomy)  
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

<a name="Astronomy.Horizon"></a>

### Astronomy.Horizon(date, location, ra, dec) ⇒ <code>Astronomy.HorizontalCoordinates</code>
Given a date and time, a geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial object,
returns horizontal coordinates (azimuth and altitude angles) for that object
as seen by that observer.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>Number</code> \| [<code>Time</code>](#Astronomy.Time) |  |
| location | [<code>Observer</code>](#Astronomy.Observer) | The location of the observer for which to find horizontal coordinates. |
| ra | <code>Number</code> | Right ascension in sidereal hours of the celestial object,       referred to the mean equinox of date for the J2000 epoch. |
| dec | <code>Number</code> | Declination in degrees of the celestial object,       referred to the mean equator of date for the J2000 epoch.      Positive values are north of the celestial equator and negative values are south. |

<a name="Astronomy.MakeObserver"></a>

### Astronomy.MakeObserver(latitude_degrees, longitude_degrees, height_in_meters)
Creates an [Observer](#Astronomy.Observer) object.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| latitude_degrees | <code>Number</code> | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude_degrees | <code>Number</code> | The observer's geographic longitude in degrees east of the prime meridian       passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height_in_meters | <code>Number</code> | The observer's elevation above mean sea level, expressed in meters.      If omitted, the elevation is assumed to be 0 meters. |

<a name="Astronomy.SunPosition"></a>

### Astronomy.SunPosition()
Returns apparent geocentric true ecliptic coordinates of date for the Sun.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
<a name="Astronomy.SkyPos"></a>

### Astronomy.SkyPos()
Returns equatorial coordinates (right ascension and declination)
in two different systems: J2000 and true-equator-of-date.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
<a name="Astronomy.Ecliptic"></a>

### Astronomy.Ecliptic()
Given J2000 equatorial Cartesian coordinates, 
returns J2000 ecliptic latitude, longitude, and cartesian coordinates.
You can call [Astronomy.GeoVector](Astronomy.GeoVector) and use its (x, y, z) return values
to pass into this function.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
<a name="Astronomy.GeoMoon"></a>

### Astronomy.GeoMoon(The) ⇒ <code>Astronomy.Vector</code>
Calculates the geocentric Cartesian coordinates for the Moon in the J2000 equatorial system.
Based on the Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
which in turn derives from E. W. Brown's lunar theories.
Adapted from Turbo Pascal code from the book 
<a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a> 
by Montenbruck and Pfleger.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| The | <code>Date</code> \| <code>Number</code> \| [<code>Time</code>](#Astronomy.Time) | date and time for which to calculate the Moon's geocentric position. |

<a name="Astronomy.Illumination"></a>

### Astronomy.Illumination(body, date) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
Calculates the phase angle, visual maginitude, 
and other values relating to the body's illumination
at the given date and time, as seen from the Earth.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the celestial body being observed.      Not allowed to be <code>"Earth"</code>. |
| date | <code>Date</code> \| <code>Number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time for which to calculate the illumination data for the given body. |

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
| body | <code>string</code> | either "Mercury" or "Venus" |
| startDate | <code>Date</code> | the date and time after which to search for the next maximum elongation event |

<a name="Astronomy.SearchPeakMagnitude"></a>

### Astronomy.SearchPeakMagnitude(body, startDate) ⇒ [<code>IlluminationInfo</code>](#Astronomy.IlluminationInfo)
Searches for the date and time Venus will next appear brightest as seen from the Earth.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | Currently only <code>"Venus"</code> is supported.      Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from Earth,      so peak magnitude events have little practical value for this planet.      The Moon reaches peak magnitude very close to full moon, which can be found using       [Astronomy.SearchMoonQuarter](Astronomy.SearchMoonQuarter) or [Astronomy.SearchMoonPhase](Astronomy.SearchMoonPhase).      The other planets reach peak magnitude very close to opposition,       which can be found using [Astronomy.SearchRelativeLongitude](Astronomy.SearchRelativeLongitude). |
| startDate | <code>Date</code> \| <code>Number</code> \| [<code>Time</code>](#Astronomy.Time) | The date and time after which to find the next peak magnitude event. |

