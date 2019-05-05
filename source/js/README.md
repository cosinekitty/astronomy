# API Reference
<a name="Astronomy"></a>

## Astronomy : <code>object</code>
**Kind**: global namespace  

* [Astronomy](#Astronomy) : <code>object</code>
    * [.Time](#Astronomy.Time)
        * [new Time(date)](#new_Astronomy.Time_new)
        * [.toString()](#Astronomy.Time+toString) ⇒ <code>string</code>
        * [.AddDays()](#Astronomy.Time+AddDays) ⇒ [<code>Time</code>](#Astronomy.Time)
    * [.ElongationEvent](#Astronomy.ElongationEvent)
        * [new ElongationEvent()](#new_Astronomy.ElongationEvent_new)
    * [.Bodies](#Astronomy.Bodies) : <code>Array.&lt;string&gt;</code>
    * [.MakeTime(date)](#Astronomy.MakeTime) ⇒ [<code>Time</code>](#Astronomy.Time)
    * [.Horizon()](#Astronomy.Horizon)
    * [.Elongation(body)](#Astronomy.Elongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
    * [.SearchMaxElongation(body, startDate)](#Astronomy.SearchMaxElongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)

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
<a name="Astronomy.ElongationEvent"></a>

### Astronomy.ElongationEvent
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**See**: Astronomy.Elongation  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>Time</code>](#Astronomy.Time) | When the event occurs. |
| visibility | <code>string</code> | Either "morning" or "evening", indicating when the body is most easily seen. |
| elongation | <code>number</code> | The angle in degrees, as seen from the center of the Earth,       of the apparent separation between the body and the Sun.      This angle is measured in 3D space and is not projected onto the ecliptic plane. |
| relative_longitude | <code>number</code> | The angle in degrees, as seen from the Sun, between the      observed body and the Earth. This value is always between      0 and 180. More precisely, relative_longitude is the absolute      value of the difference between the heliocentric ecliptic longitudes of      the centers of the observed body and the Earth. |

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

### Astronomy.Horizon()
Given a date and time, a geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial object,
returns horizontal coordinates (azimuth and altitude angles) for that object
as seen by that observer.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  
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
| body | <code>string</code> | The name of the observed body. Not allowed to be "Earth". |

<a name="Astronomy.SearchMaxElongation"></a>

### Astronomy.SearchMaxElongation(body, startDate) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
Searches for the next maximum elongation event for Mercury or Venus 
that occurs after the given start date. Calling with other values 
of 'body' will result in an exception. 
Maximum elongation occurs when the body has the greatest
angular separation from the Sun, as seen from the Earth.
Returns an object containing the date and time of the next
maximum elongation, the elongation in degrees, and whether
the body is visible in the morning or evening.

**Kind**: static method of [<code>Astronomy</code>](#Astronomy)  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | either "Mercury" or "Venus" |
| startDate | <code>Date</code> | the date and time after which to search for the next maximum elongation event |

