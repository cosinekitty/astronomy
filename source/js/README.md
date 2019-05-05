# API Reference
<a name="Astronomy"></a>

## Astronomy : <code>object</code>
**Kind**: global namespace  

* [Astronomy](#Astronomy) : <code>object</code>
    * [.ElongationEvent](#Astronomy.ElongationEvent)
        * [new ElongationEvent()](#new_Astronomy.ElongationEvent_new)
    * [.Elongation(body)](#Astronomy.Elongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)
    * [.SearchMaxElongation(body, startDate)](#Astronomy.SearchMaxElongation) ⇒ [<code>ElongationEvent</code>](#Astronomy.ElongationEvent)

<a name="Astronomy.ElongationEvent"></a>

### Astronomy.ElongationEvent
**Kind**: static class of [<code>Astronomy</code>](#Astronomy)  
**See**: Astronomy.Elongation  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | <code>Time</code> | When the event occurs. |
| visibility | <code>string</code> | Either "morning" or "evening", indicating when the body is most easily seen. |
| elongation | <code>number</code> | The angle in degrees, as seen from the center of the Earth,       of the apparent separation between the body and the Sun.      This angle is measured in 3D space and is not projected onto the ecliptic plane. |
| relative_longitude | <code>number</code> | The angle in degrees, as seen from the Sun, between the      observed body and the Earth. This value is always between      0 and 180. More precisely, relative_longitude is the absolute      value of the difference between the heliocentric ecliptic longitudes of      the centers of observed body and the Earth. |

<a name="new_Astronomy.ElongationEvent_new"></a>

#### new ElongationEvent()
Represents the visibility of a planet or the Moon relative to the Sun.
Includes angular separation from the Sun and whether visibility is
best in the morning or the evening.

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

