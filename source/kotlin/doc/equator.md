//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[equator](equator.md)

# equator

fun [equator](equator.md)(body: [Body](-body/index.md), time: [Time](-time/index.md), observer: [Observer](-observer/index.md), equdate: [EquatorEpoch](-equator-epoch/index.md), aberration: [Aberration](-aberration/index.md)): [Equatorial](-equatorial/index.md)

Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.

Calculates topocentric equatorial coordinates in one of two different systems: J2000 or true-equator-of-date, depending on the value of the equdate parameter. Equatorial coordinates include right ascension, declination, and distance in astronomical units.

This function corrects for light travel time: it adjusts the apparent location of the observed body based on how long it takes for light to travel from the body to the Earth.

This function corrects for *topocentric parallax*, meaning that it adjusts for the angular shift depending on where the observer is located on the Earth. This is most significant for the Moon, because it is so close to the Earth. However, parallax corection has a small effect on the apparent positions of other bodies.

Correction for aberration is optional, using the aberration parameter.

#### Return

Topocentric equatorial coordinates of the celestial body.

## Parameters

| | |
|---|---|
| body | The celestial body to be observed. Not allowed to be [Body.Earth]. |
| time | The date and time at which the observation takes place. |
| observer | A location on or near the surface of the Earth. |
| equdate | Selects the date of the Earth's equator in which to express the equatorial coordinates. |
| aberration | Selects whether or not to correct for aberration. |
