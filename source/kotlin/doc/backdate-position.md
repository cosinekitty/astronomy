//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[backdatePosition](backdate-position.md)

# backdatePosition

fun [backdatePosition](backdate-position.md)(time: [Time](-time/index.md), observerBody: [Body](-body/index.md), targetBody: [Body](-body/index.md), aberration: [Aberration](-aberration/index.md)): [Vector](-vector/index.md)

Solve for light travel time correction of apparent position.

When observing a distant object, for example Jupiter as seen from Earth, the amount of time it takes for light to travel from the object to the observer can significantly affect the object's apparent position.

This function solves the light travel time correction for the apparent relative position vector of a target body as seen by an observer body at a given observation time.

For geocentric calculations, #geoVector also includes light travel time correction, but the time t embedded in its returned vector refers to the observation time, not the backdated time that light left the observed body. Thus backdatePosition provides direct access to the light departure time for callers that need it.

For a more generalized light travel correction solver, see [correctLightTravel](correct-light-travel.md).

#### Return

The position vector at the solved backdated time. Its t field holds the time that light left the observed body to arrive at the observer at the observation time.

## Parameters

| | |
|---|---|
| time | The time of observation. |
| observerBody | The body to be used as the observation location. |
| targetBody | The body to be observed. |
| aberration | Aberration.Corrected to correct for aberration, or Aberration.None to leave uncorrected. |
