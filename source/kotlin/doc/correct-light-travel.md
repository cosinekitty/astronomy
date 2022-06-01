//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[correctLightTravel](correct-light-travel.md)

# correctLightTravel

fun [correctLightTravel](correct-light-travel.md)(func: [PositionFunction](-position-function/index.md), time: [Time](-time/index.md)): [Vector](-vector/index.md)

Solve for light travel time of a vector function.

When observing a distant object, for example Jupiter as seen from Earth, the amount of time it takes for light to travel from the object to the observer can significantly affect the object's apparent position. This function is a generic solver that figures out how long in the past light must have left the observed object to reach the observer at the specified observation time. It uses [PositionFunction](-position-function/index.md) to express an arbitrary position vector as a function of time.

This function repeatedly calls func.Position, passing a series of time estimates in the past. Then func.Position must return a relative state vector between the observer and the target. correctLightTravel keeps calling func.Position with more and more refined estimates of the time light must have left the target to arrive at the observer.

For common use cases, it is simpler to use [backdatePosition](backdate-position.md) for calculating the light travel time correction of one body observing another body.

#### Return

The position vector at the solved backdated time. The t field holds the time that light left the observed body to arrive at the observer at the observation time.

## Parameters

| | |
|---|---|
| func | An arbitrary position vector as a function of time. |
| time | The observation time for which to solve for light travel delay. |
