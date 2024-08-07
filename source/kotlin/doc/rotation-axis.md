//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationAxis](rotation-axis.md)

# rotationAxis

fun [rotationAxis](rotation-axis.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [AxisInfo](-axis-info/index.md)

Calculates information about a body's rotation axis at a given time.

Calculates the orientation of a body's rotation axis, along with the rotation angle of its prime meridian, at a given moment in time.

This function uses formulas standardized by the IAU Working Group on Cartographics and Rotational Elements 2015 report, as described in the following document:

https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf

See [AxisInfo](-axis-info/index.md) for more detailed information.

#### Return

North pole orientation and body spin angle.

## Parameters

| | |
|---|---|
| body | One of the following values: [Body.Sun](-body/-sun/index.md), [Body.Moon](-body/-moon/index.md), [Body.Mercury](-body/-mercury/index.md), [Body.Venus](-body/-venus/index.md), [Body.Earth](-body/-earth/index.md), [Body.Mars](-body/-mars/index.md), [Body.Jupiter](-body/-jupiter/index.md), [Body.Saturn](-body/-saturn/index.md), [Body.Uranus](-body/-uranus/index.md), [Body.Neptune](-body/-neptune/index.md), [Body.Pluto](-body/-pluto/index.md). |
| time | The time at which to calculate the body's rotation axis. |
