//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationAxis](rotation-axis.md)

# rotationAxis

[jvm]\
fun [rotationAxis](rotation-axis.md)(body: [Body](../-body/index.md), time: [AstroTime](../-astro-time/index.md)): [AxisInfo](../-axis-info/index.md)

Calculates information about a body's rotation axis at a given time.

Calculates the orientation of a body's rotation axis, along with the rotation angle of its prime meridian, at a given moment in time.

This function uses formulas standardized by the IAU Working Group on Cartographics and Rotational Elements 2015 report, as described in the following document:

https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf

See [AxisInfo](../-axis-info/index.md) for more detailed information.

#### Return

North pole orientation and body spin angle.

## Parameters

jvm

| | |
|---|---|
| body | One of the following values:     `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`,     `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. |
| time | The time at which to calculate the body's rotation axis. |
