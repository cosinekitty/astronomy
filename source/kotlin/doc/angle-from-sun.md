//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[angleFromSun](angle-from-sun.md)

# angleFromSun

[jvm]\
fun [angleFromSun](angle-from-sun.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Returns the angle between the given body and the Sun, as seen from the Earth.

This function calculates the angular separation between the given body and the Sun, as seen from the center of the Earth. This angle is helpful for determining how easy it is to see the body away from the glare of the Sun.

#### Return

The angle in degrees between the Sun and the specified body as seen from the center of the Earth.

## Parameters

jvm

| | |
|---|---|
| body | The celestial body whose angle from the Sun is to be measured.     Not allowed to be `Body.Earth`. |
| time | The time at which the observation is made. |
