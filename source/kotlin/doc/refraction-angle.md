//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[refractionAngle](refraction-angle.md)

# refractionAngle

[jvm]\
fun [refractionAngle](refraction-angle.md)(refraction: [Refraction](-refraction/index.md), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.

Given an altitude angle and a refraction option, calculates the amount of "lift" caused by atmospheric refraction. This is the number of degrees higher in the sky an object appears due to the lensing of the Earth's atmosphere.

## Parameters

jvm

| | |
|---|---|
| refraction | The option selecting which refraction correction to use.     If `Refraction.Normal`, uses a well-behaved refraction model that works well for     all valid values (-90 to +90) of `altitude`.     If `Refraction.JplHor`, this function returns a compatible value with the JPL Horizons tool.     If any other value, including `Refraction.None`, this function returns 0. |
| altitude | An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90. |
