//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[inverseRefractionAngle](inverse-refraction-angle.md)

# inverseRefractionAngle

fun [inverseRefractionAngle](inverse-refraction-angle.md)(refraction: [Refraction](-refraction/index.md), bentAltitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Calculates the inverse of an atmospheric refraction angle.

Given an observed altitude angle that includes atmospheric refraction, calculates the negative angular correction to obtain the unrefracted altitude. This is useful for cases where observed horizontal coordinates are to be converted to another orientation system, but refraction first must be removed from the observed position.

## Parameters

| | |
|---|---|
| refraction | The option selecting which refraction correction to use. |
| bentAltitude | The apparent altitude that includes atmospheric refraction. |
|  | The angular adjustment in degrees to be added to the     altitude angle to remove atmospheric lensing.     This will be less than or equal to zero. |
