//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[StateVector](index.md)

# StateVector

[jvm]\
data class [StateVector](index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vx: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vy: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vz: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](../-astro-time/index.md))

## Constructors

| | |
|---|---|
| [StateVector](-state-vector.md) | [jvm]<br>fun [StateVector](-state-vector.md)(pos: [AstroVector](../-astro-vector/index.md), vel: [AstroVector](../-astro-vector/index.md), time: [AstroTime](../-astro-time/index.md))<br>Combines a position vector and a velocity vector into a single state vector. |

## Functions

| Name | Summary |
|---|---|
| [position](position.md) | [jvm]<br>fun [position](position.md)(): [AstroVector](../-astro-vector/index.md)<br>Returns the position vector associated with this state vector. |
| [velocity](velocity.md) | [jvm]<br>fun [velocity](velocity.md)(): [AstroVector](../-astro-vector/index.md)<br>Returns the velocity vector associated with this state vector. |

## Properties

| Name | Summary |
|---|---|
| [t](t.md) | [jvm]<br>val [t](t.md): [AstroTime](../-astro-time/index.md)<br>The date and time at which this vector is valid. |
| [vx](vx.md) | [jvm]<br>val [vx](vx.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian velocity x-component expressed in AU/day. |
| [vy](vy.md) | [jvm]<br>val [vy](vy.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian velocity y-component expressed in AU/day. |
| [vz](vz.md) | [jvm]<br>val [vz](vz.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian velocity z-component expressed in AU/day. |
| [x](x.md) | [jvm]<br>val [x](x.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian position x-coordinate expressed in AU. |
| [y](y.md) | [jvm]<br>val [y](y.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian position y-coordinate expressed in AU. |
| [z](z.md) | [jvm]<br>val [z](z.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian position z-coordinate expressed in AU. |
