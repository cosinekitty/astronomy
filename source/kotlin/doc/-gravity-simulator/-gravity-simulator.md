//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GravitySimulator](index.md)/[GravitySimulator](-gravity-simulator.md)

# GravitySimulator

fun [GravitySimulator](-gravity-simulator.md)(originBody: [Body](../-body/index.md), time: [Time](../-time/index.md), bodyStates: [List](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin.collections/-list/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;)

Creates a gravity simulation object.

## Parameters

| | |
|---|---|
| originBody | Specifies the origin of the reference frame. All position vectors and velocity vectors will use originBody as the origin of the coordinate system. This origin applies to all the input vectors provided in the bodyStates parameter of this function, along with all output vectors returned by [GravitySimulator.update](update.md). Most callers will want to provide one of the following: [Body.Sun](../-body/-sun/index.md) for heliocentric coordinates, [Body.SSB](../-body/-s-s-b/index.md) for solar system barycentric coordinates, or [Body.Earth](../-body/-earth/index.md) for geocentric coordinates. Note that the gravity simulator does not correct for light travel time; all state vectors are tied to a Newtonian "instantaneous" time. |
| time | The initial time at which to start the simulation. |
| bodyStates | An array of initial state vectors (positions and velocities) of the small bodies to be simulated. The caller must know the positions and velocities of the small bodies at an initial moment in time. Their positions and velocities are expressed with respect to originBody, using equatorial J2000 orientation (EQJ). Positions are expressed in astronomical units (AU). Velocities are expressed in AU/day. All the times embedded within the state vectors must be exactly equal to time, or this constructor will throw an exception. |
