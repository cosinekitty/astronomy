//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GravitySimulator](index.md)

# GravitySimulator

class [GravitySimulator](index.md)

A simulation of zero or more small bodies moving through the Solar System.

A gravity simulator object simulates a series of incremental time steps, calculating the movement of the Sun and planets around the Solar System Barycenter (SSB). It calculates the resulting gravitational forces on an arbitrary list of small bodies provided by the caller.

## Constructors

| | |
|---|---|
| [GravitySimulator](-gravity-simulator.md)<br>fun [GravitySimulator](-gravity-simulator.md)(originBody: [Body](../-body/index.md), time: [Time](../-time/index.md), bodyStates: [List](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.collections/-list/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;)<br>Creates a gravity simulation object. |

## Functions

| Name | Summary |
|---|---|
| [solarSystemBodyState](solar-system-body-state.md)<br>fun [solarSystemBodyState](solar-system-body-state.md)(body: [Body](../-body/index.md)): [StateVector](../-state-vector/index.md)<br>Get the position and velocity of a Solar System body included in the simulation. |
| [swap](swap.md)<br>fun [swap](swap.md)()<br>Exchange the current time step with the previous time step. |
| [update](update.md)<br>fun [update](update.md)(time: [Time](../-time/index.md)): [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;<br>Advances a gravity simulation by a small time step. |

## Properties

| Name | Summary |
|---|---|
| [originBody](origin-body.md)<br>val [originBody](origin-body.md): [Body](../-body/index.md)<br>The origin of the reference frame. See constructor for more info. |
