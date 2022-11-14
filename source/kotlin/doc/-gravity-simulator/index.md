//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GravitySimulator](index.md)

# GravitySimulator

class [GravitySimulator](index.md)

A simulation of zero or more small bodies moving through the Solar System.

This class calculates the movement of arbitrary small bodies, such as asteroids or comets, that move through the Solar System. It does so by calculating the gravitational forces on the bodies from the Sun and planets. The user of this class supplies a list of initial positions and velocities for the small bodies. Then the class can update the positions and velocities over small time steps. The gravity simulator also provides access to the positions and velocities of the Sun and planets used in the simulation.

## Constructors

| | |
|---|---|
| [GravitySimulator](-gravity-simulator.md)<br>fun [GravitySimulator](-gravity-simulator.md)(originBody: [Body](../-body/index.md), time: [Time](../-time/index.md), bodyStates: List&lt;[StateVector](../-state-vector/index.md)&gt;)<br>Creates a gravity simulation object. |

## Functions

| Name | Summary |
|---|---|
| [solarSystemBodyState](solar-system-body-state.md)<br>fun [solarSystemBodyState](solar-system-body-state.md)(body: [Body](../-body/index.md)): [StateVector](../-state-vector/index.md)<br>Get the position and velocity of a Solar System body included in the simulation. |
| [swap](swap.md)<br>fun [swap](swap.md)()<br>Exchange the current time step with the previous time step. |
| [time](time.md)<br>fun [time](time.md)(): [Time](../-time/index.md)<br>Returns the time of the current simulation step. |
| [update](update.md)<br>fun [update](update.md)(time: [Time](../-time/index.md)): Array&lt;[StateVector](../-state-vector/index.md)&gt;<br>Advances the gravity simulation by a small time step. |

## Properties

| Name | Summary |
|---|---|
| [originBody](origin-body.md)<br>val [originBody](origin-body.md): [Body](../-body/index.md)<br>The origin of the reference frame. See constructor for more info. |
