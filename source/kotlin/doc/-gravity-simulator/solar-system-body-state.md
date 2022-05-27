//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GravitySimulator](index.md)/[solarSystemBodyState](solar-system-body-state.md)

# solarSystemBodyState

fun [solarSystemBodyState](solar-system-body-state.md)(body: [Body](../-body/index.md)): [StateVector](../-state-vector/index.md)

Get the position and velocity of a Solar System body included in the simulation.

In order to simulate the movement of small bodies through the Solar System, the simulator needs to calculate the state vectors for the Sun and planets.

If an application wants to know the positions of one or more of the planets in addition to the small bodies, this function provides a way to obtain their state vectors. This is provided for the sake of efficiency, to avoid redundant calculations.

The state vector is returned relative to the position and velocity of the originBody parameter that was passed to this object's constructor.

## Parameters

| | |
|---|---|
| body | The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune. |
