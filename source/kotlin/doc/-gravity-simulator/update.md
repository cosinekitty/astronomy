//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GravitySimulator](index.md)/[update](update.md)

# update

fun [update](update.md)(time: [Time](../-time/index.md)): [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;

Advances the gravity simulation by a small time step.

Updates the simulation of the user-supplied small bodies to the time indicated by the time parameter. Retuns an updated array of state vectors for the small bodies. The positions and velocities in the returned array are referenced to the originBody that was used to construct this simulator.

## Parameters

| | |
|---|---|
| time | A time that is a small increment away from the current simulation time. It is up to the developer to figure out an appropriate time increment. Depending on the trajectories, a smaller or larger increment may be needed for the desired accuracy. Some experimentation may be needed. Generally, bodies that stay in the outer Solar System and move slowly can use larger time steps.  Bodies that pass into the inner Solar System and move faster will need a smaller time step to maintain accuracy. The time value may be after or before the current simulation time to move forward or backward in time. |
