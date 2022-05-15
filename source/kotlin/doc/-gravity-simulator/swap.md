//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GravitySimulator](index.md)/[swap](swap.md)

# swap

fun [swap](swap.md)()

Exchange the current time step with the previous time step.

Sometimes it is helpful to "explore" various times near a given simulation time step, while repeatedly returning to the original time step. For example, when backdating a position for light travel time, the caller may wish to repeatedly try different amounts of backdating. When the backdating solver has converged, the caller wants to leave the simulation in its original state.

This function allows a single "undo" of a simulation, and does so very efficiently.

Usually this function will be called immediately after a matching call to [GravitySimulator.update](update.md). It has the effect of rolling back the most recent update. If called twice in a row, it reverts the swap and thus has no net effect.

The constructor initializes the current state and previous state to be identical. Both states represent the time parameter that was passed into the constructor. Therefore, swap will have no effect from the caller's point of view when passed a simulator that has not yet been updated by a call to [GravitySimulator.update](update.md).
