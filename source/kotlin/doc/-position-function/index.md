//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[PositionFunction](index.md)

# PositionFunction

fun interface [PositionFunction](index.md)

A function for which to solve a light-travel time problem.

The function [correctLightTravel](../correct-light-travel.md) solves a generalized problem of deducing how far in the past light must have left a target object to be seen by an observer at a specified time. This interface expresses an arbitrary position vector as function of time that is passed to [correctLightTravel](../correct-light-travel.md).

## Functions

| Name | Summary |
|---|---|
| [position](position.md)<br>abstract fun [position](position.md)(time: [Time](../-time/index.md)): [Vector](../-vector/index.md)<br>Returns a relative position vector for a given time. |
