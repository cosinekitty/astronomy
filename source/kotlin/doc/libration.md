//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[libration](libration.md)

# libration

fun [libration](libration.md)(time: [Time](-time/index.md)): [LibrationInfo](-libration-info/index.md)

Calculates the Moon's libration angles at a given moment in time.

Libration is an observed back-and-forth wobble of the portion of the Moon visible from the Earth. It is caused by the imperfect tidal locking of the Moon's fixed rotation rate, compared to its variable angular speed of orbit around the Earth.

This function calculates a pair of perpendicular libration angles, one representing rotation of the Moon in eclitpic longitude elon, the other in ecliptic latitude elat, both relative to the Moon's mean Earth-facing position.

This function also returns the geocentric position of the Moon expressed in ecliptic longitude mlon, ecliptic latitude mlat, the distance dist_km between the centers of the Earth and Moon expressed in kilometers, and the apparent angular diameter of the Moon diam_deg.

#### Return

The Moon's ecliptic position and libration angles as seen from the Earth.

## Parameters

jvm

| | |
|---|---|
| time | The date and time for which to calculate lunar libration. |
