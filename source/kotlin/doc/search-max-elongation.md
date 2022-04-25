//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchMaxElongation](search-max-elongation.md)

# searchMaxElongation

fun [searchMaxElongation](search-max-elongation.md)(body: [Body](-body/index.md), startTime: [Time](-time/index.md)): [ElongationInfo](-elongation-info/index.md)

Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.

Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is. Mercury especially is almost always impossible to see because it gets lost in the Sun's glare. The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through a telescope without atmospheric interference, are when these planets reach maximum elongation. These are events where the planets reach the maximum angle from the Sun as seen from the Earth.

This function solves for those times, reporting the next maximum elongation event's date and time, the elongation value itself, the relative longitude with the Sun, and whether the planet is best observed in the morning or evening. See [elongation](elongation.md) for more details about the returned structure.

## Parameters

jvm

| | |
|---|---|
| body | Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception.     To find the best viewing opportunites for planets farther from the Sun than the Earth is (Mars through Pluto)     use [searchRelativeLongitude] to find the next opposition event. |
| startTime | The date and time at which to begin the search. The maximum elongation event found will always     be the first one that occurs after this date and time. |
