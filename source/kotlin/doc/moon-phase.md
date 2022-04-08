//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[moonPhase](moon-phase.md)

# moonPhase

[jvm]\
fun [moonPhase](moon-phase.md)(time: [AstroTime](-astro-time/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Returns the Moon's phase as an angle from 0 to 360 degrees.

This function determines the phase of the Moon using its apparent ecliptic longitude relative to the Sun, as seen from the center of the Earth. Certain values of the angle have conventional definitions:

- 
   0 = new moon
- 
   90 = first quarter
- 
   180 = full moon
- 
   270 = third quarter

#### Return

The angle as described above, a value in the range 0..360 degrees.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the observation. |
