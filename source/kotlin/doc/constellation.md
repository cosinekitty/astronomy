//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[constellation](constellation.md)

# constellation

fun [constellation](constellation.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [ConstellationInfo](-constellation-info/index.md)

Determines the constellation that contains the given point in the sky.

Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the constellation that contains that point.

#### Return

A structure that contains the 3-letter abbreviation and full name of the constellation that contains the given (ra,dec), along with the converted B1875 (ra,dec) for that point.

## Parameters

jvm

| | |
|---|---|
| ra | The right ascension (RA) of a point in the sky, using the J2000 equatorial system (EQJ). |
| dec | The declination (DEC) of a point in the sky, using the J2000 equatorial system (EQJ). |
