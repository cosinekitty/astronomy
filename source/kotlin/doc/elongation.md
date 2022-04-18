//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[elongation](elongation.md)

# elongation

[jvm]\
fun [elongation](elongation.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [ElongationInfo](-elongation-info/index.md)

Determines visibility of a celestial body relative to the Sun, as seen from the Earth.

This function returns an [ElongationInfo](-elongation-info/index.md) object, which provides the following information about the given celestial body at the given time:

- 
   visibility is an enumerated type that specifies whether the body is more easily seen     in the morning before sunrise, or in the evening after sunset.
- 
   elongation is the angle in degrees between two vectors: one from the center of the Earth to the     center of the Sun, the other from the center of the Earth to the center of the specified body.     This angle indicates how far away the body is from the glare of the Sun.     The elongation angle is always in the range 0, 180.
- 
   eclipticSeparation is the absolute value of the difference between the body's ecliptic longitude and the Sun's ecliptic longitude, both as seen from the center of the Earth. This angle measures around the plane of the Earth's orbit, and ignores how far above or below that plane the body is. The ecliptic separation is measured in degrees and is always in the range 0, 180.

## Parameters

jvm

| | |
|---|---|
| body | The celestial body whose visibility is to be calculated. |
| time | The date and time of the observation. |
