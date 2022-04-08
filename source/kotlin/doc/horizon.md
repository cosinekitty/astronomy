//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[horizon](horizon.md)

# horizon

[jvm]\
fun [horizon](horizon.md)(time: [AstroTime](-astro-time/index.md), observer: [Observer](-observer/index.md), ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), refraction: [Refraction](-refraction/index.md)): [Topocentric](-topocentric/index.md)

Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.

Given a date and time, the geographic location of an observer on the Earth, and equatorial coordinates (right ascension and declination) of a celestial body, this function returns horizontal coordinates (azimuth and altitude angles) for the body relative to the horizon at the geographic location.

The right ascension ra and declination dec passed in must be *equator of date* coordinates, based on the Earth's true equator at the date and time of the observation. Otherwise the resulting horizontal coordinates will be inaccurate. Equator of date coordinates can be obtained by calling [equator](equator.md), passing in [EquatorEpoch.OfDate](-equator-epoch/-of-date/index.md) as its equdate parameter. It is also recommended to enable aberration correction by passing in [Aberration.Corrected](-aberration/-corrected/index.md) as the aberration parameter.

This function optionally corrects for atmospheric refraction. For most uses, it is recommended to pass [Refraction.Normal](-refraction/-normal/index.md) in the refraction parameter to correct for optical lensing of the Earth's atmosphere that causes objects to appear somewhat higher above the horizon than they actually are. However, callers may choose to avoid this correction by passing in [Refraction.None](-refraction/-none/index.md). If refraction correction is enabled, the azimuth, altitude, right ascension, and declination in the [Topocentric](-topocentric/index.md) object returned by this function will all be corrected for refraction. If refraction is disabled, none of these four coordinates will be corrected; in that case, the right ascension and declination in the returned structure will be numerically identical to the respective ra and dec values passed in.

#### Return

The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the observation. |
| observer | The geographic location of the observer. |
| ra | The right ascension of the body in sidereal hours. See remarks above for more details. |
| dec | The declination of the body in degrees. See remarks above for more details. |
| refraction | Selects whether to correct for atmospheric refraction, and if so, which model to use.     The recommended value for most uses is `Refraction.Normal`.     See remarks above for more details. |
