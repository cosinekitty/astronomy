//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[seasons](seasons.md)

# seasons

fun [seasons](seasons.md)(year: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)): [SeasonsInfo](-seasons-info/index.md)

Finds both equinoxes and both solstices for a given calendar year.

The changes of seasons are defined by solstices and equinoxes. Given a calendar year number, this function calculates the March and September equinoxes and the June and December solstices.

The equinoxes are the moments twice each year when the plane of the Earth's equator passes through the center of the Sun. In other words, the Sun's declination is zero at both equinoxes. The March equinox defines the beginning of spring in the northern hemisphere and the beginning of autumn in the southern hemisphere. The September equinox defines the beginning of autumn in the northern hemisphere and the beginning of spring in the southern hemisphere.

The solstices are the moments twice each year when one of the Earth's poles is most tilted toward the Sun. More precisely, the Sun's declination reaches its minimum value at the December solstice, which defines the beginning of winter in the northern hemisphere and the beginning of summer in the southern hemisphere. The Sun's declination reaches its maximum value at the June solstice, which defines the beginning of summer in the northern hemisphere and the beginning of winter in the southern hemisphere.

#### Return

A [SeasonsInfo](-seasons-info/index.md) object that contains four [Time](-time/index.md) values: the March and September equinoxes and the June and December solstices.

## Parameters

| | |
|---|---|
| year | The calendar year number for which to calculate equinoxes and solstices. The value may be any integer, but only the years 1800 through 2100 have been validated for accuracy: unit testing against data from the United States Naval Observatory confirms that all equinoxes and solstices for that range of years are within 2 minutes of the correct time. |
