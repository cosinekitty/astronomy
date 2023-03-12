//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[atmosphere](atmosphere.md)

# atmosphere

fun [atmosphere](atmosphere.md)(elevationMeters: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [AtmosphereInfo](-atmosphere-info/index.md)

Calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.

This function calculates idealized values of pressure, temperature, and density using the U.S. Standard Atmosphere (1976) model.

1. 
   COESA, U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, DC, 1976.
2. 
   Jursa, A. S., Ed., Handbook of Geophysics and the Space Environment, Air Force Geophysics Laboratory, 1985. See: https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf

## Parameters

| | |
|---|---|
| elevationMeters | The elevation above sea level at which to calculate atmospheric variables. Must be in the range -500 to +100000, or an exception will occur. |
