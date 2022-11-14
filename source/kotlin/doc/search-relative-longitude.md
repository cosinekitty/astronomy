//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchRelativeLongitude](search-relative-longitude.md)

# searchRelativeLongitude

fun [searchRelativeLongitude](search-relative-longitude.md)(body: [Body](-body/index.md), targetRelativeLongitude: Double, startTime: [Time](-time/index.md)): [Time](-time/index.md)

Searches for the time when the Earth and another planet are separated by a specified angle in ecliptic longitude, as seen from the Sun.

A relative longitude is the angle between two bodies measured in the plane of the Earth's orbit (the ecliptic plane). The distance of the bodies above or below the ecliptic plane is ignored. If you imagine the shadow of the body cast onto the ecliptic plane, and the angle measured around that plane from one body to the other in the direction the planets orbit the Sun, you will get an angle somewhere between 0 and 360 degrees. This is the relative longitude.

Given a planet other than the Earth in body and a time to start the search in startTime, this function searches for the next time that the relative longitude measured from the planet to the Earth is targetRelLon.

Certain astronomical events are defined in terms of relative longitude between the Earth and another planet:

- 
   When the relative longitude is 0 degrees, it means both planets are in the same direction from the Sun. For planets that orbit closer to the Sun (Mercury and Venus), this is known as *inferior conjunction*, a time when the other planet becomes very difficult to see because of being lost in the Sun's glare. (The only exception is in the rare event of a transit, when we see the silhouette of the planet passing between the Earth and the Sun.)
- 
   When the relative longitude is 0 degrees and the other planet orbits farther from the Sun, this is known as *opposition*.  Opposition is when the planet is closest to the Earth, and also when it is visible for most of the night, so it is considered the best time to observe the planet.
- 
   When the relative longitude is 180 degrees, it means the other planet is on the opposite side of the Sun from the Earth. This is called *superior conjunction*. Like inferior conjunction, the planet is very difficult to see from the Earth. Superior conjunction is possible for any planet other than the Earth.

#### Return

The time of the first relative longitude event that occurs after startTime.

## Parameters

| | |
|---|---|
| body | A planet other than the Earth. Any other body will cause an exception. |
| targetRelativeLongitude | The desired relative longitude, expressed in degrees. Must be in the range [0, 360). |
| startTime | The date and time at which to begin the search. |
