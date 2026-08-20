# Astronomy Engine examples in Java

The class file [Main.java](src/main/java/io/github/cosinekitty/astronomy/demo/Main.java)
shows examples of how to use the Kotlin Multiplatform version of Astronomy Engine from a Java program.

The demo program is split into separate source files, one for each topic, as listed below.
Each topic is chosen by a command line parameter.

## Building and running

```
./gradlew run --args="seasons 2026"
./gradlew run --args="moonphase"
./gradlew run --args="positions 29 -81"
```

Run without arguments to see usage text:

```
./gradlew run
```

This requires the library at `../../source/kotlin` to be available (it is referenced
via Gradle composite build).

---

### [Constellation.java](src/main/java/io/github/cosinekitty/astronomy/demo/Constellation.java)
This demo finds what constellation the Moon is in at a given time.
It also shows how to do a binary search to find the moment in time
when the Moon moves across the border between constellations.

### [JupiterMoons.java](src/main/java/io/github/cosinekitty/astronomy/demo/JupiterMoons.java)
Calculates the coordinates of Jupiter and its four major moons
(Io, Europa, Ganymede, and Callisto) as seen from the Earth
at a given date and time. This demo illustrates how to correct
for the delay caused by the time it takes for light to reach
the Earth from the Jupiter system.

### [LunarEclipse.java](src/main/java/io/github/cosinekitty/astronomy/demo/LunarEclipse.java)
Searches for the first 10 lunar eclipses (partial or total)
that occur after the specified time. Penumbral lunar eclipses
are ignored, as these are difficult to observe in practice.

### [MoonPhase.java](src/main/java/io/github/cosinekitty/astronomy/demo/MoonPhase.java)
This example shows how to determine the Moon's current phase,
and how to predict when the next 10 quarter phases will occur.

### [Positions.java](src/main/java/io/github/cosinekitty/astronomy/demo/Positions.java)
Given an observer's geographic latitude and longitude,
and an optional date and time, this demo displays the
equatorial and horizontal coordinates of the Sun, Moon, and planets.

### [RiseSetCulm.java](src/main/java/io/github/cosinekitty/astronomy/demo/RiseSetCulm.java)
Calculates rise, set, and culmination times of the Sun and Moon.
Culmination is when a body reaches the highest
point in an observer's sky as it crosses the meridian.

### [Seasons.java](src/main/java/io/github/cosinekitty/astronomy/demo/Seasons.java)
Calculates the equinoxes and solstices for a given calendar year.
