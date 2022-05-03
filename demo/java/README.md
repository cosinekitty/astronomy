# Astronomy Engine examples in Java

The class file [Main.java](src/main/java/io/github/cosinekitty/astronomy/demo/Main.java)
shows examples of how to use the Kotlin version of Astronomy Engine from a Java program.

The demo program is split into separate source files, one for each topic, as listed below.
Each topic is chosen by a command line parameter.
To build the demo program, run this command in Linux or Mac:

```
./gradlew jar
```

On Windows, use this command to build the demo program:

```
gradlew.bat jar
```

Then use the `rundemo` script (or `rundemo.bat` file on Windows) to
run the demo program, to see usage text:

```
./rundemo
```

As an example, to run the MoonPhase demo, try this:

```
./rundemo moonphase
```

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

### [MoonPhase.java](src/main/java/io/github/cosinekitty/astronomy/demo/MoonPhase.java)
This example shows how to determine the Moon's current phase,
and how to predict when the next 10 quarter phases will occur.

### [Positions.java](src/main/java/io/github/cosinekitty/astronomy/demo/Positions.java)
Given an observer's geographic latitude and longitude,
and an optional date and time, this demo displays the
equatorial and horizontal coordinates of the Sun, Moon, and planets.

### [Seasons.java](src/main/java/io/github/cosinekitty/astronomy/demo/Seasons.java)
Calculates the equinoxes and solstices for a given calendar year.
