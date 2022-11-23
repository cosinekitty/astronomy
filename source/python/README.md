# Astronomy Engine (Python)

This is the complete programming reference for the Python version of
Astronomy Engine. Supports Python 3. Does NOT support Python 2.
See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start

[![pypi](https://img.shields.io/pypi/v/astronomy-engine)](https://pypi.org/project/astronomy-engine/)

To include Astronomy Engine in your own Python program,
you can use the [astronomy-engine](https://pypi.org/project/astronomy-engine/) package:

```
pip install astronomy-engine
```

Alternatively, you can copy the file [astronomy/astronomy.py](astronomy/astronomy.py)
into your project directory.

With either approach, add the following line toward the top of your program:

```python
import astronomy
```


To get started quickly, here are some [examples](../../demo/python/).

---

## Contents

- [Topic Index](#topics)
- [Constants](#constants)
- [Classes](#classes)
- [Enumerated Types](#enumerations)
- [Error Types](#errors)
- [Functions](#functions)

---

<a name="topics"></a>
## Topic Index

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#HelioVector) | Calculates body position vector with respect to the center of the Sun. |
| [GeoVector](#GeoVector)     | Calculates body position vector with respect to the center of the Earth. |
| [Equator](#Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Ecliptic)       | Converts J2000 equatorial coordinates to J2000 ecliptic coordinates. |
| [EclipticLongitude](#EclipticLongitude) | Calculates ecliptic longitude of a body in the J2000 system. |
| [Horizon](#Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |
| [PairLongitude](#PairLongitude) | Calculates the difference in apparent ecliptic longitude between two bodies, as seen from the Earth. |
| [BaryState](#BaryState) | Calculates the barycentric position and velocity vectors of the Sun or a planet. |

### Geographic helper functions

| Function | Description |
| -------- | ----------- |
| [ObserverVector](#ObserverVector) | Calculates a vector from the center of the Earth to an observer on the Earth's surface. |
| [VectorObserver](#VectorObserver) | Calculates the geographic coordinates for a geocentric equatorial vector. |

### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSet](#SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchAltitude](#SearchAltitude) | Finds time when a body reaches a given altitude above or below the horizon. Useful for finding civil, nautical, or astronomical twilight. |
| [SearchHourAngle](#SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonPhase](#SearchMoonPhase) | Finds the next instance of the Moon reaching a specific ecliptic longitude separation from the Sun. |
| [SearchMoonQuarter](#SearchMoonQuarter) | Finds the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#NextMoonQuarter) | Finds the next quarter moon phase after a previous one that has been found. |

### Eclipses and Transits

| Function | Description |
| -------- | ----------- |
| [SearchLunarEclipse](#SearchLunarEclipse) | Search for the first lunar eclipse after a given date. |
| [NextLunarEclipse](#NextLunarEclipse) | Continue searching for more lunar eclipses. |
| [SearchGlobalSolarEclipse](#SearchGlobalSolarEclipse) | Search for the first solar eclipse after a given date that is visible anywhere on the Earth. |
| [NextGlobalSolarEclipse](#NextGlobalSolarEclipse) | Continue searching for solar eclipses visible anywhere on the Earth. |
| [SearchLocalSolarEclipse](#SearchLocalSolarEclipse) | Search for the first solar eclipse after a given date that is visible at a particular location on the Earth. |
| [NextLocalSolarEclipse](#NextLocalSolarEclipse) | Continue searching for solar eclipses visible at a particular location on the Earth. |
| [SearchTransit](#SearchTransit) | Search for the next transit of Mercury or Venus. |
| [NextTransit](#NextTransit) | Continue searching for transits of Mercury or Venus. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#NextLunarApsis) | Given an already-found apsis, finds the next perigee or apogee of the Moon. |

### Planet perihelion and aphelion

| Function | Description |
| -------- | ----------- |
| [SearchPlanetApsis](#SearchPlanetApsis) | Finds the next perihelion or aphelion of a planet after a specified date. |
| [NextPlanetApsis](#NextPlanetApsis) | Given an already-found apsis, finds the next perihelion or aphelion of a planet. |

### Visual magnitude and elongation

| Function | Description |
| -------- | ----------- |
| [Illumination](#Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |
| [SearchPeakMagnitude](#SearchPeakMagnitude) | Searches for the date and time Venus will next appear brightest as seen from the Earth. |
| [AngleFromSun](#AngleFromSun) | Returns full angle seen from Earth between body and Sun. |
| [Elongation](#Elongation) | Calculates ecliptic longitude angle between a body and the Sun, as seen from the Earth. |
| [SearchMaxElongation](#SearchMaxElongation) | Searches for the next maximum elongation event for Mercury or Venus that occurs after the given date. |

### Oppositions and conjunctions

| Function | Description |
| -------- | ----------- |
| [SearchRelativeLongitude](#SearchRelativeLongitude) | Finds oppositions and conjunctions of planets. |

### Equinoxes, solstices, and apparent solar motion

| Function | Description |
| -------- | ----------- |
| [SearchSunLongitude](#SearchSunLongitude) | Finds the next time the Sun reaches a specified apparent ecliptic longitude in the *true equator of date* system. |
| [Seasons](#Seasons) | Finds the equinoxes and solstices for a given calendar year. |
| [SunPosition](#SunPosition) | Calculates the Sun's apparent ecliptic coordinates as seen from the Earth. |

### Coordinate transforms

The following five orientation systems are supported.
Astronomy Engine can convert a vector from any of these orientations to any of the others.
It also allows converting from a vector to spherical (angular) coordinates and back,
within a given orientation. Note the 3-letter codes for each of the orientation systems;
these are used in function and type names.

- **EQJ = Equatorial J2000**: Uses the Earth's equator on January 1, 2000, at noon UTC.
- **EQD = Equator of-date**: Uses the Earth's equator on a given date and time, adjusted for precession and nutation.
- **ECL = Ecliptic**: Uses the mean plane of the Earth's orbit around the Sun. The x-axis is referenced against the J2000 equinox.
- **HOR = Horizontal**: Uses the viewpoint of an observer at a specific location on the Earth at a given date and time.
- **GAL = Galactic**: Based on the IAU 1958 definition of galactic coordinates.

| Function | Description |
| -------- | ----------- |
| [RotateVector](#RotateVector) | Applies a rotation matrix to a vector, yielding a vector in another orientation system. |
| [InverseRotation](#InverseRotation) | Given a rotation matrix, finds the inverse rotation matrix that does the opposite transformation. |
| [CombineRotation](#CombineRotation) | Given two rotation matrices, returns a rotation matrix that combines them into a net transformation. |
| [IdentityMatrix](#IdentityMatrix) | Returns a 3x3 identity matrix, which can be used to form other rotation matrices. |
| [Pivot](#Pivot) | Transforms a rotation matrix by pivoting it around a given axis by a given angle. |
| [VectorFromSphere](#VectorFromSphere) | Converts spherical coordinates to Cartesian coordinates. |
| [SphereFromVector](#SphereFromVector) | Converts Cartesian coordinates to spherical coordinates. |
| [EquatorFromVector](#EquatorFromVector) | Given an equatorial vector, calculates equatorial angular coordinates. |
| [VectorFromHorizon](#VectorFromHorizon) | Given apparent angular horizontal coordinates, calculates horizontal vector. |
| [HorizonFromVector](#HorizonFromVector) | Given a vector in horizontal orientation, calculates horizontal angular coordinates. |
| [Rotation_EQD_EQJ](#Rotation_EQD_EQJ) | Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ). |
| [Rotation_EQD_ECL](#Rotation_EQD_ECL) | Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL). |
| [Rotation_EQD_HOR](#Rotation_EQD_HOR) | Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR). |
| [Rotation_EQJ_EQD](#Rotation_EQJ_EQD) | Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD). |
| [Rotation_EQJ_ECL](#Rotation_EQJ_ECL) | Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL). |
| [Rotation_EQJ_HOR](#Rotation_EQJ_HOR) | Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR). |
| [Rotation_ECL_EQD](#Rotation_ECL_EQD) | Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD). |
| [Rotation_ECL_EQJ](#Rotation_ECL_EQJ) | Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ). |
| [Rotation_ECL_HOR](#Rotation_ECL_HOR) | Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR). |
| [Rotation_HOR_EQD](#Rotation_HOR_EQD) | Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD). |
| [Rotation_HOR_EQJ](#Rotation_HOR_EQJ) | Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ). |
| [Rotation_HOR_ECL](#Rotation_HOR_ECL) | Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL). |
| [Rotation_EQJ_GAL](#Rotation_EQJ_GAL) | Calculates a rotation matrix from equatorial J2000 (EQJ) to galactic (GAL). |
| [Rotation_GAL_EQJ](#Rotation_GAL_EQJ) | Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ). |

### Gravitational simulation of small bodies

Astronomy Engine provides a [GravitySimulator](#GravitySimulator) class
that allows you to model the trajectories of one or more small bodies like asteroids,
comets, or coasting spacecraft. If you know an initial position vector
and velocity vector for a small body, the gravity simulator can incrementally
simulate the pull of gravity on it from the Sun and planets, to calculate its
movement through the Solar System.

---

---

<a name="constants"></a>
## Constants
The following numeric constants are exported by the `astronomy` module.
They may be of use for unit conversion.
Note: For the other supported programming languages, Astronomy Engine defines
helper constants `DEG2RAD` and `RAD2DEG` to convert between angular degrees and radians.
However, because Python defines the [angular conversion functions](https://docs.python.org/3/library/math.html#angular-conversion)
`math.degrees()` and `math.radians()`, they are not needed in the Python version.

---

<a name="AU_PER_LY"></a>
### `AU_PER_LY = 63241.07708807546`

**The number of astronomical units in one light-year.**

---

<a name="CALLISTO_RADIUS_KM"></a>
### `CALLISTO_RADIUS_KM = 2410.3`

**The mean radius of Jupiter's moon Callisto, expressed in kilometers.**

---

<a name="C_AUDAY"></a>
### `C_AUDAY = 173.1446326846693`

**The speed of light expressed in astronomical units per day.**

---

<a name="EUROPA_RADIUS_KM"></a>
### `EUROPA_RADIUS_KM = 1560.8`

**The mean radius of Jupiter's moon Europa, expressed in kilometers.**

---

<a name="GANYMEDE_RADIUS_KM"></a>
### `GANYMEDE_RADIUS_KM = 2631.2`

**The mean radius of Jupiter's moon Ganymede, expressed in kilometers.**

---

<a name="IO_RADIUS_KM"></a>
### `IO_RADIUS_KM = 1821.6`

**The mean radius of Jupiter's moon Io, expressed in kilometers.**

---

<a name="JUPITER_EQUATORIAL_RADIUS_KM"></a>
### `JUPITER_EQUATORIAL_RADIUS_KM = 71492.0`

**The equatorial radius of Jupiter, expressed in kilometers.**

---

<a name="JUPITER_MEAN_RADIUS_KM"></a>
### `JUPITER_MEAN_RADIUS_KM = 69911.0`

**The volumetric mean radius of Jupiter, expressed in kilometers.**

---

<a name="JUPITER_POLAR_RADIUS_KM"></a>
### `JUPITER_POLAR_RADIUS_KM = 66854.0`

**The polar radius of Jupiter, expressed in kilometers.**

---

<a name="KM_PER_AU"></a>
### `KM_PER_AU = 1.4959787069098932e+8`

**The number of kilometers per astronomical unit.**

---

<a name="classes"></a>
## Classes

---

<a name="Apsis"></a>
### class Apsis

**An event where a satellite is closest to or farthest from the body it orbits.**

For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
event where the orbiting body reaches its closest or farthest point from the primary body.
The closest approach is called *pericenter* and the farthest point is *apocenter*.
More specific terminology is common for particular orbiting bodies.
The Moon's closest approach to the Earth is called *perigee* and its furthest
point is called *apogee*. The closest approach of a planet to the Sun is called
*perihelion* and the furthest point is called *aphelion*.
This data structure is returned by [`SearchLunarApsis`](#SearchLunarApsis) and [`NextLunarApsis`](#NextLunarApsis)
to iterate through consecutive alternating perigees and apogees.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the apsis. |
| [`ApsisKind`](#ApsisKind) | `kind` | Whether this is a pericenter or apocenter event. |
| `float` | `dist_au` | The distance between the centers of the bodies in astronomical units. |
| `float` | `dist_km` | The distance between the centers of the bodies in kilometers. |

---

<a name="AxisInfo"></a>
### class AxisInfo

**Information about a body's rotation axis at a given time.**

This structure is returned by [`RotationAxis`](#RotationAxis) to report
the orientation of a body's rotation axis at a given moment in time.
The axis is specified by the direction in space that the body's north pole
points, using angular equatorial coordinates in the J2000 system (EQJ).
Thus `ra` is the right ascension, and `dec` is the declination, of the
body's north pole vector at the given moment in time. The north pole
of a body is defined as the pole that lies on the north side of the
[Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
regardless of the body's direction of rotation.
The `spin` field indicates the angular position of a prime meridian
arbitrarily recommended for the body by the International Astronomical
Union (IAU).
The fields `ra`, `dec`, and `spin` correspond to the variables
α0, δ0, and W, respectively, from
[Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
The field `north` is a unit vector pointing in the direction of the body's north pole.
It is expressed in the equatorial J2000 system (EQJ).

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ra` | The J2000 right ascension of the body's north pole direction, in sidereal hours. |
| `float` | `dec` | The J2000 declination of the body's north pole direction, in degrees. |
| `float` | `spin` | Rotation angle of the body's prime meridian, in degrees. |
| [`Vector`](#Vector) | `north` | A J2000 dimensionless unit vector pointing in the direction of the body's north pole. |

---

<a name="ConstellationInfo"></a>
### class ConstellationInfo

**Reports the constellation that a given celestial point lies within.**

The [`Constellation`](#Constellation) function returns a `ConstellationInfo` object
to report which constellation corresponds with a given point in the sky.
Constellations are defined with respect to the B1875 equatorial system
per IAU standard. Although the `Constellation` function requires J2000 equatorial
coordinates as input, the returned object contains converted B1875 coordinates for reference.

| Type | Attribute | Description |
| --- | --- | --- |
| `string` | `symbol` | 3-character mnemonic symbol for the constellation, e.g. "Ori". |
| `string` | `name` | Full name of constellation, e.g. "Orion". |
| `float` | `ra1875` | Right ascension expressed in B1875 coordinates. |
| `float` | `dec1875` | Declination expressed in B1875 coordinates. |

---

<a name="EclipseEvent"></a>
### class EclipseEvent

**Holds a time and the observed altitude of the Sun at that time.**

When reporting a solar eclipse observed at a specific location on the Earth
(a "local" solar eclipse), a series of events occur. In addition
to the time of each event, it is important to know the altitude of the Sun,
because each event may be invisible to the observer if the Sun is below
the horizon.
If `altitude` is negative, the event is theoretical only; it would be
visible if the Earth were transparent, but the observer cannot actually see it.
If `altitude` is positive but less than a few degrees, visibility will be impaired by
atmospheric interference (sunrise or sunset conditions).

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the event. |
| `float` | `altitude` | The angular altitude of the center of the Sun above/below the horizon, at `time`, corrected for atmospheric refraction and expressed in degrees. |

---

<a name="EclipticCoordinates"></a>
### class EclipticCoordinates

**Ecliptic angular and Cartesian coordinates.**

Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).

| Type | Attribute | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vec` | Ecliptic cartesian vector with the following components: x: in the direction of the equinox along the ecliptic plane. y: Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox. z: Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north. |
| `float` | `elat` | Latitude in degrees north (positive) or south (negative) of the ecliptic plane. |
| `float` | `elon` | Longitude in degrees around the ecliptic plane prograde from the equinox. |

---

<a name="ElongationEvent"></a>
### class ElongationEvent

**Contains information about the visibility of a celestial body at a given date and time.**

See the [`Elongation`](#Elongation) function for more detailed information about the members of this class.
See also [`SearchMaxElongation`](#SearchMaxElongation) for how to search for maximum elongation events.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |
| [`Visibility`](#Visibility) | `visibility` | Whether the body is best seen in the morning or the evening. |
| `float` | `elongation` | The angle in degrees between the body and the Sun, as seen from the Earth. |
| `float` | `ecliptic_separation` | The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth. |

---

<a name="Equatorial"></a>
### class Equatorial

**Equatorial angular coordinates**

Coordinates of a celestial body as seen from the Earth.
Can be geocentric or topocentric, depending on context.
The coordinates are oriented with respect to the Earth's
equator projected onto the sky.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ra` | Right ascension in sidereal hours. |
| `float` | `dec` | Declination in degrees. |
| `float` | `dist` | Distance to the celestial body in AU. |
| [`Vector`](#Vector) | `vec` | The equatorial coordinates in cartesian form, using AU distance units. x = direction of the March equinox, y = direction of the June solstice, z = north. |

---

<a name="GlobalSolarEclipseInfo"></a>
### class GlobalSolarEclipseInfo

**Reports the time and geographic location of the peak of a solar eclipse.**

Returned by [`SearchGlobalSolarEclipse`](#SearchGlobalSolarEclipse) or [`NextGlobalSolarEclipse`](#NextGlobalSolarEclipse)
to report information about a solar eclipse event.
The eclipse is classified as partial, annular, or total, depending on the
maximum amount of the Sun's disc obscured, as seen at the peak location
on the surface of the Earth.
The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
observes either a total or annular eclipse.
If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
fields give the geographic coordinates of the center of the Moon's shadow projected
onto the daytime side of the Earth at the instant of the eclipse's peak.
If `kind` has any other value, `latitude` and `longitude` are undefined and should
not be used.
For total or annular eclipses, the `obscuration` field holds the fraction (0, 1]
of the Sun's apparent disc area that is blocked from view by the Moon's silhouette,
as seen by an observer located at the geographic coordinates `latitude`, `longitude`
at the darkest time `peak`. The value will always be 1 for total eclipses, and less than
1 for annular eclipses.
For partial eclipses, `obscuration` holds the value `None`.
This is because there is little practical use for an obscuration value of
a partial eclipse without supplying a particular observation location.
Developers who wish to find an obscuration value for partial solar eclipses should therefore use
[`SearchLocalSolarEclipse`](#SearchLocalSolarEclipse) and provide the geographic coordinates of an observer.

| Type | Attribute | Description |
| --- | --- | --- |
| [`EclipseKind`](#EclipseKind) | `kind` | The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`. |
| `float` | `obscuration` | The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only). |
| [`Time`](#Time) | `peak` | The date and time when the solar eclipse is darkest. This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center. |
| `float` | `distance` | The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers. |
| `float` | `latitude` | The geographic latitude at the center of the peak eclipse shadow. |
| `float` | `longitude` | The geographic longitude at the center of the peak eclipse shadow. |

---

<a name="GravitySimulator"></a>
### class GravitySimulator

**A simulation of zero or more small bodies moving through the Solar System.**

This class calculates the movement of arbitrary small bodies,
such as asteroids or comets, that move through the Solar System.
It does so by calculating the gravitational forces on the bodies
from the Sun and planets. The user of this class supplies a
list of initial positions and velocities for the small bodies.
Then the class can update the positions and velocities over small
time steps.

#### member functions

<a name="GravitySimulator.__init__"></a>
### GravitySimulator.__init__(self, originBody, time, bodyStates)

**Creates a gravity simulation object.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `originBody` | Specifies the origin of the reference frame. All position vectors and velocity vectors will use `originBody` as the origin of the coordinate system. This origin applies to all the input vectors provided in the `bodyStates` parameter of this function, along with all output vectors returned by [`GravitySimulator.Update`](#GravitySimulator.Update). Most callers will want to provide one of the following: `Body.Sun` for heliocentric coordinates, `Body.SSB` for solar system barycentric coordinates, or `Body.Earth` for geocentric coordinates. Note that the gravity simulator does not correct for light travel time; all state vectors are tied to a Newtonian "instantaneous" time. |
| [`Time`](#Time) | `time` | The initial time at which to start the simulation. |
| [`StateVector`](#StateVector)`[]` | `bodyStates` | An array of zero or more initial state vectors (positions and velocities) of the small bodies to be simulated. The caller must know the positions and velocities of the small bodies at an initial moment in time. Their positions and velocities are expressed with respect to `originBody`, using equatorial J2000 orientation (EQJ). Positions are expressed in astronomical units (AU). Velocities are expressed in AU/day. All the times embedded within the state vectors must exactly match `time`, or this constructor will throw an exception. |

<a name="GravitySimulator.OriginBody"></a>
### GravitySimulator.OriginBody(self)

**The origin of the reference frame. See constructor for more info.**

**Returns**: [`Body`](#Body)

<a name="GravitySimulator.SolarSystemBodyState"></a>
### GravitySimulator.SolarSystemBodyState(self, body)

**Get the position and velocity of a Solar System body included in the simulation.**

In order to simulate the movement of small bodies through the Solar System,
the simulator needs to calculate the state vectors for the Sun and planets.
If an application wants to know the positions of one or more of the planets
in addition to the small bodies, this function provides a way to obtain
their state vectors. This is provided for the sake of efficiency, to avoid
redundant calculations.
The state vector is returned relative to the position and velocity
of the `originBody` parameter that was passed to this object's constructor.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune. |

**Returns**: [`StateVector`](#StateVector)
The state vector of the requested Solar System body.

<a name="GravitySimulator.Swap"></a>
### GravitySimulator.Swap(self)

**Exchange the current time step with the previous time step.**

Sometimes it is helpful to "explore" various times near a given
simulation time step, while repeatedly returning to the original
time step. For example, when backdating a position for light travel
time, the caller may wish to repeatedly try different amounts of
backdating. When the backdating solver has converged, the caller
wants to leave the simulation in its original state.
This function allows a single "undo" of a simulation, and does so
very efficiently.
Usually this function will be called immediately after a matching
call to [`GravitySimulator.Update`](#GravitySimulator.Update). It has the effect of rolling
back the most recent update. If called twice in a row, it reverts
the swap and thus has no net effect.
The constructor initializes the current state and previous
state to be identical. Both states represent the `time` parameter that was
passed into the constructor. Therefore, `Swap` will
have no effect from the caller's point of view when passed a simulator
that has not yet been updated by a call to [`GravitySimulator.Update`](#GravitySimulator.Update).

<a name="GravitySimulator.Time"></a>
### GravitySimulator.Time(self)

**The time represented by the current step of the gravity simulation.**

**Returns**: [`Time`](#Time)

<a name="GravitySimulator.Update"></a>
### GravitySimulator.Update(self, time)

**Advances the gravity simulation by a small time step.**

Updates the simulation of the user-supplied small bodies
to the time indicated by the `time` parameter.
Returns an array of state vectors for the simulated bodies.
The array is in the same order as the original array that
was used to construct this simulator object.
The positions and velocities in the returned array are
referenced to the `originBody` that was used to construct
this simulator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | A time that is a small increment away from the current simulation time. It is up to the developer to figure out an appropriate time increment. Depending on the trajectories, a smaller or larger increment may be needed for the desired accuracy. Some experimentation may be needed. Generally, bodies that stay in the outer Solar System and move slowly can use larger time steps. Bodies that pass into the inner Solar System and move faster will need a smaller time step to maintain accuracy. The `time` value may be after or before the current simulation time to move forward or backward in time. |

**Returns**: [`StateVector`](#StateVector)`[]`
An array of state vectors, one for each small body.

---

<a name="HorizontalCoordinates"></a>
### class HorizontalCoordinates

**Coordinates of a celestial body as seen by a topocentric observer.**

Contains horizontal and equatorial coordinates as seen by an observer
on or near the surface of the Earth (a topocentric observer).
All coordinates are optionally corrected for atmospheric refraction.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `azimuth` | The compass direction laterally around the observer's horizon, measured in degrees. North is 0 degrees, east is 90 degrees, south is 180 degrees, etc. |
| `float` | `altitude` | The angle in degrees above (positive) or below (negative) the observer's horizon. |
| `float` | `ra` | The right ascension in sidereal hours. |
| `float` | `dec` | The declination in degrees. |

---

<a name="HourAngleEvent"></a>
### class HourAngleEvent

**Information about a celestial body crossing a specific hour angle.**

Returned by the function [`SearchHourAngle`](#SearchHourAngle) to report information about
a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time when the body crosses the specified hour angle. |
| [`HorizontalCoordinates`](#HorizontalCoordinates) | `hor` | Apparent coordinates of the body at the time it crosses the specified hour angle. |

---

<a name="IlluminationInfo"></a>
### class IlluminationInfo

**Information about the brightness and illuminated shape of a celestial body.**

Returned by functions [`Illumination`](#Illumination) and [`SearchPeakMagnitude`](#SearchPeakMagnitude)
to report the visual magnitude and illuminated fraction of a celestial
body at a given date and time.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |
| `float` | `mag` | The visual magnitude of the body. Smaller values are brighter. |
| `float` | `phase_angle` | The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth. |
| `float` | `phase_fraction` | A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth. |
| `float` | `helio_dist` | The distance between the Sun and the body at the observation time, in AU. |
| `dist` | `geo_dist` | The distance between the Earth and the both at the observation time, in AU. |
| [`Vector`](#Vector) | `hc` | The body's heliocentric vector. |
| [`Vector`](#Vector) | `gc` | The body's geocentric vector. |
| `float` | `ring_tilt` | For Saturn, the tilt angle in degrees of its rings as seen from Earth. When the `ring_tilt` is very close to 0, it means the rings are edge-on as seen from observers on the Earth, and are thus very difficult to see. For bodies other than Saturn, `ring_tilt` is `None`. |

---

<a name="JupiterMoonsInfo"></a>
### class JupiterMoonsInfo

**Holds the positions and velocities of Jupiter's major 4 moons.**

The [`JupiterMoons`](#JupiterMoons) function returns an object of this type
to report position and velocity vectors for Jupiter's largest 4 moons
Io, Europa, Ganymede, and Callisto. Each position vector is relative
to the center of Jupiter. Both position and velocity are oriented in
the EQJ system (that is, using Earth's equator at the J2000 epoch).
The positions are expressed in astronomical units (AU),
and the velocities in AU/day.

| Type | Attribute | Description |
| --- | --- | --- |
| [`StateVector`](#StateVector) | `io` | The position and velocity of Jupiter's moon Io. |
| [`StateVector`](#StateVector) | `europa` | The position and velocity of Jupiter's moon Europa. |
| [`StateVector`](#StateVector) | `ganymede` | The position and velocity of Jupiter's moon Ganymede. |
| [`StateVector`](#StateVector) | `callisto` | The position and velocity of Jupiter's moon Callisto. |

---

<a name="LibrationInfo"></a>
### class LibrationInfo

**Lunar libration angles, returned by [`Libration`](#Libration).**

Contains lunar libration angles and lunar position information
for a given moment in time. See [`Libration`](#Libration) for more details.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `elat` | Sub-Earth libration ecliptic latitude angle, in degrees. |
| `float` | `elon` | Sub-Earth libration ecliptic longitude angle, in degrees. |
| `float` | `mlat` | Moon's geocentric ecliptic latitude, in degrees. |
| `float` | `mlon` | Moon's geocentric ecliptic longitude, in degrees. |
| `float` | `dist_km` | Distance between the centers of the Earth and Moon in kilometers. |
| `float` | `diam_deg` | The apparent angular diameter of the Moon as seen from the center of the Earth. |

---

<a name="LocalSolarEclipseInfo"></a>
### class LocalSolarEclipseInfo

**Information about a solar eclipse as seen by an observer at a given time and geographic location.**

Returned by [`SearchLocalSolarEclipse`](#SearchLocalSolarEclipse) or [`NextLocalSolarEclipse`](#NextLocalSolarEclipse)
to report information about a solar eclipse as seen at a given geographic location.
When a solar eclipse is found, it is classified as partial, annular, or total.
The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
A partial solar eclipse is when the Moon does not line up directly enough with the Sun
to completely block the Sun's light from reaching the observer.
An annular eclipse occurs when the Moon's disc is completely visible against the Sun
but the Moon is too far away to completely block the Sun's light; this leaves the
Sun with a ring-like appearance.
A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
Sun just right to completely block all sunlight from reaching the observer.
The `obscuration` field reports what fraction of the Sun's disc appears blocked
by the Moon when viewed by the observer at the peak eclipse time.
This is a value that ranges from 0 (no blockage) to 1 (total eclipse).
The obscuration value will be between 0 and 1 for partial eclipses and annular eclipses.
The value will be exactly 1 for total eclipses. Obscuration gives an indication
of how dark the eclipse appears.
There are 5 "event" fields, each of which contains a time and a solar altitude.
Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
The fields `partial_begin` and `partial_end` are always set, and indicate when
the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
`total_begin` and `total_end` indicate when the total/annular phase begins/ends.
When an event field is valid, the caller must also check its `altitude` field to
see whether the Sun is above the horizon at the time indicated by the `time` field.
See [`EclipseEvent`](#EclipseEvent) for more information.

| Type | Attribute | Description |
| --- | --- | --- |
| [`EclipseKind`](#EclipseKind) | `kind` | The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`. |
| `float` | `obscuration` | The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak. |
| [`EclipseEvent`](#EclipseEvent) | `partial_begin` | The time and Sun altitude at the beginning of the eclipse. |
| [`EclipseEvent`](#EclipseEvent) | `total_begin` | If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise `None`. |
| [`EclipseEvent`](#EclipseEvent) | `peak` | The time and Sun altitude when the eclipse reaches its peak. |
| [`EclipseEvent`](#EclipseEvent) | `total_end` | If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise `None`. |
| [`EclipseEvent`](#EclipseEvent) | `partial_end` | The time and Sun altitude at the end of the eclipse. |

---

<a name="LunarEclipseInfo"></a>
### class LunarEclipseInfo

**Returns information about a lunar eclipse.**

Returned by [`SearchLunarEclipse`](#SearchLunarEclipse) or [`NextLunarEclipse`](#NextLunarEclipse)
to report information about a lunar eclipse event.
When a lunar eclipse is found, it is classified as penumbral, partial, or total.
Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
Total eclipses occur when the entire Moon passes into the Earth's umbra.
The `kind` field thus holds one of the values `EclipseKind.Penumbral`, `EclipseKind.Partial`,
or `EclipseKind.Total`, depending on the kind of lunar eclipse found.
The `obscuration` field holds a value in the range [0, 1] that indicates what fraction
of the Moon's apparent disc area is covered by the Earth's umbra at the eclipse's peak.
This indicates how dark the peak eclipse appears. For penumbral eclipses, the obscuration
is 0, because the Moon does not pass through the Earth's umbra. For partial eclipses,
the obscuration is somewhere between 0 and 1. For total lunar eclipses, the obscuration is 1.
Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.
Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
of the eclipse, which is half of the amount of time the eclipse spends in each
phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
By converting from minutes to days, and subtracting/adding with `peak`, the caller
may determine the date and time of the beginning/end of each eclipse phase.

| Type | Attribute | Description |
| --- | --- | --- |
| [`EclipseKind`](#EclipseKind) | `kind` | The type of lunar eclipse found. |
| `float` | `obscuration` | The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra. |
| [`Time`](#Time) | `peak` | The time of the eclipse at its peak. |
| `float` | `sd_penum` | The semi-duration of the penumbral phase in minutes. |
| `float` | `sd_partial` | The semi-duration of the penumbral phase in minutes, or 0.0 if none. |
| `float` | `sd_total` | The semi-duration of the penumbral phase in minutes, or 0.0 if none. |

---

<a name="MoonQuarter"></a>
### class MoonQuarter

**A lunar quarter event along with its date and time.**

An object of this type represents one of the four major
lunar phases that appear on calendars:
new moon, first quarter, full moon, or third quarter.
Along with the `quarter` attribute that specifies the
type of quarter, it contains a `time` field that indicates
when the lunar quarter event happens.

| Type | Attribute | Description |
| --- | --- | --- |
| `int` | `quarter` | 0=new moon, 1=first quarter, 2=full moon, 3=third quarter. |
| [`Time`](#Time) | `time` | The date and time of the lunar quarter. |

---

<a name="NodeEventInfo"></a>
### class NodeEventInfo

**Information about an ascending or descending node of a body.**

This object is returned by [`SearchMoonNode`](#SearchMoonNode) and [`NextMoonNode`](#NextMoonNode)
to report information about the center of the Moon passing through the ecliptic plane.

| Type | Attribute | Description |
| --- | --- | --- |
| [`NodeEventKind`](#NodeEventKind) | `kind` | Whether the node is ascending (south to north) or descending (north to south). |
| [`Time`](#Time) | `time` | The time when the body passes through the ecliptic plane. |

---

<a name="Observer"></a>
### class Observer

**Represents the geographic location of an observer on the surface of the Earth.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `latitude` | Geographic latitude in degrees north of the equator. |
| `float` | `longitude` | Geographic longitude in degrees east of the prime meridian at Greenwich, England. |
| `float` | `height` | Elevation above sea level in meters. |

---

<a name="PositionFunction"></a>
### class PositionFunction

**A function for which to solve a light-travel time problem.**

This abstract class defines the contract for wrapping a
position vector as a function of time. A class derived from
`PositionFunction` must define a `Position` method that
returns a position vector for a given time.
The function [`CorrectLightTravel`](#CorrectLightTravel) solves a generalized
problem of deducing how far in the past light must have
left a target object to be seen by an observer at a
specified time. It is passed an instance of `PositionFunction`
that expresses a relative position vector function.

#### member functions

<a name="PositionFunction.Position"></a>
### PositionFunction.Position(self, time)

**Returns a relative position vector for a given time.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The time at which to evaluate a relative position vector. |

**Returns**: [`Vector`](#Vector)

---

<a name="RotationMatrix"></a>
### class RotationMatrix

Contains a rotation matrix that can be used to transform one
coordinate system into another.

| Type | Parameter | Description |
| --- | --- | --- |
| `float[3][3]` | `rot` | A normalized 3x3 rotation matrix. |

---

<a name="SeasonInfo"></a>
### class SeasonInfo

**The dates and times of changes of season for a given calendar year.**

Call [`Seasons`](#Seasons) to calculate this data structure for a given year.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `mar_equinox` | The date and time of the March equinox for the specified year. |
| [`Time`](#Time) | `jun_solstice` | The date and time of the June solstice for the specified year. |
| [`Time`](#Time) | `sep_equinox` | The date and time of the September equinox for the specified year. |
| [`Time`](#Time) | `dec_solstice` | The date and time of the December solstice for the specified year. |

---

<a name="Spherical"></a>
### class Spherical

**Holds spherical coordinates: latitude, longitude, distance.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `lat` | The latitude angle: -90..+90 degrees. |
| `float` | `lon` | The longitude angle: 0..360 degrees. |
| `float` | `dist` | Distance in AU. |

---

<a name="StateVector"></a>
### class StateVector

**A combination of a position vector, a velocity vector, and a time.**

The position (x, y, z) is measured in astronomical units (AU).
The velocity (vx, vy, vz) is measured in AU/day.
The coordinate system varies and depends on context.
The state vector also includes a time stamp.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `x` | The x-coordinate of the position, measured in AU. |
| `float` | `y` | The y-coordinate of the position, measured in AU. |
| `float` | `z` | The z-coordinate of the position, measured in AU. |
| `float` | `vx` | The x-component of the velocity, measured in AU/day. |
| `float` | `vy` | The y-component of the velocity, measured in AU/day. |
| `float` | `vz` | The z-component of the velocity, measured in AU/day. |
| [`Time`](#Time) | `t` | The date and time at which the position and velocity vectors are valid. |

#### member functions

<a name="StateVector.Position"></a>
### StateVector.Position(self)

Extracts a position vector from this state vector.

<a name="StateVector.Velocity"></a>
### StateVector.Velocity(self)

Extracts a velocity vector from this state vector.

---

<a name="Time"></a>
### class Time

**Represents a date and time used for performing astronomy calculations.**

All calculations performed by Astronomy Engine are based on
dates and times represented by `Time` objects.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ut` | UT1/UTC number of days since noon on January 1, 2000. See the `ut` attribute of this class for more details. |

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ut` | The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine. Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed. UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed. The value in `ut` is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time. Before the era of atomic timekeeping, days based on the Earth's rotation were often known as *mean solar days*. |
| `float` | `tt` | Terrestrial Time days since noon on January 1, 2000. Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html) divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments for changes in the Earth's rotation. The value in `tt` is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth. Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET). |

#### member functions

<a name="Time.AddDays"></a>
### Time.AddDays(self, days)

**Calculates the sum or difference of a [`Time`](#Time) with a specified real-valued number of days.**

Sometimes we need to adjust a given [`Time`](#Time) value by a certain amount of time.
This function adds the given real number of days in `days` to the date and time
in the calling object.
More precisely, the result's Universal Time field `ut` is exactly adjusted by `days`
and the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
using a best-fit piecewise polynomial model devised by
[Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).
The value of the calling object is not modified. This function creates a brand new
[`Time`](#Time) object and returns it.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `days` | A floating point number of days by which to adjust `time`. May be negative, 0, or positive. |

**Returns**: [`Time`](#Time)

<a name="Time.FromTerrestrialTime"></a>
### Time.FromTerrestrialTime(tt)

**Creates a [`Time`](#Time) object from a Terrestrial Time day value.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `tt` | The number of days after the J2000 epoch. |

**Returns**: [`Time`](#Time)

<a name="Time.Make"></a>
### Time.Make(year, month, day, hour, minute, second)

**Creates a [`Time`](#Time) object from a UTC calendar date and time.**

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` | The UTC year value, e.g. 2019. |
| `int` | `month` | The UTC month in the range 1..12. |
| `int` | `day` | The UTC day of the month, in the range 1..31. |
| `int` | `hour` | The UTC hour, in the range 0..23. |
| `int` | `minute` | The UTC minute, in the range 0..59. |
| `float` | `second` | The real-valued UTC second, in the range [0, 60). |

**Returns**: [`Time`](#Time)

<a name="Time.Now"></a>
### Time.Now()

**Returns the computer's current date and time in the form of a [`Time`](#Time) object.**

Uses the computer's system clock to find the current UTC date and time.
Converts that date and time to a [`Time`](#Time) value and returns the result.
Callers can pass this value to other Astronomy Engine functions to
calculate current observational conditions.

**Returns**: [`Time`](#Time)

<a name="Time.Parse"></a>
### Time.Parse(text)

**Creates a [`Time`](#Time) object from a string of the form 'yyyy-mm-ddThh:mm:ss.sssZ'**

Parses a UTC date and time from a string and returns a [`Time`](#Time) object.
Permits a subset of ISO 8601 format.
The year, month, and day are required.
Hours, minutes, seconds, and fractions of a second are optional.
If time is specified, there must be a 'T' between the date and the time
and a 'Z' at the end of the time.

| Type | Parameter | Description |
| --- | --- | --- |
| `string` | `text` | A string of the following formats: `yyyy-mm-dd` `yyyy-mm-ddThh:mmZ` `yyyy-mm-ddThh:mm:ssZ` `yyyy-mm-ddThh:mm:ss.sssZ` |

**Returns**: [`Time`](#Time)

<a name="Time.Utc"></a>
### Time.Utc(self)

**Returns the UTC date and time as a `datetime` object.**

Uses the standard [`datetime`](https://docs.python.org/3/library/datetime.html) class
to represent the date and time in this Time object.

**Returns**: `datetime`

---

<a name="TransitInfo"></a>
### class TransitInfo

**Information about a transit of Mercury or Venus, as seen from the Earth.**

Returned by [`SearchTransit`](#SearchTransit) or [`NextTransit`](#NextTransit) to report
information about a transit of Mercury or Venus.
A transit is when Mercury or Venus passes between the Sun and Earth so that
the other planet is seen in silhouette against the Sun.
The calculations are performed from the point of view of a geocentric observer.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `start` | The date and time at the beginning of the transit. This is the moment the planet first becomes visible against the Sun in its background. |
| [`Time`](#Time) | `peak` | When the planet is most aligned with the Sun, as seen from the Earth. |
| [`Time`](#Time) | `finish` | The date and time at the end of the transit. This is the moment the planet is last seen against the Sun in its background. |
| `float` | `separation` | The minimum angular separation, in arcminutes, between the centers of the Sun and the planet. This angle pertains to the time stored in `peak`. |

---

<a name="Vector"></a>
### class Vector

**A Cartesian vector with 3 space coordinates and 1 time coordinate.**

The vector's space coordinates are measured in astronomical units (AU).
The coordinate system varies and depends on context.
The vector also includes a time stamp.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `x` | The x-coordinate of the vector, measured in AU. |
| `float` | `y` | The y-coordinate of the vector, measured in AU. |
| `float` | `z` | The z-coordinate of the vector, measured in AU. |
| [`Time`](#Time) | `t` | The date and time at which the coordinate is valid. |

#### member functions

<a name="Vector.Length"></a>
### Vector.Length(self)

Returns the length of the vector in AU.

<a name="Vector.format"></a>
### Vector.format(self, coord_format)

Returns a custom format string representation of the vector.

---

<a name="enumerations"></a>
## Enumerated Types

---

<a name="ApsisKind"></a>
### enum ApsisKind

**Represents whether a satellite is at a closest or farthest point in its orbit.**

An apsis is a point in a satellite's orbit that is closest to,
or farthest from, the body it orbits (its primary).
`ApsisKind` is an enumerated type that indicates which of these
two cases applies to a particular apsis event.

| Value | Description |
| --- | --- |
| `Pericenter` | The satellite is at its closest point to its primary. |
| `Apocenter` | The satellite is at its farthest point from its primary. |
| `Invalid` | A placeholder for an undefined, unknown, or invalid apsis. |

---

<a name="Body"></a>
### enum Body

**The celestial bodies supported by Astronomy Engine calculations.**

| Value | Description |
| --- | --- |
| `Invalid` | An unknown, invalid, or undefined celestial body. |
| `Mercury` | The planet Mercury. |
| `Venus` | The planet Venus. |
| `Earth` | The planet Earth. |
| `Mars` | The planet Mars. |
| `Jupiter` | The planet Jupiter. |
| `Saturn` | The planet Saturn. |
| `Uranus` | The planet Uranus. |
| `Neptune` | The planet Neptune. |
| `Pluto` | The planet Pluto. |
| `Sun` | The Sun. |
| `Moon` | The Earth's moon. |
| `EMB` | The Earth/Moon Barycenter. |
| `SSB` | The Solar System Barycenter. |
| `Star1` | User-defined star 1. |
| `Star2` | User-defined star 2. |
| `Star3` | User-defined star 3. |
| `Star4` | User-defined star 4. |
| `Star5` | User-defined star 5. |
| `Star6` | User-defined star 6. |
| `Star7` | User-defined star 7. |
| `Star8` | User-defined star 8. |

---

<a name="Direction"></a>
### enum Direction

**Indicates whether a body is rising above or setting below the horizon.**

Specifies the direction of a rising or setting event for a body.
For example, `Direction.Rise` is used to find sunrise times,
and `Direction.Set` is used to find sunset times.

| Value | Description |
| --- | --- |
| `Rise` | First appearance of a body as it rises above the horizon. |
| `Set` | Last appearance of a body as it sinks below the horizon. |

---

<a name="EclipseKind"></a>
### enum EclipseKind

**The different kinds of lunar/solar eclipses.**

| Value | Description |
| --- | --- |
| `Invalid` | No eclipse found. |
| `Penumbral` | A penumbral lunar eclipse. (Never used for a solar eclipse.) |
| `Partial` | A partial lunar/solar eclipse. |
| `Annular` | An annular solar eclipse. (Never used for a lunar eclipse.) |
| `Total` | A total lunar/solar eclipse. |

---

<a name="NodeEventKind"></a>
### enum NodeEventKind

**Indicates whether a crossing through the ecliptic plane is ascending or descending.**

| Value | Description |
| --- | --- |
| `Invalid` | A placeholder for an invalid or undefined node. |
| `Ascending` | indicates a body passing through the ecliptic plane from south to north. |
| `Descending` | indicates a body passing through the ecliptic plane from north to south. |

---

<a name="Refraction"></a>
### enum Refraction

**Selects if/how to correct for atmospheric refraction.**

Some functions allow enabling or disabling atmospheric refraction
for the calculated apparent position of a celestial body
as seen by an observer on the surface of the Earth.

| Value | Description |
| --- | --- |
| `Airless` | No atmospheric refraction correction. |
| `Normal` | Recommended correction for standard atmospheric refraction. |
| `JplHorizons` | Used only for compatibility testing with JPL Horizons online tool. |

---

<a name="Visibility"></a>
### enum Visibility

**Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.**

| Value | Description |
| --- | --- |
| `Morning` | The body is best visible in the morning, before sunrise. |
| `Evening` | The body is best visible in the evening, after sunset. |

---

<a name="errors"></a>
## Error Types

---

<a name="BadVectorError"></a>
### BadVectorError

A vector magnitude is too small to have a direction in space.

---

<a name="DateTimeFormatError"></a>
### DateTimeFormatError

The syntax of a UTC date/time string was not valid, or it contains invalid values.

---

<a name="EarthNotAllowedError"></a>
### EarthNotAllowedError

The Earth is not allowed as the celestial body in this calculation.

---

<a name="Error"></a>
### Error

Indicates an error in an astronomical calculation.

---

<a name="InternalError"></a>
### InternalError

**An internal error occured that should be reported as a bug.**

Indicates an unexpected and unrecoverable condition occurred.
If you encounter this error using Astronomy Engine, it would be very
helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
page on GitHub. Please include a copy of the stack trace, along with a description
of how to reproduce the error. This will help improve the quality of
Astronomy Engine for everyone! (Thank you in advance from the author.)

---

<a name="InvalidBodyError"></a>
### InvalidBodyError

The celestial body is not allowed for this calculation.

---

<a name="NoConvergeError"></a>
### NoConvergeError

**A numeric solver did not converge.**

Indicates that there was a failure of a numeric solver to converge.
If you encounter this error using Astronomy Engine, it would be very
helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
page on GitHub. Please include a copy of the stack trace, along with a description
of how to reproduce the error. This will help improve the quality of
Astronomy Engine for everyone! (Thank you in advance from the author.)

---

<a name="functions"></a>
## Functions

---

<a name="AngleBetween"></a>
### AngleBetween(a, b)

**Calculates the angle in degrees between two vectors.**

Given a pair of vectors, this function returns the angle in degrees
between the two vectors in 3D space.
The angle is measured in the plane that contains both vectors.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `a` | The first of a pair of vectors between which to measure an angle. |
| [`Vector`](#Vector) | `b` | The second of a pair of vectors between which to measure an angle. |

**Returns**: `float`
The angle between the two vectors expressed in degrees.
The value is in the range [0, 180].

---

<a name="AngleFromSun"></a>
### AngleFromSun(body, time)

**Returns the angle between the given body and the Sun, as seen from the Earth.**

This function calculates the angular separation between the given body and the Sun,
as seen from the center of the Earth. This angle is helpful for determining how
easy it is to see the body away from the glare of the Sun.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose angle from the Sun is to be measured. Not allowed to be `Body.Earth`. |
| [`Time`](#Time) | `time` | The time at which the observation is made. |

**Returns**: `float`
A numeric value indicating the angle in degrees between the Sun
and the specified body as seen from the center of the Earth.

---

<a name="BackdatePosition"></a>
### BackdatePosition(time, observerBody, targetBody, aberration)

**Solve for light travel time correction of apparent position.**

When observing a distant object, for example Jupiter as seen from Earth,
the amount of time it takes for light to travel from the object to the
observer can significantly affect the object's apparent position.
This function solves the light travel time correction for the apparent
relative position vector of a target body as seen by an observer body
at a given observation time.
For geocentric calculations, [`GeoVector`](#GeoVector) also includes light
travel time correction, but the time `t` embedded in its returned vector
refers to the observation time, not the backdated time that light left
the observed body. Thus `BackdatePosition` provides direct
access to the light departure time for callers that need it.
For a more generalized light travel correction solver, see [`CorrectLightTravel`](#CorrectLightTravel).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The time of observation. |
| [`Body`](#Body) | `observerBody` | The body to be used as the observation location. |
| [`Body`](#Body) | `targetBody` | The body to be observed. |
| `bool` | `aberration` | `True` to correct for aberration, or `False` to leave uncorrected. |

**Returns**: [`Vector`](#Vector)
The position vector at the solved backdated time.
Its `t` field holds the time that light left the observed
body to arrive at the observer at the observation time.

---

<a name="BaryState"></a>
### BaryState(body, time)

**Calculates barycentric position and velocity vectors for the given body.**

Given a body and a time, calculates the barycentric position and velocity
vectors for the center of that body at that time.
The vectors are expressed in equatorial J2000 coordinates (EQJ).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose barycentric state vector is to be calculated. Supported values are `Body.Sun`, `Body.SSB`, `Body.Moon`, `Body.EMB`, and all planets: `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. |
| [`Time`](#Time) | `time` | The date and time for which to calculate position and velocity. |

**Returns**: [`StateVector`](#StateVector)
An object that contains barycentric position and velocity vectors.

---

<a name="BodyCode"></a>
### BodyCode(name)

**Finds the Body enumeration value, given the name of a body.**

```
>>> astronomy.BodyCode('Mars')
<Body.Mars: 3>
```

| Type | Parameter | Description |
| --- | --- | --- |
| `str` | `name` | The common English name of a supported celestial body. |

**Returns**: [`Body`](#Body)
If `name` is a valid body name, returns the enumeration
value associated with that body.
Otherwise, returns `Body.Invalid`.

---

<a name="CombineRotation"></a>
### CombineRotation(a, b)

**Creates a rotation based on applying one rotation followed by another.**

Given two rotation matrices, returns a combined rotation matrix that is
equivalent to rotating based on the first matrix, followed by the second.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `a` | The first rotation to apply. |
| [`RotationMatrix`](#RotationMatrix) | `b` | The second rotation to apply. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
The combined rotation matrix.

---

<a name="Constellation"></a>
### Constellation(ra, dec)

**Determines the constellation that contains the given point in the sky.**

Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
constellation that contains that point.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ra` | The right ascension (RA) of a point in the sky, using the J2000 equatorial system. |
| `float` | `dec` | The declination (DEC) of a point in the sky, using the J2000 equatorial system. |

**Returns**: [`ConstellationInfo`](#ConstellationInfo)
A structure that contains the 3-letter abbreviation and full name
of the constellation that contains the given (ra,dec), along with
the converted B1875 (ra,dec) for that point.

---

<a name="CorrectLightTravel"></a>
### CorrectLightTravel(func, time)

**Solve for light travel time of a vector function.**

When observing a distant object, for example Jupiter as seen from Earth,
the amount of time it takes for light to travel from the object to the
observer can significantly affect the object's apparent position.
This function is a generic solver that figures out how long in the
past light must have left the observed object to reach the observer
at the specified observation time. It uses [`PositionFunction`](#PositionFunction)
to express an arbitrary position vector as a function of time.
This function repeatedly calls `func.Position`, passing a series of time
estimates in the past. Then `func.Position` must return a relative state vector between
the observer and the target. `CorrectLightTravel` keeps calling
`func.Position` with more and more refined estimates of the time light must have
left the target to arrive at the observer.
For common use cases, it is simpler to use [`BackdatePosition`](#BackdatePosition)
for calculating the light travel time correction of one body observing another body.
For geocentric calculations, [`GeoVector`](#GeoVector) also backdates the returned
position vector for light travel time, only it returns the observation time in
the returned vector's `t` field rather than the backdated time.
time : Time
    The observation time for which to solve for light travel delay.

| Type | Parameter | Description |
| --- | --- | --- |
| [`PositionFunction`](#PositionFunction) | `func` | An arbitrary position vector as a function of time. |

**Returns**: [`Vector`](#Vector)
The position vector at the solved backdated time.
The `t` field holds the time that light left the observed
body to arrive at the observer at the observation time.

---

<a name="DefineStar"></a>
### DefineStar(body, ra, dec, distanceLightYears)

**Assign equatorial coordinates to a user-defined star.**

Some Astronomy Engine functions allow their `body` parameter to
be a user-defined fixed point in the sky, loosely called a "star".
This function assigns a right ascension, declination, and distance
to one of the eight user-defined stars `Body.Star1`..`Body.Star8`.
Stars are not valid until defined. Once defined, they retain their
definition until re-defined by another call to `DefineStar`.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | One of the eight user-defined star identifiers: `Body.Star1`, `Body.Star2`, ..., `Body.Star8`. |
| `float` | `ra` | The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of sidereal hours, and must be within the half-open range [0, 24). |
| `float` | `dec` | The declination to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of degrees north (positive) or south (negative) of the J2000 equator, and must be within the closed range [-90, +90]. |
| `float` | `distanceLightYears` | The distance between the star and the Sun, expressed in light-years. This value is used to calculate the tiny parallax shift as seen by an observer on Earth. If you don't know the distance to the star, using a large value like 1000 will generally work well. The minimum allowed distance is 1 light-year, which is required to provide certain internal optimizations. |

---

<a name="DeltaT_EspenakMeeus"></a>
### DeltaT_EspenakMeeus(ut)

**The default Delta T function used by Astronomy Engine.**

Espenak and Meeus use a series of piecewise polynomials to
approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
This is the default Delta T function used by Astronomy Engine.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ut` | The floating point number of days since noon UTC on January 1, 2000. |

**Returns**: `float`
The estimated difference TT-UT on the given date, expressed in seconds.

---

<a name="Ecliptic"></a>
### Ecliptic(equ)

**Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.**

Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
which are relative to the plane of the Earth's orbit around the Sun.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Equatorial`](#Equatorial) | `equ` | Equatorial coordinates in the J2000 frame of reference. |

**Returns**: [`EclipticCoordinates`](#EclipticCoordinates)
Ecliptic coordinates in the J2000 frame of reference.

---

<a name="EclipticGeoMoon"></a>
### EclipticGeoMoon(time)

**Calculates spherical ecliptic geocentric position of the Moon.**

Given a time of observation, calculates the Moon's geocentric position
in ecliptic spherical coordinates. Provides the ecliptic latitude and
longitude in degrees, and the geocentric distance in astronomical units (AU).
The ecliptic longitude is measured relative to the equinox of date.
This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
It is adapted from Turbo Pascal code from the book
[Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
by Montenbruck and Pfleger.
To calculate an equatorial J2000 vector instead, use [`GeoMoon`](#GeoMoon).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Moon's position. |

**Returns**: [`Spherical`](#Spherical)
The Moon's position as a distance, ecliptic latitude, and ecliptic longitude.

---

<a name="EclipticLongitude"></a>
### EclipticLongitude(body, time)

**Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.**

This function calculates the angle around the plane of the Earth's orbit
of a celestial body, as seen from the center of the Sun.
The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body other than the Sun. |
| [`Time`](#Time) | `time` | The date and time at which the body's ecliptic longitude is to be calculated. |

**Returns**: `float`
An angular value in degrees indicating the ecliptic longitude of the body.

---

<a name="Elongation"></a>
### Elongation(body, time)

**Determines visibility of a celestial body relative to the Sun, as seen from the Earth.**

This function returns an [`ElongationEvent`](#ElongationEvent) object, which provides the following
information about the given celestial body at the given time:
- `visibility` is an enumerated type that specifies whether the body is more
  easily seen in the morning before sunrise, or in the evening after sunset.
- `elongation` is the angle in degrees between two vectors: one from the center
  of the Earth to the center of the Sun, the other from the center of the Earth
  to the center of the specified body. This angle indicates how far away the body
  is from the glare of the Sun. The elongation angle is always in the range [0, 180].
- `ecliptic_separation` is the absolute value of the difference between the body's
  ecliptic longitude and the Sun's ecliptic longitude, both as seen from the center
  of the Earth. This angle measures around the plane of the Earth's orbit, and ignores
  how far above or below that plane the body is.
  The ecliptic separation is measured in degrees and is always in the range [0, 180].

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose visibility is to be calculated. |
| [`Time`](#Time) | `time` | The date and time of the observation. |

**Returns**: [`ElongationEvent`](#ElongationEvent)

---

<a name="Equator"></a>
### Equator(body, time, observer, ofdate, aberration)

**Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.**

Calculates topocentric equatorial coordinates in one of two different systems:
J2000 or true-equator-of-date, depending on the value of the `ofdate` parameter.
Equatorial coordinates include right ascension, declination, and distance in astronomical units.
This function corrects for light travel time: it adjusts the apparent location
of the observed body based on how long it takes for light to travel from the body to the Earth.
This function corrects for *topocentric parallax*, meaning that it adjusts for the
angular shift depending on where the observer is located on the Earth. This is most
significant for the Moon, because it is so close to the Earth. However, parallax corection
has a small effect on the apparent positions of other bodies.
Correction for aberration is optional, using the `aberration` parameter.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body to be observed. Not allowed to be `Body.Earth`. |
| [`Time`](#Time) | `time` | The date and time at which the observation takes place. |
| [`Observer`](#Observer) | `observer` | A location on or near the surface of the Earth. |
| `bool` | `ofdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. If `True`, returns coordinates using the equator and equinox of date. If `False`, returns coordinates converted to the J2000 system. |
| `bool` | `aberration` | If `True`, corrects for aberration of light based on the motion of the Earth with respect to the heliocentric origin. If `False`, does not correct for aberration. |

**Returns**: [`Equatorial`](#Equatorial)
Equatorial coordinates in the specified frame of reference.

---

<a name="EquatorFromVector"></a>
### EquatorFromVector(vec)

**Given an equatorial vector, calculates equatorial angular coordinates.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vec` | A vector in an equatorial coordinate system. |

**Returns**: [`Equatorial`](#Equatorial)
Angular coordinates expressed in the same equatorial system as `vec`.

---

<a name="GeoEmbState"></a>
### GeoEmbState(time)

**Calculates the geocentric position and velocity of the Earth/Moon barycenter.**

Given a time of observation, calculates the geocentric position and velocity vectors
of the Earth/Moon barycenter (EMB).
The position (x, y, z) components are expressed in AU (astronomical units).
The velocity (vx, vy, vz) components are expressed in AU/day.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the EMB's geocentric state. |

**Returns**: [`StateVector`](#StateVector)
The EMB's position and velocity vectors in J2000 equatorial coordinates.

---

<a name="GeoMoon"></a>
### GeoMoon(time)

**Calculates equatorial geocentric position of the Moon at a given time.**

Given a time of observation, calculates the Moon's position as a vector.
The vector gives the location of the Moon's center relative to the Earth's center
with x-, y-, and z-components measured in astronomical units.
The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
In Astronomy Engine, this orientation is called EQJ.
This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
It is adapted from Turbo Pascal code from the book
[Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
by Montenbruck and Pfleger.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Moon's position. |

**Returns**: [`Vector`](#Vector)
The Moon's position as a vector in J2000 Cartesian equatorial coordinates (EQJ).

---

<a name="GeoMoonState"></a>
### GeoMoonState(time)

**Calculates equatorial geocentric position and velocity of the Moon at a given time.**

Given a time of observation, calculates the Moon's position and velocity vectors.
The position and velocity are of the Moon's center relative to the Earth's center.
The position (x, y, z) components are expressed in AU (astronomical units).
The velocity (vx, vy, vz) components are expressed in AU/day.
The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
In Astronomy Engine, this orientation is called EQJ.
If you need the Moon's position only, and not its velocity,
it is much more efficient to use [`GeoMoon`](#GeoMoon) instead.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Moon's position and velocity. |

**Returns**: [`StateVector`](#StateVector)
The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).

---

<a name="GeoVector"></a>
### GeoVector(body, time, aberration)

**Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Earth as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.
If given an invalid value for `body`, this function will raise an exception.
Unlike [`HelioVector`](#HelioVector), this function corrects for light travel time.
This means the position of the body is "back-dated" by the amount of time it takes
light to travel from that body to an observer on the Earth.
Also, the position can optionally be corrected for
[aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
causing the apparent direction of the body to be shifted due to transverse
movement of the Earth with respect to the rays of light coming from that body.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets. |
| [`Time`](#Time) | `time` | The date and time for which to calculate the position. |
| `bool` | `aberration` | A boolean value indicating whether to correct for aberration. |

**Returns**: [`Vector`](#Vector)
A geocentric position vector of the center of the given body.

---

<a name="HelioDistance"></a>
### HelioDistance(body, time)

**Calculates the distance between a body and the Sun at a given time.**

Given a date and time, this function calculates the distance between
the center of `body` and the center of the Sun.
For the planets Mercury through Neptune, this function is significantly
more efficient than calling [`HelioVector`](#HelioVector) followed by taking the length
of the resulting vector.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric distance: the Sun, Moon, any of the planets, or a user-defined star. |
| [`Time`](#Time) | `time` | The date and time for which to calculate the heliocentric distance. |

**Returns**: `float`
The heliocentric distance in AU.

---

<a name="HelioState"></a>
### HelioState(body, time)

**Calculates heliocentric position and velocity vectors for the given body.**

Given a body and a time, calculates the position and velocity
vectors for the center of that body at that time, relative to the center of the Sun.
The vectors are expressed in equatorial J2000 coordinates (EQJ).
If you need the position vector only, it is more efficient to call [`HelioVector`](#HelioVector).
The Sun's center is a non-inertial frame of reference. In other words, the Sun
experiences acceleration due to gravitational forces, mostly from the larger
planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
kinetic energy, or other quantities that require a non-accelerating frame
of reference, consider using [`BaryState`](#BaryState) instead.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose heliocentric state vector is to be calculated. Supported values are `Body.Sun`, `Body.SSB`, `Body.Moon`, `Body.EMB`, and all planets: `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. Also allowed to be a user-defined star created by [`DefineStar`](#DefineStar). |
| [`Time`](#Time) | `time` | The date and time for which to calculate position and velocity. |

**Returns**: [`StateVector`](#StateVector)
An object that contains heliocentric position and velocity vectors.

---

<a name="HelioVector"></a>
### HelioVector(body, time)

**Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Sun as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.
The position is not corrected for light travel time or aberration.
This is different from the behavior of [`GeoVector`](#GeoVector).
If given an invalid value for `body`, this function raises an exception.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose heliocentric position is to be calculated: The Sun, Moon, EMB, SSB, or any of the planets. Also allowed to be a user-defined star created by [`DefineStar`](#DefineStar). |
| [`Time`](#Time) | `time` | The time at which to calculate the heliocentric position. |

**Returns**: [`Vector`](#Vector)
A heliocentric position vector of the center of the given body
at the given time.

---

<a name="Horizon"></a>
### Horizon(time, observer, ra, dec, refraction)

**Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.**

Given a date and time, the geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial body,
this function returns horizontal coordinates (azimuth and altitude angles) for the body
relative to the horizon at the geographic location.
The right ascension `ra` and declination `dec` passed in must be *equator of date*
coordinates, based on the Earth's true equator at the date and time of the observation.
Otherwise the resulting horizontal coordinates will be inaccurate.
Equator of date coordinates can be obtained by calling [`Equator`](#Equator), passing in
`True` as its `ofdate` parameter. It is also recommended to enable
aberration correction by passing in `True` for the `aberration` parameter.
This function optionally corrects for atmospheric refraction.
For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
correct for optical lensing of the Earth's atmosphere that causes objects
to appear somewhat higher above the horizon than they actually are.
However, callers may choose to avoid this correction by passing in `Refraction.Airless`.
If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
in the [`HorizontalCoordinates`](#HorizontalCoordinates) object returned by this function will all be corrected for refraction.
If refraction is disabled, none of these four coordinates will be corrected; in that case,
the right ascension and declination in the returned object will be numerically identical
to the respective `ra` and `dec` values passed in.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to find horizontal coordinates. |
| [`Observer`](#Observer) | `observer` | The location of the observer for which to find horizontal coordinates. |
| `float` | `ra` | Right ascension in sidereal hours of the celestial object, referred to the mean equinox of date for the J2000 epoch. |
| `float` | `dec` | Declination in degrees of the celestial object, referred to the mean equator of date for the J2000 epoch. Positive values are north of the celestial equator and negative values are south of it. |
| [`Refraction`](#Refraction) | `refraction` | The option for selecting whether to correct for atmospheric lensing. If `Refraction.Normal`, a well-behaved refraction model is used. If `Refraction.Airless`, no refraction correct is performed. `Refraction.JplHorizons` is used only for compatibility testing with the JPL Horizons online tool. |

**Returns**: [`HorizontalCoordinates`](#HorizontalCoordinates)
The horizontal coordinates (altitude and azimuth), along with
equatorial coordinates (right ascension and declination), all
optionally corrected for atmospheric refraction. See remarks above
for more details.

---

<a name="HorizonFromVector"></a>
### HorizonFromVector(vector, refraction)

**Converts Cartesian coordinates to horizontal coordinates.**

Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
*IMPORTANT:* This function differs from `SphereFromVector` in two ways:
- `SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
  from north (e.g., west = +90), but this function represents a clockwise rotation
  (e.g., east = +90). The difference is because `SphereFromVector` is intended
  to preserve the vector "right-hand rule", while this function defines azimuth in a more
  traditional way as used in navigation and cartography.
- This function optionally corrects for atmospheric refraction, while `SphereFromVector` does not.
The returned object contains the azimuth in `lon`.
It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
The altitude is stored in `lat`.
The distance to the observed object is stored in `dist`,
and is expressed in astronomical units (AU).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vector` | Cartesian vector to be converted to horizontal angular coordinates. |
| [`Refraction`](#Refraction) | `refraction` | See comments in the [`RefractionAngle`](#RefractionAngle) function. |

---

<a name="IdentityMatrix"></a>
### IdentityMatrix()

**Creates an identity rotation matrix.**

Returns a rotation matrix that has no effect on orientation.
This matrix can be the starting point for other operations,
such as using a series of calls to [`Pivot`](#Pivot) to
create a custom rotation matrix.

**Returns**: [`RotationMatrix`](#RotationMatrix)
The identity rotation matrix.

---

<a name="Illumination"></a>
### Illumination(body, time)

**Finds visual magnitude, phase angle, and other illumination information about a celestial body.**

This function calculates information about how bright a celestial body appears from the Earth,
reported as visual magnitude, which is a smaller (or even negative) number for brighter objects,
and a larger number for dimmer objects.
For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
of the body appears illuminated as seen from the Earth. For example, when the phase angle is
near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
of 90 degrees means the body appears "half full".
For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
so it doesn't have a phase angle.
When the body is Saturn, the returned object contains a field `ring_tilt` that holds
the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
0 for all bodies other than Saturn.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, or any planet other than the Earth. |
| [`Time`](#Time) | `time` | The date and time of the observation. |

**Returns**: [`IlluminationInfo`](#IlluminationInfo)

---

<a name="InverseRefractionAngle"></a>
### InverseRefractionAngle(refraction, bent_altitude)

**Calculates the inverse of an atmospheric refraction angle.**

Given an observed altitude angle that includes atmospheric refraction,
calculates the negative angular correction to obtain the unrefracted
altitude. This is useful for cases where observed horizontal
coordinates are to be converted to another orientation system,
but refraction first must be removed from the observed position.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Refraction`](#Refraction) | `refraction` | `Refraction.Normal` - corrects for atmospheric refraction (recommended). `Refraction.Airless` - no correction is performed. `Refraction.JplHorizons` - For JPL Horizons compatibility testing only. |
| `float` | `bent_altitude` | The apparent altitude that includes atmospheric refraction. |

**Returns**: `float`
The angular adjustment in degrees, to be added to the
altitude angle to correct for atmospheric lensing.
This will be less than or equal to zero.

---

<a name="InverseRotation"></a>
### InverseRotation(rotation)

**Calculates the inverse of a rotation matrix.**

Given a rotation matrix that performs some coordinate transform,
this function returns the matrix that reverses that transform.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | The rotation matrix to be inverted. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
The inverse rotation matrix.

---

<a name="JupiterMoons"></a>
### JupiterMoons(time)

**Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.**

Calculates position and velocity vectors for Jupiter's moons
Io, Europa, Ganymede, and Callisto, at the given date and time.
The vectors are jovicentric (relative to the center of Jupiter).
Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
The position components are expressed in astronomical units (AU), and the
velocity components are in AU/day.
To convert to heliocentric vectors, call [`HelioVector`](#HelioVector)
with `Body.Jupiter` to get Jupiter's heliocentric position, then
add the jovicentric vectors. Likewise, you can call [`GeoVector`](#GeoVector)
to convert to geocentric vectors.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate Jupiter's moons. |

**Returns**: [`JupiterMoonsInfo`](#JupiterMoonsInfo)
The positions and velocities of Jupiter's 4 largest moons.

---

<a name="LagrangePoint"></a>
### LagrangePoint(point, time, major_body, minor_body)

**Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.**

Given a more massive "major" body and a much less massive "minor" body,
calculates one of the five Lagrange points in relation to the minor body's
orbit around the major body. The parameter `point` is an integer that
selects the Lagrange point as follows:
1 = the Lagrange point between the major body and minor body.
2 = the Lagrange point on the far side of the minor body.
3 = the Lagrange point on the far side of the major body.
4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
5 = the Lagrange point 60 degrees behind the minor body's orbital position.
The function returns the state vector for the selected Lagrange point
in equatorial J2000 coordinates (EQJ), with respect to the center of the
major body.
To calculate Sun/Earth Lagrange points, pass in `Body.Sun` for `major_body`
and `Body.EMB` (Earth/Moon barycenter) for `minor_body`.
For Lagrange points of the Sun and any other planet, pass in just that planet
(e.g. `Body.Jupiter`) for `minor_body`.
To calculate Earth/Moon Lagrange points, pass in `Body.Earth` and `Body.Moon`
for the major and minor bodies respectively.
In some cases, it may be more efficient to call [`LagrangePointFast`](#LagrangePointFast),
especially when the state vectors have already been calculated, or are needed
for some other purpose.

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `point` | An integer 1..5 that selects which of the Lagrange points to calculate. |
| [`Time`](#Time) | `time` | The time for which the Lagrange point is to be calculated. |
| [`Body`](#Body) | `major_body` | The more massive of the co-orbiting bodies: `Body.Sun` or `Body.Earth`. |
| [`Body`](#Body) | `minor_body` | The less massive of the co-orbiting bodies. See main remarks. |

**Returns**: [`StateVector`](#StateVector)
The position and velocity of the selected Lagrange point with respect to the major body's center.

---

<a name="LagrangePointFast"></a>
### LagrangePointFast(point, major_state, major_mass, minor_state, minor_mass)

**Calculates one of the 5 Lagrange points from body masses and state vectors.**

Given a more massive "major" body and a much less massive "minor" body,
calculates one of the five Lagrange points in relation to the minor body's
orbit around the major body. The parameter `point` is an integer that
selects the Lagrange point as follows:
1 = the Lagrange point between the major body and minor body.
2 = the Lagrange point on the far side of the minor body.
3 = the Lagrange point on the far side of the major body.
4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
5 = the Lagrange point 60 degrees behind the minor body's orbital position.
The caller passes in the state vector and mass for both bodies.
The state vectors can be in any orientation and frame of reference.
The body masses are expressed as GM products, where G = the universal
gravitation constant and M = the body's mass. Thus the units for
`major_mass` and `minor_mass` must be au^3/day^2.
Use [`MassProduct`](#MassProduct) to obtain GM values for various solar system bodies.
The function returns the state vector for the selected Lagrange point
using the same orientation as the state vector parameters `major_state` and `minor_state`,
and the position and velocity components are with respect to the major body's center.
Consider calling [`LagrangePoint`](#LagrangePoint), instead of this function, for simpler usage in most cases.

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `point` | An integer 1..5 that selects which of the Lagrange points to calculate. |
| [`StateVector`](#StateVector) | `major_state` | The state vector of the major (more massive) of the pair of bodies. |
| `float` | `major_mass` | The mass product GM of the major body. |
| [`StateVector`](#StateVector) | `minor_state` | The state vector of the minor (less massive) of the pair of bodies. |
| `float` | `minor_mass` | The mass product GM of the minor body. |

**Returns**: [`StateVector`](#StateVector)
The position and velocity of the selected Lagrange point with respect to the major body's center.

---

<a name="Libration"></a>
### Libration(time)

**Calculates the Moon's libration angles at a given moment in time.**

Libration is an observed back-and-forth wobble of the portion of the
Moon visible from the Earth. It is caused by the imperfect tidal locking
of the Moon's fixed rotation rate, compared to its variable angular speed
of orbit around the Earth.
This function calculates a pair of perpendicular libration angles,
one representing rotation of the Moon in eclitpic longitude `elon`, the other
in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.
This function also returns the geocentric position of the Moon
expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
and the apparent angular diameter of the Moon `diam_deg`.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Moon's libration angles. |

**Returns**: [`LibrationInfo`](#LibrationInfo)

---

<a name="MassProduct"></a>
### MassProduct(body)

**Returns the product of mass and universal gravitational constant of a Solar System body.**

For problems involving the gravitational interactions of Solar System bodies,
it is helpful to know the product GM, where G = the universal gravitational constant
and M = the mass of the body. In practice, GM is known to a higher precision than
either G or M alone, and thus using the product results in the most accurate results.
This function returns the product GM in the units au^3/day^2.
The values come from page 10 of a
[JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The body for which to find the GM product. Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet. Any other value will cause an exception to be thrown. |

**Returns**: `float`
The mass product of the given body in au^3/day^2.

---

<a name="MoonPhase"></a>
### MoonPhase(time)

**Returns the Moon's phase as an angle from 0 to 360 degrees.**

This function determines the phase of the Moon using its apparent
ecliptic longitude relative to the Sun, as seen from the center of the Earth.
Certain values of the angle have conventional definitions:
- 0 = new moon
- 90 = first quarter
- 180 = full moon
- 270 = third quarter

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |

**Returns**: `float`

---

<a name="NextGlobalSolarEclipse"></a>
### NextGlobalSolarEclipse(prevEclipseTime)

**Searches for the next global solar eclipse in a series.**

After using [`SearchGlobalSolarEclipse`](#SearchGlobalSolarEclipse) to find the first solar eclipse
in a series, you can call this function to find the next consecutive solar eclipse.
Pass in the `peak` value from the [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo) returned by the
previous call to `SearchGlobalSolarEclipse` or `NextGlobalSolarEclipse`
to find the next solar eclipse.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `prevEclipseTime` | A date and time near a new moon. Solar eclipse search will start at the next new moon. |

**Returns**: [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo)

---

<a name="NextLocalSolarEclipse"></a>
### NextLocalSolarEclipse(prevEclipseTime, observer)

**Searches for the next local solar eclipse in a series.**

After using [`SearchLocalSolarEclipse`](#SearchLocalSolarEclipse) to find the first solar eclipse
in a series, you can call this function to find the next consecutive solar eclipse.
Pass in the `peak` value from the [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo) returned by the
previous call to `SearchLocalSolarEclipse` or `NextLocalSolarEclipse`
to find the next solar eclipse.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `prevEclipseTime` | A date and time near a new moon. Solar eclipse search will start at the next new moon. |
| [`Observer`](#Observer) | `observer` | The geographic location of the observer. |

**Returns**: [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo)

---

<a name="NextLunarApsis"></a>
### NextLunarApsis(apsis)

**Finds the next lunar perigee or apogee in a series.**

This function requires an [`Apsis`](#Apsis) value obtained from a call to
[`SearchLunarApsis`](#SearchLunarApsis) or `NextLunarApsis`.
Given an apogee event, this function finds the next perigee event,
and vice versa.
See [`SearchLunarApsis`](#SearchLunarApsis) for more details.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Apsis`](#Apsis) | `apsis` |  |

**Returns**: [`Apsis`](#Apsis)

---

<a name="NextLunarEclipse"></a>
### NextLunarEclipse(prevEclipseTime)

**Searches for the next lunar eclipse in a series.**

 After using [`SearchLunarEclipse`](#SearchLunarEclipse) to find the first lunar eclipse
 in a series, you can call this function to find the next consecutive lunar eclipse.
 Pass in the `peak` value from the [`LunarEclipseInfo`](#LunarEclipseInfo) returned by the
 previous call to `SearchLunarEclipse` or `NextLunarEclipse`
 to find the next lunar eclipse.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `prevEclipseTime` | A date and time near a full moon. Lunar eclipse search will start at the next full moon. |

**Returns**: [`LunarEclipseInfo`](#LunarEclipseInfo)

---

<a name="NextMoonNode"></a>
### NextMoonNode(prevNode)

**Searches for the next time when the Moon's center crosses through the ecliptic plane.**

Call [`SearchMoonNode`](#SearchMoonNode) to find the first of a series of nodes.
Then call `NextMoonNode` to find as many more consecutive nodes as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`NodeEventInfo`](#NodeEventInfo) | `prevNode` | The previous node find from calling [`SearchMoonNode`](#SearchMoonNode) or `NextMoonNode`. |

**Returns**: [`NodeEventInfo`](#NodeEventInfo)

---

<a name="NextMoonQuarter"></a>
### NextMoonQuarter(mq)

**Continues searching for lunar quarters from a previous search.**

After calling [`SearchMoonQuarter`](#SearchMoonQuarter), this function can be called
one or more times to continue finding consecutive lunar quarters.
This function finds the next consecutive moon quarter event after
the one passed in as the parameter `mq`.

| Type | Parameter | Description |
| --- | --- | --- |
| [`MoonQuarter`](#MoonQuarter) | `mq` | A value returned by a prior call to [`SearchMoonQuarter`](#SearchMoonQuarter) or [`NextMoonQuarter`](#NextMoonQuarter). |

**Returns**: [`MoonQuarter`](#MoonQuarter)

---

<a name="NextPlanetApsis"></a>
### NextPlanetApsis(body, apsis)

**Finds the next planetary perihelion or aphelion event in a series.**

This function requires an [`Apsis`](#Apsis) value obtained from a call
to [`SearchPlanetApsis`](#SearchPlanetApsis) or `NextPlanetApsis`.
Given an aphelion event, this function finds the next perihelion event, and vice versa.
See [`SearchPlanetApsis`](#SearchPlanetApsis) for more details.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find the next perihelion/aphelion event. Not allowed to be `Body.Sun` or `Body.Moon`. Must match the body passed into the call that produced the `apsis` parameter. |
| [`Apsis`](#Apsis) | `apsis` | An apsis event obtained from a call to [`SearchPlanetApsis`](#SearchPlanetApsis) or `NextPlanetApsis`. |

**Returns**: [`Apsis`](#Apsis)

---

<a name="NextTransit"></a>
### NextTransit(body, prevTransitTime)

**Searches for another transit of Mercury or Venus.**

After calling [`SearchTransit`](#SearchTransit) to find a transit of Mercury or Venus,
this function finds the next transit after that.
Keep calling this function as many times as you want to keep finding more transits.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`. |
| [`Time`](#Time) | `prevTransitTime` | A date and time near the previous transit. |

**Returns**: [`TransitInfo`](#TransitInfo)

---

<a name="ObserverGravity"></a>
### ObserverGravity(latitude, height)

**Calculates the gravitational acceleration experienced by an observer on the Earth.**

This function implements the WGS 84 Ellipsoidal Gravity Formula.
The result is a combination of inward gravitational acceleration
with outward centrifugal acceleration, as experienced by an observer
in the Earth's rotating frame of reference.
The resulting value increases toward the Earth's poles and decreases
toward the equator, consistent with changes of the weight measured
by a spring scale of a fixed mass moved to different latitudes and heights
on the Earth.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `latitude` | The latitude of the observer in degrees north or south of the equator. By formula symmetry, positive latitudes give the same answer as negative latitudes, so the sign does not matter. |
| `float` | `height` | The height above the sea level geoid in meters. No range checking is done; however, accuracy is only valid in the range 0 to 100000 meters. |

**Returns**: `float`
The effective gravitational acceleration expressed in meters per second squared [m/s^2].

---

<a name="ObserverState"></a>
### ObserverState(time, observer, ofdate)

**Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.**

This function calculates position and velocity vectors of an observer
on or near the surface of the Earth, expressed in equatorial
coordinates. It takes into account the rotation of the Earth at the given
time, along with the given latitude, longitude, and elevation of the observer.
The caller may pass `ofdate` as `True` to return coordinates relative to the Earth's
equator at the specified time, or `False` to use the J2000 equator.
The returned position vector has components expressed in astronomical units (AU).
To convert to kilometers, multiply the `x`, `y`, and `z` values by
the constant value [`KM_PER_AU`](#KM_PER_AU).
The returned velocity vector has components expressed in AU/day.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the observer's position and velocity vectors. |
| [`Observer`](#Observer) | `observer` | The geographic location of a point on or near the surface of the Earth. |
| `bool` | `ofdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. The caller may pass `False` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the `time` parameter. Or the caller may pass `True` to use the Earth's equator at `time` as the orientation. |

**Returns**: [`StateVector`](#StateVector)
An equatorial position vector and velocity vector relative to the center of the Earth.

---

<a name="ObserverVector"></a>
### ObserverVector(time, observer, ofdate)

**Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.**

This function calculates a vector from the center of the Earth to
a point on or near the surface of the Earth, expressed in equatorial
coordinates. It takes into account the rotation of the Earth at the given
time, along with the given latitude, longitude, and elevation of the observer.
The caller may pass `ofdate` as `True` to return coordinates relative to the Earth's
equator at the specified time, or `False` to use the J2000 equator.
The returned vector has components expressed in astronomical units (AU).
To convert to kilometers, multiply the `x`, `y`, and `z` values by
the constant value [`KM_PER_AU`](#KM_PER_AU).
The inverse of this function is also available: [`VectorObserver`](#VectorObserver).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the observer's position vector. |
| [`Observer`](#Observer) | `observer` | The geographic location of a point on or near the surface of the Earth. |
| `bool` | `ofdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. The caller may pass `False` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the `time` parameter. Or the caller may pass `True` to use the Earth's equator at `time` as the orientation. |

**Returns**: [`Vector`](#Vector)
An equatorial vector from the center of the Earth to the specified location
on (or near) the Earth's surface.

---

<a name="PairLongitude"></a>
### PairLongitude(body1, body2, time)

**Returns one body's ecliptic longitude with respect to another, as seen from the Earth.**

This function determines where one body appears around the ecliptic plane
(the plane of the Earth's orbit around the Sun) as seen from the Earth,
relative to the another body's apparent position.
The function returns an angle in the half-open range [0, 360) degrees.
The value is the ecliptic longitude of `body1` relative to the ecliptic
longitude of `body2`.
The angle is 0 when the two bodies are at the same ecliptic longitude
as seen from the Earth. The angle increases in the prograde direction
(the direction that the planets orbit the Sun and the Moon orbits the Earth).
When the angle is 180 degrees, it means the two bodies appear on opposite sides
of the sky for an Earthly observer.
Neither `body1` nor `body2` is allowed to be `Body.Earth`.
If this happens, the function throws an exception.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body1` | The first body, whose longitude is to be found relative to the second body. |
| [`Body`](#Body) | `body2` | The second body, relative to which the longitude of the first body is to be found. |
| [`Time`](#Time) | `time` | The date and time of the observation. |

**Returns**: `float`
An angle in degrees in the range [0, 360).

---

<a name="Pivot"></a>
### Pivot(rotation, axis, angle)

**Re-orients a rotation matrix by pivoting it by an angle around one of its axes.**

Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
this function pivots the rotation matrix by that angle around that coordinate axis.
For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
of a telescope camera pointed at a given body, you can use `Pivot` twice:
(1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
western axis by the body's altitude angle. The resulting rotation matrix will then
reorient ECL coordinates to the orientation of your telescope camera.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | The input rotation matrix. |
| `int` | `axis` | An integer that selects which coordinate axis to rotate around: 0 = x, 1 = y, 2 = z. Any other value will cause an exception. |
| `float` | `angle` | An angle in degrees indicating the amount of rotation around the specified axis. Positive angles indicate rotation counterclockwise as seen from the positive direction along that axis, looking towards the origin point of the orientation system. Any finite number of degrees is allowed, but best precision will result from keeping `angle` in the range [-360, +360]. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A pivoted matrix object.

---

<a name="PlanetOrbitalPeriod"></a>
### PlanetOrbitalPeriod(body)

**Returns the average number of days it takes for a planet to orbit the Sun.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | One of the planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, or Pluto. |

**Returns**: `float`
The mean orbital period of the body in days.

---

<a name="RefractionAngle"></a>
### RefractionAngle(refraction, altitude)

**Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.**

Given an altitude angle and a refraction option, calculates
the amount of "lift" caused by atmospheric refraction.
This is the number of degrees higher in the sky an object appears
due to lensing of the Earth's atmosphere.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Refraction`](#Refraction) | `refraction` | The option for selecting whether to correct for atmospheric lensing. If `Refraction.Normal`, a well-behaved refraction model is used. If `Refraction.Airless`, no refraction correct is performed. `Refraction.JplHorizons` is used only for compatibility testing with the JPL Horizons online tool. Any other value raises an exception. |
| `float` | `altitude` | The number of degrees above (positive) or below (negative) the horizon an object is, before being corrected for refraction. |

**Returns**: `float`
The number of additional degrees of altitude an object appears
to have, due to atmospheric refraction, depending on the
option selected by the `refraction` parameter.

---

<a name="RotateState"></a>
### RotateState(rotation, state)

**Applies a rotation to a state vector, yielding a rotated state vector.**

This function transforms a state vector in one orientation to a
state vector in another orientation. Both the position and velocity
vectors are rotated the same way.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | A rotation matrix that specifies how the orientation of the vector is to be changed. |
| [`StateVector`](#StateVector) | `state` | The state vector whose orientation is to be changed. |

**Returns**: [`StateVector`](#StateVector)
A state vector in the orientation specified by `rotation`.

---

<a name="RotateVector"></a>
### RotateVector(rotation, vector)

**Applies a rotation to a vector, yielding a rotated vector.**

This function transforms a vector in one orientation to a vector
in another orientation.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | A rotation matrix that specifies how the orientation of the vector is to be changed. |
| [`Vector`](#Vector) | `vector` | The vector whose orientation is to be changed. |

**Returns**: [`Vector`](#Vector)
A vector in the orientation specified by `rotation`.

---

<a name="RotationAxis"></a>
### RotationAxis(body, time)

**Calculates information about a body's rotation axis at a given time.**

Calculates the orientation of a body's rotation axis, along with
the rotation angle of its prime meridian, at a given moment in time.
This function uses formulas standardized by the IAU Working Group
on Cartographics and Rotational Elements 2015 report, as described
in the following document:
https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf
See [`AxisInfo`](#AxisInfo) for more detailed information.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | One of the following values: `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. |
| [`Time`](#Time) | `time` | The time at which to calculate the body's rotation axis. |

**Returns**: [`AxisInfo`](#AxisInfo)
The body's north pole direction and angle of its prime meridian.

---

<a name="Rotation_ECL_EQD"></a>
### Rotation_ECL_EQD(time)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of date.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the desired equator. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts ECL to EQD.

---

<a name="Rotation_ECL_EQJ"></a>
### Rotation_ECL_EQJ()

**Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQJ = equatorial system, using equator at J2000 epoch.

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts ECL to EQJ.

---

<a name="Rotation_ECL_HOR"></a>
### Rotation_ECL_HOR(time, observer)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: HOR = horizontal system.
Use [`HorizonFromVector`](#HorizonFromVector) to convert the return value
to a traditional altitude/azimuth pair.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the desired horizontal orientation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts ECL to HOR at `time` and for `observer`.
The components of the horizontal vector are:
x = north, y = west, z = zenith (straight up from the observer).
These components are chosen so that the "right-hand rule" works for the vector
and so that north represents the direction where azimuth = 0.

---

<a name="Rotation_EQD_ECL"></a>
### Rotation_EQD_ECL(time)

**Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of date.
Target: ECL = ecliptic system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the source equator. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQD to ECL.

---

<a name="Rotation_EQD_EQJ"></a>
### Rotation_EQD_EQJ(time)

**Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: EQJ = equatorial system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator defines the source orientation. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQD at `time` to EQJ.

---

<a name="Rotation_EQD_HOR"></a>
### Rotation_EQD_HOR(time, observer)

**Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: HOR = horizontal system.
Use [`HorizonFromVector`](#HorizonFromVector) to convert the return value
to a traditional altitude/azimuth pair.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator applies. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's location. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQD to HOR at `time` and for `observer`.
The components of the horizontal vector are:
x = north, y = west, z = zenith (straight up from the observer).
These components are chosen so that the "right-hand rule" works for the vector
and so that north represents the direction where azimuth = 0.

---

<a name="Rotation_EQJ_ECL"></a>
### Rotation_EQJ_ECL()

**Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: ECL = ecliptic system, using equator at J2000 epoch.

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to ECL.

---

<a name="Rotation_EQJ_EQD"></a>
### Rotation_EQJ_EQD(time)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator defines the target orientation. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to EQD at `time`.

---

<a name="Rotation_EQJ_GAL"></a>
### Rotation_EQJ_GAL()

**Calculates a rotation matrix from equatorial J2000 (EQJ) to galactic (GAL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: GAL = galactic system (IAU 1958 definition).

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to GAL.

---

<a name="Rotation_EQJ_HOR"></a>
### Rotation_EQJ_HOR(time, observer)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: HOR = horizontal system.
Use [`HorizonFromVector`](#HorizonFromVector) to convert the return value to
a traditional altitude/azimuth pair.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the desired horizontal orientation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
The components of the horizontal vector are:
x = north, y = west, z = zenith (straight up from the observer).
These components are chosen so that the "right-hand rule" works for the vector
and so that north represents the direction where azimuth = 0.

---

<a name="Rotation_GAL_EQJ"></a>
### Rotation_GAL_EQJ()

**Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: GAL = galactic system (IAU 1958 definition).
Target: EQJ = equatorial system, using the equator at the J2000 epoch.

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts GAL to EQJ.

---

<a name="Rotation_HOR_ECL"></a>
### Rotation_HOR_ECL(time, observer)

**Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system.
Target: ECL = ecliptic system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the horizontal observation. |
| [`Observer`](#Observer) | `observer` | The location of the horizontal observer. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts HOR to ECL.

---

<a name="Rotation_HOR_EQD"></a>
### Rotation_HOR_EQD(time, observer)

**Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQD = equatorial system, using equator of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator applies. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts HOR to EQD at `time` and for `observer`.

---

<a name="Rotation_HOR_EQJ"></a>
### Rotation_HOR_EQJ(time, observer)

**Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQJ = equatorial system, using equator at the J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that define's the observer's horizon. |

**Returns**: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts HOR to EQJ at `time` and for `observer`.

---

<a name="Search"></a>
### Search(func, context, t1, t2, dt_tolerance_seconds)

**Searches for a time at which a function's value increases through zero.**

Certain astronomy calculations involve finding a time when an event occurs.
Often such events can be defined as the root of a function:
the time at which the function's value becomes zero.
`Search` finds the *ascending root* of a function: the time at which
the function's value becomes zero while having a positive slope. That is, as time increases,
the function transitions from a negative value, through zero at a specific moment,
to a positive value later. The goal of the search is to find that specific moment.
The search function is specified by two parameters: `func` and `context`.
The `func` parameter is a function itself that accepts a time
and a context containing any other arguments needed to evaluate the function.
The `context` parameter supplies that context for the given search.
As an example, a caller may wish to find the moment a celestial body reaches a certain
ecliptic longitude. In that case, the caller might create a type (class, tuple, whatever)
that contains a [`Body`](#Body) member to specify the body and a numeric value to hold the target longitude.
A different function might use a completely different context type.
Every time it is called, `func` returns a `float` value or it raises an exception.
If `func` raises an exception, the search immediately fails and the exception is
propagated back to the caller.
Otherwise, the search proceeds until it either finds the ascending root or fails for some reason.
The search calls `func` repeatedly to rapidly narrow in on any ascending
root within the time window specified by `t1` and `t2`. The search never
reports a solution outside this time window.
`Search` uses a combination of bisection and quadratic interpolation
to minimize the number of function calls. However, it is critical that the
supplied time window be small enough that there cannot be more than one root
(ascedning or descending) within it; otherwise the search can fail.
Beyond that, it helps to make the time window as small as possible, ideally
such that the function itself resembles a smooth parabolic curve within that window.
If an ascending root is not found, or more than one root
(ascending and/or descending) exists within the window `t1`..`t2`,
`Search` will return `None` to indicate a normal search failure.
If the search does not converge within 20 iterations, it will raise
an [`Error`](#Error) exception.

| Type | Parameter | Description |
| --- | --- | --- |
| `function(context, Time)` | `func` | A function that takes an arbitrary context parameter and a [`Time`](#Time) parameter. Returns a float value.  See remarks above for more details. |
| `object` | `context` | An arbitrary data structure needed to be passed to the function `func` every time it is called. |
| `float` | `t1` | The lower time bound of the search window. See remarks above for more details. |
| `float` | `t2` | The upper time bound of the search window. See remarks above for more details. |
| `float` | `dt_tolerance_seconds` | Specifies an amount of time in seconds within which a bounded ascending root is considered accurate enough to stop. A typical value is 1 second. |

**Returns**: [`Time`](#Time) or `None`
If the search is successful, returns a #Time object that is within
`dt_tolerance_seconds` of an ascending root.
In this case, the returned time value will always be within the
inclusive range [`t1`, `t2`].
If there is no ascending root, or there is more than one ascending root,
the function returns `None`.

---

<a name="SearchAltitude"></a>
### SearchAltitude(body, observer, direction, startTime, limitDays, altitude)

**Finds the next time a body reaches a given altitude.**

Finds when the given body ascends or descends through a given
altitude angle, as seen by an observer at the specified location on the Earth.
By using the appropriate combination of `direction` and `altitude` parameters,
this function can be used to find when civil, nautical, or astronomical twilight
begins (dawn) or ends (dusk).
Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
the horizon. To find civil dawn, pass `Direction.Rise` for `direction` and -6 for `altitude`.
Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
To find civil dusk, pass `Direction.Set` for `direction` and -6 for `altitude`.
Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.
Astronomical twilight uses -18 degrees as the `altitude` value.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, or any planet other than the Earth. |
| [`Observer`](#Observer) | `observer` | The location where observation takes place. |
| [`Direction`](#Direction) | `direction` | Either `Direction.Rise` to find an ascending altitude event or `Direction.Set` to find a descending altitude event. |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |
| `float` | `limitDays` | Limits how many days to search for the body reaching the altitude angle, and defines the direction in time to search. When `limitDays` is positive, the search is performed into the future, after `startTime`. When negative, the search is performed into the past, before `startTime`. To limit the search to the same day, you can use a value of 1 day. In cases where you want to find the altitude event no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |
| `float` | `altitude` | The desired altitude angle of the body's center above (positive) or below (negative) the observer's local horizon, expressed in degrees. Must be in the range [-90, +90]. |

**Returns**: [`Time`](#Time) or `None`
If the altitude event time is found within the specified time window,
this function returns that time. Otherwise, it returns `None`.

---

<a name="SearchGlobalSolarEclipse"></a>
### SearchGlobalSolarEclipse(startTime)

**Searches for a solar eclipse visible anywhere on the Earth's surface.**

This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo) for more information.
To find a series of solar eclipses, call this function once,
then keep calling [`NextGlobalSolarEclipse`](#NextGlobalSolarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time for starting the search for a solar eclipse. |

**Returns**: [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo)

---

<a name="SearchHourAngle"></a>
### SearchHourAngle(body, observer, hourAngle, startTime, direction=1)

**Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.**

The *hour angle* of a celestial body indicates its position in the sky with respect
to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
the number of hours that have passed since the most recent time that the body has culminated,
or reached its highest point.
This function searches for the next time a celestial body reaches the given hour angle
after the date and time specified by `startTime`.
To find when a body culminates, pass 0 for `hourAngle`.
To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.
Note that, especially close to the Earth's poles, a body as seen on a given day
may always be above the horizon or always below the horizon, so the caller cannot
assume that a culminating object is visible nor that an object is below the horizon
at its minimum altitude.
On success, the function reports the date and time, along with the horizontal coordinates
of the body at that time, as seen by the given observer.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body, which can the Sun, the Moon, or any planet other than the Earth. |
| [`Observer`](#Observer) | `observer` | Indicates a location on or near the surface of the Earth where the observer is located. |
| `float` | `hourAngle` | An hour angle value in the range [0.0, 24.0) indicating the number of sidereal hours after the body's most recent culmination. |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |
| `int` | `direction` | The direction in time to perform the search: a positive value searches forward in time, a negative value searches backward in time. The function throws an exception if `direction` is zero. |

**Returns**: [`HourAngleEvent`](#HourAngleEvent)

---

<a name="SearchLocalSolarEclipse"></a>
### SearchLocalSolarEclipse(startTime, observer)

Searches for a solar eclipse visible at a specific location on the Earth's surface.
This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo) for more information.
To find a series of solar eclipses, call this function once,
then keep calling [`NextLocalSolarEclipse`](#NextLocalSolarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.
IMPORTANT: An eclipse reported by this function might be partly or
completely invisible to the observer due to the time of day.
See [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo) for more information about this topic.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time for starting the search for a solar eclipse. |
| [`Observer`](#Observer) | `observer` | The geographic location of the observer. |

**Returns**: [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo)

---

<a name="SearchLunarApsis"></a>
### SearchLunarApsis(startTime)

**Finds the time of the first lunar apogee or perigee after the given time.**

Given a date and time to start the search in `startTime`, this function finds
the next date and time that the center of the Moon reaches the closest or
farthest point in its orbit with respect to the center of the Earth, whichever
comes first after `startTime`.  The return value (of type [`Apsis`](#Apsis)) also
contains an indicator of whether the event is apogee or perigee.
The closest point is called *perigee* and the farthest point is called *apogee*.
The word *apsis* refers to either event.
To iterate through consecutive alternating perigee and apogee events,
call [`SearchLunarApsis`](#SearchLunarApsis) once, then use the return value to call [`NextLunarApsis`](#NextLunarApsis).
After that, keep feeding the previous return value from `NextLunarApsis` into
another call of `NextLunarApsis` as many times as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time at which to start searching for the next perigee or apogee. |

**Returns**: [`Apsis`](#Apsis)

---

<a name="SearchLunarEclipse"></a>
### SearchLunarEclipse(startTime)

**Searches for a lunar eclipse.**

This function finds the first lunar eclipse that occurs after `startTime`.
A lunar eclipse may be penumbral, partial, or total.
See [`LunarEclipseInfo`](#LunarEclipseInfo) for more information.
To find a series of lunar eclipses, call this function once,
then keep calling [`NextLunarEclipse`](#NextLunarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time for starting the search for a lunar eclipse. |

**Returns**: [`LunarEclipseInfo`](#LunarEclipseInfo)

---

<a name="SearchMaxElongation"></a>
### SearchMaxElongation(body, startTime)

**Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.**

Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
a telescope without atmospheric interference, are when these planets reach maximum elongation.
These are events where the planets reach the maximum angle from the Sun as seen from the Earth.
This function solves for those times, reporting the next maximum elongation event's date and time,
the elongation value itself, the relative longitude with the Sun, and whether the planet is best
observed in the morning or evening. See [`ElongationEvent`](#ElongationEvent) for more details about the returned object.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception. To find the best viewing opportunities for planets farther from the Sun than the Earth is (Mars through Pluto), use [`SearchRelativeLongitude`](#SearchRelativeLongitude) to find the next opposition event. |
| [`Time`](#Time) | `startTime` | The date and time at which to begin the search. The maximum elongation event found will always be the first one that occurs after this date and time. |

**Returns**: [`ElongationEvent`](#ElongationEvent)

---

<a name="SearchMoonNode"></a>
### SearchMoonNode(startTime)

**Searches for a time when the Moon's center crosses through the ecliptic plane.**

Searches for the first ascending or descending node of the Moon after `startTime`.
An ascending node is when the Moon's center passes through the ecliptic plane
(the plane of the Earth's orbit around the Sun) from south to north.
A descending node is when the Moon's center passes through the ecliptic plane
from north to south. Nodes indicate possible times of solar or lunar eclipses,
if the Moon also happens to be in the correct phase (new or full, respectively).
Call `SearchMoonNode` to find the first of a series of nodes.
Then call [`NextMoonNode`](#NextMoonNode) to find as many more consecutive nodes as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time for starting the search for an ascending or descending node of the Moon. |

**Returns**: [`NodeEventInfo`](#NodeEventInfo)

---

<a name="SearchMoonPhase"></a>
### SearchMoonPhase(targetLon, startTime, limitDays)

**Searches for the time that the Moon reaches a specified phase.**

Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
longitude with respect to the Sun's geocentric ecliptic longitude.
When the Moon and the Sun have the same longitude, that is defined as a new moon.
When their longitudes are 180 degrees apart, that is defined as a full moon.
This function searches for any value of the lunar phase expressed as an
angle in degrees in the range [0, 360).
If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
it is much easier to call the functions [`SearchMoonQuarter`](#SearchMoonQuarter) and [`NextMoonQuarter`](#NextMoonQuarter).
This function is useful for finding general phase angles outside those four quarters.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `targetLon` | The difference in geocentric longitude between the Sun and Moon that specifies the lunar phase being sought. This can be any value in the range [0, 360).  Certain values have conventional names: 0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter. |
| [`Time`](#Time) | `startTime` | The beginning of the time window in which to search for the Moon reaching the specified phase. |
| `float` | `limitDays` | The number of days away from `startTime` that limits the time window for the search. If the value is negative, the search is performed into the past from `startTime`. Otherwise, the search is performed into the future from `startTime`. |

**Returns**: [`Time`](#Time) or `None`

---

<a name="SearchMoonQuarter"></a>
### SearchMoonQuarter(startTime)

**Finds the first lunar quarter after the specified date and time.**

A lunar quarter is one of the following four lunar phase events:
new moon, first quarter, full moon, third quarter.
This function finds the lunar quarter that happens soonest
after the specified date and time.
To continue iterating through consecutive lunar quarters, call this function once,
followed by calls to [`NextMoonQuarter`](#NextMoonQuarter) as many times as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |

**Returns**: [`MoonQuarter`](#MoonQuarter)

---

<a name="SearchPeakMagnitude"></a>
### SearchPeakMagnitude(body, startTime)

**Searches for the date and time Venus will next appear brightest as seen from the Earth.**

This function searches for the date and time Venus appears brightest as seen from the Earth.
Currently only Venus is supported for the `body` parameter, though this could change in the future.
Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see
from the Earth, so peak magnitude events have little practical value for that planet.
Planets other than Venus and Mercury reach peak magnitude at opposition, which can
be found using [`SearchRelativeLongitude`](#SearchRelativeLongitude).
The Moon reaches peak magnitude at full moon, which can be found using
[`SearchMoonQuarter`](#SearchMoonQuarter) or [`SearchMoonPhase`](#SearchMoonPhase).
The Sun reaches peak magnitude at perihelion, which occurs each year in January.
However, the difference is minor and has little practical value.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | Currently only `Body.Venus` is allowed. Any other value results in an exception. See remarks above for more details. |
| [`Time`](#Time) | `startTime` | The date and time to start searching for the next peak magnitude event. |

**Returns**: [`IlluminationInfo`](#IlluminationInfo)

---

<a name="SearchPlanetApsis"></a>
### SearchPlanetApsis(body, startTime)

**Finds the next planet perihelion or aphelion, after a given time.**

Given a date and time to start the search in `startTime`, this function finds the
next date and time that the center of the specified planet reaches the closest or farthest point
in its orbit with respect to the center of the Sun, whichever comes first after `startTime`.
The closest point is called *perihelion* and the farthest point is called *aphelion*.
The word *apsis* refers to either event.
To iterate through consecutive alternating perihelion and aphelion events,
call `SearchPlanetApsis` once, then use the return value to call [`NextPlanetApsis`](#NextPlanetApsis).
After that, keep feeding the previous return value from `NextPlanetApsis`
into another call of `NextPlanetApsis` as many times as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find the next perihelion/aphelion event. Not allowed to be `Body.Sun` or `Body.Moon`. |
| [`Time`](#Time) | `startTime` | The date and time at which to start searching for the next perihelion or aphelion. |

**Returns**: [`Apsis`](#Apsis)

---

<a name="SearchRelativeLongitude"></a>
### SearchRelativeLongitude(body, targetRelLon, startTime)

**Searches for when the Earth and another planet are separated by a certain ecliptic longitude.**

Searches for the time when the Earth and another planet are separated by a specified angle
in ecliptic longitude, as seen from the Sun.
A relative longitude is the angle between two bodies measured in the plane of the
Earth's orbit (the ecliptic plane). The distance of the bodies above or below the ecliptic
plane is ignored. If you imagine the shadow of the body cast onto the ecliptic plane,
and the angle measured around that plane from one body to the other in the direction
the planets orbit the Sun, you will get an angle somewhere between 0 and 360 degrees.
This is the relative longitude.
Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
this function searches for the next time that the relative longitude measured from the
planet to the Earth is `targetRelLon`.
Certain astronomical events are defined in terms of relative longitude between
the Earth and another planet:
- When the relative longitude is 0 degrees, it means both planets are in the same
  direction from the Sun. For planets that orbit closer to the Sun (Mercury and Venus),
  this is known as *inferior conjunction*, a time when the other planet becomes very
  difficult to see because of being lost in the Sun's glare.
  (The only exception is in the rare event of a transit, when we see the silhouette
  of the planet passing between the Earth and the Sun.)
- When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
  this is known as *opposition*. Opposition is when the planet is closest to the Earth,
  and also when it is visible for most of the night, so it is considered the best time
  to observe the planet.
- When the relative longitude is 180 degrees, it means the other planet is on the opposite
  side of the Sun from the Earth.  This is called *superior conjunction*.  Like inferior
  conjunction, the planet is very difficult to see from the Earth.
  Superior conjunction is possible for any planet other than the Earth.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A planet other than the Earth. If `body` is not a planet, or if it is `Body.Earth`, an error occurs. |
| `float` | `targetRelLon` | The desired relative longitude, expressed in degrees. Must be in the range [0, 360). |
| [`Time`](#Time) | `startTime` | The date and time at which to begin the search. |

**Returns**: [`Time`](#Time)
The date and time of the relative longitude event.

---

<a name="SearchRiseSet"></a>
### SearchRiseSet(body, observer, direction, startTime, limitDays)

**Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.**

This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
Rise time is when the body first starts to be visible above the horizon.
For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
Set time is the moment when the body appears to vanish below the horizon.
This function corrects for typical atmospheric refraction, which causes celestial
bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).
Note that rise or set may not occur in every 24 hour period.
For example, near the Earth's poles, there are long periods of time where
the Sun stays below the horizon, never rising.
Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
significant amount during each rotation of the Earth.
Therefore callers must not assume that the function will always succeed.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, or any planet other than the Earth. |
| [`Observer`](#Observer) | `observer` | The location where observation takes place. |
| [`Direction`](#Direction) | `direction` | Either `Direction.Rise` to find a rise time or `Direction.Set` to find a set time. |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |
| `float` | `limitDays` | Limits how many days to search for a rise or set time, and defines the direction in time to search. When `limitDays` is positive, the search is performed into the future, after `startTime`. When negative, the search is performed into the past, before `startTime`. To limit a rise or set time to the same day, you can use a value of 1 day. In cases where you want to find the next rise or set time no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |

**Returns**: [`Time`](#Time) or `None`
If the rise or set time is found within the specified time window,
this function returns that time. Otherwise, it returns `None`.

---

<a name="SearchSunLongitude"></a>
### SearchSunLongitude(targetLon, startTime, limitDays)

**Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.**

This function finds the moment in time, if any exists in the given time window,
that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.
This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call [`Seasons`](#Seasons)
to calculate all equinoxes and solstices for a given calendar year.
The function searches the window of time specified by `startTime` and `startTime+limitDays`.
The search will return `None` if the Sun never reaches the longitude `targetLon` or
if the window is so large that the longitude ranges more than 180 degrees within it.
It is recommended to keep the window smaller than 10 days when possible.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `targetLon` | The desired ecliptic longitude in degrees, relative to the true equinox of date. This may be any value in the range [0, 360), although certain values have conventional meanings: 0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice. |
| [`Time`](#Time) | `startTime` | The date and time for starting the search for the desired longitude event. |
| `float` | `limitDays` | The real-valued number of days, which when added to `startTime`, limits the range of time over which the search looks. It is recommended to keep this value between 1 and 10 days. See remarks above for more details. |

**Returns**: [`Time`](#Time) or `None`

---

<a name="SearchTransit"></a>
### SearchTransit(body, startTime)

**Searches for the first transit of Mercury or Venus after a given date.**

Finds the first transit of Mercury or Venus after a specified date.
A transit is when an inferior planet passes between the Sun and the Earth
so that the silhouette of the planet is visible against the Sun in the background.
To continue the search, pass the `finish` time in the returned structure to
[`NextTransit`](#NextTransit).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`. |
| [`Time`](#Time) | `startTime` | The date and time for starting the search for a transit. |

**Returns**: [`TransitInfo`](#TransitInfo)

---

<a name="Seasons"></a>
### Seasons(year)

**Finds both equinoxes and both solstices for a given calendar year.**

The changes of seasons are defined by solstices and equinoxes.
Given a calendar year number, this function calculates the
March and September equinoxes and the June and December solstices.
The equinoxes are the moments twice each year when the plane of the
Earth's equator passes through the center of the Sun. In other words,
the Sun's declination is zero at both equinoxes.
The March equinox defines the beginning of spring in the northern hemisphere
and the beginning of autumn in the southern hemisphere.
The September equinox defines the beginning of autumn in the northern hemisphere
and the beginning of spring in the southern hemisphere.
The solstices are the moments twice each year when one of the Earth's poles
is most tilted toward the Sun. More precisely, the Sun's declination reaches
its minimum value at the December solstice, which defines the beginning of
winter in the northern hemisphere and the beginning of summer in the southern
hemisphere. The Sun's declination reaches its maximum value at the June solstice,
which defines the beginning of summer in the northern hemisphere and the beginning
of winter in the southern hemisphere.

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` | The calendar year number for which to calculate equinoxes and solstices. The value may be any integer, but only the years 1800 through 2100 have been validated for accuracy: unit testing against data from the United States Naval Observatory confirms that all equinoxes and solstices for that range of years are within 2 minutes of the correct time. |

**Returns**: [`SeasonInfo`](#SeasonInfo)

---

<a name="SiderealTime"></a>
### SiderealTime(time)

**Calculates Greenwich Apparent Sidereal Time (GAST).**

Given a date and time, this function calculates the rotation of the
Earth, represented by the equatorial angle of the Greenwich prime meridian
with respect to distant stars (not the Sun, which moves relative to background
stars by almost one degree per day).
This angle is called Greenwich Apparent Sidereal Time (GAST).
GAST is measured in sidereal hours in the half-open range [0, 24).
When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
corrected at that time for precession and nutation of the Earth's axis.
In this context, the "equinox" is the direction in space where the Earth's
orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
at the location on the Earth's orbit of the (seasonal) March equinox.
As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
then starts over at 0.
To convert to degrees, multiply the return value by 15.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to find GAST. As an optimization, this function caches the sidereal time value in `time`, unless it has already been cached, in which case the cached value is reused. |

**Returns**: `float`
GAST expressed in sidereal hours.

---

<a name="SphereFromVector"></a>
### SphereFromVector(vector)

**Converts Cartesian coordinates to spherical coordinates.**

Given a Cartesian vector, returns latitude, longitude, and distance.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vector` | Cartesian vector to be converted to spherical coordinates. |

**Returns**: [`Spherical`](#Spherical)
Spherical coordinates that are equivalent to the given vector.

---

<a name="SunPosition"></a>
### SunPosition(time)

**Calculates geocentric ecliptic coordinates for the Sun.**

This function calculates the position of the Sun as seen from the Earth.
The returned value includes both Cartesian and spherical coordinates.
The x-coordinate and longitude values in the returned object are based
on the *true equinox of date*: one of two points in the sky where the instantaneous
plane of the Earth's equator at the given date and time (the *equatorial plane*)
intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
By convention, the apparent location of the Sun at the March equinox is chosen
as the longitude origin and x-axis direction, instead of the one for September.
`SunPosition` corrects for precession and nutation of the Earth's axis
in order to obtain the exact equatorial plane at the given time.
This function can be used for calculating changes of seasons: equinoxes and solstices.
In fact, the function [`Seasons`](#Seasons) does use this function for that purpose.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Sun's position. |

**Returns**: [`EclipticCoordinates`](#EclipticCoordinates)
The ecliptic coordinates of the Sun using the Earth's true equator of date.

---

<a name="VectorFromHorizon"></a>
### VectorFromHorizon(sphere, time, refraction)

**Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Spherical`](#Spherical) | `sphere` | A structure that contains apparent horizontal coordinates: `lat` holds the refracted azimuth angle, `lon` holds the azimuth in degrees clockwise from north, and `dist` holds the distance from the observer to the object in AU. |
| [`Time`](#Time) | `time` | The date and time of the observation. This is needed because the returned vector object requires a valid time value when passed to certain other functions. |
| [`Refraction`](#Refraction) | `refraction` | See remarks in function [`RefractionAngle`](#RefractionAngle). |

**Returns**: [`Vector`](#Vector)
A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).

---

<a name="VectorFromSphere"></a>
### VectorFromSphere(sphere, time)

**Converts spherical coordinates to Cartesian coordinates.**

Given spherical coordinates and a time at which they are valid,
returns a vector of Cartesian coordinates. The returned value
includes the time, as required by all `Time` objects.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Spherical`](#Spherical) | `sphere` | Spherical coordinates to be converted. |
| [`Time`](#Time) | `time` | The time that should be included in the returned vector. |

**Returns**: [`Vector`](#Vector)
The vector form of the supplied spherical coordinates.

---

<a name="VectorObserver"></a>
### VectorObserver(vector, ofdate)

**Calculates the geographic location corresponding to an equatorial vector.**

This is the inverse function of [`ObserverVector`](#ObserverVector).
Given a geocentric equatorial vector, it returns the geographic
latitude, longitude, and elevation for that vector.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vector` | The geocentric equatorial position vector for which to find geographic coordinates. The components are expressed in astronomical units (AU). The time `vector.t` determines the Earth's rotation. |
| `bool` | `ofdate` | Selects the date of the Earth's equator in which `vector` is expressed. The caller may pass `False` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the the time `vector.t`. Or the caller may pass `True` to use the Earth's equator at `vector.t` as the orientation. |

**Returns**: [`Observer`](#Observer)
The geographic latitude, longitude, and elevation above sea level
that corresponds to the given equatorial vector.

