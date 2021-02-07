## Classes

<dl>
<dt><a href="#AstroTime">AstroTime</a></dt>
<dd></dd>
<dt><a href="#Vector">Vector</a></dt>
<dd><p>Holds the Cartesian coordinates of a vector in 3D space,
along with the time at which the vector is valid.</p>
</dd>
<dt><a href="#Spherical">Spherical</a></dt>
<dd><p>Holds spherical coordinates: latitude, longitude, distance.</p>
</dd>
<dt><a href="#EquatorialCoordinates">EquatorialCoordinates</a></dt>
<dd><p>Holds right ascension, declination, and distance of a celestial object.</p>
</dd>
<dt><a href="#RotationMatrix">RotationMatrix</a></dt>
<dd><p>Contains a rotation matrix that can be used to transform one coordinate system to another.</p>
</dd>
<dt><a href="#HorizontalCoordinates">HorizontalCoordinates</a></dt>
<dd><p>Holds azimuth (compass direction) and altitude (angle above/below the horizon)
of a celestial object as seen by an observer at a particular location on the Earth&#39;s surface.
Also holds right ascension and declination of the same object.
All of these coordinates are optionally adjusted for atmospheric refraction;
therefore the right ascension and declination values may not exactly match
those found inside a corresponding <a href="#EquatorialCoordinates">EquatorialCoordinates</a> object.</p>
</dd>
<dt><a href="#EclipticCoordinates">EclipticCoordinates</a></dt>
<dd><p>Holds ecliptic coordinates of a celestial body.
The origin and date of the coordinate system may vary depending on the caller&#39;s usage.
In general, ecliptic coordinates are measured with respect to the mean plane of the Earth&#39;s
orbit around the Sun.
Includes Cartesian coordinates <code>(ex, ey, ez)</code> measured in
<a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU)
and spherical coordinates <code>(elon, elat)</code> measured in degrees.</p>
</dd>
<dt><a href="#Observer">Observer</a></dt>
<dd><p>Represents the geographic location of an observer on the surface of the Earth.</p>
</dd>
<dt><a href="#IlluminationInfo">IlluminationInfo</a></dt>
<dd><p>Contains information about the apparent brightness and sunlit phase of a celestial object.</p>
</dd>
<dt><a href="#MoonQuarter">MoonQuarter</a></dt>
<dd><p>Represents a quarter lunar phase, along with when it occurs.</p>
</dd>
<dt><a href="#HourAngleEvent">HourAngleEvent</a></dt>
<dd><p>Returns information about an occurrence of a celestial body
reaching a given hour angle as seen by an observer at a given
location on the surface of the Earth.</p>
</dd>
<dt><a href="#SeasonInfo">SeasonInfo</a></dt>
<dd><p>Represents the dates and times of the two solstices
and the two equinoxes in a given calendar year.
These four events define the changing of the seasons on the Earth.</p>
</dd>
<dt><a href="#ElongationEvent">ElongationEvent</a></dt>
<dd><p>Represents the angular separation of a body from the Sun as seen from the Earth
and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.</p>
</dd>
<dt><a href="#Apsis">Apsis</a></dt>
<dd><p>Represents a closest or farthest point in a body&#39;s orbit around its primary.
For a planet orbiting the Sun, this is a perihelion or aphelion, respectively.
For the Moon orbiting the Earth, this is a perigee or apogee, respectively.</p>
</dd>
<dt><a href="#ConstellationInfo">ConstellationInfo</a></dt>
<dd><p>Reports the constellation that a given celestial point lies within.</p>
</dd>
<dt><a href="#LunarEclipseInfo">LunarEclipseInfo</a></dt>
<dd><p>Returns information about a lunar eclipse.</p>
<p>Returned by <a href="#SearchLunarEclipse">SearchLunarEclipse</a> or <a href="#NextLunarEclipse">NextLunarEclipse</a>
to report information about a lunar eclipse event.
When a lunar eclipse is found, it is classified as penumbral, partial, or total.
Penumbral eclipses are difficult to observe, because the moon is only slightly dimmed
by the Earth&#39;s penumbra; no part of the Moon touches the Earth&#39;s umbra.
Partial eclipses occur when part, but not all, of the Moon touches the Earth&#39;s umbra.
Total eclipses occur when the entire Moon passes into the Earth&#39;s umbra.</p>
<p>The <code>kind</code> field thus holds one of the strings <code>&quot;penumbral&quot;</code>, <code>&quot;partial&quot;</code>,
or <code>&quot;total&quot;</code>, depending on the kind of lunar eclipse found.</p>
<p>Field <code>peak</code> holds the date and time of the peak of the eclipse, when it is at its peak.</p>
<p>Fields <code>sd_penum</code>, <code>sd_partial</code>, and <code>sd_total</code> hold the semi-duration of each phase
of the eclipse, which is half of the amount of time the eclipse spends in each
phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
By converting from minutes to days, and subtracting/adding with <code>peak</code>, the caller
may determine the date and time of the beginning/end of each eclipse phase.</p>
</dd>
<dt><a href="#GlobalSolarEclipseInfo">GlobalSolarEclipseInfo</a></dt>
<dd><p>Reports the time and geographic location of the peak of a solar eclipse.</p>
<pre><code>Returned by [SearchGlobalSolarEclipse](#SearchGlobalSolarEclipse) or [NextGlobalSolarEclipse](#NextGlobalSolarEclipse)
to report information about a solar eclipse event.

Field `peak` holds the date and time of the peak of the eclipse, defined as
the instant when the axis of the Moon&#39;s shadow cone passes closest to the Earth&#39;s center.

The eclipse is classified as partial, annular, or total, depending on the
maximum amount of the Sun&#39;s disc obscured, as seen at the peak location
on the surface of the Earth.

The `kind` field thus holds one of the strings `&quot;partial&quot;`, `&quot;annular&quot;`, or `&quot;total&quot;`.
A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
An annular eclipse is like a total eclipse, but the Moon is too far from the Earth&#39;s surface
to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
A partial eclipse is when the Moon blocks part of the Sun&#39;s disc, but nobody on the Earth
observes either a total or annular eclipse.

If `kind` is `&quot;total&quot;` or `&quot;annular&quot;`, the `latitude` and `longitude`
fields give the geographic coordinates of the center of the Moon&#39;s shadow projected
onto the daytime side of the Earth at the instant of the eclipse&#39;s peak.
If `kind` has any other value, `latitude` and `longitude` are undefined and should
not be used.
</code></pre>
</dd>
<dt><a href="#EclipseEvent">EclipseEvent</a></dt>
<dd></dd>
<dt><a href="#LocalSolarEclipseInfo">LocalSolarEclipseInfo</a></dt>
<dd></dd>
<dt><a href="#TransitInfo">TransitInfo</a></dt>
<dd></dd>
</dl>

## Constants

<dl>
<dt><a href="#Bodies">Bodies</a> : <code>Array.&lt;string&gt;</code></dt>
<dd><p>An array of strings, each a name of a supported astronomical body.
     Not all bodies are valid for all functions, but any string not in this
     list is not supported at all.</p>
</dd>
</dl>

## Functions

<dl>
<dt><a href="#AngleBetween">AngleBetween(a, b)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates the angle in degrees between two vectors.
The angle is measured in the plane that contains both vectors.</p>
</dd>
<dt><a href="#TerrestrialTime">TerrestrialTime(ut)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates Terrestrial Time (TT) from Universal Time (UT).</p>
</dd>
<dt><a href="#MakeTime">MakeTime(date)</a> ⇒ <code><a href="#AstroTime">AstroTime</a></code></dt>
<dd><p>Given a Date object or a number days since noon (12:00) on January 1, 2000 (UTC),
this function creates an <a href="#AstroTime">AstroTime</a> object.
Given an <a href="#AstroTime">AstroTime</a> object, returns the same object unmodified.
Use of this function is not required for any of the other exposed functions in this library,
because they all guarantee converting date/time parameters to AstroTime
as needed. However, it may be convenient for callers who need to understand
the difference between UTC and TT (Terrestrial Time). In some use cases,
converting once to AstroTime format and passing the result into multiple
function calls may be more efficient than passing in native JavaScript Date objects.</p>
</dd>
<dt><a href="#MakeSpherical">MakeSpherical(lat, lon, dist)</a> ⇒ <code><a href="#Spherical">Spherical</a></code></dt>
<dd><p>Create spherical coordinates.</p>
</dd>
<dt><a href="#MakeRotation">MakeRotation(rot)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Creates a rotation matrix that can be used to transform one coordinate system to another.</p>
</dd>
<dt><a href="#Horizon">Horizon(date, observer, ra, dec, refraction)</a> ⇒ <code><a href="#HorizontalCoordinates">HorizontalCoordinates</a></code></dt>
<dd><p>Given a date and time, a geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial body,
returns horizontal coordinates (azimuth and altitude angles) for that body
as seen by that observer. Allows optional correction for atmospheric refraction.</p>
</dd>
<dt><a href="#MakeObserver">MakeObserver(latitude_degrees, longitude_degrees, height_in_meters)</a></dt>
<dd><p>Creates an <a href="#Observer">Observer</a> object that represents a location
on the surface of the Earth from which observations are made.</p>
</dd>
<dt><a href="#SunPosition">SunPosition(date)</a> ⇒ <code><a href="#EclipticCoordinates">EclipticCoordinates</a></code></dt>
<dd><p>Returns apparent geocentric true ecliptic coordinates of date for the Sun.
<i>Geocentric</i> means coordinates as the Sun would appear to a hypothetical observer
at the center of the Earth.
<i>Ecliptic coordinates of date</i> are measured along the plane of the Earth&#39;s mean
orbit around the Sun, using the
<a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>
of the Earth as adjusted for precession and nutation of the Earth&#39;s
axis of rotation on the given date.</p>
</dd>
<dt><a href="#Equator">Equator(body, date, observer, ofdate, aberration)</a> ⇒ <code><a href="#EquatorialCoordinates">EquatorialCoordinates</a></code></dt>
<dd><p>Returns topocentric equatorial coordinates (right ascension and declination)
in one of two different systems: J2000 or true-equator-of-date.
Allows optional correction for aberration.
Always corrects for light travel time (represents the object as seen by the observer
with light traveling to the Earth at finite speed, not where the object is right now).
<i>Topocentric</i> refers to a position as seen by an observer on the surface of the Earth.
This function corrects for
<a href="https://en.wikipedia.org/wiki/Parallax">parallax</a>
of the object between a geocentric observer and a topocentric observer.
This is most significant for the Moon, because it is so close to the Earth.
However, it can have a small effect on the apparent positions of other bodies.</p>
</dd>
<dt><a href="#Ecliptic">Ecliptic(gx, gy, gz)</a> ⇒ <code><a href="#EclipticCoordinates">EclipticCoordinates</a></code></dt>
<dd><p>Given J2000 equatorial Cartesian coordinates,
returns J2000 ecliptic latitude, longitude, and cartesian coordinates.
You can call <a href="#GeoVector">GeoVector</a> and use its (x, y, z) return values
to pass into this function.</p>
</dd>
<dt><a href="#GeoMoon">GeoMoon(date)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Calculates the geocentric Cartesian coordinates for the Moon in the J2000 equatorial system.
Based on the Nautical Almanac Office&#39;s <i>Improved Lunar Ephemeris</i> of 1954,
which in turn derives from E. W. Brown&#39;s lunar theories.
Adapted from Turbo Pascal code from the book
<a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a>
by Montenbruck and Pfleger.</p>
</dd>
<dt><a href="#HelioVector">HelioVector(body, date)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Calculates heliocentric (i.e., with respect to the center of the Sun)
Cartesian coordinates in the J2000 equatorial system of a celestial
body at a specified time. The position is not corrected for light travel time or aberration.</p>
</dd>
<dt><a href="#HelioDistance">HelioDistance(body, date)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates the distance between a body and the Sun at a given time.</p>
<p>Given a date and time, this function calculates the distance between
the center of <code>body</code> and the center of the Sun.
For the planets Mercury through Neptune, this function is significantly
more efficient than calling <a href="#HelioVector">HelioVector</a> followed by taking the length
of the resulting vector.</p>
</dd>
<dt><a href="#GeoVector">GeoVector(body, date, aberration)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Calculates geocentric (i.e., with respect to the center of the Earth)
Cartesian coordinates in the J2000 equatorial system of a celestial
body at a specified time. The position is always corrected for light travel time:
this means the position of the body is &quot;back-dated&quot; based on how long it
takes light to travel from the body to an observer on the Earth.
Also, the position can optionally be corrected for aberration, an effect
causing the apparent direction of the body to be shifted based on
transverse movement of the Earth with respect to the rays of light
coming from that body.</p>
</dd>
<dt><a href="#Search">Search(func, t1, t2, options)</a> ⇒ <code>null</code> | <code><a href="#AstroTime">AstroTime</a></code></dt>
<dd><p>Search for next time <i>t</i> (such that <i>t</i> is between <code>t1</code> and <code>t2</code>)
that <code>func(t)</code> crosses from a negative value to a non-negative value.
The given function must have &quot;smooth&quot; behavior over the entire inclusive range [<code>t1</code>, <code>t2</code>],
meaning that it behaves like a continuous differentiable function.
It is not required that <code>t1</code> &lt; <code>t2</code>; <code>t1</code> &gt; <code>t2</code>
allows searching backward in time.
Note: <code>t1</code> and <code>t2</code> must be chosen such that there is no possibility
of more than one zero-crossing (ascending or descending), or it is possible
that the &quot;wrong&quot; event will be found (i.e. not the first event after t1)
or even that the function will return null, indicating that no event was found.</p>
</dd>
<dt><a href="#SearchSunLongitude">SearchSunLongitude(targetLon, dateStart, limitDays)</a> ⇒ <code><a href="#AstroTime">AstroTime</a></code> | <code>null</code></dt>
<dd><p>Searches for the moment in time when the center of the Sun reaches a given apparent
ecliptic longitude, as seen from the center of the Earth, within a given range of dates.
This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call <a href="#Seasons">Seasons</a>
to calculate equinoxes and solstices for a given calendar year.
<code>SearchSunLongitude</code> is more general in that it allows searching for arbitrary longitude values.</p>
</dd>
<dt><a href="#LongitudeFromSun">LongitudeFromSun(body, date)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates the ecliptic longitude difference
between the given body and the Sun as seen from
the Earth at a given moment in time.
The returned value ranges [0, 360) degrees.
By definition, the Earth and the Sun are both in the plane of the ecliptic.
Ignores the height of the <code>body</code> above or below the ecliptic plane;
the resulting angle is measured around the ecliptic plane for the &quot;shadow&quot;
of the body onto that plane.</p>
</dd>
<dt><a href="#AngleFromSun">AngleFromSun(body, date)</a> ⇒ <code>number</code></dt>
<dd><p>Returns the full angle seen from
the Earth, between the given body and the Sun.
Unlike <a href="#LongitudeFromSun">LongitudeFromSun</a>, this function does not
project the body&#39;s &quot;shadow&quot; onto the ecliptic;
the angle is measured in 3D space around the plane that
contains the centers of the Earth, the Sun, and <code>body</code>.</p>
</dd>
<dt><a href="#EclipticLongitude">EclipticLongitude(body, date)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates heliocentric ecliptic longitude based on the J2000 equinox.</p>
</dd>
<dt><a href="#Illumination">Illumination(body, date)</a> ⇒ <code><a href="#IlluminationInfo">IlluminationInfo</a></code></dt>
<dd><p>Calculates the phase angle, visual maginitude,
and other values relating to the body&#39;s illumination
at the given date and time, as seen from the Earth.</p>
</dd>
<dt><a href="#SearchRelativeLongitude">SearchRelativeLongitude(body, targetRelLon, startDate)</a> ⇒ <code><a href="#AstroTime">AstroTime</a></code></dt>
<dd><p>Searches for the date and time the relative ecliptic longitudes of
the specified body and the Earth, as seen from the Sun, reach a certain
difference. This function is useful for finding conjunctions and oppositions
of the planets. For the opposition of a superior planet (Mars, Jupiter, ..., Pluto),
or the inferior conjunction of an inferior planet (Mercury, Venus),
call with <code>targetRelLon</code> = 0. The 0 value indicates that both
planets are on the same ecliptic longitude line, ignoring the other planet&#39;s
distance above or below the plane of the Earth&#39;s orbit.
For superior conjunctions, call with <code>targetRelLon</code> = 180.
This means the Earth and the other planet are on opposite sides of the Sun.</p>
</dd>
<dt><a href="#MoonPhase">MoonPhase(date)</a> ⇒ <code>number</code></dt>
<dd><p>Determines the moon&#39;s phase expressed as an ecliptic longitude.</p>
</dd>
<dt><a href="#SearchMoonPhase">SearchMoonPhase(targetLon, dateStart, limitDays)</a> ⇒ <code><a href="#AstroTime">AstroTime</a></code> | <code>null</code></dt>
<dd><p>Searches for the date and time that the Moon reaches a specified phase.
Lunar phases are defined in terms of geocentric ecliptic longitudes
with respect to the Sun.  When the Moon and the Sun have the same ecliptic
longitude, that is defined as a new moon. When the two ecliptic longitudes
are 180 degrees apart, that is defined as a full moon.
To enumerate quarter lunar phases, it is simpler to call
<a href="#SearchMoonQuarter">SearchMoonQuarter</a> once, followed by repeatedly calling
<a href="#NextMoonQuarter">NextMoonQuarter</a>. <code>SearchMoonPhase</code> is only
necessary for finding other lunar phases than the usual quarter phases.</p>
</dd>
<dt><a href="#SearchMoonQuarter">SearchMoonQuarter(dateStart)</a> ⇒ <code><a href="#MoonQuarter">MoonQuarter</a></code></dt>
<dd><p>Finds the first quarter lunar phase after the specified date and time.
The quarter lunar phases are: new moon, first quarter, full moon, and third quarter.
To enumerate quarter lunar phases, call <code>SearchMoonQuarter</code> once,
then pass its return value to <a href="#NextMoonQuarter">NextMoonQuarter</a> to find the next
<code>MoonQuarter</code>. Keep calling <code>NextMoonQuarter</code> in a loop,
passing the previous return value as the argument to the next call.</p>
</dd>
<dt><a href="#NextMoonQuarter">NextMoonQuarter(mq)</a></dt>
<dd><p>Given a <a href="#MoonQuarter">MoonQuarter</a> object, finds the next consecutive
quarter lunar phase. See remarks in <a href="#SearchMoonQuarter">SearchMoonQuarter</a>
for explanation of usage.</p>
</dd>
<dt><a href="#SearchRiseSet">SearchRiseSet(body, observer, direction, dateStart, limitDays)</a> ⇒ <code><a href="#AstroTime">AstroTime</a></code> | <code>null</code></dt>
<dd><p>Finds a rise or set time for the given body as
seen by an observer at the specified location on the Earth.
Rise time is defined as the moment when the top of the body
is observed to first appear above the horizon in the east.
Set time is defined as the moment the top of the body
is observed to sink below the horizon in the west.
The times are adjusted for typical atmospheric refraction conditions.</p>
</dd>
<dt><a href="#SearchHourAngle">SearchHourAngle(body, observer, hourAngle, dateStart)</a> ⇒ <code><a href="#HourAngleEvent">HourAngleEvent</a></code></dt>
<dd><p>Finds the next time the given body is seen to reach the specified
<a href="https://en.wikipedia.org/wiki/Hour_angle">hour angle</a>
by the given observer.
Providing <code>hourAngle</code> = 0 finds the next maximum altitude event (culmination).
Providing <code>hourAngle</code> = 12 finds the next minimum altitude event.
Note that, especially close to the Earth&#39;s poles, a body as seen on a given day
may always be above the horizon or always below the horizon, so the caller cannot
assume that a culminating object is visible nor that an object is below the horizon
at its minimum altitude.</p>
</dd>
<dt><a href="#Seasons">Seasons(year)</a> ⇒ <code><a href="#SeasonInfo">SeasonInfo</a></code></dt>
<dd><p>Finds the equinoxes and solstices for a given calendar year.</p>
</dd>
<dt><a href="#Elongation">Elongation(body)</a> ⇒ <code><a href="#ElongationEvent">ElongationEvent</a></code></dt>
<dd><p>Calculates angular separation of a body from the Sun as seen from the Earth
and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.
See the return type <a href="#ElongationEvent">ElongationEvent</a> for details.</p>
<p>This function is helpful for determining how easy
it is to view a planet away from the Sun&#39;s glare on a given date.
It also determines whether the object is visible in the morning or evening;
this is more important the smaller the elongation is.
It is also used to determine how far a planet is from opposition, conjunction, or quadrature.</p>
</dd>
<dt><a href="#SearchMaxElongation">SearchMaxElongation(body, startDate)</a> ⇒ <code><a href="#ElongationEvent">ElongationEvent</a></code></dt>
<dd><p>Searches for the next maximum elongation event for Mercury or Venus
that occurs after the given start date. Calling with other values
of <code>body</code> will result in an exception.
Maximum elongation occurs when the body has the greatest
angular separation from the Sun, as seen from the Earth.
Returns an <code>ElongationEvent</code> object containing the date and time of the next
maximum elongation, the elongation in degrees, and whether
the body is visible in the morning or evening.</p>
</dd>
<dt><a href="#SearchPeakMagnitude">SearchPeakMagnitude(body, startDate)</a> ⇒ <code><a href="#IlluminationInfo">IlluminationInfo</a></code></dt>
<dd><p>Searches for the date and time Venus will next appear brightest as seen from the Earth.</p>
</dd>
<dt><a href="#SearchLunarApsis">SearchLunarApsis(startDate)</a> ⇒ <code><a href="#Apsis">Apsis</a></code></dt>
<dd><p>Finds the next perigee (closest approach) or apogee (farthest remove) of the Moon
that occurs after the specified date and time.</p>
</dd>
<dt><a href="#NextLunarApsis">NextLunarApsis(apsis)</a> ⇒ <code><a href="#Apsis">Apsis</a></code></dt>
<dd><p>Given a lunar apsis returned by an initial call to <a href="#SearchLunarApsis">SearchLunarApsis</a>,
or a previous call to <code>NextLunarApsis</code>, finds the next lunar apsis.
If the given apsis is a perigee, this function finds the next apogee, and vice versa.</p>
</dd>
<dt><a href="#SearchPlanetApsis">SearchPlanetApsis(body, startTime)</a> ⇒ <code><a href="#Apsis">Apsis</a></code></dt>
<dd><p>Finds the date and time of a planet&#39;s perihelion (closest approach to the Sun)
or aphelion (farthest distance from the Sun) after a given time.</p>
<p>Given a date and time to start the search in <code>startTime</code>, this function finds the
next date and time that the center of the specified planet reaches the closest or farthest point
in its orbit with respect to the center of the Sun, whichever comes first
after <code>startTime</code>.</p>
<p>The closest point is called <em>perihelion</em> and the farthest point is called <em>aphelion</em>.
The word <em>apsis</em> refers to either event.</p>
<p>To iterate through consecutive alternating perihelion and aphelion events,
call <code>SearchPlanetApsis</code> once, then use the return value to call
<a href="#NextPlanetApsis">NextPlanetApsis</a>. After that, keep feeding the previous return value
from <code>NextPlanetApsis</code> into another call of <code>NextPlanetApsis</code>
as many times as desired.</p>
</dd>
<dt><a href="#NextPlanetApsis">NextPlanetApsis(body, apsis)</a> ⇒ <code><a href="#Apsis">Apsis</a></code></dt>
<dd><p>Finds the next planetary perihelion or aphelion event in a series.</p>
<p>This function requires an <a href="#Apsis">Apsis</a> value obtained from a call
to <a href="#SearchPlanetApsis">SearchPlanetApsis</a> or <code>NextPlanetApsis</code>.
Given an aphelion event, this function finds the next perihelion event, and vice versa.
See <a href="#SearchPlanetApsis">SearchPlanetApsis</a> for more details.</p>
</dd>
<dt><a href="#InverseRotation">InverseRotation(rotation)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates the inverse of a rotation matrix.
Given a rotation matrix that performs some coordinate transform,
this function returns the matrix that reverses that trasnform.</p>
</dd>
<dt><a href="#CombineRotation">CombineRotation(a, b)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Creates a rotation based on applying one rotation followed by another.
Given two rotation matrices, returns a combined rotation matrix that is
equivalent to rotating based on the first matrix, followed by the second.</p>
</dd>
<dt><a href="#VectorFromSphere">VectorFromSphere(sphere, time)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Converts spherical coordinates to Cartesian coordinates.
Given spherical coordinates and a time at which they are valid,
returns a vector of Cartesian coordinates. The returned value
includes the time, as required by <code>AstroTime</code>.</p>
</dd>
<dt><a href="#VectorFromEquator">VectorFromEquator(equ, time)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Given angular equatorial coordinates in <code>equ</code>, calculates equatorial vector.</p>
</dd>
<dt><a href="#EquatorFromVector">EquatorFromVector(vec)</a> ⇒ <code><a href="#EquatorialCoordinates">EquatorialCoordinates</a></code></dt>
<dd><p>Given an equatorial vector, calculates equatorial angular coordinates.</p>
</dd>
<dt><a href="#SphereFromVector">SphereFromVector(vector)</a> ⇒ <code><a href="#Spherical">Spherical</a></code></dt>
<dd><p>Converts Cartesian coordinates to spherical coordinates.</p>
<p>Given a Cartesian vector, returns latitude, longitude, and distance.</p>
</dd>
<dt><a href="#HorizonFromVector">HorizonFromVector(vector, refraction)</a> ⇒ <code><a href="#Spherical">Spherical</a></code></dt>
<dd><p>Converts Cartesian coordinates to horizontal coordinates.</p>
<p>Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.</p>
<p><em>IMPORTANT:</em> This function differs from <a href="#SphereFromVector">SphereFromVector</a> in two ways:</p>
<ul>
<li><code>SphereFromVector</code> returns a <code>lon</code> value that represents azimuth defined counterclockwise
from north (e.g., west = +90), but this function represents a clockwise rotation
(e.g., east = +90). The difference is because <code>SphereFromVector</code> is intended
to preserve the vector &quot;right-hand rule&quot;, while this function defines azimuth in a more
traditional way as used in navigation and cartography.</li>
<li>This function optionally corrects for atmospheric refraction, while <code>SphereFromVector</code> does not.</li>
</ul>
<p>The returned object contains the azimuth in <code>lon</code>.
It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.</p>
<p>The altitude is stored in <code>lat</code>.</p>
<p>The distance to the observed object is stored in <code>dist</code>,
and is expressed in astronomical units (AU).</p>
</dd>
<dt><a href="#VectorFromHorizon">VectorFromHorizon(sphere, time, refraction)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Given apparent angular horizontal coordinates in <code>sphere</code>, calculate horizontal vector.</p>
</dd>
<dt><a href="#Refraction">Refraction(refraction, altitude)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates the amount of &quot;lift&quot; to an altitude angle caused by atmospheric refraction.</p>
<p>Given an altitude angle and a refraction option, calculates
the amount of &quot;lift&quot; caused by atmospheric refraction.
This is the number of degrees higher in the sky an object appears
due to the lensing of the Earth&#39;s atmosphere.</p>
</dd>
<dt><a href="#InverseRefraction">InverseRefraction(refraction, bent_altitude)</a> ⇒ <code>number</code></dt>
<dd><p>Calculates the inverse of an atmospheric refraction angle.</p>
<p>Given an observed altitude angle that includes atmospheric refraction,
calculate the negative angular correction to obtain the unrefracted
altitude. This is useful for cases where observed horizontal
coordinates are to be converted to another orientation system,
but refraction first must be removed from the observed position.</p>
</dd>
<dt><a href="#RotateVector">RotateVector(rotation, vector)</a> ⇒ <code><a href="#Vector">Vector</a></code></dt>
<dd><p>Applies a rotation to a vector, yielding a rotated vector.</p>
<p>This function transforms a vector in one orientation to a vector
in another orientation.</p>
</dd>
<dt><a href="#Rotation_EQJ_ECL">Rotation_EQJ_ECL()</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: ECL = ecliptic system, using equator at J2000 epoch.</p>
</dd>
<dt><a href="#Rotation_ECL_EQJ">Rotation_ECL_EQJ()</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQJ = equatorial system, using equator at J2000 epoch.</p>
</dd>
<dt><a href="#Rotation_EQJ_EQD">Rotation_EQJ_EQD(time)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of the specified date/time.</p>
</dd>
<dt><a href="#Rotation_EQD_EQJ">Rotation_EQD_EQJ(time)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: EQJ = equatorial system, using equator at J2000 epoch.</p>
</dd>
<dt><a href="#Rotation_EQD_HOR">Rotation_EQD_HOR(time, observer)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: HOR = horizontal system.</p>
<p>Use <code>HorizonFromVector</code> to convert the return value
to a traditional altitude/azimuth pair.</p>
</dd>
<dt><a href="#Rotation_HOR_EQD">Rotation_HOR_EQD(time, observer)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQD = equatorial system, using equator of the specified date/time.</p>
</dd>
<dt><a href="#Rotation_HOR_EQJ">Rotation_HOR_EQJ(time, observer)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQJ = equatorial system, using equator at the J2000 epoch.</p>
</dd>
<dt><a href="#Rotation_EQJ_HOR">Rotation_EQJ_HOR(time, observer)</a> ⇒</dt>
<dd><p>Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: HOR = horizontal system.</p>
<p>Use <a href="#HorizonFromVector">HorizonFromVector</a> to convert the return value
to a traditional altitude/azimuth pair.</p>
</dd>
<dt><a href="#Rotation_EQD_ECL">Rotation_EQD_ECL(time)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of date.
Target: ECL = ecliptic system, using equator at J2000 epoch.</p>
</dd>
<dt><a href="#Rotation_ECL_EQD">Rotation_ECL_EQD(time)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of date.</p>
</dd>
<dt><a href="#Rotation_ECL_HOR">Rotation_ECL_HOR(time, observer)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: HOR = horizontal system.</p>
<p>Use <a href="#HorizonFromVector">HorizonFromVector</a> to convert the return value
to a traditional altitude/azimuth pair.</p>
</dd>
<dt><a href="#Rotation_HOR_ECL">Rotation_HOR_ECL(time, observer)</a> ⇒ <code><a href="#RotationMatrix">RotationMatrix</a></code></dt>
<dd><p>Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).</p>
<p>This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system.
Target: ECL = ecliptic system, using equator at J2000 epoch.</p>
</dd>
<dt><a href="#Constellation">Constellation(ra, dec)</a> ⇒ <code><a href="#ConstellationInfo">ConstellationInfo</a></code></dt>
<dd><p>Determines the constellation that contains the given point in the sky.</p>
<p>Given J2000 equatorial (EQJ) coordinates of a point in the sky,
determines the constellation that contains that point.</p>
</dd>
<dt><a href="#SearchLunarEclipse">SearchLunarEclipse(date)</a> ⇒ <code><a href="#LunarEclipseInfo">LunarEclipseInfo</a></code></dt>
<dd></dd>
<dt><a href="#NextLunarEclipse">NextLunarEclipse(prevEclipseTime)</a> ⇒ <code><a href="#LunarEclipseInfo">LunarEclipseInfo</a></code></dt>
<dd></dd>
<dt><a href="#SearchGlobalSolarEclipse">SearchGlobalSolarEclipse(startTime)</a> ⇒ <code><a href="#GlobalSolarEclipseInfo">GlobalSolarEclipseInfo</a></code></dt>
<dd></dd>
<dt><a href="#NextGlobalSolarEclipse">NextGlobalSolarEclipse(prevEclipseTime)</a> ⇒ <code><a href="#GlobalSolarEclipseInfo">GlobalSolarEclipseInfo</a></code></dt>
<dd></dd>
<dt><a href="#SearchLocalSolarEclipse">SearchLocalSolarEclipse(startTime, observer)</a> ⇒ <code><a href="#LocalSolarEclipseInfo">LocalSolarEclipseInfo</a></code></dt>
<dd></dd>
<dt><a href="#NextLocalSolarEclipse">NextLocalSolarEclipse(prevEclipseTime, observer)</a> ⇒ <code><a href="#LocalSolarEclipseInfo">LocalSolarEclipseInfo</a></code></dt>
<dd></dd>
<dt><a href="#SearchTransit">SearchTransit(body, startTime)</a> ⇒ <code><a href="#TransitInfo">TransitInfo</a></code></dt>
<dd></dd>
<dt><a href="#NextTransit">NextTransit(body, prevTransitTime)</a> ⇒ <code><a href="#TransitInfo">TransitInfo</a></code></dt>
<dd></dd>
</dl>

## Typedefs

<dl>
<dt><a href="#SearchOptions">SearchOptions</a> : <code>Object</code></dt>
<dd><p>Options for the <a href="#Search">Search</a> function.</p>
</dd>
</dl>

<a name="AstroTime"></a>

## AstroTime
**Kind**: global class  
**Brief**: The date and time of an astronomical observation.

Objects of this type are used throughout the internals
of the Astronomy library, and are included in certain return objects.
The constructor is not accessible outside the Astronomy library;
outside users should call the [MakeTime](#MakeTime) function
to create an `AstroTime` object.  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> | The JavaScript Date object for the given date and time.      This Date corresponds to the numeric day value stored in the `ut` property. |
| ut | <code>number</code> | Universal Time (UT1/UTC) in fractional days since the J2000 epoch.      Universal Time represents time measured with respect to the Earth's rotation,      tracking mean solar days.      The Astronomy library approximates UT1 and UTC as being the same thing.      This gives sufficient accuracy for the precision requirements of this project. |
| tt | <code>number</code> | Terrestrial Time in fractional days since the J2000 epoch.      TT represents a continuously flowing ephemeris timescale independent of      any variations of the Earth's rotation, and is adjusted from UT      using historical and predictive models of those variations. |


* [AstroTime](#AstroTime)
    * [new AstroTime(date)](#new_AstroTime_new)
    * [.toString()](#AstroTime+toString) ⇒ <code>string</code>
    * [.AddDays(days)](#AstroTime+AddDays) ⇒ [<code>AstroTime</code>](#AstroTime)


* * *

<a name="new_AstroTime_new"></a>

### new AstroTime(date)

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> | A JavaScript Date object or a numeric UTC value expressed in J2000 days. |


* * *

<a name="AstroTime+toString"></a>

### astroTime.toString() ⇒ <code>string</code>
Formats an `AstroTime` object as an [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601)
date/time string in UTC, to millisecond resolution.
Example: `2018-08-17T17:22:04.050Z`

**Kind**: instance method of [<code>AstroTime</code>](#AstroTime)  

* * *

<a name="AstroTime+AddDays"></a>

### astroTime.AddDays(days) ⇒ [<code>AstroTime</code>](#AstroTime)
Returns a new `AstroTime` object adjusted by the floating point number of days.
Does NOT modify the original `AstroTime` object.

**Kind**: instance method of [<code>AstroTime</code>](#AstroTime)  

| Param | Type | Description |
| --- | --- | --- |
| days | <code>number</code> | The floating point number of days by which to adjust the given date and time.      Positive values adjust the date toward the future, and      negative values adjust the date toward the past. |


* * *

<a name="Vector"></a>

## Vector
Holds the Cartesian coordinates of a vector in 3D space,
along with the time at which the vector is valid.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| x | <code>number</code> | The x-coordinate expressed in astronomical units (AU). |
| y | <code>number</code> | The y-coordinate expressed in astronomical units (AU). |
| z | <code>number</code> | The z-coordinate expressed in astronomical units (AU). |
| t | [<code>AstroTime</code>](#AstroTime) | The time at which the vector is valid. |


* * *

<a name="Vector+Length"></a>

### vector.Length() ⇒ <code>number</code>
Returns the length of the vector in astronomical units (AU).

**Kind**: instance method of [<code>Vector</code>](#Vector)  

* * *

<a name="Spherical"></a>

## Spherical
Holds spherical coordinates: latitude, longitude, distance.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| lat | <code>number</code> | The latitude angle: -90..+90 degrees. |
| lon | <code>number</code> | The longitude angle: 0..360 degrees. |
| dist | <code>number</code> | Distance in AU. |


* * *

<a name="EquatorialCoordinates"></a>

## EquatorialCoordinates
Holds right ascension, declination, and distance of a celestial object.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| ra | <code>number</code> | Right ascension in sidereal hours: [0, 24). |
| dec | <code>number</code> | Declination in degrees: [-90, +90]. |
| dist | <code>number</code> | Distance to the celestial object expressed in      <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU). |


* * *

<a name="RotationMatrix"></a>

## RotationMatrix
Contains a rotation matrix that can be used to transform one coordinate system to another.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| rot | <code>Array.&lt;Array.&lt;number&gt;&gt;</code> | A normalized 3x3 rotation matrix. |


* * *

<a name="HorizontalCoordinates"></a>

## HorizontalCoordinates
Holds azimuth (compass direction) and altitude (angle above/below the horizon)
of a celestial object as seen by an observer at a particular location on the Earth's surface.
Also holds right ascension and declination of the same object.
All of these coordinates are optionally adjusted for atmospheric refraction;
therefore the right ascension and declination values may not exactly match
those found inside a corresponding [EquatorialCoordinates](#EquatorialCoordinates) object.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| azimuth | <code>number</code> | A horizontal compass direction angle in degrees measured starting at north      and increasing positively toward the east.      The value is in the range [0, 360).      North = 0, east = 90, south = 180, west = 270. |
| altitude | <code>number</code> | A vertical angle in degrees above (positive) or below (negative) the horizon.      The value is in the range [-90, +90].      The altitude angle is optionally adjusted upward due to atmospheric refraction. |
| ra | <code>number</code> | The right ascension of the celestial body in sidereal hours.      The value is in the reange [0, 24).      If `altitude` was adjusted for atmospheric reaction, `ra`      is likewise adjusted. |
| dec | <code>number</code> | The declination of of the celestial body in degrees.      The value in the range [-90, +90].      If `altitude` was adjusted for atmospheric reaction, `dec`      is likewise adjusted. |


* * *

<a name="EclipticCoordinates"></a>

## EclipticCoordinates
Holds ecliptic coordinates of a celestial body.
The origin and date of the coordinate system may vary depending on the caller's usage.
In general, ecliptic coordinates are measured with respect to the mean plane of the Earth's
orbit around the Sun.
Includes Cartesian coordinates `(ex, ey, ez)` measured in
<a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU)
and spherical coordinates `(elon, elat)` measured in degrees.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| ex | <code>number</code> | The Cartesian x-coordinate of the body in astronomical units (AU).      The x-axis is within the ecliptic plane and is oriented in the direction of the      <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>. |
| ey | <code>number</code> | The Cartesian y-coordinate of the body in astronomical units (AU).      The y-axis is within the ecliptic plane and is oriented 90 degrees      counterclockwise from the equinox, as seen from above the Sun's north pole. |
| ez | <code>number</code> | The Cartesian z-coordinate of the body in astronomical units (AU).      The z-axis is oriented perpendicular to the ecliptic plane,      along the direction of the Sun's north pole. |
| elat | <code>number</code> | The ecliptic latitude of the body in degrees.      This is the angle north or south of the ecliptic plane.      The value is in the range [-90, +90].      Positive values are north and negative values are south. |
| elon | <code>number</code> | The ecliptic longitude of the body in degrees.      This is the angle measured counterclockwise around the ecliptic plane,      as seen from above the Sun's north pole.      This is the same direction that the Earth orbits around the Sun.      The angle is measured starting at 0 from the equinox and increases      up to 360 degrees. |


* * *

<a name="Observer"></a>

## Observer
Represents the geographic location of an observer on the surface of the Earth.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| latitude | <code>number</code> | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude | <code>number</code> | The observer's geographic longitude in degrees east of the prime meridian      passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height | <code>number</code> | The observer's elevation above mean sea level, expressed in meters. |


* * *

<a name="IlluminationInfo"></a>

## IlluminationInfo
Contains information about the apparent brightness and sunlit phase of a celestial object.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time pertaining to the other calculated values in this object. |
| mag | <code>number</code> | The <a href="https://en.wikipedia.org/wiki/Apparent_magnitude">apparent visual magnitude</a> of the celestial body. |
| phase_angle | <code>number</code> | The angle in degrees as seen from the center of the celestial body between the Sun and the Earth.      The value is always in the range 0 to 180.      The phase angle provides a measure of what fraction of the body's face appears      illuminated by the Sun as seen from the Earth.      When the observed body is the Sun, the `phase` property is set to 0,      although this has no physical meaning because the Sun emits, rather than reflects, light.      When the phase is near 0 degrees, the body appears "full".      When it is 90 degrees, the body appears "half full".      And when it is 180 degrees, the body appears "new" and is very difficult to see      because it is both dim and lost in the Sun's glare as seen from the Earth. |
| phase_fraction | <code>number</code> | The fraction of the body's face that is illuminated by the Sun, as seen from the Earth.      Calculated from `phase_angle` for convenience.      This value ranges from 0 to 1. |
| helio_dist | <code>number</code> | The distance between the center of the Sun and the center of the body in      <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU). |
| geo_dist | <code>number</code> | The distance between the center of the Earth and the center of the body in AU. |
| gc | [<code>Vector</code>](#Vector) | Geocentric coordinates: the 3D vector from the center of the Earth to the center of the body.      The components are in expressed in AU and are oriented with respect to the J2000 equatorial plane. |
| hc | [<code>Vector</code>](#Vector) | Heliocentric coordinates: The 3D vector from the center of the Sun to the center of the body.      Like `gc`, `hc` is expressed in AU and oriented with respect      to the J2000 equatorial plane. |
| ring_tilt | <code>number</code> \| <code>null</code> | For Saturn, this is the angular tilt of the planet's rings in degrees away      from the line of sight from the Earth. When the value is near 0, the rings      appear edge-on from the Earth and are therefore difficult to see.      When `ring_tilt` approaches its maximum value (about 27 degrees),      the rings appear widest and brightest from the Earth.      Unlike the <a href="https://ssd.jpl.nasa.gov/horizons.cgi">JPL Horizons</a> online tool,      this library includes the effect of the ring tilt angle in the calculated value      for Saturn's visual magnitude.      For all bodies other than Saturn, the value of `ring_tilt` is `null`. |


* * *

<a name="MoonQuarter"></a>

## MoonQuarter
Represents a quarter lunar phase, along with when it occurs.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| quarter | <code>number</code> | An integer as follows:      0 = new moon,      1 = first quarter,      2 = full moon,      3 = third quarter. |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the quarter lunar phase. |


* * *

<a name="HourAngleEvent"></a>

## HourAngleEvent
Returns information about an occurrence of a celestial body
reaching a given hour angle as seen by an observer at a given
location on the surface of the Earth.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the celestial body reaching the hour angle. |
| hor | [<code>HorizontalCoordinates</code>](#HorizontalCoordinates) | Topocentric horizontal coordinates for the body      at the time indicated by the `time` property. |


* * *

<a name="SeasonInfo"></a>

## SeasonInfo
Represents the dates and times of the two solstices
and the two equinoxes in a given calendar year.
These four events define the changing of the seasons on the Earth.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| mar_equinox | [<code>AstroTime</code>](#AstroTime) | The date and time of the March equinox in the given calendar year.      This is the moment in March that the plane of the Earth's equator passes      through the center of the Sun; thus the Sun's declination      changes from a negative number to a positive number.      The March equinox defines      the beginning of spring in the northern hemisphere and      the beginning of autumn in the southern hemisphere. |
| jun_solstice | [<code>AstroTime</code>](#AstroTime) | The date and time of the June solstice in the given calendar year.      This is the moment in June that the Sun reaches its most positive      declination value.      At this moment the Earth's north pole is most tilted most toward the Sun.      The June solstice defines      the beginning of summer in the northern hemisphere and      the beginning of winter in the southern hemisphere. |
| sep_equinox | [<code>AstroTime</code>](#AstroTime) | The date and time of the September equinox in the given calendar year.      This is the moment in September that the plane of the Earth's equator passes      through the center of the Sun; thus the Sun's declination      changes from a positive number to a negative number.      The September equinox defines      the beginning of autumn in the northern hemisphere and      the beginning of spring in the southern hemisphere. |
| dec_solstice | [<code>AstroTime</code>](#AstroTime) | The date and time of the December solstice in the given calendar year.      This is the moment in December that the Sun reaches its most negative      declination value.      At this moment the Earth's south pole is tilted most toward the Sun.      The December solstice defines      the beginning of winter in the northern hemisphere and      the beginning of summer in the southern hemisphere. |


* * *

<a name="ElongationEvent"></a>

## ElongationEvent
Represents the angular separation of a body from the Sun as seen from the Earth
and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.

**Kind**: global class  
**See**: [Elongation](#Elongation)  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the observation. |
| visibility | <code>string</code> | Either `"morning"` or `"evening"`,      indicating when the body is most easily seen. |
| elongation | <code>number</code> | The angle in degrees, as seen from the center of the Earth,      of the apparent separation between the body and the Sun.      This angle is measured in 3D space and is not projected onto the ecliptic plane.      When `elongation` is less than a few degrees, the body is very      difficult to see from the Earth because it is lost in the Sun's glare.      The elongation is always in the range [0, 180]. |
| ecliptic_separation | <code>number</code> | The absolute value of the difference between the body's ecliptic longitude      and the Sun's ecliptic longitude, both as seen from the center of the Earth.      This angle measures around the plane of the Earth's orbit (the ecliptic),      and ignores how far above or below that plane the body is.      The ecliptic separation is measured in degrees and is always in the range [0, 180]. |


* * *

<a name="Apsis"></a>

## Apsis
Represents a closest or farthest point in a body's orbit around its primary.
For a planet orbiting the Sun, this is a perihelion or aphelion, respectively.
For the Moon orbiting the Earth, this is a perigee or apogee, respectively.

**Kind**: global class  
**See**

- [SearchLunarApsis](#SearchLunarApsis)
- [NextLunarApsis](#NextLunarApsis)

**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the apsis. |
| kind | <code>number</code> | For a closest approach (perigee or perihelion), `kind` is 0.      For a farthest distance event (apogee or aphelion), `kind` is 1. |
| dist_au | <code>number</code> | The distance between the centers of the two bodies in astronomical units (AU). |
| dist_km | <code>number</code> | The distance between the centers of the two bodies in kilometers. |


* * *

<a name="ConstellationInfo"></a>

## ConstellationInfo
Reports the constellation that a given celestial point lies within.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| symbol | <code>string</code> | 3-character mnemonic symbol for the constellation, e.g. "Ori". |
| name | <code>string</code> | Full name of constellation, e.g. "Orion". |
| ra1875 | <code>number</code> | Right ascension expressed in B1875 coordinates. |
| dec1875 | <code>number</code> | Declination expressed in B1875 coordinates. |


* * *

<a name="LunarEclipseInfo"></a>

## LunarEclipseInfo
Returns information about a lunar eclipse.

Returned by [SearchLunarEclipse](#SearchLunarEclipse) or [NextLunarEclipse](#NextLunarEclipse)
to report information about a lunar eclipse event.
When a lunar eclipse is found, it is classified as penumbral, partial, or total.
Penumbral eclipses are difficult to observe, because the moon is only slightly dimmed
by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
Total eclipses occur when the entire Moon passes into the Earth's umbra.

The `kind` field thus holds one of the strings `"penumbral"`, `"partial"`,
or `"total"`, depending on the kind of lunar eclipse found.

Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.

Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
of the eclipse, which is half of the amount of time the eclipse spends in each
phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
By converting from minutes to days, and subtracting/adding with `peak`, the caller
may determine the date and time of the beginning/end of each eclipse phase.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| kind | <code>string</code> | The type of lunar eclipse found. |
| peak | [<code>AstroTime</code>](#AstroTime) | The time of the eclipse at its peak. |
| sd_penum | <code>number</code> | The semi-duration of the penumbral phase in minutes. |
| sd_partial | <code>number</code> | The semi-duration of the penumbral phase in minutes, or 0.0 if none. |
| sd_total | <code>number</code> | The semi-duration of the penumbral phase in minutes, or 0.0 if none. |


* * *

<a name="GlobalSolarEclipseInfo"></a>

## GlobalSolarEclipseInfo
Reports the time and geographic location of the peak of a solar eclipse.

    Returned by [SearchGlobalSolarEclipse](#SearchGlobalSolarEclipse) or [NextGlobalSolarEclipse](#NextGlobalSolarEclipse)
    to report information about a solar eclipse event.

    Field `peak` holds the date and time of the peak of the eclipse, defined as
    the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.

    The eclipse is classified as partial, annular, or total, depending on the
    maximum amount of the Sun's disc obscured, as seen at the peak location
    on the surface of the Earth.

    The `kind` field thus holds one of the strings `"partial"`, `"annular"`, or `"total"`.
    A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
    An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
    to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
    A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
    observes either a total or annular eclipse.

    If `kind` is `"total"` or `"annular"`, the `latitude` and `longitude`
    fields give the geographic coordinates of the center of the Moon's shadow projected
    onto the daytime side of the Earth at the instant of the eclipse's peak.
    If `kind` has any other value, `latitude` and `longitude` are undefined and should
    not be used.

**Kind**: global class  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| kind | <code>string</code> | One of the following string values: `"partial"`, `"annular"`, `"total"`. |
| peak | [<code>AstroTime</code>](#AstroTime) | The date and time of the peak of the eclipse, defined as the instant         when the axis of the Moon's shadow cone passes closest to the Earth's center. |
| distance | <code>number</code> | The distance in kilometers between the axis of the Moon's shadow cone         and the center of the Earth at the time indicated by `peak`. |
| latitude | <code>undefined</code> \| <code>number</code> | If `kind` holds `"total"`, the geographic latitude in degrees         where the center of the Moon's shadow falls on the Earth at the         time indicated by `peak`; otherwise, `latitude` holds `undefined`. |
| longitude | <code>undefined</code> \| <code>number</code> | If `kind` holds `"total"`, the geographic longitude in degrees         where the center of the Moon's shadow falls on the Earth at the         time indicated by `peak`; otherwise, `longitude` holds `undefined`. |


* * *

<a name="EclipseEvent"></a>

## EclipseEvent
**Kind**: global class  
**Brief**: Holds a time and the observed altitude of the Sun at that time.

When reporting a solar eclipse observed at a specific location on the Earth
(a "local" solar eclipse), a series of events occur. In addition
to the time of each event, it is important to know the altitude of the Sun,
because each event may be invisible to the observer if the Sun is below
the horizon (i.e. it at night).

If `altitude` is negative, the event is theoretical only; it would be
visible if the Earth were transparent, but the observer cannot actually see it.
If `altitude` is positive but less than a few degrees, visibility will be impaired by
atmospheric interference (sunrise or sunset conditions).  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the event. |
| altitude | <code>number</code> | The angular altitude of the center of the Sun above/below the horizon, at `time`,      corrected for atmospheric refraction and expressed in degrees. |


* * *

<a name="LocalSolarEclipseInfo"></a>

## LocalSolarEclipseInfo
**Kind**: global class  
**Brief**: Information about a solar eclipse as seen by an observer at a given time and geographic location.

Returned by [SearchLocalSolarEclipse](#SearchLocalSolarEclipse) or [NextLocalSolarEclipse](#NextLocalSolarEclipse)
to report information about a solar eclipse as seen at a given geographic location.

When a solar eclipse is found, it is classified by setting `kind`
to `"partial"`, `"annular"`, or `"total"`.
A partial solar eclipse is when the Moon does not line up directly enough with the Sun
to completely block the Sun's light from reaching the observer.
An annular eclipse occurs when the Moon's disc is completely visible against the Sun
but the Moon is too far away to completely block the Sun's light; this leaves the
Sun with a ring-like appearance.
A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
Sun just right to completely block all sunlight from reaching the observer.

There are 5 "event" fields, each of which contains a time and a solar altitude.
Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
The fields `partial_begin` and `partial_end` are always set, and indicate when
the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
`total_begin` and `total_end` indicate when the total/annular phase begins/ends.
When an event field is valid, the caller must also check its `altitude` field to
see whether the Sun is above the horizon at the time indicated by the `time` field.
See #EclipseEvent for more information.  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| kind | <code>string</code> | The type of solar eclipse found: `"partial"`, `"annular"`, or `"total"`. |
| partial_begin | [<code>EclipseEvent</code>](#EclipseEvent) | The time and Sun altitude at the beginning of the eclipse. |
| total_begin | [<code>EclipseEvent</code>](#EclipseEvent) | If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise undefined. |
| peak | [<code>EclipseEvent</code>](#EclipseEvent) | The time and Sun altitude when the eclipse reaches its peak. |
| total_end | [<code>EclipseEvent</code>](#EclipseEvent) | If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise undefined. |
| partial_end | [<code>EclipseEvent</code>](#EclipseEvent) | The time and Sun altitude at the end of the eclipse. |


* * *

<a name="TransitInfo"></a>

## TransitInfo
**Kind**: global class  
**Brief**: Information about a transit of Mercury or Venus, as seen from the Earth.

Returned by [SearchTransit](#SearchTransit) or [NextTransit](#NextTransit) to report
information about a transit of Mercury or Venus.
A transit is when Mercury or Venus passes between the Sun and Earth so that
the other planet is seen in silhouette against the Sun.

The calculations are performed from the point of view of a geocentric observer.  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| start | [<code>AstroTime</code>](#AstroTime) | The date and time at the beginning of the transit.      This is the moment the planet first becomes visible against the Sun in its background. |
| peak | [<code>AstroTime</code>](#AstroTime) | When the planet is most aligned with the Sun, as seen from the Earth. |
| finish | [<code>AstroTime</code>](#AstroTime) | The date and time at the end of the transit.      This is the moment the planet is last seen against the Sun in its background. |
| separation; | <code>number</code> | The minimum angular separation, in arcminutes, between the centers of the Sun and the planet.      This angle pertains to the time stored in `peak`. |


* * *

<a name="Bodies"></a>

## Bodies : <code>Array.&lt;string&gt;</code>
An array of strings, each a name of a supported astronomical body.
     Not all bodies are valid for all functions, but any string not in this
     list is not supported at all.

**Kind**: global constant  

* * *

<a name="AngleBetween"></a>

## AngleBetween(a, b) ⇒ <code>number</code>
Calculates the angle in degrees between two vectors.
The angle is measured in the plane that contains both vectors.

**Kind**: global function  
**Returns**: <code>number</code> - The angle between the two vectors expressed in degrees.
     The value is in the range [0, 180].  

| Param | Type | Description |
| --- | --- | --- |
| a | [<code>Vector</code>](#Vector) | The first of a pair of vectors between which to measure an angle. |
| b | [<code>Vector</code>](#Vector) | The second of a pair of vectors between which to measure an angle. |


* * *

<a name="TerrestrialTime"></a>

## TerrestrialTime(ut) ⇒ <code>number</code>
Calculates Terrestrial Time (TT) from Universal Time (UT).

**Kind**: global function  
**Returns**: <code>number</code> - A Terrestrial Time expressed as a floating point number of days since the 2000.0 epoch.  

| Param | Type | Description |
| --- | --- | --- |
| ut | <code>number</code> | The Universal Time expressed as a floating point number of days since the 2000.0 epoch. |


* * *

<a name="MakeTime"></a>

## MakeTime(date) ⇒ [<code>AstroTime</code>](#AstroTime)
Given a Date object or a number days since noon (12:00) on January 1, 2000 (UTC),
this function creates an [AstroTime](#AstroTime) object.
Given an [AstroTime](#AstroTime) object, returns the same object unmodified.
Use of this function is not required for any of the other exposed functions in this library,
because they all guarantee converting date/time parameters to AstroTime
as needed. However, it may be convenient for callers who need to understand
the difference between UTC and TT (Terrestrial Time). In some use cases,
converting once to AstroTime format and passing the result into multiple
function calls may be more efficient than passing in native JavaScript Date objects.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | A Date object, a number of UTC days since the J2000 epoch (noon on January 1, 2000),      or an AstroTime object. See remarks above. |


* * *

<a name="MakeSpherical"></a>

## MakeSpherical(lat, lon, dist) ⇒ [<code>Spherical</code>](#Spherical)
Create spherical coordinates.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| lat | <code>number</code> | The angular distance above or below the reference plane, in degrees. |
| lon | <code>number</code> | The angular distance around the reference plane, in degrees. |
| dist | <code>number</code> | A radial distance in AU. |


* * *

<a name="MakeRotation"></a>

## MakeRotation(rot) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Creates a rotation matrix that can be used to transform one coordinate system to another.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| rot | <code>Array.&lt;Array.&lt;number&gt;&gt;</code> | An array [3][3] of numbers. Defines a rotation matrix used to premultiply      a 3D vector to reorient it into another coordinate system. |


* * *

<a name="Horizon"></a>

## Horizon(date, observer, ra, dec, refraction) ⇒ [<code>HorizontalCoordinates</code>](#HorizontalCoordinates)
Given a date and time, a geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial body,
returns horizontal coordinates (azimuth and altitude angles) for that body
as seen by that observer. Allows optional correction for atmospheric refraction.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which to find horizontal coordinates. |
| observer | [<code>Observer</code>](#Observer) | The location of the observer for which to find horizontal coordinates. |
| ra | <code>number</code> | Right ascension in sidereal hours of the celestial object,      referred to the mean equinox of date for the J2000 epoch. |
| dec | <code>number</code> | Declination in degrees of the celestial object,      referred to the mean equator of date for the J2000 epoch.      Positive values are north of the celestial equator and negative values are south. |
| refraction | <code>string</code> | If omitted or has a false-like value (false, null, undefined, etc.)      the calculations are performed without any correction for atmospheric      refraction. If the value is the string `"normal"`,      uses the recommended refraction correction based on Meeus "Astronomical Algorithms"      with a linear taper more than 1 degree below the horizon. The linear      taper causes the refraction to linearly approach 0 as the altitude of the      body approaches the nadir (-90 degrees).      If the value is the string `"jplhor"`, uses a JPL Horizons      compatible formula. This is the same algorithm as `"normal"`,      only without linear tapering; this can result in physically impossible      altitudes of less than -90 degrees, which may cause problems for some applications.      (The `"jplhor"` option was created for unit testing against data      generated by JPL Horizons, and is otherwise not recommended for use.) |


* * *

<a name="MakeObserver"></a>

## MakeObserver(latitude_degrees, longitude_degrees, height_in_meters)
Creates an [Observer](#Observer) object that represents a location
on the surface of the Earth from which observations are made.

**Kind**: global function  

| Param | Type | Default | Description |
| --- | --- | --- | --- |
| latitude_degrees | <code>number</code> |  | The observer's geographic latitude in degrees north of the Earth's equator.      The value is negative for observers south of the equator.      Must be in the range -90 to +90. |
| longitude_degrees | <code>number</code> |  | The observer's geographic longitude in degrees east of the prime meridian      passing through Greenwich, England.      The value is negative for observers west of the prime meridian.      The value should be kept in the range -180 to +180 to minimize floating point errors. |
| height_in_meters | <code>number</code> | <code>0</code> | The observer's elevation above mean sea level, expressed in meters.      If omitted, the elevation is assumed to be 0 meters. |


* * *

<a name="SunPosition"></a>

## SunPosition(date) ⇒ [<code>EclipticCoordinates</code>](#EclipticCoordinates)
Returns apparent geocentric true ecliptic coordinates of date for the Sun.
<i>Geocentric</i> means coordinates as the Sun would appear to a hypothetical observer
at the center of the Earth.
<i>Ecliptic coordinates of date</i> are measured along the plane of the Earth's mean
orbit around the Sun, using the
<a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>
of the Earth as adjusted for precession and nutation of the Earth's
axis of rotation on the given date.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time at which to calculate the Sun's apparent location as seen from      the center of the Earth. |


* * *

<a name="Equator"></a>

## Equator(body, date, observer, ofdate, aberration) ⇒ [<code>EquatorialCoordinates</code>](#EquatorialCoordinates)
Returns topocentric equatorial coordinates (right ascension and declination)
in one of two different systems: J2000 or true-equator-of-date.
Allows optional correction for aberration.
Always corrects for light travel time (represents the object as seen by the observer
with light traveling to the Earth at finite speed, not where the object is right now).
<i>Topocentric</i> refers to a position as seen by an observer on the surface of the Earth.
This function corrects for
<a href="https://en.wikipedia.org/wiki/Parallax">parallax</a>
of the object between a geocentric observer and a topocentric observer.
This is most significant for the Moon, because it is so close to the Earth.
However, it can have a small effect on the apparent positions of other bodies.

**Kind**: global function  
**Returns**: [<code>EquatorialCoordinates</code>](#EquatorialCoordinates) - The topocentric coordinates of the body as adjusted for the given observer.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the body for which to find equatorial coordinates.      Not allowed to be `"Earth"`. |
| date | <code>Date</code> \| <code>number</code> \| <code>Time</code> | Specifies the date and time at which the body is to be observed. |
| observer | [<code>Observer</code>](#Observer) | The location on the Earth of the observer.      Call [MakeObserver](#MakeObserver) to create an observer object. |
| ofdate | <code>bool</code> | Pass `true` to return equatorial coordinates of date,      i.e. corrected for precession and nutation at the given date.      This is needed to get correct horizontal coordinates when you call      [Horizon](#Horizon).      Pass `false` to return equatorial coordinates in the J2000 system. |
| aberration | <code>bool</code> | Pass `true` to correct for      <a href="https://en.wikipedia.org/wiki/Aberration_of_light">aberration</a>,      or `false` to leave uncorrected. |


* * *

<a name="Ecliptic"></a>

## Ecliptic(gx, gy, gz) ⇒ [<code>EclipticCoordinates</code>](#EclipticCoordinates)
Given J2000 equatorial Cartesian coordinates,
returns J2000 ecliptic latitude, longitude, and cartesian coordinates.
You can call [GeoVector](#GeoVector) and use its (x, y, z) return values
to pass into this function.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| gx | <code>number</code> | The x-coordinate of a 3D vector in the J2000 equatorial coordinate system. |
| gy | <code>number</code> | The y-coordinate of a 3D vector in the J2000 equatorial coordinate system. |
| gz | <code>number</code> | The z-coordinate of a 3D vector in the J2000 equatorial coordinate system. |


* * *

<a name="GeoMoon"></a>

## GeoMoon(date) ⇒ [<code>Vector</code>](#Vector)
Calculates the geocentric Cartesian coordinates for the Moon in the J2000 equatorial system.
Based on the Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
which in turn derives from E. W. Brown's lunar theories.
Adapted from Turbo Pascal code from the book
<a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a>
by Montenbruck and Pfleger.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which to calculate the Moon's geocentric position. |


* * *

<a name="HelioVector"></a>

## HelioVector(body, date) ⇒ [<code>Vector</code>](#Vector)
Calculates heliocentric (i.e., with respect to the center of the Sun)
Cartesian coordinates in the J2000 equatorial system of a celestial
body at a specified time. The position is not corrected for light travel time or aberration.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | One of the strings      `"Sun"`, `"Moon"`, `"Mercury"`, `"Venus"`,      `"Earth"`, `"Mars"`, `"Jupiter"`, `"Saturn"`,      `"Uranus"`, `"Neptune"`, `"Pluto"`,      `"SSB"`, or `"EMB"`. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which the body's position is to be calculated. |


* * *

<a name="HelioDistance"></a>

## HelioDistance(body, date) ⇒ <code>number</code>
Calculates the distance between a body and the Sun at a given time.

Given a date and time, this function calculates the distance between
the center of `body` and the center of the Sun.
For the planets Mercury through Neptune, this function is significantly
more efficient than calling [HelioVector](#HelioVector) followed by taking the length
of the resulting vector.

**Kind**: global function  
**Returns**: <code>number</code> - The heliocentric distance in AU.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | A body for which to calculate a heliocentric distance:      the Sun, Moon, or any of the planets. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which to calculate the heliocentric distance. |


* * *

<a name="GeoVector"></a>

## GeoVector(body, date, aberration) ⇒ [<code>Vector</code>](#Vector)
Calculates geocentric (i.e., with respect to the center of the Earth)
Cartesian coordinates in the J2000 equatorial system of a celestial
body at a specified time. The position is always corrected for light travel time:
this means the position of the body is "back-dated" based on how long it
takes light to travel from the body to an observer on the Earth.
Also, the position can optionally be corrected for aberration, an effect
causing the apparent direction of the body to be shifted based on
transverse movement of the Earth with respect to the rays of light
coming from that body.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | One of the strings      `"Sun"`, `"Moon"`, `"Mercury"`, `"Venus"`,      `"Earth"`, `"Mars"`, `"Jupiter"`, `"Saturn"`,      `"Uranus"`, `"Neptune"`, or `"Pluto"`. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which the body's position is to be calculated. |
| aberration | <code>bool</code> | Pass `true` to correct for      <a href="https://en.wikipedia.org/wiki/Aberration_of_light">aberration</a>,      or `false` to leave uncorrected. |


* * *

<a name="Search"></a>

## Search(func, t1, t2, options) ⇒ <code>null</code> \| [<code>AstroTime</code>](#AstroTime)
Search for next time <i>t</i> (such that <i>t</i> is between `t1` and `t2`)
that `func(t)` crosses from a negative value to a non-negative value.
The given function must have "smooth" behavior over the entire inclusive range [`t1`, `t2`],
meaning that it behaves like a continuous differentiable function.
It is not required that `t1` &lt; `t2`; `t1` &gt; `t2`
allows searching backward in time.
Note: `t1` and `t2` must be chosen such that there is no possibility
of more than one zero-crossing (ascending or descending), or it is possible
that the "wrong" event will be found (i.e. not the first event after t1)
or even that the function will return null, indicating that no event was found.

**Kind**: global function  
**Returns**: <code>null</code> \| [<code>AstroTime</code>](#AstroTime) - If the search is successful, returns the date and time of the solution.
     If the search fails, returns null.  

| Param | Type | Description |
| --- | --- | --- |
| func | <code>ContinuousFunction</code> | The function to find an ascending zero crossing for.      The function must accept a single parameter of type [AstroTime](#AstroTime)      and return a numeric value. |
| t1 | [<code>AstroTime</code>](#AstroTime) | The lower time bound of a search window. |
| t2 | [<code>AstroTime</code>](#AstroTime) | The upper time bound of a search window. |
| options | <code>null</code> \| [<code>SearchOptions</code>](#SearchOptions) | Options that can tune the behavior of the search.      Most callers can omit this argument or pass in `null`. |


* * *

<a name="SearchSunLongitude"></a>

## SearchSunLongitude(targetLon, dateStart, limitDays) ⇒ [<code>AstroTime</code>](#AstroTime) \| <code>null</code>
Searches for the moment in time when the center of the Sun reaches a given apparent
ecliptic longitude, as seen from the center of the Earth, within a given range of dates.
This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call [Seasons](#Seasons)
to calculate equinoxes and solstices for a given calendar year.
`SearchSunLongitude` is more general in that it allows searching for arbitrary longitude values.

**Kind**: global function  
**Returns**: [<code>AstroTime</code>](#AstroTime) \| <code>null</code> - The date and time when the Sun reaches the apparent ecliptic longitude `targetLon`
     within the range of times specified by `dateStart` and `limitDays`.
     If the Sun does not reach the target longitude within the specified time range, or the
     time range is excessively wide, the return value is `null`.
     To avoid a `null` return value, the caller must pick a time window around
     the event that is within a few days but not so small that the event might fall outside the window.  

| Param | Type | Description |
| --- | --- | --- |
| targetLon | <code>number</code> | The desired ecliptic longitude of date in degrees.      This may be any value in the range [0, 360), although certain      values have conventional meanings:      When `targetLon` is 0, finds the March equinox,      which is the moment spring begins in the northern hemisphere      and the beginning of autumn in the southern hemisphere.      When `targetLon` is 180, finds the September equinox,      which is the moment autumn begins in the northern hemisphere and      spring begins in the southern hemisphere.      When `targetLon` is 90, finds the northern solstice, which is the      moment summer begins in the northern hemisphere and winter      begins in the southern hemisphere.      When `targetLon` is 270, finds the southern solstice, which is the      moment winter begins in the northern hemisphere and summer      begins in the southern hemisphere. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | A date and time known to be earlier than the desired longitude event. |
| limitDays | <code>number</code> | A floating point number of days, which when added to `dateStart`,      yields a date and time known to be after the desired longitude event. |


* * *

<a name="LongitudeFromSun"></a>

## LongitudeFromSun(body, date) ⇒ <code>number</code>
Calculates the ecliptic longitude difference
between the given body and the Sun as seen from
the Earth at a given moment in time.
The returned value ranges [0, 360) degrees.
By definition, the Earth and the Sun are both in the plane of the ecliptic.
Ignores the height of the `body` above or below the ecliptic plane;
the resulting angle is measured around the ecliptic plane for the "shadow"
of the body onto that plane.

**Kind**: global function  
**Returns**: <code>number</code> - An angle in degrees in the range [0, 360).
     Values less than 180 indicate that the body is to the east
     of the Sun as seen from the Earth; that is, the body sets after
     the Sun does and is visible in the evening sky.
     Values greater than 180 indicate that the body is to the west of
     the Sun and is visible in the morning sky.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a supported celestial body other than the Earth. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The time at which the relative longitude is to be found. |


* * *

<a name="AngleFromSun"></a>

## AngleFromSun(body, date) ⇒ <code>number</code>
Returns the full angle seen from
the Earth, between the given body and the Sun.
Unlike [LongitudeFromSun](#LongitudeFromSun), this function does not
project the body's "shadow" onto the ecliptic;
the angle is measured in 3D space around the plane that
contains the centers of the Earth, the Sun, and `body`.

**Kind**: global function  
**Returns**: <code>number</code> - An angle in degrees in the range [0, 180].  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a supported celestial body other than the Earth. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The time at which the angle from the Sun is to be found. |


* * *

<a name="EclipticLongitude"></a>

## EclipticLongitude(body, date) ⇒ <code>number</code>
Calculates heliocentric ecliptic longitude based on the J2000 equinox.

**Kind**: global function  
**Returns**: <code>number</code> - The ecliptic longitude angle of the body in degrees measured counterclockwise around the mean
     plane of the Earth's orbit, as seen from above the Sun's north pole.
     Ecliptic longitude starts at 0 at the J2000
     <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a> and
     increases in the same direction the Earth orbits the Sun.
     The returned value is always in the range [0, 360).  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a celestial body other than the Sun. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which to calculate the ecliptic longitude. |


* * *

<a name="Illumination"></a>

## Illumination(body, date) ⇒ [<code>IlluminationInfo</code>](#IlluminationInfo)
Calculates the phase angle, visual maginitude,
and other values relating to the body's illumination
at the given date and time, as seen from the Earth.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the celestial body being observed.      Not allowed to be `"Earth"`. |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which to calculate the illumination data for the given body. |


* * *

<a name="SearchRelativeLongitude"></a>

## SearchRelativeLongitude(body, targetRelLon, startDate) ⇒ [<code>AstroTime</code>](#AstroTime)
Searches for the date and time the relative ecliptic longitudes of
the specified body and the Earth, as seen from the Sun, reach a certain
difference. This function is useful for finding conjunctions and oppositions
of the planets. For the opposition of a superior planet (Mars, Jupiter, ..., Pluto),
or the inferior conjunction of an inferior planet (Mercury, Venus),
call with `targetRelLon` = 0. The 0 value indicates that both
planets are on the same ecliptic longitude line, ignoring the other planet's
distance above or below the plane of the Earth's orbit.
For superior conjunctions, call with `targetRelLon` = 180.
This means the Earth and the other planet are on opposite sides of the Sun.

**Kind**: global function  
**Returns**: [<code>AstroTime</code>](#AstroTime) - The time when the Earth and the body next reach the specified relative longitudes.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a planet other than the Earth. |
| targetRelLon | <code>number</code> | The desired angular difference in degrees between the ecliptic longitudes      of `body` and the Earth. Must be in the range (-180, +180]. |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which to find the next occurrence of the      body and the Earth reaching the desired relative longitude. |


* * *

<a name="MoonPhase"></a>

## MoonPhase(date) ⇒ <code>number</code>
Determines the moon's phase expressed as an ecliptic longitude.

**Kind**: global function  
**Returns**: <code>number</code> - A value in the range [0, 360) indicating the difference
     in ecliptic longitude between the center of the Sun and the
     center of the Moon, as seen from the center of the Earth.
     Certain longitude values have conventional meanings:

* 0 = new moon
* 90 = first quarter
* 180 = full moon
* 270 = third quarter  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for which to calculate the moon's phase. |


* * *

<a name="SearchMoonPhase"></a>

## SearchMoonPhase(targetLon, dateStart, limitDays) ⇒ [<code>AstroTime</code>](#AstroTime) \| <code>null</code>
Searches for the date and time that the Moon reaches a specified phase.
Lunar phases are defined in terms of geocentric ecliptic longitudes
with respect to the Sun.  When the Moon and the Sun have the same ecliptic
longitude, that is defined as a new moon. When the two ecliptic longitudes
are 180 degrees apart, that is defined as a full moon.
To enumerate quarter lunar phases, it is simpler to call
[SearchMoonQuarter](#SearchMoonQuarter) once, followed by repeatedly calling
[NextMoonQuarter](#NextMoonQuarter). `SearchMoonPhase` is only
necessary for finding other lunar phases than the usual quarter phases.

**Kind**: global function  
**Returns**: [<code>AstroTime</code>](#AstroTime) \| <code>null</code> - If the specified lunar phase occurs after `dateStart`
     and before `limitDays` days after `dateStart`,
     this function returns the date and time of the first such occurrence.
     Otherwise, it returns `null`.  

| Param | Type | Description |
| --- | --- | --- |
| targetLon | <code>number</code> | The difference in geocentric ecliptic longitude between the Sun and Moon      that specifies the lunar phase being sought. This can be any value      in the range [0, 360). Here are some helpful examples:      0 = new moon,      90 = first quarter,      180 = full moon,      270 = third quarter. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The beginning of the window of time in which to search. |
| limitDays | <code>number</code> | The floating point number of days after `dateStart`      that limits the window of time in which to search. |


* * *

<a name="SearchMoonQuarter"></a>

## SearchMoonQuarter(dateStart) ⇒ [<code>MoonQuarter</code>](#MoonQuarter)
Finds the first quarter lunar phase after the specified date and time.
The quarter lunar phases are: new moon, first quarter, full moon, and third quarter.
To enumerate quarter lunar phases, call `SearchMoonQuarter` once,
then pass its return value to [NextMoonQuarter](#NextMoonQuarter) to find the next
`MoonQuarter`. Keep calling `NextMoonQuarter` in a loop,
passing the previous return value as the argument to the next call.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which to find the first quarter lunar phase. |


* * *

<a name="NextMoonQuarter"></a>

## NextMoonQuarter(mq)
Given a [MoonQuarter](#MoonQuarter) object, finds the next consecutive
quarter lunar phase. See remarks in [SearchMoonQuarter](#SearchMoonQuarter)
for explanation of usage.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| mq | [<code>MoonQuarter</code>](#MoonQuarter) | The return value of a prior call to [MoonQuarter](#MoonQuarter) or `NextMoonQuarter`. |


* * *

<a name="SearchRiseSet"></a>

## SearchRiseSet(body, observer, direction, dateStart, limitDays) ⇒ [<code>AstroTime</code>](#AstroTime) \| <code>null</code>
Finds a rise or set time for the given body as
seen by an observer at the specified location on the Earth.
Rise time is defined as the moment when the top of the body
is observed to first appear above the horizon in the east.
Set time is defined as the moment the top of the body
is observed to sink below the horizon in the west.
The times are adjusted for typical atmospheric refraction conditions.

**Kind**: global function  
**Returns**: [<code>AstroTime</code>](#AstroTime) \| <code>null</code> - The date and time of the rise or set event, or null if no such event
     occurs within the specified time window.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the body to find the rise or set time for. |
| observer | [<code>Observer</code>](#Observer) | Specifies the geographic coordinates and elevation above sea level of the observer.      Call [MakeObserver](#MakeObserver) to create an observer object. |
| direction | <code>number</code> | Either +1 to find rise time or -1 to find set time.      Any other value will cause an exception to be thrown. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which the specified rise or set time is to be found. |
| limitDays | <code>number</code> | The fractional number of days after `dateStart` that limits      when the rise or set time is to be found. |


* * *

<a name="SearchHourAngle"></a>

## SearchHourAngle(body, observer, hourAngle, dateStart) ⇒ [<code>HourAngleEvent</code>](#HourAngleEvent)
Finds the next time the given body is seen to reach the specified
<a href="https://en.wikipedia.org/wiki/Hour_angle">hour angle</a>
by the given observer.
Providing `hourAngle` = 0 finds the next maximum altitude event (culmination).
Providing `hourAngle` = 12 finds the next minimum altitude event.
Note that, especially close to the Earth's poles, a body as seen on a given day
may always be above the horizon or always below the horizon, so the caller cannot
assume that a culminating object is visible nor that an object is below the horizon
at its minimum altitude.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of a celestial body other than the Earth. |
| observer | [<code>Observer</code>](#Observer) | Specifies the geographic coordinates and elevation above sea level of the observer.      Call [MakeObserver](#MakeObserver) to create an observer object. |
| hourAngle | <code>number</code> | The hour angle expressed in      <a href="https://en.wikipedia.org/wiki/Sidereal_time">sidereal</a>      hours for which the caller seeks to find the body attain.      The value must be in the range [0, 24).      The hour angle represents the number of sidereal hours that have      elapsed since the most recent time the body crossed the observer's local      <a href="https://en.wikipedia.org/wiki/Meridian_(astronomy)">meridian</a>.      This specifying `hourAngle` = 0 finds the moment in time      the body reaches the highest angular altitude in a given sidereal day. |
| dateStart | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which the desired hour angle crossing event      is to be found. |


* * *

<a name="Seasons"></a>

## Seasons(year) ⇒ [<code>SeasonInfo</code>](#SeasonInfo)
Finds the equinoxes and solstices for a given calendar year.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| year | <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The integer value or `AstroTime` object that specifies      the UTC calendar year for which to find equinoxes and solstices. |


* * *

<a name="Elongation"></a>

## Elongation(body) ⇒ [<code>ElongationEvent</code>](#ElongationEvent)
Calculates angular separation of a body from the Sun as seen from the Earth
and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.
See the return type [ElongationEvent](#ElongationEvent) for details.

This function is helpful for determining how easy
it is to view a planet away from the Sun's glare on a given date.
It also determines whether the object is visible in the morning or evening;
this is more important the smaller the elongation is.
It is also used to determine how far a planet is from opposition, conjunction, or quadrature.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The name of the observed body. Not allowed to be `"Earth"`. |


* * *

<a name="SearchMaxElongation"></a>

## SearchMaxElongation(body, startDate) ⇒ [<code>ElongationEvent</code>](#ElongationEvent)
Searches for the next maximum elongation event for Mercury or Venus
that occurs after the given start date. Calling with other values
of `body` will result in an exception.
Maximum elongation occurs when the body has the greatest
angular separation from the Sun, as seen from the Earth.
Returns an `ElongationEvent` object containing the date and time of the next
maximum elongation, the elongation in degrees, and whether
the body is visible in the morning or evening.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | Either `"Mercury"` or `"Venus"`. |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which to search for the next maximum elongation event. |


* * *

<a name="SearchPeakMagnitude"></a>

## SearchPeakMagnitude(body, startDate) ⇒ [<code>IlluminationInfo</code>](#IlluminationInfo)
Searches for the date and time Venus will next appear brightest as seen from the Earth.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | Currently only `"Venus"` is supported.      Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from Earth,      so peak magnitude events have little practical value for that planet.      The Moon reaches peak magnitude very close to full moon, which can be found using      [SearchMoonQuarter](#SearchMoonQuarter) or [SearchMoonPhase](#SearchMoonPhase).      The other planets reach peak magnitude very close to opposition,      which can be found using [SearchRelativeLongitude](#SearchRelativeLongitude). |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which to find the next peak magnitude event. |


* * *

<a name="SearchLunarApsis"></a>

## SearchLunarApsis(startDate) ⇒ [<code>Apsis</code>](#Apsis)
Finds the next perigee (closest approach) or apogee (farthest remove) of the Moon
that occurs after the specified date and time.

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| startDate | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time after which to find the next perigee or apogee. |


* * *

<a name="NextLunarApsis"></a>

## NextLunarApsis(apsis) ⇒ [<code>Apsis</code>](#Apsis)
Given a lunar apsis returned by an initial call to [SearchLunarApsis](#SearchLunarApsis),
or a previous call to `NextLunarApsis`, finds the next lunar apsis.
If the given apsis is a perigee, this function finds the next apogee, and vice versa.

**Kind**: global function  
**Returns**: [<code>Apsis</code>](#Apsis) - The successor apogee for the given perigee, or the successor perigee for the given apogee.  

| Param | Type | Description |
| --- | --- | --- |
| apsis | [<code>Apsis</code>](#Apsis) | A lunar perigee or apogee event. |


* * *

<a name="SearchPlanetApsis"></a>

## SearchPlanetApsis(body, startTime) ⇒ [<code>Apsis</code>](#Apsis)
Finds the date and time of a planet's perihelion (closest approach to the Sun)
or aphelion (farthest distance from the Sun) after a given time.

Given a date and time to start the search in `startTime`, this function finds the
next date and time that the center of the specified planet reaches the closest or farthest point
in its orbit with respect to the center of the Sun, whichever comes first
after `startTime`.

The closest point is called *perihelion* and the farthest point is called *aphelion*.
The word *apsis* refers to either event.

To iterate through consecutive alternating perihelion and aphelion events,
call `SearchPlanetApsis` once, then use the return value to call
[NextPlanetApsis](#NextPlanetApsis). After that, keep feeding the previous return value
from `NextPlanetApsis` into another call of `NextPlanetApsis`
as many times as desired.

**Kind**: global function  
**Returns**: [<code>Apsis</code>](#Apsis) - The next perihelion or aphelion that occurs after `startTime`.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The planet for which to find the next perihelion/aphelion event.      Not allowed to be `"Sun"` or `"Moon"`. |
| startTime | [<code>AstroTime</code>](#AstroTime) | The date and time at which to start searching for the next perihelion or aphelion. |


* * *

<a name="NextPlanetApsis"></a>

## NextPlanetApsis(body, apsis) ⇒ [<code>Apsis</code>](#Apsis)
Finds the next planetary perihelion or aphelion event in a series.

This function requires an [Apsis](#Apsis) value obtained from a call
to [SearchPlanetApsis](#SearchPlanetApsis) or `NextPlanetApsis`.
Given an aphelion event, this function finds the next perihelion event, and vice versa.
See [SearchPlanetApsis](#SearchPlanetApsis) for more details.

**Kind**: global function  
**Returns**: [<code>Apsis</code>](#Apsis) - Same as the return value for [SearchPlanetApsis](#SearchPlanetApsis).  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The planet for which to find the next perihelion/aphelion event.      Not allowed to be `"Sun"` or `"Moon"`.      Must match the body passed into the call that produced the `apsis` parameter. |
| apsis | [<code>Apsis</code>](#Apsis) | An apsis event obtained from a call to [SearchPlanetApsis](#SearchPlanetApsis) or `NextPlanetApsis`. |


* * *

<a name="InverseRotation"></a>

## InverseRotation(rotation) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates the inverse of a rotation matrix.
Given a rotation matrix that performs some coordinate transform,
this function returns the matrix that reverses that trasnform.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - The inverse rotation matrix.  

| Param | Type | Description |
| --- | --- | --- |
| rotation | [<code>RotationMatrix</code>](#RotationMatrix) | The rotation matrix to be inverted. |


* * *

<a name="CombineRotation"></a>

## CombineRotation(a, b) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Creates a rotation based on applying one rotation followed by another.
Given two rotation matrices, returns a combined rotation matrix that is
equivalent to rotating based on the first matrix, followed by the second.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - The combined rotation matrix.  

| Param | Type | Description |
| --- | --- | --- |
| a | [<code>RotationMatrix</code>](#RotationMatrix) | The first rotation to apply. |
| b | [<code>RotationMatrix</code>](#RotationMatrix) | The second rotation to apply. |


* * *

<a name="VectorFromSphere"></a>

## VectorFromSphere(sphere, time) ⇒ [<code>Vector</code>](#Vector)
Converts spherical coordinates to Cartesian coordinates.
Given spherical coordinates and a time at which they are valid,
returns a vector of Cartesian coordinates. The returned value
includes the time, as required by `AstroTime`.

**Kind**: global function  
**Returns**: [<code>Vector</code>](#Vector) - The vector form of the supplied spherical coordinates.  

| Param | Type | Description |
| --- | --- | --- |
| sphere | [<code>Spherical</code>](#Spherical) | Spherical coordinates to be converted. |
| time | [<code>AstroTime</code>](#AstroTime) | The time that should be included in the returned vector. |


* * *

<a name="VectorFromEquator"></a>

## VectorFromEquator(equ, time) ⇒ [<code>Vector</code>](#Vector)
Given angular equatorial coordinates in `equ`, calculates equatorial vector.

**Kind**: global function  
**Returns**: [<code>Vector</code>](#Vector) - A vector in the equatorial system.  

| Param | Type | Description |
| --- | --- | --- |
| equ | [<code>EquatorialCoordinates</code>](#EquatorialCoordinates) | An object that contains angular equatorial coordinates to be converted to a vector. |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the observation. This is needed because the returned      vector object requires a valid time value when passed to certain other functions. |


* * *

<a name="EquatorFromVector"></a>

## EquatorFromVector(vec) ⇒ [<code>EquatorialCoordinates</code>](#EquatorialCoordinates)
Given an equatorial vector, calculates equatorial angular coordinates.

**Kind**: global function  
**Returns**: [<code>EquatorialCoordinates</code>](#EquatorialCoordinates) - Angular coordinates expressed in the same equatorial system as `vec`.  

| Param | Type | Description |
| --- | --- | --- |
| vec | [<code>Vector</code>](#Vector) | A vector in an equatorial coordinate system. |


* * *

<a name="SphereFromVector"></a>

## SphereFromVector(vector) ⇒ [<code>Spherical</code>](#Spherical)
Converts Cartesian coordinates to spherical coordinates.

Given a Cartesian vector, returns latitude, longitude, and distance.

**Kind**: global function  
**Returns**: [<code>Spherical</code>](#Spherical) - Spherical coordinates that are equivalent to the given vector.  

| Param | Type | Description |
| --- | --- | --- |
| vector | [<code>Vector</code>](#Vector) | Cartesian vector to be converted to spherical coordinates. |


* * *

<a name="HorizonFromVector"></a>

## HorizonFromVector(vector, refraction) ⇒ [<code>Spherical</code>](#Spherical)
Converts Cartesian coordinates to horizontal coordinates.

Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.

*IMPORTANT:* This function differs from [SphereFromVector](#SphereFromVector) in two ways:
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

**Kind**: global function  

| Param | Type | Description |
| --- | --- | --- |
| vector | [<code>Vector</code>](#Vector) | Cartesian vector to be converted to horizontal coordinates. |
| refraction | <code>string</code> | `"normal"`: correct altitude for atmospheric refraction (recommended).      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.      `null`: no atmospheric refraction correction is performed. |


* * *

<a name="VectorFromHorizon"></a>

## VectorFromHorizon(sphere, time, refraction) ⇒ [<code>Vector</code>](#Vector)
Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.

**Kind**: global function  
**Returns**: [<code>Vector</code>](#Vector) - A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).  

| Param | Type | Description |
| --- | --- | --- |
| sphere | [<code>Spherical</code>](#Spherical) | A structure that contains apparent horizontal coordinates:      `lat` holds the refracted azimuth angle,      `lon` holds the azimuth in degrees clockwise from north,      and `dist` holds the distance from the observer to the object in AU. |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the observation. This is needed because the returned      vector object requires a valid time value when passed to certain other functions. |
| refraction | <code>string</code> | `"normal"`: correct altitude for atmospheric refraction (recommended).      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.      `null`: no atmospheric refraction correction is performed. |


* * *

<a name="Refraction"></a>

## Refraction(refraction, altitude) ⇒ <code>number</code>
Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.

Given an altitude angle and a refraction option, calculates
the amount of "lift" caused by atmospheric refraction.
This is the number of degrees higher in the sky an object appears
due to the lensing of the Earth's atmosphere.

**Kind**: global function  
**Returns**: <code>number</code> - The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing.  

| Param | Type | Description |
| --- | --- | --- |
| refraction | <code>string</code> | `"normal"`: correct altitude for atmospheric refraction (recommended).      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.      `null`: no atmospheric refraction correction is performed. |
| altitude | <code>number</code> | An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90. |


* * *

<a name="InverseRefraction"></a>

## InverseRefraction(refraction, bent_altitude) ⇒ <code>number</code>
Calculates the inverse of an atmospheric refraction angle.

Given an observed altitude angle that includes atmospheric refraction,
calculate the negative angular correction to obtain the unrefracted
altitude. This is useful for cases where observed horizontal
coordinates are to be converted to another orientation system,
but refraction first must be removed from the observed position.

**Kind**: global function  
**Returns**: <code>number</code> - The angular adjustment in degrees to be added to the
     altitude angle to correct for atmospheric lensing.
     This will be less than or equal to zero.  

| Param | Type | Description |
| --- | --- | --- |
| refraction | <code>string</code> | `"normal"`: correct altitude for atmospheric refraction (recommended).      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.      `null`: no atmospheric refraction correction is performed. |
| bent_altitude | <code>number</code> | The apparent altitude that includes atmospheric refraction. |


* * *

<a name="RotateVector"></a>

## RotateVector(rotation, vector) ⇒ [<code>Vector</code>](#Vector)
Applies a rotation to a vector, yielding a rotated vector.

This function transforms a vector in one orientation to a vector
in another orientation.

**Kind**: global function  
**Returns**: [<code>Vector</code>](#Vector) - A vector in the orientation specified by `rotation`.  

| Param | Type | Description |
| --- | --- | --- |
| rotation | [<code>RotationMatrix</code>](#RotationMatrix) | A rotation matrix that specifies how the orientation of the vector is to be changed. |
| vector | [<code>Vector</code>](#Vector) | The vector whose orientation is to be changed. |


* * *

<a name="Rotation_EQJ_ECL"></a>

## Rotation\_EQJ\_ECL() ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: ECL = ecliptic system, using equator at J2000 epoch.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts EQJ to ECL.  

* * *

<a name="Rotation_ECL_EQJ"></a>

## Rotation\_ECL\_EQJ() ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQJ = equatorial system, using equator at J2000 epoch.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts ECL to EQJ.  

* * *

<a name="Rotation_EQJ_EQD"></a>

## Rotation\_EQJ\_EQD(time) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of the specified date/time.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts EQJ to EQD at `time`.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time at which the Earth's equator defines the target orientation. |


* * *

<a name="Rotation_EQD_EQJ"></a>

## Rotation\_EQD\_EQJ(time) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: EQJ = equatorial system, using equator at J2000 epoch.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts EQD at `time` to EQJ.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time at which the Earth's equator defines the source orientation. |


* * *

<a name="Rotation_EQD_HOR"></a>

## Rotation\_EQD\_HOR(time, observer) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: HOR = horizontal system.

Use `HorizonFromVector` to convert the return value
to a traditional altitude/azimuth pair.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts EQD to HOR at `time` and for `observer`.
     The components of the horizontal vector are:
     x = north, y = west, z = zenith (straight up from the observer).
     These components are chosen so that the "right-hand rule" works for the vector
     and so that north represents the direction where azimuth = 0.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time at which the Earth's equator applies. |
| observer | [<code>Observer</code>](#Observer) | A location near the Earth's mean sea level that defines the observer's horizon. |


* * *

<a name="Rotation_HOR_EQD"></a>

## Rotation\_HOR\_EQD(time, observer) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQD = equatorial system, using equator of the specified date/time.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts HOR to EQD at `time` and for `observer`.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time at which the Earth's equator applies. |
| observer | [<code>Observer</code>](#Observer) | A location near the Earth's mean sea level that defines the observer's horizon. |


* * *

<a name="Rotation_HOR_EQJ"></a>

## Rotation\_HOR\_EQJ(time, observer) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQJ = equatorial system, using equator at the J2000 epoch.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts HOR to EQD at `time` and for `observer`.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the observation. |
| observer | [<code>Observer</code>](#Observer) | A location near the Earth's mean sea level that defines the observer's horizon. |


* * *

<a name="Rotation_EQJ_HOR"></a>

## Rotation\_EQJ\_HOR(time, observer) ⇒
Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: HOR = horizontal system.

Use [HorizonFromVector](#HorizonFromVector) to convert the return value
to a traditional altitude/azimuth pair.

**Kind**: global function  
**Returns**: A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
     The components of the horizontal vector are:
     x = north, y = west, z = zenith (straight up from the observer).
     These components are chosen so that the "right-hand rule" works for the vector
     and so that north represents the direction where azimuth = 0.  

| Param | Description |
| --- | --- |
| time | The date and time of the desired horizontal orientation. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |


* * *

<a name="Rotation_EQD_ECL"></a>

## Rotation\_EQD\_ECL(time) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of date.
Target: ECL = ecliptic system, using equator at J2000 epoch.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts EQD to ECL.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the source equator. |


* * *

<a name="Rotation_ECL_EQD"></a>

## Rotation\_ECL\_EQD(time) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of date.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts ECL to EQD.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the desired equator. |


* * *

<a name="Rotation_ECL_HOR"></a>

## Rotation\_ECL\_HOR(time, observer) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: HOR = horizontal system.

Use [HorizonFromVector](#HorizonFromVector) to convert the return value
to a traditional altitude/azimuth pair.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts ECL to HOR at `time` and for `observer`.
     The components of the horizontal vector are:
     x = north, y = west, z = zenith (straight up from the observer).
     These components are chosen so that the "right-hand rule" works for the vector
     and so that north represents the direction where azimuth = 0.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the desired horizontal orientation. |
| observer | [<code>Observer</code>](#Observer) | A location near the Earth's mean sea level that defines the observer's horizon. |


* * *

<a name="Rotation_HOR_ECL"></a>

## Rotation\_HOR\_ECL(time, observer) ⇒ [<code>RotationMatrix</code>](#RotationMatrix)
Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system.
Target: ECL = ecliptic system, using equator at J2000 epoch.

**Kind**: global function  
**Returns**: [<code>RotationMatrix</code>](#RotationMatrix) - A rotation matrix that converts HOR to ECL.  

| Param | Type | Description |
| --- | --- | --- |
| time | [<code>AstroTime</code>](#AstroTime) | The date and time of the horizontal observation. |
| observer | [<code>Observer</code>](#Observer) | The location of the horizontal observer. |


* * *

<a name="Constellation"></a>

## Constellation(ra, dec) ⇒ [<code>ConstellationInfo</code>](#ConstellationInfo)
Determines the constellation that contains the given point in the sky.

Given J2000 equatorial (EQJ) coordinates of a point in the sky,
determines the constellation that contains that point.

**Kind**: global function  
**Returns**: [<code>ConstellationInfo</code>](#ConstellationInfo) - An object that contains the 3-letter abbreviation and full name
     of the constellation that contains the given (ra,dec), along with
     the converted B1875 (ra,dec) for that point.  

| Param | Type | Description |
| --- | --- | --- |
| ra | <code>number</code> | The right ascension (RA) of a point in the sky, using the J2000 equatorial system. |
| dec | <code>number</code> | The declination (DEC) of a point in the sky, using the J2000 equatorial system. |


* * *

<a name="SearchLunarEclipse"></a>

## SearchLunarEclipse(date) ⇒ [<code>LunarEclipseInfo</code>](#LunarEclipseInfo)
**Kind**: global function  
**Brief**: Searches for a lunar eclipse.

This function finds the first lunar eclipse that occurs after `startTime`.
A lunar eclipse may be penumbral, partial, or total.
See [LunarEclipseInfo](#LunarEclipseInfo) for more information.
To find a series of lunar eclipses, call this function once,
then keep calling [NextLunarEclipse](#NextLunarEclipse) as many times as desired,
passing in the `center` value returned from the previous call.  

| Param | Type | Description |
| --- | --- | --- |
| date | <code>Date</code> \| <code>number</code> \| [<code>AstroTime</code>](#AstroTime) | The date and time for starting the search for a lunar eclipse. |


* * *

<a name="NextLunarEclipse"></a>

## NextLunarEclipse(prevEclipseTime) ⇒ [<code>LunarEclipseInfo</code>](#LunarEclipseInfo)
**Kind**: global function  
**Brief**: Searches for the next lunar eclipse in a series.

After using [SearchLunarEclipse](#SearchLunarEclipse) to find the first lunar eclipse
in a series, you can call this function to find the next consecutive lunar eclipse.
Pass in the `center` value from the [LunarEclipseInfo](#LunarEclipseInfo) returned by the
previous call to `SearchLunarEclipse` or `NextLunarEclipse`
to find the next lunar eclipse.  

| Param | Type | Description |
| --- | --- | --- |
| prevEclipseTime | [<code>AstroTime</code>](#AstroTime) | A date and time near a full moon. Lunar eclipse search will start at the next full moon. |


* * *

<a name="SearchGlobalSolarEclipse"></a>

## SearchGlobalSolarEclipse(startTime) ⇒ [<code>GlobalSolarEclipseInfo</code>](#GlobalSolarEclipseInfo)
**Kind**: global function  
**Brief**: Searches for a solar eclipse visible anywhere on the Earth's surface.

This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [GlobalSolarEclipseInfo](#GlobalSolarEclipseInfo) for more information.
To find a series of solar eclipses, call this function once,
then keep calling [NextGlobalSolarEclipse](#NextGlobalSolarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.  

| Param | Type | Description |
| --- | --- | --- |
| startTime | [<code>AstroTime</code>](#AstroTime) | The date and time for starting the search for a solar eclipse. |


* * *

<a name="NextGlobalSolarEclipse"></a>

## NextGlobalSolarEclipse(prevEclipseTime) ⇒ [<code>GlobalSolarEclipseInfo</code>](#GlobalSolarEclipseInfo)
**Kind**: global function  
**Brief**: Searches for the next global solar eclipse in a series.

After using [SearchGlobalSolarEclipse](#SearchGlobalSolarEclipse) to find the first solar eclipse
in a series, you can call this function to find the next consecutive solar eclipse.
Pass in the `peak` value from the [GlobalSolarEclipseInfo](#GlobalSolarEclipseInfo) returned by the
previous call to `SearchGlobalSolarEclipse` or `NextGlobalSolarEclipse`
to find the next solar eclipse.  

| Param | Type | Description |
| --- | --- | --- |
| prevEclipseTime | [<code>AstroTime</code>](#AstroTime) | A date and time near a new moon. Solar eclipse search will start at the next new moon. |


* * *

<a name="SearchLocalSolarEclipse"></a>

## SearchLocalSolarEclipse(startTime, observer) ⇒ [<code>LocalSolarEclipseInfo</code>](#LocalSolarEclipseInfo)
**Kind**: global function  
**Brief**: Searches for a solar eclipse visible at a specific location on the Earth's surface.

This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [LocalSolarEclipseInfo](#LocalSolarEclipseInfo) for more information.

To find a series of solar eclipses, call this function once,
then keep calling [NextLocalSolarEclipse](#NextLocalSolarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.

IMPORTANT: An eclipse reported by this function might be partly or
completely invisible to the observer due to the time of day.
See [LocalSolarEclipseInfo](#LocalSolarEclipseInfo) for more information about this topic.  

| Param | Type | Description |
| --- | --- | --- |
| startTime | [<code>AstroTime</code>](#AstroTime) | The date and time for starting the search for a solar eclipse. |
| observer | [<code>Observer</code>](#Observer) | The geographic location of the observer. |


* * *

<a name="NextLocalSolarEclipse"></a>

## NextLocalSolarEclipse(prevEclipseTime, observer) ⇒ [<code>LocalSolarEclipseInfo</code>](#LocalSolarEclipseInfo)
**Kind**: global function  
**Brief**: Searches for the next local solar eclipse in a series.

After using [SearchLocalSolarEclipse](#SearchLocalSolarEclipse) to find the first solar eclipse
in a series, you can call this function to find the next consecutive solar eclipse.
Pass in the `peak` value from the [LocalSolarEclipseInfo](#LocalSolarEclipseInfo) returned by the
previous call to `SearchLocalSolarEclipse` or `NextLocalSolarEclipse`
to find the next solar eclipse.
This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [LocalSolarEclipseInfo](#LocalSolarEclipseInfo) for more information.  

| Param | Type | Description |
| --- | --- | --- |
| prevEclipseTime | [<code>AstroTime</code>](#AstroTime) | The date and time for starting the search for a solar eclipse. |
| observer | [<code>Observer</code>](#Observer) | The geographic location of the observer. |


* * *

<a name="SearchTransit"></a>

## SearchTransit(body, startTime) ⇒ [<code>TransitInfo</code>](#TransitInfo)
**Kind**: global function  
**Brief**: Searches for the first transit of Mercury or Venus after a given date.

Finds the first transit of Mercury or Venus after a specified date.
A transit is when an inferior planet passes between the Sun and the Earth
so that the silhouette of the planet is visible against the Sun in the background.
To continue the search, pass the `finish` time in the returned structure to
[NextTransit](#NextTransit).  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The planet whose transit is to be found. Must be `"Mercury"` or `"Venus"`. |
| startTime | [<code>AstroTime</code>](#AstroTime) | The date and time for starting the search for a transit. |


* * *

<a name="NextTransit"></a>

## NextTransit(body, prevTransitTime) ⇒ [<code>TransitInfo</code>](#TransitInfo)
**Kind**: global function  
**Brief**: Searches for another transit of Mercury or Venus.

After calling [SearchTransit](#SearchTransit) to find a transit of Mercury or Venus,
this function finds the next transit after that.
Keep calling this function as many times as you want to keep finding more transits.  

| Param | Type | Description |
| --- | --- | --- |
| body | <code>string</code> | The planet whose transit is to be found. Must be `"Mercury"` or `"Venus"`. |
| prevTransitTime | [<code>AstroTime</code>](#AstroTime) | A date and time near the previous transit. |


* * *

<a name="SearchOptions"></a>

## SearchOptions : <code>Object</code>
Options for the [Search](#Search) function.

**Kind**: global typedef  
**Properties**

| Name | Type | Description |
| --- | --- | --- |
| dt_tolerance_seconds | <code>number</code> \| <code>null</code> | The number of seconds for a time window smaller than which the search      is considered successful.  Using too large a tolerance can result in      an inaccurate time estimate.  Using too small a tolerance can cause      excessive computation, or can even cause the search to fail because of      limited floating-point resolution.  Defaults to 1 second. |
| init_f1 | <code>number</code> \| <code>null</code> | As an optimization, if the caller of [Search](#Search)      has already calculated the value of the function being searched (the parameter `func`)      at the time coordinate `t1`, it can pass in that value as `init_f1`.      For very expensive calculations, this can measurably improve performance. |
| init_f2 | <code>number</code> \| <code>null</code> | The same as `init_f1`, except this is the optional initial value of `func(t2)`      instead of `func(t1)`. |


* * *

