unit Astronomy;

{$mode objfpc}{$H+}

interface

uses
  Classes
  { you can add units after this };

const
  LibraryName = 'astronomy'; // should refer the astronomy.dll file

  C_AUDAY : Double                      = 173.1446326846693;
  KM_PER_AU : Double                    = 1.4959787069098932e+8;
  AU_PER_LY : Double                    = 63241.07708807546;
  DEG2RAD : Double                      = 0.017453292519943296;
  HOUR2RAD : Double                     = 0.2617993877991494365;
  RAD2DEG : Double                      = 57.295779513082321;
  RAD2HOUR : Double                     = 3.819718634205488;
  SUN_RADIUS_KM : Double                = 695700.0;
  MERCURY_EQUATORIAL_RADIUS_KM : Double = 2440.5;
  MERCURY_POLAR_RADIUS_KM : Double      = 2438.3;
  VENUS_RADIUS_KM : Double              = 6051.8;
  EARTH_EQUATORIAL_RADIUS_KM : Double   = 6378.1366;
  EARTH_FLATTENING : Double             = 0.996647180302104;
  EARTH_POLAR_RADIUS_KM : Double        = 6378.1366 * 0.996647180302104;
  MOON_EQUATORIAL_RADIUS_KM : Double    = 1738.1;
  MOON_POLAR_RADIUS_KM : Double         = 1736.0;
  MARS_EQUATORIAL_RADIUS_KM : Double    = 3396.2;
  MARS_POLAR_RADIUS_KM : Double         = 3376.2;
  JUPITER_EQUATORIAL_RADIUS_KM : Double = 71492.0;
  JUPITER_POLAR_RADIUS_KM : Double      = 66854.0;
  JUPITER_MEAN_RADIUS_KM : Double       = 69911.0;
  IO_RADIUS_KM : Double                 = 1821.6;
  EUROPA_RADIUS_KM : Double             = 1560.8;
  GANYMEDE_RADIUS_KM : Double           = 2631.2;
  CALLISTO_RADIUS_KM : Double           = 2410.3;
  SATURN_EQUATORIAL_RADIUS_KM : Double  = 60268.0;
  SATURN_POLAR_RADIUS_KM : Double       = 54364.0;
  URANUS_EQUATORIAL_RADIUS_KM : Double  = 25559.0;
  URANUS_POLAR_RADIUS_KM : Double       = 24973.0;
  NEPTUNE_EQUATORIAL_RADIUS_KM : Double = 24764.0;
  NEPTUNE_POLAR_RADIUS_KM : Double      = 24341.0;
  PLUTO_RADIUS_KM : Double              = 1188.3;

type
  // Status of the last operation.
  AstroStatus = (
    ASTRO_SUCCESS            = 0,  // The operation was successful.
    ASTRO_NOT_INITIALIZED    = 1,  // A placeholder that can be used for data that is not yet initialized.
    ASTRO_INVALID_BODY       = 2,  // The celestial body was not valid. Different sets of bodies are supported depending on the function.
    ASTRO_NO_CONVERGE        = 3,  // A numeric solver failed to converge. This should not happen unless there is a bug in Astronomy Engine.
    ASTRO_BAD_TIME           = 4,  // The provided date/time is outside the range allowed by this function.
    ASTRO_BAD_VECTOR         = 5,  // Vector magnitude is too small to be normalized into a unit vector.
    ASTRO_SEARCH_FAILURE     = 6,  // Search was not able to find an ascending root crossing of the function in the specified time interval.
    ASTRO_EARTH_NOT_ALLOWED  = 7,  // The Earth cannot be treated as a celestial body seen from an observer on the Earth itself.
    ASTRO_NO_MOON_QUARTER    = 8,  // No lunar quarter occurs inside the specified time range.
    ASTRO_WRONG_MOON_QUARTER = 9,  // Internal error: Astronomy_NextMoonQuarter found the wrong moon quarter.
    ASTRO_INTERNAL_ERROR     = 10, // A self-check failed inside the code somewhere, indicating a bug needs to be fixed.
    ASTRO_INVALID_PARAMETER  = 11, // A parameter value passed to a function was not valid.
    ASTRO_FAIL_APSIS         = 12, // Special-case logic for finding Neptune/Pluto apsis failed.
    ASTRO_BUFFER_TOO_SMALL   = 13, // A provided buffer's size is too small to receive the requested data.
    ASTRO_OUT_OF_MEMORY      = 14, // An attempt to allocate memory failed.
    ASTRO_INCONSISTENT_TIMES = 15  // The provided initial state vectors did not have matching times.
  );

  // A celestial body.
  AstroBody = (
    BODY_INVALID = -1,  // An invalid or undefined celestial body.
    BODY_MERCURY = 0,   // Mercury
    BODY_VENUS   = 1,   // Venus
    BODY_EARTH   = 2,   // Earth
    BODY_MARS    = 3,   // Mars
    BODY_JUPITER = 4,   // Jupiter
    BODY_SATURN  = 5,   // Saturn
    BODY_URANUS  = 6,   // Uranus
    BODY_NEPTUNE = 7,   // Neptune
    BODY_PLUTO   = 8,   // Pluto
    BODY_SUN     = 9,   // Sun
    BODY_MOON    = 10,  // Moon
    BODY_EMB     = 11,  // Earth/Moon Barycenter
    BODY_SSB     = 12,  // Solar System Barycenter
    BODY_STAR1   = 101, // user-defined star #1
    BODY_STAR2   = 102, // user-defined star #2
    BODY_STAR3   = 103, // user-defined star #3
    BODY_STAR4   = 104, // user-defined star #4
    BODY_STAR5   = 105, // user-defined star #5
    BODY_STAR6   = 106, // user-defined star #6
    BODY_STAR7   = 107, // user-defined star #7
    BODY_STAR8   = 108  // user-defined star #8
  );

  // Selects whether to correct for atmospheric refraction, and if so, how.
  AstroRefraction = (
    REFRACTION_NONE = 0,   // No atmospheric refraction correction (airless).
    REFRACTION_NORMAL = 1, // Recommended correction for standard atmospheric refraction.
    REFRACTION_JPLHOR = 2  // Used only for compatibility testing with JPL Horizons online tool.
  );

  // Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.
  AstroVisibility = (
    VISIBLE_MORNING = 0,    // The body is best visible in the morning, before sunrise.
    VISIBLE_EVENING = 1     // The body is best visible in the evening, after sunset.
  );

  // The type of apsis: pericenter (closest approach) or apocenter (farthest distance).
  AstroApsisKind = (
    APSIS_PERICENTER = 0,   // The body is at its closest approach to the object it orbits.
    APSIS_APOCENTER = 1,    // The body is at its farthest distance from the object it orbits.
    APSIS_INVALID = 2       // Undefined or invalid apsis.
  );

  // The different kinds of lunar/solar eclipses.
  AstroEclipseKind = (
    ECLIPSE_NONE = 0,       // No eclipse found.
    ECLIPSE_PENUMBRAL = 1,  // A penumbral lunar eclipse. (Never used for a solar eclipse.)
    ECLIPSE_PARTIAL = 2,    // A partial lunar/solar eclipse.
    ECLIPSE_ANNULAR = 3,    // An annular solar eclipse. (Never used for a lunar eclipse.)
    ECLIPSE_TOTAL = 4       // A total lunar/solar eclipse.
  );

  // Aberration calculation options.
  AstroAberration = (
    ABERRATION = 0,     // Request correction for aberration.
    NO_ABERRATION = 1   // Do not correct for aberration.
  );

  // Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.
  AstroEquatorDate = (
    EQUATOR_J2000 = 0,      // Represent equatorial coordinates in the J2000 epoch.
    EQUATOR_OF_DATE = 1     // Represent equatorial coordinates using the Earth's equator at the given date and time.
  );

  // Selects whether to search for a rise time or a set time.
  AstroDirection = (
    DIRECTION_SET  = -1,  // Search for the time a body finishes sinking below the horizon.
    DIRECTION_RISE = 1    // Search for the time a body begins to rise above the horizon.
  );

  // Selects the output format of the function #Astronomy_FormatTime.
  AstroTimeFormat = (
    TIME_FORMAT_DAY = 0,    // Truncate to UTC calendar date only, e.g. `2020-12-31`. Buffer size must be at least 11 characters.
    TIME_FORMAT_MINUTE = 1, // Round to nearest UTC minute, e.g. `2020-12-31T15:47Z`. Buffer size must be at least 18 characters.
    TIME_FORMAT_SECOND = 2, // Round to nearest UTC second, e.g. `2020-12-31T15:47:32Z`. Buffer size must be at least 21 characters.
    TIME_FORMAT_MILLI = 3   // Round to nearest UTC millisecond, e.g. `2020-12-31T15:47:32.397Z`. Buffer size must be at least 25 characters.
  );

  // Indicates whether a crossing through the ecliptic plane is ascending or descending.
  AstroNodeKind = (
    DESCENDING_NODE = -1, // The body passes through the ecliptic plane from north to south.
    INVALID_NODE    = 0,  // Placeholder value for a missing or invalid node.
    ASCENDING_NODE  = 1   // The body passes through the ecliptic plane from south to north.
  );

  // A date and time used for astronomical calculations.
  AstroTime = packed record
    UT: Double;  // UT1/UTC number of days since noon on January 1, 2000.
    TT: Double;  // Terrestrial Time days since noon on January 1, 2000.
    Psi: Double; // For internal use only. Used to optimize Earth tilt calculations.
    Eps: Double; // For internal use only.  Used to optimize Earth tilt calculations.
    ST: Double;  // For internal use only.  Lazy-caches sidereal time (Earth rotation).
  end;

  // Pointer to #AstroTime record.
  PAstroTime = ^AstroTime;

  // A calendar date and time expressed in UTC.
  AstroUTC = packed record
    Year: Int32;    // The year value, e.g. 2019.
    Month: Int32;   // The month value: 1=January, 2=February, ..., 12=December.
    Day: Int32;     // The day of the month in the range 1..31.
    Hour: Int32;    // The hour of the day in the range 0..23.
    Minute: Int32;  // The minute of the hour in the range 0..59.
    Second: Double; // The floating point number of seconds in the range [0,60).
  end;

  // A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
  AstroVector = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    X: Double;           // The Cartesian x-coordinate of the vector in AU.
    Y: Double;           // The Cartesian y-coordinate of the vector in AU.
    Z: Double;           // The Cartesian z-coordinate of the vector in AU.
    T: AstroTime;        // The date and time at which this vector is valid.
  end;

  // Pointer to #AstroVector record
  PAstroVector = ^AstroVector;

  // A state vector that contains a position (AU) and velocity (AU/day).
  AstroStateVector = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    X: Double;           // The Cartesian position x-coordinate of the vector in AU.
    Y: Double;           // The Cartesian position y-coordinate of the vector in AU.
    Z: Double;           // The Cartesian position z-coordinate of the vector in AU.
    VX: Double;          // The Cartesian velocity x-coordinate of the vector in AU/day.
    VY: Double;          // The Cartesian velocity y-coordinate of the vector in AU/day.
    VZ: Double;          // The Cartesian velocity z-coordinate of the vector in AU/day.
    T: AstroTime;        // The date and time at which this state vector is valid.
  end;

  // Pointer to `AstroStateVector` record.
  PAstroStateVector = ^AstroStateVector;

  // Spherical coordinates: latitude, longitude, distance.
  AstroSpherical = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Lat: Double;         // The latitude angle: -90..+90 degrees.
    Lon: Double;         // The longitude angle: 0..360 degrees.
    Dist: Double;        // Distance in AU.
  end;

  // An angular value expressed in degrees.
  AstroAngleResult = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Angle: Double;       // An angle expressed in degrees.
  end;

  // The location of an observer on (or near) the surface of the Earth.
  AstroObserver = packed record
    Latitude: Double;  // Geographic latitude in degrees north (positive) or south (negative) of the equator.
    Longitude: Double; // Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.
    Height: Double;    // The height above (positive) or below (negative) sea level, expressed in meters.
  end;

  // Equatorial angular and cartesian coordinates.
  AstroEquatorial = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    RA: Double;          // right ascension in sidereal hours.
    Dec: Double;         // declination in degrees
    Dist: Double;        // distance to the celestial body in AU.
    Vec: AstroVector;    // equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north.
  end;

  // Ecliptic angular and Cartesian coordinates.
  AstroEcliptic = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Vec: AstroVector;    // Cartesian ecliptic vector: x=equinox, y=90 degrees prograde in ecliptic plane, z=northward perpendicular to ecliptic.
    ELat: Double;        // Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
    ELon: Double;        // Longitude in degrees around the ecliptic plane prograde from the equinox.
  end;

  // Coordinates of a celestial body as seen by a topocentric observer.
  AstroHorizon = packed record
    Azimuth: Double;  // Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West.
    Altitude: Double; // Angle in degrees above (positive) or below (negative) the observer's horizon.
    RA: Double;       // Right ascension in sidereal hours.
    Dec: Double;      // Declination in degrees.
  end;

  // Contains a rotation matrix that can be used to transform one coordinate system to another.
  AstroRotation = packed record
    Status: AstroStatus;              // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Rot: array[1..3, 1..3] of Double; // A normalized 3x3 rotation matrix.
  end;

  // Information about idealized atmospheric variables at a given elevation.
  AstroAtmosphere = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Pressure: Double;    // Atmospheric pressure in pascals
    Temperature: Double; // Atmospheric temperature in kelvins
    Density: Double;     // Atmospheric density relative to sea level
  end;

  // The result of a search for an astronomical event.
  AstroSearchResult = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Time: AstroTime;     // The time at which a searched-for event occurs.
  end;

{ The dates and times of changes of season for a given calendar year.
  Call #Astronomy_Seasons to calculate this data structure for a given year. }
  AstroSeasons = packed record
    Status: AstroStatus;    // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    MarEquinox: AstroTime;  // The date and time of the March equinox for the specified year.
    JunSolstice: AstroTime; // The date and time of the June soltice for the specified year.
    SepEquinox: AstroTime;  // The date and time of the September equinox for the specified year.
    DecSolstice: AstroTime; // The date and time of the December solstice for the specified year.
  end;

  // A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.
  AstroMoonQuarter = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Quarter: Int32;      // 0=new moon, 1=first quarter, 2=full moon, 3=third quarter.
    Time: AstroTime;     // The date and time of the lunar quarter.
  end;

  // A real value returned by a function whose ascending root is to be found.
  AstroFuncResult = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Value: Double;       // The value returned by a function whose ascending root is to be found.
  end;

  // Contains information about the visibility of a celestial body at a given date and time.
  AstroElongation = packed record
    Status: AstroStatus;         // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Time: AstroTime;             // The date and time of the observation.
    Visibility: AstroVisibility; // Whether the body is best seen in the morning or the evening.
    Elongation: Double;          // The angle in degrees between the body and the Sun, as seen from the Earth.
    EclipticSeparation: Double;  // The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.
  end;

  // Information about a celestial body crossing a specific hour angle.
  AstroHourAngle = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Time: AstroTime;     // The date and time when the body crosses the specified hour angle.
    Hor: AstroHorizon;   // Apparent coordinates of the body at the time it crosses the specified hour angle.
  end;

  // Information about the brightness and illuminated shape of a celestial body.
  AstroIllum = packed record
    Status: AstroStatus;   // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Time: AstroTime;       // The date and time of the observation.
    Mag: Double;           // The visual magnitude of the body. Smaller values are brighter.
    PhaseAngle: Double;    // The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth.
    PhaseFraction: Double; // A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth.
    HelioDist: Double;     // The distance between the Sun and the body at the observation time.
    RingTilt: Double;      // For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0.
  end;

  // An apsis event: pericenter (closest approach) or apocenter (farthest distance).
  AstroApsis = packed record
    Status: AstroStatus;  // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Time: AstroTime;      // The date and time of the apsis.
    Kind: AstroApsisKind; // Whether this is a pericenter or apocenter event.
    DistAU: Double;       // The distance between the centers of the bodies in astronomical units.
    DistKM: Double;       // The distance between the centers of the bodies in kilometers.
  end;

  // Information about a lunar eclipse.
  AstroLunarEclipse = packed record
    Status: AstroStatus;    // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Kind: AstroEclipseKind; // The type of lunar eclipse found.
    Obscuration: Double;    // The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra.
    Peak: AstroTime;        // The time of the eclipse at its peak.
    SDPenum: Double;        // The semi-duration of the penumbral phase in minutes.
    SDPartial: Double;      // The semi-duration of the partial phase in minutes, or 0.0 if none.
    SDTotal: Double;        // The semi-duration of the total phase in minutes, or 0.0 if none.
  end;

  // Reports the time and geographic location of the peak of a solar eclipse.
  AstroGlobalSolarEclipse = packed record
    Status: AstroStatus;    // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Kind: AstroEclipseKind; // The type of solar eclipse found.
    Obscuration: Double;    // The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only).
    Peak: AstroTime;        // The date and time when the solar eclipse is darkest. This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.
    Distance: Double;       // The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.
    Latitude: Double;       // The geographic latitude at the center of the peak eclipse shadow.
    Longitude: Double;      // The geographic longitude at the center of the peak eclipse shadow.
  end;

  // Holds a time and the observed altitude of the Sun at that time.
  AstroEclipseEvent = packed record
    Time: AstroTime;  // The date and time of the event.
    Altitude: Double; // The angular altitude of the center of the Sun above/below the horizon, at `time`, corrected for atmospheric refraction and expressed in degrees.
  end;

  // Information about a solar eclipse as seen by an observer at a given time and geographic location.
  AstroLocalSolarEclipse = packed record
    Status: AstroStatus;             // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Kind: AstroEclipseKind;          // The type of solar eclipse found: `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`.
    Obscuration: Double;             // The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak.
    PartialBegin: AstroEclipseEvent; // The time and Sun altitude at the beginning of the eclipse.
    TotalBegin: AstroEclipseEvent;   // If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise invalid.
    Peak: AstroEclipseEvent;         // The time and Sun altitude when the eclipse reaches its peak.
    TotalEnd: AstroEclipseEvent;     // If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise invalid.
    PartialEnd: AstroEclipseEvent;   // The time and Sun altitude at the end of the eclipse.
  end;

  // Information about a transit of Mercury or Venus, as seen from the Earth.
  AstroTransit = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Start: AstroTime;    // Date and time at the beginning of the transit.
    Peak: AstroTime;     // Date and time of the peak of the transit.
    Finish: AstroTime;   // Date and time at the end of the transit.
    Separation: Double;  // Angular separation in arcminutes between the centers of the Sun and the planet at time `peak`.
  end;

  // Reports the constellation that a given celestial point lies within.
  AstroConstellation = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Symbol: PChar;       // 3-character mnemonic symbol for the constellation, e.g. "Ori".
    Name: PChar;         // Full name of constellation, e.g. "Orion".
    RA1875: Double;      // Right ascension expressed in B1875 coordinates.
    Dec1875: Double;     // Declination expressed in B1875 coordinates.
  end;

  // Lunar libration angles, returned by #Astronomy_Libration.
  AstroLibration = packed record
    ELat: Double;    // Sub-Earth libration ecliptic latitude angle, in degrees.
    ELon: Double;    // Sub-Earth libration ecliptic longitude angle, in degrees.
    MLat: Double;    // Moon's geocentric ecliptic latitude, in degrees.
    MLon: Double;    // Moon's geocentric ecliptic longitude, in degrees.
    DistKM: Double;  // Distance between the centers of the Earth and Moon in kilometers.
    DiamDeg: Double; // The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth.
  end;

  // Information about a body's rotation axis at a given time.
  AstroAxis = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    RA: Double;          // The J2000 right ascension of the body's north pole direction, in sidereal hours.
    Dec: Double;         // The J2000 declination of the body's north pole direction, in degrees.
    Spin: Double;        // Rotation angle of the body's prime meridian, in degrees.
    North: AstroVector;  // A J2000 dimensionless unit vector pointing in the direction of the body's north pole.
  end;

  // Holds the positions and velocities of Jupiter's major 4 moons.
  AstroJupiterMoons = packed record
    IO: AstroStateVector;       // Jovicentric position and velocity of Io.
    Europa: AstroStateVector;   // Jovicentric position and velocity of Europa.
    Ganymede: AstroStateVector; // Jovicentric position and velocity of Ganymede.
    Callisto: AstroStateVector; // Jovicentric position and velocity of Callisto.
  end;

  // Information about an ascending or descending node of a body.
  AstroNodeEvent = packed record
    Status: AstroStatus; // `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.
    Time: AstroTime;     // The time when the body passes through the ecliptic plane.
    Kind: AstroNodeKind; // Either `ASCENDING_NODE` or `DESCENDING_NODE`, depending on the direction of the ecliptic plane crossing.
  end;

  // A data type used for managing simulation of the gravitational forces on a small body.
  AstroGravSim = packed record
  end;

  // Pointer to `AstroGravSim` record.
  PAstroGravSim = ^AstroGravSim;

  // A pointer to a function that is to be passed as a callback to #Astronomy_Search.
  PAstroSearchFunc = function(Context: Pointer; Time: AstroTime): AstroFuncResult; cdecl;

  // A pointer to a function that calculates Delta T.
  PAstroDeltaTFunc = function(UT: Double): Double; cdecl;

  // A pointer to a function for which to solve a light-travel time problem.
  PAstroPositionFunc = function(Context: Pointer; Time: AstroTime): AstroVector; cdecl;

  // The default Delta T function used by Astronomy Engine.
  function Astronomy_DeltaT_EspenakMeeus(UT: Double): Double; cdecl; external LibraryName;

  // A Delta T function that approximates the one used by the JPL Horizons tool
  function Astronomy_DeltaT_JplHorizons(UT: Double): Double; cdecl; external LibraryName;

  // Changes the function Astronomy Engine uses to calculate Delta T.
  procedure Astronomy_SetDeltaTFunction(Func: PAstroDeltaTFunc); cdecl; external LibraryName;

  // Frees up all dynamic memory allocated by Astronomy Engine.
  procedure Astronomy_Reset(); cdecl; external LibraryName;

  // Calculates the length of the given vector. The returned length is expressed usually in AU.
  function Astronomy_VectorLength(Vector: AstroVector): Double; cdecl; external LibraryName;

  // Calculates the angle between two vectors. The returned angle is expressed in degrees.
  function Astronomy_AngleBetween(A: AstroVector; B: AstroVector): AstroAngleResult; cdecl; external LibraryName;

  // Finds the name of a celestial body.
  function Astronomy_BodyName(Body: AstroBody): PChar; cdecl; external LibraryName;

{ Returns the #AstroBody value corresponding to the given English name.
  Valid parameter values are one of: Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, EMB, SSB. }
  function Astronomy_BodyCode(Name: PChar): AstroBody; cdecl; external LibraryName;

  // Creates an observer object that represents a location on or near the surface of the Earth.
  function Astronomy_MakeObserver(Latitude: Double; Longitude: Double; Height: Double): AstroObserver; cdecl; external LibraryName;

  // Returns the computer's current date and time in the form of an #AstroTime.
  function Astronomy_CurrentTime(): AstroTime; cdecl; external LibraryName;

  // Creates an #AstroTime value from a given calendar date and time.
  function Astronomy_MakeTime(Year: Int32; Month: Int32; Day: Int32; Hour: Int32; Minute: Int32; Second: Double): AstroTime; cdecl; external LibraryName;

  // Creates an #AstroTime value from a given calendar date and time.
  function Astronomy_TimeFromUtc(UTC: AstroUTC): AstroTime; cdecl; external LibraryName;

  // Determines the calendar year, month, day and time from an #AstroTime value.
  function Astronomy_UtcFromTime(Time: AstroTime): AstroUTC; cdecl; external LibraryName;

  // Formats an #AstroTime value as an ISO 8601 string.
  function Astronomy_FormatTime(Time: AstroTime; Format: AstroTimeFormat; Text: PChar; Size: Int64): AstroStatus; cdecl; external LibraryName;

  // Converts a J2000 day value to an #AstroTime value.
  function Astronomy_TimeFromDays(UT: Double): AstroTime; cdecl; external LibraryName;

  // Converts a terrestrial time value into an #AstroTime value.
  function Astronomy_TerrestrialTime(TT: Double): AstroTime; cdecl; external LibraryName;

  // Calculates the sum or difference of an #AstroTime with a specified floating point number of days.
  function Astronomy_AddDays(Time: AstroTime; Days: Double): AstroTime; cdecl; external LibraryName;

  // Calculates Greenwich Apparent Sidereal Time (GAST) and returns in sidereal hours.
  function Astronomy_SiderealTime(Time: PAstroTime): Double; cdecl; external LibraryName;

  // Calculates the distance from a body to the Sun at a given time.
  function Astronomy_HelioDistance(Body: AstroBody; Time: AstroTime): AstroFuncResult; cdecl; external LibraryName;

  // Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.
  function Astronomy_HelioVector(Body: AstroBody; Time: AstroTime): AstroVector; cdecl; external LibraryName;

  // Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.
  function Astronomy_GeoVector(Body: AstroBody; Time: AstroTime; Aberration: AstroAberration): AstroVector; cdecl; external LibraryName;

  // Calculates equatorial geocentric position of the Moon at a given time in the J2000 equatorial system.
  function Astronomy_GeoMoon(Time: AstroTime): AstroVector; cdecl; external LibraryName;

  // Calculates spherical ecliptic geocentric postition of the Moon in the equatorial system of date.
  function Astronomy_EclipticGeoMoon(Time: AstroTime): AstroSpherical; cdecl; external LibraryName;

  // Calculates equatorial geocentric postition of the Moon at a given time in the J2000 equatorial system.
  function Astronomy_GeoMoonState(Time: AstroTime): AstroStateVector; cdecl; external LibraryName;

  // Calculates the geocentric position and velocity of the Earth/Moon barycenter in the J2000 equatorial system.
  function Astronomy_GeoEmbState(Time: AstroTime): AstroStateVector; cdecl; external LibraryName;

  // Calculates the Moon's libration angles at a given moment in time.
  function Astronomy_Libration(Time: AstroTime): AstroLibration; cdecl; external LibraryName;

  // Calculates barycentric position and velocity vectors for the given body (everything except stars).
  function Astronomy_BaryState(Body: AstroBody; Time: AstroTime): AstroStateVector; cdecl; external LibraryName;

  // Calculates heliocentric position and velocity vectors for the given body (including stars).
  function Astronomy_HelioState(Body: AstroBody; Time: AstroTime): AstroStateVector; cdecl; external LibraryName;

  // Returns the product of mass and universal gravitational constant of a Solar System body.
  function Astronomy_MassProduct(Body: AstroBody): Double; cdecl; external LibraryName;

  // Returns the average number of days it takes for a planet to orbit the Sun.
  function Astronomy_PlanetOrbitalPeriod(Body: AstroBody): Double; cdecl; external LibraryName;

  // Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.
  function Astronomy_LagrangePoint(Point: Int32; Time: AstroTime; MajorBody: AstroBody; MinorBody: AstroBody): AstroStateVector; cdecl; external LibraryName;

  // Calculates one of the 5 Lagrange points from body masses and state vectors.
  function Astronomy_LagrangePointFast(Point: Int32; MajorState: AstroStateVector; MajorMass: Double; MinorState: AstroStateVector; MinorMass: Double): AstroStateVector; cdecl; external LibraryName;

  // Calculates joviocentric positions and velocities of Jupiter's largest 4 moons.
  function Astronomy_JupiterMoons(Time: AstroTime): AstroJupiterMoons; cdecl; external LibraryName;

  // Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.
  function Astronomy_Equator(Body: AstroBody; Time: PAstroTime; Observer: AstroObserver; EquDate: AstroEquatorDate; Aberration: AstroAberration): AstroEquatorial; cdecl; external LibraryName;

  // Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.
  function Astronomy_ObserverVector(Time: PAstroTime; Observer: AstroObserver; EquDate: AstroEquatorDate): AstroVector; cdecl; external LibraryName;

  // Calculates geocentric equatorial postition and velocity of an observer on the surface of the Earth.
  function Astronomy_ObserverState(Time: PAstroTime; Observer: AstroObserver; EquDate: AstroEquatorDate): AstroStateVector; cdecl; external LibraryName;

  // Calculates the geographic location corresponding to an equatorial vector.
  function Astronomy_VectorObserver(Vector: PAstroVector; EquDate: AstroEquatorDate): AstroObserver; cdecl; external LibraryName;

  // Calculates the gravitational acceleration experienced by an observer on the Earth and returns acceleration expressed in meters per second squared.
  function Astronomy_ObserverGravity(Latitude: Double; Height: Double): Double; cdecl; external LibraryName;

  // Calculates geocentric ecliptic coordinates for the Sun.
  function Astronomy_SunPosition(Time: AstroTime): AstroEcliptic; cdecl; external LibraryName;

  // Converts a J2000 mean equator (EQJ) vector to a true ecliptic of date (ETC) vector and angles.
  function Astronomy_Ecliptic(EQJ: AstroVector): AstroEcliptic; cdecl; external LibraryName;

  // Calculates heliocentric ecliptic longitude of a body.
  function Astronomy_EclipticLongitude(Body: AstroBody; Time: AstroTime): AstroAngleResult; cdecl; external LibraryName;

  // Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.
  function Astronomy_Horizon(Time: PAstroTime; Observer: AstroObserver; RA: Double; Dec: Double; Refraction: AstroRefraction): AstroHorizon; cdecl; external LibraryName;

  // Returns the angle between the given body and the Sun, as seen from the Earth.
  function Astronomy_AngleFromSun(Body: AstroBody; Time: AstroTime): AstroAngleResult; cdecl; external LibraryName;

  // Determins visibility of a celestial body realtive to the Sun, as seen from the Earth.
  function Astronomy_Elongation(Body: AstroBody; Time: AstroTime): AstroElongation; cdecl; external LibraryName;

  // Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.
  function Astronomy_SearchMaxElongation(Body: AstroBody; StartTime: AstroTime): AstroElongation; cdecl; external LibraryName;

  // Returns one body's ecliptic longitude with respect to another, as seen from Earth.
  function Astronomy_PairLongitude(Body1: AstroBody; Body2: AstroBody; Time: AstroTime): AstroAngleResult; cdecl; external LibraryName;

  // Searches for the time when the Earth and another planet are separated by a specified angle in eclipctic longitude, as seen from the Sun.
  function Astronomy_SearchRelativeLongitude(Body: AstroBody; TargetRelLon: Double; StartTime: AstroTime): AstroSearchResult; cdecl; external LibraryName;

  // Returns the Moon's phase as an angle from 0 to 360 degrees.
  function Astronomy_MoonPhase(Time: AstroTime): AstroAngleResult; cdecl; external LibraryName;

  // Searches for the time that the Moon reaches a specified phase.
  function Astronomy_SearchMoonPhase(TargetLon: Double; StartTime: AstroTime; LimitDays: Double): AstroSearchResult; cdecl; external LibraryName;

  // Finds the first lunar quarter after the specified date and time.
  function Astronomy_SearchMoonQuarter(StartTime: AstroTime): AstroMoonQuarter; cdecl; external LibraryName;

  // Continues searching for lunar quarters from previous search.
  function Astronomy_NextMoonQuarter(MQ: AstroMoonQuarter): AstroMoonQuarter; cdecl; external LibraryName;

  // Searches for a lunar eclipse.
  function Astronomy_SearchLunarEclipse(StartTime: AstroTime): AstroLunarEclipse; cdecl; external LibraryName;

  // Searches for the next lunar eclipse in a series.
  function Astronomy_NextLunarEclipse(PrevEclipseTime: AstroTime): AstroLunarEclipse; cdecl; external LibraryName;

  // Searches for a solar eclipse visible anywhere on the Earth's surface.
  function Astronomy_SearchGlobalSolarEclipse(StartTime: AstroTime): AstroGlobalSolarEclipse; cdecl; external LibraryName;

  // Searches for the next global solar eclipse in a series.
  function Astronomy_NextGlobalSolarEclipse(PrevEclipseTime: AstroTime): AstroGlobalSolarEclipse; cdecl; external LibraryName;

  // Searches for a solar eclipse visible at specific location on the Earth's surface.
  function Astronomy_SearchLocalSolarEclipse(StartTime: AstroTime; Observer: AstroObserver): AstroLocalSolarEclipse; cdecl; external LibraryName;

  // Searches for the next local solar eclipse in a series.
  function Astronomy_NextLocalSolarEclipse(PrevEclipseTime: AstroTime; Observer: AstroObserver): AstroLocalSolarEclipse; cdecl; external LibraryName;

  // Searches for the first transit of Mercury or Venus after a given date.
  function Astronomy_SearchTransit(Body: AstroBody; StartTime: AstroTime): AstroTransit; cdecl; external LibraryName;

  // Searches for another transit of Mercury or Venus.
  function Astronomy_NextTransit(Body: AstroBody; PrevTransitTime: AstroTime): AstroTransit; cdecl; external LibraryName;

  // Searches for a time when the Moon's center crosses through the ecliptic plane.
  function Astronomy_SearchMoonNode(StartTime: AstroTime): AstroNodeEvent; cdecl; external LibraryName;

  // Searches for the next time when the Moon's center crosses through the ecliptic plane.
  function Astronomy_NextMoonNode(PrevNode: AstroNodeEvent): AstroNodeEvent; cdecl; external LibraryName;

  // Searches for a time at which a function's value increases through zero.
  function Astronomy_Search(Func: PAstroSearchFunc; Context: Pointer; T1: AstroTime; T2: AstroTime; DTToleranceSeconds: Double): AstroSearchResult; cdecl; external LibraryName;

  // Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.
  function Astronomy_SearchSunLongitude(TargetLon: Double; StartTime: AstroTime; LimitDays: Double): AstroSearchResult; cdecl; external LibraryName;

  // Searches for the time when the center of a body reaches a specified hour angle as seen by an observer on the Earth.
  function Astronomy_SearchHourAngleEx(Body: AstroBody; Observer: AstroObserver; HourAngle: Double; StartTime: AstroTime; Direction: Int32): AstroHourAngle; cdecl; external LibraryName;

  // Finds the hour angle of a body for a given observer and time.
  function Astronomy_HourAngle(Body: AstroBody; Time: PAstroTime; Observer: AstroObserver): AstroFuncResult; cdecl; external LibraryName;

  // Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.
  function Astronomy_SearchRiseSetEx(Body: AstroBody; Observer: AstroObserver; Direction: AstroDirection; StartTime: AstroTime; LimitDays: Double; MetersAboveGround: Double): AstroSearchResult; cdecl; external LibraryName;

  // Finds the next time the center of a body passes through a given altitude.
  function Astronomy_SearchAltitude(Body: AstroBody; Observer: AstroObserver; Direction: AstroDirection; StartTime: AstroTime; LimitDays: Double; Altitude: Double): AstroSearchResult; cdecl; external LibraryName;

  // Calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.
  function Astronomy_Atmosphere(ElevationMeters: Double): AstroAtmosphere; cdecl; external LibraryName;

  // Calculates information about a body's rotation axis at a given time.
  function Astronomy_RotationAxis(Body: AstroBody; Time: PAstroTime): AstroAxis; cdecl; external LibraryName;

  // Finds both equinoxes and both solstices for a given calendar year.
  function Astronomy_Seasons(Year: Int32): AstroSeasons; cdecl; external LibraryName;

  // Finds visual magnitude, phase angle, and other illumination information about a celestial body.
  function Astronomy_Illumination(Body: AstroBody; Time: AstroTime): AstroIllum; cdecl; external LibraryName;

  // Searches for the date and time Venus will next appear brightest as seen from the Earth.
  function Astronomy_SearchPeakMagnitude(Body: AstroBody; StartTime: AstroTime): AstroIllum; cdecl; external LibraryName;

  // Finds the date and time of the Moon's closest distance (perigee) or farthest distance (apogee) with the respect to the Earth.
  function Astronomy_SearchLunarApsis(StartTime: AstroTime): AstroApsis; cdecl; external LibraryName;

  // Finds the next lunar perigee or apogee event in a series.
  function Astronomy_NextLunarApsis(Apsis: AstroApsis): AstroApsis; cdecl; external LibraryName;

  // Finds the date and time of a planet's perihelion (closest approach to the Sun) or aphelion (farthest distance from the Sun) after a given time.
  function Astronomy_SearchPlanetApsis(Body: AstroBody; StartTime: AstroTime): AstroApsis; cdecl; external LibraryName;

  // Finds the next planetary perihelion or aphelion event in a series.
  function Astronomy_NextPlanetApsis(Body: AstroBody; Apsis: AstroApsis): AstroApsis; cdecl; external LibraryName;

  // Creates an identity rotation matrix.
  function Astronomy_IdentityMatrix(): AstroRotation; cdecl; external LibraryName;

  // Calculatues the inverse of a rotation matrix.
  function Astronomy_InverseRotation(Rotation: AstroRotation): AstroRotation; cdecl; external LibraryName;

  // Creates a rotation based on applying one rotation followed by another.
  function Astronomy_CombineRotation(A: AstroRotation; B: AstroRotation): AstroRotation; cdecl; external LibraryName;

  // Re-orients a rotation matrix by pivoting it by an angle around one of its axes.
  function Astronomy_Pivot(Rotation: AstroRotation; Axis: Int32; Angle: Double): AstroRotation; cdecl; external LibraryName;

  // Converts spherical coordinates to Cartesian coordinates.
  function Astronomy_VectorFromSphere(Sphere: AstroSpherical; Time: AstroTime): AstroVector; cdecl; external LibraryName;

  // Converts Cartesian coordinates to spherical coordinates.
  function Astronomy_SphereFromVector(Vector: AstroVector): AstroSpherical; cdecl; external LibraryName;

  // Given an equatorial vector, calculates equatorial angular coordinates.
  function Astronomy_EquatorFromVector(Vector: AstroVector): AstroEquatorial; cdecl; external LibraryName;

  // Given apparent angular horizontal coordinates in `sphere`, calculates horizontal vector.
  function Astronomy_VectorFromHorizon(Sphere: AstroSpherical; Time: AstroTime; Refraction: AstroRefraction): AstroVector; cdecl; external LibraryName;

  // Converts Cartesian coordinates to horizontal coordinates.
  function Astronomy_HorizonFromVector(Vector: AstroVector; Refraction: AstroRefraction): AstroSpherical; cdecl; external LibraryName;

  // Applies a rotation to a vector, yielding a rotated vector.
  function Astronomy_RotateVector(Rotation: AstroRotation; Vector: AstroVector): AstroVector; cdecl; external LibraryName;

  // Applies a rotation to a state vector, yielding a rotated vector.
  function Astronomy_RotateState(Rotation: AstroRotation; State: AstroStateVector): AstroStateVector; cdecl; external LibraryName;

  // Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean equator (EQJ).
  function Astronomy_Rotation_EQD_EQJ(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean ecliptic (ECL).
  function Astronomy_Rotation_EQD_ECL(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Returns a rotation matrix from equator of date (EQD) to true ecliptic of date (ECT).
  function Astronomy_Rotation_EQD_ECT(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
  function Astronomy_Rotation_EQD_HOR(Time: PAstroTime; Observer: AstroObserver): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean equator (EQJ) to equatorial of-date (EQD).
  function Astronomy_Rotation_EQJ_EQD(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT).
  function Astronomy_Rotation_EQJ_ECT(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL).
  function Astronomy_Rotation_EQJ_ECL(): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR).
  function Astronomy_Rotation_EQJ_HOR(Time: PAstroTime; Observer: AstroObserver): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean ecliptic (ECL) to equatorial of-date (EQD).
  function Astronomy_Rotation_ECL_EQD(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ).
  function Astronomy_Rotation_ECL_EQJ(): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR).
  function Astronomy_Rotation_ECL_HOR(Time: PAstroTime; Observer: AstroObserver): AstroRotation; cdecl; external LibraryName;

  // Calculates a rotation matrix from true ecliptic of date (ECT) to J2000 mean equator (EQJ).
  function Astronomy_Rotation_ECT_EQJ(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calulcates a rotation matrix from true ecliptic of date (ECT) to equator of date (EQD).
  function Astronomy_Rotation_ECT_EQD(Time: PAstroTime): AstroRotation; cdecl; external LibraryName;

  // Calulcates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).
  function Astronomy_Rotation_HOR_EQD(Time: PAstroTime; Observer: AstroObserver): AstroRotation; cdecl; external LibraryName;

  // Calulcates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).
  function Astronomy_Rotation_HOR_EQJ(Time: PAstroTime; Observer: AstroObserver): AstroRotation; cdecl; external LibraryName;

  // Calulcates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL).
  function Astronomy_Rotation_HOR_ECL(Time: PAstroTime; Observer: AstroObserver): AstroRotation; cdecl; external LibraryName;

  // Calulcates a rotation matrix from J2000 mean ecliptic (EQJ) to galactic (GAL).
  function Astronomy_Rotation_EQJ_GAL(): AstroRotation; cdecl; external LibraryName;

  // Calulcates a rotation matrix from galactic (GAL) to J2000 mean ecliptic (EQJ).
  function Astronomy_Rotation_GAL_EQJ(): AstroRotation; cdecl; external LibraryName;

  // Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
  function Astronomy_Refraction(Refraction: AstroRefraction; Altitude: Double): Double; cdecl; external LibraryName;

  // Calculates the inverse of an atmospheric refraction angle.
  function Astronomy_InverseRefraction(Refraction: AstroRefraction; BentAltitude: Double): Double; cdecl; external LibraryName;

  // Determines the constellation that contains the given point in the sky.
  function Astronomy_Constellation(RA: Double; Dec: Double): AstroConstellation; cdecl; external LibraryName;

{ Allocate and initialize a gravity step simulator.
  Prepares to simulate a series of incremental time steps,
  simulating the movement of zero or more small bodies through the Solar System
  acting under gravitational attraction from the Sun and planets. }
  function Astronomy_GravSimInit(SimOut: PAstroGravSim; OriginBody: AstroBody; Time: AstroTime; NumBodies: Int32; BodyStateArray: PAstroStateVector): AstroStatus; cdecl; external LibraryName;

  // Advances a gravity simulation by a small time step.
  function Astronomy_GravSimUpdate(Sim: PAstroGravSim; Time: AstroTime; NumBodies: Int32; BodyStateArray: PAstroStateVector): AstroStatus; cdecl; external LibraryName;

  // Get the position and velocity of a Solar System body included in the simulation.
  function Astronomy_GravSimBodyState(Sim: PAstroGravSim; Body: AstroBody): AstroStateVector; cdecl; external LibraryName;

  // Returns the time of the current simulation step.
  function Astronomy_GravSimTime(Sim: PAstroGravSim): AstroTime; cdecl; external LibraryName;

  // Returns the number of small bodies represented in this simulation.
  function Astronomy_GravSimNumBodies(Sim: PAstroGravSim): Int32; cdecl; external LibraryName;

  // Returns the body whose center is the coordinate origin that small bodies are referenced to.
  function Astronomy_GravSimOrigin(Sim: PAstroGravSim): AstroBody; cdecl; external LibraryName;

  // Exchange the current time step with the previous time step.
  procedure Astronomy_GravSimSwap(Sim: PAstroGravSim); cdecl; external LibraryName;

  // Releases memory allocated to a gravity simulator object.
  procedure Astronomy_GravSimFree(Sim: PAstroGravSim); cdecl; external LibraryName;

  // Solve for light travel time of a vector function.
  function Astronomy_CorrectLightTravel(Context: Pointer; Func: PAstroPositionFunc; Time: AstroTime): AstroVector; cdecl; external LibraryName;

  // Solve for light travel time correction of apparent position.
  function Astronomy_BackdatePosition(Time: AstroTime; ObserverBody: AstroBody; TargetBody: AstroBody; Aberration: AstroAberration): AstroVector; cdecl; external LibraryName;

  // Assign equatorial coordinates to a user-defined star.
  function Astronomy_DefineStar(Body: AstroBody; RA: Double; Dec: Double; DistanceLightYears: Double): AstroStatus; cdecl; external LibraryName;

implementation

end.

