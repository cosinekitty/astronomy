# Document

The file `iau_rotation_elements.pdf` is a copy of the source document [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009](https://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf), stored here as an archival reference. It shows where I got the formulas for calculating planetary rotation axes as a function of time.

# Planetary Rotation Axis Data
The `*.txt` data files were generated using the [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/app.html#/) online ephemeris tool.

Here are the settings I used to generate these data files:

- Ephemeris Type = Observer Table
- Target Body = *whichever planet the file is named for*
- Observer Location = Geocentric (but must use the Sun for generating Earth)
- Time Specification = 1970-01-01 to 2050-01-01, step = 30 days
- Table Settings:
    - check 32. "North Pole RA / DEC"
    - Reference Frame = ICRF
    - Date/Time format = calendar and Julian Day number
    - Time digits = HH:MM
    - Angle format = decimal degrees
    - Refraction model = no refraction (airless)
    - Range units = astronomical units
    - RTS flag = disabled
    - Object summary = checked
