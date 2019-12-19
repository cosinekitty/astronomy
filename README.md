# Astronomy Engine <span style="vertical-align: middle;">[![Build Status](https://travis-ci.org/cosinekitty/astronomy.svg)](https://travis-ci.org/cosinekitty/astronomy)</span> <span style="vertical-align: middle;">[![npm package](https://img.shields.io/npm/v/astronomy-engine.svg)](https://www.npmjs.com/package/astronomy-engine)</span>

### Supported Programming Languages

<table style="border-width: 0px;" cellspacing="0" cellpadding="10">
    <tr>
        <td style="text-align: center;">
            <div>C</div>
            <div><img src="source/c/c_language.svg" width="100" height="100" alt="C" /></div>
        </td>
        <td style="text-align: center;">
            <div>C#</div>
            <div><img src="source/csharp/csharp_language.svg" width="100" height="100" alt="C#" /></div>
        </td>
        <td style="text-align: center;">
            <div>Browser</div>
            <div><img src="source/js/javascript.svg" width="100" height="100" alt="JavaScript" /></div>
        </td>
        <td style="text-align: center;">
            <div>Node.js</div>
            <div><img src="source/js/nodejs.svg" width="100" height="100" alt="Node.js" /></div>
        </td>
        <td style="text-align: center;">
            <div>Python</div>
            <div><img src="source/python/python_language.svg" width="100" height="100" alt="Python" /></div>
        </td>
    </tr>
    <tr>
        <td style="text-align: center;"><a href="demo/c/">Examples</a></td>
        <td style="text-align: center;">&nbsp;</td>
        <td style="text-align: center;"><a href="demo/browser/">Examples</a></td>
        <td style="text-align: center;"><a href="demo/nodejs/">Examples</a></td>
        <td style="text-align: center;"><a href="demo/python/">Examples</a></td>
    </tr>
    <tr>
        <td style="text-align: center;"><a href="source/c/">Documentation</a></td>
        <td style="text-align: center;">&nbsp;</td>
        <td style="text-align: center;"><a href="source/js/">Documentation</a></td>
        <td style="text-align: center;"><a href="source/js/">Documentation</a></td>
        <td style="text-align: center;"><a href="source/python/">Documentation</a></td>
    </tr>
</table>

### Overview

The Astronomy Engine is a suite of open source libraries for calculating positions of
the Sun, Moon, and planets, and for predicting interesting events like oppositions,
conjunctions, rise and set times, lunar phases, and more.

It supports several popular programming langauges with a consistent API.
Function and type names are uniform across all the supported languages.

The Astronomy Engine is designed to be small, fast, and accurate to within &plusmn;1 arcminute.
It is based on the authoritative and well-tested models
[VSOP87](https://en.wikipedia.org/wiki/VSOP_(planets))
and
[NOVAS C 3.1](https://aa.usno.navy.mil/software/novas/novas_c/novasc_info.php).

These libraries are rigorously unit-tested against NOVAS,
[JPL Horizons](https://ssd.jpl.nasa.gov/horizons.cgi),
and other reliable sources of ephemeris data.
Calculations are also verified to be identical among all the supported programming languages.

### Features

- Provides calculations for the Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, and Pluto.

- Calculates for any calendar date and time between the years 1700 and 2200.

- Provides heliocentric and geocentric Cartesian vectors of all the above bodies.

- Determines apparent horizon-based positions for an observer anywhere on the Earth,
  given that observer's latitude, longitude, and elevation in meters.
  Optionally corrects for atmospheric refraction.

- Calculates rise, set, and culmination times of Sun, Moon, and planets.

- Finds date and time of Moon phases: new, first quarter, full, third quarter
  (or anywhere in between as expressed in degrees of ecliptic longitude).

- Predicts lunar apogee and perigee dates, times, and distances.

- Predicts date and time of equinoxes and solstices for a given calendar year.

- Determines apparent visual magnitudes of all the supported celestial bodies.

- Predicts dates of planetary conjunctions and oppositions.

- Predicts dates of Venus' peak visual magnitude.

- Predicts dates of maximum elongation for Mercury and Venus.

- Converts angular and vector coordinates among the following orientations:
    - Equatorial J2000
    - Equatorial equator-of-date
    - Ecliptic J2000
    - Horizontal
