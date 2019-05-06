### Overview

[![Build Status](https://travis-ci.org/cosinekitty/astronomy.svg)](https://travis-ci.org/cosinekitty/astronomy)

A suite of open source libraries for calculating positions of the Sun, Moon, and planets.

This code is designed to be small, fast, and accurate to within &plusmn;1 arcminute.
It is based on the authoritative and well-tested models
[VSOP87](https://en.wikipedia.org/wiki/VSOP_(planets))
and 
[NOVAS C 3.1](https://aa.usno.navy.mil/software/novas/novas_c/novasc_info.php).

Rigorously unit-tested against NOVAS, [JPL Horizons](https://ssd.jpl.nasa.gov/horizons.cgi),
and other reliable sources of ephemeris data.

### Features

- Provides calculations for the Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, and Pluto.

- Calculates for any calendar date and time between the years 1600 and 2200.

- Provides heliocentric and geocentric Cartesian vectors of all the above bodies.

- Determines apparent horizon-based positions for an observer anywhere on the Earth, 
  given that observer's latitude, longitude, and elevation in meters. 
  Optionally corrects for atmospheric refraction.

- Rise, set and culmination times of Sun, Moon, and planets.

- Date and time of Moon phases: new, first quarter, full, third quarter 
  (or anywhere in between as expressed in degrees of ecliptic longitude).

- Finds equinoxes and solstices for a given calendar year.

- Finds apparent visual magnitudes of all the supported celestial bodies.

- Predicts dates of planetary conjunctions and oppositions.

- Predicts dates of Venus' peak visual magnitude.

- Predicts dates of maximum elongation for Mercury and Venus.

# Supported Languages

### JavaScript

[API Reference](source/js/README.md)

### C/C++

(Coming soon.)

### Go

(Coming soon.)

### Python

(Coming soon.)
