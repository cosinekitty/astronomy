# Astronomy Engine examples in Go

This directory contains example programs that demonstrate how to use
the [Astronomy Engine](https://github.com/cosinekitty/astronomy) Go package.
The examples are pure Go and use only the standard library; there is no CGO
and no external dependency beyond the astronomy package itself.

## Building

```bash
go build ./...
```

## Running the demos

All demos are subcommands of a single program:

```bash
go run . positions  +45.6 -90.7 2018-11-30T17:55:07.234Z
go run . riseset    45 -90 2019-06-14T21:45:25Z
go run . moonphase  2019-06-15T09:15:32.987Z
go run . seasons    2020
go run . culminate  +30 -90 2015-02-28T00:00:00Z
go run . lunar_eclipse 1988-01-01T00:00:00Z
go run . solar_time +38.88 -77.03 2023-02-12T17:00:00Z
go run . gravity    29.0 1000.0
go run . galactic   38.92056 -77.0658 22.793498 197.070510 2025-04-06T00:00:00Z
go run . horizon    +25.5 -85.3 2016-12-25T12:30:45Z
```

## What each demo shows

| Demo | Description |
|------|-------------|
| `positions` | Equatorial (J2000) and horizontal coordinates of the Sun, Moon, and planets for a given observer and time. |
| `riseset` | Sunrise, sunset, moonrise, and moonset times. |
| `moonphase` | The Moon's current phase and the next eight quarter phases. |
| `seasons` | The March/September equinoxes and June/December solstices of a year. |
| `culminate` | The next time each body reaches its highest point (crosses the meridian). |
| `lunar_eclipse` | The next 10 lunar eclipses after a given date. |
| `solar_time` | Uses the Sun's hour angle to find true solar noon. |
| `gravity` | The WGS 84 acceleration of gravity at a latitude and elevation. |
| `galactic` | Converts galactic coordinates to horizontal coordinates. |
| `horizon` | Rotates the ecliptic plane into the observer's horizontal frame. |

## Testing

Run all demos to confirm they work:

```bash
./demotest
```
