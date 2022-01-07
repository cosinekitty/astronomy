#!/usr/bin/env python3
#
#    camera.py  -  by Don Cross - 2021-03-27
#
#    Example Python program for Astronomy Engine:
#    https://github.com/cosinekitty/astronomy
#
#    Suppose you want to photograph the Moon,
#    and you want to know what it will look like in the photo.
#    Given a location on the Earth, and a date/time,
#    this program calculates the orientation of the sunlit
#    side of the Moon with respect to the top of your
#    photo image. It assumes the camera faces directly
#    toward the Moon's azimuth and tilts upward to its
#    altitude angle above the horizon.
#
#    To execute, run the command:
#
#    python3 camera.py latitude longitude [yyyy-mm-ddThh:mm:ssZ]
#
import sys
import math
import astronomy
from astro_demo_common import ParseArgs

def Camera(observer, time):
    tolerance = 1.0e-15

    # Calculate the topocentric equatorial coordinates of date for the Moon.
    # Assume aberration does not matter because the Moon is so close and has such a small relative velocity.
    moon_equ = astronomy.Equator(astronomy.Body.Moon, time, observer, True, False)

    # Also calculate the Sun's topocentric position in the same coordinate system.
    sun_equ = astronomy.Equator(astronomy.Body.Sun, time, observer, True, False)

    # Get the Moon's horizontal coordinates, so we know how much to pivot azimuth and altitude.
    moon_hor = astronomy.Horizon(time, observer, moon_equ.ra, moon_equ.dec, astronomy.Refraction.Airless)
    print('Moon horizontal position: azimuth = {:0.3f}, altitude = {:0.3f}'.format(moon_hor.azimuth, moon_hor.altitude))

    # Get the rotation matrix that converts equatorial to horizontal coordintes for this place and time.
    rot = astronomy.Rotation_EQD_HOR(time, observer)

    # Modify the rotation matrix in two steps:
    # First, rotate the orientation so we are facing the Moon's azimuth.
    # We do this by pivoting around the zenith axis.
    # Horizontal axes are: 0 = north, 1 = west, 2 = zenith.
    # Tricky: because the pivot angle increases counterclockwise, and azimuth
    # increases clockwise, we undo the azimuth by adding the positive value.
    rot = astronomy.Pivot(rot, 2, moon_hor.azimuth)

    # Second, pivot around the leftward axis to bring the Moon to the camera's altitude level.
    # From the point of view of the leftward axis, looking toward the camera,
    # adding the angle is the correct sense for subtracting the altitude.
    rot = astronomy.Pivot(rot, 1, moon_hor.altitude)

    # As a sanity check, apply this rotation to the Moon's equatorial (EQD) coordinates and verify x=0, y=0.
    vec = astronomy.RotateVector(rot, moon_equ.vec)

    # Convert to unit vector.
    radius = vec.Length()
    vec.x /= radius
    vec.y /= radius
    vec.z /= radius
    print('Moon check: x={:0.6f}, y={:0.6f}, z={:0.6f}'.format(vec.x, abs(vec.y), abs(vec.z)))
    if abs(vec.x - 1.0) > tolerance:
        print("Excessive error in moon check (x).")
        return 1

    if abs(vec.y) > tolerance:
        print("Excessive error in moon check (y).")
        return 1

    if abs(vec.z) > tolerance:
        print("Excessive error in moon check (z).")
        return 1

    # Apply the same rotation to the Sun's equatorial vector.
    # The x- and y-coordinates now tell us which side appears sunlit in the camera!

    vec = astronomy.RotateVector(rot, sun_equ.vec)

    # Don't bother normalizing the Sun vector, because in AU it will be close to unit anyway.
    print('Sun vector: {}'.format(vec))

    # Calculate the tilt angle of the sunlit side, as seen by the camera.
    # The x-axis is now pointing directly at the object, z is up in the camera image, y is to the left.
    tilt = math.degrees(math.atan2(vec.z, vec.y))
    print('Tilt angle of sunlit side of the Moon = {:0.3f} degrees counterclockwise from up.'.format(tilt))
    illum = astronomy.Illumination(astronomy.Body.Moon, time)
    print('Moon magnitude = {:0.2f}, phase angle = {:0.2f} degrees.'.format(illum.mag, illum.phase_angle))
    angle = astronomy.AngleFromSun(astronomy.Body.Moon, time)
    print('Angle between Moon and Sun as seen from Earth = {:0.2f} degrees.'.format(angle))
    return 0

if __name__ == '__main__':
    observer, time = ParseArgs(sys.argv)
    sys.exit(Camera(observer, time))
