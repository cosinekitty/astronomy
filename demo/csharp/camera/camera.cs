using System;
using demo_helper;
using CosineKitty;

//    camera.cs  -  by Don Cross - 2021-03-27
//
//    Example C# program for Astronomy Engine:
//    https://github.com/cosinekitty/astronomy
//
//    Suppose you want to photograph the Moon,
//    and you want to know what it will look like in the photo.
//    Given a location on the Earth, and a date/time,
//    this program calculates the orientation of the sunlit
//    side of the Moon with respect to the top of your
//    photo image. It assumes the camera faces directly
//    toward the Moon's azimuth and tilts upward to its
//    altitude angle above the horizon.

namespace camera
{
    class Program
    {
        static int Main(string[] args)
        {
            DemoHelper.ParseArgs("camera", args, out Observer observer, out AstroTime time);
            return CameraImage(observer, time);
        }

        static int CameraImage(Observer observer, AstroTime time)
        {
            const double tolerance = 1.0e-15;

            // Calculate the topocentric equatorial coordinates of date for the Moon.
            // Assume aberration does not matter because the Moon is so close and has such a small relative velocity.
            Equatorial moon_equ = Astronomy.Equator(Body.Moon, time, observer, EquatorEpoch.OfDate, Aberration.None);

            // Also calculate the Sun's topocentric position in the same coordinate system.
            Equatorial sun_equ = Astronomy.Equator(Body.Sun, time, observer, EquatorEpoch.OfDate, Aberration.None);

            // Get the Moon's horizontal coordinates, so we know how much to pivot azimuth and altitude.
            Topocentric moon_hor = Astronomy.Horizon(time, observer, moon_equ.ra, moon_equ.dec, Refraction.None);
            Console.WriteLine($"Moon horizontal position: azimuth = {moon_hor.azimuth:F3}, altitude = {moon_hor.altitude:F3}");

            // Get the rotation matrix that converts equatorial to horizontal coordintes for this place and time.
            RotationMatrix rot = Astronomy.Rotation_EQD_HOR(time, observer);

            // Modify the rotation matrix in two steps:
            // First, rotate the orientation so we are facing the Moon's azimuth.
            // We do this by pivoting around the zenith axis.
            // Horizontal axes are: 0 = north, 1 = west, 2 = zenith.
            // Tricky: because the pivot angle increases counterclockwise, and azimuth
            // increases clockwise, we undo the azimuth by adding the positive value.
            rot = Astronomy.Pivot(rot, 2, moon_hor.azimuth);

            // Second, pivot around the leftward axis to bring the Moon to the camera's altitude level.
            // From the point of view of the leftward axis, looking toward the camera,
            // adding the angle is the correct sense for subtracting the altitude.
            rot = Astronomy.Pivot(rot, 1, moon_hor.altitude);

            // As a sanity check, apply this rotation to the Moon's equatorial (EQD) coordinates and verify x=0, y=0.
            AstroVector vec = Astronomy.RotateVector(rot, moon_equ.vec);

            // Convert to unit vector.
            double radius = vec.Length();
            vec.x /= radius;
            vec.y /= radius;
            vec.z /= radius;
            Console.WriteLine($"Moon check: x = {vec.x}, y = {vec.y}, z = {vec.z}");
            if (!double.IsFinite(vec.x) || Math.Abs(vec.x - 1.0) > tolerance)
            {
                Console.WriteLine("Excessive error in moon check (x).");
                return 1;
            }

            if (!double.IsFinite(vec.y) || Math.Abs(vec.y) > tolerance)
            {
                Console.WriteLine("Excessive error in moon check (y).");
                return 1;
            }

            if (!double.IsFinite(vec.z) || Math.Abs(vec.z) > tolerance)
            {
                Console.WriteLine("Excessive error in moon check (z).");
                return 1;
            }

            // Apply the same rotation to the Sun's equatorial vector.
            // The x- and y-coordinates now tell us which side appears sunlit in the camera!

            vec = Astronomy.RotateVector(rot, sun_equ.vec);

            // Don't bother normalizing the Sun vector, because in AU it will be close to unit anyway.
            Console.WriteLine($"Sun vector: x = {vec.x:F6}, y = {vec.y:F6}, z = {vec.z:F6}");

            // Calculate the tilt angle of the sunlit side, as seen by the camera.
            // The x-axis is now pointing directly at the object, z is up in the camera image, y is to the left.
            double tilt = Astronomy.RAD2DEG * Math.Atan2(vec.z, vec.y);
            Console.WriteLine($"Tilt angle of sunlit side of the Moon = {tilt:F3} degrees counterclockwise from up.");

            IllumInfo illum = Astronomy.Illumination(Body.Moon, time);

            Console.WriteLine($"Moon magnitude = {illum.mag:F2}, phase angle = {illum.phase_angle:F2} degrees.");

            double angle = Astronomy.AngleFromSun(Body.Moon, time);

            Console.WriteLine($"Angle between Moon and Sun as seen from Earth = {angle:F2} degrees.");
            return 0;
        }
    }
}
