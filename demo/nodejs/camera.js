/*
    camera.js  -  by Don Cross - 2021-03-26

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Given an observer's location on the Earth and a date/time,
    calculates the angle of the sunlit side of the Moon as
    seen through a camera aimed at it.

    To execute, run the command:
    node camera latitude longitude [date]
*/

const Astronomy = require('./astronomy.js');

function ParseNumber(text, name) {
    const x = Number(text);
    if (!Number.isFinite(x)) {
        console.error(`ERROR: Not a valid numeric value for ${name}: "${text}"`);
        process.exit(1);
    }
    return x;
}

function ParseDate(text) {
    const d = new Date(text);
    if (!Number.isFinite(d.getTime())) {
        console.error(`ERROR: Not a valid date: "${text}"`);
        process.exit(1);
    }
    return d;
}

function Camera(observer, time) {
    const tolerance = 1.0e-15;
    const RAD2DEG = 57.295779513082321;

    // Calculate the topocentric equatorial coordinates of date for the Moon.
    // Assume aberration does not matter because the Moon is so close and has such a small relative velocity.
    const moon_equ = Astronomy.Equator(Astronomy.Body.Moon, time, observer, true, false);

    // Also calculate the Sun's topocentric position in the same coordinate system.
    const sun_equ = Astronomy.Equator(Astronomy.Body.Sun, time, observer, true, false);

    // Get the Moon's horizontal coordinates, so we know how much to pivot azimuth and altitude.
    const moon_hor = Astronomy.Horizon(time, observer, moon_equ.ra, moon_equ.dec, false);
    console.log(`Moon horizontal position: azimuth = ${moon_hor.azimuth.toFixed(3)}, altitude = ${moon_hor.altitude.toFixed(3)}`);

    // Get the rotation matrix that converts equatorial to horizontal coordintes for this place and time.
    let rot = Astronomy.Rotation_EQD_HOR(time, observer);

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
    let vec = Astronomy.RotateVector(rot, moon_equ.vec);

    // Convert to unit vector.
    const radius = vec.Length();
    vec.x /= radius;
    vec.y /= radius;
    vec.z /= radius;
    console.log(`Moon check: x = ${vec.x.toFixed(6)}, y = ${Math.abs(vec.y).toFixed(6)}, z = ${Math.abs(vec.z).toFixed(6)}`);
    if (!Number.isFinite(vec.x) || Math.abs(vec.x - 1.0) > tolerance) {
        console.error("Excessive error in moon check (x).");
        return 1;
    }

    if (!Number.isFinite(vec.y) || Math.abs(vec.y) > tolerance) {
        console.error("Excessive error in moon check (y).");
        return 1;
    }

    if (!Number.isFinite(vec.z) || Math.abs(vec.z) > tolerance) {
        console.error("Excessive error in moon check (z).");
        return 1;
    }

    // Apply the same rotation to the Sun's equatorial vector.
    // The x- and y-coordinates now tell us which side appears sunlit in the camera!

    vec = Astronomy.RotateVector(rot, sun_equ.vec);

    // Don't bother normalizing the Sun vector, because in AU it will be close to unit anyway.
    console.log(`Sun vector: x = ${vec.x.toFixed(6)}, y = ${vec.y.toFixed(6)}, z = ${vec.z.toFixed(6)}`);

    // Calculate the tilt angle of the sunlit side, as seen by the camera.
    // The x-axis is now pointing directly at the object, z is up in the camera image, y is to the left.
    const tilt = RAD2DEG * Math.atan2(vec.z, vec.y);
    console.log(`Tilt angle of sunlit side of the Moon = ${tilt.toFixed(3)} degrees counterclockwise from up.`);

    const illum = Astronomy.Illumination(Astronomy.Body.Moon, time);

    console.log(`Moon magnitude = ${illum.mag.toFixed(2)}, phase angle = ${illum.phase_angle.toFixed(2)} degrees.`);

    const angle = Astronomy.AngleFromSun(Astronomy.Body.Moon, time);

    console.log(`Angle between Moon and Sun as seen from Earth = ${angle.toFixed(2)} degrees.`);
}

function Demo() {
    if (process.argv.length === 4 || process.argv.length === 5) {
        const latitude  = ParseNumber(process.argv[2]);
        const longitude = ParseNumber(process.argv[3]);
        const observer = new Astronomy.Observer(latitude, longitude, 0);
        const time = Astronomy.MakeTime((process.argv.length === 5) ? ParseDate(process.argv[4]) : new Date());
        Camera(observer, time);
        process.exit(0);
    } else {
        console.log('USAGE: node camera latitude longitude [date]');
        process.exit(1);
    }
}

Demo();
