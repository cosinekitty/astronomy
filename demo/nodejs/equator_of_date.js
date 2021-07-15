/*
    equator_of_date.js  -  by Don Cross  -  2021-07-06

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy
*/

const UsageText = `
USAGE:  node equator_of_date.js [a|n] ra dec [yyyy-mm-ddThh:mm:ssZ]

Converts J2000 equatorial coordinates to
equator-of-date coordinates.

[a|n] = aberration correction / no aberration correction
ra  = J2000 Right ascension in sidereal hours (0 .. 24).
dec = J2000 Declination in degrees (-90 .. +90).
yyyy-mm-ddThh:mm:ssZ = Optional date and time in UTC.
(If omitted, the current date and time are used.)

This program prints out the right ascension and declination
of the same point in the sky, but expressed in the Earth's
equator at the given date and time.
`;

const Astronomy = require('./astronomy.js');

function ParseNumber(name, text, minValue, maxValue) {
    const x = Number(text);
    if (!Number.isFinite(x) || (x < minValue) || (x > maxValue)) {
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

function Demo() {
    if (process.argv.length < 5 || process.argv.length > 6) {
        console.log(UsageText);
        process.exit(1);
    } else {
        // Parse the command line arguments.
        const correct_aberration = (process.argv[2] === 'a');
        const ra = ParseNumber("RA", process.argv[3], 0, 24);
        const dec = ParseNumber("DEC", process.argv[4], -90, +90);
        const date = (process.argv.length > 5) ? ParseDate(process.argv[5]) : new Date();
        const time = Astronomy.MakeTime(date);
        console.log(`time = ${time}`);

        // Create a rotation matrix that converts J2000 equatorial (EQJ)
        // orientation to equator-of-date (EQD) orientation.
        const rot = Astronomy.Rotation_EQJ_EQD(time);

        // Convert the spherical angular EQJ coordinates to vector.
        // Multiply ra by 15 to convert sidereal hours to degrees.
        // Scale the vector by the speed of light so that we can optionally
        // apply an aberration correction below.
        const eqj_sphere = new Astronomy.Spherical(dec, 15*ra, Astronomy.C_AUDAY);
        const eqj_vec = Astronomy.VectorFromSphere(eqj_sphere, time);

        if (correct_aberration) {
            // Use non-relativistic approximation: add barycentric Earth velocity vector
            // to the light ray vector. The direction of the light vector points toward the star,
            // which is opposite to the direction light actually travels.
            // The result is the aberration-corrected apparent position of the start in EQJ.
            const eqj_earth = Astronomy.BaryState(Astronomy.Body.Earth, time);
            eqj_vec.x += eqj_earth.vx;
            eqj_vec.y += eqj_earth.vy;
            eqj_vec.z += eqj_earth.vz;
        }

        // Use the rotation matrix to re-orient the EQJ vector to a EQD vector.
        const eqd_vec = Astronomy.RotateVector(rot, eqj_vec);

        // Convert the EQD vector back to spherical angular coordinates.
        const eqd_sphere = Astronomy.SphereFromVector(eqd_vec);

        // Print out the converted angular coordinates.
        const eqd_ra = eqd_sphere.lon / 15;     // convert degrees to sidereal hours
        const eqd_dec = eqd_sphere.lat;
        console.log(`Equator-of-date coordinates: RA=${eqd_ra.toFixed(6)}, DEC=${eqd_dec.toFixed(6)}`);

        // Success!
        process.exit(0);
    }
}

Demo();
