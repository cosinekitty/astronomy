/*
    vernal_point.js  -  by Don Cross  -  2023-02-20

    Calculate the movement of the vernal equinox point
    between two different moments in time.
    Given two times, calculate how much the equinox
    point moves from the first time to the second time,
    relative to true ecliptic coordinates (ECT) expressed
    at the first time.
*/

const Astronomy = require('./astronomy.js');


function VernalPointLongitudeChange(time1, time2) {
    console.log(`time1 = ${time1}`);
    console.log(`time2 = ${time2}`);

    // Create a vector pointing toward the vernal point at time2.
    const vec2 = new Astronomy.Vector(1, 0, 0, time2);

    // Find the rotation matrix that converts true ecliptic of date (ECT)
    // coordinates from the second time to the first time.
    // We accomplish this in two rotations: ECT(t2) --> EQJ --> ECT(t1).
    const rot = Astronomy.CombineRotation(
        Astronomy.Rotation_ECT_EQJ(time2),
        Astronomy.Rotation_EQJ_ECT(time1)
    );

    // Apply the rotation matrix to `vec2` to obtain `vec1`: the
    // second time's vernal point expressed in the first time's ecliptic system.
    const vec1 = Astronomy.RotateVector(rot, vec2);

    // Convert ecliptic direction from a vector to angles.
    const sphere = Astronomy.SphereFromVector(vec1);

    return (sphere.lon > 180) ? (360 - sphere.lon) : sphere.lon;
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
    if (process.argv.length === 4) {
        const time1 = Astronomy.MakeTime(ParseDate(process.argv[2]));
        const time2 = Astronomy.MakeTime(ParseDate(process.argv[3]));
        const longitudeChange = VernalPointLongitudeChange(time1, time2);
        console.log(`The vernal point's ecliptic longitude changed by ${longitudeChange.toFixed(4)} degrees.`);
        process.exit(0);
    } else {
        console.log('USAGE: node vernal_point time1 time2');
        console.log('where the times are in the format yyyy-mm-dd or yyyy-mm-ddThh:mm:ssZ');
        process.exit(1);
    }
}


Demo();
