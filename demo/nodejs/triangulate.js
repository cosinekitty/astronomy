/*
    triangulate.js  -  by Don Cross - 2021-06-22

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy
*/

const UsageText = `
USAGE:  node triangulate.js  lat1 lon1 elv1 az1 alt1  lat2 lon2 elv2 az2 alt2

Calculate the best-fit location of a point as observed
from two different locations on or near the Earth's surface.

lat1, lat2 = Geographic latitudes in degrees north of the equator.
lon1, lon2 = Geographic longitudes in degrees east of the prime meridian.
elv1, elv2 = Elevations above sea level in meters.
az1,  az2  = Azimuths toward observed object in degrees clockwise from north.
alt1, alt2 = Altitude angles toward observed object in degrees above horizon.

This program extrapolates lines in the given directions from the two
geographic locations and finds the location in space where they
come closest to intersecting. It then prints out the coordinates
of that triangulation point, along with the error radius in meters.
`;

const Astronomy = require('./astronomy.js');

function ParseNumber(name, text) {
    const x = Number(text);
    if (!Number.isFinite(x)) {
        console.error(`ERROR: Not a valid numeric value for ${name}: "${text}"`);
        process.exit(1);
    }
    return x;
}


function DotProduct(a, b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}


function AddScale(sa, va, sb, vb) {
    return new Astronomy.Vector(
        sa*va.x + sb*vb.x,
        sa*va.y + sb*vb.y,
        sa*va.z + sb*vb.z,
        va.t
    );
}


function DirectionVector(time, observer, altitude, azimuth) {
    // Convert horizontal angles to a horizontal unit vector.
    const hor = new Astronomy.Spherical(altitude, azimuth, 1.0);
    const hvec = Astronomy.VectorFromHorizon(hor, time, null);

    // Find the rotation matrix that converts horizontal vectors to equatorial vectors.
    const rot = Astronomy.Rotation_HOR_EQD(time, observer);

    // Rotate the horizontal (HOR) vector to an equator-of-date (EQD) vector.
    const evec = Astronomy.RotateVector(rot, hvec);

    return evec;
}


function Intersect(pos1, dir1, pos2, dir2) {
    const F = DotProduct(dir1, dir2);
    const amb = AddScale(+1, pos1, -1, pos2);    // amb = pos1 - pos2
    const E = DotProduct(dir1, amb);
    const G = DotProduct(dir2, amb);
    const denom = 1 - F*F;
    if (denom == 0.0) {
        console.error('ERROR: Cannot solve because directions are parallel.');
        process.exit(1);
    }

    const u = (F*G - E) / denom;
    const v = G + F*u;
    if (u < 0.0 || v < 0.0) {
        console.error('ERROR: Lines of sight do not converge.');
        process.exit(1);
    }

    const a = AddScale(1, pos1, u, dir1);     //  a = pos1 + u*dir1
    const b = AddScale(1, pos2, v, dir2);     //  b = pos2 + v*dir2
    const c = AddScale(0.5, a, 0.5, b);       //  c = (a+b)/2
    const miss = AddScale(+1, a, -1, b);      //  miss = a-b

    const dist = (Astronomy.KM_PER_AU * 1000 / 2) * miss.Length();   // error radius in meters
    const obs = Astronomy.VectorObserver(c, true);

    console.log(`Solution: lat = ${obs.latitude.toFixed(6)}, lon = ${obs.longitude.toFixed(6)}, elv = ${obs.height.toFixed(3)} meters; error = ${dist.toFixed(3)} meters`);
}


function Demo() {
    if (process.argv.length === 12) {
        // Validate and parse command line arguments.
        const lat1 = ParseNumber("lat1", process.argv[ 2]);
        const lon1 = ParseNumber("lon1", process.argv[ 3]);
        const elv1 = ParseNumber("elv1", process.argv[ 4]);
        const  az1 = ParseNumber("az1",  process.argv[ 5]);
        const alt1 = ParseNumber("alt1", process.argv[ 6]);
        const lat2 = ParseNumber("lat2", process.argv[ 7]);
        const lon2 = ParseNumber("lon2", process.argv[ 8]);
        const elv2 = ParseNumber("elv2", process.argv[ 9]);
        const  az2 = ParseNumber("az2",  process.argv[10]);
        const alt2 = ParseNumber("alt2", process.argv[11]);

        const obs1 = new Astronomy.Observer(lat1, lon1, elv1);
        const obs2 = new Astronomy.Observer(lat2, lon2, elv2);

        // Use an arbitrary but consistent time for the Earth's rotation.
        const time = Astronomy.MakeTime(0.0);

        // Convert geographic coordinates of the observers to vectors.
        const pos1 = Astronomy.ObserverVector(time, obs1, true);
        const pos2 = Astronomy.ObserverVector(time, obs2, true);

        // Convert horizontal coordinates into unit direction vectors.
        const dir1 = DirectionVector(time, obs1, alt1, az1);
        const dir2 = DirectionVector(time, obs2, alt2, az2);

        // Find the closest point between the skew lines.
        Intersect(pos1, dir1, pos2, dir2);
        process.exit(0);
    } else {
        console.log(UsageText);
        process.exit(1);
    }
}

Demo();
