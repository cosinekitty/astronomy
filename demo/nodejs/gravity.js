/*
    gravity.js  -  by Don Cross  -  2021-07-19

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy
*/

const UsageText = `
    USAGE:

    node gravity.js latitude height

    Calculates the gravitational acceleration experienced
    by an observer on the surface of the Earth at the specified
    latitude (degrees north of the equator) and height
    (meters of above sea level).
    The output is the gravitational acceleration in m/s^2.
`;

const Astronomy = require('./astronomy.js');

function ParseNumber(name, text, minValue, maxValue) {
    const x = Number(text);
    if (!Number.isFinite(x) || (x < minValue) || (x > maxValue)) {
        console.error(`ERROR: Not a valid numeric value for ${name}: "${text}".`);
        console.error(`Must be in the range ${minValue} .. ${maxValue}.`);
        process.exit(1);
    }
    return x;
}

function Format(value, width, precision) {
    let s = value.toFixed(precision);
    while (s.length < width) {
        s = ' ' + s;
    }
    return s;
}

function Demo() {
    if (process.argv.length !== 4) {
        console.log(UsageText);
        process.exit(1);
    } else {
        const latitude = ParseNumber('latitude', process.argv[2], -90, +90);
        const height = ParseNumber('height', process.argv[3], 0, 100000);
        const gravity = Astronomy.ObserverGravity(latitude, height);
        console.log(`latitude = ${Format(latitude, 8, 4)},  height = ${Format(height, 6, 0)},  gravity = ${Format(gravity, 8, 6)}`);
        process.exit(0);
    }
}

Demo();
