/*
    solar_time.js  -  by Don Cross - 2023-02-12

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the true solar time for
    a given observer and UTC time.

    To execute, run the command:

        node solar_time.js latitude longitude [date]

    where

        latitude = geographic latitude of the observer (-90 to +90).
        longitude = geographic longitude of the observer (-180 to +180).
        date = optional date and time string.

    If date is omitted, this program uses the computer's current date and time.
    If date is present, date is any string that Node.js can parse as a date and time,
    for example the ISO 8601 UTC format "yyyy-mm-ddThh:mm:ssZ".
*/

const Astronomy = require('./astronomy.js');

function f(x, n) {
    let s = x.toFixed(0);
    while (s.length < n)
        s = '0' + s;
    return s;
}

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

function Demo() {
    if (process.argv.length === 4 || process.argv.length === 5) {
        const latitude  = ParseNumber(process.argv[2], 'latitude');
        const longitude = ParseNumber(process.argv[3], 'longitude');
        const observer = new Astronomy.Observer(latitude, longitude, 0);
        const date = (process.argv.length === 5) ? ParseDate(process.argv[4]) : new Date();

        const hourAngle = Astronomy.HourAngle(Astronomy.Body.Sun, date, observer);
        const solarTimeHours = (hourAngle + 12) % 24;

        let milli = Math.round(solarTimeHours * 3.6e+6);
        let second = 0 | (milli / 1000);
        milli %= 1000;
        let minute = 0 | (second / 60);
        second %= 60;
        let hour = 0 | (minute / 60);
        minute %= 60;
        hour %= 24;

        console.log(`True solar time = ${solarTimeHours.toFixed(4)} hours (${f(hour,2)}:${f(minute,2)}:${f(second,2)}.${f(milli,3)})`);
        process.exit(0);
    }
    console.log('USAGE: node positions.js latitude longitude [date]');
    process.exit(1);
}

Demo();
