/*
    positions.js  -  by Don Cross - 2019-06-12

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the equatorial and horizontal coordinates
    of the Sun, Moon, and planets, as seen by an observer at a specified
    location on the Earth.

    To execute, run the command:
    
        node positions.js latitude longitude [date]

    where

        latitude = geographic latitude of the observer (-90 to +90).
        longitude = geographic longitude of the observer (-180 to +180).
        date = optional date and time string.

    If date is omitted, this program uses the computer's current date and time.
    If date is present, date is any string that Node.js can parse as a date and time,
    for example the ISO 8601 UTC format "yyyy-mm-ddThh:mm:ssZ".
*/

const Astronomy = require('../../source/js/astronomy.js');      // adjust path as needed for your system

function Format(x) {
    return x.toFixed(2).padStart(8);
}

function Demo() {
    if (process.argv.length === 4 || process.argv.length === 5) {
        const latitude = parseFloat(process.argv[2]);
        const longitude = parseFloat(process.argv[3]);
        const observer = Astronomy.MakeObserver(latitude, longitude, 0);
        const date = (process.argv.length === 5) ? new Date(process.argv[4]) : new Date();
        console.log(`UTC date = ${date.toISOString()}`);
        console.log('');
        console.log(`${'BODY'.padEnd(8)} ${'RA'.padStart(8)} ${'DEC'.padStart(8)} ${'AZ'.padStart(8)} ${'ALT'.padStart(8)}`);
        for (let body of ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']) {
            let equ_2000 = Astronomy.Equator(body, date, observer, false, true);
            let equ_ofdate = Astronomy.Equator(body, date, observer, true, true);
            let hor = Astronomy.Horizon(date, observer, equ_ofdate.ra, equ_ofdate.dec, 'normal');
            console.log(`${body.padEnd(8)} ${Format(equ_2000.ra)} ${Format(equ_2000.dec)} ${Format(hor.azimuth)} ${Format(hor.altitude)}`);
        }
        process.exit(0);
    }
    console.log('USAGE: node positions.js latitude longitude [date]');
    process.exit(1);
}

Demo();
