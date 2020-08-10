/*
    riseset.js  -  by Don Cross - 2019-06-15

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the time of the next
    sunrise, sunset, moonrise, and moonset.

    To execute, run the command:
    node riseset latitude longitude [date]
*/

const Astronomy = require('astronomy.js');

function DisplayEvent(name, evt) {
    let text = evt ? evt.date.toISOString() : '';
    console.log(name.padEnd(8) + ' : ' + text);
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
        const latitude  = ParseNumber(process.argv[2]);
        const longitude = ParseNumber(process.argv[3]);
        const observer = Astronomy.MakeObserver(latitude, longitude, 0);
        const date = (process.argv.length === 5) ? ParseDate(process.argv[4]) : new Date();
        let sunrise  = Astronomy.SearchRiseSet('Sun',  observer, +1, date, 300);
        let sunset   = Astronomy.SearchRiseSet('Sun',  observer, -1, date, 300);
        let moonrise = Astronomy.SearchRiseSet('Moon', observer, +1, date, 300);
        let moonset  = Astronomy.SearchRiseSet('Moon', observer, -1, date, 300);
        console.log('search   : ' + date.toISOString());
        DisplayEvent('sunrise',  sunrise);
        DisplayEvent('sunset',   sunset);
        DisplayEvent('moonrise', moonrise);
        DisplayEvent('moonset',  moonset);
        process.exit(0);
    } else {
        console.log('USAGE: node riseset.js latitude longitude [date]');
        process.exit(1);
    }
}

Demo();
