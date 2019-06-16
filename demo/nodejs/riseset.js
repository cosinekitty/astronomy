/*
    riseset.js  -  by Don Cross - 2019-06-15

    Example Node.js program for Astronomy Engine:
    https://cosinekitty.github.io/astronomy/

    This program calculates the time of the next
    sunrise, sunset, moonrise, and moonset.

    To execute, run the command:
    node riseset latitude longitude [date]
*/

const Astronomy = require('../../source/js/astronomy.js');      // adjust path as needed for your system

function DisplayEvent(name, evt) {
    let text = evt ? evt.date.toISOString() : '';
    console.log(name.padEnd(8) + ' : ' + text);
}

function Demo() {
    if (process.argv.length === 4 || process.argv.length === 5) {
        const latitude = parseFloat(process.argv[2]);
        const longitude = parseFloat(process.argv[3]);
        const observer = Astronomy.MakeObserver(latitude, longitude, 0);
        const date = (process.argv.length === 5) ? new Date(process.argv[4]) : new Date();
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
