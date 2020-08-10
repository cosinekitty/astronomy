/*
    seasons.js  -  by Don Cross - 2019-06-16

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the time of the next
    sunrise, sunset, moonrise, and moonset.

    To execute, run the command:
    node riseset latitude longitude [date]
*/

const Astronomy = require('./astronomy.js');

function DisplayEvent(name, evt) {
    let text = evt ? evt.date.toISOString() : '';
    console.log(name.padEnd(17) + ' : ' + text);
}

function Demo() {
    if (process.argv.length === 3) {
        let year = parseInt(process.argv[2]);
        if (!Number.isSafeInteger(year)) {
            console.log(`ERROR: Not a valid year: "${process.argv[2]}"`);
            process.exit(1);
        }
        let seasons = Astronomy.Seasons(year);
        DisplayEvent('March equinox', seasons.mar_equinox);
        DisplayEvent('June solstice', seasons.jun_solstice);
        DisplayEvent('September equinox', seasons.sep_equinox);
        DisplayEvent('December solstice', seasons.dec_solstice);
        process.exit(0);
    } else {
        console.log('USAGE: node seasons.js year');
        process.exit(1);
    }
}

Demo();
