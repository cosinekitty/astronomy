/*
    culminate.js  -  by Don Cross - 2019-06-17

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This example program shows how to calculate the time
    the Sun, Moon, and planets will next reach their highest point in the sky
    as seen by an observer at a given location on the Earth.
    This is called culmination, and is found by finding when
    each object's "hour angle" is 0.

    To execute, run the command:
    node culminate latitude longitude [date]
*/

const Astronomy = require('../../source/js/astronomy.js');

function DisplayEvent(name, evt) {
    let text ;
    if (evt) {
        text = evt.time.date.toISOString() + '  altitude=' + evt.hor.altitude.toFixed(2).padStart(6);
        text +=  '  azimuth=' + evt.hor.azimuth.toFixed(2).padStart(7);
    } else {
        text = '(not found)';
    }
    console.log(name.padEnd(8) + ' : ' + text);
}

function Demo() {
    if (process.argv.length === 4 || process.argv.length === 5) {
        const latitude = parseFloat(process.argv[2]);
        const longitude = parseFloat(process.argv[3]);
        const observer = Astronomy.MakeObserver(latitude, longitude, 0);
        const date = (process.argv.length === 5) ? new Date(process.argv[4]) : new Date();
        console.log('search   : ' + date.toISOString());

        for (let body of Astronomy.Bodies) {
            if (body !== 'Earth') {
                let culm = Astronomy.SearchHourAngle(body, observer, 0, date);
                DisplayEvent(body, culm);
            }
        }

        process.exit(0);
    } else {
        console.log('USAGE: node culminate.js latitude longitude [date]');
        process.exit(1);
    }
}

Demo();
