/*
    culminate.js  -  by Don Cross - 2019-06-17

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This example program shows how to calculate the time
    the Sun, Moon, and planets will next reach their highest point in the sky
    as seen by an observer at a given location on the Earth.
    This is called culmination, and is found by finding when
    each body's "hour angle" is 0.

    Having an hour angle of 0 is another way of saying that the body is
    crossing the meridian, the imaginary semicircle in the sky that passes
    from due north on the horizon, through the zenith (straight up),
    toward due south on the horizon. At this moment the body appears to
    have an azimuth of either 180 degrees (due south) or 0 (due north).

    To execute, run the command:
    node culminate latitude longitude [date]
*/

const Astronomy = require('astronomy.js');

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
        const latitude  = ParseNumber(process.argv[2]);
        const longitude = ParseNumber(process.argv[3]);
        const observer = Astronomy.MakeObserver(latitude, longitude, 0);
        const date = (process.argv.length === 5) ? ParseDate(process.argv[4]) : new Date();
        console.log('search   : ' + date.toISOString());

        for (let body of Astronomy.Bodies) {
            if (body !== 'Earth' && body !== 'EMB' && body !== 'SSB') {
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
