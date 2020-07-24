/*
    lunar_eclipse.js  -  by Don Cross - 2020-05-17

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Searches for the next 10 partial/total lunar eclipses after
    the current date, or a date specified on the command line.

    To execute, run the command:
    node lunar_eclipse [date]
*/

const Astronomy = require('../../source/js/astronomy.js');

function Pad(s, w) {
    s = s.toFixed(0);
    while (s.length < w) {
        s = '0' + s;
    }
    return s;
}

function FormatDate(t) {
    const date = t.date;
    var year = Pad(date.getUTCFullYear(), 4);
    var month = Pad(1 + date.getUTCMonth(), 2);
    var day = Pad(date.getUTCDate(), 2);
    var hour = Pad(date.getUTCHours(), 2);
    var minute = Pad(date.getUTCMinutes(), 2);
    var svalue = date.getUTCSeconds() + (date.getUTCMilliseconds() / 1000);
    var second = Pad(Math.round(svalue), 2);
    return `${year}-${month}-${day} ${hour}:${minute}:${second} UTC`;
}

function PrintEclipse(e) {
    // Calculate beginning/ending of different phases
    // of an eclipse by subtracting/adding the peak time
    // with the number of minutes indicated by the "semi-duration"
    // fields sd_partial and sd_total.
    const MINUTES_PER_DAY = 24 * 60;

    const p1 = e.peak.AddDays(-e.sd_partial / MINUTES_PER_DAY);
    console.log(`${FormatDate(p1)} - Partial eclipse begins.`);

    if (e.sd_total > 0) {
        const t1 = e.peak.AddDays(-e.sd_total / MINUTES_PER_DAY);
        console.log(`${FormatDate(t1)} - Total eclipse begins.`);
    }

    console.log(`${FormatDate(e.peak)} - Peak of ${e.kind} eclipse.`);

    if (e.sd_total > 0) {
        const t2 = e.peak.AddDays(+e.sd_total / MINUTES_PER_DAY);
        console.log(`${FormatDate(t2)} - Total eclipse ends.`);
    }

    const p2 = e.peak.AddDays(+e.sd_partial / MINUTES_PER_DAY);
    console.log(`${FormatDate(p2)} - Partial eclipse ends.`);
    console.log('');
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
    const date = (process.argv.length === 3) ? ParseDate(process.argv[2]) : new Date();
    let count = 0;
    let eclipse = Astronomy.SearchLunarEclipse(date);
    for(;;) {
        if (eclipse.kind !== 'penumbral') {
            PrintEclipse(eclipse);
            if (++count === 10) {
                break;
            }
        }
        eclipse = Astronomy.NextLunarEclipse(eclipse.peak);
    }
}

Demo();
process.exit(0);
