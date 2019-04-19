'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function Fail(message) {
    console.log(`FATAL(rise_set_test.js): ${message}`);
    process.exit(1);
}

function LoadTestData(filename) {
    // Moon  150 -45 2050-03-07T19:13Z s
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.trimRight().split('\n');
    let data = [];
    for (let row of lines) {        
        let token = row.split(/\s+/g);
        data.push({
            body: token[0],
            lon: parseFloat(token[1]),
            lat: parseFloat(token[2]),
            date: new Date(token[3]),
            direction: { r:+1, s:-1 }[token[4]] || Fail(`Invalid event code ${token[4]}`)
        });
    }
    return data;
}

function Test() {
    const data = LoadTestData('riseset/riseset.txt');

    // The test data is sorted by body, then geographic location, then date/time.

    let body;
    let observer;
    let r_search_date, r_date;
    let s_search_date, s_date;    
    let a_date, b_date, a_dir, b_dir;
    for (let evt of data) {
        if (!observer || observer.latitude !== evt.lat || observer.longitude !== evt.lon || body !== evt.body) {
            // Every time we see a new geographic location, start a new iteration
            // of finding all rise/set times for that UTC calendar year.
            body = evt.body;
            observer = Astronomy.MakeObserver(evt.lat, evt.lon, 0);
            r_search_date = s_search_date = new Date(Date.UTC(evt.date.getUTCFullYear(), 0, 1));
            b_date = null;
        }

        if (b_date) {
            // recycle the second event from the previous iteration as the first event
            a_date = b_date;
            a_dir = b_dir;
            b_date = null;
        } else {
            r_date = Astronomy.SearchRiseSet(body, observer, +1, r_search_date, 366) ||
                Fail(`Did not find ${body} rise after ${r_search_date.toISOString()}`);

            s_date = Astronomy.SearchRiseSet(body, observer, -1, s_search_date, 366) ||
                Fail(`Did not find ${body} set after ${s_search_date.toISOString()}`);

            // Expect the current event to match the earlier of the found dates.
            if (r_date.tt < s_date.tt) {
                a_date = r_date;
                b_date = s_date;
                a_dir = +1;
                b_dir = -1;
            } else {
                a_date = s_date;
                b_date = r_date;
                a_dir = -1;
                b_dir = +1;
            }

            r_search_date = r_date.AddDays(0.01).date;
            s_search_date = s_date.AddDays(0.01).date;
        }

        if (a_dir !== evt.direction) {
            Fail(`Expected ${body} dir=${evt.direction} at ${evt.date.toISOString()} but found ${a_dir} ${a_date.toString()}`);
        }
    }
}

Test();
process.exit(0);