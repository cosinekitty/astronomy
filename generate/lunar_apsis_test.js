/*
    lunar_apsis_test.js

    Exercises finding of lunar apogee and perigee using test data derived from:
    http://astropixels.com/ephemeris/moon/moonperap2001.html
*/

'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function Test() {
    const filename = 'apsides/moon.txt';
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.split(/\r?\n/);

    const time_before = new Date();
    let evt = Astronomy.SearchLunarApsis(new Date(Date.UTC(2001, 0, 1)));

    let lnum = 0;
    let count = 0;
    let max_minute_error = 0;
    let max_dist_error = 0;
    for (let line of lines) {
        ++lnum;

        let token = line.split(/\s+/);
        if (token.length !== 3)
            continue;

        let kind = parseInt(token[0]);
        let date = new Date(token[1]);
        let dist = parseInt(token[2]);

        if (evt.kind !== kind) {
            console.log('line = ', line);
            throw `${filename} line ${lnum}: Expected apsis type ${kind}, found ${evt.kind}`;
        }

        let diff_minutes = Math.abs(evt.time.date - date) / (1000 * 60);
        if (diff_minutes > 35) {
            throw `${filename} line ${lnum}: Excessive time error: ${diff_minutes} minutes`;
        }
        max_minute_error = Math.max(max_minute_error, diff_minutes);

        let diff_dist = Math.abs(evt.dist_km - dist);
        if (diff_dist > 25) {
            throw `${filename} line ${lnum}: Excessive distance error: ${diff_dist} km`;
        }
        max_dist_error = Math.max(max_dist_error, diff_dist);

        ++count;
        evt = Astronomy.NextLunarApsis(evt);
    }
    const time_after = new Date();
    const elapsed = (time_after - time_before) / 1000;

    console.log(`lunar_apsis_test: verified ${count} lines, max time error = ${max_minute_error.toFixed(3)} minutes, max dist error = ${max_dist_error.toFixed(3)} km.`);

    if (count !== 2651)
        throw 'FATAL: Did not process the expected number of data rows!';

    const perf = Astronomy.GetPerformanceMetrics();

    console.log(`lunar_apsis_test PERFORMANCE: time=${elapsed.toFixed(3)}, iter=${perf.lunar_apsis_iter}, calcmoon=${perf.calcmoon}, calcmoon/call=${perf.calcmoon/perf.lunar_apsis_calls}`);
}

Test();
process.exit(0);
