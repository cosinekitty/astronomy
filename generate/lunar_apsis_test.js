'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function Test() {
    const filename = 'apsides/moon.txt';
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.split(/\n/);

    let evt = Astronomy.SearchLunarApsis(new Date(Date.UTC(2001, 0, 1)));

    let lnum = 0;
    let count = 0;
    for (let line of lines) {
        ++lnum;

        let token = line.split(/\s+/);
        if (token.length !== 3)
            continue;

        let kind = parseInt(token[0]);
        let date = new Date(token[1]);
        let dist = parseInt(token[2]);

        if (evt.apsisType !== kind) {
            console.log('line = ', line);
            throw `${filename} line ${lnum}: Expected apsis type ${kind}, found ${evt.apsisType}`;
        }

        let diff_minutes = Math.abs(evt.time.date - date) / (1000 * 60);
        if (diff_minutes > 35) {
            throw `${filename} line ${lnum}: Excessive time error: ${diff_minutes} minutes`;
        }

        let diff_dist = Math.abs(evt.dist_km - dist);
        if (diff_dist > 25) {
            throw `${filename} line ${lnum}: Excessive distance error: ${diff_dist} km`;
        }

        ++count;
        evt = Astronomy.SearchLunarApsis(evt.time.AddDays(1));
    }

    console.log(`lunar_apsis_test: successfully verified ${count} lines.`);
}

Test();
process.exit(0);
