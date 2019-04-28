'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function Fail(message) {
    console.log(`FATAL(mag_test.js): ${message}`);
    process.exit(1);
}

function LoadMagnitudeData(filename) {
    const text = fs.readFileSync(filename, 'utf8');
    const lines = text.split(/[\r\n]+/);
    let lnum = 0;
    let rows = [];
    for (let line of lines) {
        ++lnum;

        //  Date__(UT)__HR:MN      APmag  S-brt            delta      deldot    S-T-O
        // [ 2016-Mar-28 00:00      -3.83   0.94 1.60152932868679   6.0989077  25.8529]
        let m = line.match(/^\s(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$/);
        if (m) {
            const year = parseInt(m[1]);
            const month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].indexOf(m[2]);
            const day = parseInt(m[3]);
            const hour = parseInt(m[4]);
            const minute = parseInt(m[5]);

            const item = {
                lnum: lnum,
                date: new Date(Date.UTC(year, month, day, hour, minute)),
                mag: parseFloat(m[6]),
                geo_dist: parseFloat(m[8]),
                phase_angle: parseFloat(m[10])
            };

            rows.push(item);
        }
    }

    console.log(`${filename} : ${rows.length} rows`);

    return {
        filename: filename,
        rows: rows
    };
}

function CheckMagnitudeData(body, data) {
    let diff_lo, diff_hi;
    for (let item of data.rows) {
        let illum = Astronomy.Illumination(body, item.date);
        let diff = illum.mag - item.mag;
        if (diff_lo === undefined) {
            diff_lo = diff_hi = diff;
        } else {
            diff_lo = Math.min(diff_lo, diff);
            diff_hi = Math.max(diff_hi, diff);
        }
    }
    console.log(`${body}: diff_lo=${diff_lo}, diff_hi=${diff_hi}`);
    if (Math.abs(diff_lo) < 0.01 && Math.abs(diff_hi) < 0.01) {
        console.log(`${body}: PASS`);
    } else {
        console.log(`${body}: FAIL`);
    }
}

for (let body of Astronomy.Bodies) {
    if (body !== 'Earth') {
        const data = LoadMagnitudeData(`magnitude/${body}.txt`);
        CheckMagnitudeData(body, data);
    }
}