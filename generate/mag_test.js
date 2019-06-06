'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

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

        // [ Date__(UT)__HR:MN      APmag  S-brt               r        rdot            delta      deldot    S-T-O]
        // [ 2023-Mar-30 00:00      -4.01   1.17  0.719092953368  -0.1186373 1.20453495004726 -11.0204917  55.9004]
        let m = line.match(/^\s(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$/);
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
                sbrt: parseFloat(m[7]),
                helio_dist: parseFloat(m[8]),
                helio_radvel: parseFloat(m[9]),
                geo_dist: parseFloat(m[10]),
                geo_radvel: parseFloat(m[11]),
                phase_angle: parseFloat(m[12])
            };

            rows.push(item);
        }
    }

    //console.log(`${filename} : ${rows.length} rows`);

    return {
        filename: filename,
        rows: rows
    };
}

function CheckMagnitudeData(body, data) {
    let diff_lo, diff_hi;
    let sum_squared_diff = 0;
    for (let item of data.rows) {
        let illum = Astronomy.Illumination(body, item.date);
        let diff = illum.mag - item.mag;
        sum_squared_diff += diff*diff;
        if (diff_lo === undefined) {
            diff_lo = diff_hi = diff;
        } else {
            diff_lo = Math.min(diff_lo, diff);
            diff_hi = Math.max(diff_hi, diff);
        }
    }
    let rms = Math.sqrt(sum_squared_diff / data.rows.length);
    const limit = 0.012;
    const pass = (Math.abs(diff_lo) < limit && Math.abs(diff_hi) < limit);
    console.log(`${body.padEnd(8)} ${pass?"    ":"FAIL"}  diff_lo=${diff_lo.toFixed(4).padStart(8)}, diff_hi=${diff_hi.toFixed(4).padStart(8)}, rms=${rms.toFixed(4).padStart(8)}`);
    return pass;
}

function TestSaturn() {
    // JPL Horizons does not include Saturn's rings in its magnitude models.
    // I still don't have authoritative test data for Saturn's magnitude.
    // For now, I just test for consistency with Paul Schlyter's formulas at:
    // http://www.stjarnhimlen.se/comp/ppcomp.html#15

    let success = true;

    const data = [
        { date: '1972-01-01T00:00Z', mag: -0.31904865,  tilt: +24.50061220 },
        { date: '1980-01-01T00:00Z', mag: +0.85213663,  tilt:  -1.85761461 },
        { date: '2009-09-04T00:00Z', mag: +1.01626809,  tilt:  +0.08380716 },
        { date: '2017-06-15T00:00Z', mag: -0.12318790,  tilt: -26.60871409 },
        { date: '2019-05-01T00:00Z', mag: +0.32954097,  tilt: -23.53880802 },
        { date: '2025-09-25T00:00Z', mag: +0.51286575,  tilt:  +1.52327932 },
        { date: '2032-05-15T00:00Z', mag: -0.04652109,  tilt: +26.95717765 }
    ];

    for (let item of data) {
        let illum = Astronomy.Illumination('Saturn', new Date(item.date));
        console.log(`Saturn: date=${illum.time.date.toISOString()}  mag=${illum.mag.toFixed(8).padStart(12)}  ring_tilt=${illum.ring_tilt.toFixed(8).padStart(12)}`);
        const mag_diff = Math.abs(illum.mag - item.mag);
        if (mag_diff > 1.0e-8) {
            console.log(`ERROR: Excessive magnitude error ${mag_diff}`);
            success = false;
        }
        const tilt_diff = Math.abs(illum.ring_tilt - item.tilt);
        if (tilt_diff > 1.0e-8) {
            console.log(`ERROR: Excessive ring tilt error ${tilt_diff}`);
            success = false;
        }
    }

    return success;
}

function TestMaxMag(filename, body) {
    // Test that we can find maximum magnitude events for Venus within
    // ranges found using JPL Horizons ephemeris data that has been
    // pre-processed by magnitude/findmax.py.

    console.log('TestMaxMag: entering');
    const text = fs.readFileSync(filename, 'utf8');
    const lines = text.trim().split(/[\r\n]+/);
    let date = new Date(Date.UTC(2000, 0, 1));
    let max_diff = 0;
    for (let line of lines) {
        //console.log(`TestMaxMag: searching after ${date.toISOString()}`);
        let token = line.split(/\s+/);
        let date1 = new Date(token[0]);
        let date2 = new Date(token[1]);
        let evt = Astronomy.SearchPeakMagnitude(body, date);
        if (evt.time.date < date1 || evt.time.date > date2)
            throw `Event time ${evt.time.toString()} is outside the range ${date1.toISOString()} .. ${date2.toISOString()}`;

        // How close are we to the center date?
        let date_center = new Date((date1.getTime() + date2.getTime())/2);
        let diff_hours = Math.abs(evt.time.date - date_center) / (1000 * 3600);
        if (diff_hours > 7.1)
            throw `Excessive diff_hours = ${diff_hours} from center date ${date_center.toISOString()}`;

        max_diff = Math.max(max_diff, diff_hours);
        date = date2;
    }
    console.log(`TestMaxMag: ${lines.length} events, max error = ${max_diff.toFixed(3)} hours.`);
    return true;
}

function Test() {
    let all_passed = true;
    for (let body of Astronomy.Bodies) {
        if (body !== 'Earth' && body !== 'Saturn') {
            const data = LoadMagnitudeData(`magnitude/${body}.txt`);
            if (!CheckMagnitudeData(body, data))
                all_passed = false;
        }
    }

    if (!TestSaturn())
        all_passed = false;

    if (!TestMaxMag('magnitude/maxmag_Venus.txt', 'Venus'))
        all_passed = false;

    all_passed || Fail('Found excessive error in at least one test.');
}

Test();
console.log('mag_test: success');
process.exit(0);