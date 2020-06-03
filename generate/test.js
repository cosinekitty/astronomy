'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');
let Verbose = false;

function AstroCheck() {
    var date = Astronomy.MakeTime(new Date('1700-01-01T00:00:00Z'));
    var stop = Astronomy.MakeTime(new Date('2200-01-01T00:00:00Z'));
    var body, pos, hor, dt, j2000, ofdate, time;
    const observer = Astronomy.MakeObserver(29, -81, 10);

    console.log(`o ${observer.latitude.toFixed(6)} ${observer.longitude.toFixed(6)} ${observer.height.toFixed(6)}`);

    dt = (10 + Math.PI/100);       // 10.03141592... days; exercise different times of day
    while (date.tt < stop.tt) {
        time = Astronomy.MakeTime(date);

        for (body of Astronomy.Bodies) {
            if (body !== 'Moon') {
                pos = Astronomy.HelioVector(body, date);
                console.log(`v ${body} ${pos.t.tt.toFixed(16)} ${pos.x.toFixed(16)} ${pos.y.toFixed(16)} ${pos.z.toFixed(16)}`);

                if (body !== 'Earth' && body !== 'EMB' && body !== 'SSB') {
                    j2000 = Astronomy.Equator(body, date, observer, false, false);
                    ofdate = Astronomy.Equator(body, date, observer, true, true);
                    hor = Astronomy.Horizon(date, observer, ofdate.ra, ofdate.dec);
                    console.log(`s ${body} ${time.tt.toFixed(16)} ${time.ut.toFixed(16)} ${j2000.ra.toFixed(16)} ${j2000.dec.toFixed(16)} ${j2000.dist.toFixed(16)} ${hor.azimuth.toFixed(16)} ${hor.altitude.toFixed(16)}`);
                }
            }
        }
        pos = Astronomy.GeoMoon(date);
        console.log(`v GM ${pos.t.tt.toFixed(16)} ${pos.x.toFixed(16)} ${pos.y.toFixed(16)} ${pos.z.toFixed(16)}`);

        j2000 = Astronomy.Equator('Moon', date, observer, false, false);
        ofdate = Astronomy.Equator('Moon', date, observer, true, true);
        hor = Astronomy.Horizon(date, observer, ofdate.ra, ofdate.dec);
        console.log(`s GM ${time.tt.toFixed(16)} ${time.ut.toFixed(16)} ${j2000.ra.toFixed(16)} ${j2000.dec.toFixed(16)} ${j2000.dist.toFixed(16)} ${hor.azimuth.toFixed(16)} ${hor.altitude.toFixed(16)}`);

        date = date.AddDays(dt);
    }

    return 0;
}


function MoonPhase() {
    function LoadMoonPhaseData(filename) {
        // Load known moon phase times from US Naval Observatory.
        const text = fs.readFileSync(filename, {encoding:'utf8'});
        const lines = text.trimRight().split('\n');
        let data = [];
        for (let row of lines) {
            let token = row.split(' ');
            data.push({quarter:parseInt(token[0]), date:new Date(token[1])});
        }
        return data;
    }

    function TestLongitudes(data) {
        // Using known moon phase times from US Naval Obs
        let max_arcmin = 0;
        for (let row of data) {
            let elong = Astronomy.MoonPhase(row.date);
            let expected_elong = 90 * row.quarter;
            let degree_error = Math.abs(elong - expected_elong);
            if (degree_error > 180) degree_error = 360 - degree_error;
            let arcmin = 60 * degree_error;
            max_arcmin = Math.max(max_arcmin, arcmin);
        }
        console.log(`JS TestLongitudes: nrows = ${data.length}, max_arcmin = ${max_arcmin}`);
        if (max_arcmin > 1.0) {
            console.log(`JS TestLongitudes: EXCESSIVE ANGULAR ERROR`);
            return 1;
        }
        return 0;
    }

    function SearchYear(year, data, index) {
        const millis_per_minute = 60*1000;
        const threshold_minutes = 2;    // max tolerable prediction error in minutes
        let date = new Date(Date.UTC(year, 0, 1));
        let maxdiff = 0;
        let count = 0;
        let mq = Astronomy.SearchMoonQuarter(date);
        while (index < data.length && data[index].date.getUTCFullYear() === year) {

            // Verify that we found anything at all.
            if (!mq) {
                console.log(`JS SearchYear: could not find next moon quarter after ${date.toISOString()}`);
                return 1;
            }

            // Verify that we found the correct quarter.
            if (data[index].quarter !== mq.quarter) {
                console.log(`JS SearchYear: predicted quarter ${mq.quarter} but correct is ${data[index].quarter} for ${date.toISOString()}`);
                return 1;
            }

            // Verify that the date and time we found is very close to the correct answer.
            // Calculate the discrepancy in minutes.
            // This is appropriate because the "correct" answers are only given to the minute.
            let diff = Math.abs(mq.time.date - data[index].date) / millis_per_minute;
            if (diff > threshold_minutes) {
                console.log(`JS SearchYear: EXCESSIVE ERROR = ${diff.toFixed(3)} minutes, correct=${data[index].date.toISOString()}, calculated=${mq.time.toString()}`);
                return 1;
            }
            maxdiff = Math.max(maxdiff, diff);

            ++index;
            ++count;
            mq = Astronomy.NextMoonQuarter(mq);
            date = mq.time.date;
        }
        if (Verbose) console.log(`JS SearchYear(${year}): count=${count}, maxdiff=${maxdiff.toFixed(3)}`);
        return 0;
    }

    function TestSearch(data) {
        // Search each year finding each quarter moon phase and confirm
        // they match up with the correct answers stored in the 'data' array.
        let index = 0;  // index into 'data' for the current year we are interested in.
        for (let year=1800; year <= 2100; year += 10) {
            while (data[index].date.getUTCFullYear() < year) ++index;
            let error = SearchYear(year, data, index);
            if (error) return error;
        }
        return 0;
    }

    function Statistics(data) {
        // Figure out mean, min, and max differences between various moon quarters.
        const stats = [null, null, null, null];
        let prev = null;
        for (let curr of data) {
            if (prev && curr.date.getUTCFullYear() === prev.date.getUTCFullYear()) {
                let dt = (curr.date - prev.date) / (24 * 3600 * 1000);
                let s = stats[prev.quarter];
                if (s) {
                    s.min = Math.min(s.min, dt);
                    s.max = Math.max(s.max, dt);
                    s.sum += dt;
                    ++s.count;
                } else {
                    stats[prev.quarter] = {
                        min: dt,
                        max: dt,
                        sum: dt,
                        count: 1
                    };
                }
            }
            prev = curr;
        }

        for (let q=0; q < 4; ++q) {
            let s = stats[q];
            if (s) {
                console.log(`JS Statistics: q=${q} min=${s.min.toFixed(3)} avg=${(s.sum/s.count).toFixed(3)} max=${s.max.toFixed(3)} span=${(s.max-s.min).toFixed(3)}`);
            }
        }
    }

    const TestData = LoadMoonPhaseData('moonphase/moonphases.txt');
    Statistics(TestData);
    return TestLongitudes(TestData) || TestSearch(TestData);
}


function LunarApsis() {
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

    console.log(`lunar_apsis_test.js: verified ${count} lines, max time error = ${max_minute_error.toFixed(3)} minutes, max dist error = ${max_dist_error.toFixed(3)} km.`);

    if (count !== 2651)
        throw 'FATAL: Did not process the expected number of data rows!';

    console.log(`LunarApsis PERFORMANCE: time=${elapsed.toFixed(3)}`);
    return 0;
}


function LunarEclipse() {
    Astronomy.CalcMoonCount = 0;
    const filename = 'eclipse/lunar_eclipse.txt';
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.trim().split(/\r?\n/);
    const diff_limit = 2.0;
    let lnum = 0;
    let skip_count = 0;
    let max_diff_minutes = 0;
    let sum_diff_minutes = 0;
    let diff_count = 0;
    let eclipse = Astronomy.SearchLunarEclipse(new Date(Date.UTC(1701, 0)));
    for (let line of lines) {
        ++lnum;
        if (line.length < 17) {
            console.error(`JS TestLunarEclipse(${filename} line ${lnum}): line is too short.`);
            return 1;
        }
        const time_text = line.substr(0, 17);
        const peak_time = Astronomy.MakeTime(new Date(time_text));
        const token = line.substr(17).trim().split(/\s+/);
        if (token.length !== 2) {
            console.error(`JS TestLunarEclipse(${filename} line ${lnum}): invalid number of tokens.`);
            return 1;
        }
        const partial_minutes = parseFloat(token[0]);
        const total_minutes = parseFloat(token[1]);

        let valid = false;
        switch (eclipse.kind) {
        case 'penumbral':
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial == 0.0) && (eclipse.sd_total == 0.0);
            break;

        case 'partial':
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total == 0.0);
            break;

        case 'total':
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total > 0.0);
            break;

        default:
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): invalid eclipse kind ${eclipse.kind}.`);
            return 1;
        }

        if (!valid) {
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): invalid semidurations for kind ${kind}: penum=${sd_penum}, partial=${sd_partial}, total=${sd_total}.`);
            return 1;
        }

        // Check eclipse peak.
        let diff_days = eclipse.peak.ut - peak_time.ut;

        // Tolerate missing penumbral eclipses - skip to next input line without calculating next eclipse.
        if (partial_minutes == 0.0 && diff_days > 20.0) {
            ++skip_count;
            continue;
        }

        let diff_minutes = (24.0 * 60.0) * Math.abs(diff_days);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit) {
            console.error(`JS LunarEclipseTest expected peak: ${peak_time}`);
            console.error(`JS LunarEclipseTest found    peak: ${eclipse.peak}`);
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): EXCESSIVE peak time error = ${diff_minutes} minutes (${diff_days} days).`);
            return 1;
        }

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* check partial eclipse duration */

        diff_minutes = Math.abs(partial_minutes - eclipse.sd_partial);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit) {
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): EXCESSIVE partial eclipse semiduration error: ${diff_minutes} minutes.`);
            return 1;
        }

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* check total eclipse duration */

        diff_minutes = Math.abs(total_minutes - eclipse.sd_total);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit) {
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): EXCESSIVE total eclipse semiduration error: ${diff_minutes} minutes.`);
            return 1;
        }

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* calculate for next iteration */

        eclipse = Astronomy.NextLunarEclipse(eclipse.peak);
    }
    console.log(`JS LunarEclipseTest: PASS (verified ${lnum}, skipped ${skip_count}, max_diff_minutes = ${max_diff_minutes}, avg_diff_minutes = ${sum_diff_minutes / diff_count}, moon calcs = ${Astronomy.CalcMoonCount})`);
    return 0;
}


function main() {
    let args = process.argv.slice(2);
    if (args.length > 0 && args[0] === '-v') {
        Verbose = true;
        args = args.slice(1);
    }

    if (args.length === 1 && args[0] === 'astro_check') {
        return AstroCheck();
    }

    if (args.length === 1 && args[0] === 'moon_phase') {
        return MoonPhase();
    }

    if (args.length === 1 && args[0] === 'lunar_apsis') {
        return LunarApsis();
    }

    if (args.length === 1 && args[0] === 'lunar_eclipse') {
        return LunarEclipse();
    }

    console.log('test.js: Invalid command line arguments.');
    return 1;
}


process.exit(main());
