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
    if (TestLongitudes(TestData)) return 1;
    if (TestSearch(TestData)) return 1;
    console.log('JS MoonPhase: PASS');
    return 0;
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

    console.log(`JS LunarApsis PASS: time=${elapsed.toFixed(3)}`);
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


function PlanetApsis() {
    const Planet = {
        Mercury: { OrbitalPeriod:    87.969 },
        Venus:   { OrbitalPeriod:   224.701 },
        Earth:   { OrbitalPeriod:   365.256 },
        Mars:    { OrbitalPeriod:   686.980 },
        Jupiter: { OrbitalPeriod:  4332.589 },
        Saturn:  { OrbitalPeriod: 10759.22  },
        Uranus:  { OrbitalPeriod: 30685.4   },
        Neptune: { OrbitalPeriod: 60189.0   },
        Pluto:   { OrbitalPeriod: 90560.0   }
    };
    const degree_threshold = 0.1
    const start_time = Astronomy.MakeTime(new Date('1700-01-01T00:00:00Z'));
    let found_bad_planet = false;
    let pindex = 0;
    for (let body of ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']) {
        let count = 1;
        const period = Planet[body].OrbitalPeriod;
        let min_interval = -1.0;
        let max_interval = -1.0;
        let max_diff_days = 0.0;
        let max_dist_ratio = 0.0;
        let apsis = Astronomy.SearchPlanetApsis(body, start_time);
        const filename = `apsides/apsis_${pindex}.txt`;
        const text = fs.readFileSync(filename, {encoding:'utf8'});
        const lines = text.split(/\r?\n/);
        for (const line of lines) {
            if (line.trim() == '') {
                continue;
            }
            const token = line.split(/\s+/);
            if (token.length != 3) {
                throw `${filename} line ${count}: Invalid data format: ${token.length} tokens.`;
            }
            const expected_kind = parseInt(token[0]);
            const expected_time = Astronomy.MakeTime(new Date(token[1]));
            const expected_distance = parseFloat(token[2]);
            if (expected_kind !== apsis.kind) {
                throw `${filename} line ${count}: WRONG APSIS KIND: expected ${expected_kind}, found ${apsis.kind}`;
            }
            const diff_days = Math.abs(expected_time.tt - apsis.time.tt);
            max_diff_days = Math.max(max_diff_days, diff_days);
            const diff_degrees = (diff_days / period) * 360;
            if (diff_degrees > degree_threshold) {
                found_bad_planet = true;
            }
            const diff_dist_ratio = Math.abs(expected_distance - apsis.dist_au) / expected_distance;
            max_dist_ratio = Math.max(max_dist_ratio, diff_dist_ratio);
            if (diff_dist_ratio > 1.0e-4) {
                throw `${filename} line ${count}: distance ratio ${diff_dist_ratio} is too large.`;
            }

            // Calculate the next apsis
            const prev_time = apsis.time;
            try {
                apsis = Astronomy.NextPlanetApsis(body, apsis);
            } catch (e) {
                if (body === 'Pluto') {
                    // It is OK for Pluto to fail due to the time being out-of-bounds.
                    break;
                }
                throw e;
            }
            ++count;
            const interval = apsis.time.tt - prev_time.tt;
            if (min_interval < 0.0) {
                min_interval = max_interval = interval;
            } else {
                min_interval = Math.min(min_interval, interval);
                max_interval = Math.max(max_interval, interval);
            }
        }
        if (count < 2) {
            throw `Failed to find apsides for ${body}`;
        }
        if (Verbose) console.log(`JS PlanetApsis: ${count} apsides for ${body} -- intervals: min=${min_interval}, max=${max_interval}, ratio=${max_interval/min_interval}; max day=${max_diff_days}, degrees=${(max_diff_days / period) * 360}, dist ratio=${max_dist_ratio}`);
        ++pindex;
    }
    if (found_bad_planet) {
        throw `APSIS FAIL: Planet(s) exceeded angular threshold (${degree_threshold} degrees).`;
    }
    return 0;
}


function Elongation() {
    function SearchElongTest() {
        // Max elongation data obtained from:
        // http://www.skycaramba.com/greatest_elongations.shtml
        TestMaxElong('Mercury', '2010-01-17T05:22Z', '2010-01-27T05:22Z', 24.80, 'morning');
        TestMaxElong('Mercury', '2010-05-16T02:15Z', '2010-05-26T02:15Z', 25.10, 'morning');
        TestMaxElong('Mercury', '2010-09-09T17:24Z', '2010-09-19T17:24Z', 17.90, 'morning');
        TestMaxElong('Mercury', '2010-12-30T14:33Z', '2011-01-09T14:33Z', 23.30, 'morning');
        TestMaxElong('Mercury', '2011-04-27T19:03Z', '2011-05-07T19:03Z', 26.60, 'morning');
        TestMaxElong('Mercury', '2011-08-24T05:52Z', '2011-09-03T05:52Z', 18.10, 'morning');
        TestMaxElong('Mercury', '2011-12-13T02:56Z', '2011-12-23T02:56Z', 21.80, 'morning');
        TestMaxElong('Mercury', '2012-04-08T17:22Z', '2012-04-18T17:22Z', 27.50, 'morning');
        TestMaxElong('Mercury', '2012-08-06T12:04Z', '2012-08-16T12:04Z', 18.70, 'morning');
        TestMaxElong('Mercury', '2012-11-24T22:55Z', '2012-12-04T22:55Z', 20.60, 'morning');
        TestMaxElong('Mercury', '2013-03-21T22:02Z', '2013-03-31T22:02Z', 27.80, 'morning');
        TestMaxElong('Mercury', '2013-07-20T08:51Z', '2013-07-30T08:51Z', 19.60, 'morning');
        TestMaxElong('Mercury', '2013-11-08T02:28Z', '2013-11-18T02:28Z', 19.50, 'morning');
        TestMaxElong('Mercury', '2014-03-04T06:38Z', '2014-03-14T06:38Z', 27.60, 'morning');
        TestMaxElong('Mercury', '2014-07-02T18:22Z', '2014-07-12T18:22Z', 20.90, 'morning');
        TestMaxElong('Mercury', '2014-10-22T12:36Z', '2014-11-01T12:36Z', 18.70, 'morning');
        TestMaxElong('Mercury', '2015-02-14T16:20Z', '2015-02-24T16:20Z', 26.70, 'morning');
        TestMaxElong('Mercury', '2015-06-14T17:10Z', '2015-06-24T17:10Z', 22.50, 'morning');
        TestMaxElong('Mercury', '2015-10-06T03:20Z', '2015-10-16T03:20Z', 18.10, 'morning');
        TestMaxElong('Mercury', '2016-01-28T01:22Z', '2016-02-07T01:22Z', 25.60, 'morning');
        TestMaxElong('Mercury', '2016-05-26T08:45Z', '2016-06-05T08:45Z', 24.20, 'morning');
        TestMaxElong('Mercury', '2016-09-18T19:27Z', '2016-09-28T19:27Z', 17.90, 'morning');
        TestMaxElong('Mercury', '2017-01-09T09:42Z', '2017-01-19T09:42Z', 24.10, 'morning');
        TestMaxElong('Mercury', '2017-05-07T23:19Z', '2017-05-17T23:19Z', 25.80, 'morning');
        TestMaxElong('Mercury', '2017-09-02T10:14Z', '2017-09-12T10:14Z', 17.90, 'morning');
        TestMaxElong('Mercury', '2017-12-22T19:48Z', '2018-01-01T19:48Z', 22.70, 'morning');
        TestMaxElong('Mercury', '2018-04-19T18:17Z', '2018-04-29T18:17Z', 27.00, 'morning');
        TestMaxElong('Mercury', '2018-08-16T20:35Z', '2018-08-26T20:35Z', 18.30, 'morning');
        TestMaxElong('Mercury', '2018-12-05T11:34Z', '2018-12-15T11:34Z', 21.30, 'morning');
        TestMaxElong('Mercury', '2019-04-01T19:40Z', '2019-04-11T19:40Z', 27.70, 'morning');
        TestMaxElong('Mercury', '2019-07-30T23:08Z', '2019-08-09T23:08Z', 19.00, 'morning');
        TestMaxElong('Mercury', '2019-11-18T10:31Z', '2019-11-28T10:31Z', 20.10, 'morning');
        TestMaxElong('Mercury', '2010-03-29T23:32Z', '2010-04-08T23:32Z', 19.40, 'evening');
        TestMaxElong('Mercury', '2010-07-28T01:03Z', '2010-08-07T01:03Z', 27.40, 'evening');
        TestMaxElong('Mercury', '2010-11-21T15:42Z', '2010-12-01T15:42Z', 21.50, 'evening');
        TestMaxElong('Mercury', '2011-03-13T01:07Z', '2011-03-23T01:07Z', 18.60, 'evening');
        TestMaxElong('Mercury', '2011-07-10T04:56Z', '2011-07-20T04:56Z', 26.80, 'evening');
        TestMaxElong('Mercury', '2011-11-04T08:40Z', '2011-11-14T08:40Z', 22.70, 'evening');
        TestMaxElong('Mercury', '2012-02-24T09:39Z', '2012-03-05T09:39Z', 18.20, 'evening');
        TestMaxElong('Mercury', '2012-06-21T02:00Z', '2012-07-01T02:00Z', 25.70, 'evening');
        TestMaxElong('Mercury', '2012-10-16T21:59Z', '2012-10-26T21:59Z', 24.10, 'evening');
        TestMaxElong('Mercury', '2013-02-06T21:24Z', '2013-02-16T21:24Z', 18.10, 'evening');
        TestMaxElong('Mercury', '2013-06-02T16:45Z', '2013-06-12T16:45Z', 24.30, 'evening');
        TestMaxElong('Mercury', '2013-09-29T09:59Z', '2013-10-09T09:59Z', 25.30, 'evening');
        TestMaxElong('Mercury', '2014-01-21T10:00Z', '2014-01-31T10:00Z', 18.40, 'evening');
        TestMaxElong('Mercury', '2014-05-15T07:06Z', '2014-05-25T07:06Z', 22.70, 'evening');
        TestMaxElong('Mercury', '2014-09-11T22:20Z', '2014-09-21T22:20Z', 26.40, 'evening');
        TestMaxElong('Mercury', '2015-01-04T20:26Z', '2015-01-14T20:26Z', 18.90, 'evening');
        TestMaxElong('Mercury', '2015-04-27T04:46Z', '2015-05-07T04:46Z', 21.20, 'evening');
        TestMaxElong('Mercury', '2015-08-25T10:20Z', '2015-09-04T10:20Z', 27.10, 'evening');
        TestMaxElong('Mercury', '2015-12-19T03:11Z', '2015-12-29T03:11Z', 19.70, 'evening');
        TestMaxElong('Mercury', '2016-04-08T14:00Z', '2016-04-18T14:00Z', 19.90, 'evening');
        TestMaxElong('Mercury', '2016-08-06T21:24Z', '2016-08-16T21:24Z', 27.40, 'evening');
        TestMaxElong('Mercury', '2016-12-01T04:36Z', '2016-12-11T04:36Z', 20.80, 'evening');
        TestMaxElong('Mercury', '2017-03-22T10:24Z', '2017-04-01T10:24Z', 19.00, 'evening');
        TestMaxElong('Mercury', '2017-07-20T04:34Z', '2017-07-30T04:34Z', 27.20, 'evening');
        TestMaxElong('Mercury', '2017-11-14T00:32Z', '2017-11-24T00:32Z', 22.00, 'evening');
        TestMaxElong('Mercury', '2018-03-05T15:07Z', '2018-03-15T15:07Z', 18.40, 'evening');
        TestMaxElong('Mercury', '2018-07-02T05:24Z', '2018-07-12T05:24Z', 26.40, 'evening');
        TestMaxElong('Mercury', '2018-10-27T15:25Z', '2018-11-06T15:25Z', 23.30, 'evening');
        TestMaxElong('Mercury', '2019-02-17T01:23Z', '2019-02-27T01:23Z', 18.10, 'evening');
        TestMaxElong('Mercury', '2019-06-13T23:14Z', '2019-06-23T23:14Z', 25.20, 'evening');
        TestMaxElong('Mercury', '2019-10-10T04:00Z', '2019-10-20T04:00Z', 24.60, 'evening');
        TestMaxElong('Venus', '2010-12-29T15:57Z', '2011-01-08T15:57Z', 47.00, 'morning');
        TestMaxElong('Venus', '2012-08-05T08:59Z', '2012-08-15T08:59Z', 45.80, 'morning');
        TestMaxElong('Venus', '2014-03-12T19:25Z', '2014-03-22T19:25Z', 46.60, 'morning');
        TestMaxElong('Venus', '2015-10-16T06:57Z', '2015-10-26T06:57Z', 46.40, 'morning');
        TestMaxElong('Venus', '2017-05-24T13:09Z', '2017-06-03T13:09Z', 45.90, 'morning');
        TestMaxElong('Venus', '2018-12-27T04:24Z', '2019-01-06T04:24Z', 47.00, 'morning');
        TestMaxElong('Venus', '2010-08-10T03:19Z', '2010-08-20T03:19Z', 46.00, 'evening');
        TestMaxElong('Venus', '2012-03-17T08:03Z', '2012-03-27T08:03Z', 46.00, 'evening');
        TestMaxElong('Venus', '2013-10-22T08:00Z', '2013-11-01T08:00Z', 47.10, 'evening');
        TestMaxElong('Venus', '2015-05-27T18:46Z', '2015-06-06T18:46Z', 45.40, 'evening');
        TestMaxElong('Venus', '2017-01-02T13:19Z', '2017-01-12T13:19Z', 47.10, 'evening');
        TestMaxElong('Venus', '2018-08-07T17:02Z', '2018-08-17T17:02Z', 45.90, 'evening');
    }

    function LoadData(filename) {
        const text = fs.readFileSync(filename, {encoding:'utf8'});
        const lines = text.trimRight().split('\n');
        let data = [];
        for (let row of lines) {
            let token = row.split(/\s+/);
            data.push({date:new Date(token[0]), body:token[1]});
        }
        return data;
    }

    function TestFile(filename, startSearchYear, targetRelLon) {
        const data = LoadData(filename);
        const startDate = new Date(Date.UTC(startSearchYear, 0, 1));
        for (let item of data) {
            let time = Astronomy.SearchRelativeLongitude(item.body, targetRelLon, startDate);
            let diff_minutes = (time.date - item.date) / 60000;
            if (Verbose) console.log(`JS ${item.body}: error = ${diff_minutes.toFixed(3)} minutes`);
            if (Math.abs(diff_minutes) > 15)
                throw `!!! Excessive error for body ${item.body}`;
        }
    }

    function TestPlanet(outFileName, body, startYear, stopYear, zeroLonEventName) {
        let rlon = 0;
        let date = new Date(Date.UTC(startYear, 0, 1));
        let stopDate = new Date(Date.UTC(stopYear, 0, 1));
        let text = '';
        let count = 0;
        let prev_time, min_diff, max_diff, sum_diff=0;

        while (date < stopDate) {
            let event = (rlon === 0) ? zeroLonEventName : 'sup';
            let evt_time = Astronomy.SearchRelativeLongitude(body, rlon, date);
            if (prev_time) {
                // Check for consistent intervals.
                // Mainly I don't want to accidentally skip over an event!
                let day_diff = evt_time.tt - prev_time.tt;
                if (min_diff === undefined) {
                    min_diff = max_diff = day_diff;
                } else {
                    min_diff = Math.min(min_diff, day_diff);
                    max_diff = Math.max(max_diff, day_diff);
                }
                sum_diff += day_diff;
            }
            let geo = Astronomy.GeoVector(body, evt_time);
            let dist = Math.sqrt(geo.x*geo.x + geo.y*geo.y + geo.z*geo.z);
            text += `e ${body} ${event} ${evt_time.tt} ${dist}\n`;
            rlon = 180 - rlon;
            date = evt_time.date;
            ++count;
            prev_time = evt_time;
        }

        fs.writeFileSync(outFileName, text);

        const ratio = max_diff / min_diff;
        if (Verbose) console.log(`JS TestPlanet(${body}): ${count} events, interval min=${min_diff.toFixed(1)}, max=${max_diff.toFixed(1)}, avg=${(sum_diff/count).toFixed(1)}, ratio=${ratio.toFixed(3)}`);

        let thresh = {Mercury:1.65, Mars:1.30}[body] || 1.07;
        if (ratio > thresh)
            throw `TestPlanet: Excessive event interval ratio for ${body} = ${ratio}`;
    }

    function TestMaxElong(body, startText, verifyText, verifyAngle, verifyVisibility) {
        let startDate = new Date(startText);
        let verifyDate = new Date(verifyText);
        let evt = Astronomy.SearchMaxElongation(body, startDate);

        let hour_diff = (verifyDate - evt.time.date) / (1000 * 3600);
        let arcmin_diff = 60.0 * Math.abs(evt.elongation - verifyAngle);
        if (Verbose) console.log(`JS TestMaxElong: ${body.padStart(8)} ${evt.visibility.padStart(8)} elong=${evt.elongation.toFixed(2).padStart(5)} (${arcmin_diff.toFixed(2).padStart(4)} arcmin)  ${evt.time.toString()} (err ${hour_diff.toFixed(2).padStart(5)} hours)`);

        if (evt.visibility !== verifyVisibility)
            throw `TestMaxElong: expected visibility ${verifyVisibility}, but found ${evt.visibility}`;

        if (arcmin_diff > 4.0)
            throw `TestMaxElong: excessive angular error = ${angle_diff} arcmin`;

        if (Math.abs(hour_diff) > 0.603)
            throw `TestMaxElong: excessive hour error = ${hour_diff}`;
    }

    TestFile('longitude/opposition_2018.txt', 2018, 0);

    for (let body of ['Mercury', 'Venus'])
        TestPlanet(`temp/js_longitude_${body}.txt`, body, 1700, 2200, 'inf');

    for (let body of ['Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
        TestPlanet(`temp/js_longitude_${body}.txt`, body, 1700, 2200, 'opp');

    SearchElongTest();

    console.log('JS Elongation: SUCCESS');
    return 0;
}


function Seasons() {
    function LoadTestData(filename) {
        // Moon  150 -45 2050-03-07T19:13Z s
        const text = fs.readFileSync(filename, {encoding:'utf8'});
        const lines = text.trimRight().split('\n');
        let data = [];
        let lnum = 0;
        let minByMonth = [];
        let maxByMonth = [];
        for (let row of lines) {
            let token = row.split(/\s+/g);
            let item = {
                lnum: ++lnum,
                date: new Date(token[0]),
                name: token[1]  // Perihelion, Equinox, Solstice, Aphelion
            };
            data.push(item);

            if (item.name === 'Equinox' || item.name == 'Solstice') {
                let month = 1 + item.date.getUTCMonth();
                let format = item.date.toISOString().substring(8);
                if (!minByMonth[month] || format < minByMonth[month])
                    minByMonth[month] = format;
                if (!maxByMonth[month] || format > maxByMonth[month])
                    maxByMonth[month] = format;
            }
        }

        if (Verbose) {
            console.log(`JS Seasons LoadTestData: count = ${data.length}`);
            for (let month of [3, 6, 9, 12]) {
                console.log(`Month ${month}: earliest ${minByMonth[month]}, latest ${maxByMonth[month]}`);
            }
        }
        return data;
    }
    const data = LoadTestData('seasons/seasons.txt');
    let current_year;
    let seasons;
    let calc_date;
    let min_diff, max_diff, sum_diff = 0, count = 0;
    let month_max_diff = [];
    for (let item of data) {
        let year = item.date.getUTCFullYear();
        if (current_year !== year) {
            current_year = year;
            seasons = Astronomy.Seasons(year);
        }
        calc_date = null;
        let month = 1 + item.date.getUTCMonth();
        switch (item.name) {
        case 'Equinox':
            switch (month) {
            case 3:
                calc_date = seasons.mar_equinox.date;
                break;
            case 9:
                calc_date = seasons.sep_equinox.date;
                break;
            default:
                throw `ERROR: Invalid equinox date in test data: ${item.date.toISOString()}`;
            }
            break;
        case 'Solstice':
            switch (month) {
            case 6:
                calc_date = seasons.jun_solstice.date;
                break;
            case 12:
                calc_date = seasons.dec_solstice.date;
                break;
            default:
                throw `ERROR: Invalid solstice date in test data: ${item.date.toISOString()}`;
            }
            break;
        default:
            continue;   // ignore the other kinds of events for now
        }
        if (!calc_date)
            throw `ERROR: Missing calc_date for test date ${item.date.toISOString()}`;
        let diff_minutes = (calc_date - item.date) / 60000;
        if (Math.abs(diff_minutes) > 2.37) {
            throw `ERROR: Excessive error in season calculation: ${diff_minutes.toFixed(3)} minutes`;
        }

        if (min_diff === undefined) {
            min_diff = max_diff = diff_minutes;
        } else {
            min_diff = Math.min(min_diff, diff_minutes);
            max_diff = Math.max(max_diff, diff_minutes);
        }
        if (month_max_diff[month] === undefined) {
            month_max_diff[month] = Math.abs(diff_minutes);
        } else {
            month_max_diff[month] = Math.max(month_max_diff[month], Math.abs(diff_minutes));
        }
        sum_diff += diff_minutes;
        ++count;
    }

    if (Verbose) {
        console.log(`JS Seasons: n=${count}, minute errors: min=${min_diff.toFixed(3)}, avg=${(sum_diff/count).toFixed(3)}, max=${max_diff.toFixed(3)}`);
        for (let month of [3, 6, 9, 12]) {
            console.log(`JS Seasons: max diff by month ${month} = ${month_max_diff[month].toFixed(3)}`);
        }
    }
    console.log('JS Seasons: PASS');
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

    if (args.length === 1 && args[0] === 'planet_apsis') {
        return PlanetApsis();
    }

    if (args.length === 1 && args[0] == 'elongation') {
        return Elongation();
    }

    if (args.length === 1 && args[0] == 'seasons') {
        return Seasons();
    }

    console.log('test.js: Invalid command line arguments.');
    return 1;
}


process.exit(main());
