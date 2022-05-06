'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');
let Verbose = false;

function Fail(message) {
    console.trace();
    console.log(`FATAL(test.js): ${message}`);
    process.exit(1);
}

function v(x) {
    // Verify that 'x' is really a number.

    if (!Number.isFinite(x)) {
        console.trace();
        throw `Not a finite number: ${x}`;
    }

    return x;
}

function int(s) {
    return v(parseInt(s));
}

function float(s) {
    return v(parseFloat(s));
}

function abs(x) {
    return Math.abs(v(x));
}

function max(a, b) {
    return Math.max(v(a), v(b));
}

function min(a, b) {
    return Math.min(v(a), v(b));
}

function sin(x) {
    return Math.sin(v(x));
}

function cos(x) {
    return Math.cos(v(x));
}

function sqrt(x) {
    return v(Math.sqrt(v(x)));
}

function ReadLines(filename) {
    // Load the entire text of the given file.
    const text = fs.readFileSync(filename, {encoding:'utf8'});

    // Split it into lines, handling the different line endings
    // in Linux and Windows.
    // FIXFIXFIX: This might not work on Mac OS or other operating systems.
    const lines = text.trimEnd().split(/\r?\n/);
    return lines;
}

function SelectJupiterMoon(jm, mindex) {
    return [jm.io, jm.europa, jm.ganymede, jm.callisto][mindex] ||
        Fail(`SelectJupiterMoon: invalid mindex = ${mindex}`);
}

function AstroCheck() {
    var date = Astronomy.MakeTime(new Date('1700-01-01T00:00:00Z'));
    var stop = Astronomy.MakeTime(new Date('2200-01-01T00:00:00Z'));
    var body, pos, hor, dt, j2000, ofdate, time;
    const observer = new Astronomy.Observer(29, -81, 10);

    console.log(`o ${observer.latitude.toFixed(6)} ${observer.longitude.toFixed(6)} ${observer.height.toFixed(6)}`);

    dt = (10 + Math.PI/100);       // 10.03141592... days; exercise different times of day
    while (date.tt < stop.tt) {
        time = Astronomy.MakeTime(date);

        for (body in Astronomy.Body) {
            if (body !== 'Moon') {
                pos = Astronomy.HelioVector(body, date);
                console.log(`v ${body} ${pos.t.tt.toExponential(18)} ${pos.x.toExponential(18)} ${pos.y.toExponential(18)} ${pos.z.toExponential(18)}`);

                if (body !== 'Earth' && body !== 'EMB' && body !== 'SSB') {
                    j2000 = Astronomy.Equator(body, date, observer, false, false);
                    ofdate = Astronomy.Equator(body, date, observer, true, true);
                    hor = Astronomy.Horizon(date, observer, ofdate.ra, ofdate.dec);
                    console.log(`s ${body} ${time.tt.toExponential(18)} ${time.ut.toExponential(18)} ${j2000.ra.toExponential(18)} ${j2000.dec.toExponential(18)} ${j2000.dist.toExponential(18)} ${hor.azimuth.toExponential(18)} ${hor.altitude.toExponential(18)}`);
                }
            }
        }
        pos = Astronomy.GeoMoon(date);
        console.log(`v GM ${pos.t.tt.toExponential(18)} ${pos.x.toExponential(18)} ${pos.y.toExponential(18)} ${pos.z.toExponential(18)}`);

        j2000 = Astronomy.Equator('Moon', date, observer, false, false);
        ofdate = Astronomy.Equator('Moon', date, observer, true, true);
        hor = Astronomy.Horizon(date, observer, ofdate.ra, ofdate.dec);
        console.log(`s GM ${time.tt.toExponential(18)} ${time.ut.toExponential(18)} ${j2000.ra.toExponential(18)} ${j2000.dec.toExponential(18)} ${j2000.dist.toExponential(18)} ${hor.azimuth.toExponential(18)} ${hor.altitude.toExponential(18)}`);

        const jm = Astronomy.JupiterMoons(time);
        for (let mindex = 0; mindex < 4; ++mindex) {
            const moon = SelectJupiterMoon(jm, mindex);
            console.log(`j ${mindex} ${time.tt.toExponential(18)} ${time.ut.toExponential(18)} ${moon.x.toExponential(18)} ${moon.y.toExponential(18)} ${moon.z.toExponential(18)} ${moon.vx.toExponential(18)} ${moon.vy.toExponential(18)} ${moon.vz.toExponential(18)}`);
        }

        date = date.AddDays(dt);
    }

    return 0;
}


function MoonPhase() {
    function LoadMoonPhaseData(filename) {
        // Load known moon phase times from US Naval Observatory.
        const lines = ReadLines(filename);
        let data = [];
        for (let row of lines) {
            let token = row.split(' ');
            data.push({quarter:int(token[0]), date:new Date(token[1])});
        }
        return data;
    }

    function TestLongitudes(data) {
        let max_arcmin = 0;
        for (let row of data) {
            let elong = v(Astronomy.MoonPhase(row.date));
            let expected_elong = 90 * v(row.quarter);
            let degree_error = abs(elong - expected_elong);
            if (degree_error > 180) degree_error = 360 - degree_error;
            let arcmin = 60 * degree_error;
            max_arcmin = max(max_arcmin, arcmin);
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
        const threshold_minutes = 1.5;    // max tolerable prediction error in minutes
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
            let diff = abs(mq.time.date - data[index].date) / millis_per_minute;
            if (diff > threshold_minutes) {
                console.log(`JS SearchYear: EXCESSIVE ERROR = ${diff.toFixed(3)} minutes, correct=${data[index].date.toISOString()}, calculated=${mq.time.toString()}`);
                return 1;
            }
            maxdiff = max(maxdiff, diff);

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
                    s.min = min(s.min, dt);
                    s.max = max(s.max, dt);
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
                console.log(`JS MoonPhase statistics: q=${q} min=${s.min.toFixed(3)} avg=${(s.sum/s.count).toFixed(3)} max=${s.max.toFixed(3)} span=${(s.max-s.min).toFixed(3)}`);
            }
        }
    }

    const TestData = LoadMoonPhaseData('moonphase/moonphases.txt');
    if (Verbose) {
        Statistics(TestData);
    }
    if (TestLongitudes(TestData)) return 1;
    if (TestSearch(TestData)) return 1;
    console.log('JS MoonPhase: PASS');
    return 0;
}


function LunarApsis() {
    const filename = 'apsides/moon.txt';
    const lines = ReadLines(filename);

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

        let kind = int(token[0]);
        let date = new Date(token[1]);
        let dist = int(token[2]);

        if (evt.kind !== kind) {
            console.log('line = ', line);
            throw `${filename} line ${lnum}: Expected apsis type ${kind}, found ${evt.kind}`;
        }

        let diff_minutes = abs(evt.time.date - date) / (1000 * 60);
        if (diff_minutes > 35) {
            throw `${filename} line ${lnum}: Excessive time error: ${diff_minutes} minutes`;
        }
        max_minute_error = max(max_minute_error, diff_minutes);

        let diff_dist = abs(evt.dist_km - dist);
        if (diff_dist > 25) {
            throw `${filename} line ${lnum}: Excessive distance error: ${diff_dist} km`;
        }
        max_dist_error = max(max_dist_error, diff_dist);

        ++count;
        evt = Astronomy.NextLunarApsis(evt);
    }
    const time_after = new Date();
    const elapsed = (time_after - time_before) / 1000;

    if (Verbose) {
        console.log(`JS LunarApsis: verified ${count} lines, max time error = ${max_minute_error.toFixed(3)} minutes, max dist error = ${max_dist_error.toFixed(3)} km.`);
    }

    if (count !== 2651)
        throw 'FATAL: Did not process the expected number of data rows!';

    console.log(`JS LunarApsis PASS: time=${elapsed.toFixed(3)}`);
    return 0;
}


function LunarEclipseIssue78() {
    // https://github.com/cosinekitty/astronomy/issues/78

    let eclipse = Astronomy.SearchLunarEclipse(new Date(Date.UTC(2020, 11, 19)));
    const expected_peak = new Date('2021-05-26T11:18:42Z');  // https://www.timeanddate.com/eclipse/lunar/2021-may-26
    const dt = (expected_peak - eclipse.peak.date) / 1000;
    if (abs(dt) > 40.0)
        throw `LunarEclipseIssue78: Excessive prediction error = ${dt} seconds.`;
    if (eclipse.kind !== Astronomy.EclipseKind.Total)
        throw `Expected total eclipse; found: ${eclipse.kind}`;
    console.log(`JS LunarEclipseIssue78: PASS`);
    return 0;
}


function LunarEclipse() {
    Astronomy.CalcMoonCount = 0;
    const filename = 'eclipse/lunar_eclipse.txt';
    const lines = ReadLines(filename);
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
        const partial_minutes = float(token[0]);
        const total_minutes = float(token[1]);

        let valid = false;
        switch (eclipse.kind) {
        case Astronomy.EclipseKind.Penumbral:
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial == 0.0) && (eclipse.sd_total == 0.0);
            break;

        case Astronomy.EclipseKind.Partial:
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total == 0.0);
            break;

        case Astronomy.EclipseKind.Total:
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
        let diff_days = v(eclipse.peak.ut) - v(peak_time.ut);

        // Tolerate missing penumbral eclipses - skip to next input line without calculating next eclipse.
        if (partial_minutes == 0.0 && diff_days > 20.0) {
            ++skip_count;
            continue;
        }

        let diff_minutes = (24.0 * 60.0) * abs(diff_days);
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

        diff_minutes = abs(partial_minutes - eclipse.sd_partial);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit) {
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): EXCESSIVE partial eclipse semiduration error: ${diff_minutes} minutes.`);
            return 1;
        }

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* check total eclipse duration */

        diff_minutes = abs(total_minutes - eclipse.sd_total);
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

function VectorFromAngles(lat, lon) {
    const latrad = Astronomy.DEG2RAD * lat;
    const lonrad = Astronomy.DEG2RAD * lon;
    const coslat = cos(latrad);
    return [
        cos(lonrad) * coslat,
        sin(lonrad) * coslat,
        sin(latrad)
    ];
}


function AngleDiff(alat, alon, blat, blon) {
    const a = VectorFromAngles(alat, alon);
    const b = VectorFromAngles(blat, blon);
    const dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    if (dot <= -1.0) {
        return 180.0;
    }
    if (dot >= +1.0) {
        return 0.0;
    }
    return v(Astronomy.RAD2DEG * Math.acos(dot));
}


function GlobalSolarEclipse() {
    const expected_count = 1180;
    const filename = 'eclipse/solar_eclipse.txt';
    const lines = ReadLines(filename);
    let max_minutes = 0.0;
    let max_angle = 0.0;
    let skip_count = 0;
    let eclipse = Astronomy.SearchGlobalSolarEclipse(new Date(Date.UTC(1701, 0)));
    let lnum = 0;
    for (let line of lines) {
        ++lnum;
        // 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
        let token = line.trim().split(/\s+/);
        if (token.length !== 5) {
            console.error(`JS GlobalSolarEclipse(${filename} line ${lnum}): invalid token count = ${token.length}`);
            return 1;
        }
        const peak = Astronomy.MakeTime(new Date(token[0]));
        const typeChar = token[2];
        const lat = float(token[3]);
        const lon = float(token[4]);
        const expected_kind = {
            'P': Astronomy.EclipseKind.Partial,
            'A': Astronomy.EclipseKind.Annular,
            'T': Astronomy.EclipseKind.Total,
            'H': Astronomy.EclipseKind.Total
        }[typeChar];

        let diff_days = eclipse.peak.tt - peak.tt;
        // Sometimes we find marginal eclipses that aren't listed in the test data.
        // Ignore them if the distance between the Sun/Moon shadow axis and the Earth's center is large.
        while (diff_days < -25.0 && eclipse.distance > 9000.0) {
            ++skip_count;
            eclipse = Astronomy.NextGlobalSolarEclipse(eclipse.peak);
            diff_days = eclipse.peak.ut - peak.ut;
        }

        // Validate the eclipse prediction.
        const diff_minutes = (24 * 60) * abs(diff_days);
        if (diff_minutes > 6.93) {
            console.error(`JS GlobalSolarEclipse(${filename} line ${lnum}): EXCESSIVE TIME ERROR = ${diff_minutes} minutes`);
            return 1;
        }

        if (diff_minutes > max_minutes) {
            max_minutes = diff_minutes;
        }

        // Validate the eclipse kind, but only when it is not a "glancing" eclipse.
        if ((eclipse.distance < 6360) && (eclipse.kind != expected_kind)) {
            console.error(`JS GlobalSolarEclipse(${filename} line ${lnum}): WRONG ECLIPSE KIND: expected ${expected_kind}, found ${eclipse.kind}`);
            return 1;
        }

        if (eclipse.kind === Astronomy.EclipseKind.Total || eclipse.kind === Astronomy.EclipseKind.Annular) {
            // When the distance between the Moon's shadow ray and the Earth's center is beyond 6100 km,
            // it creates a glancing blow whose geographic coordinates are excessively sensitive to
            // slight changes in the ray. Therefore, it is unreasonable to count large errors there.
            if (eclipse.distance < 6100.0) {
                const diff_angle = AngleDiff(lat, lon, eclipse.latitude, eclipse.longitude);
                if (diff_angle > 0.247) {
                    console.error(`JS GlobalSolarEclipse(${filename} line ${lnum}): EXCESSIVE GEOGRAPHIC LOCATION ERROR = ${diff_angle} degrees`);
                    console.log(`   calculated lat=${eclipse.latitude}, lon=${eclipse.longitude}; expected ${lat}, ${lon}`);
                    return 1;
                }
                if (diff_angle > max_angle) {
                    max_angle = diff_angle;
                }
            }
        }

        eclipse = Astronomy.NextGlobalSolarEclipse(eclipse.peak);
    }

    if (lnum != expected_count) {
        console.error(`JS GlobalSolarEclipse: WRONG LINE COUNT = ${lnum}, expected ${expected_count}`);
        return 1;
    }

    if (skip_count > 2) {
        console.error(`JS GlobalSolarEclipse: EXCESSIVE SKIP COUNT = ${skip_count}`);
        return 1;
    }

    console.log(`JS GlobalSolarEclipse: PASS (${lnum} verified, ${skip_count} skipped, max minutes = ${max_minutes}, max angle = ${max_angle})`);
    return 0;
}


function LocalSolarEclipse1() {
    // Re-use the test data for global solar eclipses, only feed the given coordinates
    // into the local solar eclipse predictor as the observer's location.
    // In each case, start the search 20 days before the expected eclipse.
    // Then verify that the peak time and eclipse type is correct in each case.

    const expected_count = 1180;
    const filename = 'eclipse/solar_eclipse.txt';
    const lines = ReadLines(filename);
    let max_minutes = 0.0;
    let skip_count = 0;
    let lnum = 0;
    for (let line of lines) {
        ++lnum;
        // 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
        let token = line.trim().split(/\s+/);
        if (token.length !== 5) {
            console.error(`JS LocalSolarEclipse1(${filename} line ${lnum}): invalid token count = ${token.length}`);
            return 1;
        }
        const peak = Astronomy.MakeTime(new Date(token[0]));
        //const typeChar = token[2];
        const lat = float(token[3]);
        const lon = float(token[4]);
        const observer = new Astronomy.Observer(lat, lon, 0);

        // Start the search 20 days before we know the eclipse should peak.
        const search_start = peak.AddDays(-20);
        const eclipse = Astronomy.SearchLocalSolarEclipse(search_start, observer);

        // Validate the predicted eclipse peak time.
        const diff_days = v(eclipse.peak.time.tt) - v(peak.tt);
        if (diff_days > 20) {
            ++skip_count;
            continue;
        }

        const diff_minutes = (24 * 60) * abs(diff_days);
        if (diff_minutes > 7.14) {
            console.error(`JS LocalSolarEclipse1(${filename} line ${lnum}): EXCESSIVE TIME ERROR = ${diff_minutes} minutes`);
            return 1;
        }

        if (diff_minutes > max_minutes) {
            max_minutes = diff_minutes;
        }

        // Confirm that we can call NextLocalSolarEclipse and it doesn't throw an exception.
        // Could add more stringent checking here later.
        const next = Astronomy.NextLocalSolarEclipse(eclipse.peak.time, observer);
        v(next.peak.time.tt);
    }

    if (lnum != expected_count) {
        console.error(`JS LocalSolarEclipse1: WRONG LINE COUNT = ${lnum}, expected ${expected_count}`);
        return 1;
    }

    if (skip_count > 6) {
        console.error(`JS LocalSolarEclipse1: EXCESSIVE SKIP COUNT = ${skip_count}`);
        return 1;
    }

    console.log(`JS LocalSolarEclipse1: PASS (${lnum} verified, ${skip_count} skipped, max minutes = ${max_minutes})`);
    return 0;
}


function TrimLine(line) {
    // Treat '#' as a comment character.
    const poundIndex = line.indexOf('#');
    if (poundIndex >= 0) {
        line = line.substr(0, poundIndex);
    }
    line = line.trim();
    return line;
}


function ParseEvent(time_str, alt_str, required) {
    if (required) {
        return {
            time: new Astronomy.MakeTime(new Date(time_str)),
            altitude: float(alt_str)
        };
    }
    if (time_str !== '-') {
        throw `Expected event time to be "-" but found "${time_str}"`;
    }
    return null;
}


function LocalSolarEclipse2() {
    // Test ability to calculate local solar eclipse conditions away from
    // the peak position on the Earth.

    const filename = 'eclipse/local_solar_eclipse.txt';
    const lines = ReadLines(filename);
    let lnum = 0;
    let verify_count = 0;
    let max_minutes = 0.0;
    let max_degrees = 0.0;

    function CheckEvent(calc, expect) {
        const diff_minutes = (24 * 60) * abs(expect.time.ut - calc.time.ut);
        if (diff_minutes > max_minutes) {
            max_minutes = diff_minutes;
        }
        if (diff_minutes > 1.0) {
            throw `CheckEvent(${filename} line ${lnum}): EXCESSIVE TIME ERROR: ${diff_minutes} minutes.`;
        }
        const diff_alt = abs(expect.altitude - calc.altitude);
        if (diff_alt > max_degrees) {
            max_degrees = diff_alt;
        }
        if (diff_alt > 0.5) {
            throw `CheckEvent(${filename} line ${lnum}): EXCESSIVE ALTITUDE ERROR: ${diff_alt} degrees.`;
        }
    }

    for (let line of lines) {
        ++lnum;
        line = TrimLine(line);
        if (line === '') continue;
        const token = line.split(/\s+/);
        if (token.length !== 13) {
            console.error(`JS LocalSolarEclipse2(${filename} line ${lnum}): Incorrect token count = ${token.length}`);
            return 1;
        }
        const latitude = float(token[0]);
        const longitude = float(token[1]);
        const observer = new Astronomy.Observer(latitude, longitude, 0);
        const typeChar = token[2];
        const expected_kind = {
            'P': Astronomy.EclipseKind.Partial,
            'A': Astronomy.EclipseKind.Annular,
            'T': Astronomy.EclipseKind.Total,
            'H': Astronomy.EclipseKind.Total
        }[typeChar];
        const p1    = ParseEvent(token[3],  token[4],   true);
        const t1    = ParseEvent(token[5],  token[6],   (typeChar !== 'P'));
        const peak  = ParseEvent(token[7],  token[8],   true);
        const t2    = ParseEvent(token[9],  token[10],  (typeChar !== 'P'));
        const p2    = ParseEvent(token[11], token[12],  true);

        const search_time = p1.time.AddDays(-20);
        const eclipse = Astronomy.SearchLocalSolarEclipse(search_time, observer);
        if (eclipse.kind !== expected_kind) {
            console.error(`JS LocalSolarEclipse2(${filename} line ${lnum}): expected eclipse kind "${expected_kind}" but found "${eclipse.kind}".`);
            return 1;
        }
        CheckEvent(eclipse.peak, peak);
        CheckEvent(eclipse.partial_begin, p1);
        CheckEvent(eclipse.partial_end, p2);
        if (typeChar != 'P') {
            CheckEvent(eclipse.total_begin, t1);
            CheckEvent(eclipse.total_end, t2);
        }
        ++verify_count;
    }

    console.log(`JS LocalSolarEclipse2: PASS (${verify_count} verified, max_minutes = ${max_minutes}, max_degrees = ${max_degrees})`);
    return 0;
}


function LocalSolarEclipse() {
    if (0 !== LocalSolarEclipse1())
        return 1;

    if (0 !== LocalSolarEclipse2())
        return 1;

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
    const start_time = Astronomy.MakeTime(new Date('1700-01-01T00:00:00Z'));
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
        const lines = ReadLines(filename);
        for (const line of lines) {
            if (line.trim() == '') {
                continue;
            }
            const token = line.split(/\s+/);
            if (token.length != 3) {
                throw `${filename} line ${count}: Invalid data format: ${token.length} tokens.`;
            }
            const expected_kind = int(token[0]);
            const expected_time = Astronomy.MakeTime(new Date(token[1]));
            const expected_distance = float(token[2]);
            if (expected_kind !== apsis.kind) {
                throw `${filename} line ${count}: WRONG APSIS KIND: expected ${expected_kind}, found ${apsis.kind}`;
            }
            const diff_days = abs(expected_time.tt - apsis.time.tt);
            max_diff_days = max(max_diff_days, diff_days);
            const diff_degrees = (diff_days / period) * 360;
            const degree_threshold = 0.1;
            if (diff_degrees > 0.1) {
                throw `APSIS FAIL: ${body} exceeded angular threshold (${diff_degrees} vs ${degree_threshold} degrees).`;
            }
            const diff_dist_ratio = abs(expected_distance - apsis.dist_au) / expected_distance;
            max_dist_ratio = max(max_dist_ratio, diff_dist_ratio);
            if (diff_dist_ratio > 1.05e-4) {
                throw `${filename} line ${count}: distance ratio ${diff_dist_ratio} is too large.`;
            }

            // Calculate the next apsis
            const prev_time = apsis.time;
            apsis = Astronomy.NextPlanetApsis(body, apsis);
            ++count;
            const interval = apsis.time.tt - prev_time.tt;
            if (min_interval < 0.0) {
                min_interval = max_interval = interval;
            } else {
                min_interval = min(min_interval, interval);
                max_interval = max(max_interval, interval);
            }
        }
        if (count < 2) {
            throw `Failed to find apsides for ${body}`;
        }
        if (Verbose) console.log(`JS PlanetApsis: ${count} apsides for ${body} -- intervals: min=${min_interval}, max=${max_interval}, ratio=${max_interval/min_interval}; max day=${max_diff_days}, degrees=${(max_diff_days / period) * 360}, dist ratio=${max_dist_ratio}`);
        ++pindex;
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
        const lines = ReadLines(filename);
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
            if (abs(diff_minutes) > 6.8)
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
                    min_diff = min(min_diff, day_diff);
                    max_diff = max(max_diff, day_diff);
                }
                sum_diff += day_diff;
            }
            let geo = Astronomy.GeoVector(body, evt_time, false);
            let dist = sqrt(geo.x*geo.x + geo.y*geo.y + geo.z*geo.z);
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
        let arcmin_diff = 60.0 * abs(evt.elongation - verifyAngle);
        if (Verbose) console.log(`JS TestMaxElong: ${body.padStart(8)} ${evt.visibility.padStart(8)} elong=${evt.elongation.toFixed(2).padStart(5)} (${arcmin_diff.toFixed(2).padStart(4)} arcmin)  ${evt.time.toString()} (err ${hour_diff.toFixed(2).padStart(5)} hours)`);

        if (evt.visibility !== verifyVisibility)
            throw `TestMaxElong: expected visibility ${verifyVisibility}, but found ${evt.visibility}`;

        if (arcmin_diff > 3.4)
            throw `TestMaxElong: excessive angular error = ${angle_diff} arcmin`;

        if (abs(hour_diff) > 0.6)
            throw `TestMaxElong: excessive hour error = ${hour_diff}`;
    }

    TestFile('longitude/opposition_2018.txt', 2018, 0);

    for (let body of ['Mercury', 'Venus'])
        TestPlanet(`temp/js_longitude_${body}.txt`, body, 1700, 2200, 'inf');

    for (let body of ['Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
        TestPlanet(`temp/js_longitude_${body}.txt`, body, 1700, 2200, 'opp');

    SearchElongTest();

    console.log('JS Elongation: PASS');
    return 0;
}


function Seasons() {
    function LoadTestData(filename) {
        // Moon  150 -45 2050-03-07T19:13Z s
        const lines = ReadLines(filename);
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
        if (abs(diff_minutes) > 2.37) {
            throw `ERROR: Excessive error in season calculation: ${diff_minutes.toFixed(3)} minutes`;
        }

        if (min_diff === undefined) {
            min_diff = max_diff = diff_minutes;
        } else {
            min_diff = min(min_diff, diff_minutes);
            max_diff = max(max_diff, diff_minutes);
        }
        if (month_max_diff[month] === undefined) {
            month_max_diff[month] = abs(diff_minutes);
        } else {
            month_max_diff[month] = max(month_max_diff[month], abs(diff_minutes));
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
    return 0;
}


function SeasonsIssue187() {
    // This is a regression test for:
    // https://github.com/cosinekitty/astronomy/issues/187
    // For years far from the present, the seasons search was sometimes failing.

    for (let year = -2000; year <= +9999; ++year) {
        try {
            Astronomy.Seasons(year);
        } catch (e) {
            console.error(`JS SeasonsIssue187: FAIL (year = ${year}): `, e);
            return 1;
        }
    }

    console.log('JS SeasonsIssue187: PASS');
    return 0;
}


function RiseSet() {
    function LoadTestData(filename) {
        // Moon  150 -45 2050-03-07T19:13Z s
        const lines = ReadLines(filename);
        let data = [];
        let lnum = 0;
        for (let row of lines) {
            let token = row.split(/\s+/g);
            data.push({
                lnum: ++lnum,
                body: token[0],
                lon: float(token[1]),
                lat: float(token[2]),
                date: new Date(token[3]),
                direction: { r:+1, s:-1 }[token[4]] || Fail(`Invalid event code ${token[4]}`)
            });
        }
        return data;
    }
    const data = LoadTestData('riseset/riseset.txt');

    // The test data is sorted by body, then geographic location, then date/time.

    const before_date = new Date();

    let body;
    let observer;
    let r_search_date, r_date;
    let s_search_date, s_date;
    let a_date, b_date, a_dir, b_dir;
    let sum_minutes = 0;
    let max_minutes = 0;
    for (let evt of data) {
        if (!observer || observer.latitude !== evt.lat || observer.longitude !== evt.lon || body !== evt.body) {
            // Every time we see a new geographic location, start a new iteration
            // of finding all rise/set times for that UTC calendar year.
            body = evt.body;
            observer = new Astronomy.Observer(evt.lat, evt.lon, 0);
            r_search_date = s_search_date = new Date(Date.UTC(evt.date.getUTCFullYear(), 0, 1));
            b_date = null;
            if (Verbose) {
                console.log(`JS RiseSet: ${body} ${evt.lat} ${evt.lon}`);
            }
        }

        if (b_date) {
            // recycle the second event from the previous iteration as the first event
            a_date = b_date;
            a_dir = b_dir;
            b_date = null;
        } else {
            r_date = Astronomy.SearchRiseSet(body, observer, +1, r_search_date, 366) ||
                Fail(`JS RiseSet: Did not find ${body} rise after ${r_search_date.toISOString()}`);

            s_date = Astronomy.SearchRiseSet(body, observer, -1, s_search_date, 366) ||
                Fail(`JS RiseSet: Did not find ${body} set after ${s_search_date.toISOString()}`);

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
            Fail(`[line ${evt.lnum}] Expected ${body} dir=${evt.direction} at ${evt.date.toISOString()} but found ${a_dir} ${a_date.toString()}`);
        }

        let error_minutes = abs(a_date.date - evt.date) / 60000;
        sum_minutes += error_minutes * error_minutes;
        if (error_minutes > max_minutes) {
            max_minutes = error_minutes;
            if (Verbose) {
                console.log(`Line ${evt.lnum} : error = ${error_minutes.toFixed(4)}`);
            }
        }
        if (error_minutes > 0.57) {
            console.log(`Expected ${evt.date.toISOString()}`);
            console.log(`Found    ${a_date.toString()}`);
            Fail("Excessive prediction time error.");
        }
    }

    const after_date = new Date();
    const elapsed_seconds = (after_date - before_date) / 1000;

    console.log(`JS RiseSet PASS: elapsed=${elapsed_seconds.toFixed(3)}, error in minutes: rms=${sqrt(sum_minutes/data.length).toFixed(4)}, max=${max_minutes.toFixed(4)}`);
    return 0;
}


function Rotation() {
    function CompareMatrices(caller, a, b, tolerance) {
        for (let i=0; i<3; ++i) {
            for (let j=0; j<3; ++j) {
                const diff = abs(a.rot[i][j] - b.rot[i][j]);
                if (diff > tolerance) {
                    throw `ERROR(${caller}): matrix[${i}][${j}] = ${a.rot[i][j]}, expected ${b.rot[i][j]}, diff ${diff}`;
                }
            }
        }
    }

    function CompareVectors(caller, a, b, tolerance) {
        let diff;

        diff = abs(a.x - b.x);
        if (diff > tolerance) {
            throw `ERROR(${caller}): vector x = ${a.x}, expected ${b.x}, diff ${diff}`;
        }

        diff = abs(a.y - b.y);
        if (diff > tolerance) {
            throw `ERROR(${caller}): vector y = ${a.y}, expected ${b.y}, diff ${diff}`;
        }

        diff = abs(a.z - b.z);
        if (diff > tolerance) {
            throw `ERROR(${caller}): vector z = ${a.z}, expected ${b.z}, diff ${diff}`;
        }
    }

    function Rotation_MatrixInverse() {
        const a = Astronomy.MakeRotation([
            [1, 4, 7],
            [2, 5, 8],
            [3, 6, 9]
        ]);

        const v = Astronomy.MakeRotation([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]
        ]);

        const b = Astronomy.InverseRotation(a);
        CompareMatrices('Rotation_MatrixInverse', b, v, 0);
    }

    function Rotation_MatrixMultiply() {
        const a = Astronomy.MakeRotation([
            [1, 4, 7],
            [2, 5, 8],
            [3, 6, 9]
        ]);

        const b = Astronomy.MakeRotation([
            [10, 13, 16],
            [11, 14, 17],
            [12, 15, 18]
        ]);

        const v = Astronomy.MakeRotation([
            [84, 201, 318],
            [90, 216, 342],
            [96, 231, 366]
        ]);

        const c = Astronomy.CombineRotation(b, a);
        CompareMatrices('Rotation_MatrixMultiply', c, v, 0);
    }

    function VectorDiff(a, b) {
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const dz = a.z - b.z;
        return sqrt(dx*dx + dy*dy + dz*dz);
    }

    function Test_EQJ_ECL() {
        const r = Astronomy.Rotation_EQJ_ECL();

        /* Calculate heliocentric Earth position at a test time. */
        const time = Astronomy.MakeTime(new Date('2019-12-08T19:39:15Z'));
        const ev = Astronomy.HelioVector('Earth', time);

        /* Use the existing Astronomy.Ecliptic() to calculate ecliptic vector and angles. */
        const ecl = Astronomy.Ecliptic(ev);
        if (Verbose) console.log(`JS Test_EQJ_ECL ecl = (${ecl.vec.x}, ${ecl.vec.y}, ${ecl.vec.z})`);

        /* Now compute the same vector via rotation matrix. */
        const ee = Astronomy.RotateVector(r, ev);
        const dx = ee.x - ecl.vec.x;
        const dy = ee.y - ecl.vec.y;
        const dz = ee.z - ecl.vec.z;
        const diff = sqrt(dx*dx + dy*dy + dz*dz);
        if (Verbose) console.log(`JS Test_EQJ_ECL ee = (${ee.x}, ${ee.y}, ${ee.z}); diff = ${diff}`);
        if (diff > 2.0e-15)
            throw 'Test_EQJ_ECL: EXCESSIVE VECTOR ERROR';

        /* Reverse the test: go from ecliptic back to equatorial. */
        const ir = Astronomy.Rotation_ECL_EQJ();
        const et = Astronomy.RotateVector(ir, ee);
        const idiff = VectorDiff(et, ev);
        if (Verbose) console.log(`JS Test_EQJ_ECL ev diff = ${idiff}`);
        if (idiff > 2.3e-16)
            throw 'Test_EQJ_ECL: EXCESSIVE REVERSE ROTATION ERROR';
    }

    function Test_GAL_EQJ_NOVAS(filename) {
        const THRESHOLD_SECONDS = 8.8;
        const rot = Astronomy.Rotation_EQJ_GAL();
        const lines = ReadLines(filename);
        const time = new Astronomy.AstroTime(0);    // placeholder time - value does not matter
        let lnum = 0;
        let max_diff = 0;
        for (let row of lines) {
            ++lnum;
            let token = row.trim().split(/\s+/);
            if (token.length !== 4)
                throw `Test_GAL_EQJ_NOVAS(${filename} line ${lnum}): expected 4 tokens, found ${token.length}`;
            const ra = float(token[0]);
            const dec = float(token[1]);
            const glon = float(token[2]);
            const glat = float(token[3]);

            const eqj_sphere = new Astronomy.Spherical(dec, 15*ra, 1);
            const eqj_vec = Astronomy.VectorFromSphere(eqj_sphere, time);
            const gal_vec = Astronomy.RotateVector(rot, eqj_vec);
            const gal_sphere = Astronomy.SphereFromVector(gal_vec);
            const dlat = v(gal_sphere.lat - glat);
            const dlon = cos(Astronomy.DEG2RAD * glat) * v(gal_sphere.lon - glon);
            const diff = 3600 * sqrt(dlon*dlon + dlat*dlat);
            if (diff > THRESHOLD_SECONDS)
                throw `Test_GAL_EQJ_NOVAS(${filename} line ${lnum}): EXCESSIVE ERROR = ${diff.toFixed(3)} arcseconds.`;
            if (diff > max_diff)
                max_diff = diff;
        }
        if (Verbose) console.log(`JS Test_GAL_EQJ_NOVAS: PASS. max_diff = ${max_diff.toFixed(3)} arcseconds.`);
        return 0;
    }

    function Test_EQJ_EQD(body) {
        /* Verify conversion of equatorial J2000 to equatorial of-date, and back. */
        /* Use established functions to calculate spherical coordinates for the body, in both EQJ and EQD. */
        const time = Astronomy.MakeTime(new Date('2019-12-08T20:50:00Z'));
        const observer = new Astronomy.Observer(+35, -85, 0);
        const eq2000 = Astronomy.Equator(body, time, observer, false, true);
        const eqdate = Astronomy.Equator(body, time, observer, true, true);

        /* Convert EQJ spherical coordinates to vector. */
        const v2000 = eq2000.vec;

        /* Find rotation matrix. */
        const r = Astronomy.Rotation_EQJ_EQD(time);

        /* Rotate EQJ vector to EQD vector. */
        const vdate = Astronomy.RotateVector(r, v2000);

        /* Convert vector back to angular equatorial coordinates. */
        let equcheck = Astronomy.EquatorFromVector(vdate);

        /* Compare the result with the eqdate. */
        const ra_diff = abs(equcheck.ra - eqdate.ra);
        const dec_diff = abs(equcheck.dec - eqdate.dec);
        const dist_diff = abs(equcheck.dist - eqdate.dist);
        if (Verbose) console.log(`JS Test_EQJ_EQD: ${body} ra=${eqdate.ra}, dec=${eqdate.dec}, dist=${eqdate.dist}, ra_diff=${ra_diff}, dec_diff=${dec_diff}, dist_diff=${dist_diff}`);
        if (ra_diff > 1.0e-14 || dec_diff > 1.0e-14 || dist_diff > 4.0e-15)
            throw 'Test_EQJ_EQD: EXCESSIVE ERROR';

        /* Perform the inverse conversion back to equatorial J2000 coordinates. */
        const ir = Astronomy.Rotation_EQD_EQJ(time);
        const t2000 = Astronomy.RotateVector(ir, vdate);
        const diff = VectorDiff(t2000, v2000);
        if (Verbose) console.log(`JS Test_EQJ_EQD: ${body} inverse diff = ${diff}`);
        if (diff > 5.0e-15)
            throw 'Test_EQJ_EQD: EXCESSIVE INVERSE ERROR';
    }

    function Test_EQD_HOR(body) {
        /* Use existing functions to calculate horizontal coordinates of the body for the time+observer. */
        const time = Astronomy.MakeTime(new Date('1970-12-13T05:15:00Z'));
        const observer = new Astronomy.Observer(-37, +45, 0);
        const eqd = Astronomy.Equator(body, time, observer, true, true);
        if (Verbose) console.log(`JS Test_EQD_HOR ${body}: OFDATE ra=${eqd.ra}, dec=${eqd.dec}`);
        const hor = Astronomy.Horizon(time, observer, eqd.ra, eqd.dec, 'normal');

        /* Calculate the position of the body as an equatorial vector of date. */
        const vec_eqd = eqd.vec;

        /* Calculate rotation matrix to convert equatorial J2000 vector to horizontal vector. */
        const rot = Astronomy.Rotation_EQD_HOR(time, observer);

        /* Rotate the equator of date vector to a horizontal vector. */
        const vec_hor = Astronomy.RotateVector(rot, vec_eqd);

        /* Convert the horizontal vector to horizontal angular coordinates. */
        const xsphere = Astronomy.HorizonFromVector(vec_hor, 'normal');
        const diff_alt = abs(xsphere.lat - hor.altitude);
        const diff_az = abs(xsphere.lon - hor.azimuth);

        if (Verbose) console.log(`JS Test_EQD_HOR ${body}: trusted alt=${hor.altitude}, az=${hor.azimuth}; test alt=${xsphere.lat}, az=${xsphere.lon}; diff_alt=${diff_alt}, diff_az=${diff_az}`);
        if (diff_alt > 4.3e-14 || diff_az > 1.2e-13)
            throw 'Test_EQD_HOR: EXCESSIVE HORIZONTAL ERROR.';

        /* Confirm that we can convert back to horizontal vector. */
        const check_hor = Astronomy.VectorFromHorizon(xsphere, time, 'normal');
        let diff = VectorDiff(check_hor, vec_hor);
        if (Verbose) console.log(`JS Test_EQD_HOR ${body}: horizontal recovery: diff = ${diff}`);
        if (diff > 2.0e-15)
            throw 'Test_EQD_HOR: EXCESSIVE ERROR IN HORIZONTAL RECOVERY.';

        /* Verify the inverse translation from horizontal vector to equatorial of-date vector. */
        const irot = Astronomy.Rotation_HOR_EQD(time, observer);
        const check_eqd = Astronomy.RotateVector(irot, vec_hor);
        diff = VectorDiff(check_eqd, vec_eqd);
        if (Verbose) console.log(`JS Test_EQD_HOR ${body}: OFDATE inverse rotation diff = ${diff}`);
        if (diff > 2.1e-15)
            throw 'Test_EQD_HOR: EXCESSIVE OFDATE INVERSE HORIZONTAL ERROR.';

        /* Exercise HOR to EQJ translation. */
        const eqj = Astronomy.Equator(body, time, observer, false, true);
        const vec_eqj = eqj.vec;
        const yrot = Astronomy.Rotation_HOR_EQJ(time, observer);
        const check_eqj = Astronomy.RotateVector(yrot, vec_hor);
        diff = VectorDiff(check_eqj, vec_eqj);
        if (Verbose) console.log(`JS Test_EQD_HOR ${body}: J2000 inverse rotation diff = ${diff}`);
        if (diff > 6.0e-15)
            throw 'Test_EQD_HOR: EXCESSIVE J2000 INVERSE HORIZONTAL ERROR.';

        /* Verify the inverse translation: EQJ to HOR. */
        const zrot = Astronomy.Rotation_EQJ_HOR(time, observer);
        const another_hor = Astronomy.RotateVector(zrot, vec_eqj);
        diff = VectorDiff(another_hor, vec_hor);
        if (Verbose) console.log(`JS Test_EQD_HOR ${body}: EQJ inverse rotation diff = ${diff}`);
        if (diff > 3.0e-15)
            throw 'Test_EQD_HOR: EXCESSIVE EQJ INVERSE HORIZONTAL ERROR.';
    }

    const IdentityMatrix = Astronomy.IdentityMatrix();

    function CheckInverse(aname, bname, arot, brot) {
        const crot = Astronomy.CombineRotation(arot, brot);
        CompareMatrices(`CheckInverse(${aname},${bname})`, crot, IdentityMatrix, 2.0e-15);
    }

    function CheckCycle(cyclename, arot, brot, crot) {
        const xrot = Astronomy.CombineRotation(arot, brot);
        const irot = Astronomy.InverseRotation(xrot);
        CompareMatrices(cyclename, crot, irot, 2.0e-15);
    }

    function Test_RotRoundTrip() {
        const time = Astronomy.MakeTime(new Date('2067-05-30T14:45:00Z'));
        const observer = new Astronomy.Observer(+28, -82, 0);

        /*
            In each round trip, calculate a forward rotation and a backward rotation.
            Verify the two are inverse matrices.
        */

        /* Round trip #1: EQJ <==> EQD. */
        const eqj_eqd = Astronomy.Rotation_EQJ_EQD(time);
        const eqd_eqj = Astronomy.Rotation_EQD_EQJ(time);
        CheckInverse('eqj_eqd', 'eqd_eqj', eqj_eqd, eqd_eqj);

        /* Round trip #2: EQJ <==> ECL. */
        const eqj_ecl = Astronomy.Rotation_EQJ_ECL();
        const ecl_eqj = Astronomy.Rotation_ECL_EQJ();
        CheckInverse('eqj_ecl', 'ecl_eqj', eqj_ecl, ecl_eqj);

        /* Round trip #3: EQJ <==> HOR. */
        const eqj_hor = Astronomy.Rotation_EQJ_HOR(time, observer);
        const hor_eqj = Astronomy.Rotation_HOR_EQJ(time, observer);
        CheckInverse('eqj_hor', 'hor_eqj', eqj_hor, hor_eqj);

        /* Round trip #4: EQD <==> HOR. */
        const eqd_hor = Astronomy.Rotation_EQD_HOR(time, observer);
        const hor_eqd = Astronomy.Rotation_HOR_EQD(time, observer);
        CheckInverse('eqd_hor', 'hor_eqd', eqd_hor, hor_eqd);

        /* Round trip #5: EQD <==> ECL. */
        const eqd_ecl = Astronomy.Rotation_EQD_ECL(time);
        const ecl_eqd = Astronomy.Rotation_ECL_EQD(time);
        CheckInverse('eqd_ecl', 'ecl_eqd', eqd_ecl, ecl_eqd);

        /* Round trip #6: HOR <==> ECL. */
        const hor_ecl = Astronomy.Rotation_HOR_ECL(time, observer);
        const ecl_hor = Astronomy.Rotation_ECL_HOR(time, observer);
        CheckInverse('hor_ecl', 'ecl_hor', hor_ecl, ecl_hor);

        /*
            Verify that combining different sequences of rotations result
            in the expected combination.
            For example, (EQJ ==> HOR ==> ECL) must be the same matrix as (EQJ ==> ECL).
            Each of these is a "triangle" of relationships between 3 orientations.
            There are 4 possible ways to pick 3 orientations from the 4 to form a triangle.
            Because we have just proved that each transformation is reversible,
            we only need to verify the triangle in one cyclic direction.
        */
       CheckCycle('eqj_ecl, ecl_eqd, eqd_eqj', eqj_ecl, ecl_eqd, eqd_eqj);     /* excluded corner = HOR */
       CheckCycle('eqj_hor, hor_ecl, ecl_eqj', eqj_hor, hor_ecl, ecl_eqj);     /* excluded corner = EQD */
       CheckCycle('eqj_hor, hor_eqd, eqd_eqj', eqj_hor, hor_eqd, eqd_eqj);     /* excluded corner = ECL */
       CheckCycle('ecl_eqd, eqd_hor, hor_ecl', ecl_eqd, eqd_hor, hor_ecl);     /* excluded corner = EQJ */

       if (Verbose) console.log('JS Test_RotRoundTrip: PASS');
    }

    function Rotation_Pivot() {
        const tolerance = 1.0e-15;

        /* Test #1 */

        /* Start with an identity matrix. */
        const ident = Astronomy.IdentityMatrix();

        /* Pivot 90 degrees counterclockwise around the z-axis. */
        let r = Astronomy.Pivot(ident, 2, +90.0);

        /* Put the expected answer in 'a'. */
        const a = Astronomy.MakeRotation([
            [ 0, +1,  0],
            [-1,  0,  0],
            [ 0,  0, +1],
        ]);

        /* Compare actual 'r' with expected 'a'. */
        CompareMatrices('Rotation_Pivot #1', r, a, tolerance);

        /* Test #2. */

        /* Pivot again, -30 degrees around the x-axis. */
        r = Astronomy.Pivot(r, 0, -30.0);

        /* Pivot a third time, 180 degrees around the y-axis. */
        r = Astronomy.Pivot(r, 1, +180.0);

        /* Use the 'r' matrix to rotate a vector. */
        const v1 = new Astronomy.Vector(1, 2, 3, Astronomy.MakeTime(0));

        const v2 = Astronomy.RotateVector(r, v1);

        /* Initialize the expected vector 've'. */
        const ve = new Astronomy.Vector(+2.0, +2.3660254037844390, -2.0980762113533156, v1.t);

        CompareVectors('Rotation_Pivot #2', v2, ve, tolerance);

        if (Verbose) console.log('JS Rotation_Pivot: PASS');
    }


    Rotation_MatrixInverse();
    Rotation_MatrixMultiply();
    Rotation_Pivot();
    Test_EQJ_ECL();
    Test_GAL_EQJ_NOVAS('temp/galeqj.txt');

    Test_EQJ_EQD('Mercury');
    Test_EQJ_EQD('Venus');
    Test_EQJ_EQD('Mars');
    Test_EQJ_EQD('Jupiter');
    Test_EQJ_EQD('Saturn');

    Test_EQD_HOR('Mercury');
    Test_EQD_HOR('Venus');
    Test_EQD_HOR('Mars');
    Test_EQD_HOR('Jupiter');
    Test_EQD_HOR('Saturn');

    Test_RotRoundTrip();

    console.log('JS Rotation: PASS');
    return 0;
}


function Refraction() {
    for (let alt = -90.1; alt <= +90.1; alt += 0.001) {
        const refr = Astronomy.Refraction('normal', alt);
        const corrected = alt + refr;
        const inv_refr = Astronomy.InverseRefraction('normal', corrected);
        const check_alt = corrected + inv_refr;
        const diff = abs(check_alt - alt);
        if (diff > 2.0e-14)
            throw `JS Refraction: alt=${alt}, refr=${refr}, diff=${diff}`;
    }

    console.log('JS Refraction: PASS');
    return 0;
}


function MonthNumber(mtext) {
    const index = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].indexOf(mtext);
    if (index < 0)
        throw `Invalid month text "${mtext}"`;
    return 1 + index;
}


function Magnitude() {
    function LoadMagnitudeData(filename) {
        const lines = ReadLines(filename);
        let lnum = 0;
        let rows = [];
        for (let line of lines) {
            ++lnum;

            // [ Date__(UT)__HR:MN      APmag  S-brt               r        rdot            delta      deldot    S-T-O]
            // [ 2023-Mar-30 00:00      -4.01   1.17  0.719092953368  -0.1186373 1.20453495004726 -11.0204917  55.9004]
            let m = line.match(/^\s(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$/);
            if (m) {
                const year = int(m[1]);
                const month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].indexOf(m[2]);
                const day = int(m[3]);
                const hour = int(m[4]);
                const minute = int(m[5]);

                const item = {
                    lnum: lnum,
                    date: new Date(Date.UTC(year, month, day, hour, minute)),
                    mag: float(m[6]),
                    sbrt: float(m[7]),
                    helio_dist: float(m[8]),
                    helio_radvel: float(m[9]),
                    geo_dist: float(m[10]),
                    geo_radvel: float(m[11]),
                    phase_angle: float(m[12])
                };

                rows.push(item);
            }
        }

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
                diff_lo = min(diff_lo, diff);
                diff_hi = max(diff_hi, diff);
            }
        }
        let rms = sqrt(sum_squared_diff / data.rows.length);
        const limit = 0.012;
        const pass = (abs(diff_lo) < limit && abs(diff_hi) < limit);
        if (!pass || Verbose) console.log(`JS ${body.padEnd(8)} ${pass?"    ":"FAIL"}  diff_lo=${diff_lo.toFixed(4).padStart(8)}, diff_hi=${diff_hi.toFixed(4).padStart(8)}, rms=${rms.toFixed(4).padStart(8)}`);
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
            if (Verbose) console.log(`JS Saturn: date=${illum.time.date.toISOString()}  mag=${illum.mag.toFixed(8).padStart(12)}  ring_tilt=${illum.ring_tilt.toFixed(8).padStart(12)}`);
            const mag_diff = abs(illum.mag - item.mag);
            if (mag_diff > 1.0e-4) {
                console.log(`ERROR: Excessive magnitude error ${mag_diff}`);
                success = false;
            }
            const tilt_diff = abs(illum.ring_tilt - item.tilt);
            if (tilt_diff > 3.0e-5) {
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

        const lines = ReadLines(filename);
        let date = new Date(Date.UTC(2000, 0, 1));
        let max_diff = 0;
        for (let line of lines) {
            let token = line.split(/\s+/);
            let date1 = new Date(token[0]);
            let date2 = new Date(token[1]);
            let evt = Astronomy.SearchPeakMagnitude(body, date);
            if (evt.time.date < date1 || evt.time.date > date2)
                throw `Event time ${evt.time.toString()} is outside the range ${date1.toISOString()} .. ${date2.toISOString()}`;

            // How close are we to the center date?
            let date_center = new Date((date1.getTime() + date2.getTime())/2);
            let diff_hours = abs(evt.time.date - date_center) / (1000 * 3600);
            if (diff_hours > 7.1)
                throw `Excessive diff_hours = ${diff_hours} from center date ${date_center.toISOString()}`;

            max_diff = max(max_diff, diff_hours);
            date = date2;
        }
        if (Verbose) console.log(`JS TestMaxMag: ${lines.length} events, max error = ${max_diff.toFixed(3)} hours.`);
        return true;
    }

    let all_passed = true;
    const bodies = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Uranus', 'Neptune', 'Pluto'];
    for (let body of bodies) {
        const data = LoadMagnitudeData(`magnitude/${body}.txt`);
        if (!CheckMagnitudeData(body, data))
            all_passed = false;
    }

    if (!TestSaturn())
        all_passed = false;

    if (!TestMaxMag('magnitude/maxmag_Venus.txt', 'Venus'))
        all_passed = false;

    all_passed || Fail('Found excessive error in at least one test.');
    console.log('JS Magnitude: PASS');
    return 0;
}


function Constellation() {
    const filename = 'constellation/test_input.txt';
    const lines = ReadLines(filename);
    let lnum = 0;
    let failcount = 0;
    for (let row of lines) {
        ++lnum;
        let token = row.trim().split(/\s+/);
        if (token.length !== 4) {
            Fail(`Bad data in ${filename} line ${lnum}: found ${token.length} tokens.`);
        }
        const id = int(token[0]);
        const ra = float(token[1]);
        const dec = float(token[2]);
        const symbol = token[3];
        if (!/^[A-Z][A-Za-z]{2}$/.test(symbol)) {
            Fail(`Invalid symbol "${symbol}" in ${filename} line ${lnum}`);
        }
        const constel = Astronomy.Constellation(ra, dec);
        if (constel.symbol !== symbol) {
            ++failcount;
            console.error(`Star ${id}: expected ${symbol}, found ${constel.symbol} at B1875 RA=${constel.ra1875}, DEC=${constel.dec1875}`);
        }
    }
    if (failcount > 0) {
        Fail(`JS Constellation: ${failcount} failures.`);
    }
    console.log(`JS Constellation: PASS (verified ${lines.length})`);
    return 0;
}


function TransitFile(body, filename, limit_minutes, limit_sep) {
    const lines = ReadLines(filename);
    let lnum = 0;
    let max_minutes = 0;
    let max_sep = 0;
    let transit = Astronomy.SearchTransit(body, Astronomy.MakeTime(new Date(Date.UTC(1600, 0))));
    for (let row of lines) {
        ++lnum;
        let token = row.trim().split(/\s+/);

        /* 22:17 1881-11-08T00:57Z 03:38  3.8633 */
        if (token.length !== 4) {
            Fail(`JS TransitFile(${filename} line ${lnum}): bad data format.`);
        }
        const textp = token[1];
        const text1 = textp.substr(0, 11) + token[0] + 'Z';
        const text2 = textp.substr(0, 11) + token[2] + 'Z';
        let timep = Astronomy.MakeTime(new Date(textp));
        let time1 = Astronomy.MakeTime(new Date(text1));
        let time2 = Astronomy.MakeTime(new Date(text2));
        const separation = float(token[3]);

        // If the start time is after the peak time, it really starts on the previous day.
        if (time1.ut > timep.ut)
            time1 = time1.AddDays(-1.0);

        // If the finish time is before the peak time, it really starts on the following day.
        if (time2.ut < timep.ut)
            time2 = time2.AddDays(+1.0);

        const diff_start  = (24.0 * 60.0) * abs(time1.ut - transit.start.ut );
        const diff_peak   = (24.0 * 60.0) * abs(timep.ut - transit.peak.ut  );
        const diff_finish = (24.0 * 60.0) * abs(time2.ut - transit.finish.ut);
        const diff_sep = abs(separation - transit.separation);

        max_minutes = max(max_minutes, diff_start);
        max_minutes = max(max_minutes, diff_peak);
        max_minutes = max(max_minutes, diff_finish);
        if (max_minutes > limit_minutes) {
            Fail(`JS TransitFile(${filename} line ${lnum}): EXCESSIVE TIME ERROR = ${max_minutes} minutes.`);
        }

        max_sep = max(max_sep, diff_sep);
        if (max_sep > limit_sep) {
            Fail(`JS TransitFile(${filename} line ${lnum}): EXCESSIVE SEPARATION ERROR = ${max_sep} arcminutes.`);
        }

        transit = Astronomy.NextTransit(body, transit.finish);
    }

    console.log(`JS TransitFile(${filename}): PASS - verified ${lnum}, max minutes = ${max_minutes}, max sep arcmin = ${max_sep}`);
    return 0;
}


function Transit() {
    if (0 !== TransitFile('Mercury', 'eclipse/mercury.txt', 10.710, 0.2121)) {
        return 1;
    }

    if (0 !== TransitFile('Venus', 'eclipse/venus.txt', 9.109, 0.6772)) {
        return 1;
    }

    return 0;
}


function PlutoCheckDate(ut, arcmin_tolerance, x, y, z) {
    const time = Astronomy.MakeTime(ut);
    const timeText = time.toString();
    if (Verbose) console.log(`JS PlutoCheck: ${timeText} = ${time.ut} UT = ${time.tt} TT`);
    const vector = Astronomy.HelioVector('Pluto', time);
    const dx = v(vector.x - x);
    const dy = v(vector.y - y);
    const dz = v(vector.z - z);
    const diff = sqrt(dx*dx + dy*dy + dz*dz);
    const dist = (sqrt(x*x + y*y + z*z) - 1.0);       /* worst-case distance between Pluto and Earth */
    const arcmin = (diff / dist) * (180.0 * 60.0 / Math.PI);
    if (Verbose) console.log(`JS PlutoCheck: calc pos = [${vector.x}, ${vector.y}, ${vector.z}]`);
    if (Verbose) console.log(`JS PlutoCheck: ref  pos = [${x}, ${y}, ${z}]`);
    if (Verbose) console.log(`JS PlutoCheck: del  pos = [${vector.x - x}, ${vector.y - y}, ${vector.z - z}]`);
    if (Verbose) console.log(`JS PlutoCheck: diff = ${diff} AU, ${arcmin} arcmin`);
    if (Verbose) console.log("");
    if (v(arcmin) > arcmin_tolerance) {
        console.error("JS PlutoCheck: EXCESSIVE ERROR");
        return 1;
    }
    return 0;
}


function PlutoCheck() {
    if (0 != PlutoCheckDate(  +18250.0,  0.089, +37.4377303523676090, -10.2466292454075898, -14.4773101310875809)) return 1;
    if (0 != PlutoCheckDate( -856493.0,  4.067, +23.4292113199166252, +42.1452685817740829,  +6.0580908436642940)) return 1;
    if (0 != PlutoCheckDate( +435633.0,  0.016, -27.3178902095231813, +18.5887022581070305, +14.0493896259306936)) return 1;
    if (0 != PlutoCheckDate(       0.0,  8.e-9,  -9.8753673425269000, -27.9789270580402771,  -5.7537127596369588)) return 1;
    if (0 != PlutoCheckDate( +800916.0,  2.286, -29.5266052645301365, +12.0554287322176474, +12.6878484911631091)) return 1;
    console.log("JS PlutoCheck: PASS");
    return 0;
}


function GeoidTestCase(time, observer, ofdate) {
    let topo_moon = Astronomy.Equator(Astronomy.Body.Moon, time, observer, ofdate, false);
    let surface = Astronomy.ObserverVector(time, observer, ofdate);
    let geo_moon = Astronomy.GeoVector(Astronomy.Body.Moon, time, false);

    if (ofdate) {
        // GeoVector() returns J2000 coordinates. Convert to equator-of-date coordinates.
        const rot = Astronomy.Rotation_EQJ_EQD(time);
        geo_moon = Astronomy.RotateVector(rot, geo_moon);
    }

    const dx = Astronomy.KM_PER_AU * v((geo_moon.x - surface.x) - topo_moon.vec.x);
    const dy = Astronomy.KM_PER_AU * v((geo_moon.y - surface.y) - topo_moon.vec.y);
    const dz = Astronomy.KM_PER_AU * v((geo_moon.z - surface.z) - topo_moon.vec.z);
    const diff = sqrt(dx*dx + dy*dy + dz*dz);
    if (Verbose) console.log(`JS GeoidTestCase: ofdate=${ofdate}, time=${time.toISOString()}, lat=${observer.latitude}, lon=${observer.longitude}, ht=${observer.height}, surface=(${Astronomy.KM_PER_AU * surface.x}, ${Astronomy.KM_PER_AU * surface.y}, ${Astronomy.KM_PER_AU * surface.z}), diff = ${diff} km`);

    // Require 1 millimeter accuracy! (one millionth of a kilometer).
    if (diff > 1.0e-6) {
        console.error('JS GeoidTestCase: EXCESSIVE POSITION ERROR.');
        return 1;
    }

    // Verify that we can convert the surface vector back to an observer.
    const vobs = Astronomy.VectorObserver(surface, ofdate);
    const lat_diff = abs(vobs.latitude - observer.latitude);
    let lon_diff;

    /* Longitude is meaningless at the poles, so don't bother checking it there. */
    if (-89.99 <= observer.latitude && observer.latitude <= +89.99) {
        lon_diff = abs(vobs.longitude - observer.longitude);
        if (lon_diff > 180.0)
            lon_diff = 360.0 - lon_diff;
        lon_diff = abs(lon_diff * Math.cos(Astronomy.DEG2RAD * observer.latitude));
        if (lon_diff > 1.0e-6) {
            console.error(`JS GeoidTestCase: EXCESSIVE longitude check error = ${lon_diff}`);
            return 1;
        }
    } else {
        lon_diff = 0.0;
    }

    const h_diff = abs(vobs.height - observer.height);
    if (Verbose) console.log(`JS GeoidTestCase: vobs=(lat=${vobs.latitude}, lon=${vobs.longitude}, height=${vobs.height}), lat_diff=${lat_diff}, lon_diff=${lon_diff}, h_diff=${h_diff}`);

    if (lat_diff > 1.0e-6) {
        console.error(`JS GeoidTestCase: EXCESSIVE latitude check error = ${lat_diff}`);
        return 1;
    }

    if (h_diff > 0.001) {
        console.error(`JS GeoidTestCase: EXCESSIVE height check error = ${h_diff}`);
        return 1;
    }

    return 0;
}


function Geoid() {
    const time_list = [
        new Date('1066-09-27T18:00:00Z'),
        new Date('1970-12-13T15:42:00Z'),
        new Date('1970-12-13T15:43:00Z'),
        new Date('2015-03-05T02:15:45Z')
    ];

    const observer_list = [
        new Astronomy.Observer( +1.5,   +2.7,    7.4),
        new Astronomy.Observer(-53.7, +141.7, +100.0),
        new Astronomy.Observer(+30.0,  -85.2,  -50.0),
        new Astronomy.Observer(+90.0,  +45.0,  -50.0),
        new Astronomy.Observer(-90.0, -180.0,    0.0)
    ];

    // Test a variety of times and locations, in both supported orientation systems.

    let observer, time;
    for (observer of observer_list) {
        for (time of time_list) {
            if (0 != GeoidTestCase(time, observer, false))
                return 1;
            if (0 != GeoidTestCase(time, observer, true))
                return 1;
        }
    }

    // More exhaustive tests for a single time value across many different geographic coordinates.

    time = new Date('2021-06-20T15:08:00Z');
    for (let lat = -90; lat <= +90; lat += 1) {
        for (let lon = -175; lon <= +180; lon += 5) {
            observer = new Astronomy.Observer(lat, lon, 0.0);
            if (0 != GeoidTestCase(time, observer, true))
                return 1;
        }
    }

    console.log('JS GeoidTest: PASS');
    return 0;
}


function JupiterMoons_CheckJpl(mindex, tt, pos, vel) {
    const pos_tolerance = 9.0e-4;
    const vel_tolerance = 9.0e-4;
    const time = Astronomy.AstroTime.FromTerrestrialTime(tt);
    const jm = Astronomy.JupiterMoons(time);
    const moon = SelectJupiterMoon(jm, mindex);

    let dx = v(pos[0] - moon.x);
    let dy = v(pos[1] - moon.y);
    let dz = v(pos[2] - moon.z);
    let mag = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    const pos_diff = sqrt(dx*dx + dy*dy + dz*dz) / mag;
    if (pos_diff > pos_tolerance) {
        console.error(`JS JupiterMoons_CheckJpl(mindex=${mindex}, tt=${tt}): excessive position error ${pos_diff}`);
        return 1;
    }

    dx = v(vel[0] - moon.vx);
    dy = v(vel[1] - moon.vy);
    dz = v(vel[2] - moon.vz);
    mag = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    const vel_diff = sqrt(dx*dx + dy*dy + dz*dz) / mag;
    if (vel_diff > vel_tolerance) {
        console.error(`JS JupiterMoons_CheckJpl(mindex=${mindex}, tt=${tt}): excessive velocity error ${vel_diff}`);
        return 1;
    }

    if (Verbose) console.debug(`JS JupiterMoons_CheckJpl: mindex=${mindex}, tt=${tt}, pos_diff=${pos_diff}, vel_diff=${vel_diff}`);
    return 0;
}


function JupiterMoons() {
    for (let mindex = 0; mindex < 4; ++mindex) {
        const filename = `jupiter_moons/horizons/jm${mindex}.txt`;
        const lines = ReadLines(filename);
        let lnum = 0;
        let found = false;
        let part = -1;
        const expected_count = 5001;
        let count = 0;
        let tt, match, pos, vel;
        for (let line of lines) {
            ++lnum;
            if (!found) {
                if (line == '$$SOE') {
                    found = true;
                    part = 0;
                } else if (line.startsWith('Revised:')) {
                    if (line.length !== 79) {
                        console.error(`JS JupiterMoons(${filename} line ${lnum}): unexpected line length.`);
                        return 1;
                    }
                    console.log(line.substr(76));
                    const check_mindex = int(line.substr(76)) - 501;
                    if (mindex !== check_mindex) {
                        console.error(`JS JupiterMoons(${filename} line ${lnum}): moon index does not match: check=${check_mindex}, mindex=${mindex}.`);
                        return 1;
                    }
                }
            } else if (line == '$$EOE') {
                break;
            } else {
                switch (part) {
                case 0:
                    // 2446545.000000000 = A.D. 1986-Apr-24 12:00:00.0000 TDB
                    tt = float(line.split()[0]) - 2451545.0;    // convert JD to J2000 TT
                    break;

                case 1:
                    // X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                    match = /\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)/.exec(line);
                    if (!match) {
                        console.error(`JS JupiterMoons(${filename} line ${lnum}): cannot parse position vector.`);
                        return 1;
                    }
                    pos = [ float(match[1]), float(match[2]), float(match[3]) ];
                    break;

                case 2:
                    // VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                    match = /\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)/.exec(line);
                    if (!match) {
                        console.error(`JS JupiterMoons(${filename} line ${lnum}): cannot parse velocity vector.`);
                        return 1;
                    }
                    vel = [ float(match[1]), float(match[2]), float(match[3]) ];
                    if (JupiterMoons_CheckJpl(mindex, tt, pos, vel)) {
                        console.error(`JS JupiterMoons(${filename} line ${lnum}): FAILED VERIFICATION.`);
                        return 1;
                    }
                    ++count;
                    break;

                default:
                    console.error(`JS JupiterMoons(${filename} line ${lnum}): unexpected part = ${part}`);
                    return 1;
                }
                part = (part + 1) % 3;
            }
        }
        if (count !== expected_count) {
            console.error(`JS JupiterMoons: Expected ${expected_count} test cases, but found ${count}`);
            return 1;
        }
    }
    console.log(`JS JupiterMoons: PASS`);
    return 0;
}


function Issue103() {
    // https://github.com/cosinekitty/astronomy/issues/103

    const observer = new Astronomy.Observer(29, -81, 10);
    const ut = -51279.9420508868643083;
    const time = Astronomy.MakeTime(ut);
    const ofdate = Astronomy.Equator(Astronomy.Body.Moon, time, observer, true, true);

    console.log(`tt  = ${time.tt.toFixed(16)}`);
    console.log(`ra  = ${ofdate.ra.toFixed(16)}`);
    console.log(`dec = ${ofdate.dec.toFixed(16)}`);
    const hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, false);
    console.log(`az  = ${hor.azimuth.toFixed(16)}`);
    console.log(`alt = ${hor.altitude.toFixed(16)}`);

    return 0;
}


function AberrationTest() {
    const THRESHOLD_SECONDS = 0.4;
    const filename = 'equatorial/Mars_j2000_ofdate_aberration.txt';
    const lines = ReadLines(filename);
    let lnum = 0;
    let found_begin = false;
    let max_diff_seconds = 0;
    let count = 0;
    for (let line of lines) {
        ++lnum;
        if (!found_begin) {
            if (line === '$$SOE')
                found_begin = true;
        } else if (line === '$$EOE') {
            break;
        } else {
            // 2459371.500000000 *   118.566080210  22.210647456 118.874086738  22.155784122
            const ut = float(line.trim().split(/\s+/)[0]) - 2451545.0;    // convert JD to J2000 UT day value
            const time = Astronomy.MakeTime(ut);
            const tokens = line.substr(22).trim().split(/\s+/);
            const jra = float(tokens[0]);
            const jdec = float(tokens[1]);
            const dra = float(tokens[2]);
            const ddec = float(tokens[3]);

            // Create spherical coordinates with magnitude = speed of light.
            // This helps us perform the aberration correction later.
            const eqj_sphere = new Astronomy.Spherical(jdec, jra, Astronomy.C_AUDAY);

            // Convert EQJ angular coordinates (jra, jdec) to an EQJ vector.
            // This vector represents a ray of light, only it is travelling from the observer toward the star.
            // The backwards direction makes the math simpler.
            const eqj_vec = Astronomy.VectorFromSphere(eqj_sphere, time);

            // Calculate the Earth's barycentric velocity vector in EQJ coordinates.
            const eqj_earth = Astronomy.BaryState(Astronomy.Body.Earth, time);

            // Use non-relativistic approximation: add light vector to Earth velocity vector.
            // This gives aberration-corrected apparent position of the star in EQJ.
            eqj_vec.x += eqj_earth.vx;
            eqj_vec.y += eqj_earth.vy;
            eqj_vec.z += eqj_earth.vz;

            // Calculate the rotation matrix that converts J2000 coordinates to of-date coordinates.
            const rot = Astronomy.Rotation_EQJ_EQD(time);

            // Use the rotation matrix to re-orient the EQJ vector to a EQD vector.
            const eqd_vec = Astronomy.RotateVector(rot, eqj_vec);

            // Convert the EQD vector back to spherical angular coordinates.
            const eqd_sphere = Astronomy.SphereFromVector(eqd_vec);

            // Calculate the differences in RA and DEC between expected and calculated values.
            const factor = cos(eqd_sphere.lat * Astronomy.DEG2RAD);     // RA errors are less important toward the poles.
            const xra = factor * abs(eqd_sphere.lon - dra);
            const xdec = abs(eqd_sphere.lat - ddec);
            const diff_seconds = 3600 * sqrt(xra*xra + xdec*xdec);
            if (Verbose) console.debug(`JS AberrationTest(${filename} line ${lnum}): xra=${xra.toFixed(6)} deg, xdec=${xdec.toFixed(6)} deg, diff_seconds=${diff_seconds.toFixed(3)}`);
            if (diff_seconds > THRESHOLD_SECONDS) {
                console.error(`JS AberrationTest(${filename} line ${lnum}): EXCESSIVE ANGULAR ERROR = ${diff_seconds.toFixed(3)} seconds.`);
                return 1;
            }
            if (diff_seconds > max_diff_seconds)
                max_diff_seconds = diff_seconds;
            ++count;
        }
    }
    console.log(`JS AberrationTest(${filename}): PASS - Tested ${count} cases. max_diff_seconds = ${max_diff_seconds.toFixed(3)}`);
    return 0;
}


function StateVectorDiff(relative, vec, x, y, z) {
    const dx = v(vec[0] - x);
    const dy = v(vec[1] - y);
    const dz = v(vec[2] - z);
    let diff_squared = dx*dx + dy*dy + dz*dz;
    if (relative)
        diff_squared /= (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    return sqrt(diff_squared);
}


function VerifyState(func, score, filename, lnum, time, pos, vel, r_thresh, v_thresh) {
    const state = func.Eval(time);

    const rdiff = StateVectorDiff((r_thresh > 0.0), pos, state.x, state.y, state.z);
    if (rdiff > score.max_rdiff)
        score.max_rdiff = rdiff;

    const vdiff = StateVectorDiff((v_thresh > 0.0), vel, state.vx, state.vy, state.vz);
    if (vdiff > score.max_vdiff)
        score.max_vdiff = vdiff;

    if (rdiff > Math.abs(r_thresh)) {
        console.error(`JS VerifyState(${filename} line ${lnum}): EXCESSIVE POSITION ERROR = ${rdiff.toExponential(3)}`);
        return 1;
    }

    if (vdiff > Math.abs(v_thresh)) {
        console.error(`JS VerifyState(${filename} line ${lnum}): EXCESSIVE VELOCITY ERROR = ${vdiff.toExponential(3)}`);
        return 1;
    }

    return 0;
}


function VerifyStateBody(func, filename, r_thresh, v_thresh) {
    const lines = ReadLines(filename);
    let lnum = 0;
    let found = false;
    let part = -1;
    let count = 0;
    let tt, match, pos, vel, time;
    let score = { max_rdiff:0, max_vdiff:0 };
    for (let line of lines) {
        ++lnum;
        if (!found) {
            if (line == '$$SOE') {
                found = true;
                part = 0;
            }
        } else if (line == '$$EOE') {
            break;
        } else {
            switch (part) {
            case 0:
                // 2446545.000000000 = A.D. 1986-Apr-24 12:00:00.0000 TDB
                tt = float(line.split()[0]) - 2451545.0;    // convert JD to J2000 TT
                time = Astronomy.AstroTime.FromTerrestrialTime(tt);
                break;

            case 1:
                // X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                match = /\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)/.exec(line);
                if (!match) {
                    console.error(`JS VerifyStateBody(${filename} line ${lnum}): cannot parse position vector.`);
                    return 1;
                }
                pos = [ float(match[1]), float(match[2]), float(match[3]) ];
                break;

            case 2:
                // VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                match = /\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)/.exec(line);
                if (!match) {
                    console.error(`JS VerifyStateBody(${filename} line ${lnum}): cannot parse velocity vector.`);
                    return 1;
                }
                vel = [ float(match[1]), float(match[2]), float(match[3]) ];
                if (VerifyState(func, score, filename, lnum, time, pos, vel, r_thresh, v_thresh))
                    return 1;
                ++count;
                break;

            default:
                console.error(`JS VerifyStateBody(${filename} line ${lnum}): unexpected part = ${part}`);
                return 1;
            }
            part = (part + 1) % 3;
        }
    }

    if (Verbose) console.debug(`JS VerifyStateBody(${filename}): PASS - Tested ${count} cases. max rdiff=${score.max_rdiff.toExponential(3)}, vdiff=${score.max_vdiff.toExponential(3)}`);
    return 0;
}


// Constants for use inside unit tests only; they doesn't make sense for public consumption.
const Body_GeoMoon = -100;
const Body_Geo_EMB = -101;

class BaryStateFunc {
    constructor(body) {
        this.body = body;
    }

    Eval(time) {
        if (this.body === Body_GeoMoon)
            return Astronomy.GeoMoonState(time);

        if (this.body === Body_Geo_EMB)
            return Astronomy.GeoEmbState(time);

        return Astronomy.BaryState(this.body, time);
    }
}


function BaryStateTest() {
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Sun),     'barystate/Sun.txt',     -1.224e-05, -1.134e-07)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Mercury), 'barystate/Mercury.txt',  1.672e-04,  2.698e-04)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Venus),   'barystate/Venus.txt',    4.123e-05,  4.308e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Earth),   'barystate/Earth.txt',    2.296e-05,  6.359e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Mars),    'barystate/Mars.txt',     3.107e-05,  5.550e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Jupiter), 'barystate/Jupiter.txt',  7.389e-05,  2.471e-04)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Saturn),  'barystate/Saturn.txt',   1.067e-04,  3.220e-04)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Uranus),  'barystate/Uranus.txt',   9.035e-05,  2.519e-04)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Neptune), 'barystate/Neptune.txt',  9.838e-05,  4.446e-04)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Pluto),   'barystate/Pluto.txt',    4.259e-05,  7.827e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.Moon),    "barystate/Moon.txt",     2.354e-05,  6.604e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Astronomy.Body.EMB),     "barystate/EMB.txt",      2.353e-05,  6.511e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Body_GeoMoon),           "barystate/GeoMoon.txt",  4.086e-05,  5.347e-05)) return 1;
    if (VerifyStateBody(new BaryStateFunc(Body_Geo_EMB),           "barystate/GeoEMB.txt",   4.076e-05,  5.335e-05)) return 1;
    console.log('JS BaryStateTest: PASS');
    return 0;
}


class HelioStateFunc {
    constructor(body) {
        this.body = body;
    }

    Eval(time) {
        return Astronomy.HelioState(this.body, time);
    }
}


function HelioStateTest() {
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.SSB),     'heliostate/SSB.txt',     -1.209e-05, -1.125e-07)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Mercury), 'heliostate/Mercury.txt',  1.481e-04,  2.756e-04)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Venus),   'heliostate/Venus.txt',    3.528e-05,  4.485e-05)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Earth),   'heliostate/Earth.txt',    1.476e-05,  6.105e-05)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Mars),    'heliostate/Mars.txt',     3.154e-05,  5.603e-05)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Jupiter), 'heliostate/Jupiter.txt',  7.455e-05,  2.562e-04)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Saturn),  'heliostate/Saturn.txt',   1.066e-04,  3.150e-04)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Uranus),  'heliostate/Uranus.txt',   9.034e-05,  2.712e-04)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Neptune), 'heliostate/Neptune.txt',  9.834e-05,  4.534e-04)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Pluto),   'heliostate/Pluto.txt',    4.271e-05,  1.198e-04)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.Moon),    'heliostate/Moon.txt',     1.477e-05,  6.195e-05)) return 1;
    if (VerifyStateBody(new HelioStateFunc(Astronomy.Body.EMB),     'heliostate/EMB.txt',      1.476e-05,  6.106e-05)) return 1;
    console.log('JS HelioStateTest: PASS');
    return 0;
}


class TopoStateFunc {
    constructor(body) {
        this.body = body;
    }

    Eval(time) {
        const observer = new Astronomy.Observer(30.0, -80.0, 1000.0);

        let observer_state = Astronomy.ObserverState(time, observer, false);
        let state;
        if (this.body == Body_Geo_EMB) {
            state = Astronomy.GeoEmbState(time);
        } else if (this.body == Astronomy.Body.Earth) {
            state = new Astronomy.StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time);
        } else {
            throw `JS TopoStateFunction: unsupported body ${this.body}`;
        }

        state.x  -= observer_state.x;
        state.y  -= observer_state.y;
        state.z  -= observer_state.z;
        state.vx -= observer_state.vx;
        state.vy -= observer_state.vy;
        state.vz -= observer_state.vz;

        return state;
    }
}

function TopoStateTest() {
    if (VerifyStateBody(new TopoStateFunc(Astronomy.Body.Earth),  "topostate/Earth_N30_W80_1000m.txt",  2.108e-04, 2.430e-04)) return 1;
    if (VerifyStateBody(new TopoStateFunc(Body_Geo_EMB),          "topostate/EMB_N30_W80_1000m.txt",    7.195e-04, 2.497e-04)) return 1;
    console.log("JS TopoStateTest: PASS");
    return 0;
}



function TwilightTest() {
    const tolerance_seconds = 60.0;
    const filename = 'riseset/twilight.txt';
    const lines = ReadLines(filename);
    let lnum = 0;
    let max_diff = 0.0;

    const name = [
        "astronomical dawn",
        "nautical dawn",
        "civil dawn",
        "civil dusk",
        "nautical dusk",
        "astronomical dusk"
    ];

    for (let line of lines) {
        ++lnum;
        const tokens = line.split(/\s+/);
        if (tokens.length !== 9) {
            console.error(`JS TwilightTest: FAIL(${filename} line ${lnum}): incorrect number of tokens = ${tokens.length}.`);
            return 1;
        }
        const lat = float(tokens[0]);
        const lon = float(tokens[1]);
        const observer = new Astronomy.Observer(lat, lon, 0);

        const searchDate = Astronomy.MakeTime(new Date(tokens[2]));
        const correctTimes = [
            Astronomy.MakeTime(new Date(tokens[3])),      // astronomical dawn
            Astronomy.MakeTime(new Date(tokens[4])),      // nautical dawn
            Astronomy.MakeTime(new Date(tokens[5])),      // civil dawn
            Astronomy.MakeTime(new Date(tokens[6])),      // civil dusk
            Astronomy.MakeTime(new Date(tokens[7])),      // nautical dusk
            Astronomy.MakeTime(new Date(tokens[8]))       // astronomical dusk
        ];

        const calcTimes = [
            Astronomy.SearchAltitude(Astronomy.Body.Sun, observer, +1, searchDate, 1, -18),   // astronomical dawn
            Astronomy.SearchAltitude(Astronomy.Body.Sun, observer, +1, searchDate, 1, -12),   // nautical dawn
            Astronomy.SearchAltitude(Astronomy.Body.Sun, observer, +1, searchDate, 1, -6),    // civil dawn
            Astronomy.SearchAltitude(Astronomy.Body.Sun, observer, -1, searchDate, 1, -6),    // civil dusk
            Astronomy.SearchAltitude(Astronomy.Body.Sun, observer, -1, searchDate, 1, -12),   // nautical dusk
            Astronomy.SearchAltitude(Astronomy.Body.Sun, observer, -1, searchDate, 1, -18),   // astronomical dusk
        ];

        for (let i = 0; i < correctTimes.length; ++i) {
            const correct = correctTimes[i];
            const calc = calcTimes[i];
            const diff = 86400 * abs(calc.ut - correct.ut);
            if (diff > tolerance_seconds) {
                console.error(`JS TwilightTest(${filename} line ${lnum}): EXCESSIVE ERROR = ${diff} seconds for ${name[i]}`);
                console.error(`Expected ${correct} but calculated ${calc}`);
                return 1;
            }
            if (diff > max_diff)
                max_diff = diff;
        }
    }
    console.log(`JS TwilightTest: PASS (${lnum} test cases, max error = ${max_diff.toFixed(3)} seconds)`);
    return 0;
}


function Libration(filename) {
    const lines = ReadLines(filename);
    let max_diff_elon = 0.0;
    let max_diff_elat = 0.0;
    let max_diff_distance = 0.0;
    let max_diff_diam = 0.0;
    let max_eclip_lon = -900.0;
    let count = 0;
    let lnum = 0;
    for (let line of lines) {
        ++lnum;
        if (lnum === 1) {
            if (line !== "   Date       Time    Phase    Age    Diam    Dist     RA        Dec      Slon      Slat     Elon     Elat   AxisA") {
                console.error(`JS Libration(${filename} line ${lnum}): unexpected header line.`);
                return 1;
            }
        } else {
            const token = line.split(/\s+/);
            if (token.length !== 16) {
                console.error(`JS Libration: FAIL(${filename} line ${lnum}): incorrect number of tokens = ${token.length}.`);
                return 1;
            }

            const day = int(token[0]);
            const month = MonthNumber(token[1]);
            const year = int(token[2]);
            const hmtoken = token[3].split(':');
            if (hmtoken.length !== 2) {
                Console.WriteLine(`JS Libration(${filename} line ${lnum}): expected hh:mm but found '${token[3]}'`);
                return 1;
            }
            const hour = int(hmtoken[0]);
            const minute = int(hmtoken[1]);
            const time = Astronomy.MakeTime(new Date(Date.UTC(year, month-1, day, hour, minute)));

            const diam = float(token[7]) / 3600.0;
            const dist = float(token[8]);
            const elon = float(token[13]);
            const elat = float(token[14]);

            const lib = Astronomy.Libration(time);

            const diff_elon = 60.0 * abs(lib.elon - elon);
            if (diff_elon > max_diff_elon)
                max_diff_elon = diff_elon;

            const diff_elat = 60.0 * abs(lib.elat - elat);
            if (diff_elat > max_diff_elat)
                max_diff_elat = diff_elat;

            const diff_distance = abs(lib.dist_km - dist);
            if (diff_distance > max_diff_distance)
                max_diff_distance = diff_distance;

            const diff_diam = abs(lib.diam_deg - diam);
            if (diff_diam > max_diff_diam)
                max_diff_diam = diff_diam;

            if (lib.mlon > max_eclip_lon)
                max_eclip_lon = lib.mlon;

            if (diff_elon > 0.1304) {
                console.error(`JS Libration(${filename} line ${lnum}): EXCESSIVE diff_elon = ${diff_elon} arcmin`);
                return 1;
            }

            if (diff_elat > 1.6476) {
                console.error(`JS Libration(${filename} line ${lnum}): EXCESSIVE diff_elat = ${diff_elat} arcmin`);
                return 1;
            }

            if (diff_distance > 54.377) {
                console.error(`JS Libration(${filename} line ${lnum}): EXCESSIVE diff_distance = ${diff_distance} km`);
                return 1;
            }

            if (diff_diam > 0.00009) {
                console.error(`JS Libration(${filename}): EXCESSIVE diff_diam = ${diff_diam} degrees.`);
                return 1;
            }
            ++count;
        }
    }
    if (max_eclip_lon < 359.0 || max_eclip_lon > 360.0) {
        console.error(`JS Libration(${filename}): INVALID max ecliptic longitude = ${max_eclip_lon.toFixed(3)} degrees.`);
        return 1;
    }
    console.log(`JS Libration(${filename}): PASS (${count} test cases, max_diff_elon = ${max_diff_elon} arcmin, max_diff_elat = ${max_diff_elat} arcmin, max_diff_distance = ${max_diff_distance} km, max_diff_diam = ${max_diff_diam} deg)`);
    return 0;
}


function LibrationTest() {
    return (
        Libration("libration/mooninfo_2020.txt") ||
        Libration("libration/mooninfo_2021.txt") ||
        Libration("libration/mooninfo_2022.txt")
    );
}


function AxisTest() {
    if (0 !== AxisTestBody(Astronomy.Body.Sun,      "axis/Sun.txt",       0.0))        return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Mercury,  "axis/Mercury.txt",   0.074340))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Venus,    "axis/Venus.txt",     0.0))        return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Earth,    "axis/Earth.txt",     0.000591))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Moon,     "axis/Moon.txt",      0.264845))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Mars,     "axis/Mars.txt",      0.075323))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Jupiter,  "axis/Jupiter.txt",   0.000324))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Saturn,   "axis/Saturn.txt",    0.000304))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Uranus,   "axis/Uranus.txt",    0.0))        return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Neptune,  "axis/Neptune.txt",   0.000464))   return 1;
    if (0 !== AxisTestBody(Astronomy.Body.Pluto,    "axis/Pluto.txt",     0.0))        return 1;
    console.log("JS AxisTest: PASS");
    return 0;
}


function AxisTestBody(body, filename, arcmin_tolerance) {
    const lines = ReadLines(filename);
    let max_arcmin = 0;
    let lnum = 0;
    let count = 0;
    let found_data = false;
    for (let line of lines) {
        ++lnum;
        if (!found_data) {
            if (line === '$$SOE')
                found_data = true;
        } else {
            if (line === '$$EOE')
                break;

            const token = line.trim().split(/\s+/);
            // [ '1970-Jan-01', '00:00', '2440587.500000000', '281.01954', '61.41577' ]
            if (token.length !== 5) {
                console.error(`JS AxisBodyTest(${filename} line ${lnum}): expected 5 tokens, found ${tokens.length}`);
                return 1;
            }
            const jd = float(token[2]);
            const ra = float(token[3]);
            const dec = float(token[4]);
            const time = Astronomy.MakeTime(jd - 2451545.0);
            const axis = Astronomy.RotationAxis(body, time);
            const sphere = new Astronomy.Spherical(dec, ra, 1);
            const north = Astronomy.VectorFromSphere(sphere, time);
            const arcmin = 60 * Astronomy.AngleBetween(north, axis.north);
            if (arcmin > max_arcmin)
                max_arcmin = arcmin;

            ++count;
        }
    }
    if (Verbose) console.debug(`JS AxisTestBody(${body}): ${count} test cases, max arcmin error = ${max_arcmin}.`);
    if (max_arcmin > arcmin_tolerance) {
        console.log(`JS AxisTestBody(${body}): EXCESSIVE ERROR = ${max_arcmin} arcmin.`);
        return 1;
    }
    return 0;
}


function MoonNodes() {
    const filename = 'moon_nodes/moon_nodes.txt';
    const lines = ReadLines(filename);
    let lnum = 0;
    let prev_kind = '?';
    let max_angle = 0;
    let max_minutes = 0;
    let node;
    for (let line of lines) {
        ++lnum;
        const token = line.trim().split(/\s+/);
        if (token.length !== 4) {
            console.error(`JS MoonNodes(${filename} line ${lnum}): Unexpected number of tokens = ${tokens.length}`);
            return 1;
        }
        const kind = token[0];
        if (kind !== 'A' && kind !== 'D') {
            console.error(`JS MoonNodes(${filename} line ${lnum}): Invalid node kind '${kind}'`);
            return 1;
        }
        if (kind === prev_kind) {
            console.error(`JS MoonNodes(${filename} line ${lnum}): Duplicate node kind.`);
            return 1;
        }

        // Convert file's EQD Moon position angles to vector form.
        const time = Astronomy.MakeTime(new Date(token[1]));
        const ra = float(token[2]);
        const dec = float(token[3]);
        const sphere = new Astronomy.Spherical(dec, 15*ra, 1);
        const vec_test = Astronomy.VectorFromSphere(sphere, time);

        // Calculate EQD coordinates of the Moon. Verify against input file.
        const vec_eqj = Astronomy.GeoMoon(time);
        const rot = Astronomy.Rotation_EQJ_EQD(time);
        const vec_eqd = Astronomy.RotateVector(rot, vec_eqj);
        const angle = Astronomy.AngleBetween(vec_test, vec_eqd);
        const diff_angle = 60 * abs(angle);
        if (diff_angle > max_angle)
            max_angle = diff_angle;
        if (diff_angle > 1.54) {
            console.error(`JS MoonNodes(${filename} line ${lnum}): EXCESSIVE equatorial error = ${diff_angle.toFixed(3)} arcmin.`);
            return 1;
        }

        // Test the Astronomy Engine moon node searcher.
        if (lnum === 1) {
            // The very first time, so search for the first node in the series.
            // Back up a few days to make sure we really are finding it ourselves.
            const earlier = time.AddDays(-6.5472);      // back up by a weird amount of time
            node = Astronomy.SearchMoonNode(earlier);
        } else {
            node = Astronomy.NextMoonNode(node);
        }

        // Verify the ecliptic latitude is very close to zero at the alleged node.
        const ecl = Astronomy.EclipticGeoMoon(node.time);
        const diff_lat = 60 * abs(ecl.lat);
        if (diff_lat > 8.1e-4) {
            console.error(`JS MoonNodes(${filename} line ${lnum}): found node has excessive latitude = ${diff_lat} arcmin.`);
            return 1;
        }

        // Verify the time agrees with Espenak's time to within a few minutes.
        const diff_minutes = (24 * 60) * abs(node.time.tt - time.tt);
        if (diff_minutes > max_minutes)
            max_minutes = diff_minutes;

        // Verify the kind of node matches what Espenak says (ascending or descending).
        if (kind === 'A' && node.kind !== Astronomy.NodeEventKind.Ascending) {
            console.error(`JS MoonNodes(${filename} line ${lnum}): did not find ascending node as expected.`);
            return 1;
        }
        if (kind === 'D' && node.kind !== Astronomy.NodeEventKind.Descending) {
            console.error(`JS MoonNodes(${filename} line ${lnum}): did not find descending node as expected.`);
            return 1;
        }

        // Prepare for the next iteration.
        prev_kind = kind;
    }
    if (max_minutes > 3.681) {
        console.error(`JS MoonNodes: EXCESSIVE time prediction error = ${max_minutes.toFixed(3)} minutes.`);
        return 1;
    }
    console.log(`JS MoonNodes: PASS (${lnum} nodes, max equ error = ${max_angle.toFixed(3)} arcmin, max time error = ${max_minutes.toFixed(3)} minutes.)`);
    return 0;
}


class LagrangeFunc {
    constructor(point, major_body, minor_body) {
        this.point = point;
        this.major_body = major_body;
        this.minor_body = minor_body;
    }

    Eval(time) {
        return Astronomy.LagrangePoint(this.point, time, this.major_body, this.minor_body);
    }
}


function VerifyStateLagrange(major_body, minor_body, point, filename, r_thresh, v_thresh) {
    const func = new LagrangeFunc(point, major_body, minor_body);
    return VerifyStateBody(func, filename, r_thresh, v_thresh);
}


function LagrangeTest() {
    // Test Sun/EMB Lagrange points.
    if (0 != VerifyStateLagrange(Astronomy.Body.Sun, Astronomy.Body.EMB, 1, "lagrange/semb_L1.txt",   1.33e-5, 6.13e-5)) return 1;
    if (0 != VerifyStateLagrange(Astronomy.Body.Sun, Astronomy.Body.EMB, 2, "lagrange/semb_L2.txt",   1.33e-5, 6.13e-5)) return 1;
    if (0 != VerifyStateLagrange(Astronomy.Body.Sun, Astronomy.Body.EMB, 4, "lagrange/semb_L4.txt",   3.75e-5, 5.28e-5)) return 1;
    if (0 != VerifyStateLagrange(Astronomy.Body.Sun, Astronomy.Body.EMB, 5, "lagrange/semb_L5.txt",   3.75e-5, 5.28e-5)) return 1;

    // Test Earth/Moon Lagrange points.
    if (0 != VerifyStateLagrange(Astronomy.Body.Earth, Astronomy.Body.Moon, 1, "lagrange/em_L1.txt",  3.79e-5, 5.06e-5)) return 1;
    if (0 != VerifyStateLagrange(Astronomy.Body.Earth, Astronomy.Body.Moon, 2, "lagrange/em_L2.txt",  3.79e-5, 5.06e-5)) return 1;
    if (0 != VerifyStateLagrange(Astronomy.Body.Earth, Astronomy.Body.Moon, 4, "lagrange/em_L4.txt",  3.79e-5, 1.59e-3)) return 1;
    if (0 != VerifyStateLagrange(Astronomy.Body.Earth, Astronomy.Body.Moon, 5, "lagrange/em_L5.txt",  3.79e-5, 1.59e-3)) return 1;

    console.log("JS LagrangeTest: PASS");
    return 0;
}


function SiderealTimeTest() {
    const date = new Date('2022-03-15T21:50:00Z');
    const gast = Astronomy.SiderealTime(date);
    const correct = 9.398368460418821;
    const diff = abs(gast - correct);
    console.log(`JS SiderealTimeTest: gast=${gast.toFixed(15)}, correct=${correct.toFixed(15)}, diff=${diff.toExponential(3)}`);
    if (diff > 1.0e-15) {
        console.error('JS SiderealTimeTest: FAIL - excessive error.');
        return 1;
    }
    console.log('JS SiderealTimeTest: PASS');
    return 0;
}


const UnitTests = {
    aberration:             AberrationTest,
    axis:                   AxisTest,
    barystate:              BaryStateTest,
    constellation:          Constellation,
    elongation:             Elongation,
    geoid:                  Geoid,
    global_solar_eclipse:   GlobalSolarEclipse,
    heliostate:             HelioStateTest,
    issue_103:              Issue103,
    jupiter_moons:          JupiterMoons,
    lagrange:               LagrangeTest,
    libration:              LibrationTest,
    local_solar_eclipse:    LocalSolarEclipse,
    lunar_apsis:            LunarApsis,
    lunar_eclipse:          LunarEclipse,
    lunar_eclipse_78:       LunarEclipseIssue78,
    magnitude:              Magnitude,
    moon_nodes:             MoonNodes,
    moon_phase:             MoonPhase,
    planet_apsis:           PlanetApsis,
    pluto:                  PlutoCheck,
    refraction:             Refraction,
    rise_set:               RiseSet,
    rotation:               Rotation,
    seasons:                Seasons,
    seasons187:             SeasonsIssue187,
    sidereal:               SiderealTimeTest,
    topostate:              TopoStateTest,
    transit:                Transit,
    twilight:               TwilightTest,
};


function TestAll() {
    for (let name in UnitTests) {
        if (UnitTests[name]()) {
            return 1;
        }
    }
    console.log('JS ALL PASS');
}


function main() {
    let args = process.argv.slice(2);
    if (args.length > 0 && args[0] === '-v') {
        Verbose = true;
        args = args.slice(1);
    }

    if (args.length === 1) {
        const name = args[0];
        if (name === 'all') {
            return TestAll();
        }
        if (name === 'astro_check') {
            // This is a special case because the output is redirected and parsed.
            // It needs to be invoked separately, and cannot emit any extraneous output.
            return AstroCheck();
        }
        const func = UnitTests[name];
        if (typeof func !== 'function') {
            console.log(`test.js: Unknown unit test "${name}"`);
            return 1;
        }
        return func();
    }

    console.log('test.js: Invalid command line arguments.');
    return 1;
}


process.exit(main());
