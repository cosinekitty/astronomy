'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function TestLunarEclipse() {
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

        // Check eclipse center.
        let diff_days = eclipse.center.ut - peak_time.ut;

        // Tolerate missing penumbral eclipses - skip to next input line without calculating next eclipse.
        if (partial_minutes == 0.0 && diff_days > 20.0) {
            ++skip_count;
            continue;
        }

        let diff_minutes = (24.0 * 60.0) * Math.abs(diff_days);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit) {
            console.error(`JS LunarEclipseTest expected center: ${peak_time}`);
            console.error(`JS LunarEclipseTest found    center: ${eclipse.center}`);
            console.error(`JS LunarEclipseTest(${filename} line ${lnum}): EXCESSIVE center time error = ${diff_minutes} minutes (${diff_days} days).`);
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

        eclipse = Astronomy.NextLunarEclipse(eclipse.center);
    }
    console.log(`JS LunarEclipseTest: PASS (verified ${lnum}, skipped ${skip_count}, max_diff_minutes = ${max_diff_minutes}, avg_diff_minutes = ${sum_diff_minutes / diff_count}, moon calcs = ${Astronomy.CalcMoonCount})`);
    return 0;
}

process.exit(TestLunarEclipse());
