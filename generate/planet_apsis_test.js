/*
    planet_apsis_test.js

    Exercises finding planet perihelion/aphelion.
*/

'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');
const DebugMode = (process.argv.length > 2 && process.argv[2] === '-d');

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

function Test() {
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
        if (DebugMode) console.log(`JS PlanetApsis: ${count} apsides for ${body} -- intervals: min=${min_interval}, max=${max_interval}, ratio=${max_interval/min_interval}; max day=${max_diff_days}, degrees=${(max_diff_days / period) * 360}, dist ratio=${max_dist_ratio}`);
        ++pindex;
    }
    if (found_bad_planet) {
        throw `APSIS FAIL: Planet(s) exceeded angular threshold (${degree_threshold} degrees).`;
    }
}

Test();
console.log('planet_apsis_test.js: PASS');
process.exit(0);
