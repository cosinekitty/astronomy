'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function Fail(message) {
    console.log(`FATAL(seasons_test.js): ${message}`);
    process.exit(1);
}

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

    console.log(`LoadTestData: count = ${data.length}`);
    for (let month of [3, 6, 9, 12]) {
        console.log(`Month ${month}: earliest ${minByMonth[month]}, latest ${maxByMonth[month]}`);
    }
    return data;
}

function Test() {
    const data = LoadTestData('seasons/seasons.txt');
    let index = 0;
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
        if (Math.abs(diff_minutes) > 1.7) {
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

    console.log(`seasons_test: n=${count}, minute errors: min=${min_diff.toFixed(3)}, avg=${(sum_diff/count).toFixed(3)}, max=${max_diff.toFixed(3)}`);
    for (let month of [3, 6, 9, 12]) {
        console.log(`seasons_test: max diff by month ${month} = ${month_max_diff[month].toFixed(3)}`);
    }
}

console.log(`seasons_test: Beginning.`);
Test();
console.log(`seasons_test: Finished.`);
process.exit(0);