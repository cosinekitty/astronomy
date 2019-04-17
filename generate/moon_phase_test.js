const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function LoadMoonPhaseData(filename) {
    // Load known moon phase times from US Naval Observatory.
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.trimRight().split('\n');
    let data = [];
    for (let row of lines) {        
        let token = row.split(' ');
        data.push({phase:parseInt(token[0]), date:new Date(token[1])});
    }
    return data;
}

function TestLongitudes(data) {
    // Using known moon phase times from US Naval Obs
    let max_arcmin = 0;
    for (let row of data) {
        let gm = Astronomy.GeoVector('Moon', row.date);
        const moon_eclip = Astronomy.Ecliptic(gm.x, gm.y, gm.z);

        let sun = Astronomy.GeoVector('Sun', row.date);
        const sun_eclip = Astronomy.Ecliptic(sun.x, sun.y, sun.z);

        let elong = moon_eclip.elon - sun_eclip.elon;
        if (elong < 0) elong += 360;
        let expected_elong = 90 * row.phase;
        let degree_error = Math.abs(elong - expected_elong);
        if (degree_error > 180) degree_error = 360 - degree_error;
        let arcmin = 60 * degree_error;
        max_arcmin = Math.max(max_arcmin, arcmin);
        //console.log(`${row.phase} ${row.date.toISOString()} ${elong} ${arcmin}`);
    }
    console.log(`TestLongitudes: nrows = ${data.length}, max_arcmin = ${max_arcmin}`);
    if (max_arcmin > 1.0) {
        console.log(`TestLongitudes: EXCESSIVE ANGULAR ERROR`);
        return 1;
    }
    return 0;
}

function TestSearch(data) {
    return 0;
}

const TestData = LoadMoonPhaseData('moonphase/moonphases.txt')
process.exit(TestLongitudes(TestData) || TestSearch(TestData));
