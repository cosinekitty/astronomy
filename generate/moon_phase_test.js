const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function TestMoonPhases() {
    const text = fs.readFileSync('moonphase/moonphases.txt', {encoding:'utf8'});
    const lines = text.trimRight().split('\n');
    for (let row of lines) {        
        let token = row.split(' ');
        let phase = parseInt(token[0]);
        let date = new Date(token[1]);
        let gm = Astronomy.GeoVector('Moon', date);
        const moon_eclip = Astronomy.Ecliptic(gm.x, gm.y, gm.z);

        let sun = Astronomy.GeoVector('Sun', date);
        const sun_eclip = Astronomy.Ecliptic(sun.x, sun.y, sun.z);

        let elong = moon_eclip.elon - sun_eclip.elon;
        if (elong < 0) elong += 360;
        console.log(`${phase} ${date.toISOString()} ${elong}`);
    }
}

TestMoonPhases();