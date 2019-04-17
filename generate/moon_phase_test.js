const Astronomy = require('../source/js/astronomy.js');

function Test1() {
    const when = new Date('2019-04-17T00:00:00Z');
    const gm = Astronomy.GeoVector('Moon', when);
    const eclip = Astronomy.Ecliptic(gm.x, gm.y, gm.z);
    console.log(eclip);
}

Test1();