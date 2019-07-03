'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function Run() {
    const observer = Astronomy.MakeObserver(29, -81, 10);
    const time = Astronomy.MakeTime(-93692.7685882873047376);
    console.log(`time.ut     = ${time.ut.toFixed(16)}`);
    console.log(`time.tt     = ${time.tt.toFixed(16)}`);
    const body = 'Sun';
    const j2000 = Astronomy.Equator(body, time, observer, false, false);
    console.log(`j2000  ra   = ${j2000.ra.toFixed(16)}`);
    console.log(`j2000  dec  = ${j2000.dec.toFixed(16)}`);
    const ofdate = Astronomy.Equator(body, time, observer, true, true);
    console.log(`ofdate ra   = ${ofdate.ra.toFixed(16)}`);
    console.log(`ofdate dec  = ${ofdate.dec.toFixed(16)}`);
    const hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, false);
    console.log(`azimuth     = ${hor.azimuth.toFixed(16)}`);
    console.log(`altitude    = ${hor.altitude.toFixed(16)}`);
}

Run();
