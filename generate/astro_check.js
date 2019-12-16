'use strict';
var Astronomy = require('../source/js/astronomy.min.js');

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

            if (body !== 'Earth') {
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
