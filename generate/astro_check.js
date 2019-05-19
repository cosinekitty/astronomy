var Astronomy = require('../source/js/astronomy.js');

var date = new Date('1700-01-01T00:00:00Z');
var stop = new Date('2200-01-01T00:00:00Z');
var body, pos, sky, hor, dt;
const observer = Astronomy.MakeObserver(29, -81, 10);

console.log(`o ${observer.latitude.toFixed(6)} ${observer.longitude.toFixed(6)} ${observer.height.toFixed(6)}`);

dt = (24*3600*1000) * (10 + Math.PI/100);       // 10.03141592... days; exercise different times of day
while (date < stop) {
    for (body of Astronomy.Bodies) {
        if (body !== 'Moon') {
            pos = Astronomy.HelioVector(body, date);
            console.log(`v ${body} ${pos.t.tt.toFixed(16)} ${pos.x.toFixed(16)} ${pos.y.toFixed(16)} ${pos.z.toFixed(16)}`);
    
            if (body !== 'Earth') {
                pos = Astronomy.GeoVector(body, date);
                sky = Astronomy.SkyPos(pos, observer);
                hor = Astronomy.Horizon(sky.t, observer, sky.ofdate.ra, sky.ofdate.dec);
                console.log(`s ${body} ${pos.t.tt.toFixed(16)} ${pos.t.ut.toFixed(16)} ${sky.j2000.ra.toFixed(16)} ${sky.j2000.dec.toFixed(16)} ${sky.j2000.dist.toFixed(16)} ${hor.azimuth.toFixed(16)} ${hor.altitude.toFixed(16)}`);
            }
        }
    }
    pos = Astronomy.GeoMoon(date);
    console.log(`v GM ${pos.t.tt.toFixed(16)} ${pos.x.toFixed(16)} ${pos.y.toFixed(16)} ${pos.z.toFixed(16)}`);

    sky = Astronomy.SkyPos(pos, observer);
    hor = Astronomy.Horizon(sky.t, observer, sky.ofdate.ra, sky.ofdate.dec);
    console.log(`s GM ${pos.t.tt.toFixed(16)} ${pos.t.ut.toFixed(16)} ${sky.j2000.ra.toFixed(16)} ${sky.j2000.dec.toFixed(16)} ${sky.j2000.dist.toFixed(16)} ${hor.azimuth.toFixed(16)} ${hor.altitude.toFixed(16)}`);

    date = new Date(date.getTime() + dt);
}
