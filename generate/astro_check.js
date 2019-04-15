var Astronomy = require('../source/js/astronomy.js');

var date = new Date('1700-01-01T00:00:00Z');
var stop = new Date('2200-01-01T00:00:00Z');
var body, pos, sky, hor, dt;
const observer = Astronomy.MakeObserver(29, -81, 10);

console.log(`o ${observer.latitude} ${observer.longitude} ${observer.height}`);

dt = (24*3600*1000) * (10 + Math.PI/100);       // 10.03141592... days; exercise different times of day
while (date < stop) {
    for (body of Astronomy.Bodies) {
        if (body !== 'Moon') {
            pos = Astronomy.HelioVector(body, date);
            console.log(`v ${body} ${pos.t.tt} ${pos.x} ${pos.y} ${pos.z}`);
    
            if (body !== 'Earth') {
                pos = Astronomy.GeoVector(body, date);
                sky = Astronomy.SkyPos(pos, observer);
                hor = Astronomy.Horizon(sky.t, observer, sky.ofdate.ra, sky.ofdate.dec);
                console.log(`s ${body} ${pos.t.tt} ${pos.t.ut} ${sky.j2000.ra} ${sky.j2000.dec} ${sky.j2000.dist} ${hor.azimuth} ${hor.altitude}`);
            }
        }
    }
    pos = Astronomy.GeoMoon(date);
    console.log(`v GM ${pos.t.tt} ${pos.x} ${pos.y} ${pos.z}`);

    sky = Astronomy.SkyPos(pos, observer);
    hor = Astronomy.Horizon(sky.t, observer, sky.ofdate.ra, sky.ofdate.dec);
    console.log(`s GM ${pos.t.tt} ${pos.t.ut} ${sky.j2000.ra} ${sky.j2000.dec} ${sky.j2000.dist} ${hor.azimuth} ${hor.altitude}`);

    date = new Date(date.getTime() + dt);
}
