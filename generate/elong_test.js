'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function LoadData(filename) {
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.trimRight().split('\n');
    let data = [];
    for (let row of lines) {        
        let token = row.split(/\s+/);
        data.push({date:new Date(token[0]), body:token[1]});
    }
    return data;
}

function TestFile(filename, startSearchYear, targetRelLon) {
    const data = LoadData(filename);
    const startDate = new Date(Date.UTC(startSearchYear, 0, 1));
    for (let item of data) {
        let time = Astronomy.SearchRelativeLongitude(item.body, targetRelLon, startDate);
        let diff_minutes = (time.date - item.date) / 60000;
        console.log(`${item.body}: error = ${diff_minutes.toFixed(3)} minutes`);
        if (Math.abs(diff_minutes) > 15)
            throw `!!! Excessive error for body ${item.body}`;
    }
}

function TestSuperiorPlanet(outFileName, body, startYear, stopYear) {
    // e Jupiter opp <tt> <au>
    // e Jupiter sup <tt> <au>
    let rlon = 0;   // start with opposition, then alternate between conjunction and opposition
    let date = new Date(Date.UTC(startYear, 0, 1));
    let stopDate = new Date(Date.UTC(stopYear, 0, 1));
    let text = '';
    let count = 0;

    while (date < stopDate) {
        let event = (rlon === 0) ? 'opp' : 'sup';
        let evt_time = Astronomy.SearchRelativeLongitude(body, rlon, date);
        let geo = Astronomy.GeoVector(body, evt_time);
        let dist = Math.sqrt(geo.x*geo.x + geo.y*geo.y + geo.z*geo.z);
        text += `e ${body} ${event} ${evt_time.tt} ${dist}\n`;
        rlon = 180 - rlon;
        date = evt_time.date;
        ++count;
    }

    fs.writeFileSync(outFileName, text);
    console.log(`TestSuperiorPlanet(${body}): wrote ${count} events to file ${outFileName}`);
}

console.log('elong_test.js: Starting');
TestFile('longitude/opposition_2018.txt', 2018, 0);
for (let body of ['Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
    TestSuperiorPlanet(`temp/longitude_${body}.txt`, body, 1700, 2200);

console.log('elong_test.js: SUCCESS')
process.exit(0);