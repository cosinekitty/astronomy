'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function SearchElongTest() {
    // Max elongation data obtained from:
    // http://www.skycaramba.com/greatest_elongations.shtml
    TestMaxElong('Mercury', '2010-01-17T05:22Z', '2010-01-27T05:22Z', 24.80, 'morning');
    TestMaxElong('Mercury', '2010-05-16T02:15Z', '2010-05-26T02:15Z', 25.10, 'morning');
    TestMaxElong('Mercury', '2010-09-09T17:24Z', '2010-09-19T17:24Z', 17.90, 'morning');
    TestMaxElong('Mercury', '2010-12-30T14:33Z', '2011-01-09T14:33Z', 23.30, 'morning');
    TestMaxElong('Mercury', '2011-04-27T19:03Z', '2011-05-07T19:03Z', 26.60, 'morning');
    TestMaxElong('Mercury', '2011-08-24T05:52Z', '2011-09-03T05:52Z', 18.10, 'morning');
    TestMaxElong('Mercury', '2011-12-13T02:56Z', '2011-12-23T02:56Z', 21.80, 'morning');
    TestMaxElong('Mercury', '2012-04-08T17:22Z', '2012-04-18T17:22Z', 27.50, 'morning');
    TestMaxElong('Mercury', '2012-08-06T12:04Z', '2012-08-16T12:04Z', 18.70, 'morning');
    TestMaxElong('Mercury', '2012-11-24T22:55Z', '2012-12-04T22:55Z', 20.60, 'morning');
    TestMaxElong('Mercury', '2013-03-21T22:02Z', '2013-03-31T22:02Z', 27.80, 'morning');
    TestMaxElong('Mercury', '2013-07-20T08:51Z', '2013-07-30T08:51Z', 19.60, 'morning');
    TestMaxElong('Mercury', '2013-11-08T02:28Z', '2013-11-18T02:28Z', 19.50, 'morning');
    TestMaxElong('Mercury', '2014-03-04T06:38Z', '2014-03-14T06:38Z', 27.60, 'morning');
    TestMaxElong('Mercury', '2014-07-02T18:22Z', '2014-07-12T18:22Z', 20.90, 'morning');
    TestMaxElong('Mercury', '2014-10-22T12:36Z', '2014-11-01T12:36Z', 18.70, 'morning');
    TestMaxElong('Mercury', '2015-02-14T16:20Z', '2015-02-24T16:20Z', 26.70, 'morning');
    TestMaxElong('Mercury', '2015-06-14T17:10Z', '2015-06-24T17:10Z', 22.50, 'morning');
    TestMaxElong('Mercury', '2015-10-06T03:20Z', '2015-10-16T03:20Z', 18.10, 'morning');
    TestMaxElong('Mercury', '2016-01-28T01:22Z', '2016-02-07T01:22Z', 25.60, 'morning');
    TestMaxElong('Mercury', '2016-05-26T08:45Z', '2016-06-05T08:45Z', 24.20, 'morning');
    TestMaxElong('Mercury', '2016-09-18T19:27Z', '2016-09-28T19:27Z', 17.90, 'morning');
    TestMaxElong('Mercury', '2017-01-09T09:42Z', '2017-01-19T09:42Z', 24.10, 'morning');
    TestMaxElong('Mercury', '2017-05-07T23:19Z', '2017-05-17T23:19Z', 25.80, 'morning');
    TestMaxElong('Mercury', '2017-09-02T10:14Z', '2017-09-12T10:14Z', 17.90, 'morning');
    TestMaxElong('Mercury', '2017-12-22T19:48Z', '2018-01-01T19:48Z', 22.70, 'morning');
    TestMaxElong('Mercury', '2018-04-19T18:17Z', '2018-04-29T18:17Z', 27.00, 'morning');
    TestMaxElong('Mercury', '2018-08-16T20:35Z', '2018-08-26T20:35Z', 18.30, 'morning');
    TestMaxElong('Mercury', '2018-12-05T11:34Z', '2018-12-15T11:34Z', 21.30, 'morning');
    TestMaxElong('Mercury', '2019-04-01T19:40Z', '2019-04-11T19:40Z', 27.70, 'morning');
    TestMaxElong('Mercury', '2019-07-30T23:08Z', '2019-08-09T23:08Z', 19.00, 'morning');
    TestMaxElong('Mercury', '2019-11-18T10:31Z', '2019-11-28T10:31Z', 20.10, 'morning');
    TestMaxElong('Mercury', '2010-03-29T23:32Z', '2010-04-08T23:32Z', 19.40, 'evening');
    TestMaxElong('Mercury', '2010-07-28T01:03Z', '2010-08-07T01:03Z', 27.40, 'evening');
    TestMaxElong('Mercury', '2010-11-21T15:42Z', '2010-12-01T15:42Z', 21.50, 'evening');
    TestMaxElong('Mercury', '2011-03-13T01:07Z', '2011-03-23T01:07Z', 18.60, 'evening');
    TestMaxElong('Mercury', '2011-07-10T04:56Z', '2011-07-20T04:56Z', 26.80, 'evening');
    TestMaxElong('Mercury', '2011-11-04T08:40Z', '2011-11-14T08:40Z', 22.70, 'evening');
    TestMaxElong('Mercury', '2012-02-24T09:39Z', '2012-03-05T09:39Z', 18.20, 'evening');
    TestMaxElong('Mercury', '2012-06-21T02:00Z', '2012-07-01T02:00Z', 25.70, 'evening');
    TestMaxElong('Mercury', '2012-10-16T21:59Z', '2012-10-26T21:59Z', 24.10, 'evening');
    TestMaxElong('Mercury', '2013-02-06T21:24Z', '2013-02-16T21:24Z', 18.10, 'evening');
    TestMaxElong('Mercury', '2013-06-02T16:45Z', '2013-06-12T16:45Z', 24.30, 'evening');
    TestMaxElong('Mercury', '2013-09-29T09:59Z', '2013-10-09T09:59Z', 25.30, 'evening');
    TestMaxElong('Mercury', '2014-01-21T10:00Z', '2014-01-31T10:00Z', 18.40, 'evening');
    TestMaxElong('Mercury', '2014-05-15T07:06Z', '2014-05-25T07:06Z', 22.70, 'evening');
    TestMaxElong('Mercury', '2014-09-11T22:20Z', '2014-09-21T22:20Z', 26.40, 'evening');
    TestMaxElong('Mercury', '2015-01-04T20:26Z', '2015-01-14T20:26Z', 18.90, 'evening');
    TestMaxElong('Mercury', '2015-04-27T04:46Z', '2015-05-07T04:46Z', 21.20, 'evening');
    TestMaxElong('Mercury', '2015-08-25T10:20Z', '2015-09-04T10:20Z', 27.10, 'evening');
    TestMaxElong('Mercury', '2015-12-19T03:11Z', '2015-12-29T03:11Z', 19.70, 'evening');
    TestMaxElong('Mercury', '2016-04-08T14:00Z', '2016-04-18T14:00Z', 19.90, 'evening');
    TestMaxElong('Mercury', '2016-08-06T21:24Z', '2016-08-16T21:24Z', 27.40, 'evening');
    TestMaxElong('Mercury', '2016-12-01T04:36Z', '2016-12-11T04:36Z', 20.80, 'evening');
    TestMaxElong('Mercury', '2017-03-22T10:24Z', '2017-04-01T10:24Z', 19.00, 'evening');
    TestMaxElong('Mercury', '2017-07-20T04:34Z', '2017-07-30T04:34Z', 27.20, 'evening');
    TestMaxElong('Mercury', '2017-11-14T00:32Z', '2017-11-24T00:32Z', 22.00, 'evening');
    TestMaxElong('Mercury', '2018-03-05T15:07Z', '2018-03-15T15:07Z', 18.40, 'evening');
    TestMaxElong('Mercury', '2018-07-02T05:24Z', '2018-07-12T05:24Z', 26.40, 'evening');
    TestMaxElong('Mercury', '2018-10-27T15:25Z', '2018-11-06T15:25Z', 23.30, 'evening');
    TestMaxElong('Mercury', '2019-02-17T01:23Z', '2019-02-27T01:23Z', 18.10, 'evening');
    TestMaxElong('Mercury', '2019-06-13T23:14Z', '2019-06-23T23:14Z', 25.20, 'evening');
    TestMaxElong('Mercury', '2019-10-10T04:00Z', '2019-10-20T04:00Z', 24.60, 'evening');
    TestMaxElong('Venus', '2010-12-29T15:57Z', '2011-01-08T15:57Z', 47.00, 'morning');
    TestMaxElong('Venus', '2012-08-05T08:59Z', '2012-08-15T08:59Z', 45.80, 'morning');
    TestMaxElong('Venus', '2014-03-12T19:25Z', '2014-03-22T19:25Z', 46.60, 'morning');
    TestMaxElong('Venus', '2015-10-16T06:57Z', '2015-10-26T06:57Z', 46.40, 'morning');
    TestMaxElong('Venus', '2017-05-24T13:09Z', '2017-06-03T13:09Z', 45.90, 'morning');
    TestMaxElong('Venus', '2018-12-27T04:24Z', '2019-01-06T04:24Z', 47.00, 'morning');
    TestMaxElong('Venus', '2010-08-10T03:19Z', '2010-08-20T03:19Z', 46.00, 'evening');
    TestMaxElong('Venus', '2012-03-17T08:03Z', '2012-03-27T08:03Z', 46.00, 'evening');
    TestMaxElong('Venus', '2013-10-22T08:00Z', '2013-11-01T08:00Z', 47.10, 'evening');
    TestMaxElong('Venus', '2015-05-27T18:46Z', '2015-06-06T18:46Z', 45.40, 'evening');
    TestMaxElong('Venus', '2017-01-02T13:19Z', '2017-01-12T13:19Z', 47.10, 'evening');
    TestMaxElong('Venus', '2018-08-07T17:02Z', '2018-08-17T17:02Z', 45.90, 'evening');
}

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
        Astronomy.ResetPerformanceMetrics();
        let time = Astronomy.SearchRelativeLongitude(item.body, targetRelLon, startDate);
        const metrics = Astronomy.GetPerformanceMetrics();
        let diff_minutes = (time.date - item.date) / 60000;
        console.log(`${item.body}: error = ${diff_minutes.toFixed(3)} minutes, iterations = ${metrics.longitude_iter}`);
        if (Math.abs(diff_minutes) > 15)
            throw `!!! Excessive error for body ${item.body}`;
    }
}

function TestPlanet(outFileName, body, startYear, stopYear, zeroLonEventName) {
    let rlon = 0;
    let date = new Date(Date.UTC(startYear, 0, 1));
    let stopDate = new Date(Date.UTC(stopYear, 0, 1));
    let text = '';
    let count = 0;
    let prev_time, min_diff, max_diff, sum_diff=0;

    Astronomy.ResetPerformanceMetrics();
    while (date < stopDate) {
        let event = (rlon === 0) ? zeroLonEventName : 'sup';
        let evt_time = Astronomy.SearchRelativeLongitude(body, rlon, date);
        if (prev_time) {
            // Check for consistent intervals.
            // Mainly I don't want to accidentally skip over an event!
            let day_diff = evt_time.tt - prev_time.tt;
            if (min_diff === undefined) {
                min_diff = max_diff = day_diff;
            } else {
                min_diff = Math.min(min_diff, day_diff);
                max_diff = Math.max(max_diff, day_diff);
            }
            sum_diff += day_diff;
        }
        let geo = Astronomy.GeoVector(body, evt_time);
        let dist = Math.sqrt(geo.x*geo.x + geo.y*geo.y + geo.z*geo.z);
        text += `e ${body} ${event} ${evt_time.tt} ${dist}\n`;
        rlon = 180 - rlon;
        date = evt_time.date;
        ++count;
        prev_time = evt_time;
    }
    const metrics = Astronomy.GetPerformanceMetrics();

    fs.writeFileSync(outFileName, text);

    const ratio = max_diff / min_diff;
    const iter_per_call = metrics.longitude_iter / metrics.longitude_search;
    console.log(`TestPlanet(${body}): ${count} events, ${iter_per_call.toFixed(3)} iter/call, interval min=${min_diff.toFixed(1)}, max=${max_diff.toFixed(1)}, avg=${(sum_diff/count).toFixed(1)}, ratio=${ratio.toFixed(3)}`);

    let thresh = {Mercury:1.65, Mars:1.30}[body] || 1.07;
    if (ratio > thresh)
        throw `TestPlanet: Excessive event interval ratio for ${body} = ${ratio}`;
}

function TestMaxElong(body, startText, verifyText, verifyAngle, verifyVisibility) {
    let startDate = new Date(startText);
    let verifyDate = new Date(verifyText);
    let evt = Astronomy.SearchMaxElongation(body, startDate);

    let hour_diff = (verifyDate - evt.time.date) / (1000 * 3600);
    let arcmin_diff = 60.0 * Math.abs(evt.elongation - verifyAngle);
    console.log(`TestMaxElong: ${body.padStart(8)} ${evt.visibility.padStart(8)} elong=${evt.elongation.toFixed(2).padStart(5)} (${arcmin_diff.toFixed(2).padStart(4)} arcmin)  ${evt.time.toString()} (err ${hour_diff.toFixed(2).padStart(5)} hours)`);

    if (evt.visibility !== verifyVisibility)
        throw `TestMaxElong: expected visibility ${verifyVisibility}, but found ${evt.visibility}`;

    if (arcmin_diff > 4.0)
        throw `TestMaxElong: excessive angular error = ${angle_diff} arcmin`;

    if (Math.abs(hour_diff) > 0.603)
        throw `TestMaxElong: excessive hour error = ${hour_diff}`;
}

console.log('elong_test.js: Starting');
TestFile('longitude/opposition_2018.txt', 2018, 0);

for (let body of ['Mercury', 'Venus'])
    TestPlanet(`temp/js_longitude_${body}.txt`, body, 1700, 2200, 'inf');

for (let body of ['Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
    TestPlanet(`temp/js_longitude_${body}.txt`, body, 1700, 2200, 'opp');

SearchElongTest();

console.log('elong_test.js: SUCCESS');
process.exit(0);
