/*
    moonphase.js  -  by Don Cross - 2019-05-13

    Example Node.js program for Astronomy Engine:
    https://cosinekitty.github.io/astronomy/

    To execute, run the command:
    node moonphase
*/

const Astronomy = require('../../source/js/astronomy.js');      // adjust path as needed for your system

function Pad(s, w) {
    s = s.toFixed(0);
    while (s.length < w) {
        s = '0' + s;
    }
    return s;
}

function FormatDate(date) {
    var year = Pad(date.getFullYear(), 4);
    var month = Pad(1 + date.getMonth(), 2);
    var day = Pad(date.getDate(), 2);
    var hour = Pad(date.getHours(), 2);
    var minute = Pad(date.getMinutes(), 2);
    var second = Pad(date.getSeconds(), 2);
    return `${year}-${month}-${day} ${hour}:${minute}:${second}`;
}

function Demo() {
    const now = new Date();

    // Calculate the Moon's current phase angle, 
    // which ranges from 0 to 360 degrees.
    //   0 degrees = new moon,
    //  90 degrees = first quarter,
    // 180 degrees = full moon,
    // 270 degrees = third quarter.
    const phase = Astronomy.MoonPhase(now);     // https://cosinekitty.github.io/astronomy/source/js/#Astronomy.MoonPhase
    console.log(`${FormatDate(now)} : Moon's phase angle = ${phase.toFixed(6)} degrees.`);        
    console.log('');

    // Now we predict when the next 10 lunar quarter phases will happen.
    console.log('The next 10 lunar quarters are:');
    const QuarterName = ['New Moon', 'First Quarter', 'Full Moon', 'Third Quarter'];
    let mq;
    for (let i=0; i < 10; ++i) {
        if (mq === undefined) {
            // The first time around the for loop, we search forward
            // from the current date and time to find the next quarter
            // phase, whatever it might be.
            // https://cosinekitty.github.io/astronomy/source/js/#Astronomy.SearchMoonQuarter
            mq = Astronomy.SearchMoonQuarter(now);

            // The object now stored in 'mq' is documented here:
            // https://cosinekitty.github.io/astronomy/source/js/#Astronomy.MoonQuarter
        } else {
            // Use the previous moon quarter information to find the next quarter phase event.
            // https://cosinekitty.github.io/astronomy/source/js/#Astronomy.NextMoonQuarter
            mq = Astronomy.NextMoonQuarter(mq);
        }
        console.log(`${FormatDate(mq.time.date)} : ${QuarterName[mq.quarter]}`);
    }
}

Demo();
