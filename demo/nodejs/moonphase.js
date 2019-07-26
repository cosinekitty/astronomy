/*
    moonphase.js  -  by Don Cross - 2019-05-13

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the Moon's phase for a given date and time,
    or for the computer's current date and time if none is given.
    It also finds the dates and times of the subsequent 10 quarter phase changes.

    To execute, run the command:
    node moonphase [date]
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
    var year = Pad(date.getUTCFullYear(), 4);
    var month = Pad(1 + date.getUTCMonth(), 2);
    var day = Pad(date.getUTCDate(), 2);
    var hour = Pad(date.getUTCHours(), 2);
    var minute = Pad(date.getUTCMinutes(), 2);
    var svalue = date.getUTCSeconds() + (date.getUTCMilliseconds() / 1000);
    var second = Pad(Math.round(svalue), 2);
    return `${year}-${month}-${day} ${hour}:${minute}:${second} UTC`;
}

function Demo() {
    const date = (process.argv.length === 3) ? new Date(process.argv[2]) : new Date();

    // Calculate the Moon's current phase angle, 
    // which ranges from 0 to 360 degrees.
    //   0 degrees = new moon,
    //  90 degrees = first quarter,
    // 180 degrees = full moon,
    // 270 degrees = third quarter.
    const phase = Astronomy.MoonPhase(date);
    console.log(`${FormatDate(date)} : Moon's phase angle = ${phase.toFixed(6)} degrees.`);        
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
            mq = Astronomy.SearchMoonQuarter(date);
        } else {
            // Use the previous moon quarter information to find the next quarter phase event.
            mq = Astronomy.NextMoonQuarter(mq);
        }
        console.log(`${FormatDate(mq.time.date)} : ${QuarterName[mq.quarter]}`);
    }
}

Demo();
process.exit(0);
