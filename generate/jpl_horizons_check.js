'use strict';
/*
    jpl_horizons_check.js

    Compares ephemeris data from the online JPL Horizons tool
    with calculations performed by the JavaScript astronomy library.

    -----------------------------------------------------------------------------

    MIT License

    Copyright (c) 2019 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

const fs = require('fs');
const Astronomy = require('../source/js/astronomy.js');

function ProcessRow(context, row) {
    // [ 1900-Jan-01 00:00 A   286.68085 -23.49769 285.80012 -23.36851 250.7166 -13.8716]
    // [ 2019-Aug-25 00:00 C   156.34214  11.05463 156.94501  11.15274 282.1727   1.0678]
    let m = row.match(/^\s(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s[\* ACN][ mrts]\s*(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)/);
    if (!m) {
        throw `Invalid line ${context.lnum} in file ${context.inFileName}: ${row}`;
    }
    const year = parseInt(m[1]);
    const month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].indexOf(m[2]);
    const day = parseInt(m[3]);
    const hour = parseInt(m[4]);
    const minute = parseInt(m[5]);
    const date = new Date(Date.UTC(year, month, day, hour, minute));
    //console.log(date);
}

function ProcessFile(inFileName) {
    const text = fs.readFileSync(inFileName, 'utf8');
    const lines = text.split('\n');
    var inHeader = true;
    var context = {
        observer: null,
        body: null,
        lnum: 0,
        inFileName: inFileName,
        maxArcminError_MetricEquatorial: 0,
        maxArcminError_ApparentEquatorial: 0,
        maxArcminError_Horizontal: 0
    };
    for (var row of lines) {
        ++context.lnum;
        if (inHeader) {
            // Center geodetic : 279.000000,29.0000000,0.0100000 {E-lon(deg),Lat(deg),Alt(km)}
            let m = row.match(/^Center geodetic : ([^,]+),([^,]+),(\S+)/);
            if (m) {
                let longitude = parseFloat(m[1]);
                if (longitude > 180) {
                    longitude -= 360;
                }
                let latitude = parseFloat(m[2]);
                let elevation = 1000 * parseFloat(m[3]);
                context.observer = Astronomy.MakeObserver(latitude, longitude, elevation);
                console.debug(`Observer: lat=${latitude}, lon=${longitude}, elev=${elevation}`);
            }

            // Target body name: Mars (499)                      {source: mar097}
            m = row.match(/^Target body name: ([A-Za-z]+)/);
            if (m) {
                context.body = m[1];
                console.debug(`Body: ${context.body}`);
            }

            if (row === '$$SOE') {
                inHeader = false;
                if (!context.observer) {
                    throw `Missing observer data in file: ${inFileName}`;
                }
                if (!context.body) {
                    throw `Missing target body in file: ${inFileName}`;
                }
            }
        } else if (row === '$$EOE') {
            console.debug(`SUCCESS: ${inFileName}`);
            return;
        } else {
            ProcessRow(context, row);
        }
    }
    throw `Missing $$EOE token in file: ${inFileName}`;
}

if (process.argv.length !== 3) {
    console.log(`USAGE:  node ${process.argv[1]} infile`);
    process.exit(1);
}

ProcessFile(process.argv[2]);
process.exit(0);
