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

const DEG2RAD = 0.017453292519943296;

function ParseRightAscension(degText) {
    const hours = +(parseFloat(degText) / 15).toFixed(6);
    if (hours < 0 || hours >= 24) {
        throw `Invalid RA degrees=${degText} ==> hours=${hours}`;
    }
    return hours;
}

function ParseDeclination(text) {
    const dec = parseFloat(text);
    if (dec < -90 || dec > +90) {
        throw `Invalid DEC text "${text}"`;
    }
    return dec;
}

function ParseAzimuth(text) {
    const az = parseFloat(text);
    if (az < 0 || az >= 360) {
        throw `Invalid azimuth text "${text}"`;
    }
    return az;
}

function ParseAltitude(text) {
    const alt = parseFloat(text);
    if (alt < -90 || alt > +90) {
        throw `Invalid altitude text "${text}"`;
    }
    return alt;
}

function Fatal(context, message) {
    throw `FATAL(${context.fn} : ${context.lnum}): ${message}`;
}

function CompareEquatorial(context, type, jpl_ra, jpl_dec, calc_ra, calc_dec) {
    let ra_error = Math.abs(calc_ra - jpl_ra);
    if (ra_error > 12) {
        ra_error = 24 - ra_error;
    }
    ra_error *= 15 * 60 * Math.cos(jpl_dec * DEG2RAD);    // scale to arcminutes at this declination

    let dec_error = 60 * (calc_dec - jpl_dec);
    let arcmin = Math.sqrt(ra_error*ra_error + dec_error*dec_error);
    if (arcmin > context.arcmin_threshold) {
        console.log(`JPL ra=${jpl_ra}, dec=${jpl_dec}; CALC ra=${calc_ra}, dec=${calc_dec}; deltas=${ra_error}, ${dec_error}`);
        Fatal(context, `Excessive ${type} equatorial error = ${arcmin} arcmin`);
    }
    return arcmin;
}

function CompareHorizontal(context, jpl_az, jpl_alt, calc_az, calc_alt) {
    let az_error = Math.abs(calc_az - jpl_az);
    if (az_error > 180) {
        az_error = 360 - az_error;
    }
    az_error *= 60 * Math.cos(jpl_alt * DEG2RAD);   // scale to arcminutes at this altitude

    let alt_error = 60 * (calc_alt - jpl_alt);
    let arcmin = Math.sqrt(alt_error*alt_error + az_error*az_error);
    if (arcmin > context.arcmin_threshold) {
        Fatal(context, `Excessive horizontal error = ${arcmin} arcmin`);
    }
    return arcmin;
}

function ProcessRow(context, row) {
    // [ 1900-Jan-01 00:00 A   286.68085 -23.49769 285.80012 -23.36851 250.7166 -13.8716]
    // [ 2019-Aug-25 00:00 C   156.34214  11.05463 156.94501  11.15274 282.1727   1.0678]
    let m = row.match(/^\s(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s[\* ACN][ mrts]\s*(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)/);
    if (!m) {
        throw `Invalid line ${context.lnum} in file ${context.fn}: ${row}`;
    }
    const year = parseInt(m[1]);
    const month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].indexOf(m[2]);
    const day = parseInt(m[3]);
    const hour = parseInt(m[4]);
    const minute = parseInt(m[5]);
    const date = new Date(Date.UTC(year, month, day, hour, minute));
    
    const jpl = {
        m_ra:  ParseRightAscension(m[6]),
        m_dec: ParseDeclination(m[7]),
        a_ra:  ParseRightAscension(m[8]),
        a_dec: ParseDeclination(m[9]),
        az:    ParseAzimuth(m[10]),
        alt:   ParseAltitude(m[11])
    };

    const pos = Astronomy.GeoVector(context.body, date);
    const sky = Astronomy.SkyPos(pos, context.observer);
    const hor = Astronomy.Horizon(sky.t, context.observer, sky.ofdate.ra, sky.ofdate.dec);

    let arcmin = CompareEquatorial(context, 'metric', jpl.m_ra, jpl.m_dec, sky.j2000.ra, sky.j2000.dec);
    context.maxArcminError_MetricEquatorial = Math.max(context.maxArcminError_MetricEquatorial, arcmin);

    arcmin = CompareEquatorial(context, 'apparent', jpl.a_ra, jpl.a_dec, sky.ofdate.ra, sky.ofdate.dec);
    context.maxArcminError_ApparentEquatorial = Math.max(context.maxArcminError_ApparentEquatorial, arcmin);

    arcmin = CompareHorizontal(context, jpl.az, jpl.alt, hor.azimuth, hor.altitude);
    context.maxArcminError_Horizontal = Math.max(context.maxArcminError_Horizontal, arcmin);
}

function PrintSummary(context) {
    console.log(`${context.fn} : m_eq=${context.maxArcminError_MetricEquatorial.toFixed(3)}, a_eq=${context.maxArcminError_ApparentEquatorial.toFixed(3)}, hor=${context.maxArcminError_Horizontal.toFixed(3)}`);
}

function ProcessFile(inFileName) {
    const text = fs.readFileSync(inFileName, 'utf8');
    const lines = text.split('\n');
    var inHeader = true;
    var context = {
        observer: null,
        body: null,
        lnum: 0,
        fn: inFileName,
        arcmin_threshold: 1,
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
                //console.debug(`Observer: lat=${latitude}, lon=${longitude}, elev=${elevation}`);
            }

            // Target body name: Mars (499)                      {source: mar097}
            m = row.match(/^Target body name: ([A-Za-z]+)/);
            if (m) {
                context.body = m[1];
                //console.debug(`Body: ${context.body}`);
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
            return PrintSummary(context);
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
