'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function Fail(message) {
    console.log(`FATAL(mag_test.js): ${message}`);
    process.exit(1);
}

function CompareMatrices(caller, a, b, tolerance) {
    for (let i=0; i<3; ++i) {
        for (let j=0; j<3; ++j) {
            const diff = Math.abs(a.rot[i][j] - b.rot[i][j]);
            if (diff > tolerance) {
                throw `ERROR(${caller}): matrix[${i}][${j}] = ${a.rot[i][j]}, expected ${b.rot[i][j]}, diff ${diff}`;
            }
        }
    }
}

function Rotation_MatrixInverse() {
    const a = {rot:[
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ]};

    const v = {rot:[
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ]};

    const b = Astronomy.InverseRotation(a);
    CompareMatrices('Rotation_MatrixInverse', b, v, 0);
}

function Rotation_MatrixMultiply() {
    const a = {rot:[
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ]};

    const b = {rot:[
        [10, 13, 16],
        [11, 14, 17],
        [12, 15, 18]
    ]};

    const v = {rot:[
        [84, 201, 318],
        [90, 216, 342],
        [96, 231, 366]
    ]};

    const c = Astronomy.CombineRotation(b, a);
    CompareMatrices('Rotation_MatrixMultiply', c, v, 0);
}

function VectorDiff(a, b) {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    const dz = a.z - b.z;
    return Math.sqrt(dx*dx + dy*dy + dz*dz);
}

function Test_EQJ_ECL() {
    const r = Astronomy.Rotation_EQJ_ECL();

    /* Calculate heliocentric Earth position at a test time. */
    const time = Astronomy.MakeTime(new Date('2019-12-08T19:39:15Z'));
    const ev = Astronomy.HelioVector('Earth', time);

    /* Use the existing Astronomy.Ecliptic() to calculate ecliptic vector and angles. */
    const ecl = Astronomy.Ecliptic(ev.x, ev.y, ev.z);
    console.log(`Test_EQJ_ECL ecl = (${ecl.ex}, ${ecl.ey}, ${ecl.ez})`);

    /* Now compute the same vector via rotation matrix. */
    const ee = Astronomy.RotateVector(r, ev);
    const dx = ee.x - ecl.ex;
    const dy = ee.y - ecl.ey;
    const dz = ee.z - ecl.ez;
    const diff = Math.sqrt(dx*dx + dy*dy + dz*dz);
    console.log(`Test_EQJ_ECL ee = (${ee.x}, ${ee.y}, ${ee.z}); diff = ${diff}`);
    if (diff > 1.0e-16)
        throw 'Test_EQJ_ECL: EXCESSIVE VECTOR ERROR';

    /* Reverse the test: go from ecliptic back to equatorial. */
    const ir = Astronomy.Rotation_ECL_EQJ();
    const et = Astronomy.RotateVector(ir, ee);
    const idiff = VectorDiff(et, ev);
    console.log(`Test_EQJ_ECL ev diff = ${idiff}`);
    if (idiff > 1.0e-16)
        throw 'Test_EQJ_ECL: EXCESSIVE REVERSE ROTATION ERROR';
}

function Test_EQJ_EQD(body) {
    /* Verify conversion of equatorial J2000 to equatorial of-date, and back. */
    /* Use established functions to calculate spherical coordinates for the body, in both EQJ and EQD. */
    const time = Astronomy.MakeTime(new Date('2019-12-08T20:50:00Z'));
    const observer = Astronomy.MakeObserver(+35, -85, 0);
    const eq2000 = Astronomy.Equator(body, time, observer, false, true);
    const eqdate = Astronomy.Equator(body, time, observer, true, true);

    /* Convert EQJ spherical coordinates to vector. */
    const sphere = {lat:eq2000.dec, lon:15*eq2000.ra, dist:eq2000.dist};
    const v2000 = Astronomy.VectorFromSphere(sphere, time);

    /* Find rotation matrix. */
    const r = Astronomy.Rotation_EQJ_EQD(time);

    /* Rotate EQJ vector to EQD vector. */
    const vdate = Astronomy.RotateVector(r, v2000);

    /* Convert vector back to spherical coordinates. */
    let xsphere = Astronomy.SphereFromVector(vdate);

    /* Compare the result with the eqdate. */
    const ra_diff = Math.abs((xsphere.lon / 15) - eqdate.ra);
    const dec_diff = Math.abs(xsphere.lat - eqdate.dec);
    const dist_diff = Math.abs(xsphere.dist - eqdate.dist);
    console.log(`Test_EQJ_EQD: ${body} ra=${eqdate.ra}, dec=${eqdate.dec}, dist=${eqdate.dist}, ra_diff=${ra_diff}, dec_diff=${dec_diff}, dist_diff=${dist_diff}`);
    if (ra_diff > 1.0e-14 || dec_diff > 1.0e-14 || dist_diff > 4.0e-15)
        throw 'Test_EQJ_EQD: EXCESSIVE ERROR';

    /* Perform the inverse conversion back to equatorial J2000 coordinates. */
    const ir = Astronomy.Rotation_EQD_EQJ(time);
    const t2000 = Astronomy.RotateVector(ir, vdate);
    const diff = VectorDiff(t2000, v2000);
    console.log(`Test_EQJ_EQD: ${body} inverse diff = ${diff}`);
    if (diff > 3.0e-15)
        throw 'Test_EQJ_EQD: EXCESSIVE INVERSE ERROR';
}

function RotationTest() {
    Rotation_MatrixInverse();
    Rotation_MatrixMultiply();
    Test_EQJ_ECL();
    Test_EQJ_EQD('Mercury');
    Test_EQJ_EQD('Venus');
    Test_EQJ_EQD('Mars');
    Test_EQJ_EQD('Jupiter');
    Test_EQJ_EQD('Saturn');
}

RotationTest();
console.log('rotation.js: success');
process.exit(0);
