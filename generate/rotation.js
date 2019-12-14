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
    const a = Astronomy.MakeRotation([
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ]);

    const v = Astronomy.MakeRotation([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ]);

    const b = Astronomy.InverseRotation(a);
    CompareMatrices('Rotation_MatrixInverse', b, v, 0);
}

function Rotation_MatrixMultiply() {
    const a = Astronomy.MakeRotation([
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ]);

    const b = Astronomy.MakeRotation([
        [10, 13, 16],
        [11, 14, 17],
        [12, 15, 18]
    ]);

    const v = Astronomy.MakeRotation([
        [84, 201, 318],
        [90, 216, 342],
        [96, 231, 366]
    ]);

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
    const v2000 = Astronomy.VectorFromEquator(eq2000, time);

    /* Find rotation matrix. */
    const r = Astronomy.Rotation_EQJ_EQD(time);

    /* Rotate EQJ vector to EQD vector. */
    const vdate = Astronomy.RotateVector(r, v2000);

    /* Convert vector back to angular equatorial coordinates. */
    let equcheck = Astronomy.EquatorFromVector(vdate);

    /* Compare the result with the eqdate. */
    const ra_diff = Math.abs(equcheck.ra - eqdate.ra);
    const dec_diff = Math.abs(equcheck.dec - eqdate.dec);
    const dist_diff = Math.abs(equcheck.dist - eqdate.dist);
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

function Test_EQD_HOR(body) {
    /* Use existing functions to calculate horizontal coordinates of the body for the time+observer. */
    const time = Astronomy.MakeTime(new Date('1970-12-13T05:15:00Z'));
    const observer = Astronomy.MakeObserver(-37, +45, 0);
    const eqd = Astronomy.Equator(body, time, observer, true, true);
    console.log(`Test_EQD_HOR ${body}: OFDATE ra=${eqd.ra}, dec=${eqd.dec}`);
    const hor = Astronomy.Horizon(time, observer, eqd.ra, eqd.dec, 'normal');

    /* Calculate the position of the body as an equatorial vector of date. */
    const vec_eqd = Astronomy.VectorFromEquator(eqd, time);

    /* Calculate rotation matrix to convert equatorial J2000 vector to horizontal vector. */
    const rot = Astronomy.Rotation_EQD_HOR(time, observer);

    /* Rotate the equator of date vector to a horizontal vector. */
    const vec_hor = Astronomy.RotateVector(rot, vec_eqd);

    /* Convert the horizontal vector to horizontal angular coordinates. */
    const xsphere = Astronomy.HorizonFromVector(vec_hor, 'normal');
    const diff_alt = Math.abs(xsphere.lat - hor.altitude);
    const diff_az = Math.abs(xsphere.lon - hor.azimuth);

    console.log(`Test_EQD_HOR ${body}: trusted alt=${hor.altitude}, az=${hor.azimuth}; test alt=${xsphere.lat}, az=${xsphere.lon}; diff_alt=${diff_alt}, diff_az=${diff_az}`);
    if (diff_alt > 4.0e-14 || diff_az > 1.0e-13)
        throw 'Test_EQD_HOR: EXCESSIVE HORIZONTAL ERROR.';

    /* Confirm that we can convert back to horizontal vector. */
    const check_hor = Astronomy.VectorFromHorizon(xsphere, time, 'normal');
    let diff = VectorDiff(check_hor, vec_hor);
    console.log(`Test_EQD_HOR ${body}: horizontal recovery: diff = ${diff}`);
    if (diff > 1.0e-15)
        throw 'Test_EQD_HOR: EXCESSIVE ERROR IN HORIZONTAL RECOVERY.';

    /* Verify the inverse translation from horizontal vector to equatorial of-date vector. */
    const irot = Astronomy.Rotation_HOR_EQD(time, observer);
    const check_eqd = Astronomy.RotateVector(irot, vec_hor);
    diff = VectorDiff(check_eqd, vec_eqd);
    console.log(`Test_EQD_HOR ${body}: OFDATE inverse rotation diff = ${diff}`);
    if (diff > 1.0e-15)
        throw 'Test_EQD_HOR: EXCESSIVE OFDATE INVERSE HORIZONTAL ERROR.';

    /* Exercise HOR to EQJ translation. */
    const eqj = Astronomy.Equator(body, time, observer, false, true);
    const vec_eqj = Astronomy.VectorFromEquator(eqj, time);
    const yrot = Astronomy.Rotation_HOR_EQJ(time, observer);
    const check_eqj = Astronomy.RotateVector(yrot, vec_hor);
    diff = VectorDiff(check_eqj, vec_eqj);
    console.log(`Test_EQD_HOR ${body}: J2000 inverse rotation diff = ${diff}`);
    if (diff > 3.0e-15)
        throw 'Test_EQD_HOR: EXCESSIVE J2000 INVERSE HORIZONTAL ERROR.';

    /* Verify the inverse translation: EQJ to HOR. */
    const zrot = Astronomy.Rotation_EQJ_HOR(time, observer);
    const another_hor = Astronomy.RotateVector(zrot, vec_eqj);
    diff = VectorDiff(another_hor, vec_hor);
    console.log(`Test_EQD_HOR ${body}: EQJ inverse rotation diff = ${diff}`);
    if (diff > 2.0e-15)
        throw 'Test_EQD_HOR: EXCESSIVE EQJ INVERSE HORIZONTAL ERROR.';
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

    Test_EQD_HOR('Mercury');
    Test_EQD_HOR('Venus');
    Test_EQD_HOR('Mars');
    Test_EQD_HOR('Jupiter');
    Test_EQD_HOR('Saturn');
}

RotationTest();
console.log('rotation.js: success');
process.exit(0);
