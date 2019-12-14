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

function RotationTest() {
    Rotation_MatrixInverse();
    Rotation_MatrixMultiply();
}

RotationTest();
console.log('rotation.js: success');
process.exit(0);
