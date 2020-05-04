'use strict';
const fs = require('fs');
const Astronomy = require('../source/js/astronomy.min.js');

function Fail(message) {
    console.error(`FATAL(constellation.js): ${message}`);
    process.exit(1);
}

function ConstellationTest() {
    const filename = 'constellation/test_input.txt';
    const text = fs.readFileSync(filename, {encoding:'utf8'});
    const lines = text.trimRight().split('\n');
    let lnum = 0;
    let failcount = 0;
    for (let row of lines) {
        ++lnum;
        let token = row.trim().split(/\s+/);
        if (token.length !== 4) {
            Fail(`Bad data in ${filename} line ${lnum}: found ${token.length} tokens.`);
        }
        const id = parseInt(token[0]);
        const ra = parseFloat(token[1]);
        const dec = parseFloat(token[2]);
        const symbol = token[3];
        if (!/^[A-Z][A-Za-z]{2}$/.test(symbol)) {
            Fail(`Invalid symbol "${symbol}" in ${filename} line ${lnum}`);
        }
        const constel = Astronomy.Constellation(ra, dec);
        if (constel.symbol !== symbol) {
            ++failcount;
            console.error(`Star ${id}: expected ${symbol}, found ${constel.symbol} at B1875 RA=${constel.ra1875}, DEC=${constel.dec1875}`);
        }
    }
    if (failcount > 0) {
        Fail(`constellation.js: ${failcount} failures.`);
    }
    console.log(`constellation.js: PASS (verified ${lines.length})`);
}

ConstellationTest();
process.exit(0);
