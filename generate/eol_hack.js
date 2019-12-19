/*
    eol_hack.js  -  Don Cross

    Fixes line endings in a text file to be CR LF pairs.
    This is a hack so that there is no diff noise
    caused by running jsdoc2md from Windows.

    See open bug at:
    https://github.com/jsdoc2md/jsdoc-to-markdown/issues/112
*/
'use strict';
const fs = require('fs');

function FixFile(filename) {
    const inText = fs.readFileSync(filename, 'utf8');
    const outText = inText.replace(/\r?\n|\r\n?/g, '\r\n');

    if (inText === outText) {
        console.log(`eol_hack.js: Leaving file as-is: ${filename}`);
    } else {
        fs.writeFileSync(filename, outText, 'utf8');
        console.log(`eol_hack.js: Rewrote line endings: ${filename}`);
    }
}

if (process.argv.length !== 3) {
    console.log('USAGE: eol_hack.js filename');
    process.exit(1);
}

FixFile(process.argv[2]);
process.exit(0);
