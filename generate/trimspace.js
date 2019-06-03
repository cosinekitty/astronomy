/*
    trimspace.js  -  Don Cross

    Trim all trailing whitespace from every line in a text file.
    This is to help me automatically clean up source code before
    feeding it through doxygen, because that has caused issues
    in Windows builds.
*/
'use strict';
const fs = require('fs');

function FixFile(filename) {
    const inText = fs.readFileSync(filename, 'utf8');
    const outText = inText.replace(/[ \t]+(\r\n?)/g, '$1');
    
    if (inText === outText) {
        console.log(`trimspace.js: Leaving file as-is: ${filename}`);
    } else {
        fs.writeFileSync(filename, outText, 'utf8');        
        console.log(`trimspace.js: Removed trailing whitespace from file: ${filename}`);
    }
}

if (process.argv.length !== 3) {
    console.log('USAGE: trimspace.js filename');
    process.exit(1);
}

FixFile(process.argv[2]);
process.exit(0);
