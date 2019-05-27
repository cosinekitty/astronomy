/*
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
const xml2js = require('xml2js');

class Parm {
    constructor(name, type) {
        this.name = name;
        this.type = type;
    }
}

class FuncInfo {
    constructor(name, type, argtext, parmlist) {
        this.name = name;
        this.type = type;
        this.argtext = argtext;
        this.parmlist = parmlist;
    }
}

class Transformer {
    constructor(sectlist) {
        for (let sect of sectlist) {
            switch (sect.$.kind) {
            case 'define':
                console.log(`hydrogen: processing ${sect.memberdef.length} defines`);
                break;
            case 'enum':
                console.log(`hydrogen: processing ${sect.memberdef.length} enums`);
                break;
            case 'typedef':
                console.log(`hydrogen: processing ${sect.memberdef.length} typedefs`);
                break;
            case 'func':
                console.log(`hydrogen: processing ${sect.memberdef.length} functions`);
                break;
            default:
                console.log(`hydrogen: ignoring "${sect.$.kind}"`);
                break;
            }
        }
    }

    Render() {
        return '';
    }
}

function run(headerXmlFileName, markdownFileName) {
    const headerXml = fs.readFileSync(headerXmlFileName);
    const parser = new xml2js.Parser();
    parser.parseString(headerXml, function(err, result) {
        const xform = new Transformer(result.doxygen.compounddef[0].sectiondef);
        const markdown = xform.Render();
        fs.writeFileSync(markdownFileName, markdown);
    });
}

if (process.argv.length === 3) {
    run(process.argv[2], 'README.md');
    process.exit(0);
} else {
    console.log('USAGE: node hydrogen.js header.xml');
    process.exit(1);
}
