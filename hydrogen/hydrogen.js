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

function GenerateMarkdown(mdFileName, header, source) {
    let text = '# Astronomy Engine\n';
    text += '(Documentation coming soon...)\n';

    fs.writeFileSync(mdFileName, text);
}

function run(headerXmlFileName, sourceXmlFileName, markdownFileName) {
    const headerXml = fs.readFileSync(headerXmlFileName);
    const sourceXml = fs.readFileSync(sourceXmlFileName);
    const headerParser = new xml2js.Parser();
    headerParser.parseString(headerXml, function(err, hresult) {
        const sourceParser = new xml2js.Parser();
        sourceParser.parseString(sourceXml, function(err, cresult) {
            GenerateMarkdown(markdownFileName, hresult, cresult);
        });
    });
}

if (process.argv.length === 4) {
    run(process.argv[2], process.argv[3], 'README.md');
    process.exit(0);
} else {
    console.log('USAGE: node hydrogen.js header.xml source.xml');
    process.exit(1);
}
