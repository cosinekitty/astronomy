/*
    parse_moon_phases.js
*/
const fs = require('fs');

listing = '';

for (let filename of fs.readdirSync('.')) {
    let m = filename.match(/^(\d{4})\.json$/);
    if (m) {
        let year = parseInt(m[1]);
        //console.log(year);
        let text = fs.readFileSync(filename, {encoding:'utf8'});
        let data = JSON.parse(text);
        if (data.error !== false) {
            console.log(`ERROR: found error flag (or flag is missing) in file ${filename}`);
            process.exit(1);
        }
        if (data.year !== year) {
            console.log(`ERROR: data.year=${data.year} in file ${filename}`);
            process.exit(1);
        }
        if (!(data.phasedata instanceof Array)) {
            console.log(`ERROR: Missing phasedata array in ${filename}`);
            process.exit(1);
        }
        for (let item of data.phasedata) {
            let when = new Date(`${item.date} ${item.time}Z`);
            let phase = ['New Moon', 'First Quarter', 'Full Moon', 'Last Quarter'].indexOf(item.phase);
            if (phase < 0) {
                console.log(`ERROR: unknown lunar phase name "${item.phase}" in file ${filename}`);
                process.exit(1);
            }
            listing += `${phase} ${when.toISOString()}\n`;
        }
    }
}

fs.writeFileSync('moonphases.txt', listing);
process.exit(0);
