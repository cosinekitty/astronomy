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
        for (let item of data.data) {
            listing += `${item.year}-${item.month}-${item.day}T${item.time}Z ${item.phenom}\n`;
        }
    }
}

fs.writeFileSync('seasons.txt', listing);
process.exit(0);
