'use strict';

const fs = require('fs');
const filepath = '../website/src/assets/documentation.json';

fs.readFile(filepath, (err, data) => {
    if (err) throw err;

    const docs = JSON.stringify(
        JSON.parse(data),
        (k,v) => k === 'path' || k === 'files' ? undefined : v,
        4
    );

    fs.writeFileSync(filepath, docs);
});
