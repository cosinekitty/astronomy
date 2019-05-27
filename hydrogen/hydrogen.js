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

function Look(m, key) {
    for (let e of m.$$) {
        if (e['#name'] === key) {
            return e;
        }
    }
    return null;
}

function Find(m, key) {
    let x = Look(m, key);
    if (x === null) {
        throw `Find: could not find key "${key}" in:\n${JSON.stringify(m,null,2)}`;
    }
    return x;
}

function MemberId(m) {
    return m.$.id;
}

class Item {
    constructor(m) {
        this.id = MemberId(m);
        this.name = Find(m, 'name');
        this.detail = Look(m, 'detaileddescription');
        this.brief = Look(m, 'briefdescription');
    }

    static Flat(x) {
        if (typeof x === 'string') {
            return x;
        } 
    
        if (x instanceof Array) {
            let s = '';
            for (let e of x) {
                s += Item.Flat(e);
            }
            return s;
        }     
        
        if (typeof x === 'object') {
            if (x.$$) {
                return Item.Flat(x.$$);
            }
    
            if (x._) {
                return Item.Flat(x._);
            }

            return '';
        }
    
        throw `Item.Flat: don't know how to convert: ${x}`;    
    }
}

class Define extends Item {
    constructor(m) {
        super(m);
        this.init = Find(m, 'initializer');
        console.log(`    #define ${Item.Flat(this.name)} ${Item.Flat(this.init)}  // ${Item.Flat(this.detail)} // ${this.id}`);
    }
}

class EnumInfo extends Item {
    constructor(m) {
        super(m);
        this.enumValueList = []
        for (let e of m.$$) {
            switch (e['#name']) {
            case 'enumvalue':
                this.enumValueList.push(e);
                break;
            }
        }
        console.log(`   enum ${Item.Flat(this.name)} : ${this.enumValueList.length} values.`);
    }
}

class TypeDef extends Item {
    constructor(m) {
        super(m);
        this.definition = Find(m, 'definition');
    }
}

class Parm extends Item {
    constructor(m) {
        super(m);
    }
}

class FuncInfo extends Item {
    constructor(m) {
        super(m);
    }
}

class Transformer {
    constructor(sectlist) {
        for (let sect of sectlist) {
            switch (sect.$.kind) {
            case 'define':
                this.defines = this.ParseDefineList(sect.memberdef);
                break;
            case 'enum':
                this.enums = this.ParseEnumList(sect.memberdef);
                break;
            case 'typedef':
                this.typedefs = this.ParseTypedefList(sect.memberdef);
                break;
            case 'func':
                this.funcs = this.ParseFuncList(sect.memberdef);
                break;
            default:
                console.log(`hydrogen: ignoring "${sect.$.kind}"`);
                break;
            }
        }
    }

    ParseDefineList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} defines`);
        let dlist = [];
        for (let m of mlist) {
            dlist.push(new Define(m));
        }
        return dlist;
    }

    ParseEnumList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} enums`);
        let elist = [];
        for (let m of mlist) {
            elist.push(new EnumInfo(m))
        }
        return elist;
    }

    ParseTypedefList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} typedefs`);
        let tlist = [];
        for (let m of mlist) {
            tlist.push(new TypeDef(m));
        }
        return tlist;
    }

    ParseFuncList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} funcs`);
        let flist = [];
        for (let m of mlist) {
            flist.push(new FuncInfo(m));
        }
        return flist;
    }

    Render() {
        return '';
    }
}

function run(headerXmlFileName, markdownFileName) {
    const headerXml = fs.readFileSync(headerXmlFileName);
    const parser = new xml2js.Parser({
        explicitChildren: true,
        preserveChildrenOrder: true,
        charsAsChildren: true
    });

    parser.parseString(headerXml, function(err, result) {
        const sectlist = result.doxygen.compounddef[0].sectiondef;
        //const dump = JSON.stringify(sectlist, null, 2);
        //fs.writeFileSync('dump.json', dump);
        const xform = new Transformer(sectlist);
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
