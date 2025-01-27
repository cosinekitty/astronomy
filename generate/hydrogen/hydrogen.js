/*
    MIT License

    Copyright (c) 2019-2025 Don Cross <cosinekitty@gmail.com>

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

function Find(m, key, altkey) {
    let x = Look(m, key);
    if (x === null) {
        if (altkey) {
            x = Look(m, altkey);
        }
        if (x === null) {
            throw `Find: could not find key "${key}" in:\n${JSON.stringify(m,null,2)}`;
        }
    }
    return x;
}

function MemberId(m) {
    return m.$.id;
}

class Item {
    constructor(m) {
        this.id = MemberId(m);
        this.name = Find(m, 'name', 'compoundname');
        this.detail = Look(m, 'detaileddescription');
        this.brief = Look(m, 'briefdescription');
        this.section = Look(m, 'sectiondef');
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

    static Search(x, key, kind) {
        let queue = [x];
        while (queue.length > 0) {
            let e = queue.shift();
            if (e['#name'] === key) {
                if (!kind || e.$.kind === kind) {
                    return e;
                }
            }

            if (e.$$ instanceof Array) {
                for (let c of e.$$) {
                    queue.push(c);
                }
            }
        }
        return null;
    }

    static Clean(s) {
        if (s && typeof s === 'string') {
            return s.replace(/\s+/g, ' ');
        }
        return '';
    }

    MarkdownPrefix() {
        const name = Item.Flat(this.name);
        let md = `\n\n---\n\n<a name="${name}"></a>\n`;
        return md;
    }

    MdText(x, allow_paragraphs) {
        let md = '';
        let para_sep = allow_paragraphs ? '\n\n' : ' ';
        if (x && x.$$) {
            for (let y of x.$$) {
                let n = y['#name'];
                switch (n) {
                case 'para':
                    md += para_sep + this.MdText(y, allow_paragraphs);
                    break;

                case '__text__':
                    md += y._;
                    break;

                case 'plusmn':
                case 'Delta':
                    md += `&${n};`;
                    break;

                case 'computeroutput':
                    md += '`';
                    md += this.MdText(y, allow_paragraphs);
                    md += '`';
                    break;

                case 'emphasis':
                    md += '*';
                    md += this.MdText(y, allow_paragraphs);
                    md += '*';
                    break;

                case 'ulink':
                    md += '[' + this.MdText(y, allow_paragraphs) + '](' + y.$.url + ')';
                    break;

                case 'ref':
                    md += '[`' + y._ + '`](#' + y._ + ')';
                    break;

                case 'parameterlist':
                    break;  // parameter list rendered separately

                case 'simplesect':
                    if (y.$.kind === 'return') {
                        // documentation about return value of a function.
                        md += para_sep + '**Returns:** ' + this.MdText(y, false) + para_sep;
                    } else {
                        console.log(`MdText: interesting simplesect kind = ${y.$.kind}`);
                    }
                    break;

                case 'itemizedlist':
                    if (!allow_paragraphs)
                        throw 'Cannot include an itemized list inside a table cell.';

                    md += para_sep + this.MdText(y, true);
                    break;

                case 'listitem':
                    if (!allow_paragraphs)
                        throw 'Cannot include a list item inside a table cell.';

                    md += '- ' + this.MdText(y, true).trimLeft() + '\n';
                    break;

                default:
                    console.log(JSON.stringify(y, null, 2));
                    console.log(`MdText: unknown element name: [${n}]`);
                    break;
                }
            }
        }

        if (!allow_paragraphs) {
            md = Item.Clean(md);
        }
        return md;
    }

    MdDescription(brief, detail, allow_paragraphs, require_brief) {
        // I could prevent gh-pages (kramdown) escaping html inside tables:
        // https://stackoverflow.com/questions/47262698/inline-html-is-escaped-by-jekyll
        // But then it doesn't render correctly on github.com!  Ya just can't win.

        let md = '';

        let btext = this.MdText(brief);
        if (btext) {
            md += '**' + btext.trim() + '** ';
            if (allow_paragraphs)
                md += '\n\n';
        } else if (require_brief) {
            throw `Missing brief text for ${require_brief}`;
        }

        let dtext = this.MdText(detail, allow_paragraphs);
        if (dtext) {
            md += dtext;
        }

        return md;
    }

    static MdType(type) {
        let text = Item.Flat(type);
        let md;
        if (text.indexOf('astro_') === 0) {
            // Create a link to our custom type.
            // But watch out for pointer types like 'astro_time_t *'.
            // For those, we need to extract and link to just the type itself.
            const m = /([a-z_]+)(.*)/.exec(text);
            const symbol = m[1];
            const tail = m[2];
            if (tail.length > 0) {
                // Verbose HTML syntax gets the job done best.
                md = `<code><a href="#${symbol}">${symbol}</a>${tail}</code>`;
            } else {
                md = '[`' + symbol + '`](#' + symbol + ')';
            }
        } else {
            // assume built-in type that we can't link to
            md = '`' + text + '`';
        }
        return md;
    }
}

class Define extends Item {
    constructor(m) {
        super(m);
        this.init = Find(m, 'initializer');
    }

    Markdown() {
        let md = this.MarkdownPrefix();
        let name = Item.Flat(this.name);
        md += '### `' + name + '`\n\n';
        md += this.MdDescription(this.brief, this.detail, true, '#define ' + name);
        let defn = '';
        for (let term of this.init['$$']) {
            defn += term._;
        }
        md += "\n\n```C\n#define " + name + "  " + defn + "\n```\n\n";
        return md;
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
    }

    Markdown() {
        let name = Item.Flat(this.name);
        let md = this.MarkdownPrefix();
        md += '### `' + name + '`\n\n';
        md += this.MdDescription(this.brief, this.detail, true, 'enum ' + name);
        if (this.enumValueList instanceof Array && this.enumValueList.length > 0) {
            md += '\n\n| Enum Value | Description |\n';
            md += '| --- | --- |\n';
            for (let e of this.enumValueList) {
                let ename = Item.Flat(Find(e, 'name'));
                let ebrief = Look(e, 'briefdescription');
                let edetail = Look(e, 'detaileddescription');
                let desc = this.MdDescription(ebrief, edetail, false);
                md += '| `' + ename + '` | ' + desc + ' |\n';
            }
            md += '\n';
        }
        return md;
    }
}

class TypeDef extends Item {
    constructor(m) {
        super(m);
        this.definition = Find(m, 'definition');
    }

    Markdown() {
        let name = Item.Flat(this.name);

        let defn = Item.Flat(this.definition);

        // Normalize the definition text by squashing multiple contiguous
        // spaces into a single space.
        // This is to get around an inconsistency in the output generated by
        // doxygen version 1.8.13 (Debian Buster) and 1.9.1 (Debian Bullseye).
        // -`typedef astro_func_result_t(*  astro_search_func_t) (void *context, astro_time_t time);`
        // +`typedef astro_func_result_t(* astro_search_func_t) (void *context, astro_time_t time);`        
        defn = defn.replace(/ {2,}/g, ' ');

        let md = this.MarkdownPrefix();
        md += '### `' + name + '`\n\n';
        md += '`' + defn + ';`\n\n';
        md += this.MdDescription(this.brief, this.detail, true, 'typedef ' + name);
        return md;
    }
}

class FuncInfo extends Item {
    constructor(m) {
        super(m);
        this.paramlist = this.ParamList(m);     // actual declared parameters and types in code (not doxygen comments)
        this.paramdocs = Item.Search(this.detail, 'parameterlist', 'param');        // documentation about parameters in doxygen comments
        this.rettype = Find(m, 'type');
    }

    ParamList(m) {
        var plist = [];
        if (m.$$ instanceof Array) {
            for (let e of m.$$) {
                if (e['#name'] === 'param') {
                    let type = Item.Flat(Find(e, 'type'));
                    if (type !== 'void') {      // special case: func(void) has no parameters.
                        let name = Item.Flat(Find(e, 'declname'));
                        plist.push({type:type, name:name});
                    }
                }
            }
        }
        return plist;
    }

    Markdown() {
        let name = Item.Flat(this.name);
        let md = this.MarkdownPrefix();
        md += '### ' + name;
        md += this.MdParamNameList();
        md += ' &#8658; ' + Item.MdType(this.rettype) + '\n\n';
        md += this.MdDescription(this.brief, this.detail, true, 'function ' + name);

        if (this.paramdocs) {
            md += '\n\n| Type | Parameter | Description |\n';
            md += '| --- | --- | --- |\n';
            for (let p of this.paramdocs.$$) {
                if (p['#name'] === 'parameteritem') {
                    let nameElem = Item.Search(p, 'parametername');
                    let name = nameElem ? Item.Flat(nameElem) : '';
                    let descElem = Item.Search(p, 'parameterdescription');
                    let desc = descElem ? this.MdDescription(null, descElem, false) : '';
                    let type = this.MdTypeForName(name);
                    md += '| ' + type + ' | `' + name + '` | ' + desc + ' | \n';
                }
            }
            md += '\n\n';
        }

        return md;
    }

    MdTypeForName(name) {
        if (name && typeof name === 'string') {
            for (let p of this.paramlist) {
                if (p.name === name) {
                    return Item.MdType(p.type);
                }
            }
            throw `Cannot find parameter called ${name} in function ${Item.Flat(this.name)}`;
        }
        throw `Invalid/missing parameter name in function ${Item.Flat(this.name)}`;
    }

    MdParamNameList() {
        // Generate a comma-delimited list of the parameter names.
        let md = '(';
        if (this.paramlist) {
            let count = 0;
            for (let p of this.paramlist) {
                if (++count > 1)
                    md += ', ';
                md += p.name;
            }
        }
        md += ')';
        return md;
    }
}

class StructInfo extends Item {
    constructor(m) {
        super(m);
    }

    Markdown() {
        let name = Item.Flat(this.name);
        let md = this.MarkdownPrefix();
        md += '### `' + name + '`\n\n';
        md += this.MdDescription(this.brief, this.detail, true);
        md += '\n\n';

        if (!this.section) {
            throw `struct "${name}" has no inner members.`;
        }

        md += '| Type | Member | Description |\n';
        md += '| ---- | ------ | ----------- |\n';
        for (let member of this.section.$$) {
            if (member.$.kind === 'variable' && member.$.prot === 'public') {
                md += this.MdMember(member);
            }
        }
        return md;
    }

    MdMember(member) {
        let name = Find(member, 'name');
        let type = Find(member, 'type');
        let brief = Look(member, 'briefdescription');
        let detail = Look(member, 'detaileddescription');
        let ntext = Item.Flat(name);
        let ttext = Item.MdType(type);
        let md = '| ' + ttext + ' | `' + ntext + '` | ' + this.MdDescription(brief, detail, false) + ' |\n';
        return md;
    }
}

class Transformer {
    constructor() {
        this.structs = [];
    }

    AddStruct(m) {
        this.structs.push(new StructInfo(m));
    }

    AddSectList(sectlist) {
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

    static SortList(list) {
        return list.sort((a, b) => {
            let fa = Item.Flat(a.name);
            let fb = Item.Flat(b.name);
            if (fa < fb)
                return -1;
            if (fa > fb)
                return +1;
            return 0;
        });
    }

    ParseDefineList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} defines`);
        let dlist = [];
        for (let m of mlist) {
            dlist.push(new Define(m));
        }
        return Transformer.SortList(dlist);
    }

    ParseEnumList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} enums`);
        let elist = [];
        for (let m of mlist) {
            elist.push(new EnumInfo(m))
        }
        return Transformer.SortList(elist);
    }

    ParseTypedefList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} typedefs`);
        let tlist = [];
        for (let m of mlist) {
            tlist.push(new TypeDef(m));
        }
        return Transformer.SortList(tlist);
    }

    ParseFuncList(mlist) {
        console.log(`hydrogen: processing ${mlist.length} funcs`);
        let flist = [];
        for (let m of mlist) {
            flist.push(new FuncInfo(m));
        }
        return Transformer.SortList(flist);
    }

    Markdown() {
        let md = '';

        md += '\n<a name="functions"></a>\n';
        md += '## Functions\n\n';
        for (let f of this.funcs) {
            md += f.Markdown();
        }

        md += '\n<a name="constants"></a>\n';
        md += '## Constants\n\n';
        for (let d of this.defines) {
            md += d.Markdown();
        }

        md += '\n<a name="enums"></a>\n';
        md += '## Enumerated Types\n\n';
        for (let e of this.enums) {
            md += e.Markdown();
        }

        md += '\n<a name="structs"></a>\n';
        md += '## Structures\n\n';
        for (let s of this.structs) {
            md += s.Markdown();
        }

        md += '\n<a name="typedefs"></a>\n';
        md += '## Type Definitions\n\n';
        for (let t of this.typedefs) {
            md += t.Markdown();
        }

        return md;
    }
}

function run(inPrefixFileName, inXmlFileName, outMarkdownFileName) {
    const headerXml = fs.readFileSync(inXmlFileName);
    const parser = new xml2js.Parser({
        explicitChildren: true,
        preserveChildrenOrder: true,
        charsAsChildren: true
    });

    parser.parseString(headerXml, function(err, result) {
        const xform = new Transformer();
        for (let compound of result.doxygen.compounddef) {
            if (compound.$.kind === 'struct') {
                xform.AddStruct(compound);
            } else if (compound.$.kind === 'file') {
                if (compound.$.id === 'astronomy_8h') {
                    xform.AddSectList(compound.sectiondef);
                }
            }
        }

        // https://stackoverflow.com/questions/47262698/inline-html-is-escaped-by-jekyll
        let md = fs.readFileSync(inPrefixFileName, {encoding: 'utf8'});
        md += xform.Markdown();
        fs.writeFileSync(outMarkdownFileName, md);
    });
}

if (process.argv.length === 4) {
    run(process.argv[2], process.argv[3], 'README.md');
    process.exit(0);
} else {
    console.log('USAGE: node hydrogen.js prefix.md header.xml');
    process.exit(1);
}
