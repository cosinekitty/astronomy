#!/usr/bin/env python3
import sys
import os
import re
import importlib
import inspect
import enum

def PrintUsage():
    print("""
USAGE:  pydown.py prefix.md infile.py outfile.md
""")
    return 1

def Fail(message):
    print('FATAL(pydown):', message)
    sys.exit(1)

def LoadModule(inPythonFileName):
    dir = os.path.dirname(inPythonFileName)
    if not dir:
        dir = os.getcwd()
    sys.path.append(dir)
    modname = os.path.basename(inPythonFileName)
    if modname.endswith('.py'):
        modname = modname[:-3]  # chop off '.py'
    module = importlib.import_module(modname)
    return module

def HtmlEscape(text):
    text = text.replace('&', '&amp;')
    text = text.replace('<', '&lt;')
    text = text.replace('>', '&gt;')
    return text

def SymbolLink(name):
    # Special case: Search() and related functions have return type = "Time or `None`"
    m = re.match(r'^\s*([a-zA-Z0-9_]+)\s+or\s+`None`\s*$', name)
    if m:
        return SymbolLink(m.group(1)) + ' or `None`'

    if 'a' <= name[0] <= 'z':
        # Assume built-in Python identifier, so do not link
        return '`{0}`'.format(name)

    # Other links look like `StateVector[]`. We need to link to StateVector, but exclude the [].
    m = re.match(r'^\s*([a-zA-Z0-9_]+(\.[a-zA-Z0-9_]+)*)([^a-zA-Z0-9_\s]+)', name)
    if m:
        return SymbolLink(m.group(1)) + '`' + m.group(3) + '`'

    # [`astro_time_t`](#astro_time_t)
    return '[`{0}`](#{0})'.format(name)

def FixText(s):
    # Expand "#Body" to "[`Body`](#Body)".
    # Tricky: also need to find "#GravitySimulator.Update",
    # but NOT "#GravitySimulator. Blah blah...".
    return re.sub(r'#([A-Z][A-Za-z0-9_]*(\.[A-Z][A-Za-z0-9_]*)*)', r'[`\1`](#\1)', s)

class ParmInfo:
    def __init__(self, name, type):
        self.name = name
        self.type = type
        self.description = ''

    def AppendDescriptionLine(self, line):
        self.description += line.strip() + ' '

class DocInfo:
    def __init__(self, doc):
        self.description = ''
        self.parameters = []
        self.attributes = []
        self.enumValues = []
        self.returnType = None
        self.returns = ''

        lines = doc.split('\n')

        # First line is boldfaced if followed by blank line.
        if len(lines) >= 2 and lines[0].strip() != '' and lines[1].strip() == '':
            self.summary = lines[0]
            lines = lines[2:]
        else:
            self.summary = ''

        currentAttr = None
        currentParm = None
        mode = ''
        for line in lines:
            if re.match(r'^\-+$', line):
                continue
            if line in ['Parameters', 'Returns', 'Example', 'Examples', 'Attributes', 'Values']:
                mode = line
                continue
            if line.strip() == '':
                if mode == 'code':
                    self.description += '```\n'
                mode = ''
                continue
            if mode == 'Parameters':
                currentParm = self.ProcessParmAttrLine(line, currentParm, self.parameters)
            elif mode == 'Attributes':
                currentAttr = self.ProcessParmAttrLine(line, currentAttr, self.attributes)
            elif mode == 'Returns':
                if line.startswith(' '):
                    self.returns += line.strip() + '\n'
                else:
                    self.returnType = line.strip()
            elif mode == 'Example' or mode == 'Examples':
                pass
            elif mode == 'Values':
                self.ProcessEnumValue(line)
            elif mode == '':
                if re.match(r'^\s*>>>', line):
                    mode = 'code'
                    self.description += '```\n'
                self.description += line + '\n'
            elif mode == 'code':
                self.description += line + '\n'
            else:
                raise Exception('Unknown mode = "{}"'.format(mode))
        if mode == 'code':
            self.description += '```\n'

    def ProcessEnumValue(self, line):
        m = re.match(r'^\s*([A-Za-z][A-Za-z0-9_]+)\s*:\s*(.*)$', line)
        if not m:
            raise Exception('Invalid enum documentation: "{}"'.format(line))
        pair = (m.group(1), m.group(2).strip())
        self.enumValues.append(pair)

    def ProcessParmAttrLine(self, line, item, itemlist):
        if line.startswith(' '):
            # The first line of description, or another line of description.
            item.AppendDescriptionLine(line)
        else:
            # name : type
            token = line.split(':')
            if len(token) != 2:
                raise Exception('Expected name:type but found: "{}"'.format(line))
            item = ParmInfo(token[0].strip(), token[1].strip())
            itemlist.append(item)
        return item

    def Table(self, itemlist, tag):
        md = ''
        if itemlist:
            md += '| Type | {} | Description |\n'.format(tag)
            md += '| --- | --- | --- |\n'
            for p in itemlist:
                if not p.type:
                    raise Exception('Symbol "{}" has missing type declaration.'.format(p.name))
                md += '| {} | {} | {} |\n'.format(SymbolLink(p.type), '`' + p.name + '`', FixText(p.description.strip()))
            md += '\n'
        return md

    def EnumTable(self):
        md = ''
        if self.enumValues:
            md += '| Value | Description |\n'
            md += '| --- | --- |\n'
            for (name, desc) in self.enumValues:
                md += '| {} | {} |\n'.format('`' + name + '`', desc)
        return md

    def Markdown(self):
        md = '\n'
        if self.summary:
            md += '**' + FixText(self.summary) + '**\n\n'
        if self.description:
            md += FixText(self.description) + '\n\n'
        md += self.Table(self.parameters, 'Parameter')
        md += self.Table(self.attributes, 'Attribute')
        md += self.EnumTable()
        if self.returns or self.returnType:
            md += '\n**Returns**'
            if self.returnType:
                md += ': ' + SymbolLink(self.returnType)
            md += '\n'
            md += self.returns + '\n'
            md += '\n'
        md += '\n'
        return md

    def VerifyEnum(self, members):
        defs = set(name for (name, _) in self.enumValues)
        if defs != members:
            print('Actual enums: [' + ', '.join(members) + ']')
            print('Documented enums: [' + ', '.join(defs) + ']')
            raise Exception('Documented enums do not match actual enums.')

def MdSignature(sig):
    # Convert the type signature from inspect.signature() into a string.
    text = str(sig)
    # Now that we have type hints, we get inconsistent return values from
    # inspect.signature() for functions that return either a given type value or None.
    # On older Pythons we see "Optional[astronomy.Time]".
    # On newer Pythons we see "Union[astronomy.Time, NoneType]".
    # This causes unit test failures on GitHub Actions when I check in changes.
    # I prefer the older syntax because it's clearer what my intention is.
    text = re.sub(r'Union\[([A-Za-z_][A-Za-z_0-9\.]*),\s*NoneType\]', r'Optional[\1]', text)
    # Escape square brackets as in Optional[X] to Optional\[X\] for Markdown safety.
    text = text.replace('[', r'\[').replace(']', r'\]')
    # Convert 'astronomy.X' to a link to the type X.
    text = re.sub(r'\bastronomy\.([A-Za-z_][A-Za-z_0-9]*)', lambda m : SymbolLink(m.group(1)), text)
    # Replace quoted forward declarations 'X' with links to the type X.
    text = re.sub(r"'([A-Za-z_][A-Za-z_0-9]*)'", lambda m : SymbolLink(m.group(1)), text)
    # Escape characters as needed for Markdown/HTML.
    text = HtmlEscape(text)
    # Replace clumsy '->' symbols with HTML right-arrow character.
    text = text.replace('-&gt;', '&#8594;')
    return text

def MdFunction(func, parent=None):
    md = ''
    doc = inspect.getdoc(func)
    if doc:
        if doc.startswith('Initialize self.'):
            # Special case: skip trivial constructors that have no documentation.
            return ''
        sig = inspect.signature(func)
        md += '\n'
        if parent:
            name = parent.__name__ + '.' + func.__name__
        else:
            name = func.__name__
            md += '---\n'
            md += '\n'
        md += '<a name="{}"></a>\n'.format(name)
        md += '### ' + name + MdSignature(sig) + '\n'
        info = DocInfo(doc)
        md += info.Markdown()
        md += '\n'
    else:
        Fail('No documentation for function ' + func.__name__)
    return md

def MdClass(c):
    md = ''
    doc = inspect.getdoc(c)
    if doc:
        # Do not document the placeholder type `Any`
        if c.__name__ == 'Any':
            return ''
        md += '\n'
        md += '---\n'
        md += '\n'
        md += '<a name="{}"></a>\n'.format(c.__name__)
        md += '### class ' + c.__name__ + '\n'
        info = DocInfo(doc)
        md += info.Markdown()
        md += '\n'
        func_md = ''
        for name, obj in inspect.getmembers(c):
            if name == '__init__':
                func_md += MdFunction(obj, parent=c)
        for name, obj in inspect.getmembers(c):
            if not name.startswith('_'):
                func_md += MdFunction(obj, parent=c)
        if func_md:
            md += '#### member functions\n\n' + func_md
    else:
        Fail('No documentation for class ' + c.__name__)
    return md

def MdEnumType(c):
    md = ''
    doc = inspect.getdoc(c)
    if doc:
        md += '\n'
        md += '---\n'
        md += '\n'
        md += '<a name="{}"></a>\n'.format(c.__name__)
        md += '### enum ' + c.__name__ + '\n'
        info = DocInfo(doc)
        info.VerifyEnum(set(c.__members__))
        md += info.Markdown()
        md += '\n'
    else:
        Fail('No documentation for enumeration class ' + c.__name__)
    return md

def MdErrType(c):
    md = ''
    doc = inspect.getdoc(c)
    if doc:
        md += '\n'
        md += '---\n'
        md += '\n'
        md += '<a name="{}"></a>\n'.format(c.__name__)
        md += '### ' + c.__name__ + '\n'
        info = DocInfo(doc)
        md += info.Markdown()
        md += '\n'
    else:
        Fail('No documentation for exception class ' + c.__name__)
    return md

def Markdown(module, const_md, const_set):
    md = ''
    funclist = []
    classlist = []
    enumlist = []
    errlist = []
    for name, obj in inspect.getmembers(module):
        if not name.startswith('_'):
            if inspect.isfunction(obj):
                funclist.append(obj)
            elif inspect.isclass(obj):
                if issubclass(obj, enum.Enum):
                    enumlist.append(obj)
                elif issubclass(obj, Exception):
                    errlist.append(obj)
                else:
                    classlist.append(obj)
            elif inspect.ismodule(obj):
                pass # ignore other modules pulled in
            else:
                # Assume this is a global constant. Fail if not documented
                # using my custom "#<const>" documentation text.
                if name not in const_set:
                    Fail('Undocumented symbol: ' + name)


    md += '---\n'
    md += '\n'
    md += '<a name="constants"></a>\n'
    md += '## Constants\n'
    md += 'The following numeric constants are exported by the `astronomy` module.\n'
    md += 'They may be of use for unit conversion.\n'
    md += 'Note: For the other supported programming languages, Astronomy Engine defines\n'
    md += 'helper constants `DEG2RAD` and `RAD2DEG` to convert between angular degrees and radians.\n'
    md += 'However, because Python defines the [angular conversion functions](https://docs.python.org/3/library/math.html#angular-conversion)\n'
    md += '`math.degrees()` and `math.radians()`, they are not needed in the Python version.\n'
    md += '\n'
    md += const_md

    md += '---\n'
    md += '\n'
    md += '<a name="classes"></a>\n'
    md += '## Classes\n'
    md += '\n'
    for c in classlist:
        md += MdClass(c)

    md += '---\n'
    md += '\n'
    md += '<a name="enumerations"></a>\n'
    md += '## Enumerated Types\n'
    md += '\n'
    for c in enumlist:
        md += MdEnumType(c)

    md += '---\n'
    md += '\n'
    md += '<a name="errors"></a>\n'
    md += '## Error Types\n'
    md += '\n'
    for c in errlist:
        md += MdErrType(c)

    md += '---\n'
    md += '\n'
    md += '<a name="functions"></a>\n'
    md += '## Functions\n'
    md += '\n'
    for func in funclist:
        md += MdFunction(func)

    # Remove extraneous blank lines.
    # We never need more than 2 consecutive newline characters.
    md = re.sub('\n{3,}', '\n\n', md)

    return md


def ConstantsMd(inPythonFileName):
    documentedSymbolSet = set()
    md = ''

    clist = []
    with open(inPythonFileName) as infile:
        for line in infile:
            # Consider symbols like Union defined in lines like the following:
            # "from typing import List, Optional, Union, Callable"
            tokens = line.split()
            if len(tokens) > 3 and tokens[0] == 'from' and tokens[2] == 'import':
                tokens = re.findall(r'[A-Za-z_][A-Za-z0-9_]*', line)
                if len(tokens) > 3 and tokens[0] == 'from' and tokens[2] == 'import':
                    for name in tokens[3:]:
                        documentedSymbolSet.add(name)
                    continue

            # Consider symbols like KM_PER_AU defined in lines like the following:
            # "KM_PER_AU = 1.4959787069098932e+8   #<const> The number of kilometers per astronomical unit."
            parts = line.split('#<const>')
            if len(parts) == 2:
                code = parts[0].strip()
                doc = parts[1].strip()
                tokens = code.split()
                if len(tokens) >= 3 and tokens[1] == '=':
                    # Reformat the code to remove extraneous spaces.
                    codeText = ' '.join(tokens).strip()
                    symbol = tokens[0]
                    clist.append((symbol, codeText, doc))
                    documentedSymbolSet.add(symbol)
                    continue

    for (symbol, code, doc) in sorted(clist):
        md += '\n---\n\n'
        md += '<a name="{}"></a>\n'.format(symbol)
        md += '### `{}`\n\n'.format(code)
        md += '**{}**\n\n'.format(doc)

    return md, documentedSymbolSet


def main():
    if len(sys.argv) != 4:
        return PrintUsage()
    prefixFileName = sys.argv[1]
    inPythonFileName = sys.argv[2]
    outMarkdownFileName = sys.argv[3]
    # Delete output file before we begin.
    # That way, if anything goes wrong, it won't exist,
    # and thus the error becomes conspicuous to scripts/tools.
    if os.access(outMarkdownFileName, os.F_OK):
        os.remove(outMarkdownFileName)

    # Load the prefix text.
    with open(prefixFileName, 'rt') as infile:
        prefix = infile.read()

    module = LoadModule(inPythonFileName)
    const_md, const_set = ConstantsMd(inPythonFileName)
    md = Markdown(module, const_md, const_set)

    with open(outMarkdownFileName, 'wt', encoding='utf-8') as outfile:
        outfile.write(prefix)
        outfile.write(md)

    return 0

if __name__ == '__main__':
    sys.exit(main())
