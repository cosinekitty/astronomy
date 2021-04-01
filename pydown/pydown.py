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
    # [`astro_time_t`](#astro_time_t)
    return '[`{0}`](#{0})'.format(name)

def FixText(s):
    # Expand "#Body" to "[`Body`](#Body)".
    return re.sub(r'#([A-Z][A-Za-z0-9_]*)', r'[`\1`](#\1)', s)

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
            md += '### Returns'
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
    text = str(sig)
    text = HtmlEscape(text)
    return text

def MdFunction(func, parent=None):
    md = ''
    doc = inspect.getdoc(func)
    if doc:
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
        md += '\n'
        md += '---\n'
        md += '\n'
        md += '<a name="{}"></a>\n'.format(c.__name__)
        md += '### class ' + c.__name__ + '\n'
        info = DocInfo(doc)
        md += info.Markdown()
        md += '\n'

        firstMemberFunc = True
        for name, obj in inspect.getmembers(c):
            if not name.startswith('_'):
                if firstMemberFunc:
                    firstMemberFunc = False
                    md += '#### member functions\n\n'
                md += MdFunction(obj, parent=c)
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

def Markdown(module, const_md):
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
                print('pydown.py WARNING: ignoring', name)

    md += '---\n'
    md += '\n'
    md += '<a name="constants"></a>\n'
    md += '## Constants\n'
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


def ParseConstants(inPythonFileName):
    md = ''
    with open(inPythonFileName) as infile:
        for line in infile:
            parts = line.split('#<const>')
            if len(parts) == 2:
                code = parts[0].strip()
                doc = parts[1].strip()
                tokens = code.split()
                if len(tokens) >= 3 and tokens[1] == '=':
                    symbol = tokens[0]
                    md += '\n---\n\n'
                    md += '<a name="{}"></a>\n'.format(symbol)
                    md += '### `{}`\n\n'.format(code)
                    md += '**{}**\n\n'.format(doc)
    return md


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
    const_md = ParseConstants(inPythonFileName)
    md = Markdown(module, const_md)

    with open(outMarkdownFileName, 'wt') as outfile:
        outfile.write(prefix)
        outfile.write(md)

    return 0

if __name__ == '__main__':
    sys.exit(main())
