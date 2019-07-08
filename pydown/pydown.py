#!/usr/bin/env python3
import sys
import os
import importlib
import inspect

def PrintUsage():
    print("""
USAGE:  pydown.py infile.py outfile.md
""")
    return 1

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

def Markdown(module):
    funclist = []
    classlist = []
    for name, obj in inspect.getmembers(module):
        if not name.startswith('_'):
            if inspect.isfunction(obj):
                funclist.append(obj)
            elif inspect.isclass(obj):
                classlist.append(obj)
            elif inspect.ismodule(obj):
                pass # ignore other modules pulled in
            else:
                print('other', name)
    return ''

def main():
    if len(sys.argv) != 3:
        return PrintUsage()
    inPythonFileName = sys.argv[1]
    outMarkdownFileName = sys.argv[2]
    # Delete output file before we begin.
    # That way, if anything goes wrong, it won't exist,
    # and thus the error becomes conspicuous to scripts/tools.
    if os.access(outMarkdownFileName, os.F_OK):
        os.remove(outMarkdownFileName)
    module = LoadModule(inPythonFileName)
    md = Markdown(module)
    with open(outMarkdownFileName, 'wt') as outfile:
        outfile.write(md)
    return 0

if __name__ == '__main__':
    sys.exit(main())
