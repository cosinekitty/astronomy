#!/usr/bin/env python3
import sys
import os
import importlib

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

def main():
    if len(sys.argv) != 3:
        return PrintUsage()
    inPythonFileName = sys.argv[1]
    outMarkdownFileName = sys.argv[2]
    module = LoadModule(inPythonFileName)
    return 0

if __name__ == '__main__':
    sys.exit(main())
