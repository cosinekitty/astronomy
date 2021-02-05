#!/usr/bin/env python3
import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('USAGE: apply_exports_hack.py js_filename')
        sys.exit(1)
    filename = sys.argv[1]
    with open(filename, 'rt') as infile:
        text = infile.read()

    find = "})(this.Astronomy = {});"
    repl = "})(typeof exports==='undefined' ? (this.Astronomy={}) : exports);"

    if text.find(find) < 0:
        print('ERROR: could not find target text in file: {}'.format(filename))
        sys.exit(1)

    text = text.replace(find, repl)

    with open(filename, 'wt') as outfile:
        outfile.write(text)

    print('apply_exports_hack.py: patched file {}'.format(filename))
    sys.exit(0)
