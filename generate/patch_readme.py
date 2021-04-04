#!/usr/bin/env python3
import sys
import os
import re

def PatchReadme(readmeFileName, jsFileName):
    nbytes = os.stat(jsFileName).st_size
    if nbytes >= 100000:
        print('ERROR(patch_readme.py): The size of {} has grown to {} bytes. This is too large!'.format(jsFileName, nbytes))
        return 1

    marker = '<!--MINIFIED_SIZE-->'
    with open(readmeFileName, 'rt') as infile:
        text = infile.read()

    pattern = marker + '[0-9]*'
    repl = marker + str(nbytes)
    updated = re.sub(pattern, repl, text)
    if updated != text:
        with open(readmeFileName, 'wt') as outfile:
            outfile.write(updated)
        print('patch_readme.py: Updated {} with js size reported as {} bytes.'.format(readmeFileName, nbytes))
    else:
        print('patch_readme.py: No changes needed in {}.'.format(readmeFileName))
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('USAGE: patch_readme.py readme_file js_file')
        sys.exit(1)
    readmeFileName = sys.argv[1]
    jsFileName = sys.argv[2]
    sys.exit(PatchReadme(readmeFileName, jsFileName))
