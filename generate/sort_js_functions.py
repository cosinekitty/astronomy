#!/usr/bin/env python3
#
#   sort_js_functions.py - Don Cross <cosinekitty@gmail.com>
#
#   I could not find any easy way to ask jsdoc2md to sort
#   the list of functions alphabetically in the generated markdown
#   output, so I'm doing it myself as a post-processing step.
#
import sys

def SortFunctionNames(intext):
    # Strategy: Split file into blobs separated by "* * *".
    # There is a group of non-function blobs, followed by
    # function blobs, followed again by non-function blobs.
    # Function blobs contain "**Kind**: global function".
    # Generate text with only the function blobs sorted by function name.
    delim = '\n\n* * *\n\n'
    blobs = intext.split(delim)
    foundFunctionSection = False
    blobs_before_funcs = []
    func_blobs = []
    blobs_after_funcs = []
    for b in blobs:
        if '**Kind**: global function' in b:
            foundFunctionSection = True
            func_blobs.append(b)
        else:
            if foundFunctionSection:
                blobs_after_funcs.append(b)
            else:
                blobs_before_funcs.append(b)

    # We can directly sort the markdown output, because the first
    # variable text that appears in each function blob is '<a name="FunctionName"></a>'.
    # If this changes some day, replace func_blobs with a list of tuples: (func_name, blob), then sort the tuple list.
    func_blobs.sort()
    return delim.join(blobs_before_funcs + func_blobs + blobs_after_funcs)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('USAGE: sort_js_functions.py filename.md')
        sys.exit(1)
    filename = sys.argv[1]
    with open(filename, 'rt') as infile:
        intext = infile.read()
    outtext = SortFunctionNames(intext)
    with open(filename, 'wt') as outfile:
        outfile.write(outtext)
    sys.exit(0)
