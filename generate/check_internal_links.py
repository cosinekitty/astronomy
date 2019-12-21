#!/usr/bin/env python3
import sys
import re

def FindBrokenLinks(inFileName):
    with open(inFileName, 'rt') as infile:
        text = infile.read()
    # Search for all link names, of the form: (#Some.Name)
    linkSet = set(m.group(1) for m in re.finditer(r'\(#([A-Za-z0-9_\.]+)\)', text))
    # Search for all anchor names, of the form: <a name="Some.Name">
    anchorSet = set(m.group(1) for m in re.finditer(r'<\s*a\s+name\s*=\s*"([A-Za-z0-9_\.]+)"\s*>', text))
    # Find all link names for which there is no matching anchor name.
    return sorted(linkSet - anchorSet)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('USAGE: check_internal_links.py infile.md')
        sys.exit(1)
    filename = sys.argv[1]
    badLinks = FindBrokenLinks(filename)
    if len(badLinks) > 0:
        print('ERROR(check_internal_links.py): The following {} links are bad in file {}'.format(len(badLinks), filename))
        for name in badLinks:
            print(name)
        sys.exit(1)
    print('Links are OK in: {}'.format(filename))
    sys.exit(0)
