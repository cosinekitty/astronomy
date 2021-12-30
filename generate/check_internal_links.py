#!/usr/bin/env python3
import sys
import re
import os

def FindBrokenLinks(text):
    # Search for all link names, of the form: (#Some.Name)
    linkSet = set(m.group(1) for m in re.finditer(r'\(#([A-Za-z0-9_\.]+)\)', text))
    # Search for all anchor names, of the form: <a name="Some.Name">
    anchorSet = set(m.group(1) for m in re.finditer(r'<\s*a\s+name\s*=\s*"([A-Za-z0-9_\.]+)"\s*>', text))
    # Find all link names for which there is no matching anchor name.
    return sorted(linkSet - anchorSet)

def FindBogusLinks(text):
    # Search for bogus links of the form [Symbol](Symbol).
    bogusSet = set()
    for m in re.finditer(r'\[([A-Za-z0-9_\.]+)\]\(([A-Za-z0-9_\.]+)\)', text):
        if m.group(1) == m.group(2):
            bogusSet.add(m.group(1))
    return bogusSet

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('USAGE: check_internal_links.py infile.md')
        sys.exit(1)

    rc = 0
    filename = os.path.realpath(sys.argv[1])
    with open(filename, 'rt') as infile:
        text = infile.read()

    badLinks = FindBrokenLinks(text)
    if len(badLinks) > 0:
        rc = 1
        print('ERROR(check_internal_links.py): The following {} links are bad in file {}'.format(len(badLinks), filename))
        for name in badLinks:
            print(name)

    bogusLinks = FindBogusLinks(text)
    if len(bogusLinks) > 0:
        rc = 1
        print('ERROR(check_internal_links.py): The following {} links are of the form "(Symbol)[Symbol]" in file {}'.format(len(bogusLinks), filename))
        for name in bogusLinks:
            print(name)

    if rc == 0:
        print('Links are OK in: {}'.format(filename))

    sys.exit(rc)
