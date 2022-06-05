#!/usr/bin/env python3
import sys

def Patch(version, filename, prefix, suffix):
    with open(filename, 'rt') as infile:
        text = infile.read()
    pindex = text.find(prefix)
    if pindex < 0:
        print('patch_version_numbers.py: ERROR: Cannot find prefix "{}" in file: {}'.format(prefix, filename))
        return 1
    pindex += len(prefix)
    sindex = text.find(suffix, pindex)
    if sindex < 0:
        print('patch_version_numbers.py: ERROR: Cannot find suffix "{}" in file: {}'.format(suffix, filename))
        return 1
    updated = text[:pindex] + version + text[sindex:]
    if text != updated:
        print('patch_version_numbers.py: UPDATING: {}'.format(filename))
        with open(filename, 'wt') as outfile:
            outfile.write(updated)
    return 0

def PatchVersionNumbers():
    with open('version.txt', 'rt') as infile:
        version = infile.read().strip()
    print('patch_version_numbers.py: Version = {}'.format(version))
    return (
        Patch(version, '../README.md', '2b-v', '-blue') |
        Patch(version, '../source/js/package.json', '"version": "', '"') |
        Patch(version, '../source/csharp/astronomy.csproj', '<PackageVersion>', '</PackageVersion>') |
        Patch(version, '../source/python/setup.py', "version='", "'") |
        Patch(version, '../source/kotlin/build.gradle.kts', 'version = "', '"')
    )

if __name__ == '__main__':
    sys.exit(PatchVersionNumbers())
