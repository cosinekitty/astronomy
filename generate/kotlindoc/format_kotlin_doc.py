#!/usr/bin/env python3
#
#   format_kotlin_doc.py  -  Don Cross <cosinekitty@gmail.com>
#
#   A custom utility for post-processing Markdown documentation
#   for the Kotlin version of Astronomy Engine.
#   There are a few things about the `dokkaGfm` tool I don't like.
#   This script is a hack to work around them.
#
import os
import re
import shutil


SymbolsWithUnwantedArgs = [
    'SSB',
    'EMB',
    'Moon',
    'Sun',
    'Mercury',
    'Venus',
    'Earth',
    'Mars',
    'Jupiter',
    'Saturn',
    'Uranus',
    'Neptune',
    'Pluto'
]


def RemoveJvmTags(text):
    fix = text.replace(' | [jvm]', '')
    fix = fix.replace('[jvm]\\\n', '')
    return fix


def RemoveEnumProperties(text):
    # Remove the entire 'Properties' section if it consists of `name` and `ordinal` only.
    # 1. These are unhelpful noise. They exist for all enum classes, and don't tell us anything.
    # 2. There is a bug in Dokka that makes them link to the wrong class.
    #    I reported this: https://github.com/Kotlin/dokka/issues/2464
    # In the meantime, I strip them out entirely.
    # Search for '## Properties'.
    # If all the listed items consist of `name` and `ordinal` only, then
    # chop off the entire section.
    prefix = '## Properties'
    sectionIndex = text.find(prefix)
    if sectionIndex >= 0:
        tail = text[sectionIndex + len(prefix):]
        nameSet = set(re.findall(r'\n\|\s+\[(\S+)\]', tail))
        if nameSet == set(['name', 'ordinal']):
            # Truncate the markdown text right before '## Properties'
            return text[:sectionIndex]
    # Return the entire text, unmodified.
    return text


def RemovePrivateConstructors(text):
    # enum class Body contains a private constructor where members
    # can initialize their internal properties.
    # This is leaking into the dokkaGfm output.
    # I submitted the following issue about this:
    # https://github.com/Kotlin/dokka/issues/2468
    # Example:
    # | [Saturn](-saturn/index.md) | [jvm]<br>[Saturn](-saturn/index.md)(SATURN_GM, 10759.22, VsopModel(vsopLonSaturn, vsopLatSaturn, vsopRadSaturn))<br>The planet Saturn. |
    # should be converted to:
    # | [Saturn](-saturn/index.md) | [jvm]<br>The planet Saturn. |
    fix = text
    for sym in SymbolsWithUnwantedArgs:
        prefix = '<br>[' + sym + ']'
        front = fix.find(prefix)
        if front >= 0:
            # Chop out the entire section that contains the unwanted constructor call.
            back = fix.find('<br>', front + len(prefix))
            if back > front:
                fix = fix[:front] + fix[back:]

        # Another case: the individual enum member has its own index.md file.
        # This has an unwanted constructor call.
        if ('# ' + sym + '\n') in fix:
            prefix = '\n[' + sym + '](index.md)('
            front = fix.find(prefix)
            if front >= 0:
                suffix = '\n'
                back = fix.find(suffix, front + len(prefix))
                if back > front:
                    fix = fix[:front] + fix[back + len(suffix):]

    return fix


def ReverseEnumEntries(text):
    # Enumeration entries are listed in backwards order.
    # I reported this as:
    # https://github.com/Kotlin/dokka/issues/2466
    # Work around this issue for now...
    fix = text
    prefix = '## Entries\n\n| | |\n|---|---|\n'
    middleIndex = fix.find(prefix)
    if middleIndex >= 0:
        middleIndex += len(prefix)
        front = fix[:middleIndex]
        middle = fix[middleIndex:]
        backIndex = middle.find('\n## ')
        if backIndex > 0:
            backIndex += 1
            back = middle[backIndex:]
            middle = middle[:backIndex]
        else:
            back = ''
        rows = middle.strip().split('\n')
        rows.reverse()
        fix = front + '\n'.join(rows) + '\n\n' + back
    return fix


def FixMarkdown(text):
    fix = RemoveEnumProperties(text)
    fix = RemoveJvmTags(fix)
    fix = RemovePrivateConstructors(fix)
    fix = ReverseEnumEntries(fix)
    return fix


# A typical line looks like this:
# | [Direction](-direction/index.md) | [jvm]<br>enum [Direction](-direction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Direction](-direction/index.md)&gt; <br>Selects whether to search for a rising event or a setting event for a celestial body. |
# Find all the links like (-direction/index.md) and prefix them with 'doc/':
# [Direction](doc/-direction/index.md)
# This is probably one of the most horrible regular expressions
# I have ever crafted. One tricky bit is the use of [^\):].
# Apart from looking like a surreal emoji, the ':' prevents
# substitution of external hyperlinks like
# https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html
LinkRegex = re.compile(r'\[([^\]]+)\]\(([^\):]+)\)')     # [Direction](-direction/index.md)

if __name__ == '__main__':
    prefixFileName = 'kotlin_prefix.md'
    inDocDir = '../../source/kotlin/build/dokka/gfm/astronomy/io.github.cosinekitty.astronomy'
    rootFileName = os.path.join(inDocDir, 'index.md')
    outDir = '../../source/kotlin'
    outDocDir = os.path.join(outDir, 'doc')
    outFileName = os.path.join(outDir, 'README.md')
    with open(outFileName,'wt') as outfile:
        # Copy the entire prefix file.
        with open(prefixFileName, 'rt') as infile:
            outfile.write(infile.read())

        # Slurp in the interesting parts of the documentation root.
        with open(rootFileName,'rt') as infile:
            # Skip unwanted junk at the top.
            # Start paying attention after the package declaration.
            skipJunk = True
            for line in infile:
                if skipJunk:
                    if line.startswith('# Package '):
                        skipJunk = False
                else:
                    fix = LinkRegex.sub(r'[\1](doc/\2)', line)
                    fix = RemoveJvmTags(fix)
                    outfile.write(fix)

    # Nuke the existing 'doc' directory.
    shutil.rmtree(outDocDir, ignore_errors = True)

    # Recursively process files, reading from the raw files generated by Dokka,
    # and writing to the target directory to be hosted on GitHub.
    for dirName, childDirs, files in os.walk(inDocDir):
        # Calculate the target relative path.
        targetDir = os.path.join(outDocDir, dirName[1+len(inDocDir):])
        os.makedirs(targetDir, exist_ok = True)
        for fn in files:
            if fn.endswith('.md'):
                inFileName = os.path.join(dirName, fn)
                with open(inFileName, 'rt') as infile:
                    text = infile.read()
                outFileName = os.path.join(targetDir, fn)
                with open(outFileName, 'wt') as outfile:
                    fix = FixMarkdown(text)
                    outfile.write(fix)
