#!/usr/bin/env python3
import os
import re
import shutil

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
    docRootDir = '../../source/kotlin/build/dokka/gfm/astronomy/io.github.cosinekitty.astronomy'
    rootFileName = os.path.join(docRootDir, 'index.md')
    outDir = '../../source/kotlin'
    docDir = os.path.join(outDir, 'doc')
    outFileName = os.path.join(outDir, 'README.md')
    with open(outFileName,'wt') as outfile:
        # Copy the entire prefix file.
        with open(prefixFileName, 'rt') as infile:
            outfile.write(infile.read())

        # Slurp in the interesting parts of the documentation root.
        with open(rootFileName,'rt') as infile:
            # Skip unwanted junk at the top.
            # Start paying attention at the package declaration.
            skipJunk = True
            for line in infile:
                if skipJunk:
                    if line.startswith('# Package '):
                        skipJunk = False
                else:
                    fix = LinkRegex.sub(r'[\1](doc/\2)', line)
                    outfile.write(fix)

    # Nuke the existing 'doc' directory.
    shutil.rmtree(docDir, ignore_errors = True)

    # Recursively copy files from the document root to the target directory.
    shutil.copytree(docRootDir, docDir)
