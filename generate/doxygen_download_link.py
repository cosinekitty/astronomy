#!/usr/bin/env python3
#
#   Utility for finding the current link for downloading a
#   zip file with Doxygen binaries for Windows.
#   Doxygen conveniently provides us the latest version of Windows
#   binaries for us to download on their website.
#   Not so conveniently, they periodically rename the zip file
#   when they release a new version number, which broke my GitHub Actions
#   scripts. So I created this Python script to scrape the
#   current download link from their website.
#
import sys
import re
import urllib.request

def DoxygenLink():
    url = 'https://www.doxygen.nl/download.html'
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64)',
        'Accept': 'text/html',
    }
    request = urllib.request.Request(url, None, headers)
    with urllib.request.urlopen(request) as response:
        text = response.read().decode('utf-8')
    # Search for the download link.
    mlist = re.findall(r'<a href="(https://www.doxygen.nl/files/doxygen-[0-9\.]+.windows.x64.bin.zip)">', text)
    if len(mlist) != 1:
        return 1
    print(mlist[0])
    return 0

if __name__ == '__main__':
    sys.exit(DoxygenLink())
