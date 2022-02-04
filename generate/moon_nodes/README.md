# Test data for moon node calculations

The original test data comes from Fred Espenak's
[Node Passages of the Moon Table](http://astropixels.com/ephemeris/moon/moonnodes2001.html).
I saved the HTML of that page and hand-edited it to create the file `espenak_nodes.txt`.
Then I wrote the script `parse_moon_nodes.py` to transform the data into something
easier for my test programs to process. The result is `moon_nodes.txt`.
