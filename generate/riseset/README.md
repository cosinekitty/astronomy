# Moon and Sun rise/set test data

Data originally obtained from:
https://aa.usno.navy.mil/data/docs/RS_OneYear.php

That link is now broken. The calculator moved to:
https://aa.usno.navy.mil/data/RS_OneYear

To add more test cases:

1. Visit the second link above. Enter year, latitude, longitude.
2. Generate a year's worth of rise/set data.
3. View page source.
4. Copy and paste just the preformatted calendar section.
5. Save as an html file.
6. Run the script ./parse_riseset.py
7. Verify that the new data has been inserted into riseset.txt.
