#!/usr/bin/env python3
import sys
sys.path.append('../source/python')
import astronomy


def Test_AstroTime():
    expected_ut = 6910.270978506945
    expected_tt = 6910.271779431480
    time = astronomy.MakeTime(2018, 12, 2, 18, 30, 12.543)
    diff = time.ut - expected_ut
    if abs(diff) > 1.0e-12:
        print('Test_AstroTime: excessive UT error {}'.format(diff))
        sys.exit(1)

    diff = time.tt - expected_tt
    if abs(diff) > 1.0e-12:
        print('Test_AstroTime: excessive TT error {}'.format(diff))
        sys.exit(1)


if len(sys.argv) == 2:
    if sys.argv[1] == 'time':
        Test_AstroTime()
        sys.exit(0)

print('test.py: Invalid command line arguments.')
sys.exit(1)
