#!/usr/bin/env python3
import sys
import math
import re
sys.path.append('../source/python')
import astronomy

#-----------------------------------------------------------------------------------------------------------

def Test_AstroTime():
    expected_ut = 6910.270978506945
    expected_tt = 6910.271779431480
    time = astronomy.Time.Make(2018, 12, 2, 18, 30, 12.543)
    diff = time.ut - expected_ut
    if abs(diff) > 1.0e-12:
        print('Test_AstroTime: excessive UT error {}'.format(diff))
        sys.exit(1)
    diff = time.tt - expected_tt
    if abs(diff) > 1.0e-12:
        print('Test_AstroTime: excessive TT error {}'.format(diff))
        sys.exit(1)
    s = str(time.Utc())
    if s != '2018-12-02 18:30:12.543000':
        print('Test_AstroTime: Utc() returned incorrect string "{}"'.format(s))
        sys.exit(1)
    time = astronomy.Time.Make(2018, 12, 31, 23, 59, 59.9994)
    s = str(time)
    if s != '2018-12-31T23:59:59.999Z':
        print('Test_AstroTime: expected 2018-12-31T23:59:59.999Z but found {}'.format(s))
        sys.exit(1)
    time = astronomy.Time.Make(2018, 12, 31, 23, 59, 59.9995)
    s = str(time)
    if s != '2019-01-01T00:00:00.000Z':
        print('Test_AstroTime: expected 2019-01-01T00:00:00.000Z but found {}'.format(s))
        sys.exit(1)
    print('Current time =', astronomy.Time.Now())

#-----------------------------------------------------------------------------------------------------------

def Test_GeoMoon():
    time = astronomy.Time.Make(2019, 6, 24, 15, 45, 37)
    vec = astronomy.GeoMoon(time)
    print('Test_GeoMoon: vec = {:0.16f}, {:0.16f}, {:0.16f}'.format(vec.x, vec.y, vec.z))
    # Correct values obtained from C version of GeoMoon calculation
    cx, cy, cz = 0.002674036155459549, -0.0001531716308218381, -0.0003150201604895409
    dx, dy, dz = vec.x - cx, vec.y - cy, vec.z - cz
    diff = math.sqrt(dx*dx + dy*dy + dz*dz)
    print('Test_GeoMoon: diff = {}'.format(diff))    
    if diff > 4.34e-19:
        print('Test_GeoMoon: EXCESSIVE ERROR')
        sys.exit(1)

#-----------------------------------------------------------------------------------------------------------

def Test_AstroCheck(printflag):
    time = astronomy.Time.Make(1700, 1, 1, 0, 0, 0)
    stop = astronomy.Time.Make(2200, 1, 1, 0, 0, 0)
    observer = astronomy.Observer(29, -81, 10)
    if printflag:
        print('o {:0.6f} {:0.6f} {:0.6f}'.format(observer.latitude, observer.longitude, observer.height))
    dt = 10 + math.pi/100
    bodylist = [
        astronomy.BODY_SUN, astronomy.BODY_MOON, astronomy.BODY_MERCURY, astronomy.BODY_VENUS, 
        astronomy.BODY_EARTH, astronomy.BODY_MARS, astronomy.BODY_JUPITER, astronomy.BODY_SATURN, 
        astronomy.BODY_URANUS, astronomy.BODY_NEPTUNE, astronomy.BODY_PLUTO
    ]

    while time.tt < stop.tt:
        for body in bodylist:
            name = astronomy.BodyName[body]
            if body != astronomy.BODY_MOON:
                pos = astronomy.HelioVector(body, time)                
                if printflag:
                    print('v {} {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(name, pos.t.tt, pos.x, pos.y, pos.z))
                if body != astronomy.BODY_EARTH:
                    j2000 = astronomy.Equator(body, time, observer, False, False)
                    ofdate = astronomy.Equator(body, time, observer, True, True)
                    hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.REFRACTION_NONE)
                    if printflag:
                        print('s {} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(name, time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude))
        pos = astronomy.GeoMoon(time)
        if printflag:
            print('v GM {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(pos.t.tt, pos.x, pos.y, pos.z))
        j2000 = astronomy.Equator(astronomy.BODY_MOON, time, observer, False, False)
        ofdate = astronomy.Equator(astronomy.BODY_MOON, time, observer, True, True)
        hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.REFRACTION_NONE)
        if printflag:
            print('s GM {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude))
        time = time.AddDays(dt)

#-----------------------------------------------------------------------------------------------------------

def Test_Seasons(filename):
    with open(filename, 'rt') as infile:
        lnum = 0
        current_year = 0
        mar_count = sep_count = jun_count = dec_count = 0
        max_minutes = 0.0
        for line in infile:
            lnum += 1
            line = line.strip()
            m = re.match(r'^(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([A-Za-z]+)$', line)
            if not m:
                print('Test_Seasons: Invalid data on line {} of file {}'.format(lnum, filename))
                return 1
            year = int(m.group(1))
            month = int(m.group(2))
            day = int(m.group(3))
            hour = int(m.group(4))
            minute = int(m.group(5))
            name = m.group(6)
            if year != current_year:
                current_year = current_year
                seasons = astronomy.Seasons(year)
            correct_time = astronomy.Time.Make(year, month, day, hour, minute, 0)
            if name == 'Equinox':
                if month == 3:
                    calc_time = seasons.mar_equinox
                    mar_count += 1
                elif month == 9:
                    calc_time = seasons.sep_equinox
                    sep_count += 1
                else:
                    print('Test_Seasons: {} line {}: Invalid equinox date in test data'.format(filename, lnum))
                    return 1
            elif name == 'Solstice':
                if month == 6:
                    calc_time = seasons.jun_solstice
                    jun_count += 1
                elif month == 12:
                    calc_time = seasons.dec_solstice
                    dec_count += 1
                else:
                    print('Test_Seasons: {} line {}: Invalid solstice date in test data'.format(filename, lnum))
                    return 1
            elif name == 'Aphelion':
                continue # not yet calculated
            elif name == 'Perihelion':
                continue # not yet calculated
            else:
                print('Test_Seasons: {} line {}: unknown event type {}'.format(filename, lnum, name))
                return 1

            # Verify that the calculated time matches the correct time for this event.
            diff_minutes = (24.0 * 60.0) * abs(calc_time.tt - correct_time.tt)
            if diff_minutes > max_minutes:
                max_minutes = diff_minutes
            
            if diff_minutes > 1.7:
                print('Test_Seasons: {} line {}: excessive error ({}): {} minutes.'.format(filename, lnum, name, diff_minutes))
                return 1
    print('Test_Seasons: verified {} lines from file {} : max error minutes = {:0.3f}'.format(lnum, filename, max_minutes))
    print('Test_Seasons: Event counts: mar={}, jun={}, sep={}, dec={}'.format(mar_count, jun_count, sep_count, dec_count))
    return 0

#-----------------------------------------------------------------------------------------------------------

def Test_MoonPhase(filename):
    threshold_seconds = 120.0       # max tolerable prediction error in seconds
    max_arcmin = 0.0
    maxdiff = 0.0
    quarter_count = 0

    with open(filename, 'rt') as infile:
        lnum = 0
        prev_year = 0
        for line in infile:
            lnum += 1
            line = line.strip()
            m = re.match(r'^([0-3]) (\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+\.\d+)Z$', line)
            if not m:
                print('Test_MoonPhase: invalid data format in {} line {}'.format(filename, lnum))
                return 1

            quarter = int(m.group(1))
            year = int(m.group(2))
            month = int(m.group(3))
            day = int(m.group(4))
            hour = int(m.group(5))
            minute = int(m.group(6))
            second = float(m.group(7))

            expected_elong = 90.0 * quarter
            expected_time = astronomy.Time.Make(year, month, day, hour, minute, second)
            angle = astronomy.MoonPhase(expected_time)
            degree_error = abs(angle - expected_elong)
            if degree_error > 180.0:
                degree_error = 360.0 - degree_error
            arcmin = 60.0 * degree_error
            if arcmin > 1.0:
                print('Test_MoonPhase({} line {}): EXCESSIVE ANGULAR ERROR: {} arcmin'.format(filename, lnum, arcmin))
                return 1
            max_arcmin = max(max_arcmin, arcmin)

            if year != prev_year:
                prev_year = year
                # The test data contains a single year's worth of data for every 10 years.
                # Every time we see the year value change, it breaks continuity of the phases.
                # Start the search over again.
                start_time = astronomy.Time.Make(year, 1, 1, 0, 0, 0.0)
                mq = astronomy.SearchMoonQuarter(start_time)
            else:
                # Yet another lunar quarter in the same year.
                expected_quarter = (1 + mq.quarter) % 4
                mq = astronomy.NextMoonQuarter(mq)
                # Expect the next consecutive quarter.
                if expected_quarter != mq.quarter:
                    print('Test_MoonPhase({} line {}): SearchMoonQuarter returned quarter {}, but expected {}.'.format(filename, lnum, mq.quarter, expected_quarter))
                    return 1

            quarter_count += 1

            # Make sure the time matches what we expect.
            diff_seconds = abs(mq.time.tt - expected_time.tt) * (24.0 * 3600.0)
            if diff_seconds > threshold_seconds:
                print('Test_MoonPhase({} line {}): excessive time error {:0.3f} seconds.'.format(filename, lnum, diff_seconds))
                return 1

            maxdiff = max(maxdiff, diff_seconds)

    print('Test_MoonPhase: passed {} lines for file {} : max_arcmin = {:0.6f}, maxdiff = {:0.3f} seconds, {} quarters.'
        .format(lnum, filename, max_arcmin, maxdiff, quarter_count))
    return 0

#-----------------------------------------------------------------------------------------------------------

if len(sys.argv) == 2:
    if sys.argv[1] == 'time':
        Test_AstroTime()
        sys.exit(0)

    if sys.argv[1] == 'moon':
        Test_GeoMoon()
        sys.exit(0)

    if sys.argv[1] == 'astro_check' or sys.argv[1] == 'astro_profile':
        Test_AstroCheck(sys.argv[1] == 'astro_check')
        sys.exit(0)

if len(sys.argv) == 3:
    if sys.argv[1] == 'seasons':
        sys.exit(Test_Seasons(sys.argv[2]))

    if sys.argv[1] == 'moonphase':
        sys.exit(Test_MoonPhase(sys.argv[2]))

print('test.py: Invalid command line arguments.')
sys.exit(1)
