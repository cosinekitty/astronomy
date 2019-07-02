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

def TestElongFile(filename, targetRelLon):
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            line = line.strip()
            m = re.match(r'^(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z ([A-Za-z]+)$', line)
            if not m:
                print('TestElongFile({} line {}): invalid data format'.format(filename, lnum))
                return 1
            year = int(m.group(1))
            month = int(m.group(2))
            day = int(m.group(3))
            hour = int(m.group(4))
            minute = int(m.group(5))
            name = m.group(6)
            body = astronomy.BodyCode(name)
            if body < 0:
                print('TestElongFile({} line {}): invalid body name "{}"'.format(filename, lnum, name))
                return 1
            search_time = astronomy.Time.Make(year, 1, 1, 0, 0, 0)
            expected_time = astronomy.Time.Make(year, month, day, hour, minute, 0)
            found_time = astronomy.SearchRelativeLongitude(body, targetRelLon, search_time)
            if found_time is None:
                print('TestElongFile({} line {}): SearchRelativeLongitude failed.'.format(filename, lnum))
                return 1
            diff_minutes = (24.0 * 60.0) * (found_time.tt - expected_time.tt)
            print('TestElongFile: {:<7s} error = {:6.3} minutes'.format(name, diff_minutes))
            if abs(diff_minutes) > 15.0:
                print('TestElongFile({} line {}): EXCESSIVE ERROR.'.format(filename, lnum))
                return 1
    print('TestElongFile: passed {} rows of data'.format(lnum))
    return 0

def TestPlanetLongitudes(body, outFileName, zeroLonEventName):
    startYear = 1700
    stopYear = 2200
    rlon = 0.0
    sum_diff = 0.0
    count = 0
    name = astronomy.BodyName[body]
    with open(outFileName, 'wt') as outfile:
        time = astronomy.Time.Make(startYear, 1, 1, 0, 0, 0)
        stopTime = astronomy.Time.Make(stopYear, 1, 1, 0, 0, 0)
        while time.tt < stopTime.tt:
            count += 1
            event = zeroLonEventName if rlon == 0.0 else 'sup'
            found_time = astronomy.SearchRelativeLongitude(body, rlon, time)
            if found_time is None:
                print('TestPlanetLongitudes({}): SearchRelativeLongitudes failed'.format(name))
                return 1
            if count >= 2:
                # Check for consistent intervals.
                # Mainly I don't want to skip over an event!
                day_diff = found_time.tt - time.tt
                sum_diff += day_diff
                if count == 2:
                    min_diff = max_diff = day_diff
                else:
                    min_diff = min(min_diff, day_diff)
                    max_diff = max(max_diff, day_diff)
            geo = astronomy.GeoVector(body, found_time, True)
            dist = geo.Length()
            outfile.write('e {} {} {:0.16f} {:0.16f}\n'.format(name, event, found_time.tt, dist))
            # Search for the opposite longitude vent next time.
            time = found_time
            rlon = 180.0 - rlon
    if body == astronomy.BODY_MERCURY:
        thresh = 1.65
    elif body == astronomy.BODY_MARS:
        thresh = 1.30
    else:
        thresh = 1.07
    ratio = max_diff / min_diff
    print('TestPlanetLongitudes({:<7s}): {:5d} events, ratio={:5.3f}, file: {}'.format(name, count, ratio, outFileName))
    if ratio > thresh:
        print('TestPlanetLongitudes({}): EXCESSIVE EVENT INTERVAL RATIO'.format(name))
        return 1
    return 0

ElongTestData = [
    # Max elongation data obtained from:
    # http://www.skycaramba.com/greatest_elongations.shtml
    ( astronomy.BODY_MERCURY, "2010-01-17T05:22Z", "2010-01-27T05:22Z", 24.80, 'morning' ),
    ( astronomy.BODY_MERCURY, "2010-05-16T02:15Z", "2010-05-26T02:15Z", 25.10, 'morning' ),
    ( astronomy.BODY_MERCURY, "2010-09-09T17:24Z", "2010-09-19T17:24Z", 17.90, 'morning' ),
    ( astronomy.BODY_MERCURY, "2010-12-30T14:33Z", "2011-01-09T14:33Z", 23.30, 'morning' ),
    ( astronomy.BODY_MERCURY, "2011-04-27T19:03Z", "2011-05-07T19:03Z", 26.60, 'morning' ),
    ( astronomy.BODY_MERCURY, "2011-08-24T05:52Z", "2011-09-03T05:52Z", 18.10, 'morning' ),
    ( astronomy.BODY_MERCURY, "2011-12-13T02:56Z", "2011-12-23T02:56Z", 21.80, 'morning' ),
    ( astronomy.BODY_MERCURY, "2012-04-08T17:22Z", "2012-04-18T17:22Z", 27.50, 'morning' ),
    ( astronomy.BODY_MERCURY, "2012-08-06T12:04Z", "2012-08-16T12:04Z", 18.70, 'morning' ),
    ( astronomy.BODY_MERCURY, "2012-11-24T22:55Z", "2012-12-04T22:55Z", 20.60, 'morning' ),
    ( astronomy.BODY_MERCURY, "2013-03-21T22:02Z", "2013-03-31T22:02Z", 27.80, 'morning' ),
    ( astronomy.BODY_MERCURY, "2013-07-20T08:51Z", "2013-07-30T08:51Z", 19.60, 'morning' ),
    ( astronomy.BODY_MERCURY, "2013-11-08T02:28Z", "2013-11-18T02:28Z", 19.50, 'morning' ),
    ( astronomy.BODY_MERCURY, "2014-03-04T06:38Z", "2014-03-14T06:38Z", 27.60, 'morning' ),
    ( astronomy.BODY_MERCURY, "2014-07-02T18:22Z", "2014-07-12T18:22Z", 20.90, 'morning' ),
    ( astronomy.BODY_MERCURY, "2014-10-22T12:36Z", "2014-11-01T12:36Z", 18.70, 'morning' ),
    ( astronomy.BODY_MERCURY, "2015-02-14T16:20Z", "2015-02-24T16:20Z", 26.70, 'morning' ),
    ( astronomy.BODY_MERCURY, "2015-06-14T17:10Z", "2015-06-24T17:10Z", 22.50, 'morning' ),
    ( astronomy.BODY_MERCURY, "2015-10-06T03:20Z", "2015-10-16T03:20Z", 18.10, 'morning' ),
    ( astronomy.BODY_MERCURY, "2016-01-28T01:22Z", "2016-02-07T01:22Z", 25.60, 'morning' ),
    ( astronomy.BODY_MERCURY, "2016-05-26T08:45Z", "2016-06-05T08:45Z", 24.20, 'morning' ),
    ( astronomy.BODY_MERCURY, "2016-09-18T19:27Z", "2016-09-28T19:27Z", 17.90, 'morning' ),
    ( astronomy.BODY_MERCURY, "2017-01-09T09:42Z", "2017-01-19T09:42Z", 24.10, 'morning' ),
    ( astronomy.BODY_MERCURY, "2017-05-07T23:19Z", "2017-05-17T23:19Z", 25.80, 'morning' ),
    ( astronomy.BODY_MERCURY, "2017-09-02T10:14Z", "2017-09-12T10:14Z", 17.90, 'morning' ),
    ( astronomy.BODY_MERCURY, "2017-12-22T19:48Z", "2018-01-01T19:48Z", 22.70, 'morning' ),
    ( astronomy.BODY_MERCURY, "2018-04-19T18:17Z", "2018-04-29T18:17Z", 27.00, 'morning' ),
    ( astronomy.BODY_MERCURY, "2018-08-16T20:35Z", "2018-08-26T20:35Z", 18.30, 'morning' ),
    ( astronomy.BODY_MERCURY, "2018-12-05T11:34Z", "2018-12-15T11:34Z", 21.30, 'morning' ),
    ( astronomy.BODY_MERCURY, "2019-04-01T19:40Z", "2019-04-11T19:40Z", 27.70, 'morning' ),
    ( astronomy.BODY_MERCURY, "2019-07-30T23:08Z", "2019-08-09T23:08Z", 19.00, 'morning' ),
    ( astronomy.BODY_MERCURY, "2019-11-18T10:31Z", "2019-11-28T10:31Z", 20.10, 'morning' ),
    ( astronomy.BODY_MERCURY, "2010-03-29T23:32Z", "2010-04-08T23:32Z", 19.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2010-07-28T01:03Z", "2010-08-07T01:03Z", 27.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2010-11-21T15:42Z", "2010-12-01T15:42Z", 21.50, 'evening' ),
    ( astronomy.BODY_MERCURY, "2011-03-13T01:07Z", "2011-03-23T01:07Z", 18.60, 'evening' ),
    ( astronomy.BODY_MERCURY, "2011-07-10T04:56Z", "2011-07-20T04:56Z", 26.80, 'evening' ),
    ( astronomy.BODY_MERCURY, "2011-11-04T08:40Z", "2011-11-14T08:40Z", 22.70, 'evening' ),
    ( astronomy.BODY_MERCURY, "2012-02-24T09:39Z", "2012-03-05T09:39Z", 18.20, 'evening' ),
    ( astronomy.BODY_MERCURY, "2012-06-21T02:00Z", "2012-07-01T02:00Z", 25.70, 'evening' ),
    ( astronomy.BODY_MERCURY, "2012-10-16T21:59Z", "2012-10-26T21:59Z", 24.10, 'evening' ),
    ( astronomy.BODY_MERCURY, "2013-02-06T21:24Z", "2013-02-16T21:24Z", 18.10, 'evening' ),
    ( astronomy.BODY_MERCURY, "2013-06-02T16:45Z", "2013-06-12T16:45Z", 24.30, 'evening' ),
    ( astronomy.BODY_MERCURY, "2013-09-29T09:59Z", "2013-10-09T09:59Z", 25.30, 'evening' ),
    ( astronomy.BODY_MERCURY, "2014-01-21T10:00Z", "2014-01-31T10:00Z", 18.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2014-05-15T07:06Z", "2014-05-25T07:06Z", 22.70, 'evening' ),
    ( astronomy.BODY_MERCURY, "2014-09-11T22:20Z", "2014-09-21T22:20Z", 26.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2015-01-04T20:26Z", "2015-01-14T20:26Z", 18.90, 'evening' ),
    ( astronomy.BODY_MERCURY, "2015-04-27T04:46Z", "2015-05-07T04:46Z", 21.20, 'evening' ),
    ( astronomy.BODY_MERCURY, "2015-08-25T10:20Z", "2015-09-04T10:20Z", 27.10, 'evening' ),
    ( astronomy.BODY_MERCURY, "2015-12-19T03:11Z", "2015-12-29T03:11Z", 19.70, 'evening' ),
    ( astronomy.BODY_MERCURY, "2016-04-08T14:00Z", "2016-04-18T14:00Z", 19.90, 'evening' ),
    ( astronomy.BODY_MERCURY, "2016-08-06T21:24Z", "2016-08-16T21:24Z", 27.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2016-12-01T04:36Z", "2016-12-11T04:36Z", 20.80, 'evening' ),
    ( astronomy.BODY_MERCURY, "2017-03-22T10:24Z", "2017-04-01T10:24Z", 19.00, 'evening' ),
    ( astronomy.BODY_MERCURY, "2017-07-20T04:34Z", "2017-07-30T04:34Z", 27.20, 'evening' ),
    ( astronomy.BODY_MERCURY, "2017-11-14T00:32Z", "2017-11-24T00:32Z", 22.00, 'evening' ),
    ( astronomy.BODY_MERCURY, "2018-03-05T15:07Z", "2018-03-15T15:07Z", 18.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2018-07-02T05:24Z", "2018-07-12T05:24Z", 26.40, 'evening' ),
    ( astronomy.BODY_MERCURY, "2018-10-27T15:25Z", "2018-11-06T15:25Z", 23.30, 'evening' ),
    ( astronomy.BODY_MERCURY, "2019-02-17T01:23Z", "2019-02-27T01:23Z", 18.10, 'evening' ),
    ( astronomy.BODY_MERCURY, "2019-06-13T23:14Z", "2019-06-23T23:14Z", 25.20, 'evening' ),
    ( astronomy.BODY_MERCURY, "2019-10-10T04:00Z", "2019-10-20T04:00Z", 24.60, 'evening' ),
    ( astronomy.BODY_VENUS,   "2010-12-29T15:57Z", "2011-01-08T15:57Z", 47.00, 'morning' ),
    ( astronomy.BODY_VENUS,   "2012-08-05T08:59Z", "2012-08-15T08:59Z", 45.80, 'morning' ),
    ( astronomy.BODY_VENUS,   "2014-03-12T19:25Z", "2014-03-22T19:25Z", 46.60, 'morning' ),
    ( astronomy.BODY_VENUS,   "2015-10-16T06:57Z", "2015-10-26T06:57Z", 46.40, 'morning' ),
    ( astronomy.BODY_VENUS,   "2017-05-24T13:09Z", "2017-06-03T13:09Z", 45.90, 'morning' ),
    ( astronomy.BODY_VENUS,   "2018-12-27T04:24Z", "2019-01-06T04:24Z", 47.00, 'morning' ),
    ( astronomy.BODY_VENUS,   "2010-08-10T03:19Z", "2010-08-20T03:19Z", 46.00, 'evening' ),
    ( astronomy.BODY_VENUS,   "2012-03-17T08:03Z", "2012-03-27T08:03Z", 46.00, 'evening' ),
    ( astronomy.BODY_VENUS,   "2013-10-22T08:00Z", "2013-11-01T08:00Z", 47.10, 'evening' ),
    ( astronomy.BODY_VENUS,   "2015-05-27T18:46Z", "2015-06-06T18:46Z", 45.40, 'evening' ),
    ( astronomy.BODY_VENUS,   "2017-01-02T13:19Z", "2017-01-12T13:19Z", 47.10, 'evening' ),
    ( astronomy.BODY_VENUS,   "2018-08-07T17:02Z", "2018-08-17T17:02Z", 45.90, 'evening' )
]

def ParseDate(text):
    m = re.match(r'^(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z$', text)
    if not m:
        print('ParseDate: invalid date text "{}"'.format(text))
        raise Exception('Bad elongation test data')
    year = int(m.group(1))
    month = int(m.group(2))
    day = int(m.group(3))
    hour = int(m.group(4))
    minute = int(m.group(5))
    return astronomy.Time.Make(year, month, day, hour, minute, 0)

def TestMaxElong(body, searchText, eventText, angle, visiblity):
    name = astronomy.BodyName[body]
    searchTime = ParseDate(searchText)
    eventTime = ParseDate(eventText)
    evt = astronomy.SearchMaxElongation(body, searchTime)
    if evt is None:
        print('TestMaxElong({} {}): SearchMaxElongation failed.'.format(name, searchText))
        return 1
    hour_diff = 24.0 * abs(evt.time.tt - eventTime.tt)
    arcmin_diff = 60.0 * abs(evt.elongation - angle)
    print('TestMaxElong: {:<7s} {:<7s} elong={:5.2f} ({:4.2f} arcmin, {:5.3f} hours)'.format(name, visiblity, evt.elongation, arcmin_diff, hour_diff))
    if hour_diff > 0.6:
        print('TestMaxElong({} {}): EXCESSIVE HOUR ERROR.'.format(name, searchText))
        return 1
    if arcmin_diff > 3.4:
        print('TestMaxElong({} {}): EXCESSIVE ARCMIN ERROR.'.format(name, searchText))
        return 1
    return 0

def SearchElongTest():
    for (body, searchText, eventText, angle, visibility) in ElongTestData:
        if 0 != TestMaxElong(body, searchText, eventText, angle, visibility):
            return 1
    return 0

def Test_Elongation():
    if 0 != TestElongFile('longitude/opposition_2018.txt', 0.0): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_MERCURY, "temp/py_longitude_Mercury.txt", "inf"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_VENUS,   "temp/py_longitude_Venus.txt",   "inf"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_MARS,    "temp/py_longitude_Mars.txt",    "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_JUPITER, "temp/py_longitude_Jupiter.txt", "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_SATURN,  "temp/py_longitude_Saturn.txt",  "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_URANUS,  "temp/py_longitude_Uranus.txt",  "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_NEPTUNE, "temp/py_longitude_Neptune.txt", "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.BODY_PLUTO,   "temp/py_longitude_Pluto.txt",   "opp"): return 1
    if 0 != SearchElongTest(): return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def ParseJplHorizonsDateTime(line):
    m = re.match(r'^\s*(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s+(.*)$', line)
    if not m:
        return None, None
    year = int(m.group(1))
    month = 1 + ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].index(m.group(2))
    day = int(m.group(3))
    hour = int(m.group(4))
    minute = int(m.group(5))
    rest = m.group(6)
    time = astronomy.Time.Make(year, month, day, hour, minute, 0)
    return time, rest

def CheckMagnitudeData(body, filename):
    limit = 0.012
    sum_squared_diff = 0.0
    with open(filename, 'rt') as infile:
        count = lnum = 0
        for line in infile:
            lnum += 1
            line = line.strip()
            (time, rest) = ParseJplHorizonsDateTime(line)
            if (time is not None) and (rest is not None) and not ('n.a.' in rest):
                data = [float(t) for t in rest.split()]
                if len(data) != 7:
                    print('CheckMagnitudeData({} line {}): invalid data format'.format(filename, lnum))
                    return 1
                (mag, sbrt, dist, rdot, delta, deldot, phase_angle) = data
                illum = astronomy.Illumination(body, time)
                diff = illum.mag - mag
                if abs(diff) > limit:
                    print('CheckMagnitudeData({} line {}): EXCESSIVE ERROR: correct mag={}, calc mag={}'.format(filename, lnum, mag, illum.mag))
                    return 1
                sum_squared_diff += diff * diff
                if count == 0:
                    diff_lo = diff_hi = diff
                else:
                    diff_lo = min(diff_lo, diff)
                    diff_hi = max(diff_hi, diff)
                count += 1

        if count == 0:
            print('CheckMagnitudeData: Did not find any data in file: {}'.format(filename))
            return 1
    rms = math.sqrt(sum_squared_diff / count)
    print('CheckMagnitudeData: {:<21s} {:5d} rows diff_lo={:0.4f} diff_hi={:0.4f} rms={:0.4f}'.format(filename, count, diff_lo, diff_hi, rms))
    return 0

def CheckSaturn():
    # JPL Horizons does not include Saturn's rings in its magnitude models.
    # I still don't have authoritative test data for Saturn's magnitude.
    # For now, I just test for consistency with Paul Schlyter's formulas at:
    # http://www.stjarnhimlen.se/comp/ppcomp.html#15
    data = [
        ( "1972-01-01T00:00Z", -0.31904865,  +24.50061220 ),
        ( "1980-01-01T00:00Z", +0.85213663,   -1.85761461 ),
        ( "2009-09-04T00:00Z", +1.01626809,   +0.08380716 ),
        ( "2017-06-15T00:00Z", -0.12318790,  -26.60871409 ),
        ( "2019-05-01T00:00Z", +0.32954097,  -23.53880802 ),
        ( "2025-09-25T00:00Z", +0.51286575,   +1.52327932 ),
        ( "2032-05-15T00:00Z", -0.04652109,  +26.95717765 )
    ]
    for (dtext, mag, tilt) in data:
        time = ParseDate(dtext)
        illum = astronomy.Illumination(astronomy.BODY_SATURN, time)
        print('Saturn: date={}  calc mag={:12.8f}  ring_tilt={:12.8f}'.format(dtext, illum.mag, illum.ring_tilt))
        mag_diff = abs(illum.mag - mag)
        if mag_diff > 1.0e-8:
            print('CheckSaturn: Excessive magnitude error {}'.format(mag_diff))
            return 1
        tilt_diff = abs(illum.ring_tilt - tilt)
        if (tilt_diff > 1.0e-8):
            print('CheckSaturn: Excessive ring tilt error {}'.format(tilt_diff))
            return 1
    return 0

def TestMaxMag(body, filename):
    # Example of input data:
    #
    # 2001-02-21T08:00Z 2001-02-27T08:00Z 23.17 19.53 -4.84 
    #
    # JPL Horizons test data has limited floating point precision in the magnitude values.
    # There is a pair of dates for the beginning and end of the max magnitude period,
    # given the limited precision. We pick the point halfway between as the supposed max magnitude time.
    with open(filename, 'rt') as infile:
        lnum = 0
        search_time = astronomy.Time.Make(2001, 1, 1, 0, 0, 0)
        for line in infile:
            lnum += 1
            line = line.strip()
            tokenlist = line.split()
            if len(tokenlist) == 5:
                time1 = ParseDate(tokenlist[0])
                time2 = ParseDate(tokenlist[1])
                if time1 and time2:
                    center_time = time1.AddDays(0.5*(time2.ut - time1.ut))
                    correct_mag = float(tokenlist[4])
                    illum = astronomy.SearchPeakMagnitude(body, search_time)
                    mag_diff = abs(illum.mag - correct_mag)
                    hours_diff = 24.0 * abs(illum.time.ut - center_time.ut)
                    print('TestMaxMag: mag_diff={:0.3f}, hours_diff={:0.3f}'.format(mag_diff, hours_diff))
                    if hours_diff > 7.1:
                        print('TestMaxMag({} line {}): EXCESSIVE TIME DIFFERENCE.'.format(filename, lnum))
                        return 1
                    if mag_diff > 0.005:
                        print('TestMaxMag({} line {}): EXCESSIVE MAGNITUDE DIFFERENCE.'.format(filename, lnum))
                        return 1
                    search_time = time2
    print('TestMaxMag: processed {} lines from file {}'.format(lnum, filename))
    return 0


def Test_Magnitude():
    nfailed = 0
    nfailed += CheckMagnitudeData(astronomy.BODY_SUN,     'magnitude/Sun.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_MOON,    'magnitude/Moon.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_MERCURY, 'magnitude/Mercury.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_VENUS,   'magnitude/Venus.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_MARS,    'magnitude/Mars.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_JUPITER, 'magnitude/Jupiter.txt')
    nfailed += CheckSaturn()
    nfailed += CheckMagnitudeData(astronomy.BODY_URANUS,  'magnitude/Uranus.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_NEPTUNE, 'magnitude/Neptune.txt')
    nfailed += CheckMagnitudeData(astronomy.BODY_PLUTO,   'magnitude/Pluto.txt')
    nfailed += TestMaxMag(astronomy.BODY_VENUS, 'magnitude/maxmag_Venus.txt')
    if nfailed > 0:
        print('Test_Magnitude: failed {} test(s).'.format(nfailed))
        return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def Test_RiseSet(filename):
    sum_minutes = 0.0
    max_minutes = 0.0
    nudge_days = 0.01
    observer = None
    current_body = None
    a_dir = 0
    b_dir = 0
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            line = line.strip()
            # Moon  103 -61 1944-01-02T17:08Z s
            # Moon  103 -61 1944-01-03T05:47Z r
            m = re.match(r'^([A-Za-z]+)\s+(-?[0-9\.]+)\s+(-?[0-9\.]+)\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([sr])$', line)
            if not m:
                print('Test_RiseSet({} line {}): invalid data format'.format(filename, lnum))
                return 1
            name = m.group(1)
            longitude = float(m.group(2))
            latitude = float(m.group(3))
            year = int(m.group(4))
            month = int(m.group(5))
            day = int(m.group(6))
            hour = int(m.group(7))
            minute = int(m.group(8))
            kind = m.group(9)
            correct_time = astronomy.Time.Make(year, month, day, hour, minute, 0)
            direction = astronomy.DIRECTION_RISE if kind == 'r' else astronomy.DIRECTION_SET
            body = astronomy.BodyCode(name)
            if body < 0:
                print('Test_RiseSet({} line {}): invalid body name "{}"'.format(filename, lnum, name))
                return 1

            # Every time we see a new geographic location, start a new iteration
            # of finding all rise/set times for that UTC calendar year.
            if (observer is None) or (observer.latitude != latitude) or (observer.longitude != longitude) or (current_body != body):
                current_body = body
                observer = astronomy.Observer(latitude, longitude, 0)
                r_search_date = s_search_date = astronomy.Time.Make(year, 1, 1, 0, 0, 0)
                b_evt = None
                print('Test_RiseSet: {:<7s} lat={:0.1f} lon={:0.1f}'.format(name, latitude, longitude))

            if b_evt is not None:
                # Recycle the second event from the previous iteration as the first event.
                a_evt = b_evt
                a_dir = b_dir
                b_evt = None
            else:
                r_evt = astronomy.SearchRiseSet(body, observer, astronomy.DIRECTION_RISE, r_search_date, 366.0)
                if r_evt is None:
                    print('Test_RiseSet({} line {}): rise search failed'.format(filename, lnum))
                    return 1
                s_evt = astronomy.SearchRiseSet(body, observer, astronomy.DIRECTION_SET, s_search_date, 366.0)
                if s_evt is None:
                    print('Test_RiseSet({} line {}): set search failed'.format(filename, lnum))
                    return 1
                # Expect the current event to match the earlier of the found times.
                if r_evt.tt < s_evt.tt:
                    a_evt = r_evt
                    b_evt = s_evt
                    a_dir = astronomy.DIRECTION_RISE
                    b_dir = astronomy.DIRECTION_SET
                else:
                    a_evt = s_evt
                    b_evt = r_evt
                    a_dir = astronomy.DIRECTION_SET
                    b_dir = astronomy.DIRECTION_RISE                
                # Nudge the event times forward a tiny amount.
                r_search_date = r_evt.AddDays(nudge_days)
                s_search_date = s_evt.AddDays(nudge_days)

            if a_dir != direction:
                print('Test_RiseSet({} line {}): expected dir={} but found {}'.format(filename, lnum, a_dir, direction))
                return 1

            error_minutes = (24.0 * 60.0) * abs(a_evt.tt - correct_time.tt)
            sum_minutes += error_minutes ** 2
            max_minutes = max(max_minutes, error_minutes)
            if error_minutes > 0.56:
                print('Test_RiseSet({} line {}): excessive prediction time error = {} minutes.'.format(filename, lnum, error_minutes))
                print('    correct = {}, calculated = {}'.format(correct_time, a_evt))
                return 1

    rms_minutes = math.sqrt(sum_minutes / lnum)
    print('Test_RiseSet: passed {} lines: time errors in minutes: rms={:0.4f}, max={:0.4f}'.format(lnum, rms_minutes, max_minutes))
    return 0

#-----------------------------------------------------------------------------------------------------------

def LunarApsis(filename):
    max_minutes = 0.0
    max_km = 0.0
    with open(filename, 'rt') as infile:
        start_time = astronomy.Time.Make(2001, 1, 1, 0, 0, 0)
        lnum = 0
        for line in infile:
            lnum += 1
            if lnum == 1:
                apsis = astronomy.SearchLunarApsis(start_time)
            else:
                apsis = astronomy.NextLunarApsis(apsis)
            tokenlist = line.split()
            if len(tokenlist) != 3:
                print('LunarApsis({} line {}): invalid data format'.format(filename, lnum))
                return 1
            correct_time = ParseDate(tokenlist[1])
            if not correct_time:
                print('LunarApsis({} line {}): invalid time'.format(filename, lnum))
                return 1
            kind = int(tokenlist[0])
            dist_km = float(tokenlist[2])
            diff_minutes = (24.0 * 60.0) * abs(apsis.time.ut - correct_time.ut)
            diff_km = abs(apsis.dist_km - dist_km)
            if diff_minutes > 35.0:
                print('LunarApsis({} line {}): Excessive time error = {} minutes.'.format(filename, lnum, diff_minutes))
                return 1
            if diff_km > 25.0:
                print('LunarApsis({} line {}): Excessive distance error = {} km.'.format(filename, lnum, diff_km))
                return 1
            max_minutes = max(max_minutes, diff_minutes)
            max_km = max(max_km, diff_km)
    print('LunarApsis: found {} events, max time error = {:0.3f} minutes, max distance error = {:0.3f} km.'.format(lnum, max_minutes, max_km))
    return 0
            

def Test_Apsis():
    if 0 != LunarApsis('apsides/moon.txt'):
        return 1
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
    if sys.argv[1] == 'elongation':
        sys.exit(Test_Elongation())
    if sys.argv[1] == 'magnitude':
        sys.exit(Test_Magnitude())
    if sys.argv[1] == 'apsis':
        sys.exit(Test_Apsis())

if len(sys.argv) == 3:
    if sys.argv[1] == 'seasons':
        sys.exit(Test_Seasons(sys.argv[2]))
    if sys.argv[1] == 'moonphase':
        sys.exit(Test_MoonPhase(sys.argv[2]))
    if sys.argv[1] == 'riseset':
        sys.exit(Test_RiseSet(sys.argv[2]))

print('test.py: Invalid command line arguments.')
sys.exit(1)
