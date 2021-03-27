#!/usr/bin/env python3
import sys
import math
import re
import os
sys.path.append('../source/python')
import astronomy

#-----------------------------------------------------------------------------------------------------------

Verbose = False

def Debug(text):
    if Verbose:
        print(text)

def v(x):
    # Verify that a number is really numeric
    if not isinstance(x, (int, float)):
        raise Exception('Not a numeric type: {}'.format(x))
    if not math.isfinite(x):
        raise Exception('Not a finite numeric value: {}'.format(x))
    return x

def vabs(x):
    return abs(v(x))

def vmax(a, b):
    return max(v(a), v(b))

def vmin(a, b):
    return min(v(a), v(b))

def sqrt(x):
    return v(math.sqrt(v(x)))

def AssertGoodTime(text, correct):
    time = astronomy.Time.Parse(text)
    check = str(time)
    if check != correct:
        print('Python AssertGoodTime FAILURE: parsed "{}", got "{}", expected "{}"'.format(text, check, correct))
        sys.exit(1)
    Debug('PY AssertGoodTime: "{}" OK'.format(text))

def AssertBadTime(text):
    try:
        astronomy.Time.Parse(text)
    except astronomy.DateTimeFormatError:
        Debug('PY AssertBadTime: "{}" OK'.format(text))
    else:
        print('PY AssertBadTime FAILURE: should not have parsed "{}"'.format(text))
        sys.exit(1)

def AstroTime():
    expected_ut = 6910.270978506945
    expected_tt = 6910.271800214368
    time = astronomy.Time.Make(2018, 12, 2, 18, 30, 12.543)
    diff = time.ut - expected_ut
    if vabs(diff) > 1.0e-12:
        print('PY AstroTime: excessive UT error {}'.format(diff))
        sys.exit(1)
    diff = time.tt - expected_tt
    if vabs(diff) > 1.0e-12:
        print('PY AstroTime: excessive TT error {}'.format(diff))
        sys.exit(1)
    s = str(time.Utc())
    if s != '2018-12-02 18:30:12.543000':
        print('PY AstroTime: Utc() returned incorrect string "{}"'.format(s))
        sys.exit(1)
    time = astronomy.Time.Make(2018, 12, 31, 23, 59, 59.9994)
    s = str(time)
    if s != '2018-12-31T23:59:59.999Z':
        print('PY AstroTime: expected 2018-12-31T23:59:59.999Z but found {}'.format(s))
        sys.exit(1)
    time = astronomy.Time.Make(2018, 12, 31, 23, 59, 59.9995)
    s = str(time)
    if s != '2019-01-01T00:00:00.000Z':
        print('PY AstroTime: expected 2019-01-01T00:00:00.000Z but found {}'.format(s))
        sys.exit(1)
    print('PY Current time =', astronomy.Time.Now())
    AssertGoodTime('2015-12-31', '2015-12-31T00:00:00.000Z')
    AssertGoodTime('2015-12-31T23:45Z', '2015-12-31T23:45:00.000Z')
    AssertGoodTime('2015-01-02T23:45:17Z', '2015-01-02T23:45:17.000Z')
    AssertGoodTime('1971-03-17T03:30:55.976Z', '1971-03-17T03:30:55.976Z')
    AssertBadTime('')
    AssertBadTime('1971-13-01')
    AssertBadTime('1971-12-32')
    AssertBadTime('1971-12-31T24:00:00Z')
    AssertBadTime('1971-12-31T23:60:00Z')
    AssertBadTime('1971-12-31T23:00:60Z')
    AssertBadTime('1971-03-17T03:30:55.976')
    return 0

#-----------------------------------------------------------------------------------------------------------

def GeoMoon():
    time = astronomy.Time.Make(2019, 6, 24, 15, 45, 37)
    vec = astronomy.GeoMoon(time)
    print('PY GeoMoon: vec = {:0.16f}, {:0.16f}, {:0.16f}'.format(vec.x, vec.y, vec.z))
    # Correct values obtained from C version of GeoMoon calculation
    cx, cy, cz = +0.002674037026701135, -0.0001531610316600666, -0.0003150159927069429
    dx, dy, dz = vec.x - cx, vec.y - cy, vec.z - cz
    diff = sqrt(dx*dx + dy*dy + dz*dz)
    print('PY GeoMoon: diff = {}'.format(diff))
    if diff > 4.34e-19:
        print('PY GeoMoon: EXCESSIVE ERROR')
        return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def AstroCheck(printflag):
    time = astronomy.Time.Make(1700, 1, 1, 0, 0, 0)
    stop = astronomy.Time.Make(2200, 1, 1, 0, 0, 0)
    observer = astronomy.Observer(29, -81, 10)
    if printflag:
        print('o {:0.6f} {:0.6f} {:0.6f}'.format(observer.latitude, observer.longitude, observer.height))
    dt = 10 + math.pi/100
    bodylist = [
        astronomy.Body.Sun, astronomy.Body.Moon, astronomy.Body.Mercury, astronomy.Body.Venus,
        astronomy.Body.Earth, astronomy.Body.Mars, astronomy.Body.Jupiter, astronomy.Body.Saturn,
        astronomy.Body.Uranus, astronomy.Body.Neptune, astronomy.Body.Pluto,
        astronomy.Body.SSB, astronomy.Body.EMB
    ]

    while time.tt < stop.tt:
        for body in bodylist:
            name = body.name
            if body != astronomy.Body.Moon:
                pos = astronomy.HelioVector(body, time)
                if printflag:
                    print('v {} {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(name, pos.t.tt, pos.x, pos.y, pos.z))
                if body != astronomy.Body.Earth and body != astronomy.Body.EMB and body != astronomy.Body.SSB:
                    j2000 = astronomy.Equator(body, time, observer, False, False)
                    ofdate = astronomy.Equator(body, time, observer, True, True)
                    hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.Refraction.Airless)
                    if printflag:
                        print('s {} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(name, time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude))
        pos = astronomy.GeoMoon(time)
        if printflag:
            print('v GM {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(pos.t.tt, pos.x, pos.y, pos.z))
        j2000 = astronomy.Equator(astronomy.Body.Moon, time, observer, False, False)
        ofdate = astronomy.Equator(astronomy.Body.Moon, time, observer, True, True)
        hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.Refraction.Airless)
        if printflag:
            print('s GM {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f} {:0.16f}'.format(time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude))
        time = time.AddDays(dt)
    return 0

#-----------------------------------------------------------------------------------------------------------

def Seasons(filename = 'seasons/seasons.txt'):
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
                print('PY Seasons: Invalid data on line {} of file {}'.format(lnum, filename))
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
                    print('PY Seasons: {} line {}: Invalid equinox date in test data'.format(filename, lnum))
                    return 1
            elif name == 'Solstice':
                if month == 6:
                    calc_time = seasons.jun_solstice
                    jun_count += 1
                elif month == 12:
                    calc_time = seasons.dec_solstice
                    dec_count += 1
                else:
                    print('PY Seasons: {} line {}: Invalid solstice date in test data'.format(filename, lnum))
                    return 1
            elif name == 'Aphelion':
                continue # not yet calculated
            elif name == 'Perihelion':
                continue # not yet calculated
            else:
                print('PY Seasons: {} line {}: unknown event type {}'.format(filename, lnum, name))
                return 1

            # Verify that the calculated time matches the correct time for this event.
            diff_minutes = (24.0 * 60.0) * vabs(calc_time.tt - correct_time.tt)
            if diff_minutes > max_minutes:
                max_minutes = diff_minutes

            if diff_minutes > 2.37:
                print('PY Seasons: {} line {}: excessive error ({}): {} minutes.'.format(filename, lnum, name, diff_minutes))
                return 1
    print('PY Seasons: verified {} lines from file {} : max error minutes = {:0.3f}'.format(lnum, filename, max_minutes))
    print('PY Seasons: Event counts: mar={}, jun={}, sep={}, dec={}'.format(mar_count, jun_count, sep_count, dec_count))
    return 0

#-----------------------------------------------------------------------------------------------------------

def MoonPhase(filename = 'moonphase/moonphases.txt'):
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
                print('PY MoonPhase: invalid data format in {} line {}'.format(filename, lnum))
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
            degree_error = vabs(angle - expected_elong)
            if degree_error > 180.0:
                degree_error = 360.0 - degree_error
            arcmin = 60.0 * degree_error
            if arcmin > 1.0:
                print('PY MoonPhase({} line {}): EXCESSIVE ANGULAR ERROR: {} arcmin'.format(filename, lnum, arcmin))
                return 1
            max_arcmin = vmax(max_arcmin, arcmin)

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
                    print('PY MoonPhase({} line {}): SearchMoonQuarter returned quarter {}, but expected {}.'.format(filename, lnum, mq.quarter, expected_quarter))
                    return 1

            quarter_count += 1

            # Make sure the time matches what we expect.
            diff_seconds = vabs(mq.time.tt - expected_time.tt) * (24.0 * 3600.0)
            if diff_seconds > threshold_seconds:
                print('PY MoonPhase({} line {}): excessive time error {:0.3f} seconds.'.format(filename, lnum, diff_seconds))
                return 1

            maxdiff = vmax(maxdiff, diff_seconds)

    print('PY MoonPhase: passed {} lines for file {} : max_arcmin = {:0.6f}, maxdiff = {:0.3f} seconds, {} quarters.'
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
                print('PY TestElongFile({} line {}): invalid data format'.format(filename, lnum))
                return 1
            year = int(m.group(1))
            month = int(m.group(2))
            day = int(m.group(3))
            hour = int(m.group(4))
            minute = int(m.group(5))
            name = m.group(6)
            body = astronomy.BodyCode(name)
            if body.value == astronomy.Body.Invalid:
                print('PY TestElongFile({} line {}): invalid body name "{}"'.format(filename, lnum, name))
                return 1
            search_time = astronomy.Time.Make(year, 1, 1, 0, 0, 0)
            expected_time = astronomy.Time.Make(year, month, day, hour, minute, 0)
            found_time = astronomy.SearchRelativeLongitude(body, targetRelLon, search_time)
            if found_time is None:
                print('PY TestElongFile({} line {}): SearchRelativeLongitude failed.'.format(filename, lnum))
                return 1
            diff_minutes = (24.0 * 60.0) * (found_time.tt - expected_time.tt)
            Debug('PY TestElongFile: {:<7s} error = {:6.3} minutes'.format(name, diff_minutes))
            if vabs(diff_minutes) > 15.0:
                print('PY TestElongFile({} line {}): EXCESSIVE ERROR.'.format(filename, lnum))
                return 1
    print('PY TestElongFile: passed {} rows of data'.format(lnum))
    return 0

def TestPlanetLongitudes(body, outFileName, zeroLonEventName):
    startYear = 1700
    stopYear = 2200
    rlon = 0.0
    sum_diff = 0.0
    count = 0
    name = body.name
    with open(outFileName, 'wt') as outfile:
        time = astronomy.Time.Make(startYear, 1, 1, 0, 0, 0)
        stopTime = astronomy.Time.Make(stopYear, 1, 1, 0, 0, 0)
        while time.tt < stopTime.tt:
            count += 1
            event = zeroLonEventName if rlon == 0.0 else 'sup'
            found_time = astronomy.SearchRelativeLongitude(body, rlon, time)
            if found_time is None:
                print('PY TestPlanetLongitudes({}): SearchRelativeLongitudes failed'.format(name))
                return 1
            if count >= 2:
                # Check for consistent intervals.
                # Mainly I don't want to skip over an event!
                day_diff = found_time.tt - time.tt
                sum_diff += day_diff
                if count == 2:
                    min_diff = max_diff = day_diff
                else:
                    min_diff = vmin(min_diff, day_diff)
                    max_diff = vmax(max_diff, day_diff)
            geo = astronomy.GeoVector(body, found_time, True)
            dist = geo.Length()
            outfile.write('e {} {} {:0.16f} {:0.16f}\n'.format(name, event, found_time.tt, dist))
            # Search for the opposite longitude vent next time.
            time = found_time
            rlon = 180.0 - rlon
    if body == astronomy.Body.Mercury:
        thresh = 1.65
    elif body == astronomy.Body.Mars:
        thresh = 1.30
    else:
        thresh = 1.07
    ratio = max_diff / min_diff
    Debug('PY TestPlanetLongitudes({:<7s}): {:5d} events, ratio={:5.3f}, file: {}'.format(name, count, ratio, outFileName))
    if ratio > thresh:
        print('PY TestPlanetLongitudes({}): EXCESSIVE EVENT INTERVAL RATIO'.format(name))
        return 1
    return 0

ElongTestData = [
    # Max elongation data obtained from:
    # http://www.skycaramba.com/greatest_elongations.shtml
    ( astronomy.Body.Mercury, "2010-01-17T05:22Z", "2010-01-27T05:22Z", 24.80, 'morning' ),
    ( astronomy.Body.Mercury, "2010-05-16T02:15Z", "2010-05-26T02:15Z", 25.10, 'morning' ),
    ( astronomy.Body.Mercury, "2010-09-09T17:24Z", "2010-09-19T17:24Z", 17.90, 'morning' ),
    ( astronomy.Body.Mercury, "2010-12-30T14:33Z", "2011-01-09T14:33Z", 23.30, 'morning' ),
    ( astronomy.Body.Mercury, "2011-04-27T19:03Z", "2011-05-07T19:03Z", 26.60, 'morning' ),
    ( astronomy.Body.Mercury, "2011-08-24T05:52Z", "2011-09-03T05:52Z", 18.10, 'morning' ),
    ( astronomy.Body.Mercury, "2011-12-13T02:56Z", "2011-12-23T02:56Z", 21.80, 'morning' ),
    ( astronomy.Body.Mercury, "2012-04-08T17:22Z", "2012-04-18T17:22Z", 27.50, 'morning' ),
    ( astronomy.Body.Mercury, "2012-08-06T12:04Z", "2012-08-16T12:04Z", 18.70, 'morning' ),
    ( astronomy.Body.Mercury, "2012-11-24T22:55Z", "2012-12-04T22:55Z", 20.60, 'morning' ),
    ( astronomy.Body.Mercury, "2013-03-21T22:02Z", "2013-03-31T22:02Z", 27.80, 'morning' ),
    ( astronomy.Body.Mercury, "2013-07-20T08:51Z", "2013-07-30T08:51Z", 19.60, 'morning' ),
    ( astronomy.Body.Mercury, "2013-11-08T02:28Z", "2013-11-18T02:28Z", 19.50, 'morning' ),
    ( astronomy.Body.Mercury, "2014-03-04T06:38Z", "2014-03-14T06:38Z", 27.60, 'morning' ),
    ( astronomy.Body.Mercury, "2014-07-02T18:22Z", "2014-07-12T18:22Z", 20.90, 'morning' ),
    ( astronomy.Body.Mercury, "2014-10-22T12:36Z", "2014-11-01T12:36Z", 18.70, 'morning' ),
    ( astronomy.Body.Mercury, "2015-02-14T16:20Z", "2015-02-24T16:20Z", 26.70, 'morning' ),
    ( astronomy.Body.Mercury, "2015-06-14T17:10Z", "2015-06-24T17:10Z", 22.50, 'morning' ),
    ( astronomy.Body.Mercury, "2015-10-06T03:20Z", "2015-10-16T03:20Z", 18.10, 'morning' ),
    ( astronomy.Body.Mercury, "2016-01-28T01:22Z", "2016-02-07T01:22Z", 25.60, 'morning' ),
    ( astronomy.Body.Mercury, "2016-05-26T08:45Z", "2016-06-05T08:45Z", 24.20, 'morning' ),
    ( astronomy.Body.Mercury, "2016-09-18T19:27Z", "2016-09-28T19:27Z", 17.90, 'morning' ),
    ( astronomy.Body.Mercury, "2017-01-09T09:42Z", "2017-01-19T09:42Z", 24.10, 'morning' ),
    ( astronomy.Body.Mercury, "2017-05-07T23:19Z", "2017-05-17T23:19Z", 25.80, 'morning' ),
    ( astronomy.Body.Mercury, "2017-09-02T10:14Z", "2017-09-12T10:14Z", 17.90, 'morning' ),
    ( astronomy.Body.Mercury, "2017-12-22T19:48Z", "2018-01-01T19:48Z", 22.70, 'morning' ),
    ( astronomy.Body.Mercury, "2018-04-19T18:17Z", "2018-04-29T18:17Z", 27.00, 'morning' ),
    ( astronomy.Body.Mercury, "2018-08-16T20:35Z", "2018-08-26T20:35Z", 18.30, 'morning' ),
    ( astronomy.Body.Mercury, "2018-12-05T11:34Z", "2018-12-15T11:34Z", 21.30, 'morning' ),
    ( astronomy.Body.Mercury, "2019-04-01T19:40Z", "2019-04-11T19:40Z", 27.70, 'morning' ),
    ( astronomy.Body.Mercury, "2019-07-30T23:08Z", "2019-08-09T23:08Z", 19.00, 'morning' ),
    ( astronomy.Body.Mercury, "2019-11-18T10:31Z", "2019-11-28T10:31Z", 20.10, 'morning' ),
    ( astronomy.Body.Mercury, "2010-03-29T23:32Z", "2010-04-08T23:32Z", 19.40, 'evening' ),
    ( astronomy.Body.Mercury, "2010-07-28T01:03Z", "2010-08-07T01:03Z", 27.40, 'evening' ),
    ( astronomy.Body.Mercury, "2010-11-21T15:42Z", "2010-12-01T15:42Z", 21.50, 'evening' ),
    ( astronomy.Body.Mercury, "2011-03-13T01:07Z", "2011-03-23T01:07Z", 18.60, 'evening' ),
    ( astronomy.Body.Mercury, "2011-07-10T04:56Z", "2011-07-20T04:56Z", 26.80, 'evening' ),
    ( astronomy.Body.Mercury, "2011-11-04T08:40Z", "2011-11-14T08:40Z", 22.70, 'evening' ),
    ( astronomy.Body.Mercury, "2012-02-24T09:39Z", "2012-03-05T09:39Z", 18.20, 'evening' ),
    ( astronomy.Body.Mercury, "2012-06-21T02:00Z", "2012-07-01T02:00Z", 25.70, 'evening' ),
    ( astronomy.Body.Mercury, "2012-10-16T21:59Z", "2012-10-26T21:59Z", 24.10, 'evening' ),
    ( astronomy.Body.Mercury, "2013-02-06T21:24Z", "2013-02-16T21:24Z", 18.10, 'evening' ),
    ( astronomy.Body.Mercury, "2013-06-02T16:45Z", "2013-06-12T16:45Z", 24.30, 'evening' ),
    ( astronomy.Body.Mercury, "2013-09-29T09:59Z", "2013-10-09T09:59Z", 25.30, 'evening' ),
    ( astronomy.Body.Mercury, "2014-01-21T10:00Z", "2014-01-31T10:00Z", 18.40, 'evening' ),
    ( astronomy.Body.Mercury, "2014-05-15T07:06Z", "2014-05-25T07:06Z", 22.70, 'evening' ),
    ( astronomy.Body.Mercury, "2014-09-11T22:20Z", "2014-09-21T22:20Z", 26.40, 'evening' ),
    ( astronomy.Body.Mercury, "2015-01-04T20:26Z", "2015-01-14T20:26Z", 18.90, 'evening' ),
    ( astronomy.Body.Mercury, "2015-04-27T04:46Z", "2015-05-07T04:46Z", 21.20, 'evening' ),
    ( astronomy.Body.Mercury, "2015-08-25T10:20Z", "2015-09-04T10:20Z", 27.10, 'evening' ),
    ( astronomy.Body.Mercury, "2015-12-19T03:11Z", "2015-12-29T03:11Z", 19.70, 'evening' ),
    ( astronomy.Body.Mercury, "2016-04-08T14:00Z", "2016-04-18T14:00Z", 19.90, 'evening' ),
    ( astronomy.Body.Mercury, "2016-08-06T21:24Z", "2016-08-16T21:24Z", 27.40, 'evening' ),
    ( astronomy.Body.Mercury, "2016-12-01T04:36Z", "2016-12-11T04:36Z", 20.80, 'evening' ),
    ( astronomy.Body.Mercury, "2017-03-22T10:24Z", "2017-04-01T10:24Z", 19.00, 'evening' ),
    ( astronomy.Body.Mercury, "2017-07-20T04:34Z", "2017-07-30T04:34Z", 27.20, 'evening' ),
    ( astronomy.Body.Mercury, "2017-11-14T00:32Z", "2017-11-24T00:32Z", 22.00, 'evening' ),
    ( astronomy.Body.Mercury, "2018-03-05T15:07Z", "2018-03-15T15:07Z", 18.40, 'evening' ),
    ( astronomy.Body.Mercury, "2018-07-02T05:24Z", "2018-07-12T05:24Z", 26.40, 'evening' ),
    ( astronomy.Body.Mercury, "2018-10-27T15:25Z", "2018-11-06T15:25Z", 23.30, 'evening' ),
    ( astronomy.Body.Mercury, "2019-02-17T01:23Z", "2019-02-27T01:23Z", 18.10, 'evening' ),
    ( astronomy.Body.Mercury, "2019-06-13T23:14Z", "2019-06-23T23:14Z", 25.20, 'evening' ),
    ( astronomy.Body.Mercury, "2019-10-10T04:00Z", "2019-10-20T04:00Z", 24.60, 'evening' ),
    ( astronomy.Body.Venus,   "2010-12-29T15:57Z", "2011-01-08T15:57Z", 47.00, 'morning' ),
    ( astronomy.Body.Venus,   "2012-08-05T08:59Z", "2012-08-15T08:59Z", 45.80, 'morning' ),
    ( astronomy.Body.Venus,   "2014-03-12T19:25Z", "2014-03-22T19:25Z", 46.60, 'morning' ),
    ( astronomy.Body.Venus,   "2015-10-16T06:57Z", "2015-10-26T06:57Z", 46.40, 'morning' ),
    ( astronomy.Body.Venus,   "2017-05-24T13:09Z", "2017-06-03T13:09Z", 45.90, 'morning' ),
    ( astronomy.Body.Venus,   "2018-12-27T04:24Z", "2019-01-06T04:24Z", 47.00, 'morning' ),
    ( astronomy.Body.Venus,   "2010-08-10T03:19Z", "2010-08-20T03:19Z", 46.00, 'evening' ),
    ( astronomy.Body.Venus,   "2012-03-17T08:03Z", "2012-03-27T08:03Z", 46.00, 'evening' ),
    ( astronomy.Body.Venus,   "2013-10-22T08:00Z", "2013-11-01T08:00Z", 47.10, 'evening' ),
    ( astronomy.Body.Venus,   "2015-05-27T18:46Z", "2015-06-06T18:46Z", 45.40, 'evening' ),
    ( astronomy.Body.Venus,   "2017-01-02T13:19Z", "2017-01-12T13:19Z", 47.10, 'evening' ),
    ( astronomy.Body.Venus,   "2018-08-07T17:02Z", "2018-08-17T17:02Z", 45.90, 'evening' )
]

def TestMaxElong(body, searchText, eventText, angle, visibility):
    name = body.name
    searchTime = astronomy.Time.Parse(searchText)
    eventTime = astronomy.Time.Parse(eventText)
    evt = astronomy.SearchMaxElongation(body, searchTime)
    if evt is None:
        print('PY TestMaxElong({} {}): SearchMaxElongation failed.'.format(name, searchText))
        return 1
    if evt.visibility != visibility:
        print('PY TestMaxElong({} {}): SearchMaxElongation returned visibility {}, but expected {}'.format(name, searchText, evt.visibility.name, visibility.name))
        return 1
    hour_diff = 24.0 * vabs(evt.time.tt - eventTime.tt)
    arcmin_diff = 60.0 * vabs(evt.elongation - angle)
    Debug('PY TestMaxElong: {:<7s} {:<7s} elong={:5.2f} ({:4.2f} arcmin, {:5.3f} hours)'.format(name, visibility.name, evt.elongation, arcmin_diff, hour_diff))
    if hour_diff > 0.603:
        print('PY TestMaxElong({} {}): EXCESSIVE HOUR ERROR.'.format(name, searchText))
        return 1
    if arcmin_diff > 3.4:
        print('PY TestMaxElong({} {}): EXCESSIVE ARCMIN ERROR.'.format(name, searchText))
        return 1
    return 0

def SearchElongTest():
    for (body, searchText, eventText, angle, visibility) in ElongTestData:
        if 0 != TestMaxElong(body, searchText, eventText, angle, astronomy.Visibility[visibility.title()]):
            return 1
    return 0

def Elongation():
    if 0 != TestElongFile('longitude/opposition_2018.txt', 0.0): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Mercury, "temp/py_longitude_Mercury.txt", "inf"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Venus,   "temp/py_longitude_Venus.txt",   "inf"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Mars,    "temp/py_longitude_Mars.txt",    "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Jupiter, "temp/py_longitude_Jupiter.txt", "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Saturn,  "temp/py_longitude_Saturn.txt",  "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Uranus,  "temp/py_longitude_Uranus.txt",  "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Neptune, "temp/py_longitude_Neptune.txt", "opp"): return 1
    if 0 != TestPlanetLongitudes(astronomy.Body.Pluto,   "temp/py_longitude_Pluto.txt",   "opp"): return 1
    if 0 != SearchElongTest(): return 1
    print('PY Elongation: PASS')
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
                    print('PY CheckMagnitudeData({} line {}): invalid data format'.format(filename, lnum))
                    return 1
                (mag, sbrt, dist, rdot, delta, deldot, phase_angle) = data
                illum = astronomy.Illumination(body, time)
                diff = illum.mag - mag
                if vabs(diff) > limit:
                    print('PY CheckMagnitudeData({} line {}): EXCESSIVE ERROR: correct mag={}, calc mag={}'.format(filename, lnum, mag, illum.mag))
                    return 1
                sum_squared_diff += diff * diff
                if count == 0:
                    diff_lo = diff_hi = diff
                else:
                    diff_lo = vmin(diff_lo, diff)
                    diff_hi = vmax(diff_hi, diff)
                count += 1

        if count == 0:
            print('PY CheckMagnitudeData: Did not find any data in file: {}'.format(filename))
            return 1
    rms = sqrt(sum_squared_diff / count)
    Debug('PY CheckMagnitudeData: {:<21s} {:5d} rows diff_lo={:0.4f} diff_hi={:0.4f} rms={:0.4f}'.format(filename, count, diff_lo, diff_hi, rms))
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
    error = 0
    for (dtext, mag, tilt) in data:
        time = astronomy.Time.Parse(dtext)
        illum = astronomy.Illumination(astronomy.Body.Saturn, time)
        Debug('PY Saturn: date={}  calc mag={:12.8f}  ring_tilt={:12.8f}'.format(dtext, illum.mag, illum.ring_tilt))
        mag_diff = vabs(illum.mag - mag)
        if mag_diff > 1.0e-4:
            print('PY CheckSaturn: Excessive magnitude error {}'.format(mag_diff))
            error = 1
        tilt_diff = vabs(illum.ring_tilt - tilt)
        if (tilt_diff > 3.0e-5):
            print('PY CheckSaturn: Excessive ring tilt error {}'.format(tilt_diff))
            error = 1
    return error

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
                time1 = astronomy.Time.Parse(tokenlist[0])
                time2 = astronomy.Time.Parse(tokenlist[1])
                if time1 and time2:
                    center_time = time1.AddDays(0.5*(time2.ut - time1.ut))
                    correct_mag = float(tokenlist[4])
                    illum = astronomy.SearchPeakMagnitude(body, search_time)
                    mag_diff = vabs(illum.mag - correct_mag)
                    hours_diff = 24.0 * vabs(illum.time.ut - center_time.ut)
                    Debug('PY TestMaxMag: mag_diff={:0.3f}, hours_diff={:0.3f}'.format(mag_diff, hours_diff))
                    if hours_diff > 7.1:
                        print('PY TestMaxMag({} line {}): EXCESSIVE TIME DIFFERENCE.'.format(filename, lnum))
                        return 1
                    if mag_diff > 0.005:
                        print('PY TestMaxMag({} line {}): EXCESSIVE MAGNITUDE DIFFERENCE.'.format(filename, lnum))
                        return 1
                    search_time = time2
    Debug('PY TestMaxMag: processed {} lines from file {}'.format(lnum, filename))
    return 0


def Magnitude():
    nfailed = 0
    nfailed += CheckMagnitudeData(astronomy.Body.Sun,     'magnitude/Sun.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Moon,    'magnitude/Moon.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Mercury, 'magnitude/Mercury.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Venus,   'magnitude/Venus.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Mars,    'magnitude/Mars.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Jupiter, 'magnitude/Jupiter.txt')
    nfailed += CheckSaturn()
    nfailed += CheckMagnitudeData(astronomy.Body.Uranus,  'magnitude/Uranus.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Neptune, 'magnitude/Neptune.txt')
    nfailed += CheckMagnitudeData(astronomy.Body.Pluto,   'magnitude/Pluto.txt')
    nfailed += TestMaxMag(astronomy.Body.Venus, 'magnitude/maxmag_Venus.txt')
    if nfailed == 0:
        print('PY Magnitude: PASS')
    else:
        print('PY Magnitude: failed {} test(s).'.format(nfailed))
        return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def RiseSet(filename = 'riseset/riseset.txt'):
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
                print('PY RiseSet({} line {}): invalid data format'.format(filename, lnum))
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
            direction = astronomy.Direction.Rise if kind == 'r' else astronomy.Direction.Set
            body = astronomy.BodyCode(name)
            if body == astronomy.Body.Invalid:
                print('PY RiseSet({} line {}): invalid body name "{}"'.format(filename, lnum, name))
                return 1

            # Every time we see a new geographic location, start a new iteration
            # of finding all rise/set times for that UTC calendar year.
            if (observer is None) or (observer.latitude != latitude) or (observer.longitude != longitude) or (current_body != body):
                current_body = body
                observer = astronomy.Observer(latitude, longitude, 0)
                r_search_date = s_search_date = astronomy.Time.Make(year, 1, 1, 0, 0, 0)
                b_evt = None
                Debug('PY RiseSet: {:<7s} lat={:0.1f} lon={:0.1f}'.format(name, latitude, longitude))

            if b_evt is not None:
                # Recycle the second event from the previous iteration as the first event.
                a_evt = b_evt
                a_dir = b_dir
                b_evt = None
            else:
                r_evt = astronomy.SearchRiseSet(body, observer, astronomy.Direction.Rise, r_search_date, 366.0)
                if r_evt is None:
                    print('PY RiseSet({} line {}): rise search failed'.format(filename, lnum))
                    return 1
                s_evt = astronomy.SearchRiseSet(body, observer, astronomy.Direction.Set, s_search_date, 366.0)
                if s_evt is None:
                    print('PY RiseSet({} line {}): set search failed'.format(filename, lnum))
                    return 1
                # Expect the current event to match the earlier of the found times.
                if r_evt.tt < s_evt.tt:
                    a_evt = r_evt
                    b_evt = s_evt
                    a_dir = astronomy.Direction.Rise
                    b_dir = astronomy.Direction.Set
                else:
                    a_evt = s_evt
                    b_evt = r_evt
                    a_dir = astronomy.Direction.Set
                    b_dir = astronomy.Direction.Rise
                # Nudge the event times forward a tiny amount.
                r_search_date = r_evt.AddDays(nudge_days)
                s_search_date = s_evt.AddDays(nudge_days)

            if a_dir != direction:
                print('PY RiseSet({} line {}): expected dir={} but found {}'.format(filename, lnum, a_dir, direction))
                return 1

            error_minutes = (24.0 * 60.0) * vabs(a_evt.tt - correct_time.tt)
            sum_minutes += error_minutes ** 2
            max_minutes = vmax(max_minutes, error_minutes)
            if error_minutes > 0.57:
                print('PY RiseSet({} line {}): excessive prediction time error = {} minutes.'.format(filename, lnum, error_minutes))
                print('    correct = {}, calculated = {}'.format(correct_time, a_evt))
                return 1

    rms_minutes = sqrt(sum_minutes / lnum)
    print('PY RiseSet: passed {} lines: time errors in minutes: rms={:0.4f}, max={:0.4f}'.format(lnum, rms_minutes, max_minutes))
    return 0

#-----------------------------------------------------------------------------------------------------------

def LunarApsis(filename = 'apsides/moon.txt'):
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
                print('PY LunarApsis({} line {}): invalid data format'.format(filename, lnum))
                return 1
            correct_time = astronomy.Time.Parse(tokenlist[1])
            if not correct_time:
                print('PY LunarApsis({} line {}): invalid time'.format(filename, lnum))
                return 1
            kind = astronomy.ApsisKind(int(tokenlist[0]))
            if apsis.kind != kind:
                print('PY LunarApsis({} line {}): Expected kind {} but found {}'.format(filename, lnum, kind, apsis.kind))
                return 1
            dist_km = float(tokenlist[2])
            diff_minutes = (24.0 * 60.0) * vabs(apsis.time.ut - correct_time.ut)
            diff_km = vabs(apsis.dist_km - dist_km)
            if diff_minutes > 35.0:
                print('PY LunarApsis({} line {}): Excessive time error = {} minutes.'.format(filename, lnum, diff_minutes))
                return 1
            if diff_km > 25.0:
                print('PY LunarApsis({} line {}): Excessive distance error = {} km.'.format(filename, lnum, diff_km))
                return 1
            max_minutes = vmax(max_minutes, diff_minutes)
            max_km = vmax(max_km, diff_km)
    print('PY LunarApsis: found {} events, max time error = {:0.3f} minutes, max distance error = {:0.3f} km.'.format(lnum, max_minutes, max_km))
    return 0

#-----------------------------------------------------------------------------------------------------------

def CompareMatrices(caller, a, b, tolerance):
    for i in range(3):
        for j in range(3):
            diff = vabs(a.rot[i][j] - b.rot[i][j])
            if diff > tolerance:
                print('PY CompareMatrices ERROR({}): matrix[{}][{}] = {}, expected {}, diff {}'.format(caller, i, j, a.rot[i][j], b.rot[i][j], diff))
                sys.exit(1)


def CompareVectors(caller, a, b, tolerance):
    diff = vabs(a.x - b.x)
    if diff > tolerance:
        print('PY CompareVectors ERROR({}): vector x = {}, expected {}, diff {}'.format(caller, a.x, b.x, diff))
        sys.exit(1)

    diff = vabs(a.y - b.y)
    if diff > tolerance:
        print('PY CompareVectors ERROR({}): vector y = {}, expected {}, diff {}'.format(caller, a.y, b.y, diff))
        sys.exit(1)

    diff = vabs(a.z - b.z)
    if diff > tolerance:
        print('PY CompareVectors ERROR({}): vector z = {}, expected {}, diff {}'.format(caller, a.z, b.z, diff))
        sys.exit(1)


def Rotation_MatrixInverse():
    a = astronomy.RotationMatrix([
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ])
    v = astronomy.RotationMatrix([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ])
    b = astronomy.InverseRotation(a)
    CompareMatrices('Rotation_MatrixInverse', b, v, 0)


def Rotation_MatrixMultiply():
    a = astronomy.RotationMatrix([
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ])

    b = astronomy.RotationMatrix([
        [10, 13, 16],
        [11, 14, 17],
        [12, 15, 18]
    ])

    v = astronomy.RotationMatrix([
        [84, 201, 318],
        [90, 216, 342],
        [96, 231, 366]
    ])

    c = astronomy.CombineRotation(b, a)
    CompareMatrices('Rotation_MatrixMultiply', c, v, 0)


def VectorDiff(a, b):
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return sqrt(dx*dx + dy*dy + dz*dz)


def Test_EQJ_ECL():
    r = astronomy.Rotation_EQJ_ECL()
    # Calculate heliocentric Earth position at a test time.
    time = astronomy.Time.Make(2019, 12, 8, 19, 39, 15)
    ev = astronomy.HelioVector(astronomy.Body.Earth, time)

    # Use the existing astronomy.Ecliptic() to calculate ecliptic vector and angles.
    ecl = astronomy.Ecliptic(ev)
    Debug('PY Test_EQJ_ECL ecl = ({}, {}, {})'.format(ecl.ex, ecl.ey, ecl.ez))

    # Now compute the same vector via rotation matrix.
    ee = astronomy.RotateVector(r, ev)
    dx = ee.x - ecl.ex
    dy = ee.y - ecl.ey
    dz = ee.z - ecl.ez
    diff = sqrt(dx*dx + dy*dy + dz*dz)
    Debug('PY Test_EQJ_ECL ee = ({}, {}, {}); diff = {}'.format(ee.x, ee.y, ee.z, diff))
    if diff > 1.0e-16:
        print('PY Test_EQJ_ECL: EXCESSIVE VECTOR ERROR')
        sys.exit(1)

    # Reverse the test: go from ecliptic back to equatorial.
    ir = astronomy.Rotation_ECL_EQJ()
    et = astronomy.RotateVector(ir, ee)
    idiff = VectorDiff(et, ev)
    Debug('PY Test_EQJ_ECL ev diff = {}'.format(idiff))
    if idiff > 2.0e-16:
        print('PY Test_EQJ_ECL: EXCESSIVE REVERSE ROTATION ERROR')
        sys.exit(1)


def Test_EQJ_EQD(body):
    # Verify conversion of equatorial J2000 to equatorial of-date, and back.
    # Use established functions to calculate spherical coordinates for the body, in both EQJ and EQD.
    time = astronomy.Time.Make(2019, 12, 8, 20, 50, 0)
    observer = astronomy.Observer(+35, -85, 0)
    eq2000 = astronomy.Equator(body, time, observer, False, True)
    eqdate = astronomy.Equator(body, time, observer, True, True)

    # Convert EQJ spherical coordinates to vector.
    v2000 = astronomy.VectorFromEquator(eq2000, time)

    # Find rotation matrix.
    r = astronomy.Rotation_EQJ_EQD(time)

    # Rotate EQJ vector to EQD vector.
    vdate = astronomy.RotateVector(r, v2000)

    # Convert vector back to angular equatorial coordinates.
    equcheck = astronomy.EquatorFromVector(vdate)

    # Compare the result with the eqdate.
    ra_diff = vabs(equcheck.ra - eqdate.ra)
    dec_diff = vabs(equcheck.dec - eqdate.dec)
    dist_diff = vabs(equcheck.dist - eqdate.dist)
    Debug('PY Test_EQJ_EQD: {} ra={}, dec={}, dist={}, ra_diff={}, dec_diff={}, dist_diff={}'.format(
        body.name, eqdate.ra, eqdate.dec, eqdate.dist, ra_diff, dec_diff, dist_diff
    ))
    if ra_diff > 1.0e-14 or dec_diff > 1.0e-14 or dist_diff > 4.0e-15:
        print('PY Test_EQJ_EQD: EXCESSIVE ERROR')
        sys.exit(1)

    # Perform the inverse conversion back to equatorial J2000 coordinates.
    ir = astronomy.Rotation_EQD_EQJ(time)
    t2000 = astronomy.RotateVector(ir, vdate)
    diff = VectorDiff(t2000, v2000)
    Debug('PY Test_EQJ_EQD: {} inverse diff = {}'.format(body.name, diff))
    if diff > 5.0e-15:
        print('PY Test_EQJ_EQD: EXCESSIVE INVERSE ERROR')
        sys.exit(1)


def Test_EQD_HOR(body):
    # Use existing functions to calculate horizontal coordinates of the body for the time+observer.
    time = astronomy.Time.Make(1970, 12, 13, 5, 15, 0)
    observer = astronomy.Observer(-37, +45, 0)
    eqd = astronomy.Equator(body, time, observer, True, True)
    Debug('PY Test_EQD_HOR {}: OFDATE ra={}, dec={}'.format(body.name, eqd.ra, eqd.dec))
    hor = astronomy.Horizon(time, observer, eqd.ra, eqd.dec, astronomy.Refraction.Normal)

    # Calculate the position of the body as an equatorial vector of date.
    vec_eqd = astronomy.VectorFromEquator(eqd, time)

    # Calculate rotation matrix to convert equatorial J2000 vector to horizontal vector.
    rot = astronomy.Rotation_EQD_HOR(time, observer)

    # Rotate the equator of date vector to a horizontal vector.
    vec_hor = astronomy.RotateVector(rot, vec_eqd)

    # Convert the horizontal vector to horizontal angular coordinates.
    xsphere = astronomy.HorizonFromVector(vec_hor, astronomy.Refraction.Normal)
    diff_alt = vabs(xsphere.lat - hor.altitude)
    diff_az = vabs(xsphere.lon - hor.azimuth)

    Debug('PY Test_EQD_HOR {}: trusted alt={}, az={}; test alt={}, az={}; diff_alt={}, diff_az={}'.format(
        body.name, hor.altitude, hor.azimuth, xsphere.lat, xsphere.lon, diff_alt, diff_az))

    if diff_alt > 4.0e-14 or diff_az > 1.0e-13:
        print('PY Test_EQD_HOR: EXCESSIVE HORIZONTAL ERROR.')
        sys.exit(1)

    # Confirm that we can convert back to horizontal vector.
    check_hor = astronomy.VectorFromHorizon(xsphere, time, astronomy.Refraction.Normal)
    diff = VectorDiff(check_hor, vec_hor)
    Debug('PY Test_EQD_HOR {}: horizontal recovery: diff = {}'.format(body.name, diff))
    if diff > 3.0e-15:
        print('PY Test_EQD_HOR: EXCESSIVE ERROR IN HORIZONTAL RECOVERY.')
        sys.exit(1)

    # Verify the inverse translation from horizontal vector to equatorial of-date vector.
    irot = astronomy.Rotation_HOR_EQD(time, observer)
    check_eqd = astronomy.RotateVector(irot, vec_hor)
    diff = VectorDiff(check_eqd, vec_eqd)
    Debug('PY Test_EQD_HOR {}: OFDATE inverse rotation diff = {}'.format(body.name, diff))
    if diff > 2.7e-15:
        print('PY Test_EQD_HOR: EXCESSIVE OFDATE INVERSE HORIZONTAL ERROR.')
        sys.exit(1)

    # Exercise HOR to EQJ translation.
    eqj = astronomy.Equator(body, time, observer, False, True)
    vec_eqj = astronomy.VectorFromEquator(eqj, time)
    yrot = astronomy.Rotation_HOR_EQJ(time, observer)
    check_eqj = astronomy.RotateVector(yrot, vec_hor)
    diff = VectorDiff(check_eqj, vec_eqj)
    Debug('PY Test_EQD_HOR {}: J2000 inverse rotation diff = {}'.format(body.name, diff))
    if diff > 5.0e-15:
        print('PY Test_EQD_HOR: EXCESSIVE J2000 INVERSE HORIZONTAL ERROR.')
        sys.exit(1)

    # Verify the inverse translation: EQJ to HOR.
    zrot = astronomy.Rotation_EQJ_HOR(time, observer)
    another_hor = astronomy.RotateVector(zrot, vec_eqj)
    diff = VectorDiff(another_hor, vec_hor)
    Debug('PY Test_EQD_HOR {}: EQJ inverse rotation diff = {}'.format(body.name, diff))
    if diff > 6.0e-15:
        print('PY Test_EQD_HOR: EXCESSIVE EQJ INVERSE HORIZONTAL ERROR.')
        sys.exit(1)

IdentityMatrix = astronomy.RotationMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

def CheckInverse(aname, bname, arot, brot):
    crot = astronomy.CombineRotation(arot, brot)
    caller = 'CheckInverse({},{})'.format(aname, bname)
    CompareMatrices(caller, crot, IdentityMatrix, 2.0e-15)


def CheckCycle(cyclename, arot, brot, crot):
    xrot = astronomy.CombineRotation(arot, brot)
    irot = astronomy.InverseRotation(xrot)
    CompareMatrices(cyclename, crot, irot, 2.0e-15)


def Test_RotRoundTrip():
    # In each round trip, calculate a forward rotation and a backward rotation.
    # Verify the two are inverse matrices.
    time = astronomy.Time.Make(2067, 5, 30, 14, 45, 0)
    observer = astronomy.Observer(+28, -82, 0)

    # Round trip #1: EQJ <==> EQD.
    eqj_eqd = astronomy.Rotation_EQJ_EQD(time)
    eqd_eqj = astronomy.Rotation_EQD_EQJ(time)
    CheckInverse('eqj_eqd', 'eqd_eqj', eqj_eqd, eqd_eqj)

    # Round trip #2: EQJ <==> ECL.
    eqj_ecl = astronomy.Rotation_EQJ_ECL()
    ecl_eqj = astronomy.Rotation_ECL_EQJ()
    CheckInverse('eqj_ecl', 'ecl_eqj', eqj_ecl, ecl_eqj)

    # Round trip #3: EQJ <==> HOR.
    eqj_hor = astronomy.Rotation_EQJ_HOR(time, observer)
    hor_eqj = astronomy.Rotation_HOR_EQJ(time, observer)
    CheckInverse('eqj_hor', 'hor_eqj', eqj_hor, hor_eqj)

    # Round trip #4: EQD <==> HOR.
    eqd_hor = astronomy.Rotation_EQD_HOR(time, observer)
    hor_eqd = astronomy.Rotation_HOR_EQD(time, observer)
    CheckInverse('eqd_hor', 'hor_eqd', eqd_hor, hor_eqd)

    # Round trip #5: EQD <==> ECL.
    eqd_ecl = astronomy.Rotation_EQD_ECL(time)
    ecl_eqd = astronomy.Rotation_ECL_EQD(time)
    CheckInverse('eqd_ecl', 'ecl_eqd', eqd_ecl, ecl_eqd)

    # Round trip #6: HOR <==> ECL.
    hor_ecl = astronomy.Rotation_HOR_ECL(time, observer)
    ecl_hor = astronomy.Rotation_ECL_HOR(time, observer)
    CheckInverse('hor_ecl', 'ecl_hor', hor_ecl, ecl_hor)

    # Verify that combining different sequences of rotations result
    # in the expected combination.
    # For example, (EQJ ==> HOR ==> ECL) must be the same matrix as (EQJ ==> ECL).
    # Each of these is a "triangle" of relationships between 3 orientations.
    # There are 4 possible ways to pick 3 orientations from the 4 to form a triangle.
    # Because we have just proved that each transformation is reversible,
    # we only need to verify the triangle in one cyclic direction.
    CheckCycle('eqj_ecl, ecl_eqd, eqd_eqj', eqj_ecl, ecl_eqd, eqd_eqj)     # excluded corner = HOR
    CheckCycle('eqj_hor, hor_ecl, ecl_eqj', eqj_hor, hor_ecl, ecl_eqj)     # excluded corner = EQD
    CheckCycle('eqj_hor, hor_eqd, eqd_eqj', eqj_hor, hor_eqd, eqd_eqj)     # excluded corner = ECL
    CheckCycle('ecl_eqd, eqd_hor, hor_ecl', ecl_eqd, eqd_hor, hor_ecl)     # excluded corner = EQJ

    Debug('PY Test_RotRoundTrip: PASS')


def Rotation_Pivot():
    tolerance = 1.0e-15

    # Start with an identity matrix.
    ident = astronomy.IdentityMatrix()

    # Pivot 90 degrees counterclockwise around the z-axis.
    r = astronomy.Pivot(ident, 2, +90.0)

    # Put the expected answer in 'a'.
    a = astronomy.RotationMatrix([
        [ 0, +1,  0],
        [-1,  0,  0],
        [ 0,  0, +1],
    ])

    # Compare actual 'r' with expected 'a'.
    CompareMatrices('Rotation_Pivot #1', r, a, tolerance)

    # Pivot again, -30 degrees around the x-axis.
    r = astronomy.Pivot(r, 0, -30.0)

    # Pivot a third time, 180 degrees around the y-axis.
    r = astronomy.Pivot(r, 1, +180.0)

    # Use the 'r' matrix to rotate a vector.
    v1 = astronomy.Vector(1, 2, 3, astronomy.Time(0))
    v2 = astronomy.RotateVector(r, v1)

    # Initialize the expected vector 've'.
    ve = astronomy.Vector(+2.0, +2.3660254037844390, -2.0980762113533156, v1.t)

    CompareVectors('Rotation_Pivot #2', v2, ve, tolerance)

    Debug('PY Rotation_Pivot: PASS')



def Rotation():
    Rotation_MatrixInverse()
    Rotation_MatrixMultiply()
    Rotation_Pivot()
    Test_EQJ_ECL()
    Test_EQJ_EQD(astronomy.Body.Mercury)
    Test_EQJ_EQD(astronomy.Body.Venus)
    Test_EQJ_EQD(astronomy.Body.Mars)
    Test_EQJ_EQD(astronomy.Body.Jupiter)
    Test_EQJ_EQD(astronomy.Body.Saturn)
    Test_EQD_HOR(astronomy.Body.Mercury)
    Test_EQD_HOR(astronomy.Body.Venus)
    Test_EQD_HOR(astronomy.Body.Mars)
    Test_EQD_HOR(astronomy.Body.Jupiter)
    Test_EQD_HOR(astronomy.Body.Saturn)
    Test_RotRoundTrip()
    print('PY Rotation: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def Refraction():
    alt = -90.1
    while alt <= +90.1:
        refr = astronomy.RefractionAngle(astronomy.Refraction.Normal, alt)
        corrected = alt + refr
        inv_refr = astronomy.InverseRefractionAngle(astronomy.Refraction.Normal, corrected)
        check_alt = corrected + inv_refr
        diff = vabs(check_alt - alt)
        if diff > 2.0e-14:
            print('PY Refraction: ERROR - excessive error: alt={}, refr={}, diff={}'.format(alt, refr, diff))
            return 1
        alt += 0.001
    print('PY Refraction: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def PlanetApsis():
    start_time = astronomy.Time.Make(1700, 1, 1, 0, 0, 0)
    body = astronomy.Body.Mercury
    while body.value <= astronomy.Body.Pluto.value:
        count = 1
        period = astronomy._PlanetOrbitalPeriod[body.value]
        filename = os.path.join('apsides', 'apsis_{}.txt'.format(body.value))
        min_interval = -1.0
        max_diff_days = 0.0
        max_dist_ratio = 0.0
        apsis = astronomy.SearchPlanetApsis(body, start_time)
        with open(filename, 'rt') as infile:
            for line in infile:
                token = line.split()
                if len(token) != 3:
                    print('PY PlanetApsis({} line {}): Invalid data format: {} tokens'.format(filename, count, len(token)))
                    return 1
                expected_kind = astronomy.ApsisKind(int(token[0]))
                expected_time = astronomy.Time.Parse(token[1])
                expected_distance = float(token[2])
                if apsis.kind != expected_kind:
                    print('PY PlanetApsis({} line {}): WRONG APSIS KIND: expected {}, found {}'.format(filename, count, expected_kind, apsis.kind))
                    return 1
                diff_days = vabs(expected_time.tt - apsis.time.tt)
                max_diff_days = vmax(max_diff_days, diff_days)
                diff_degrees = (diff_days / period) * 360
                if body == astronomy.Body.Pluto:
                    degree_threshold = 0.262
                else:
                    degree_threshold = 0.1
                if diff_degrees > degree_threshold:
                    print('PY PlanetApsis: FAIL - {} exceeded angular threshold ({} vs {} degrees)'.format(body.name, diff_degrees, degree_threshold))
                    return 1
                diff_dist_ratio = vabs(expected_distance - apsis.dist_au) / expected_distance
                max_dist_ratio = vmax(max_dist_ratio, diff_dist_ratio)
                if diff_dist_ratio > 1.05e-4:
                    print('PY PlanetApsis({} line {}): distance ratio {} is too large.'.format(filename, count, diff_dist_ratio))
                    return 1

                # Calculate the next apsis.
                prev_time = apsis.time
                apsis = astronomy.NextPlanetApsis(body, apsis)
                count += 1
                interval = apsis.time.tt - prev_time.tt
                if min_interval < 0.0:
                    min_interval = max_interval = interval
                else:
                    min_interval = vmin(min_interval, interval)
                    max_interval = vmax(max_interval, interval)
            if count < 2:
                print('PY PlanetApsis: FAILED to find apsides for {}'.format(body))
                return 1
            Debug('PY PlanetApsis: {:4d} apsides for {:<9s} -- intervals: min={:0.2f}, max={:0.2f}, ratio={:0.6f}; max day={:0.3f}, degrees={:0.3f}, dist ratio={:0.6f}'.format(
                count,
                body.name,
                min_interval, max_interval, max_interval / min_interval,
                max_diff_days,
                (max_diff_days / period) * 360.0,
                max_dist_ratio
            ))
        body = astronomy.Body(body.value + 1)
    print('PY PlanetApsis: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def Constellation():
    inFileName = 'constellation/test_input.txt'
    lnum = 0
    failcount = 0
    with open(inFileName, 'rt') as infile:
        for line in infile:
            lnum += 1
            m = re.match(r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+([A-Z][a-zA-Z]{2})\s*$', line)
            if not m:
                print('PY Constellation: invalid line {} in file {}'.format(lnum, inFileName))
                return 1
            id = int(m.group(1))
            ra = float(m.group(2))
            dec = float(m.group(3))
            symbol = m.group(4)
            constel = astronomy.Constellation(ra, dec)
            if constel.symbol != symbol:
                print('Star {:6d}: expected {}, found {} at B1875 RA={:10.6f}, DEC={:10.6f}'.format(id, symbol, constel.symbol, constel.ra1875, constel.dec1875))
                failcount += 1
    if failcount > 0:
        print('PY Constellation: {} failures'.format(failcount))
        return 1
    print('PY Constellation: PASS (verified {})'.format(lnum))
    return 0

#-----------------------------------------------------------------------------------------------------------

def LunarEclipseIssue78():
    # https://github.com/cosinekitty/astronomy/issues/78

    eclipse = astronomy.SearchLunarEclipse(astronomy.Time.Make(2020, 12, 19, 0, 0, 0))
    expected_peak = astronomy.Time.Make(2021, 5, 26, 11, 18, 42)  # https://www.timeanddate.com/eclipse/lunar/2021-may-26
    dt = (expected_peak.tt - eclipse.peak.tt) * (24.0 * 3600.0)
    if vabs(dt) > 40.0:
        print('LunarEclipseIssue78: Excessive prediction error = {} seconds.'.format(dt))
        return 1
    if eclipse.kind != astronomy.EclipseKind.Total:
        print('Expected total eclipse; found: {}'.format(eclipse.kind))
        return 1
    print('PY LunarEclipseIssue78: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def LunarEclipse():
    astronomy._CalcMoonCount = 0
    filename = 'eclipse/lunar_eclipse.txt'
    with open(filename, 'rt') as infile:
        eclipse = astronomy.SearchLunarEclipse(astronomy.Time.Make(1701, 1, 1, 0, 0, 0))
        lnum = 0
        skip_count = 0
        diff_count = 0
        sum_diff_minutes = 0.0
        max_diff_minutes = 0.0
        diff_limit = 2.0
        for line in infile:
            lnum += 1
            if len(line) < 17:
                print('PY LunarEclipse({} line {}): line is too short.'.format(filename, lnum))
                return 1
            time_text = line[0:17]
            peak_time = astronomy.Time.Parse(time_text)
            token = line[17:].split()
            if len(token) != 2:
                print('PY LunarEclipse({} line {}): wrong number of tokens.'.format(filename, lnum))
                return 1
            partial_minutes = float(token[0])
            total_minutes = float(token[1])
            valid = False
            # Verify that the calculated eclipse semi-durations are consistent with the kind.
            if eclipse.kind == astronomy.EclipseKind.Penumbral:
                valid = (eclipse.sd_penum > 0.0) and (eclipse.sd_partial == 0.0) and (eclipse.sd_total == 0.0)
            elif eclipse.kind == astronomy.EclipseKind.Partial:
                valid = (eclipse.sd_penum > 0.0) and (eclipse.sd_partial > 0.0) and (eclipse.sd_total == 0.0)
            elif eclipse.kind == astronomy.EclipseKind.Total:
                valid = (eclipse.sd_penum > 0.0) and (eclipse.sd_partial > 0.0) and (eclipse.sd_total > 0.0)
            else:
                print('PY LunarEclipse({} line {}): invalid eclipse kind {}.'.format(filename, lnum, eclipse.kind))
                return 1

            if not valid:
                print('PY LunarEclipse({} line {}): invalid semidurations.'.format(filename, lnum))
                return 1

            # Check eclipse peak time.
            diff_days = eclipse.peak.ut - peak_time.ut

            # Tolerate missing penumbral eclipses - skip to next input line without calculating next eclipse.
            if partial_minutes == 0.0 and diff_days > 20.0:
                skip_count += 1
                continue

            diff_minutes = (24.0 * 60.0) * vabs(diff_days)
            sum_diff_minutes += diff_minutes
            diff_count += 1

            if diff_minutes > diff_limit:
                print("PY LunarEclipse expected center: {}".format(peak_time))
                print("PY LunarEclipse found    center: {}".format(eclipse.peak))
                print("PY LunarEclipse({} line {}): EXCESSIVE center time error = {} minutes ({} days).".format(filename, lnum, diff_minutes, diff_days))
                return 1

            if diff_minutes > max_diff_minutes:
                max_diff_minutes = diff_minutes

            # check partial eclipse duration

            diff_minutes = vabs(partial_minutes - eclipse.sd_partial)
            sum_diff_minutes += diff_minutes
            diff_count += 1

            if diff_minutes > diff_limit:
                print("PY LunarEclipse({} line {}): EXCESSIVE partial eclipse semiduration error: {} minutes".format(filename, lnum, diff_minutes))
                return 1

            if diff_minutes > max_diff_minutes:
                max_diff_minutes = diff_minutes

            # check total eclipse duration

            diff_minutes = vabs(total_minutes - eclipse.sd_total)
            sum_diff_minutes += diff_minutes
            diff_count += 1

            if diff_minutes > diff_limit:
                print("PY LunarEclipse({} line {}): EXCESSIVE total eclipse semiduration error: {} minutes".format(filename, lnum, diff_minutes))
                return 1

            if diff_minutes > max_diff_minutes:
                max_diff_minutes = diff_minutes

            # calculate for next iteration

            eclipse = astronomy.NextLunarEclipse(eclipse.peak)
    print("PY LunarEclipse: PASS (verified {}, skipped {}, max_diff_minutes = {}, avg_diff_minutes = {}, moon calcs = {})".format(lnum, skip_count, max_diff_minutes, (sum_diff_minutes / diff_count), astronomy._CalcMoonCount))
    return 0

#-----------------------------------------------------------------------------------------------------------

def VectorFromAngles(lat, lon):
    rlat = math.radians(v(lat))
    rlon = math.radians(v(lon))
    coslat = math.cos(rlat)
    return [
        math.cos(rlon) * coslat,
        math.sin(rlon) * coslat,
        math.sin(rlat)
    ]


def AngleDiff(alat, alon, blat, blon):
    a = VectorFromAngles(alat, alon)
    b = VectorFromAngles(blat, blon)
    dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    if dot <= -1.0:
        return 180.0
    if dot >= +1.0:
        return 0.0
    return v(math.degrees(math.acos(dot)))


def KindFromChar(typeChar):
    return {
        'P': astronomy.EclipseKind.Partial,
        'A': astronomy.EclipseKind.Annular,
        'T': astronomy.EclipseKind.Total,
        'H': astronomy.EclipseKind.Total,
    }[typeChar]

def GlobalSolarEclipse():
    expected_count = 1180
    max_minutes = 0.0
    max_angle = 0.0
    skip_count = 0
    eclipse = astronomy.SearchGlobalSolarEclipse(astronomy.Time.Make(1701, 1, 1, 0, 0, 0))
    filename = 'eclipse/solar_eclipse.txt'
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            # 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
            token = line.split()
            if len(token) != 5:
                print('PY GlobalSolarEclipse({} line {}): invalid token count = {}'.format(filename, lnum, len(token)))
                return 1
            peak = astronomy.Time.Parse(token[0])
            expected_kind = KindFromChar(token[2])
            lat = float(token[3])
            lon = float(token[4])

            diff_days = eclipse.peak.tt - peak.tt
            # Sometimes we find marginal eclipses that aren't listed in the test data.
            # Ignore them if the distance between the Sun/Moon shadow axis and the Earth's center is large.
            while diff_days < -25.0 and eclipse.distance > 9000.0:
                skip_count += 1
                eclipse = astronomy.NextGlobalSolarEclipse(eclipse.peak)
                diff_days = eclipse.peak.ut - peak.ut

            # Validate the eclipse prediction.
            diff_minutes = (24 * 60) * vabs(diff_days)
            if diff_minutes > 6.93:
                print('PY GlobalSolarEclipse({} line {}): EXCESSIVE TIME ERROR = {} minutes'.format(filename, lnum, diff_minutes))
                return 1

            if diff_minutes > max_minutes:
                max_minutes = diff_minutes

            # Validate the eclipse kind, but only when it is not a "glancing" eclipse.
            if (eclipse.distance < 6360) and (eclipse.kind != expected_kind):
                print('PY GlobalSolarEclipse({} line {}): WRONG ECLIPSE KIND: expected {}, found {}'.format(filename, lnum, expected_kind, eclipse.kind))
                return 1

            if eclipse.kind == astronomy.EclipseKind.Total or eclipse.kind == astronomy.EclipseKind.Annular:
                # When the distance between the Moon's shadow ray and the Earth's center is beyond 6100 km,
                # it creates a glancing blow whose geographic coordinates are excessively sensitive to
                # slight changes in the ray. Therefore, it is unreasonable to count large errors there.
                if eclipse.distance < 6100.0:
                    diff_angle = AngleDiff(lat, lon, eclipse.latitude, eclipse.longitude)
                    if diff_angle > 0.247:
                        print('PY GlobalSolarEclipse({} line {}): EXCESSIVE GEOGRAPHIC LOCATION ERROR = {} degrees'.format(filename, lnum, diff_angle))
                        return 1
                    if diff_angle > max_angle:
                        max_angle = diff_angle

            eclipse = astronomy.NextGlobalSolarEclipse(eclipse.peak)

    if lnum != expected_count:
        print('PY GlobalSolarEclipse: WRONG LINE COUNT = {}, expected {}'.format(lnum, expected_count))
        return 1

    if skip_count > 2:
        print('PY GlobalSolarEclipse: EXCESSIVE SKIP COUNT = {}'.format(skip_count))
        return 1

    print('PY GlobalSolarEclipse: PASS ({} verified, {} skipped, max minutes = {}, max angle = {})'.format(lnum, skip_count, max_minutes, max_angle))
    return 0

#-----------------------------------------------------------------------------------------------------------

def LocalSolarEclipse1():
    expected_count = 1180
    max_minutes = 0.0
    skip_count = 0
    filename = 'eclipse/solar_eclipse.txt'
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            # 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
            token = line.split()
            if len(token) != 5:
                print('PY LocalSolarEclipse1({} line {}): invalid token count = {}'.format(filename, lnum, len(token)))
                return 1
            peak = astronomy.Time.Parse(token[0])
            #typeChar = token[2]
            lat = float(token[3])
            lon = float(token[4])
            observer = astronomy.Observer(lat, lon, 0.0)

            # Start the search 20 days before we know the eclipse should peak.
            search_start = peak.AddDays(-20)
            eclipse = astronomy.SearchLocalSolarEclipse(search_start, observer)

            # Validate the predicted eclipse peak time.
            diff_days = eclipse.peak.time.tt - peak.tt
            if diff_days > 20:
                skip_count += 1
                continue

            diff_minutes = (24 * 60) * vabs(diff_days)
            if diff_minutes > 7.14:
                print('PY LocalSolarEclipse1({} line {}): EXCESSIVE TIME ERROR = {} minutes'.format(filename, lnum, diff_minutes))
                return 1

            if diff_minutes > max_minutes:
                max_minutes = diff_minutes

    if lnum != expected_count:
        print('PY LocalSolarEclipse1: WRONG LINE COUNT = {}, expected {}'.format(lnum, expected_count))
        return 1

    if skip_count > 6:
        print('PY LocalSolarEclipse1: EXCESSIVE SKIP COUNT = {}'.format(skip_count))
        return 1

    print('PY LocalSolarEclipse1: PASS ({} verified, {} skipped, max minutes = {})'.format(lnum, skip_count, max_minutes))
    return 0


def TrimLine(line):
    # Treat '#' as a comment character.
    poundIndex = line.find('#')
    if poundIndex >= 0:
        line = line[:poundIndex]
    return line.strip()


def ParseEvent(time_str, alt_str, required):
    if required:
        time = astronomy.Time.Parse(time_str)
        altitude = float(alt_str)
        return astronomy.EclipseEvent(time, altitude)
    if time_str != '-':
        raise Exception('Expected event time to be "-" but found "{}"'.format(time_str))
    return None


def LocalSolarEclipse2():
    # Test ability to calculate local solar eclipse conditions away from
    # the peak position on the Earth.

    filename = 'eclipse/local_solar_eclipse.txt'
    lnum = 0
    verify_count = 0
    max_minutes = 0.0
    max_degrees = 0.0

    def CheckEvent(calc, expect):
        nonlocal max_minutes, max_degrees
        diff_minutes = (24 * 60) * vabs(expect.time.ut - calc.time.ut)
        if diff_minutes > max_minutes:
            max_minutes = diff_minutes
        if diff_minutes > 1.0:
            raise Exception('CheckEvent({} line {}): EXCESSIVE TIME ERROR: {} minutes.'.format(filename, lnum, diff_minutes))
        diff_alt = vabs(expect.altitude - calc.altitude)
        if diff_alt > max_degrees:
            max_degrees = diff_alt
        if diff_alt > 0.5:
            raise Exception('CheckEvent({} line {}): EXCESSIVE ALTITUDE ERROR: {} degrees.'.format(filename, lnum, diff_alt))

    with open(filename, 'rt') as infile:
        for line in infile:
            lnum += 1
            line = TrimLine(line)
            if line == '':
                continue
            token = line.split()
            if len(token) != 13:
                print('PY LocalSolarEclipse2({} line {}): Incorrect token count = {}'.format(filename, lnum, len(token)))
                return 1
            latitude = float(token[0])
            longitude = float(token[1])
            observer = astronomy.Observer(latitude, longitude, 0)
            expected_kind = KindFromChar(token[2])
            is_umbral = (expected_kind != astronomy.EclipseKind.Partial)
            p1    = ParseEvent(token[3],  token[4],   True)
            t1    = ParseEvent(token[5],  token[6],   is_umbral)
            peak  = ParseEvent(token[7],  token[8],   True)
            t2    = ParseEvent(token[9],  token[10],  is_umbral)
            p2    = ParseEvent(token[11], token[12],  True)
            search_time = p1.time.AddDays(-20)
            eclipse = astronomy.SearchLocalSolarEclipse(search_time, observer)
            if eclipse.kind != expected_kind:
                print('PY LocalSolarEclipse2({} line {}): expected eclipse kind "{}" but found "{}".'.format(
                    filename, lnum, expected_kind, eclipse.kind
                ))
                return 1
            CheckEvent(eclipse.peak, peak)
            CheckEvent(eclipse.partial_begin, p1)
            CheckEvent(eclipse.partial_end, p2)
            if is_umbral:
                CheckEvent(eclipse.total_begin, t1)
                CheckEvent(eclipse.total_end, t2)
            verify_count += 1
    print('PY LocalSolarEclipse2: PASS ({} verified, max_minutes = {}, max_degrees = {})'.format(verify_count, max_minutes, max_degrees))
    return 0


def LocalSolarEclipse():
    if 0 != LocalSolarEclipse1():
        return 1
    if 0 != LocalSolarEclipse2():
        return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def TransitFile(body, filename, limit_minutes, limit_sep):
    lnum = 0
    max_minutes = 0
    max_sep = 0
    with open(filename, 'rt') as infile:
        transit = astronomy.SearchTransit(body, astronomy.Time.Make(1600, 1, 1, 0, 0, 0))
        for line in infile:
            lnum += 1
            token = line.strip().split()
            # 22:17 1881-11-08T00:57Z 03:38  3.8633
            if len(token) != 4:
                print('PY TransitFile({} line {}): bad data format.'.format(filename, lnum))
                return 1

            textp = token[1]
            text1 = textp[0:11] + token[0] + 'Z'
            text2 = textp[0:11] + token[2] + 'Z'
            timep = astronomy.Time.Parse(textp)
            time1 = astronomy.Time.Parse(text1)
            time2 = astronomy.Time.Parse(text2)
            separation = float(token[3])

            # If the start time is after the peak time, it really starts on the previous day.
            if time1.ut > timep.ut:
                time1 = time1.AddDays(-1.0)

            # If the finish time is before the peak time, it really starts on the following day.
            if time2.ut < timep.ut:
                time2 = time2.AddDays(+1.0)

            diff_start  = (24.0 * 60.0) * vabs(time1.ut - transit.start.ut )
            diff_peak   = (24.0 * 60.0) * vabs(timep.ut - transit.peak.ut  )
            diff_finish = (24.0 * 60.0) * vabs(time2.ut - transit.finish.ut)
            diff_sep = vabs(separation - transit.separation)
            max_minutes = vmax(max_minutes, diff_start)
            max_minutes = vmax(max_minutes, diff_peak)
            max_minutes = vmax(max_minutes, diff_finish)
            if max_minutes > limit_minutes:
                print('PY TransitFile({} line {}): EXCESSIVE TIME ERROR = {} minutes.'.format(filename, lnum, max_minutes))
                return 1
            max_sep = vmax(max_sep, diff_sep)
            if max_sep > limit_sep:
                print('PY TransitFile({} line {}): EXCESSIVE SEPARATION ERROR = {} arcminutes.'.format(filename, lnum, max_sep))
                return 1
            transit = astronomy.NextTransit(body, transit.finish)
    print('PY TransitFile({}): PASS - verified {}, max minutes = {}, max sep arcmin = {}'.format(filename, lnum, max_minutes, max_sep))
    return 0


def Transit():
    if 0 != TransitFile(astronomy.Body.Mercury, 'eclipse/mercury.txt', 10.710, 0.2121):
        return 1
    if 0 != TransitFile(astronomy.Body.Venus, 'eclipse/venus.txt', 9.109, 0.6772):
        return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def PlutoCheckDate(ut, arcmin_tolerance, x, y, z):
    time = astronomy.Time(ut)
    try:
        timeText = str(time)
    except OverflowError:
        timeText = "???"
    Debug('PY PlutoCheck: {} = {} UT = {} TT'.format(timeText, time.ut, time.tt))
    vector = astronomy.HelioVector(astronomy.Body.Pluto, time)
    dx = v(vector.x - x)
    dy = v(vector.y - y)
    dz = v(vector.z - z)
    diff = sqrt(dx*dx + dy*dy + dz*dz)
    dist = sqrt(x*x + y*y + z*z) - 1.0
    arcmin = (diff / dist) * (180.0 * 60.0 / math.pi)
    Debug('PY PlutoCheck: calc pos = [{}, {}, {}]'.format(vector.x, vector.y, vector.z))
    Debug('PY PlutoCheck: ref  pos = [{}, {}, {}]'.format(x, y, z))
    Debug('PY PlutoCheck: del  pos = [{}, {}, {}]'.format(vector.x - x, vector.y - y, vector.z - z))
    Debug('PY PlutoCheck: diff = {} AU, {} arcmin'.format(diff, arcmin))
    if v(arcmin) > arcmin_tolerance:
        print('PY PlutoCheck: EXCESSIVE ERROR')
        return 1
    Debug('')
    return 0


def PlutoCheck():
    if PlutoCheckDate(  +18250.0,  0.271, +37.4377303523676090, -10.2466292454075898, -14.4773101310875809): return 1
    if PlutoCheckDate(  +18250.0,  0.271, +37.4377303523676090, -10.2466292454075898, -14.4773101310875809): return 1
    if PlutoCheckDate( -856493.0,  6.636, +23.4292113199166252, +42.1452685817740829,  +6.0580908436642940): return 1
    if PlutoCheckDate( +435633.0,  0.058, -27.3178902095231813, +18.5887022581070305, +14.0493896259306936): return 1
    if PlutoCheckDate(       0.0, 4.0e-9, -9.8753673425269000,  -27.9789270580402771,  -5.7537127596369588): return 1
    if PlutoCheckDate( +800916.0,  6.705, -29.5266052645301365, +12.0554287322176474, +12.6878484911631091): return 1
    print("PY PlutoCheck: PASS")
    return 0

#-----------------------------------------------------------------------------------------------------------

UnitTests = {
    'constellation':            Constellation,
    'elongation':               Elongation,
    'global_solar_eclipse':     GlobalSolarEclipse,
    'local_solar_eclipse':      LocalSolarEclipse,
    'lunar_apsis':              LunarApsis,
    'lunar_eclipse':            LunarEclipse,
    'lunar_eclipse_78':         LunarEclipseIssue78,
    'magnitude':                Magnitude,
    'moon':                     GeoMoon,
    'moonphase':                MoonPhase,
    'planet_apsis':             PlanetApsis,
    'pluto':                    PlutoCheck,
    'refraction':               Refraction,
    'riseset':                  RiseSet,
    'rotation':                 Rotation,
    'seasons':                  Seasons,
    'time':                     AstroTime,
    'transit':                  Transit,
}

#-----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-v':
        sys.argv = sys.argv[1:]
        Verbose = True

    if len(sys.argv) == 2:
        name = sys.argv[1]
        if name in UnitTests:
            sys.exit(UnitTests[name]())
        if name in ['astro_check', 'astro_profile']:
            sys.exit(AstroCheck(sys.argv[1] == 'astro_check'))
        if name == 'all':
            for name in sorted(UnitTests.keys()):
                func = UnitTests[name]
                Debug('test.py: Starting test "{}"'.format(name))
                rc = func()
                Debug('test.py: Test "{}" returned {}'.format(name, rc))
                if rc != 0:
                    sys.exit(1)
            print('test.py: ALL PASS')
            sys.exit(0)

    print('test.py: Invalid command line arguments.')
    sys.exit(1)
