#!/usr/bin/env python3
import sys
import math
import re
import os
from itertools import chain
sys.path.append('../source/python')
import astronomy
from dontrig import xcos, xsin

#-----------------------------------------------------------------------------------------------------------

if os.getenv('GITHUB_JOB'):
    # Super weird hack: re-define _sin and _cos to have exactly reproducible values.
    print('test.py - NOTE: Detected GitHub job. Replacing trig functions...')
    astronomy._sin = xsin
    astronomy._cos = xcos

#-----------------------------------------------------------------------------------------------------------

Verbose = False
SECONDS_PER_DAY = 86400.0
MINUTES_PER_DAY = 1440.0

def Debug(text):
    if Verbose:
        print(text)

def Pass(funcname):
    print('PY {}: PASS'.format(funcname))
    return 0

def Fail(funcname, reason):
    print('PY {} FAIL: {}'.format(funcname, reason))
    return 1

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

def CalendarCase(year, month, day, hour, minute, second):
    # Convert to Astronomy Engine Time object.
    time = astronomy.Time.Make(year, month, day, hour, minute, second)
    # Convert to back calendar date tuple.
    (cyear, cmonth, cday, chour, cminute, csecond) = time.Calendar()
    if (cyear, cmonth, cday) != (year, month, day):
        return Fail('CalendarCase', 'Expected {:06d}-{:02d}-{:02d} but found {:06d}-{:02d}-{:02d}'.format(
            year, month, day,
            cyear, cmonth, cday
        ))
    expectedMillis = 1000.0*(second + 60.0*(minute + 60.0*hour))
    calcMillis = 1000.0*(csecond + 60.0*(cminute + 60.0*chour))
    diffMillis = vabs(calcMillis - expectedMillis)
    if diffMillis > 4.0:
        return Fail('CalendarCase', 'EXCESSIVE millisecond error = {:0.6f} for {:06d}-{:02d}-{:02d}'.format(
            diffMillis, year, month, day
        ))
    return 0

def AstroTime():
    expected_ut = 6910.270978506945
    expected_tt = 6910.271800214368
    time = astronomy.Time.Make(2018, 12, 2, 18, 30, 12.543)
    diff = time.ut - expected_ut
    if vabs(diff) > 1.0e-12:
        print('PY AstroTime: excessive UT error {}'.format(diff))
        return 1
    diff = time.tt - expected_tt
    if vabs(diff) > 1.0e-12:
        print('PY AstroTime: excessive TT error {}'.format(diff))
        return 1
    s = str(time.Utc())
    if s != '2018-12-02 18:30:12.543000':
        print('PY AstroTime: Utc() returned incorrect string "{}"'.format(s))
        return 1
    time = astronomy.Time.Make(2018, 12, 31, 23, 59, 59.9999)
    expected = '2018-12-31T23:59:59.999Z'
    s = str(time)
    if s != expected:
        print('PY AstroTime: expected {} but found {}'.format(expected, s))
        return 1
    print('PY Current time =', astronomy.Time.Now())
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
    # Extreme year values...
    AssertGoodTime('-4172-12-02T14:30:45.123Z', '-004172-12-02T14:30:45.123Z')
    AssertGoodTime('-4173-12-02T14:30:45.123Z', '-004173-12-02T14:30:45.123Z')
    AssertGoodTime('-4174-12-02T14:30:45.123Z', '-004174-12-02T14:30:45.123Z')
    AssertGoodTime('-4175-12-02T14:30:45.123Z', '-004175-12-02T14:30:45.123Z')
    AssertGoodTime('-4176-12-02T14:30:45.123Z', '-004176-12-02T14:30:45.123Z')
    AssertGoodTime('-2300-12-19T16:22:26.325Z', '-002300-12-19T16:22:26.325Z')
    AssertGoodTime('-2300-12-19T16:22:26.325Z', '-002300-12-19T16:22:26.325Z')
    AssertGoodTime('+12345-12-11T13:30:10.041Z', '+012345-12-11T13:30:10.040Z')
    AssertGoodTime('+12346-12-11T13:30:10.041Z', '+012346-12-11T13:30:10.040Z')
    AssertGoodTime('+12347-12-11T13:30:10.041Z', '+012347-12-11T13:30:10.040Z')
    AssertGoodTime('+12348-12-11T13:30:10.041Z', '+012348-12-11T13:30:10.040Z')
    AssertGoodTime('-123456-01-14T22:55:12.000Z', '-123456-01-14T22:55:11.999Z')
    AssertGoodTime('+123456-01-14T22:55:12.000Z', '+123456-01-14T22:55:11.999Z')
    AssertGoodTime('-999995-01-14T22:55:12.297Z', '-999995-01-14T22:55:12.297Z')
    AssertGoodTime('-999996-01-14T22:55:12.297Z', '-999996-01-14T22:55:12.297Z')
    AssertGoodTime('-999997-01-14T22:55:12.297Z', '-999997-01-14T22:55:12.297Z')
    AssertGoodTime('-999998-01-14T22:55:12.297Z', '-999998-01-14T22:55:12.297Z')
    AssertGoodTime('-999999-01-14T22:55:12.000Z', '-999999-01-14T22:55:11.998Z')
    AssertGoodTime('+999999-01-14T22:55:12.000Z', '+999999-01-14T22:55:11.998Z')

    nyears = 0
    for year in chain(range(-999999, -995999), range(-3000, 3001), range(+996000, +1000000)):
        # Check just before and after each potential leap day.
        if CalendarCase(year, 2, 28, 14, 45, 28.321):
            return 1
        if CalendarCase(year, 3,  1, 14, 45, 28.321):
            return 1
        nyears += 1
    return Pass('AstroTime({} calendar years)'.format(nyears))

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

def SelectJupiterMoon(jm, mindex):
    return [jm.io, jm.europa, jm.ganymede, jm.callisto][mindex]

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
                    print('v {} {:0.18e} {:0.18e} {:0.18e} {:0.18e}'.format(name, pos.t.tt, pos.x, pos.y, pos.z))
                if body != astronomy.Body.Earth and body != astronomy.Body.EMB and body != astronomy.Body.SSB:
                    j2000 = astronomy.Equator(body, time, observer, False, False)
                    ofdate = astronomy.Equator(body, time, observer, True, True)
                    hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.Refraction.Airless)
                    if printflag:
                        print('s {} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e}'.format(name, time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude))
        pos = astronomy.GeoMoon(time)
        if printflag:
            print('v GM {:0.18e} {:0.18e} {:0.18e} {:0.18e}'.format(pos.t.tt, pos.x, pos.y, pos.z))
        j2000 = astronomy.Equator(astronomy.Body.Moon, time, observer, False, False)
        ofdate = astronomy.Equator(astronomy.Body.Moon, time, observer, True, True)
        hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.Refraction.Airless)
        if printflag:
            print('s GM {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e}'.format(time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude))
        jm = astronomy.JupiterMoons(time)
        if printflag:
            for mindex in range(4):
                moon = SelectJupiterMoon(jm, mindex)
                print('j {:d} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e} {:0.18e}'.format(mindex, time.tt, time.ut, moon.x, moon.y, moon.z, moon.vx, moon.vy, moon.vz))
        if printflag:
            # Nutation calculations
            print('n {:0.18e} {:0.18e}'.format(time._et.dpsi, time._et.deps))
        sphere = astronomy.EclipticGeoMoon(time)
        if printflag:
            print('m {:0.18f} {:0.18f} {:0.18f}'.format(sphere.lat, sphere.lon, sphere.dist))
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
                current_year = year
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


def SeasonsIssue187():
    # This is a regression test for:
    # https://github.com/cosinekitty/astronomy/issues/187
    # For years far from the present, the seasons search was sometimes failing.
    for year in range(1, 9999, 1):
        try:
            astronomy.Seasons(year)
        except astronomy.InternalError:
            print('PY SeasonsIssue187: FAIL - internal error for year {}'.format(year))
            return 1
    print('PY SeasonsIssue187: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def MoonPhase(filename = 'moonphase/moonphases.txt'):
    threshold_seconds = 90.0       # max tolerable prediction error in seconds
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
            diff_seconds = vabs(mq.time.tt - expected_time.tt) * SECONDS_PER_DAY
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
            if vabs(diff_minutes) > 6.8:
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
    if hour_diff > 0.6:
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
    return (
        TestElongFile('longitude/opposition_2018.txt', 0.0) or
        TestPlanetLongitudes(astronomy.Body.Mercury, "temp/py_longitude_Mercury.txt", "inf") or
        TestPlanetLongitudes(astronomy.Body.Venus,   "temp/py_longitude_Venus.txt",   "inf") or
        TestPlanetLongitudes(astronomy.Body.Mars,    "temp/py_longitude_Mars.txt",    "opp") or
        TestPlanetLongitudes(astronomy.Body.Jupiter, "temp/py_longitude_Jupiter.txt", "opp") or
        TestPlanetLongitudes(astronomy.Body.Saturn,  "temp/py_longitude_Saturn.txt",  "opp") or
        TestPlanetLongitudes(astronomy.Body.Uranus,  "temp/py_longitude_Uranus.txt",  "opp") or
        TestPlanetLongitudes(astronomy.Body.Neptune, "temp/py_longitude_Neptune.txt", "opp") or
        TestPlanetLongitudes(astronomy.Body.Pluto,   "temp/py_longitude_Pluto.txt",   "opp") or
        SearchElongTest() or
        Pass('Elongation')
    )

#-----------------------------------------------------------------------------------------------------------

def MonthNumber(mtext):
    return 1 + ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].index(mtext)

def ParseJplHorizonsDateTime(line):
    m = re.match(r'^\s*(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s(\d{2}):(\d{2})\s+(.*)$', line)
    if not m:
        return None, None
    year = int(m.group(1))
    month = MonthNumber(m.group(2))
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
        ( "1972-01-01T00:00Z", -0.31725492,  +24.43386475 ),
        ( "1980-01-01T00:00Z", +0.85796177,   -1.72627324 ),
        ( "2009-09-04T00:00Z", +1.01932560,   +0.01834451 ),
        ( "2017-06-15T00:00Z", -0.12303373,  -26.60068380 ),
        ( "2019-05-01T00:00Z", +0.33124502,  -23.47173574 ),
        ( "2025-09-25T00:00Z", +0.50543708,   +1.69118986 ),
        ( "2032-05-15T00:00Z", -0.04649573,  +26.95238680 )
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

def ToggleDir(dir):
    return astronomy.Direction(-dir.value)

def RiseSetSlot(ut1, ut2, direction, observer):
    maxDiff = 0.0
    nslots = 100
    for i in range(1, nslots):
        ut = ut1 + (i / nslots)*(ut2 - ut1)
        time = astronomy.Time(ut)
        result = astronomy.SearchRiseSet(astronomy.Body.Sun, observer, direction, time, -1.0)
        if not result:
            print('PY RiseSetSlot: backward slot search failed for {} before {}'.format(direction, time))
            return 1
        diff = SECONDS_PER_DAY * vabs(result.ut - ut1)
        maxDiff = max(maxDiff, diff)
        result = astronomy.SearchRiseSet(astronomy.Body.Sun, observer, direction, time, +1.0)
        if not result:
            print('PY RiseSetSlot: forward slot search failed for {} after {}'.format(direction, time))
            return 1
        diff = SECONDS_PER_DAY * vabs(result.ut - ut2)
        maxDiff = max(maxDiff, diff)

    if maxDiff > 0.13:
        print('PY RiseSetSlot: EXCESSIVE {} slot-test discrepancy = {:0.6f} seconds.'.format(direction, maxDiff))
        return 1
    Debug('PY RiseSetSlot: {} slot-test discrepancy = {:0.6f} seconds.'.format(direction, maxDiff))
    return 0


def RiseSetReverse():
    nsamples = 5000
    nudge = 0.1
    utList = []
    observer = astronomy.Observer(30.5, -90.7, 0.0)
    dtMin = +1000.0
    dtMax = -1000.0
    maxDiff = 0.0

    # Find alternating sunrise/sunset events in forward chronological order.
    dir = astronomy.Direction.Rise
    time = astronomy.Time.Make(2022, 1, 1, 0, 0, 0)
    for i in range(nsamples):
        result = astronomy.SearchRiseSet(astronomy.Body.Sun, observer, dir, time, +1.0)
        if not result:
            print('PY RiseSetReverse: cannot find {} event after {}'.format(dir, time))
            return 1
        utList.append(result.ut)
        if i > 0:
            # Check the time between consecutive sunrise/sunset events.
            # These will vary considerably with the seasons, so just make sure we don't miss any entirely.
            dt = v(utList[i] - utList[i-1])
            dtMin = min(dtMin, dt)
            dtMax = max(dtMax, dt)
        dir = ToggleDir(dir)
        time = result.AddDays(+nudge)

    Debug('PY RiseSetReverse: dtMin={:0.6f} days, dtMax={:0.6f} days.'.format(dtMin, dtMax))
    if (dtMin < 0.411) or (dtMax > 0.589):
        print('PY RiseSetReverse: Invalid intervals between sunrise/sunset.')
        return 1

    # Perform the same search in reverse. Verify we get consistent rise/set times.
    for i in range(nsamples-1, -1, -1):
        dir = ToggleDir(dir)
        result = astronomy.SearchRiseSet(astronomy.Body.Sun, observer, dir, time, -1.0)
        if not result:
            print('PY RiseSetReverse: cannot find {] event before {}.'.format(dir, time))
            return 1
        diff = SECONDS_PER_DAY * vabs(utList[i] - result.ut)
        maxDiff = max(maxDiff, diff)
        time = result.AddDays(-nudge)

    if maxDiff > 0.1:
        print('PY RiseSetReverse: EXCESSIVE forward/backward discrepancy = {:0.6f} seconds.'.format(maxDiff))
        return 1
    Debug('PY RiseSetReverse: forward/backward discrepancy = {:0.6f} seconds.'.format(maxDiff))

    # All even indexes in utList hold sunrise times.
    # All odd indexes in utList hold sunset times.
    # Verify that forward/backward searches for consecutive sunrises/sunsets
    # resolve correctly for 100 time slots between them.
    k = (nsamples // 2) & ~1

    return (
        RiseSetSlot(utList[k], utList[k+2], astronomy.Direction.Rise, observer) or
        RiseSetSlot(utList[k+1], utList[k+3], astronomy.Direction.Set, observer) or
        Pass('RiseSetReverse')
    )

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
                print('PY RiseSet({} line {}): expected dir={} but found {}'.format(filename, lnum, direction, a_dir))
                return 1

            error_minutes = (24.0 * 60.0) * vabs(a_evt.tt - correct_time.tt)
            sum_minutes += error_minutes ** 2
            max_minutes = vmax(max_minutes, error_minutes)
            if error_minutes > 1.18:
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


def Test_GAL_EQJ_NOVAS(filename):
    THRESHOLD_SECONDS = 8.8
    rot = astronomy.Rotation_EQJ_GAL()
    time = astronomy.Time(0.0)      # placeholder time - value does not matter
    with open(filename, 'rt') as infile:
        lnum = 0
        max_diff = 0.0
        for line in infile:
            lnum += 1
            token = line.split()
            if len(token) != 4:
                print('PY Test_GAL_EQJ_NOVAS({} line {}): Wrong number of tokens.'.format(filename, lnum))
                sys.exit(1)
            ra = float(token[0])
            dec = float(token[1])
            glon = float(token[2])
            glat = float(token[3])
            eqj_sphere = astronomy.Spherical(dec, 15.0*ra, 1.0)
            eqj_vec = astronomy.VectorFromSphere(eqj_sphere, time)
            gal_vec = astronomy.RotateVector(rot, eqj_vec)
            gal_sphere = astronomy.SphereFromVector(gal_vec)
            dlat = gal_sphere.lat - glat
            dlon = math.cos(math.radians(glat)) * (gal_sphere.lon - glon)
            diff = 3600.0 * math.hypot(dlon, dlat)
            if diff > THRESHOLD_SECONDS:
                print('PY Test_GAL_EQJ_NOVAS({} line {}): EXCESSIVE ERROR = {:0.3f}'.format(filename, lnum, diff))
                sys.exit(1)
            if diff > max_diff:
                max_diff = diff
        Debug('PY Test_GAL_EQJ_NOVAS: PASS. max_diff = {:0.3f} arcseconds.'.format(max_diff))
        return 0

def Test_EQJ_EQD(body):
    # Verify conversion of equatorial J2000 to equatorial of-date, and back.
    # Use established functions to calculate spherical coordinates for the body, in both EQJ and EQD.
    time = astronomy.Time.Make(2019, 12, 8, 20, 50, 0)
    observer = astronomy.Observer(+35, -85, 0)
    eq2000 = astronomy.Equator(body, time, observer, False, True)
    eqdate = astronomy.Equator(body, time, observer, True, True)

    # Convert EQJ spherical coordinates to vector.
    v2000 = eq2000.vec

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
    vec_eqd = eqd.vec

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

    if diff_alt > 4.0e-14 or diff_az > 1.2e-13:
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
    vec_eqj = eqj.vec
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

    # Round trip #7: EQD <==> ECT
    eqd_ect = astronomy.Rotation_EQD_ECT(time)
    ect_eqd = astronomy.Rotation_ECT_EQD(time)
    CheckInverse('eqd_ect', 'ect_eqd', eqd_ect, ect_eqd)

    # Round trip #8: EQJ <==> ECT
    eqj_ect = astronomy.Rotation_EQJ_ECT(time)
    ect_eqj = astronomy.Rotation_ECT_EQJ(time)
    CheckInverse('eqj_ect', 'ect_eqj', eqj_ect, ect_eqj)

    # Verify that combining different sequences of rotations result
    # in the expected combination.
    # For example, (EQJ ==> HOR ==> ECL) must be the same matrix as (EQJ ==> ECL).
    CheckCycle('eqj_ecl, ecl_eqd, eqd_eqj', eqj_ecl, ecl_eqd, eqd_eqj)
    CheckCycle('eqj_hor, hor_ecl, ecl_eqj', eqj_hor, hor_ecl, ecl_eqj)
    CheckCycle('eqj_hor, hor_eqd, eqd_eqj', eqj_hor, hor_eqd, eqd_eqj)
    CheckCycle('ecl_eqd, eqd_hor, hor_ecl', ecl_eqd, eqd_hor, hor_ecl)
    CheckCycle('eqj_eqd, eqd_ect, ect_eqj', eqj_eqd, eqd_ect, ect_eqj)

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

def Test_EQD_ECT():
    time = astronomy.Time.Make(1900, 1, 1, 0, 0, 0.0)
    stopTime = astronomy.Time.Make(2100, 1, 1, 0, 0, 0.0)
    count = 0
    max_diff = 0.0
    while time.ut <= stopTime.ut:
        # Get Moon's geocentric position in EQJ.
        eqj = astronomy.GeoMoon(time)
        # Convert EQJ to EQD.
        eqj_eqd = astronomy.Rotation_EQJ_EQD(time)
        eqd = astronomy.RotateVector(eqj_eqd, eqj)
        # Convert EQD to ECT.
        eqd_ect = astronomy.Rotation_EQD_ECT(time)
        ect = astronomy.RotateVector(eqd_ect, eqd)
        # Independently get the Moon's spherical coordinates in ECT.
        sphere = astronomy.EclipticGeoMoon(time)
        # Convert spherical coordinates to ECT vector.
        check_ect = astronomy.VectorFromSphere(sphere, time)
        # Verify the two ECT vectors are identical, within tolerance.
        max_diff = max(max_diff, VectorDiff(ect, check_ect))
        time = time.AddDays(10.0)
        count += 1
    if max_diff > 3.743e-18:
        print('PY Test_EQD_ECT: excessive vector diff = {:0.6e} au.'.format(max_diff))
        sys.exit(1)
    Debug('PY Test_EQD_ECT: PASS: count = {}, max_diff = {:0.6e} au.'.format(count, max_diff))


def Ecliptic():
    time = astronomy.Time.Make(1900, 1, 1, 0, 0, 0.0)
    stopTime = astronomy.Time.Make(2100, 1, 1, 0, 0, 0.0)
    count = 0
    max_vec_diff = 0
    max_angle_diff = 0.0
    while time.ut <= stopTime.ut:
        # Get Moon's geocentric position in EQJ.
        eqj = astronomy.GeoMoon(time)

        # Convert EQJ to ECT.
        eclip = astronomy.Ecliptic(eqj)

        # Confirm that the ecliptic angles and ecliptic vector are consistent.
        check_sphere = astronomy.Spherical(eclip.elat, eclip.elon, eclip.vec.Length())
        check_vec = astronomy.VectorFromSphere(check_sphere, time)
        max_angle_diff = max(max_angle_diff, VectorDiff(eclip.vec, check_vec))

        # Independently get the Moon's spherical coordinates in ECT.
        sphere = astronomy.EclipticGeoMoon(time)

        # Convert spherical coordinates to ECT vector.
        check_ect = astronomy.VectorFromSphere(sphere, time)

        # Verify the two ECT vectors are identical, within tolerance.
        max_vec_diff = max(max_vec_diff, VectorDiff(eclip.vec, check_ect))

        time = time.AddDays(10.0)
        count += 1

    if max_vec_diff > 3.388e-18:
        return Fail('Ecliptic', 'EXCESSIVE VECTOR DIFF = {:0.6e} au.'.format(max_vec_diff))

    if max_angle_diff > 3.007e-18:
        return Fail('Ecliptic', 'EXCESSIVE ANGLE DIFF = {:0.6e} au.'.format(max_angle_diff))

    print('PY Ecliptic: PASS: count = {:d}, max_vec_diff = {:0.6e} au, max_angle_diff = {:0.6e} au.'.format(count, max_vec_diff, max_angle_diff))
    return 0


def Rotation():
    Rotation_MatrixInverse()
    Rotation_MatrixMultiply()
    Rotation_Pivot()
    Test_GAL_EQJ_NOVAS('temp/galeqj.txt')
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
    Test_EQD_ECT()
    Test_RotRoundTrip()
    print('PY Rotation: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def PrecisionProblem():
    # [2024-05-28] Suddenly Python calculations diverged from the other languages.
    # This is to help me narrow in on where the numeric precisions has gone astray.
    # First  file: temp\c_check.txt
    # Second file: temp\py_check.txt
    # Tolerance = 7.110e-16
    #
    #             lnum                 a_value                 b_value     factor       diff  name
    # FAIL       22073  2.4668955106186205e-01  2.4668955106091364e-01    3.25733  3.089e-12  helio_x
    #
    # v Mercury -1.019385923567141144e+05 2.466895510618620502e-01 -2.985844027934208555e-01 -1.851362492545978455e-01
    tt = -101938.59235671411
    time = astronomy.Time.FromTerrestrialTime(tt)
    #print(time)
    vec = astronomy.HelioVector(astronomy.Body.Mercury, time)
    print(vec)
    # Known correct values from C code...
    cx = +2.466895510618620502e-01
    cy = -2.985844027934208555e-01
    cz = -1.851362492545978455e-01
    factor = 3.25733
    dx = (cx - vec.x) * factor
    dy = (cy - vec.y) * factor
    dz = (cz - vec.z) * factor
    tolerance = 7.110e-16
    diff = math.sqrt(dx*dx + dy*dy + dz*dz)
    print('PY PrecisionProblem: dx={:g}, dy={:g}, dz={:g}, diff={:g}'.format(dx, dy, dz, diff))
    if diff > tolerance:
        print('PY PrecisionProblem: EXCESSIVE ERROR')
        return 1
    print('PY PrecisionProblem: PASS')
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
        period = astronomy.PlanetOrbitalPeriod(body)
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
    dt = (expected_peak.tt - eclipse.peak.tt) * SECONDS_PER_DAY
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

            # Make sure numeric data are finite numbers.
            v(eclipse.obscuration)
            v(eclipse.sd_partial)
            v(eclipse.sd_penum)
            v(eclipse.sd_total)

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
            sd_valid = False
            frac_valid = False
            # Verify that the calculated eclipse semi-durations are consistent with the kind.
            # Verify that obscurations also make sense for the kind.
            if eclipse.kind == astronomy.EclipseKind.Penumbral:
                sd_valid = (eclipse.sd_penum > 0.0) and (eclipse.sd_partial == 0.0) and (eclipse.sd_total == 0.0)
                frac_valid = (eclipse.obscuration == 0.0)
            elif eclipse.kind == astronomy.EclipseKind.Partial:
                sd_valid = (eclipse.sd_penum > 0.0) and (eclipse.sd_partial > 0.0) and (eclipse.sd_total == 0.0)
                frac_valid = (0.0 < eclipse.obscuration < 1.0)
            elif eclipse.kind == astronomy.EclipseKind.Total:
                sd_valid = (eclipse.sd_penum > 0.0) and (eclipse.sd_partial > 0.0) and (eclipse.sd_total > 0.0)
                frac_valid = (eclipse.obscuration == 1.0)
            else:
                print('PY LunarEclipse({} line {}): invalid eclipse kind {}.'.format(filename, lnum, eclipse.kind))
                return 1

            if not sd_valid:
                print('PY LunarEclipse({} line {}): invalid semidurations.'.format(filename, lnum))
                return 1

            if not frac_valid:
                print('PY LunarEclipse({} line {}): invalid obscuration {:0.8f} for eclipsekind {}.'.format(filename, lnum, eclipse.obscuration, eclipse.kind))

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
    print("PY LunarEclipse: PASS (verified {}, skipped {}, max_diff_minutes = {}, avg_diff_minutes = {})".format(lnum, skip_count, max_diff_minutes, (sum_diff_minutes / diff_count)))
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
            if diff_minutes > 7.56:
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

            # Verify the obscuration value is consistent with the eclipse kind.
            if eclipse.kind == astronomy.EclipseKind.Partial:
                if eclipse.obscuration is not None:
                    print('PY GlobalSolarEclipse({} line {}): Expected obscuration = None for partial eclipse, but found {}'.format(filename, lnum, eclipse.obscuration))
                    return 1
            elif eclipse.kind == astronomy.EclipseKind.Annular:
                if not (0.8 < v(eclipse.obscuration) < 1.0):
                    print('PY GlobalSolarEclipse({} line {}): Invalid obscuration = {:0.8f} for annular eclipse.'.format(filename, lnum, eclipse.obscuration))
                    return 1
            elif eclipse.kind == astronomy.EclipseKind.Total:
                if v(eclipse.obscuration) != 1.0:
                    print('PY GlobalSolarEclipse({} line {}): Invalid obscuration = {:0.8f} for total eclipse.'.format(filename, lnum, eclipse.obscuration))
                    return 1
            else:
                print('PY GlobalSolarEclipse({} line {}): Unhandled eclipse kind {}'.format(filename, lnum, eclipse.kind))
                return 1

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
            funcname = 'LocalSolarEclipse({} line {})'.format(filename, lnum)
            # 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
            token = line.split()
            if len(token) != 5:
                return Fail(funcname, 'invalid token count = {}'.format(len(token)))
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
            if diff_minutes > 7.737:
                return Fail(funcname, 'EXCESSIVE TIME ERROR = {} minutes'.format(diff_minutes))

            if diff_minutes > max_minutes:
                max_minutes = diff_minutes

            # Verify obscuration makes sense for this kind of eclipse.
            v(eclipse.obscuration)
            if eclipse.kind in [astronomy.EclipseKind.Annular, astronomy.EclipseKind.Partial]:
                frac_valid = (0.0 < eclipse.obscuration < 1.0)
            elif eclipse.kind == astronomy.EclipseKind.Total:
                frac_valid = (eclipse.obscuration == 1.0)
            else:
                return Fail(funcname, 'Invalid eclipse kind {}'.format(eclipse.kind))
            if not frac_valid:
                return Fail(funcname, 'Invalid eclipse obscuration {:0.8f} for {} eclipse.'.format(eclipse.obscuration, eclipse.kind))

    funcname = 'LocalSolarEclipse1({})'.format(filename)

    if lnum != expected_count:
        return Fail(funcname, 'WRONG LINE COUNT = {}, expected {}'.format(lnum, expected_count))

    if skip_count > 6:
        return Fail(funcname, 'EXCESSIVE SKIP COUNT = {}'.format(skip_count))

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
        # Ignore discrepancies for negative altitudes, because of quirky and irrelevant differences in refraction models.
        if expect.altitude >= 0.0:
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
    return (
        LocalSolarEclipse1() or
        LocalSolarEclipse2()
    )

#-----------------------------------------------------------------------------------------------------------

def GlobalAnnularCase(year, month, day, obscuration):
    # Search for the first solar eclipse that occurs after the given date.
    time = astronomy.Time.Make(year, month, day, 0, 0, 0.0)
    eclipse = astronomy.SearchGlobalSolarEclipse(time)
    funcname = 'GlobalAnnularCase({:04d}-{:02d}-{:02d})'.format(year, month, day)

    # Verify the eclipse is within 1 day after the search basis time.
    dt = v(eclipse.peak.ut - time.ut)
    if not (0.0 <= dt <= 1.0):
        return Fail(funcname, 'found eclipse {:0.4f} days after search time.'.format(dt))

    # Verify we found an annular solar eclipse.
    if eclipse.kind != astronomy.EclipseKind.Annular:
        return Fail(funcname, 'expected annular eclipse but found {}'.format(eclipse.kind))

    # Check how accurately we calculated obscuration.
    diff = v(eclipse.obscuration - obscuration)
    if abs(diff) > 0.0000904:
        return Fail(funcname, 'excessive obscuration error = {:0.8f}, expected = {:0.8f}, actual = {:0.8f}'.format(diff, obscuration, eclipse.obscuration))

    Debug('{}: obscuration error = {:11.8f}'.format(funcname, diff))
    return 0


def LocalSolarCase(year, month, day, latitude, longitude, kind, obscuration, tolerance):
    funcname = 'LocalSolarCase({:04d}-{:02d}-{:02d})'.format(year, month, day)
    time = astronomy.Time.Make(year, month, day, 0, 0, 0.0)
    observer = astronomy.Observer(latitude, longitude, 0.0)
    eclipse = astronomy.SearchLocalSolarEclipse(time, observer)
    dt = v(eclipse.peak.time.ut - time.ut)
    if not (0.0 <= dt <= 1.0):
        return Fail(funcname, 'eclipse found {:0.4f} days after search date'.format(dt))

    if eclipse.kind != kind:
        return Fail(funcname, 'expected {} eclipse, but found {}.'.format(kind, eclipse.kind))

    diff = v(eclipse.obscuration - obscuration)
    if abs(diff) > tolerance:
        return Fail(funcname, 'obscuration diff = {:0.8f}, expected = {:0.8f}, actual = {:0.8f}'.format(diff, obscuration, eclipse.obscuration))

    Debug('{}: obscuration diff = {:11.8f}'.format(funcname, diff))
    return 0


def SolarFraction():
    return (
        # Verify global solar eclipse obscurations for annular eclipses only.
        # This is because they are the only nontrivial values for global solar eclipses.
        # The trivial values are all validated exactly by GlobalSolarEclipseTest().

        GlobalAnnularCase(2023, 10, 14, 0.90638) or    # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2023Oct14Aprime.html
        GlobalAnnularCase(2024, 10,  2, 0.86975) or    # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2024Oct02Aprime.html
        GlobalAnnularCase(2027,  2,  6, 0.86139) or    # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2027Feb06Aprime.html
        GlobalAnnularCase(2028,  1, 26, 0.84787) or    # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2028Jan26Aprime.html
        GlobalAnnularCase(2030,  6,  1, 0.89163) or    # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2030Jun01Aprime.html

        # Verify obscuration values for specific locations on the Earth.
        # Local solar eclipse calculations include obscuration for all types of eclipse, not just annular and total.
        LocalSolarCase(2023, 10, 14,  11.3683,  -83.1017, astronomy.EclipseKind.Annular, 0.90638, 0.000080) or  # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2023Oct14Aprime.html
        LocalSolarCase(2023, 10, 14,  25.78,    -80.22,   astronomy.EclipseKind.Partial, 0.578,   0.000023) or  # https://aa.usno.navy.mil/calculated/eclipse/solar?eclipse=22023&lat=25.78&lon=-80.22&label=Miami%2C+FL&height=0&submit=Get+Data
        LocalSolarCase(2023, 10, 14,  30.2666,  -97.7000, astronomy.EclipseKind.Partial, 0.8867,  0.001016) or  # http://astro.ukho.gov.uk/eclipse/0332023/Austin_TX_United_States_2023Oct14.png
        LocalSolarCase(2024,  4,  8,  25.2900, -104.1383, astronomy.EclipseKind.Total,   1.0,     0.0     ) or  # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2024Apr08Tprime.html
        LocalSolarCase(2024,  4,  8,  37.76,   -122.44,   astronomy.EclipseKind.Partial, 0.340,   0.000604) or  # https://aa.usno.navy.mil/calculated/eclipse/solar?eclipse=12024&lat=37.76&lon=-122.44&label=San+Francisco%2C+CA&height=0&submit=Get+Data
        LocalSolarCase(2024, 10,  2, -21.9533, -114.5083, astronomy.EclipseKind.Annular, 0.86975, 0.000061) or  # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2024Oct02Aprime.html
        LocalSolarCase(2024, 10,  2, -33.468,   -70.636,  astronomy.EclipseKind.Partial, 0.436,   0.000980) or  # https://aa.usno.navy.mil/calculated/eclipse/solar?eclipse=22024&lat=-33.468&lon=-70.636&label=Santiago%2C+Chile&height=0&submit=Get+Data
        LocalSolarCase(2030,  6,  1,  56.525,    80.0617, astronomy.EclipseKind.Annular, 0.89163, 0.000067) or  # https://www.eclipsewise.com/solar/SEprime/2001-2100/SE2030Jun01Aprime.html
        LocalSolarCase(2030,  6,  1,  40.388,    49.914,  astronomy.EclipseKind.Partial, 0.67240, 0.000599) or  # http://xjubier.free.fr/en/site_pages/SolarEclipseCalc_Diagram.html
        LocalSolarCase(2030,  6,  1,  40.3667,   49.8333, astronomy.EclipseKind.Partial, 0.6736,  0.001464) or  # http://astro.ukho.gov.uk/eclipse/0132030/Baku_Azerbaijan_2030Jun01.png

        Pass('SolarFraction')
    )

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
    if PlutoCheckDate(  +18250.0,  0.089, +37.4377303523676090, -10.2466292454075898, -14.4773101310875809): return 1
    if PlutoCheckDate( -856493.0,  4.067, +23.4292113199166252, +42.1452685817740829,  +6.0580908436642940): return 1
    if PlutoCheckDate( +435633.0,  0.016, -27.3178902095231813, +18.5887022581070305, +14.0493896259306936): return 1
    if PlutoCheckDate(       0.0,   8e-9,  -9.8753673425269000, -27.9789270580402771,  -5.7537127596369588): return 1
    if PlutoCheckDate( +800916.0,  2.286, -29.5266052645301365, +12.0554287322176474, +12.6878484911631091): return 1
    print("PY PlutoCheck: PASS")
    return 0

#-----------------------------------------------------------------------------------------------------------

def GeoidTestCase(time, observer, ofdate):
    topo_moon = astronomy.Equator(astronomy.Body.Moon, time, observer, ofdate, False)
    surface = astronomy.ObserverVector(time, observer, ofdate)
    geo_moon = astronomy.GeoVector(astronomy.Body.Moon, time, False)

    if ofdate:
        # GeoVector() returns J2000 coordinates. Convert to equator-of-date coordinates.
        rot = astronomy.Rotation_EQJ_EQD(time)
        geo_moon = astronomy.RotateVector(rot, geo_moon)

    dx = astronomy.KM_PER_AU * v((geo_moon.x - surface.x) - topo_moon.vec.x)
    dy = astronomy.KM_PER_AU * v((geo_moon.y - surface.y) - topo_moon.vec.y)
    dz = astronomy.KM_PER_AU * v((geo_moon.z - surface.z) - topo_moon.vec.z)
    diff = sqrt(dx*dx + dy*dy + dz*dz)
    Debug('PY GeoidTestCase: ofdate={}, time={}, obs={}, surface=({}, {}, {}), diff = {} km'.format(
        ofdate,
        time,
        observer,
        astronomy.KM_PER_AU * surface.x,
        astronomy.KM_PER_AU * surface.y,
        astronomy.KM_PER_AU * surface.z,
        diff
    ))

    # Require 1 millimeter accuracy! (one millionth of a kilometer).
    if diff > 1.0e-6:
        print('PY GeoidTestCase: EXCESSIVE POSITION ERROR.')
        return 1

    # Verify that we can convert the surface vector back to an observer.
    vobs = astronomy.VectorObserver(surface, ofdate)
    lat_diff = vabs(vobs.latitude - observer.latitude)

    # Longitude is meaningless at the poles, so don't bother checking it there.
    if -89.99 <= observer.latitude <= +89.99:
        lon_diff = vabs(vobs.longitude - observer.longitude)
        if lon_diff > 180.0:
            lon_diff = 360.0 - lon_diff
        lon_diff = vabs(lon_diff * math.cos(math.degrees(observer.latitude)))
        if lon_diff > 1.0e-6:
            print('PY GeoidTestCase: EXCESSIVE longitude check error = {}'.format(lon_diff))
            return 1
    else:
        lon_diff = 0.0

    h_diff = vabs(vobs.height - observer.height)
    Debug('PY GeoidTestCase: vobs={}, lat_diff={}, lon_diff={}, h_diff={}'.format(vobs, lat_diff, lon_diff, h_diff))
    if lat_diff > 1.0e-6:
        print('PY GeoidTestCase: EXCESSIVE latitude check error = {}'.format(lat_diff))
        return 1
    if h_diff > 0.001:
        print('PY GeoidTestCase: EXCESSIVE height check error = {}'.format(h_diff))
        return 1
    return 0


def Geoid():
    time_list = [
        astronomy.Time.Parse('1066-09-27T18:00:00Z'),
        astronomy.Time.Parse('1970-12-13T15:42:00Z'),
        astronomy.Time.Parse('1970-12-13T15:43:00Z'),
        astronomy.Time.Parse('2015-03-05T02:15:45Z')
    ]

    observer_list = [
        astronomy.Observer(  0.0,    0.0,    0.0),
        astronomy.Observer( +1.5,   +2.7,    7.4),
        astronomy.Observer( -1.5,   -2.7,    7.4),
        astronomy.Observer(-53.7, +141.7,  100.0),
        astronomy.Observer(+30.0,  -85.2,  -50.0),
        astronomy.Observer(+90.0,  +45.0,  -50.0),
        astronomy.Observer(-90.0, -180.0,    0.0),
        astronomy.Observer(-89.0,  -81.0, 1234.0),
        astronomy.Observer(+89.0, -103.4,  279.8),
        astronomy.Observer(+48.2,   24.5, 2019.0),
        astronomy.Observer(+28.5,  -82.3,   -3.4)
    ]

    # Test hand-crafted locations.

    for observer in observer_list:
        for time in time_list:
            if GeoidTestCase(time, observer, False):
                return 1
            if GeoidTestCase(time, observer, True):
                return 1

    # More exhaustive tests for a single time value across many different geographic coordinates.
    # Solving for latitude is the most complicated part of VectorObserver, so
    # I test for every 1-degree increment of latitude, but with 5-degree increments for longitude.
    time = astronomy.Time.Parse('2021-06-20T15:08:00Z')
    lat = -90
    while lat <= +90:
        lon = -175
        while lon <= +180:
            observer = astronomy.Observer(lat, lon, 0.0)
            if GeoidTestCase(time, observer, True):
                return 1
            lon += 5
        lat += 1

    print('PY GeoidTest: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def JupiterMoons_CheckJpl(mindex, tt, pos, vel):
    pos_tolerance = 9.0e-4
    vel_tolerance = 9.0e-4
    time = astronomy.Time.FromTerrestrialTime(tt)
    jm = astronomy.JupiterMoons(time)
    moon = SelectJupiterMoon(jm, mindex)

    dx = v(pos[0] - moon.x)
    dy = v(pos[1] - moon.y)
    dz = v(pos[2] - moon.z)
    mag = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2])
    pos_diff = sqrt(dx*dx + dy*dy + dz*dz) / mag
    if pos_diff > pos_tolerance:
        print('PY JupiterMoons_CheckJpl(mindex={}, tt={}): excessive position error {}'.format(mindex, tt, pos_diff))
        return 1

    dx = v(vel[0] - moon.vx)
    dy = v(vel[1] - moon.vy)
    dz = v(vel[2] - moon.vz)
    mag = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2])
    vel_diff = sqrt(dx*dx + dy*dy + dz*dz) / mag
    if vel_diff > vel_tolerance:
        print('PY JupiterMoons_CheckJpl(mindex={}, tt={}): excessive velocity error {}'.format(mindex, tt, vel_diff))
        return 1

    Debug('PY JupiterMoons_CheckJpl: mindex={}, tt={}, pos_diff={}, vel_diff={}'.format(mindex, tt, pos_diff, vel_diff))
    return 0


def JupiterMoons():
    for mindex in range(4):
        filename = 'jupiter_moons/horizons/jm{}.txt'.format(mindex)
        with open(filename, 'rt') as infile:
            lnum = 0
            found = False
            part = -1
            expected_count = 5001
            count = 0
            for line in infile:
                line = line.rstrip()
                lnum += 1
                if not found:
                    if line == '$$SOE':
                        found = True
                        part = 0
                    elif line.startswith('Revised:'):
                        check_mindex = int(line[76:]) - 501
                        if mindex != check_mindex:
                            print('PY JupiterMoons({} line {}): moon index does not match: check={}, mindex={}'.format(filename, lnum, check_mindex, mindex))
                            return 1
                elif line == '$$EOE':
                    break
                else:
                    if part == 0:
                        # 2446545.000000000 = A.D. 1986-Apr-24 12:00:00.0000 TDB
                        tt = float(line.split()[0]) - 2451545.0    # convert JD to J2000 TT
                    elif part == 1:
                        # X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                        match = re.match(r'\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)', line)
                        if not match:
                            print('PY JupiterMoons({} line {}): cannot parse position vector.'.format(filename, lnum))
                            return 1
                        pos = [ float(match.group(1)), float(match.group(2)), float(match.group(3)) ]
                    else:   # part == 2
                        # VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                        match = re.match(r'\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)', line)
                        if not match:
                            print('PY JupiterMoons({} line {}): cannot parse velocity vector.'.format(filename, lnum))
                            return 1
                        vel = [ float(match.group(1)), float(match.group(2)), float(match.group(3)) ]
                        if JupiterMoons_CheckJpl(mindex, tt, pos, vel):
                            print('PY JupiterMoons({} line {}): FAILED VERIFICATION.'.format(filename, lnum))
                            return 1
                        count += 1
                    part = (part + 1) % 3
            if count != expected_count:
                print('PY JupiterMoons: expected {} test cases, but found {}'.format(expected_count, count))
                return 1

    print('PY JupiterMoons: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def Issue103():
    # https://github.com/cosinekitty/astronomy/issues/103
    observer = astronomy.Observer(29, -81, 10)
    ut = -8.817548982869034808e+04
    time = astronomy.Time(ut)
    body = astronomy.Body.Venus
    ofdate = astronomy.Equator(body, time, observer, True, True)
    hor = astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, astronomy.Refraction.Airless)
    print('tt  = {:23.16f}'.format(time.tt))
    print('az  = {:23.16f}'.format(hor.azimuth))
    print('alt = {:23.16f}'.format(hor.altitude))
    return 0

#-----------------------------------------------------------------------------------------------------------

class _bary_stats_t:
    def __init__(self):
        self.max_rdiff = 0.0
        self.max_vdiff = 0.0

def StateVectorDiff(relative, vec, x, y, z):
    dx = v(vec[0] - x)
    dy = v(vec[1] - y)
    dz = v(vec[2] - z)
    diff_squared = dx*dx + dy*dy + dz*dz
    if relative:
        diff_squared /= (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])
    return sqrt(diff_squared)

#-----------------------------------------------------------------------------------------------------------

def VerifyState(func, stats, filename, lnum, time, pos, vel, r_thresh, v_thresh):
    state = func.Eval(time)

    rdiff = StateVectorDiff((r_thresh > 0.0), pos, state.x, state.y, state.z)
    if rdiff > stats.max_rdiff:
        stats.max_rdiff = rdiff

    vdiff = StateVectorDiff((v_thresh > 0.0), vel, state.vx, state.vy, state.vz)
    if vdiff > stats.max_vdiff:
        stats.max_vdiff = vdiff

    if rdiff > abs(r_thresh):
        print('PY VerifyState({} line {}): EXCESSIVE position error = {:0.4e}'.format(filename, lnum, rdiff))
        return 1

    if vdiff > abs(v_thresh):
        print('PY VerifyState({} line {}): EXCESSIVE velocity error = {:0.4e}'.format(filename, lnum, vdiff))
        return 1

    return 0


class JplStateRecord:
    def __init__(self, lnum, state):
        self.lnum = lnum
        self.state = state


def JplHorizonsStateVectors(filename):
    with open(filename, 'rt') as infile:
        lnum = 0
        part = 0
        found_begin = False
        for line in infile:
            line = line.rstrip()
            lnum += 1
            if not found_begin:
                if line == '$$SOE':
                    found_begin = True
            elif line == '$$EOE':
                break
            else:
                if part == 0:
                    # 2446545.000000000 = A.D. 1986-Apr-24 12:00:00.0000 TDB
                    tt = float(line.split()[0]) - 2451545.0    # convert JD to J2000 TT
                    time = astronomy.Time.FromTerrestrialTime(tt)
                elif part == 1:
                    # X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                    match = re.match(r'\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)', line)
                    if not match:
                        print('PY JplHorizonsStateVectors({} line {}): cannot parse position vector.'.format(filename, lnum))
                        return 1
                    rx, ry, rz = float(match.group(1)), float(match.group(2)), float(match.group(3))
                else:   # part == 2
                    # VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                    match = re.match(r'\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)', line)
                    if not match:
                        print('PY JplHorizonsStateVectors({} line {}): cannot parse velocity vector.'.format(filename, lnum))
                        return 1
                    vx, vy, vz = float(match.group(1)), float(match.group(2)), float(match.group(3))
                    yield JplStateRecord(lnum, astronomy.StateVector(rx, ry, rz, vx, vy, vz, time))
                part = (part + 1) % 3
    return 0


def VerifyStateBody(func, filename, r_thresh, v_thresh):
    stats = _bary_stats_t()
    count = 0
    for rec in JplHorizonsStateVectors(filename):
        time = rec.state.t
        pos = [rec.state.x, rec.state.y, rec.state.z]
        vel = [rec.state.vx, rec.state.vy, rec.state.vz]
        if VerifyState(func, stats, filename, rec.lnum, time, pos, vel, r_thresh, v_thresh):
            print('PY VerifyStateBody({} line {}): FAILED VERIFICATION.'.format(filename, rec.lnum))
            return 1
        count += 1
    Debug('PY VerifyStateBody({}): PASS - Tested {} cases. max rdiff={:0.3e}, vdiff={:0.3e}'.format(filename, count, stats.max_rdiff, stats.max_vdiff))
    return 0

#-----------------------------------------------------------------------------------------------------------

# Constants for use inside unit tests only; they doesn't make sense for public consumption.
_Body_GeoMoon = -100
_Body_Geo_EMB = -101

class BaryStateFunc:
    def __init__(self, body):
        self.body = body

    def Eval(self, time):
        if self.body == _Body_GeoMoon:
            return astronomy.GeoMoonState(time)
        if self.body == _Body_Geo_EMB:
            return astronomy.GeoEmbState(time)
        return astronomy.BaryState(self.body, time)

def BaryState():
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Sun),     'barystate/Sun.txt',     -1.224e-05, -1.134e-07):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Mercury), 'barystate/Mercury.txt',  1.672e-04,  2.698e-04):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Venus),   'barystate/Venus.txt',    4.123e-05,  4.308e-05):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Earth),   'barystate/Earth.txt',    2.296e-05,  6.359e-05):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Mars),    'barystate/Mars.txt',     3.107e-05,  5.550e-05):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Jupiter), 'barystate/Jupiter.txt',  7.389e-05,  2.471e-04):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Saturn),  'barystate/Saturn.txt',   1.067e-04,  3.220e-04):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Uranus),  'barystate/Uranus.txt',   9.035e-05,  2.519e-04):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Neptune), 'barystate/Neptune.txt',  9.838e-05,  4.446e-04):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Pluto),   'barystate/Pluto.txt',    4.259e-05,  7.827e-05):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.Moon),    "barystate/Moon.txt",     2.354e-05,  6.604e-05):  return 1
    if VerifyStateBody(BaryStateFunc(astronomy.Body.EMB),     "barystate/EMB.txt",      2.353e-05,  6.511e-05):  return 1
    if VerifyStateBody(BaryStateFunc(_Body_GeoMoon),          "barystate/GeoMoon.txt",  4.086e-05,  5.347e-05):  return 1
    if VerifyStateBody(BaryStateFunc(_Body_Geo_EMB),          "barystate/GeoEMB.txt",   4.076e-05,  5.335e-05):  return 1
    print('PY BaryState: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

class HelioStateFunc:
    def __init__(self, body):
        self.body = body

    def Eval(self, time):
        return astronomy.HelioState(self.body, time)

def HelioState():
    if VerifyStateBody(HelioStateFunc(astronomy.Body.SSB),     'heliostate/SSB.txt',     -1.209e-05, -1.125e-07): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Mercury), 'heliostate/Mercury.txt',  1.481e-04,  2.756e-04): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Venus),   'heliostate/Venus.txt',    3.528e-05,  4.485e-05): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Earth),   'heliostate/Earth.txt',    1.476e-05,  6.105e-05): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Mars),    'heliostate/Mars.txt',     3.154e-05,  5.603e-05): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Jupiter), 'heliostate/Jupiter.txt',  7.455e-05,  2.562e-04): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Saturn),  'heliostate/Saturn.txt',   1.066e-04,  3.150e-04): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Uranus),  'heliostate/Uranus.txt',   9.034e-05,  2.712e-04): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Neptune), 'heliostate/Neptune.txt',  9.834e-05,  4.534e-04): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Pluto),   'heliostate/Pluto.txt',    4.271e-05,  1.198e-04): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.Moon),    'heliostate/Moon.txt',     1.477e-05,  6.195e-05): return 1
    if VerifyStateBody(HelioStateFunc(astronomy.Body.EMB),     'heliostate/EMB.txt',      1.476e-05,  6.106e-05): return 1
    print('PY HelioState: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

class TopoStateFunc:
    def __init__(self, body):
        self.body = body

    def Eval(self, time):
        observer = astronomy.Observer(30.0, -80.0, 1000.0)
        observer_state = astronomy.ObserverState(time, observer, False)
        if self.body == _Body_Geo_EMB:
            state = astronomy.GeoEmbState(time)
        elif self.body == astronomy.Body.Earth:
            state = astronomy.StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time)
        else:
            raise Exception('PY TopoStateFunction: unsupported body ' + self.body)
        state.x  -= observer_state.x
        state.y  -= observer_state.y
        state.z  -= observer_state.z
        state.vx -= observer_state.vx
        state.vy -= observer_state.vy
        state.vz -= observer_state.vz
        return state

def TopoState():
    if VerifyStateBody(TopoStateFunc(astronomy.Body.Earth), 'topostate/Earth_N30_W80_1000m.txt', 2.108e-04, 2.430e-04): return 1
    if VerifyStateBody(TopoStateFunc(_Body_Geo_EMB),        'topostate/EMB_N30_W80_1000m.txt',   7.197e-04, 2.497e-04): return 1
    print('PY TopoState: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def Aberration():
    THRESHOLD_SECONDS = 0.453
    filename = 'equatorial/Mars_j2000_ofdate_aberration.txt'
    count = 0
    with open(filename, 'rt') as infile:
        lnum = 0
        found_begin = False
        max_diff_seconds = 0.0
        for line in infile:
            lnum += 1
            line = line.rstrip()
            if not found_begin:
                if line == '$$SOE':
                    found_begin = True
            elif line == '$$EOE':
                break
            else:
                # 2459371.500000000 *   118.566080210  22.210647456 118.874086738  22.155784122
                token = line.split()
                if len(token) < 5:
                    print('PY Aberration({} line {}): not enough tokens'.format(filename, lnum))
                    return 1

                jd = float(token[0])
                jra = float(token[-4])
                jdec = float(token[-3])
                dra = float(token[-2])
                ddec = float(token[-1])

                # Convert julian day value to AstroTime.
                time = astronomy.Time(jd - 2451545.0)

                # Convert EQJ angular coordinates (jra, jdec) to an EQJ vector.
                # Make the maginitude of the vector the speed of light,
                # to prepare for aberration correction.
                eqj_sphere = astronomy.Spherical(jdec, jra, astronomy.C_AUDAY)
                eqj_vec = astronomy.VectorFromSphere(eqj_sphere, time)

                # Aberration correction: calculate the Earth's barycentric
                # velocity vector in EQJ coordinates.
                eqj_earth = astronomy.BaryState(astronomy.Body.Earth, time)

                # Use non-relativistic approximation: add light vector to Earth velocity vector.
                # This gives aberration-corrected apparent position of the start in EQJ.
                eqj_vec.x += eqj_earth.vx
                eqj_vec.y += eqj_earth.vy
                eqj_vec.z += eqj_earth.vz

                # Calculate the rotation matrix that converts J2000 coordinates to of-date coordinates.
                rot = astronomy.Rotation_EQJ_EQD(time)

                # Use the rotation matrix to re-orient the EQJ vector to an EQD vector.
                eqd_vec = astronomy.RotateVector(rot, eqj_vec)

                # Convert the EQD vector back to spherical angular coordinates.
                eqd_sphere = astronomy.SphereFromVector(eqd_vec)

                # Calculate the differences in RA and DEC between expected and calculated values.
                factor = math.cos(math.radians(v(eqd_sphere.lat)))    # RA errors are less important toward the poles.
                xra = factor * vabs(eqd_sphere.lon - dra)
                xdec = vabs(eqd_sphere.lat - ddec)
                diff_seconds = 3600.0 * sqrt(xra*xra + xdec*xdec)
                Debug('PY Aberration({} line {}): xra={:0.6f} deg, xdec={:0.6f} deg, diff_seconds={:0.3f}.'.format(filename, lnum, xra, xdec, diff_seconds))
                if diff_seconds > THRESHOLD_SECONDS:
                    print('PY Aberration({} line {}): EXCESSIVE ANGULAR ERROR = {:0.3f} seconds.'.format(filename, lnum, diff_seconds));
                    return 1

                if diff_seconds > max_diff_seconds:
                    max_diff_seconds = diff_seconds

                # We have completed one more test case.
                count += 1

    print('PY AberrationTest({}): PASS - Tested {} cases. max_diff_seconds = {:0.3f}'.format(filename, count, max_diff_seconds))
    return 0

#-----------------------------------------------------------------------------------------------------------

def Twilight():
    tolerance_seconds = 60.0
    max_diff = 0.0
    filename = 'riseset/twilight.txt'
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            tokens = line.split()
            if len(tokens) != 9:
                print('PY Twilight({} line {}): incorrect number of tokens = {}'.format(filename, lnum, len(tokens)))
                return 1
            observer = astronomy.Observer(float(tokens[0]), float(tokens[1]), 0.0)
            searchDate = astronomy.Time.Parse(tokens[2])
            correctTimes = [astronomy.Time.Parse(t) for t in tokens[3:]]
            calcTimes = [
                astronomy.SearchAltitude(astronomy.Body.Sun, observer, astronomy.Direction.Rise, searchDate, 1.0, -18.0),   # astronomical dawn
                astronomy.SearchAltitude(astronomy.Body.Sun, observer, astronomy.Direction.Rise, searchDate, 1.0, -12.0),   # nautical dawn
                astronomy.SearchAltitude(astronomy.Body.Sun, observer, astronomy.Direction.Rise, searchDate, 1.0,  -6.0),   # civil dawn
                astronomy.SearchAltitude(astronomy.Body.Sun, observer, astronomy.Direction.Set,  searchDate, 1.0,  -6.0),   # civil dusk
                astronomy.SearchAltitude(astronomy.Body.Sun, observer, astronomy.Direction.Set,  searchDate, 1.0, -12.0),   # nautical dusk
                astronomy.SearchAltitude(astronomy.Body.Sun, observer, astronomy.Direction.Set,  searchDate, 1.0, -18.0)    # astronomical dusk
            ]
            for i in range(6):
                correct = correctTimes[i]
                calc = calcTimes[i]
                diff = SECONDS_PER_DAY * vabs(calc.ut - correct.ut)
                if diff > tolerance_seconds:
                    print('PY Twilight({} line {}): EXCESSIVE ERROR = {} seconds for case {}'.format(filename, lnum, diff, i))
                    return 1
                if diff > max_diff:
                    max_diff = diff

    print('PY Twilight: PASS ({} test cases, max error = {} seconds)'.format(lnum, max_diff))
    return 0

#-----------------------------------------------------------------------------------------------------------

def LibrationFile(filename):
    max_diff_elon = 0.0
    max_diff_elat = 0.0
    max_diff_distance = 0.0
    max_diff_diam = 0.0
    max_eclip_lon = -900.0
    count = 0
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            token = line.split()
            if lnum == 1:
                if line != "   Date       Time    Phase    Age    Diam    Dist     RA        Dec      Slon      Slat     Elon     Elat   AxisA\n":
                    print('PY LibrationFile({} line {}): unexpected header line'.format(filename, lnum))
                    return 1
            else:
                if len(token) != 16:
                    print('PY LibrationFile({} line {}): expected 16 tokens, found {}'.format(filename, lnum, len(token)))
                    return 1
                day = int(token[0])
                month = MonthNumber(token[1])
                year = int(token[2])
                hmtoken = token[3].split(':')
                if len(hmtoken) != 2:
                    print('PY LibrationFile({} line {}): expected hh:mm but found "{}"'.format(filename, lnum, hmtoken))
                    return 1
                hour = int(hmtoken[0])
                minute = int(hmtoken[1])
                time = astronomy.Time.Make(year, month, day, hour, minute, 0.0)
                diam = float(token[7]) / 3600.0
                dist = float(token[8])
                elon = float(token[13])
                elat = float(token[14])
                lib = astronomy.Libration(time)

                diff_elon = 60.0 * vabs(lib.elon - elon)
                if diff_elon > max_diff_elon:
                    max_diff_elon = diff_elon

                diff_elat = 60.0 * vabs(lib.elat - elat)
                if diff_elat > max_diff_elat:
                    max_diff_elat = diff_elat

                diff_distance = vabs(lib.dist_km - dist)
                if diff_distance > max_diff_distance:
                    max_diff_distance = diff_distance

                diff_diam = vabs(lib.diam_deg - diam)
                if diff_diam > max_diff_diam:
                    max_diff_diam = diff_diam

                if lib.mlon > max_eclip_lon:
                    max_eclip_lon = lib.mlon

                if diff_elon > 0.1304:
                    print('PY LibrationFile({} line {}): EXCESSIVE diff_elon = {}'.format(filename, lnum, diff_elon))
                    return 1

                if diff_elat > 1.6476:
                    print('PY LibrationFile({} line {}): EXCESSIVE diff_elat = {}'.format(filename, lnum, diff_elat))
                    return 1

                if diff_distance > 54.377:
                    print('PY LibrationFile({} line {}): EXCESSIVE diff_distance = {}'.format(filename, lnum, diff_distance))
                    return 1

                if diff_diam > 0.00009:
                    print('PY LibrationFile({} line {}): EXCESSIVE diff_diam = {}'.format(filename, lnum, diff_diam))
                    return 1

                count += 1

    if not (359.0 < max_eclip_lon < 360.0):
        print('PY LibrationFile({}): INVALID max ecliptic longitude = {:0.3f}'.format(filename, max_eclip_lon))
        return 1

    print('PY LibrationFile({}): PASS ({} test cases, max_diff_elon = {} arcmin, max_diff_elat = {} arcmin, max_diff_distance = {} km, max_diff_diam = {} deg)'.format(
        filename, count, max_diff_elon, max_diff_elat, max_diff_distance, max_diff_diam
    ))
    return 0

def Libration():
    return (
        LibrationFile('libration/mooninfo_2020.txt') or
        LibrationFile('libration/mooninfo_2021.txt') or
        LibrationFile('libration/mooninfo_2022.txt')
    )

#-----------------------------------------------------------------------------------------------------------

def Axis():
    if AxisTestBody(astronomy.Body.Sun,      'axis/Sun.txt',       0.0)      :  return 1
    if AxisTestBody(astronomy.Body.Mercury,  'axis/Mercury.txt',   0.074340) :  return 1
    if AxisTestBody(astronomy.Body.Venus,    'axis/Venus.txt',     0.0)      :  return 1
    if AxisTestBody(astronomy.Body.Earth,    'axis/Earth.txt',     0.002032) :  return 1
    if AxisTestBody(astronomy.Body.Moon,     'axis/Moon.txt',      0.264845) :  return 1
    if AxisTestBody(astronomy.Body.Mars,     'axis/Mars.txt',      0.075323) :  return 1
    if AxisTestBody(astronomy.Body.Jupiter,  'axis/Jupiter.txt',   0.000324) :  return 1
    if AxisTestBody(astronomy.Body.Saturn,   'axis/Saturn.txt',    0.000304) :  return 1
    if AxisTestBody(astronomy.Body.Uranus,   'axis/Uranus.txt',    0.0)      :  return 1
    if AxisTestBody(astronomy.Body.Neptune,  'axis/Neptune.txt',   0.000464) :  return 1
    if AxisTestBody(astronomy.Body.Pluto,    'axis/Pluto.txt',     0.0)      :  return 1
    print('PY AxisTest: PASS')
    return 0

def AxisTestBody(body, filename, arcmin_tolerance):
    max_arcmin = 0
    lnum = 0
    count = 0
    found_data = False
    with open(filename, 'rt') as infile:
        for line in infile:
            line = line.strip()
            lnum += 1
            if not found_data:
                if line == '$$SOE':
                    found_data = True
            else:
                if line == '$$EOE':
                    break
                token = line.split()
                # [ '1970-Jan-01', '00:00', '2440587.500000000', '281.01954', '61.41577' ]
                jd = float(token[2])
                ra = float(token[3])
                dec = float(token[4])
                time = astronomy.Time(jd - 2451545.0)
                axis = astronomy.RotationAxis(body, time)
                sphere = astronomy.Spherical(dec, ra, 1.0)
                north = astronomy.VectorFromSphere(sphere, time)
                arcmin = 60.0 * astronomy.AngleBetween(north, axis.north)
                if arcmin > max_arcmin:
                    max_arcmin = arcmin
                count += 1
    Debug('PY AxisTestBody({}): {} test cases, max arcmin error = {}.'.format(body, count, max_arcmin))
    if max_arcmin > arcmin_tolerance:
        print('PY AxisTestBody({}): EXCESSIVE ERROR = {}'.format(body, max_arcmin))
        return 1
    return 0

#-----------------------------------------------------------------------------------------------------------

def MoonNodes():
    filename = 'moon_nodes/moon_nodes.txt'
    with open(filename, 'rt') as infile:
        max_angle = 0.0
        max_minutes = 0.0
        prev_kind = '?'
        lnum = 0
        for line in infile:
            line = line.strip()
            lnum += 1
            token = line.split()
            if len(token) != 4:
                print('PY MoonNodes({} line {}): syntax error'.format(filename, lnum))
                return 1
            kind = token[0]
            if kind not in 'AD':
                print('PY MoonNodes({} line {}): invalid node kind'.format(filename, lnum))
                return 1
            if kind == prev_kind:
                print('PY MoonNodes({} line {}): duplicate ascending/descending node'.format(filename, lnum))
                return 1
            time = astronomy.Time.Parse(token[1])
            ra = float(token[2])
            dec = float(token[3])
            sphere = astronomy.Spherical(dec, 15.0 * ra, 1.0)
            vec_test = astronomy.VectorFromSphere(sphere, time)

            # Calculate EQD coordinates of the Moon. Verify against input file.
            vec_eqj = astronomy.GeoMoon(time)
            rot = astronomy.Rotation_EQJ_EQD(time)
            vec_eqd = astronomy.RotateVector(rot, vec_eqj)
            angle = astronomy.AngleBetween(vec_test, vec_eqd)
            diff_angle = 60.0 * abs(angle)
            if diff_angle > max_angle:
                max_angle = diff_angle
            if diff_angle > 1.54:
                print('PY MoonNodes({} line {}): EXCESSIVE equatorial error = {:0.3f} arcmin'.format(filename, lnum, diff_angle))

            if lnum == 1:
                # The very first time, so search for the first node in the series.
                # Back up a few days to make sure we really are finding it ourselves.
                earlier = time.AddDays(-6.5472)    # back up by a weird amount of time
                node = astronomy.SearchMoonNode(earlier)
            else:
                # Use the previous node to find the next node.
                node = astronomy.NextMoonNode(node)

            # Verify the ecliptic latitude is very close to zero at the alleged node.
            ecl = astronomy.EclipticGeoMoon(node.time)
            diff_lat = 60.0 * abs(ecl.lat)
            if diff_lat > 8.1e-4:
                print('PY MoonNodes({} line {}): found node has excessive latitude = {:0.4f} arcmin.'.format(filename, lnum, diff_lat))
                return 1

            # Verify the time agrees with Espenak's time to within a few minutes.
            diff_minutes = (24.0 * 60.0) * abs(node.time.tt - time.tt)
            if diff_minutes > max_minutes:
                max_minutes = diff_minutes

            # Verify the kind of node matches what Espenak says (ascending or descending).
            if kind == 'A' and node.kind != astronomy.NodeEventKind.Ascending:
                print('PY MoonNodes({} line {}): did not find ascending node as expected.'.format(filename, lnum))
                return 1

            if kind == 'D' and node.kind != astronomy.NodeEventKind.Descending:
                print('PY MoonNodes({} line {}): did not find descending node as expected.'.format(filename, lnum))
                return 1

            prev_kind = kind
    if max_minutes > 3.681:
        print('PY MoonNodes: EXCESSIVE time prediction error = {:0.3f} minutes.'.format(max_minutes))
        return 1
    print('PY MoonNodes: PASS ({} nodes, max equ error = {:0.3f} arcmin, max time error = {:0.3f} minutes.)'.format(lnum, max_angle, max_minutes))
    return 0

#-----------------------------------------------------------------------------------------------------------

def MoonReversePhase(longitude):
    # Verify that SearchMoonPhase works both forward and backward in time.
    nphases = 5000
    utList = []
    dtMin = +1000.0
    dtMax = -1000.0

    # Search forward in time from 1800 to find consecutive phase events events.
    time = astronomy.Time.Make(1800, 1, 1, 0, 0, 0.0)
    for i in range(nphases):
        result = astronomy.SearchMoonPhase(longitude, time, +40.0)
        if result is None:
            print('PY MoonReversePhase(lon={}, i={}): failed to find event after {}'.format(longitude, i, time))
            return 1
        utList.append(result.ut)
        if i > 0:
            # Verify that consecutive events are reasonably close to the synodic period (29.5 days) apart.
            dt = v(utList[i] - utList[i-1])
            if dt < dtMin:
                dtMin = dt
            if dt > dtMax:
                dtMax = dt
        time = result.AddDays(+0.1)

    Debug('PY MoonReversePhase({}): dtMin={:0.6f} days, dtMax={:0.6f} days.'.format(longitude, dtMin, dtMax))
    if (dtMin < 29.175) or (dtMax > 29.926):
        print('PY MoonReversePhase({}): Time between consecutive events is suspicious.'.format(longitude))
        return 1

    # Do a reverse chronological search and make sure the results are consistent with the forward search.
    time = time.AddDays(20.0)
    maxDiff = 0.0
    for i in range(nphases-1, -1, -1):
        result = astronomy.SearchMoonPhase(longitude, time, -40.0)
        if result is None:
            print('PY MoonReversePhase(lon={}, i={}): failed to find event before {}'.format(longitude, i, time))
            return 1
        diff = SECONDS_PER_DAY * vabs(result.ut - utList[i])
        if diff > maxDiff:
            maxDiff = diff
        time = result.AddDays(-0.1)

    Debug('PY MoonReversePhase({}): Maximum discrepancy in reverse search = {:0.6f} seconds.'.format(longitude, maxDiff))
    if maxDiff > 0.164:
        print('PY MoonReversePhase({}): EXCESSIVE DISCREPANCY in reverse search.'.format(longitude))
        return 1

    # Pick a pair of consecutive events from the middle of the list.
    # Verify forward and backward searches work correctly from many intermediate times.
    nslots = 100
    k = nphases // 2
    ut1 = utList[k]
    ut2 = utList[k+1]
    for i in range(1, nslots):
        ut = ut1 + (i/nslots)*(ut2 - ut1)
        time = astronomy.Time(ut)
        before = astronomy.SearchMoonPhase(longitude, time, -40.0)
        if before is None:
            print('PY MoonReversePhase(lon={}, time={}): backward search failed'.format(longitude, time))
            return 1
        diff = SECONDS_PER_DAY * vabs(before.ut - ut1)
        if diff > 0.07:
            print('PY MoonReversePhase(lon={}, time={}): backward search error = {:0.4e} seconds'.format(longitude, time, diff))
            return 1
        after = astronomy.SearchMoonPhase(longitude, time, +40.0)
        if after is None:
            print('PY MoonReversePhase(lon={}, time={}): forward search failed'.format(longitude, time))
            return 1
        diff = SECONDS_PER_DAY * vabs(after.ut - ut2)
        if diff > 0.07:
            print('PY MoonReversePhase(lon={}, time={}): forward search error = {:0.4e} seconds'.format(longitude, time, diff))
            return 1

    print('PY MoonReversePhase({}): PASS'.format(longitude))
    return 0


def MoonReverse():
    return (
        MoonReversePhase(0.0) or
        MoonReversePhase(90.0) or
        MoonReversePhase(180.0) or
        MoonReversePhase(270.0)
    )

#-----------------------------------------------------------------------------------------------------------

class LagrangeFunc:
    def __init__(self, point, major_body, minor_body):
        self.point = point
        self.major_body = major_body
        self.minor_body = minor_body

    def Eval(self, time):
        return astronomy.LagrangePoint(self.point, time, self.major_body, self.minor_body)


def VerifyStateLagrange(major_body, minor_body, point, filename, r_thresh, v_thresh):
    func = LagrangeFunc(point, major_body, minor_body)
    return VerifyStateBody(func, filename, r_thresh, v_thresh)


def Lagrange():
    # Test Sun/EMB Lagrange points.
    if VerifyStateLagrange(astronomy.Body.Sun, astronomy.Body.EMB, 1, 'lagrange/semb_L1.txt',   1.33e-5, 6.13e-5): return 1
    if VerifyStateLagrange(astronomy.Body.Sun, astronomy.Body.EMB, 2, 'lagrange/semb_L2.txt',   1.33e-5, 6.13e-5): return 1
    if VerifyStateLagrange(astronomy.Body.Sun, astronomy.Body.EMB, 4, 'lagrange/semb_L4.txt',   3.75e-5, 5.28e-5): return 1
    if VerifyStateLagrange(astronomy.Body.Sun, astronomy.Body.EMB, 5, 'lagrange/semb_L5.txt',   3.75e-5, 5.28e-5): return 1

    # Test Earth/Moon Lagrange points.
    if VerifyStateLagrange(astronomy.Body.Earth, astronomy.Body.Moon, 1, 'lagrange/em_L1.txt',  3.79e-5, 5.06e-5): return 1
    if VerifyStateLagrange(astronomy.Body.Earth, astronomy.Body.Moon, 2, 'lagrange/em_L2.txt',  3.79e-5, 5.06e-5): return 1
    if VerifyStateLagrange(astronomy.Body.Earth, astronomy.Body.Moon, 4, 'lagrange/em_L4.txt',  3.79e-5, 1.59e-3): return 1
    if VerifyStateLagrange(astronomy.Body.Earth, astronomy.Body.Moon, 5, 'lagrange/em_L5.txt',  3.79e-5, 1.59e-3): return 1

    print('PY Lagrange: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def SiderealTime():
    correct = 9.3983699280076483
    time = astronomy.Time.Make(2022, 3, 15, 21, 50, 0)
    gast = astronomy.SiderealTime(time)
    diff = abs(gast - correct)
    print('PY SiderealTime: gast={:0.15f}, correct={:0.15f}, diff={:0.3e}'.format(gast, correct, diff))
    if diff > 1.0e-15:
        print('PY SiderealTime: EXCESSIVE ERROR')
        return 1
    print('PY SiderealTime: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def Repr():
    time = astronomy.Time.Make(2022, 3, 31, 21, 4, 45.123)
    if str(time) != '2022-03-31T21:04:45.123Z':
        print('PY Repr: FAIL str(time)')
        return 1

    if repr(time) != "Time('2022-03-31T21:04:45.123Z')":
        print('PY Repr: FAIL repr(time)')
        return 1

    vec = astronomy.Vector(-1.8439088914585775, 1.51657508881031, 0.8366600265340756, time)
    if repr(vec) != "Vector(-1.8439088914585775, 1.51657508881031, 0.8366600265340756, Time('2022-03-31T21:04:45.123Z'))":
        print('PY Repr: FAIL repr(vec)')
        return 1

    state = astronomy.StateVector(vec.x, vec.y, vec.z, -vec.x/3, -vec.y/3, -vec.z/3, vec.t)
    if repr(state) != "StateVector(x=-1.8439088914585775, y=1.51657508881031, z=0.8366600265340756, vx=0.6146362971528592, vy=-0.5055250296034367, vz=-0.27888667551135854, t=Time('2022-03-31T21:04:45.123Z'))":
        print('PY Repr: FAIL repr(state)')
        return 1

    observer = astronomy.Observer(32.1, 45.6, 98.765)
    if repr(observer) != 'Observer(latitude=32.1, longitude=45.6, height=98.765)':
        print('PY Repr: FAIL repr(observer)')
        return 1

    rot = astronomy.Rotation_EQJ_ECL()
    if repr(rot) != 'RotationMatrix([[1, 0, 0], [0, 0.9174821430670688, -0.3977769691083922], [0, 0.3977769691083922, 0.9174821430670688]])':
        print('PY Repr: FAIL repr(rot)')
        return 1

    sph = astronomy.Spherical(lat=-27.3, lon=85.2, dist=2.54)
    if repr(sph) != 'Spherical(lat=-27.3, lon=85.2, dist=2.54)':
        print('PY Repr: FAIL repr(sph)')
        return 1

    equ = astronomy.Equatorial(8.54, -23.753, 2.986, vec)
    if repr(equ) != "Equatorial(ra=8.54, dec=-23.753, dist=2.986, vec=Vector(-1.8439088914585775, 1.51657508881031, 0.8366600265340756, Time('2022-03-31T21:04:45.123Z')))":
        print('PY Repr: FAIL repr(equ)')
        return 1

    print('PY Repr: PASS')
    return 0

#-----------------------------------------------------------------------------------------------------------

def GravSimTest():
    Debug("")

    if 0 != GravSimEmpty("barystate/Sun.txt",      astronomy.Body.SSB, astronomy.Body.Sun,      0.0269, 1.9635): return 1
    if 0 != GravSimEmpty("barystate/Mercury.txt",  astronomy.Body.SSB, astronomy.Body.Mercury,  0.5725, 0.9332): return 1
    if 0 != GravSimEmpty("barystate/Venus.txt",    astronomy.Body.SSB, astronomy.Body.Venus,    0.1433, 0.1458): return 1
    if 0 != GravSimEmpty("barystate/Earth.txt",    astronomy.Body.SSB, astronomy.Body.Earth,    0.0651, 0.2098): return 1
    if 0 != GravSimEmpty("barystate/Mars.txt",     astronomy.Body.SSB, astronomy.Body.Mars,     0.1150, 0.1896): return 1
    if 0 != GravSimEmpty("barystate/Jupiter.txt",  astronomy.Body.SSB, astronomy.Body.Jupiter,  0.2546, 0.8831): return 1
    if 0 != GravSimEmpty("barystate/Saturn.txt",   astronomy.Body.SSB, astronomy.Body.Saturn,   0.3660, 1.0818): return 1
    if 0 != GravSimEmpty("barystate/Uranus.txt",   astronomy.Body.SSB, astronomy.Body.Uranus,   0.3107, 0.9321): return 1
    if 0 != GravSimEmpty("barystate/Neptune.txt",  astronomy.Body.SSB, astronomy.Body.Neptune,  0.3382, 1.5586): return 1
    if 0 != GravSimEmpty("heliostate/Mercury.txt", astronomy.Body.Sun, astronomy.Body.Mercury,  0.5087, 0.9473): return 1
    if 0 != GravSimEmpty("heliostate/Venus.txt",   astronomy.Body.Sun, astronomy.Body.Venus,    0.1214, 0.1543): return 1
    if 0 != GravSimEmpty("heliostate/Earth.txt",   astronomy.Body.Sun, astronomy.Body.Earth,    0.0508, 0.2099): return 1
    if 0 != GravSimEmpty("heliostate/Mars.txt",    astronomy.Body.Sun, astronomy.Body.Mars,     0.1085, 0.1927): return 1
    if 0 != GravSimEmpty("heliostate/Jupiter.txt", astronomy.Body.Sun, astronomy.Body.Jupiter,  0.2564, 0.8805): return 1
    if 0 != GravSimEmpty("heliostate/Saturn.txt",  astronomy.Body.Sun, astronomy.Body.Saturn,   0.3664, 1.0826): return 1
    if 0 != GravSimEmpty("heliostate/Uranus.txt",  astronomy.Body.Sun, astronomy.Body.Uranus,   0.3106, 0.9322): return 1
    if 0 != GravSimEmpty("heliostate/Neptune.txt", astronomy.Body.Sun, astronomy.Body.Neptune,  0.3381, 1.5584): return 1

    Debug("")
    nsteps = 20

    if 0 != GravSimFile("barystate/Ceres.txt",    astronomy.Body.SSB,   nsteps, 0.6640, 0.6226): return 1
    if 0 != GravSimFile("barystate/Pallas.txt",   astronomy.Body.SSB,   nsteps, 0.4687, 0.3474): return 1
    if 0 != GravSimFile("barystate/Vesta.txt",    astronomy.Body.SSB,   nsteps, 0.5806, 0.5462): return 1
    if 0 != GravSimFile("barystate/Juno.txt",     astronomy.Body.SSB,   nsteps, 0.6760, 0.5750): return 1
    if 0 != GravSimFile("barystate/Bennu.txt",    astronomy.Body.SSB,   nsteps, 3.7444, 2.6581): return 1
    if 0 != GravSimFile("barystate/Halley.txt",   astronomy.Body.SSB,   nsteps, 0.0539, 0.0825): return 1
    if 0 != GravSimFile("heliostate/Ceres.txt",   astronomy.Body.Sun,   nsteps, 0.0445, 0.0355): return 1
    if 0 != GravSimFile("heliostate/Pallas.txt",  astronomy.Body.Sun,   nsteps, 0.1062, 0.0854): return 1
    if 0 != GravSimFile("heliostate/Vesta.txt",   astronomy.Body.Sun,   nsteps, 0.1432, 0.1308): return 1
    if 0 != GravSimFile("heliostate/Juno.txt",    astronomy.Body.Sun,   nsteps, 0.1554, 0.1328): return 1
    if 0 != GravSimFile("geostate/Ceres.txt",     astronomy.Body.Earth, nsteps, 6.5689, 6.4797): return 1
    if 0 != GravSimFile("geostate/Pallas.txt",    astronomy.Body.Earth, nsteps, 9.3288, 7.3533): return 1
    if 0 != GravSimFile("geostate/Vesta.txt",     astronomy.Body.Earth, nsteps, 3.2980, 3.8863): return 1
    if 0 != GravSimFile("geostate/Juno.txt",      astronomy.Body.Earth, nsteps, 6.0962, 7.7147): return 1

    Debug("")
    print("PY GravSimTest: PASS")

    return 0


def GravSimEmpty(filename, origin, body, rthresh, vthresh):
    max_rdiff = 0.0
    max_vdiff = 0.0
    sim = None
    for rec in JplHorizonsStateVectors(filename):
        if sim is None:
            sim = astronomy.GravitySimulator(origin, rec.state.t, [])
        sim.Update(rec.state.t)
        calc = sim.SolarSystemBodyState(body)
        if origin == astronomy.Body.SSB and body == astronomy.Body.Sun:
            rdiff = SsbArcminPosError(rec.state, calc)
        else:
            rdiff = ArcminPosError(rec.state, calc)
        if rdiff > rthresh:
            print('PY GravSimEmpty({} line {}): excessive position error = {} arcmin.'.format(filename, rec.lnum, rdiff))
            return 1
        if rdiff > max_rdiff:
            max_rdiff = rdiff

        vdiff = ArcminVelError(rec.state, calc)
        if vdiff > vthresh:
            print('PY GravSimEmpty({} line {}): excessive velocity error = {} arcmin.'.format(filename, rec.lnum, vdiff))
            return 1
        if vdiff > max_vdiff:
            max_vdiff = vdiff
    Debug('PY GravSimEmpty({:22s}): PASS - max pos error = {:0.4f} arcmin, max vel error = {:0.4f} arcmin.'.format(filename, max_rdiff, max_vdiff))
    return 0


def GravSimFile(filename, originBody, nsteps, rthresh, vthresh):
    sim = None
    max_rdiff = 0.0
    max_vdiff = 0.0
    for rec in JplHorizonsStateVectors(filename):
        if sim is None:
            sim = astronomy.GravitySimulator(originBody, rec.state.t, [rec.state])
            time = rec.state.t
            smallBodyArray = sim.Update(time)
        else:
            tt1 = prev.state.t.tt
            tt2 = rec.state.t.tt
            dt = (tt2 - tt1) / nsteps
            for k in range(1, nsteps+1):
                time = astronomy.Time.FromTerrestrialTime(tt1 + k*dt)
                smallBodyArray = sim.Update(time)
                if len(smallBodyArray) != 1:
                    print('PY GravSimFile({} line {}): unexpected smallBodyArray.length = {}'.format(filename, rec.lnum, len(smallBodyArray)))
                    return 1
                if time.tt != sim.GetTime().tt:
                    print('PY GravSimFile({} line {}): expected {} but simulator reports {}'.format(filename, rec.lnum, time, sim.GetTime()))
                    return 1
        rdiff = ArcminPosError(rec.state, smallBodyArray[0])
        if rdiff > rthresh:
            print('PY GravSimFile({} line {}): excessive position error = {}'.format(filename, rec.lnum, rdiff))
            return 1
        if rdiff > max_rdiff:
            max_rdiff = rdiff
        vdiff = ArcminVelError(rec.state, smallBodyArray[0])
        if vdiff > vthresh:
            print('PY GravSimFile({} line {}): excessive position error = {}'.format(filename, rec.lnum, vdiff))
            return 1
        if vdiff > max_vdiff:
            max_vdiff = vdiff
        prev = rec
    Debug('PY GravSimFile({:22s}): PASS - max pos error = {:0.4f} arcmin, max vel error = {:0.4f} arcmin.'.format(filename, max_rdiff, max_vdiff))
    return 0


def SsbArcminPosError(correct, calc):
    # Scale the SSB based on 1 AU, not on its absolute magnitude, which can become very close to zero.
    dx = calc.x - correct.x
    dy = calc.y - correct.y
    dz = calc.z - correct.z
    diffSquared = dx*dx + dy*dy + dz*dz
    radians = sqrt(diffSquared)
    return 60.0 * math.degrees(radians)


def ArcminPosError(correct, calc):
    dx = calc.x - correct.x
    dy = calc.y - correct.y
    dz = calc.z - correct.z
    diffSquared = dx*dx + dy*dy + dz*dz
    magSquared = correct.x*correct.x + correct.y*correct.y + correct.z*correct.z
    radians = sqrt(diffSquared / magSquared)
    return 60.0 * math.degrees(radians)


def ArcminVelError(correct, calc):
    dx = calc.vx - correct.vx
    dy = calc.vy - correct.vy
    dz = calc.vz - correct.vz
    diffSquared = dx*dx + dy*dy + dz*dz
    magSquared = correct.vx*correct.vx + correct.vy*correct.vy + correct.vz*correct.vz
    radians = sqrt(diffSquared / magSquared)
    return 60.0 * math.degrees(radians)


#-----------------------------------------------------------------------------------------------------------

def CheckDecemberSolstice(year, expected):
    si = astronomy.Seasons(year)
    actual = str(si.dec_solstice)
    if actual != expected:
        print('PY DatesIssue250: FAIL: year {}, expected [{}], actual [{}]'.format(year, expected, actual))
        return 1
    return 0

def DatesIssue250():
    # Make sure we can handle dates outside the range supported by System.DateTime.
    # https://github.com/cosinekitty/astronomy/issues/250
    return (
        CheckDecemberSolstice( 2022, "2022-12-21T21:47:54.455Z") or
        CheckDecemberSolstice(-2300, "-002300-12-19T16:22:27.929Z") or
        CheckDecemberSolstice(12345, "+012345-12-11T13:30:10.276Z") or
        Pass('DatesIssue250')
    )

#-----------------------------------------------------------------------------------------------------------

def LunarFractionCase(year, month, day, obscuration):
    time = astronomy.Time.Make(year, month, day, 0, 0, 0.0)
    eclipse = astronomy.SearchLunarEclipse(time)
    # This should be a partial lunar eclipse.
    if eclipse.kind != astronomy.EclipseKind.Partial:
        print('PY LunarFractionCase({:04d}-{:02d}-{:02d}) FAIL: expected partial eclipse, but found {}.'.format(year, month, day, eclipse.kind))
        return 1

    # The partial eclipse should always happen within 24 hours of the given date.
    dt = v(eclipse.peak.ut - time.ut)
    if dt < 0.0 or dt > 1.0:
        print('PY LunarFractionCase({:04d}-{:02d}-{:02d}) FAIL: eclipse occurs {:0.4f} days after predicted date.'.format(year, month, day, dt))
        return 1

    diff = v(eclipse.obscuration - obscuration)
    if abs(diff) > 0.00763:
        print('PY LunarFractionCase({:04d}-{:02d}-{:02d}) FAIL: excessive obscuration diff = {:0.8f}, expected = {:0.8f}, actual = {:0.8f}'.format(year, month, day, diff, obscuration, eclipse.obscuration))
        return 1

    Debug('PY LunarFractionCase({:04d}-{:02d}-{:02d}): obscuration diff = {:11.8f}'.format(year, month, day, diff))
    return 0


def LunarFraction():
    # Verify calculation of the fraction of the Moon's disc covered by the Earth's umbra during a partial eclipse.
    # Data for this is more tedious to gather, because Espenak data does not contain it.
    # We already verify fraction=0.0 for penumbral eclipses and fraction=1.0 for total eclipses in LunarEclipseTest.
    return (
        LunarFractionCase(2010,  6, 26, 0.506) or  # https://www.timeanddate.com/eclipse/lunar/2010-june-26
        LunarFractionCase(2012,  6,  4, 0.304) or  # https://www.timeanddate.com/eclipse/lunar/2012-june-4
        LunarFractionCase(2013,  4, 25, 0.003) or  # https://www.timeanddate.com/eclipse/lunar/2013-april-25
        LunarFractionCase(2017,  8,  7, 0.169) or  # https://www.timeanddate.com/eclipse/lunar/2017-august-7
        LunarFractionCase(2019,  7, 16, 0.654) or  # https://www.timeanddate.com/eclipse/lunar/2019-july-16
        LunarFractionCase(2021, 11, 19, 0.991) or  # https://www.timeanddate.com/eclipse/lunar/2021-november-19
        LunarFractionCase(2023, 10, 28, 0.060) or  # https://www.timeanddate.com/eclipse/lunar/2023-october-28
        LunarFractionCase(2024,  9, 18, 0.035) or  # https://www.timeanddate.com/eclipse/lunar/2024-september-18
        LunarFractionCase(2026,  8, 28, 0.962) or  # https://www.timeanddate.com/eclipse/lunar/2026-august-28
        LunarFractionCase(2028,  1, 12, 0.024) or  # https://www.timeanddate.com/eclipse/lunar/2028-january-12
        LunarFractionCase(2028,  7,  6, 0.325) or  # https://www.timeanddate.com/eclipse/lunar/2028-july-6
        LunarFractionCase(2030,  6, 15, 0.464) or  # https://www.timeanddate.com/eclipse/lunar/2030-june-15
        Pass('LunarFraction')
    )

#-----------------------------------------------------------------------------------------------------------

def StarRiseSetCulmCase(starName, ra, dec, distLy, observer, year, month, day, riseHour, riseMinute, culmHour, culmMinute, setHour, setMinute):
    func = 'StarRiseSetCulmCase({})'.format(starName)

    # Calculate expected event times.
    expectedRiseTime = astronomy.Time.Make(year, month, day, riseHour, riseMinute, 0.0)
    expectedCulmTime = astronomy.Time.Make(year, month, day, culmHour, culmMinute, 0.0)
    expectedSetTime  = astronomy.Time.Make(year, month, day, setHour,  setMinute,  0.0)

    # Define a custom star object.
    astronomy.DefineStar(astronomy.Body.Star1, ra, dec, distLy)

    # Use Astronomy Engine to search for event times.
    searchTime = astronomy.Time.Make(year, month, day, 0, 0, 0.0)

    rise = astronomy.SearchRiseSet(astronomy.Body.Star1, observer, astronomy.Direction.Rise, searchTime, 1.0)
    if rise is None:
        return Fail(func, 'Star rise search failed.')

    culm = astronomy.SearchHourAngle(astronomy.Body.Star1, observer, 0.0, searchTime, +1)
    if culm is None:
        return Fail(func, 'Star culmination search failed.')

    set = astronomy.SearchRiseSet(astronomy.Body.Star1, observer, astronomy.Direction.Set, searchTime, 1.0)
    if set is None:
        return Fail(func, 'Star set search failed.')

    # Compare expected times with calculated times.
    rdiff = MINUTES_PER_DAY * vabs(expectedRiseTime.ut - rise.ut)
    cdiff = MINUTES_PER_DAY * vabs(expectedCulmTime.ut - culm.time.ut)
    sdiff = MINUTES_PER_DAY * vabs(expectedSetTime.ut - set.ut)

    Debug("{}: minutes rdiff = {:0.4f}, cdiff = {:0.4f}, sdiff = {:0.4f}".format(func, rdiff, cdiff, sdiff))

    if rdiff > 0.5: return Fail(func, "excessive rise time error = {:0.4f} minutes.".format(rdiff))
    if cdiff > 0.5: return Fail(func, "excessive culm time error = {:0.4f} minutes.".format(cdiff))
    if sdiff > 0.5: return Fail(func, "excessive set time error =  {:0.4f} minutes.".format(sdiff))
    return 0


def StarRiseSetCulm():
    observer = astronomy.Observer(+25.77, -80.19, 0.0)
    return (
        StarRiseSetCulmCase("Sirius",   6.7525, -16.7183,   8.6, observer, 2022, 11, 21,  2, 37,  8,  6, 13, 34) or
        StarRiseSetCulmCase("Sirius",   6.7525, -16.7183,   8.6, observer, 2022, 11, 25,  2, 22,  7, 50, 13, 18) or
        StarRiseSetCulmCase("Canopus",  6.3992, -52.6956, 310.0, observer, 2022, 11, 21,  4, 17,  7, 44, 11, 11) or
        StarRiseSetCulmCase("Canopus",  6.3992, -52.6956, 310.0, observer, 2022, 11, 25,  4,  1,  7, 28, 10, 56) or
        Pass("StarRiseSetCulm")
    )

#-----------------------------------------------------------------------------------------------------------

class HourAngleTester:
    def __init__(self):
        self.cases = 0
        self.maxdiff = 0.0

    def Case(self, latitude, longitude, hourAngle):
        threshold = 0.1 / 3600   # SearchHourAngle() accuracy: 0.1 seconds converted to hours
        observer = astronomy.Observer(latitude, longitude, 0)
        startTime = astronomy.Time.Make(2023, 2, 11, 0, 0, 0)
        search = astronomy.SearchHourAngle(astronomy.Body.Sun, observer, hourAngle, startTime, +1)
        calc = astronomy.HourAngle(astronomy.Body.Sun, search.time, observer)
        diff = vabs(calc - hourAngle)
        if diff > 12.0:
            diff = 24.0 - diff;
        if diff > self.maxdiff:
            self.maxdiff = diff
        self.cases += 1
        if diff > threshold:
            print('PY HourAngleCase: EXCESSIVE ERROR = {:0.6e}, calc HA = {:0.16f}, for hourAngle={:0.1f}'.format(diff, calc, hourAngle))
            return False
        Debug('PY HourAngleCase: Hour angle = {:4.1f}, longitude = {:6.1f}, diff = {:9.4e}'.format(hourAngle, longitude, diff))
        return True

    def Pass(self):
        print('PY HourAngle ({:d} cases, maxdiff = {:9.4e}): PASS'.format(self.cases, self.maxdiff))
        return 0


def HourAngle():
    tester = HourAngleTester()
    latitude = 35
    longitude = -170
    while longitude <= 180:
        hour = 0
        while hour < 24:
            if not tester.Case(latitude, longitude, hour):
                return 1
            hour += 1
        longitude += 5
    return tester.Pass()

#-----------------------------------------------------------------------------------------------------------

def Atmosphere():
    filename = 'riseset/atmosphere.csv'
    maxdiff = 0.0
    ncases = 0
    tolerance = 8.8e-11
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            line = line.strip()
            lnum += 1
            if lnum == 1:
                if line != 'elevation,temperature,pressure,density,relative_density':
                    return Fail('Atmosphere', 'Expected header line but found [{}]'.format(line))
            else:
                tokens = line.split(',')
                if len(tokens) != 5:
                    return Fail('Atmosphere({} line {})'.format(filename, lnum), 'expected 5 numeric tokens but found {}'.format(len(tokens)))
                elevation = v(float(tokens[0]))
                temperature = v(float(tokens[1]))
                pressure = v(float(tokens[2]))
                # ignore tokens[3] = absolute_density
                relative_density = v(float(tokens[4]))

                atmos = astronomy.Atmosphere(elevation)

                diff = vabs(atmos.temperature - temperature)
                maxdiff = max(maxdiff, diff)
                if diff > tolerance:
                    return Fail('Atmosphere', 'EXCESSIVE temperature difference = {}'.format(diff))

                diff = vabs(atmos.pressure - pressure)
                maxdiff = max(maxdiff, diff)
                if diff > tolerance:
                    return Fail('Atmosphere', 'EXCESSIVE pressure difference = {}'.format(diff))

                diff = vabs(atmos.density - relative_density)
                maxdiff = max(maxdiff, diff)
                if diff > tolerance:
                    return Fail('Atmosphere', 'EXCESSIVE density difference = {}'.format(diff))

                ncases += 1

    if ncases != 34:
        return Fail('Atmosphere', 'expected 34 cases but found {}'.format(ncases))

    return Pass('Atmosphere')

#-----------------------------------------------------------------------------------------------------------

def RiseSetElevationBodyCase(body, observer, direction, metersAboveGround, startTime, eventOffsetDays):
    time = astronomy.SearchRiseSet(body, observer, direction, startTime, 2.0, metersAboveGround)
    if not time:
        return Fail('RiseSetElevationBodyCase {} {}: search failed.'.format(body, direction))
    diff = v(time.ut - (startTime.ut + eventOffsetDays))
    if diff > 0.5:
        diff -= 1.0     # assume event actually takes place on the next day
    diff = vabs(MINUTES_PER_DAY * diff)     # convert signed days to absolute minutes
    if diff > 0.5:
        return Fail('RiseSetElevationBodyCase {} {}: EXCESSIVE diff = {}.'.format(body, direction, diff))
    return 0


def RiseSetElevation():
    regex = re.compile(r'^(\d+)-(\d+)-(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+):(\d+)\s+(\d+):(\d+)\s+(\d+):(\d+)\s+(\d+):(\d+)\s+(\S+)\s*$')
    filename = 'riseset/elevation.txt'
    with open(filename, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            if line.startswith('#'):
                continue
            m = regex.match(line)
            if not m:
                return Fail('RiseSetElevation({} line {})'.format(filename, lnum), 'Invalid data format')
            year = int(m.group(1))
            month = int(m.group(2))
            day = int(m.group(3))
            latitude = v(float(m.group(4)))
            longitude = v(float(m.group(5)))
            height = v(float(m.group(6)))
            metersAboveGround = v(float(m.group(7)))
            srh = int(m.group( 8))
            srm = int(m.group( 9))
            ssh = int(m.group(10))
            ssm = int(m.group(11))
            mrh = int(m.group(12))
            mrm = int(m.group(13))
            msh = int(m.group(14))
            msm = int(m.group(15))

            # Get search origin time
            time = astronomy.Time.Make(year, month, day, 0, 0, 0.0)

            # Convert scanned values into sunrise, sunset, moonrise, moonset day offsets.
            sr = (srh + srm/60.0) / 24.0
            ss = (ssh + ssm/60.0) / 24.0
            mr = (mrh + mrm/60.0) / 24.0
            ms = (msh + msm/60.0) / 24.0

            observer = astronomy.Observer(latitude, longitude, height)

            if (0 != RiseSetElevationBodyCase(astronomy.Body.Sun,  observer, astronomy.Direction.Rise, metersAboveGround, time, sr) or
                0 != RiseSetElevationBodyCase(astronomy.Body.Sun,  observer, astronomy.Direction.Set,  metersAboveGround, time, ss) or
                0 != RiseSetElevationBodyCase(astronomy.Body.Moon, observer, astronomy.Direction.Rise, metersAboveGround, time, mr) or
                0 != RiseSetElevationBodyCase(astronomy.Body.Moon, observer, astronomy.Direction.Set,  metersAboveGround, time, ms)):
                return 1


    return Pass('RiseSetElevation')

#-----------------------------------------------------------------------------------------------------------

def Trigonometry() -> int:
    tolerance = 1.0e-15
    inFileName = 'trigonometry/trig.txt'
    with open(inFileName, 'rt') as infile:
        lnum = 0
        for line in infile:
            lnum += 1
            line = line.strip()
            if line == '':
                break
            token = [float(s) for s in line.split()]
            print(token)
            if len(token) != 3:
                return Fail('Trigonometry({} line {}): incorrect number of tokens in [{}]'.format(inFileName, lnum, line))
            (deg, cosCorrect, sinCorrect) = token
            rad = math.radians(deg)
            cosCalc = xcos(rad)
            sinCalc = xsin(rad)
            cosDiff = abs(cosCalc - cosCorrect)
            sinDiff = abs(sinCalc - sinCorrect)
            if (cosDiff > tolerance) or (sinDiff > tolerance):
                return Fail(
                    'Trigonometry({:s} line {:d})'.format(inFileName, lnum),
                    'EXCESS ERROR - deg={:0.1f}, cosDiff={:g}, sinDiff={:g}'.format(deg, cosDiff, sinDiff)
                )
    return Pass('Trigonometry')

#-----------------------------------------------------------------------------------------------------------


UnitTests = {
    '_trig':                    Trigonometry,
    'aberration':               Aberration,
    'atmosphere':               Atmosphere,
    'axis':                     Axis,
    'barystate':                BaryState,
    'constellation':            Constellation,
    'dates250':                 DatesIssue250,
    'ecliptic':                 Ecliptic,
    'elongation':               Elongation,
    'geoid':                    Geoid,
    'global_solar_eclipse':     GlobalSolarEclipse,
    'gravsim':                  GravSimTest,
    'heliostate':               HelioState,
    'hour_angle':               HourAngle,
    'issue_103':                Issue103,
    'jupiter_moons':            JupiterMoons,
    'lagrange':                 Lagrange,
    'libration':                Libration,
    'local_solar_eclipse':      LocalSolarEclipse,
    'lunar_apsis':              LunarApsis,
    'lunar_eclipse':            LunarEclipse,
    'lunar_eclipse_78':         LunarEclipseIssue78,
    'lunar_fraction':           LunarFraction,
    'magnitude':                Magnitude,
    'moon':                     GeoMoon,
    'moon_nodes':               MoonNodes,
    'moon_reverse':             MoonReverse,
    'moonphase':                MoonPhase,
    'planet_apsis':             PlanetApsis,
    'pluto':                    PlutoCheck,
    'precision':                PrecisionProblem,
    'refraction':               Refraction,
    'repr':                     Repr,
    'riseset':                  RiseSet,
    'riseset_elevation':        RiseSetElevation,
    'riseset_reverse':          RiseSetReverse,
    'rotation':                 Rotation,
    'seasons':                  Seasons,
    'seasons187':               SeasonsIssue187,
    'sidereal':                 SiderealTime,
    'solar_fraction':           SolarFraction,
    'star_risesetculm':         StarRiseSetCulm,
    'time':                     AstroTime,
    'topostate':                TopoState,
    'transit':                  Transit,
    'twilight':                 Twilight,
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
