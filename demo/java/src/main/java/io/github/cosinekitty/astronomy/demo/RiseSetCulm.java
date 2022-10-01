package io.github.cosinekitty.astronomy.demo;

import java.util.ArrayList;
import java.util.Collections;

import io.github.cosinekitty.astronomy.*;

public class RiseSetCulm {
    /**
     * Display rise/set/culmination times for the Sun and Moon.
     *
     * @param observer
     * The geographic location of the observer.
     *
     * @param startTime
     * The date and time after which to search for the first
     * occurrences of rise/set/culmination times of the Sun and Moon.
     */
    public static int run(Observer observer, Time startTime) {
        // At extreme latitudes (close to the poles), there can be long
        // periods of time where the Sun/Moon never rise/set.
        // Allow the search to proceed for 2 days, which will usually find events
        // for most populated latitudes.
        // In other applications (e.g. a daily calendar),
        // the developer may wish to use a 1-day search limit.
        // Note that the Moon rises/sets about 1 hour later each day, so it is
        // common for the Moon to not have a rise or set event in a given day.
        final double dayLimit = 2.0;

        var eventList = new ArrayList<AstroEvent>();

        // Rise/set times may or may not occur within any finite search window.
        maybeAdd(eventList, "sunrise",  Astronomy.searchRiseSet(Body.Sun,  observer, Direction.Rise, startTime, dayLimit));
        maybeAdd(eventList, "sunset",   Astronomy.searchRiseSet(Body.Sun,  observer, Direction.Set,  startTime, dayLimit));
        maybeAdd(eventList, "moonrise", Astronomy.searchRiseSet(Body.Moon, observer, Direction.Rise, startTime, dayLimit));
        maybeAdd(eventList, "moonset",  Astronomy.searchRiseSet(Body.Moon, observer, Direction.Set,  startTime, dayLimit));

        // Culmination times can always be found regardless of latitude.
        eventList.add(new AstroEvent("sunculm",  Astronomy.searchHourAngle(Body.Sun,  observer, 0.0, startTime, +1)));
        eventList.add(new AstroEvent("moonculm", Astronomy.searchHourAngle(Body.Moon, observer, 0.0, startTime, +1)));

        // Sort the list chronologically.
        Collections.sort(eventList);

        // Print the events.
        for (AstroEvent e : eventList) {
            e.display();
        }

        return 0;
    }

    private static void maybeAdd(ArrayList<AstroEvent> eventList, String name, Time time) {
        if (time != null) {
            eventList.add(new AstroEvent(name, time));
        }
    }

    private static class AstroEvent implements Comparable<AstroEvent> {
        public final String name;
        public final Time time;
        public final double altitude;

        public AstroEvent(String name, Time time) {
            this.name = name;
            this.time = time;
            this.altitude = Double.NaN;     // this event has no altitude data associated with it
        }

        public AstroEvent(String name, HourAngleInfo hourAngleInfo) {
            this.name = name;
            this.time = hourAngleInfo.getTime();
            this.altitude = hourAngleInfo.getHor().getAltitude();
        }

        @Override
        public int compareTo(AstroEvent other) {
            return this.time.compareTo(other.time);
        }

        public void display() {
            System.out.print(String.format("%-8s : %s", name, time));
            if (Double.isFinite(altitude)) {
                System.out.print(String.format("   altitude = %5.2f", altitude));
            }
            System.out.println();
        }
    }
}
