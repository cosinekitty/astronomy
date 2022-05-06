package io.github.cosinekitty.astronomy.demo;
import io.github.cosinekitty.astronomy.*;

public class Constellation {
    /**
     * Find current constellation of the Moon, and its constellation changes over the next 30 days.
     *
     * @param startTime
     * The date and time to start the constellation search.
     */
    public static int run(Time startTime) {
        Time stopTime = startTime.addDays(30.0);

        // There are 12 zodiac constellations, and the moon takes
        // about 27.3 days in its sidereal period. Therefore, there
        // are roughly 2.2 days per constellation. We will sample
        // the Moon's constellation once every 0.1 days to reduce
        // the chance of missing a brief transition through a small
        // part of a constellation.
        final double dayIncrement = 0.1;
        Time t1 = startTime;
        ConstellationInfo c1 = bodyConstellation(Body.Moon, t1);
        System.out.printf(String.format("%s : The Moon starts in %s.%n", t1, c1.getName()));

        while (t1.getTt() < stopTime.getTt()) {
            Time t2 = t1.addDays(dayIncrement);
            ConstellationInfo c2 = bodyConstellation(Body.Moon, t2);
            if (sameConstellation(c1, c2)) {
                // No constellation change in this time step.
                // Try again on the next time step.
                c1 = c2;
                t1 = t2;
            } else {
                // The body moved from one constellation to another during this time step.
                // Narrow in on the exact moment of transition by doing a binary search.
                ConstellationEvent change = findConstellationChange(Body.Moon, c1, t1, t2);
                System.out.printf(String.format("%s : The Moon enters %s.%n", change.time, change.con.getName()));
                c1 = change.con;
                t1 = change.time;
            }
        }

        return 0;
    }

    private static ConstellationInfo bodyConstellation(Body body, Time time) {
        // Find a vector from the center of the Earth to the center of the body.
        Vector vec = Astronomy.geoVector(body, time, Aberration.Corrected);

        // Convert cartesian vector to spherical angular coordinates.
        Equatorial equ = vec.toEquatorial();

        // Use the right ascension and declination to find the constellation.
        return Astronomy.constellation(equ.getRa(), equ.getDec());
    }

    private static boolean sameConstellation(ConstellationInfo c1, ConstellationInfo c2) {
        return c1.getSymbol().equals(c2.getSymbol());
    }

    private static class ConstellationEvent {
        public final Time time;
        public final ConstellationInfo con;

        public ConstellationEvent(Time time, ConstellationInfo con) {
            this.time = time;
            this.con = con;
        }
    }

    private static ConstellationEvent findConstellationChange(
        Body body,
        ConstellationInfo c1,
        Time startTime,
        Time endTime
    ) {
        final double tolerance = 0.1 / Astronomy.SECONDS_PER_DAY;   // one tenth of a second, expressed in days
        Time t1 = startTime;
        Time t2 = endTime;

        // Do a binary search for when the constellation changes.
        while (true) {
            // Calculate the width of the search window.
            double dt = t2.getUt() - t1.getUt();

            // Let tx = the time halfway between t1 and t2.
            Time tx = t1.addDays(dt / 2.0);

            // What constellation is the body in at tx?
            ConstellationInfo cx = bodyConstellation(body, tx);

            // Is it in the same constellation as at time t1, or a different one?
            if (sameConstellation(cx, c1)) {
                // Still in the same constellation.
                // Narrow the search window to [tx, t2].
                t1 = tx;
            } else {
                if (dt < tolerance) {
                    // We found the time of transition to within 0.1 seconds.
                    // Always end the search inside the new constellation.
                    return new ConstellationEvent(tx, cx);
                }
                // The constellation changed some time in the range [t1, tx].
                // Keep searching.
                t2 = tx;
            }
        }
    }
}
