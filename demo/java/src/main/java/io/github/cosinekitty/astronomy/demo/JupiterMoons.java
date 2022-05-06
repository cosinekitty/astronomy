package io.github.cosinekitty.astronomy.demo;
import io.github.cosinekitty.astronomy.*;

public class JupiterMoons {
    /**
     * Print positions of Jupiter's moons as seen from Earth.
     *
     * Calculates the coordinates of Jupiter and its four major moons
     * (Io, Europa, Ganymede, and Callisto) as seen from the Earth
     * at a given date and time. This function illustrates how to correct
     * for the delay caused by the time it takes for light to reach
     * the Earth from the Jupiter system.
     *
     * Without light travel time correction, observations with a telescope
     * will differ from the calculations by plus or minus 8 minutes over
     * the span of a year, because of the amount of time it takes for light
     * to travel across the diameter of the Earth's orbit.
     *
     * @param time
     * The date and time of the observation.
     */
    public static int run(Time time) {
        // Call geoVector to calculate the geocentric position of Jupiter.
        // geoVector corrects for light travel time.
        // That means it returns a vector to where Jupiter appears to be
        // in the sky, when the light left Jupiter to travel toward the
        // Earth to arrive here at the specified time. This is different from
        // where Jupiter is at that time.

        System.out.printf("Calculations for: %s%n", time);

        Vector jv = Astronomy.geoVector(Body.Jupiter, time, Aberration.Corrected);

        // Calculate the amount of time it took light to reach the Earth from Jupiter.
        // The distance to Jupiter (AU) divided by the speed of light (AU/day) = time in days.
        double lightTravelDays = jv.length() / Astronomy.C_AUDAY;
        System.out.println();
        System.out.println(String.format("It took light %.2f minutes to reach the Earth from Jupiter.", lightTravelDays * Astronomy.MINUTES_PER_DAY));
        System.out.println();

        // The jupiterMoons function calculates positions of Jupiter's moons without
        // correcting for light travel time. Correct for light travel by backdating
        // by the given amount of light travel time.
        Time backdate = time.addDays(-lightTravelDays);

        JupiterMoonsInfo jm = Astronomy.jupiterMoons(backdate);

        // Tricky: the `+` operator for adding `Vector` will throw an exception
        // if the vectors do not have matching times. We work around this
        // by using `withTime` to clone each moon's position vector to have
        // a different time. This is a manual override to work around a safety check.

        printBody("Jupiter",  jv);
        printBody("Io",       jv.plus(jm.getIo().position().withTime(jv.getT())));
        printBody("Europa",   jv.plus(jm.getEuropa().position().withTime(jv.getT())));
        printBody("Ganymede", jv.plus(jm.getGanymede().position().withTime(jv.getT())));
        printBody("Callisto", jv.plus(jm.getCallisto().position().withTime(jv.getT())));
        System.out.println();

        return 0;
    }

    private static void printBody(String name, Vector geovec) {
        Equatorial equ = geovec.toEquatorial();
        System.out.println(String.format("%-8s   RA %10.6f   DEC %10.6f  %10.6f AU", name, equ.getRa(), equ.getDec(), equ.getDist()));
    }
}
