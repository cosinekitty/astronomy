package io.github.cosinekitty.astronomy.demo;

import java.util.stream.Stream;

import io.github.cosinekitty.astronomy.*;

public class LunarEclipse {
    /**
     * Searches for the next 10 partial/total lunar eclipses.
     *
     * @param startTime
     * The date and time after which to start searching for lunar eclipses.
     */
    public static int run(Time startTime) {
        LunarEclipseInfo e = Astronomy.searchLunarEclipse(startTime);
        Stream.iterate(e, x -> Astronomy.nextLunarEclipse(x.getPeak()))
                .filter(x -> x.getKind() != EclipseKind.Penumbral)
                .limit(10)
                .forEach(LunarEclipse::printEclipse);
        return 0;
    }

    private static void printEclipse(LunarEclipseInfo e) {
        // Calculate beginning/ending of different phases
        // of an eclipse by subtracting/adding the peak time
        // with the number of minutes indicated by the "semi-duration"
        // fields sdPartial and sdTotal.

        Time peak = e.getPeak();            // the central time of the eclipse, when the Moon is darkest
        double total = e.getSdTotal();      // the semiduration of totality in minutes, or 0 if not total
        double partial = e.getSdPartial();  // the semiduration of partiality in minutes
        Time p1 = peak.addDays(-e.getSdPartial() / Astronomy.MINUTES_PER_DAY);
        System.out.printf("%s  Partial eclipse begins.%n", p1);
        if (total > 0.0) {
            Time t1 = peak.addDays(-total / Astronomy.MINUTES_PER_DAY);
            System.out.printf("%s  Total eclipse begins.%n", t1);
        }
        System.out.printf("%s  Peak of %s eclipse.%n", peak, e.getKind().toString().toLowerCase());
        if (total > 0.0) {
            Time t2 = peak.addDays(+total / Astronomy.MINUTES_PER_DAY);
            System.out.printf("%s  Total eclipse ends.%n", t2);
        }
        Time p2 = peak.addDays(+partial / Astronomy.MINUTES_PER_DAY);
        System.out.printf("%s  Partial eclipse ends.%n", p2);
        System.out.println();   // extra blank line to separate each lunar eclipse
    }
}
