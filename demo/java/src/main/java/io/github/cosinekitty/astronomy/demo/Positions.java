package io.github.cosinekitty.astronomy.demo;
import io.github.cosinekitty.astronomy.*;

public class Positions {
    private static final Body[] bodyList = new Body[] {
        Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
        Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
    };

    /**
     * Print a table of solar system body positions in equatorial and horizontal coordinates.
     *
     * @param observer
     * The geographic location for which to calculate positions.
     *
     * @param time
     * The date and time for which to calculate positions.
     */
    public static int run(Observer observer, Time time) {
        System.out.printf("UTC date = %s%n", time);
        System.out.println();
        System.out.println("BODY           RA      DEC       AZ      ALT");
        for (Body body : bodyList) {
            Equatorial equ_2000 = Astronomy.equator(body, time, observer, EquatorEpoch.J2000, Aberration.Corrected);
            Equatorial equ_ofdate = Astronomy.equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
            Topocentric hor = Astronomy.horizon(time, observer, equ_ofdate.getRa(), equ_ofdate.getDec(), Refraction.Normal);
            System.out.printf("%-8s %8.2f %8.2f %8.2f %8.2f%n", body, equ_2000.getRa(), equ_2000.getDec(), hor.getAzimuth(), hor.getAltitude());
        }
        return 0;
    }
}
