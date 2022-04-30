package io.github.cosinekitty.astronomy.demo;
import io.github.cosinekitty.astronomy.*;

public class MoonPhase {
    public static int run(Time time) {
        // Calculate the Moon's ecliptic phase angle,
        // which ranges from 0 to 360 degrees.
        //   0 degrees = new mooon,
        //  90 degrees = first quarter,
        // 180 degrees = full moon,
        // 270 degrees = third quarter
        double phase = Astronomy.moonPhase(time);
        System.out.printf("%s : Moon's ecliptic phase angle = %1.3f degrees.%n", time, phase);

        // Calculate the fraction of the Moon's disc that appears
        // illuminated, as seen from the Earth.
        IlluminationInfo illum = Astronomy.illumination(Body.Moon, time);
        System.out.printf("%s : Moon's illuminated fraction = %1.2f%%.%n", time, 100.0 * illum.getPhaseFraction());

        // Predict when the next 10 lunar quarter phases will happen.
        System.out.println();
        System.out.println("The next 10 lunar quarters are:");
        MoonQuarterInfo mq = Astronomy.searchMoonQuarter(time);
        for (int i = 0; i < 10; ++i) {
            if (i > 0) {
                mq = Astronomy.nextMoonQuarter(mq);
            }
            System.out.printf("%s : %s%n", mq.getTime(), quarterName(mq.getQuarter()));
        }

        return 0;
    }

    private static String quarterName(int quarter) {
        switch (quarter) {
            case 0: return "New Moon";
            case 1: return "First Quarter";
            case 2: return "Full Moon";
            case 3: return "Third Quarter";
            default: return "INVALID QUARTER";
        }
    }
}
