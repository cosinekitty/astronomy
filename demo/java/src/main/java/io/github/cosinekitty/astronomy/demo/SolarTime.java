package io.github.cosinekitty.astronomy.demo;
import io.github.cosinekitty.astronomy.*;

public class SolarTime {
    /**
     * Display true solar time for a given geographic location at a given time.
     *
     * @param observer
     * The geographic location for which to calculate true solar time.
     *
     * @param time
     * The date and time for which to calculate true solar time.
     */
    public static int run(Observer observer, Time time) {
        double ha = Astronomy.hourAngle(Body.Sun, time, observer);
        double solarTimeHours = (ha + 12.0) % 24.0;
        int milli = (int) Math.round(solarTimeHours * 3.6e+6);
        int second = milli / 1000;
        milli %= 1000;
        int minute = second / 60;
        second %= 60;
        int hour = minute / 60;
        minute %= 60;
        hour %= 24;

        System.out.printf("True solar time = %7.4f hours (%02d:%02d:%02d.%03d)%n", solarTimeHours, hour, minute, second, milli);
        return 0;
    }
}
