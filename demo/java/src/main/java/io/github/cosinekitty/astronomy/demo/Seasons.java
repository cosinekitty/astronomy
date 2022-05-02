package io.github.cosinekitty.astronomy.demo;
import io.github.cosinekitty.astronomy.*;

public class Seasons {
    public static int run(int year) {
        SeasonsInfo seasons = Astronomy.seasons(year);
        System.out.printf("March equinox     : %s%n", seasons.getMarchEquinox());
        System.out.printf("June solstice     : %s%n", seasons.getJuneSolstice());
        System.out.printf("September equinox : %s%n", seasons.getSeptemberEquinox());
        System.out.printf("December solstice : %s%n", seasons.getDecemberSolstice());
        return 0;
    }
}
