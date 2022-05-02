package io.github.cosinekitty.astronomy.demo;

import java.time.Instant;
import java.time.format.DateTimeParseException;

import io.github.cosinekitty.astronomy.*;

public class Main {
    private static String usageText = String.join(System.getProperty("line.separator"),
        "Command line arguments:",
        "",
        "    jupiter_moons [yyyy-mm-ddThh:mm:ssZ]",
        "        Calculates the coordinates of Jupiter and its four major moons",
        "        (Io, Europa, Ganymede, and Callisto) as seen from the Earth",
        "        at a given date and time. This demo illustrates how to correct",
        "        for the delay caused by the time it takes for light to reach",
        "        the Earth from the Jupiter system.",
        "",
        "    moonphase [yyyy-mm-ddThh:mm:ssZ]",
        "        Calculates the Moon's ecliptic phase and illumination percentage",
        "        for a given date and time, or for the computer's current date and",
        "        time if none is given on the command line.",
        "        Also finds the dates and times of the subsequent 10 quarter phases.",
        "",
        "    positions latitude longitude [yyyy-mm-ddThh:mm:ssZ]",
        "        Displays the equatorial and horizontal coordinates of",
        "        the Sun, Moon, and planets, as seen from a given",
        "        geographic location. Uses the date and time specified on",
        "        the command line, if present. Otherwise, uses the computer's",
        "        current date and time.",
        "",
        "    seasons year",
        "        Given an integer year number, displays the solstices and equinoxes for that year.",
        "        The year must be in the range 0000..9999.",
        ""
    );

    private static class DemoException extends Exception {
        public DemoException(String message) {
            super(message);
        }
    }

    public static void main(String[] args) {
        int rc = 1;
        if (args.length == 0) {
            System.out.println(usageText);
        } else {
            try {
                String verb = args[0];
                boolean found = false;
                for (int i = 0; i < demoList.length; ++i) {
                    Demo demo = demoList[i];
                    if (demo.name.equals(verb)) {
                        found = true;
                        if (args.length >= demo.minArgs && args.length <= demo.maxArgs) {
                            rc = demo.runner.run(args);
                        } else {
                            System.out.println(usageText);
                        }
                        break;
                    }
                }
                if (!found) {
                    System.out.printf("ERROR: Unknown command '%s'.%n", verb);
                }
            } catch (DateTimeParseException e) {
                System.out.println("ERROR: Invalid date/time format on command line.");
                rc = 1;
            } catch (DemoException e) {
                System.out.printf("ERROR: %s%n", e.getMessage());
                rc = 1;
            }
        }
        System.exit(rc);
    }

    private static Time parseTime(String args[], int index) {
        long millis =
            (index >= 0 && index < args.length)
            ? Instant.parse(args[index]).toEpochMilli()
            : System.currentTimeMillis();

        return Time.fromMillisecondsSince1970(millis);
    }

    private static double parseNumber(String name, String text, double minValue, double maxValue) throws DemoException {
        try {
            double value = Double.parseDouble(text);
            if (!Double.isFinite(value) || value < minValue || value > maxValue) {
                throw new DemoException(String.format("Value is out of range for %s.", name));
            }
            return value;
        } catch (NumberFormatException e) {
            throw new DemoException(String.format("Invalid numeric format '%s' for %s.%n", text, name));
        }
    }

    private static int parseYear(String text) throws DemoException {
        try {
            int year = Integer.parseInt(text);
            if (year < 0 || year > 9999) {
                throw new DemoException("Year must be in the range 0000..9999.");
            }
            return year;
        } catch (NumberFormatException e) {
            throw new DemoException(String.format("Invalid numeric format '%s' for year.", text));
        }
    }

    private static Observer parseObserver(String[] args, int index) throws DemoException {
        double latitude = parseNumber("latitude", args[index], -90.0, +90.0);
        double longitude = parseNumber("longitude", args[index+1], -180.0, +180.0);
        return new Observer(latitude, longitude, 0.0);
    }

    private static interface DemoRunner {
        public int run(String[] args) throws DemoException;
    }

    private static class Demo {
        public final String name;
        public final int minArgs;
        public final int maxArgs;
        public final DemoRunner runner;

        public Demo(String name, int minArgs, int maxArgs, DemoRunner runner) {
            this.name = name;
            this.minArgs = minArgs;
            this.maxArgs = maxArgs;
            this.runner = runner;
        }
    }

    private static Demo[] demoList = new Demo[] {
        new Demo("jupiter_moons", 1, 2, args ->
            JupiterMoons.run(
                parseTime(args, 1)
            )
        ),
        new Demo("moonphase", 1, 2, args ->
            MoonPhase.run(
                parseTime(args, 1)
            )
        ),
        new Demo("positions", 3, 4, args ->
            Positions.run(
                parseObserver(args, 1),
                parseTime(args, 3)
            )
        ),
        new Demo("seasons", 2, 2, args ->
            Seasons.run(
                parseYear(args[1])
            )
        )
    };
}
