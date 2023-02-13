package io.github.cosinekitty.astronomy.demo;

import java.time.Instant;
import java.time.format.DateTimeParseException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import io.github.cosinekitty.astronomy.*;

public class Main {
    private static final String usageText = String.join(System.getProperty("line.separator"),
        "",
        "In all demos that include [yyyy-mm-ddThh:mm:ssZ]",
        "on the command line, that pattern indicates an optional date/time.",
        "If the date/time is specified, it is used for the calculation.",
        "If absent, the computer's current date and time is used.",
        "",
        "Command line arguments:",
        "",
        "    constellation [yyyy-mm-ddThh:mm:ssZ]",
        "        Finds what constellation the Moon is in at a given time.",
        "        Then it finds the Moon's constellation changes over the",
        "        subsequent 30 days.",
        "",
        "    jupiter_moons [yyyy-mm-ddThh:mm:ssZ]",
        "        Calculates the coordinates of Jupiter and its four major moons",
        "        (Io, Europa, Ganymede, and Callisto) as seen from the Earth",
        "        at a given date and time. This demo illustrates how to correct",
        "        for the delay caused by the time it takes for light to reach",
        "        the Earth from the Jupiter system.",
        "",
        "    lunar_eclipse [yyyy-mm-ddThh:mm:ssZ]",
        "        Searches for the first 10 lunar eclipses (partial or total)",
        "        that occur after the specified time. Penumbral lunar eclipses",
        "        are ignored, as these are difficult to observe in practice.",
        "",
        "    moonphase [yyyy-mm-ddThh:mm:ssZ]",
        "        Calculates the Moon's ecliptic phase and illumination percentage",
        "        for a given date and time. Also finds the dates and times of",
        "        the subsequent 10 quarter phases of the Moon.",
        "",
        "    positions latitude longitude [yyyy-mm-ddThh:mm:ssZ]",
        "        Displays the equatorial and horizontal coordinates of",
        "        the Sun, Moon, and planets, as seen from a given",
        "        geographic location. Uses the date and time specified on",
        "        the command line, if present. Otherwise, uses the computer's",
        "        current date and time.",
        "",
        "    riseset latitude longitude [yyyy-mm-ddThh:mm:ssZ]",
        "        Displays the times of the following 6 events:",
        "        sunrise, sun culmination, sunset,",
        "        moonrise, moon culmination, moonset.",
        "        Culmination is when a body reaches the highest",
        "        point in an observer's sky as it crosses the meridian.",
        "        The specified time is the starting point of the search.",
        "        The displayed events are those that occur first after that time.",
        "        The events are displayed in chronological order.",
        "",
        "    seasons year",
        "        Given an integer year number, displays the solstices and equinoxes for that year.",
        "        The year must be in the range 0000..9999.",
        "",
        "    solar_time latitude longitude [yyyy-mm-ddThh:mm:ssZ]",
        "        Displays the true solar time for the observer at the",
        "        given geographic coordinates.",
        ""
    );

    private static int printUsage() {
        System.out.println(usageText);
        return 1;
    }

    private static class DemoException extends Exception {
        public DemoException(String message) {
            super(message);
        }
    }

    public static void main(String[] args) {
        System.exit(runDemos(args));
    }

    private static int runDemos(String[] args) {
        if (args.length == 0) return printUsage();
        String verb = args[0];
        Optional<Demo> foundDemo = demoList.stream()
                .filter(demo -> demo.name.equals(verb))
                .findFirst();
        if (foundDemo.isEmpty()) {
            System.out.printf("ERROR: Unknown command '%s'.%n", verb);
            return 1;
        }
        Demo demo = foundDemo.get();
        if (args.length < demo.minArgs || args.length > demo.maxArgs) {
            System.out.println("ERROR: Incorrect number of command-line arguments.");
            return 1;
        }
        try {
            return demo.runner.run(args);
        } catch (DateTimeParseException e) {
            System.out.println("ERROR: Invalid date/time format on command line.");
            return 1;
        } catch (DemoException e) {
            System.out.printf("ERROR: %s%n", e.getMessage());
            return 1;
        }
    }

    private static Time parseTime(String[] args, int index) {
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

    private interface DemoRunner {
        int run(String[] args) throws DemoException;
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

    private static final List<Demo> demoList = Arrays.asList(
        new Demo("constellation", 1, 2, args ->
            Constellation.run(
                parseTime(args, 1)
            )
        ),
        new Demo("jupiter_moons", 1, 2, args ->
            JupiterMoons.run(
                parseTime(args, 1)
            )
        ),
        new Demo("lunar_eclipse", 1, 2, args ->
            LunarEclipse.run(
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
        new Demo("riseset", 3, 4, args ->
            RiseSetCulm.run(
                parseObserver(args, 1),
                parseTime(args, 3)
            )
        ),
        new Demo("seasons", 2, 2, args ->
            Seasons.run(
                parseYear(args[1])
            )
        ),
        new Demo("solar_time", 3, 4, args ->
            SolarTime.run(
                parseObserver(args, 1),
                parseTime(args, 3)
            )
        )
    );
}
