package io.github.cosinekitty.astronomy.demo;

import java.time.Instant;
import java.time.format.DateTimeParseException;

import io.github.cosinekitty.astronomy.*;

public class Main {
    private static String usageText = String.join(System.getProperty("line.separator"),
        "Command line arguments:",
        "",
        "    moonphase [yyyy-mm-ddThh:mm:ssZ]",
        "       Calculates the Moon's ecliptic phase and illumination percentage",
        "       for a given date and time, or for the computer's current date and",
        "       time if none is given on the command line.",
        "       Also finds the dates and times of the subsequent 10 quarter phases.",
        "",
        "    now [yyyy-mm-ddThh:mm:ssZ]",
        "       Display current date and time, or the time supplied on the command line.",
        "",
        "    seasons year",
        "       Given an integer year number, display the solstices and equinoxes for that year.",
        ""
    );

    private static Time parseTime(String args[], int index) {
        long millis =
            (index >= 0 && index < args.length)
            ? Instant.parse(args[index]).toEpochMilli()
            : System.currentTimeMillis();

        return Time.fromMillisecondsSince1970(millis);
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
                System.out.println("FATAL: Invalid date/time syntax.");
                rc = 1;
            }
        }
        System.exit(rc);
    }

    private static interface DemoRunner {
        public int run(String[] args);
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

    private static int printTime(String[] args) {
        System.out.println(parseTime(args, 1));
        return 0;
    }

    private static Demo[] demoList = new Demo[] {
        new Demo("moonphase", 1, 2, args -> MoonPhase.run(parseTime(args, 1))),
        new Demo("now", 1, 2, args -> printTime(args)),
        new Demo("seasons", 2, 2, args -> Seasons.run(args[1]))
    };
}
