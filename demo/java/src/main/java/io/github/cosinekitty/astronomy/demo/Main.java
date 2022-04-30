package io.github.cosinekitty.astronomy.demo;

import java.time.Instant;
import java.time.format.DateTimeParseException;
import java.util.Date;

import io.github.cosinekitty.astronomy.*;

public class Main {
    private static String UsageText = String.join(System.getProperty("line.separator"),
        "Command line arguments:",
        "",
        "    moonphase [yyyy-mm-ddThh:mm:ssZ]",
        "       Calculates the Moon's ecliptic phase and illumination percentage",
        "       for a given date and time, or for the computer's current date and",
        "       time if none is given on the command line.",
        "       Also finds the dates and times of the subsequent 10 quarter phases.",
        "",
        "    now",
        "       Display current date and time.",
        ""
    );

    private static Time parseTime(String args[], int index) {
        if (index >= args.length) {
            Date now = new Date();
            return Time.fromMillisecondsSince1970(now.getTime());
        }
        try {
            Instant instant = Instant.parse(args[index]);
            return Time.fromMillisecondsSince1970(instant.toEpochMilli());
        } catch (DateTimeParseException e) {
            System.out.print("FATAL: Invalid date/time syntax: ");
            System.out.println(args[index]);
            return null;
        }
    }

    public static void main(String[] args) {
        int rc = 1;
        if (args.length == 0) {
            System.out.println(UsageText);
        } else {
            switch (args[0]) {
                case "moonphase":
                    if (args.length <= 2) {
                        Time time = parseTime(args, 1);
                        if (time != null) {
                            rc = MoonPhase.run(time);
                        }
                    } else {
                        System.out.println(UsageText);
                    }
                    break;

                case "now":
                    Time time = Time.fromMillisecondsSince1970(System.currentTimeMillis());
                    System.out.println(time);
                    rc = 0;
                    break;

                default:
                    System.out.println("ERROR: Unknown command line argument");
                    break;
            }
        }
        System.exit(rc);
    }
}
