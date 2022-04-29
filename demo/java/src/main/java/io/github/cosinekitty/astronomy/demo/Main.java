package io.github.cosinekitty.astronomy.demo;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
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
        Pattern pattern = Pattern.compile("^(\\d{4})-(\\d{2})-(\\d{2})(T(\\d{2}):(\\d{2})(:(\\d{2}(\\.\\d+)?))?Z)$");
        Matcher matcher = pattern.matcher(args[index]);
        if (matcher.find()) {
            int year = Integer.parseInt(matcher.group(1));
            int month = Integer.parseInt(matcher.group(2));
            int day = Integer.parseInt(matcher.group(3));
            int hour = 0;
            int minute = 0;
            double second = 0.0;
            if (!matcher.group(4).isEmpty()) {
                hour = Integer.parseInt(matcher.group(5));
                minute = Integer.parseInt(matcher.group(6));
                if (!matcher.group(7).isEmpty()) {
                    second = Double.parseDouble((matcher.group(8)));
                }
            }
            return new Time(year, month, day, hour, minute, second);
        }
        System.out.print("FATAL: Invalid date/time syntax: ");
        System.out.println(args[index]);
        return null;
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
                    Date now = new Date();
                    Time time = Time.fromMillisecondsSince1970(now.getTime());
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
