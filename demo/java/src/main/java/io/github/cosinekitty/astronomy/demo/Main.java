package io.github.cosinekitty.astronomy.demo;

import java.util.Date;
import io.github.cosinekitty.astronomy.AstroTime;

public class Main {
    private static final String usageText = String.join(System.getProperty("line.separator"),
        "Command line arguments:",
        "",
        "    now",
        "    Display current date and time using AstroTime.",
        ""
    );

    public static void main(String[] args) {
        int rc = 1;
        if (args.length == 0) {
            System.out.println(usageText);
        } else {
            switch (args[0]) {
                case "now" -> {
                    System.out.println(new AstroTime(new Date()));
                    rc = 0;
                }
                default -> System.out.println("ERROR: Unknown command line argument");
            }
        }
        System.exit(rc);
    }
}
