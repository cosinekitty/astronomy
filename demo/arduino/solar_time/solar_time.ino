/*
    solar_time.c  -  by Don Cross - 2019-06-11

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Given an observer's geographic latitude and longitude,
    and an optional date and time, this program displays
    true solar time for that observer and time.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"
char printBuffer[80];

void setup() {
  Serial.begin(9600);
}

void loop() {
    int error;
    astro_observer_t observer;
    astro_time_t time;
    astro_func_result_t ha;
    double solarTimeHours;
    int hour, minute, second, milli;

    const char *argvs[] = {"22.2","11.8","1988-02-02T12:55:33Z"};
    int args = 3;

    error = ParseArgs(args, argvs, &observer, &time);
    if (error)
        Serial.println("error..... Args");

    ha = Astronomy_HourAngle(BODY_SUN, &time, observer);
    if (ha.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer,"ERROR %d in Astronomy_HourAngle().\n", ha.status);
        Serial.println(printBuffer);
    }

    solarTimeHours = fmod(ha.value + 12.0, 24.0);

    milli = (int)round(solarTimeHours * 3600000.0);
    second = milli / 1000;
    milli %= 1000;
    minute = second / 60;
    second %= 60;
    hour = minute / 60;
    minute %= 60;
    hour %= 24;

    snprintf(printBuffer, sizeof printBuffer,"True solar time = %0.4lf hours (%02d:%02d:%02d.%03d)\n", solarTimeHours, hour, minute, second, milli);
    Serial.println(printBuffer);
    delay(2000);
}
