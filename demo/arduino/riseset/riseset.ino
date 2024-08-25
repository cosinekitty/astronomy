/*
    riseset.c  -  by Don Cross - 2019-06-14

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates sunrise, sunset, moonrise, and moonset
    times for an observer at a given latitude and longitude.
*/


#include "astro_demo_common.h"
char printBuffer[80];

int PrintEvent(const char *name, astro_search_result_t evt)
{
    snprintf(printBuffer, sizeof printBuffer,"%-8s : ", name);
    Serial.print(printBuffer);
    if (evt.status == ASTRO_SUCCESS)
    {
        PrintTime(evt.time);
        Serial.println();
        return 0;
    }
    else
    {
        snprintf(printBuffer, sizeof printBuffer,"ERROR %d\n", evt.status);
        Serial.print(printBuffer);
        return 1;
    }
}

void setup() {
  Serial.begin(9600);
}

void loop() {
    astro_observer_t observer;
    astro_time_t time;
    astro_search_result_t sunrise, sunset, moonrise, moonset;

    const char *argvs[] = {"22.2","11.8","1988-02-02T12:55:33Z"};
    int args = 3;

    ParseArgs(args, argvs, &observer, &time);
        

    snprintf(printBuffer, sizeof printBuffer,"search   : ");
    Serial.print(printBuffer);
    PrintTime(time);
    Serial.println();

    sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 300.0);
    sunset   = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET,  time, 300.0);
    moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 300.0);
    moonset  = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET,  time, 300.0);

    PrintEvent("sunrise", sunrise);
        

    PrintEvent("sunset", sunset);
        

    PrintEvent("moonrise", moonrise);
        

    PrintEvent("moonset", moonset);
    
    delay(2000);

    
}
