/*
    seasons.c  -  Don Cross  -  2019-06-16

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program demonstrates how to calculate the
    equinoxes and solstices for a year.
*/

#include <stdio.h>
#include "astro_demo_common.h"
char printBuffer[80];

void Printd(const char *name, astro_time_t time)
{
    snprintf(printBuffer, sizeof printBuffer,"%-17s : ", name);
    Serial.print(printBuffer);
    PrintTime(time);
    Serial.println();
}

void setup() {
  Serial.begin(9600);
}

void loop() 
{
    int year;
    astro_seasons_t seasons;

    const char *argvs[] = {"0","1988"};
    int args = 2;
    if (args != 2)
    {
        Serial.println("USAGE: seasons year");

    }

    if (1 != sscanf(argvs[1], "%d", &year))
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Invalid year '%s'.\n", argvs[1]);
        Serial.println(printBuffer);
        
    }

    seasons = Astronomy_Seasons(year);
    if (seasons.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_Seasons() returned %d\n", seasons.status);
        Serial.println(printBuffer);
    }

    Printd("March equinox", seasons.mar_equinox);
    Printd("June solstice", seasons.jun_solstice);
    Printd("September equinox", seasons.sep_equinox);
    Printd("December solstice", seasons.dec_solstice);

    delay(2000);
}
