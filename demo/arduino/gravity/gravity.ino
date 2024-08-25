


#include <stdio.h>
#include <math.h>

#include "astronomy.h"

void setup() {
  // put your setup code here, to run once:

Serial.begin(9600);



}

void loop() {
  // put your main code here, to run repeatedly:

    const double MAX_HEIGHT_METERS = 100000.0;
    double latitude, height, gravity;
    latitude = 23.3;
    height = 1000;
    gravity = Astronomy_ObserverGravity(latitude, height);
    // printf("latitude = %8.4lf,  height = %6.0lf,  gravity = %8.6lf\n", latitude, height, gravity);
    Serial.printf("latitude = %8.4lf,  height = %6.0lf,  gravity = %8.6lf\n", latitude, height, gravity);
    // Serial.println(gravity,6);
    delay(500);
}
//--->ASTRONOMY_ENGINE_WHOLE_SECOND