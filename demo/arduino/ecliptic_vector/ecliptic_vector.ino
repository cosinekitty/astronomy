
#include "stdio.h"
#include "astro_demo_common.h"


void setup() 
{
  Serial.begin(9600);

  // setSyncProvider(syncProvider);//sets internal clock
  // if(timeStatus() != timeSet) 
  //   Serial.println("Unable to sync with the unixTimeStamp");
  // else
  //   Serial.println("unixTimeStamp has set the system time");   
        

}

void loop() 
{
      static const astro_body_t body[] = {
        BODY_SUN, BODY_MERCURY, BODY_VENUS, BODY_EARTH, BODY_MOON, BODY_MARS,
        BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO
    };
    int i;
int num_bodies;
    astro_time_t time;
    astro_rotation_t rot;
    astro_time_t time2;
  num_bodies = sizeof(body) / sizeof(body[0]);;

const char *argv[] = {"123","2023-10-07T00:38:57+03:30"};

            if (0 == ParseTime(argv[1], &time))
            {

              rot = Astronomy_Rotation_EQJ_ECL();
              // Serial.println("BODY               X           Y           Z\n");
              for (i=0; i < num_bodies; ++i)
                {
                    astro_vector_t eqj = Astronomy_HelioVector(body[i], time);
                    astro_vector_t ecl = Astronomy_RotateVector(rot, eqj);
                    const char *name = Astronomy_BodyName(body[i]);
                    if (ecl.status != ASTRO_SUCCESS)
                    {
                         Serial.println( "ERROR %d calculating vector for %s.\n");
                         Serial.println((int)ecl.status);
                         Serial.println(name);
                    }

                    Serial.printf("%-8s %11.6lf %11.6lf %11.6lf\n", name, ecl.x, ecl.y, ecl.z);
                }
            }

  delay(1000);
}


