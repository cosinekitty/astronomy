/*
    planet.cpp

    Copyright (C) 2013 by Don Cross  -  http://cosinekitty.com/raytrace

    This software is provided 'as-is', without any express or implied
    warranty. In no event will the author be held liable for any damages
    arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
       distribution.

    -------------------------------------------------------------------------

    Contains simple 3D models of solar system planets.
*/

#include "planet.h"

namespace Imager
{
    /*static*/ SolidObject* Saturn::CreateRingSystem()
    {
        struct RingData
        {
            double  innerRadiusKm;
            double  outerRadiusKm;
            double  red;
            double  green;
            double  blue;
        };

        static const RingData ringData[] =
        {
            {  92154.0,  117733.0,  196.0,  194.0,  180.0 },
            { 122405.0,  133501.0,  157.0,  160.0,  158.0 },
            { 134085.0,  136888.0,  136.0,  140.0,  142.0 }
        };

        SolidObject* ringSystem = NULL;

        static const size_t NUM_RINGS = sizeof(ringData) / sizeof(ringData[0]);
        for (size_t i=0; i < NUM_RINGS; ++i)
        {
            const Color color(ringData[i].red / 255.0, ringData[i].green / 255.0, ringData[i].blue / 255.0);

            ThinRing* ringSolid = new ThinRing(ringData[i].innerRadiusKm / MEAN_EARTH_RADIUS_KM, ringData[i].outerRadiusKm / MEAN_EARTH_RADIUS_KM);
            ringSolid->SetFullMatte(color);

            if (ringSystem)
            {
                ringSystem = new SetUnion(Vector(), ringSolid, ringSystem);
            }
            else
            {
                ringSystem = ringSolid;
            }
        }

        return ringSystem;
    }
}
