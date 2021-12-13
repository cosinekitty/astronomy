/*
    planet.h

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

    Planets of the Solar System, rendered to scale: 1 unit per Earth mean diameter.
*/

#ifndef __DDC_IMAGER_PLANET_H
#define __DDC_IMAGER_PLANET_H

#include "imager.h"

namespace Imager
{
    const double MEAN_EARTH_RADIUS_KM = 6371.0;

    class Planet: public Spheroid
    {
    public:
        Planet(
            Color _color,
            double _equatorialRadiusInKm,
            double _polarRadiusInKm)
                : Spheroid(
                    _equatorialRadiusInKm / MEAN_EARTH_RADIUS_KM,
                    _equatorialRadiusInKm / MEAN_EARTH_RADIUS_KM,
                    _polarRadiusInKm / MEAN_EARTH_RADIUS_KM)
        {
            SetFullMatte(_color);
        }
    };

    class Saturn: public SetUnion
    {
    public:
        Saturn()
            : SetUnion(
                Vector(),
                new Planet(Color(1.0, 1.0, 0.6, 0.335), 60268.0, 54364.0),
                CreateRingSystem())
        {
        }

    private:
        static SolidObject* CreateRingSystem();
    };
}

#endif // __DDC_IMAGER_PLANET_H