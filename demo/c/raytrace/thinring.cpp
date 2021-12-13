/*
    thinring.cpp

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
    Implements class ThinRing.
    A thin ring is a disc-shaped surface with zero thickness 
    that has a smaller disc-shaped hole at a common center.
*/

#include "planet.h"

namespace Imager
{
    void ThinRing::ObjectSpace_AppendAllIntersections(
        const Vector& vantage, 
        const Vector& direction, 
        IntersectionList& intersectionList) const
    {
        if (fabs(direction.z) > EPSILON)
        {
            const double u = -vantage.z / direction.z;
            if (u > EPSILON)
            {
                const double x = u*direction.x + vantage.x;
                const double y = u*direction.y + vantage.y;
                const double m = x*x + y*y;
                if ((m <= r2*r2 + EPSILON) && (r1*r1 <= m + EPSILON))
                {
                    Intersection intersection;
                    intersection.point = Vector(x, y, 0.0);
                    intersection.distanceSquared = (u * direction).MagnitudeSquared();

                    // We "cheat" a little bit in calculating the normal vector by knowing too much about the caller.
                    // This is necessary because this is not really a normal solid object, but an infinitesimally thin surface.
                    // Therefore, we provide a different normal vector depending on the supplied vantage, 
                    // such that the point of view of the observer determines which side of the surface is seen.
                    // (Doing otherwise would cause the surface to appear completely black in some cases.)
                    intersection.surfaceNormal = Vector(0.0, 0.0, (vantage.z >= 0.0) ? +1.0 : -1.0);

                    intersection.solid = this;
                    intersectionList.push_back(intersection);
                }
            }
        }
    }
}
