/*
    solid.cpp

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
    Contains common code for base class SolidObject.
*/

#include "imager.h"

namespace Imager
{
    bool SolidObject::Contains(const Vector& point) const
    {
        // FIXFIXFIX:  This function does not handle the "corner case":
        // multiple intersections found at the same point but for
        // different facets of the solid.

        if (isFullyEnclosed)
        {
            // This method assumes that the solid's surfaces fully
            // enclose a volume of space without any gaps or cracks.
            // Pick an arbitrary direction in space and count the number
            // of times we enter and exit this solid.
            const Vector direction(0.0, 0.0, 1.0);

            enclosureList.clear();
            AppendAllIntersections(point, direction, enclosureList);

            int enterCount = 0;     // number of times we enter the solid
            int exitCount  = 0;     // number of times we exit the solid

            IntersectionList::const_iterator iter = enclosureList.begin();
            IntersectionList::const_iterator end  = enclosureList.end();
            for (; iter != end; ++iter)
            {
                const Intersection& intersection = *iter;

                // Calculate the dot product of the direction with 
                // the surface normal.  
                const double dotprod = DotProduct(
                    direction, 
                    intersection.surfaceNormal);

                // If it is positive, we are exiting the solid.
                // If it is negative, we are entering the solid.  
                if (dotprod > EPSILON)
                {
                    ++exitCount;
                }
                else if (dotprod < -EPSILON)
                {
                    ++enterCount;
                }
                else
                {
                    // If the dot product is too close to zero, 
                    // something odd is going on because we 
                    // should not have found an intersection
                    // with a plane in the first place.
                    throw ImagerException("Ambiguous transition.");
                }
            }

            // If the original point is within this solid,
            // we have exited the object one more time than we entered.
            // Otherwise, we have exited and entered the same number of times.
            switch (exitCount - enterCount)
            {
            case 0:
                return false;   // point is outside the solid

            case 1:
                return true;    // point is inside the solid

            default:
                // This can happen only if the solid's surfaces
                // do not properly enclose a volume of space without
                // gaps.  Either the surfaces need to be corrected
                // so as to perfectly seal the interior volume
                // or this instance should be constructed with 
                // isFullyEnclosed initialized to false.
                throw ImagerException("Cannot determine containment.");
            }
        }
        else
        {
            // Whoever constructed this object has indicated that
            // it should not be considered to contain any points.
            return false;
        }
    }
}
