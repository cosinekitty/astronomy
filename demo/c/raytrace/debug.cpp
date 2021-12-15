/*
    debug.cpp

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

    Implements output operators and other stuff 
    useful for debugging the ray tracer.
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include "imager.h"

namespace Imager
{
    void Indent(std::ostream& output, int depth)
    {
        const int numSpaces = 4 * depth;
        for (int i=0; i < numSpaces; ++i)
        {
            output << ' ';
        }
    }

    std::ostream& operator<< (std::ostream &output, const Color& color)
    {
        using namespace std;

        output.precision(3);
        output << scientific;
        output << "(R " << color.red << ", G " << color.green << ", B " << color.blue << ")";
        return output;
    }

    std::ostream& operator<< (std::ostream& output, const Intersection& intersection)
    {
        using namespace std;

        output << "{ ";
        if (intersection.solid)
        {
            output << intersection.solid->GetTag();
        }
        if (intersection.tag)
        {
            output << ":" << intersection.tag;
        }
        output.precision(3);
        output << fixed;
        output << " d=" << intersection.distanceSquared;
        output << ", p=" << intersection.point;
        output << ", n=" << intersection.surfaceNormal;
        output << " }";
        return output;
    }

    std::ostream& operator<< (std::ostream& output, const Vector& vector)
    {
        using namespace std;

        output.precision(3);
        output << fixed;
        output << "(" << vector.x << ", " << vector.y << ", " << vector.z << ")";
        return output;
    }
}
