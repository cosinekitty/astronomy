/*
    spheroid.cpp

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
    Implements class Spheroid.
    A spheroid is like a sphere, only it may have three different 
    diameters in the x-, y-, and z-directions.
*/

#include "algebra.h"
#include "imager.h"

namespace Imager
{
    void Spheroid::ObjectSpace_AppendAllIntersections(
        const Vector& vantage, 
        const Vector& direction, 
        IntersectionList& intersectionList) const
    {
        double u[2];
        const int numSolutions = Algebra::SolveQuadraticEquation(
            b2*c2*direction.x*direction.x + a2*c2*direction.y*direction.y + a2*b2*direction.z*direction.z,
            2.0*(b2*c2*vantage.x*direction.x + a2*c2*vantage.y*direction.y + a2*b2*vantage.z*direction.z),
            b2*c2*vantage.x*vantage.x + a2*c2*vantage.y*vantage.y + a2*b2*vantage.z*vantage.z - a2*b2*c2,
            u
        );

        for (int i=0; i < numSolutions; ++i)
        {
            if (u[i] > EPSILON)
            {
                Intersection intersection;
                Vector displacement = u[i] * direction;
                intersection.distanceSquared = displacement.MagnitudeSquared();
                intersection.point = vantage + displacement;

                // The surface normal vector was calculated by expressing the spheroid as a 
                // function z(x,y) = sqrt(1 - (x/a)^2 - (y/b)^2),
                // taking partial derivatives dz/dx = (c*c*x)/(a*a*z), dz/dy = (c*c*y)/(b*b*z), 
                // and using these to calculate the vectors <1, 0, dz/dx> and <0, 1, dy,dz>.
                // The normalized cross product of these two vectors yields the surface normal vector.
                const double x = intersection.point.x;
                const double y = intersection.point.y;
                const double z = intersection.point.z;

                // But we need to handle special cases when z is very close to 0.
                if (fabs(z) <= EPSILON)
                {
                    if (fabs(x) <= EPSILON)
                    {
                        // The equation devolves to (y^2)/(b^2) = 1, or y = +/- b.
                        intersection.surfaceNormal = Vector(0.0, ((y > 0.0) ? 1.0 : -1.0), 0.0);
                    }
                    else
                    {
                        // The equation devolves to an ellipse on the xy plane : 
                        // (x^2)/(a^2) + (y^2)/(b^2) = 1.
                        intersection.surfaceNormal = Vector(-1.0, -(a2*y)/(b2*x), 0.0).UnitVector();
                    }
                }
                else
                {
                    intersection.surfaceNormal = Vector((c2*x)/(a2*z), (c2*y)/(b2*z), 1.0).UnitVector();
                }

                // Handle special cases with polarity: the polarity of the components of
                // the surface normal vector must match that of the <x,y,z> intersection point,
                // because the surface normal vector always points (roughly) away from the vantage, just
                // like any point on the surface of the spheroid does.

                if (x * intersection.surfaceNormal.x < 0.0)     // negative product means opposite polarities
                {
                    intersection.surfaceNormal.x *= -1.0;
                }

                if (y * intersection.surfaceNormal.y < 0.0)     // negative product means opposite polarities
                {
                    intersection.surfaceNormal.y *= -1.0;
                }

                if (z * intersection.surfaceNormal.z < 0.0)     // negative product means opposite polarities
                {
                    intersection.surfaceNormal.z *= -1.0;
                }

                intersection.solid = this;

                intersectionList.push_back(intersection);
            }
        }
    }
}
