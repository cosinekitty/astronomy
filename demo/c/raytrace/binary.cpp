/*
    binary.cpp

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

    Implementation of class SolidObject_BinaryOperator, which is the base class
    for operators applied to two other SolidObjects.
*/

#include "imager.h"

namespace Imager
{
    // All rotations and translations are applied to the two nested solids in tandem.

    SolidObject& SolidObject_BinaryOperator::Translate(double dx, double dy, double dz)
    {
        SolidObject::Translate(dx, dy, dz);     // chain to base class
        Left().Translate(dx, dy, dz);           // translate left solid
        Right().Translate(dx, dy, dz);          // translate right solid
        return *this;
    }

    SolidObject& SolidObject_BinaryOperator::RotateX(double angleInDegrees)    // rotates counterclockwise around center looking into axis parallel to x-axis.
    {
        const double angleInRadians = RadiansFromDegrees(angleInDegrees);
        const double a = cos(angleInRadians);
        const double b = sin(angleInRadians);

        NestedRotateX(Left(),  angleInDegrees, a, b);
        NestedRotateX(Right(), angleInDegrees, a, b);

        return *this;
    }

    SolidObject& SolidObject_BinaryOperator::RotateY(double angleInDegrees)    // rotates counterclockwise around center looking into axis parallel to y-axis.
    {
        const double angleInRadians = RadiansFromDegrees(angleInDegrees);
        const double a = cos(angleInRadians);
        const double b = sin(angleInRadians);

        NestedRotateY(Left(),  angleInDegrees, a, b);
        NestedRotateY(Right(), angleInDegrees, a, b);

        return *this;
    }

    SolidObject& SolidObject_BinaryOperator::RotateZ(double angleInDegrees)    // rotates counterclockwise around center looking into axis parallel to z-axis.
    {
        const double angleInRadians = RadiansFromDegrees(angleInDegrees);
        const double a = cos(angleInRadians);
        const double b = sin(angleInRadians);

        NestedRotateZ(Left(),  angleInDegrees, a, b);
        NestedRotateZ(Right(), angleInDegrees, a, b);

        return *this;
    }

    void SolidObject_BinaryOperator::NestedRotateX(SolidObject &nested, double angleInDegrees, double a, double b)
    {
        // Rotate the nested object about its own center.
        nested.RotateX(angleInDegrees);

        // Revolve the center of the nested object around the common center of this binary operator.
        const Vector& c = Center();
        const Vector& nc = nested.Center();
        const double dy = nc.y - c.y;
        const double dz = nc.z - c.z;
        nested.Move (nc.x, c.y + (a*dy - b*dz), c.z + (a*dz + b*dy));
    }

    void SolidObject_BinaryOperator::NestedRotateY(SolidObject &nested, double angleInDegrees, double a, double b)
    {
        // Rotate the nested object about its own center.
        nested.RotateY(angleInDegrees);

        // Revolve the center of the nested object around the common center of this binary operator.
        const Vector& c = Center();
        const Vector& nc = nested.Center();
        const double dx = nc.x - c.x;
        const double dz = nc.z - c.z;
        nested.Move (c.x + (a*dx + b*dz), nc.y, c.z + (a*dz - b*dx));
    }

    void SolidObject_BinaryOperator::NestedRotateZ(SolidObject &nested, double angleInDegrees, double a, double b)
    {
        // Rotate the nested object about its own center.
        nested.RotateZ(angleInDegrees);

        // Revolve the center of the nested object around the common center of this binary operator.
        const Vector& c = Center();
        const Vector& nc = nested.Center();
        const double dx = nc.x - c.x;
        const double dy = nc.y - c.y;
        nested.Move (c.x + (a*dx - b*dy), c.y + (a*dy + b*dx), nc.z);
    }
}
