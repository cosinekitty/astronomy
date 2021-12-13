/*
    optics.cpp

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

    Member functions for class Optics, which describes the optical properties
    of a point on the surface of a solid object.
*/

#include <cmath>
#include "imager.h"

namespace Imager
{
    void Optics::ValidateReflectionColor(const Color& color) const
    {
        // A color is valid for reflection if all its
        // components are in the range 0..1.
        // Otherwise, it is possible for multiple reflections
        // to keep "amplifying" light beyond any limit.
        if (color.red < 0.0 || color.red > 1.0)
        {
            throw ImagerException("Invalid red color component.");
        }

        if (color.green < 0.0 || color.green > 1.0)
        {
            throw ImagerException("Invalid green color component.");
        }

        if (color.blue < 0.0 || color.blue > 1.0)
        {
            throw ImagerException("Invalid blue color component.");
        }
    }

    void Optics::SetMatteColor(const Color& _matteColor)
    {
        ValidateReflectionColor(_matteColor);
        matteColor = _matteColor;
    }

    void Optics::SetGlossColor(const Color& _glossColor)
    {
        ValidateReflectionColor(_glossColor);
        glossColor = _glossColor;
    }

    void Optics::SetMatteGlossBalance(
        double glossFactor,
        const Color& rawMatteColor,
        const Color& rawGlossColor)
    {
        // Make sure the raw colors have realistic values.
        // Otherwise very high component values can defeat
        // the purpose of trying to balance reflected light
        // levels to realistic ranges.
        ValidateReflectionColor(rawMatteColor);
        ValidateReflectionColor(rawGlossColor);

        // glossFactor must be in the range 0..1.
        if (glossFactor < 0.0 || glossFactor > 1.0)
        {
            throw ImagerException("Gloss factor must be in the range 0..1");
        }

        // Scale the glossy and matte parts of reflected light
        // so that an object does not reflect more light than hits it.
        SetMatteColor((1.0 - glossFactor) * rawMatteColor);
        SetGlossColor(glossFactor * rawGlossColor);
    }

    void Optics::SetOpacity(double _opacity)
    {
        if (_opacity < 0.0 || _opacity > 1.0)
        {
            throw ImagerException("Invalid opacity.");
        }
        opacity = _opacity;
    }
}