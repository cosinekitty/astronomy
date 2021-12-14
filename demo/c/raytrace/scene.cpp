/*
    scene.cpp

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

    Implements class Scene, which renders a collection of
    SolidObjects and LightSources that illuminate them.
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include "imager.h"
#include "lodepng.h"

namespace Imager
{
    // Empties out the solidObjectList and destroys/frees
    // the SolidObjects that were in it.
    void Scene::ClearSolidObjectList()
    {
        SolidObjectList::iterator iter = solidObjectList.begin();
        SolidObjectList::iterator end  = solidObjectList.end();
        for (; iter != end; ++iter)
        {
            delete *iter;
            *iter = NULL;
        }
        solidObjectList.clear();
    }

    // A limit to how deeply in recursion CalculateLighting may go
    // before it gives up, so as to avoid call stack overflow.
    const int MAX_OPTICAL_RECURSION_DEPTH = 20;

    // A limit to how weak the red, green, or blue intensity of
    // a light ray may be after recursive calls from multiple
    // reflections and/or refractions before giving up.
    // This intensity is deemed too weak to make a significant
    // difference to the image.
    const double MIN_OPTICAL_INTENSITY = 0.001;

    inline bool IsSignificant(const Color& color)
    {
        return
            (color.red   >= MIN_OPTICAL_INTENSITY) ||
            (color.green >= MIN_OPTICAL_INTENSITY) ||
            (color.blue  >= MIN_OPTICAL_INTENSITY);
    }

    Color Scene::TraceRay(
        const Vector& vantage,
        const Vector& direction,
        double refractiveIndex,
        Color rayIntensity,
        int recursionDepth) const
    {
        Intersection intersection;
        const int numClosest = FindClosestIntersection(
            vantage,
            direction,
            intersection);

        switch (numClosest)
        {
        case 0:
            // The ray of light did not hit anything.
            // Therefore we see the background color attenuated
            // by the incoming ray intensity.
            return rayIntensity * backgroundColor;

        case 1:
            // The ray of light struck exactly one closest surface.
            // Determine the lighting using that single intersection.
            return CalculateLighting(
                intersection,
                direction,
                refractiveIndex,
                rayIntensity,
                1 + recursionDepth);

        default:
            // There is an ambiguity: more than one intersection
            // has the same minimum distance.  Caller must catch
            // this exception and have a backup plan for handling
            // this ray of light.
            throw AmbiguousIntersectionException();
        }
    }

    // Determines the color of an intersection,
    // based on illumination it receives via scattering,
    // glossy reflection, and refraction (lensing).
    Color Scene::CalculateLighting(
        const Intersection& intersection,
        const Vector& direction,
        double refractiveIndex,
        Color rayIntensity,
        int recursionDepth) const
    {
        Color colorSum(0.0, 0.0, 0.0);

#if RAYTRACE_DEBUG_POINTS
        if (activeDebugPoint)
        {
            using namespace std;

            Indent(cout, recursionDepth);
            cout << "CalculateLighting[" << recursionDepth << "] {" << endl;

            Indent(cout, 1+recursionDepth);
            cout << intersection << endl;

            Indent(cout, 1+recursionDepth);
            cout << "direction=" << direction << endl;

            Indent(cout, 1+recursionDepth);
            cout.precision(4);
            cout << "refract=" << fixed << refractiveIndex;
            cout << ", intensity=" << rayIntensity << endl;

            Indent(cout, recursionDepth);
            cout << "}" << endl;
        }
#endif

        // Check for recursion stopping conditions.
        // The first is an absolute upper limit on recursion,
        // so as to avoid stack overflow crashes and to
        // limit computation time due to recursive branching.
        if (recursionDepth <= MAX_OPTICAL_RECURSION_DEPTH)
        {
            // The second limit is checking for the ray path
            // having been partially reflected/refracted until
            // it is too weak to matter significantly for
            // determining the associated pixel's color.
            if (IsSignificant(rayIntensity))
            {
                if (intersection.solid == NULL)
                {
                    // If we get here, it means some derived class forgot to
                    // initialize intersection.solid before appending to
                    // the intersection list.
                    throw ImagerException("Undefined solid at intersection.");
                }
                const SolidObject& solid = *intersection.solid;

                // Determine the optical properties at the specified
                // point on whatever solid object the ray intersected with.
                const Optics optics = solid.SurfaceOptics(
                    intersection.point,
                    intersection.context
                );

                // Opacity of a surface point is the fraction 0..1
                // of the light ray available for matte and gloss.
                // The remainder, transparency = 1-opacity, is
                // available for refraction and refractive reflection.
                const double opacity = optics.GetOpacity();
                const double transparency = 1.0 - opacity;
                if (opacity > 0.0)
                {
                    // This object is at least a little bit opaque,
                    // so calculate the part of the color caused by
                    // matte (scattered) reflection.
                    const Color matteColor =
                        opacity *
                        optics.GetMatteColor() *
                        rayIntensity *
                        CalculateMatte(intersection);

                    colorSum += matteColor;

#if RAYTRACE_DEBUG_POINTS
                    if (activeDebugPoint)
                    {
                        using namespace std;

                        Indent(cout, recursionDepth);
                        cout << "matteColor=" << matteColor;
                        cout << ", colorSum=" << colorSum;
                        cout << endl;
                    }
#endif
                }

                double refractiveReflectionFactor = 0.0;
                if (transparency > 0.0)
                {
                    // This object is at least a little bit transparent,
                    // so calculate refraction of the ray passing through
                    // the point. The refraction calculation also tells us
                    // how much reflection was caused by the interface
                    // between the current ray medium and the medium it
                    // is now passing into.  This reflection factor will
                    // be combined with glossy reflection to determine
                    // total reflection below.
                    // Note that only the 'transparent' part of the light
                    // is available for refraction and refractive reflection.

                    colorSum += CalculateRefraction(
                        intersection,
                        direction,
                        refractiveIndex,
                        transparency * rayIntensity,
                        recursionDepth,
                        refractiveReflectionFactor  // output parameter
                    );
                }

                // There are two sources of shiny reflection
                // that need to be considered together:
                // 1. Reflection caused by refraction.
                // 2. The glossy part.

                // The refractive part causes reflection of all
                // colors equally.  Each color component is
                // diminished based on transparency (the part
                // of the ray left available to refraction in
                // the first place).
                Color reflectionColor (1.0, 1.0, 1.0);
                reflectionColor *= transparency * refractiveReflectionFactor;

                // Add in the glossy part of the reflection, which
                // can be different for red, green, and blue.
                // It is diminished to the part of the ray that
                // was not available for refraction.
                reflectionColor += opacity * optics.GetGlossColor();

                // Multiply by the accumulated intensity of the
                // ray as it has traveled around the scene.
                reflectionColor *= rayIntensity;

                if (IsSignificant(reflectionColor))
                {
                    const Color matteColor = CalculateReflection(
                        intersection,
                        direction,
                        refractiveIndex,
                        reflectionColor,
                        recursionDepth);

                    colorSum += matteColor;
                }
            }
        }

#if RAYTRACE_DEBUG_POINTS
        if (activeDebugPoint)
        {
            using namespace std;

            Indent(cout, recursionDepth);
            cout << "CalculateLighting[" << recursionDepth << "] returning ";
            cout << colorSum << endl;
        }
#endif

        return colorSum;
    }

    // Determines the contribution of the illumination of a point
    // based on matte (scatter) reflection based on light incident
    // to a point on the surface of a solid object.
    Color Scene::CalculateMatte(const Intersection& intersection) const
    {
        // Start at the location where the camera ray hit
        // a surface and trace toward all light sources.
        // Add up all the color components to create a
        // composite color value.
        Color colorSum(0.0, 0.0, 0.0);

        // Iterate through all of the light sources.
        LightSourceList::const_iterator iter = lightSourceList.begin();
        LightSourceList::const_iterator end  = lightSourceList.end();
        for (; iter != end; ++iter)
        {
            // Each time through the loop, 'source'
            // will refer to one of the light sources.
            const LightSource& source = *iter;

            // See if we can draw a line from the intersection
            // point toward the light source without hitting any surfaces.
            if (HasClearLineOfSight(intersection.point, source.location))
            {
                // Since there is nothing between this point on the object's
                // surface and the given light source, add this light source's
                // contribution based on the light's color, luminosity,
                // squared distance, and angle with the surface normal.

                // Calculate a direction vector from the intersection point
                // toward the light source point.
                const Vector direction = source.location - intersection.point;

                const double incidence = DotProduct(
                    intersection.surfaceNormal,
                    direction.UnitVector()
                );

                // If the dot product of the surface normal vector and
                // the ray toward the light source is negative, it means
                // light is hitting the surface from the inside of the object,
                // even though we thought we had a clear line of sight.
                // If the dot product is zero, it means the ray grazes
                // the very edge of the object.  Only when the dot product
                // is positive does this light source make the point brighter.
                if (incidence > 0.0)
                {
                    const double intensity =
                        incidence / direction.MagnitudeSquared();

                    colorSum += intensity * source.color;
                }
            }
        }

        return colorSum;
    }


    Color Scene::CalculateReflection(
        const Intersection& intersection,
        const Vector& incidentDir,
        double refractiveIndex,
        Color rayIntensity,
        int recursionDepth) const
    {
        // Find the direction of the reflected ray based on the incident ray
        // direction and the surface normal vector.  The reflected ray has
        // the same angle with the normal vector as the incident ray, but
        // on the opposite side of the cone centered at the normal vector
        // that sweeps out the incident angle.
        const Vector& normal = intersection.surfaceNormal;
        const double perp = 2.0 * DotProduct(incidentDir, normal);
        const Vector reflectDir = incidentDir - (perp * normal);

        // Follow the ray in the new direction from the intersection point.
        return TraceRay(
            intersection.point,
            reflectDir,
            refractiveIndex,
            rayIntensity,
            recursionDepth);
    }

    Color Scene::CalculateRefraction(
        const Intersection& intersection,
        const Vector& direction,
        double sourceRefractiveIndex,
        Color rayIntensity,
        int recursionDepth,
        double& outReflectionFactor) const
    {
        // Convert direction to a unit vector so that
        // relation between angle and dot product is simpler.
        const Vector dirUnit = direction.UnitVector();

        double cos_a1 = DotProduct(dirUnit, intersection.surfaceNormal);
        double sin_a1;
        if (cos_a1 <= -1.0)
        {
            if (cos_a1 < -1.0001)
            {
                throw ImagerException("Dot product too small.");
            }
            // The incident ray points in exactly the opposite
            // direction as the normal vector, so the ray
            // is entering the solid exactly perpendicular
            // to the surface at the intersection point.
            cos_a1 = -1.0;  // clamp to lower limit
            sin_a1 =  0.0;
        }
        else if (cos_a1 >= +1.0)
        {
            if (cos_a1 > +1.0001)
            {
                throw ImagerException("Dot product too large.");
            }
            // The incident ray points in exactly the same
            // direction as the normal vector, so the ray
            // is exiting the solid exactly perpendicular
            // to the surface at the intersection point.
            cos_a1 = +1.0;  // clamp to upper limit
            sin_a1 =  0.0;
        }
        else
        {
            // The ray is entering/exiting the solid at some
            // positive angle with respect to the normal vector.
            // We need to calculate the sine of that angle
            // using the trig identity cos^2 + sin^2 = 1.
            // The angle between any two vectors is always between
            // 0 and PI, so the sine of such an angle is never negative.
            sin_a1 = sqrt(1.0 - cos_a1*cos_a1);
        }

        // The parameter sourceRefractiveIndex passed to this function
        // tells us the refractive index of the medium the light ray
        // was passing through before striking this intersection.
        // We need to figure out what the target refractive index is,
        // i.e., the refractive index of whatever substance the ray
        // is about to pass into.  We determine this by pretending that
        // the ray continues traveling in the same direction a tiny
        // amount beyond the intersection point, then asking which
        // solid object (if any) contains that test point.
        // Ties are broken by insertion order: whichever solid was
        // inserted into the scene first that contains a point is
        // considered the winner.  If a solid is found, its refractive
        // index is used as the target refractive index; otherwise,
        // we use the scene's ambient refraction, which defaults to
        // vacuum (but that can be overridden by a call to
        // Scene::SetAmbientRefraction).

        const double SMALL_SHIFT = 0.001;
        const Vector testPoint = intersection.point + SMALL_SHIFT*dirUnit;
        const SolidObject* container = PrimaryContainer(testPoint);
        const double targetRefractiveIndex =
            (container != NULL) ?
            container->GetRefractiveIndex() :
            ambientRefraction;

        const double ratio = sourceRefractiveIndex / targetRefractiveIndex;

        // Snell's Law: the sine of the refracted ray's angle
        // with the normal is obtained by multiplying the
        // ratio of refractive indices by the sine of the
        // incident ray's angle with the normal.
        const double sin_a2 = ratio * sin_a1;

        if (sin_a2 <= -1.0 || sin_a2 >= +1.0)
        {
            // Since sin_a2 is outside the bounds -1..+1, then
            // there is no such real angle a2, which in turn
            // means that the ray experiences total internal reflection,
            // so that no refracted ray exists.
            outReflectionFactor = 1.0;      // complete reflection
            return Color(0.0, 0.0, 0.0);    // no refraction at all
        }

        // Getting here means there is at least a little bit of
        // refracted light in addition to reflected light.
        // Determine the direction of the refracted light.
        // We solve a quadratic equation to help us calculate
        // the vector direction of the refracted ray.

        double k[2];
        const int numSolutions = Algebra::SolveQuadraticEquation(
            1.0,
            2.0 * cos_a1,
            1.0 - 1.0/(ratio*ratio),
            k);

        // There are generally 2 solutions for k, but only
        // one of them is correct.  The right answer is the
        // value of k that causes the light ray to bend the
        // smallest angle when comparing the direction of the
        // refracted ray to the incident ray.  This is the
        // same as finding the hypothetical refracted ray
        // with the largest positive dot product.
        // In real refraction, the ray is always bent by less
        // than 90 degrees, so all valid dot products are
        // positive numbers.
        double maxAlignment = -0.0001;  // any negative number works as a flag
        Vector refractDir;
        for (int i=0; i < numSolutions; ++i)
        {
            Vector refractAttempt = dirUnit + k[i]*intersection.surfaceNormal;
            double alignment = DotProduct(dirUnit, refractAttempt);
            if (alignment > maxAlignment)
            {
                maxAlignment = alignment;
                refractDir = refractAttempt;
            }
        }

        if (maxAlignment <= 0.0)
        {
            // Getting here means there is something wrong with the math.
            // Either there were no solutions to the quadratic equation,
            // or all solutions caused the refracted ray to bend 90 degrees
            // or more, which is not possible.
            throw ImagerException("Refraction failure.");
        }

        // Determine the cosine of the exit angle.
        double cos_a2 = sqrt(1.0 - sin_a2*sin_a2);
        if (cos_a1 < 0.0)
        {
            // Tricky bit: the polarity of cos_a2 must
            // match that of cos_a1.
            cos_a2 = -cos_a2;
        }

        // Determine what fraction of the light is
        // reflected at the interface.  The caller
        // needs to know this for calculating total
        // reflection, so it is saved in an output parameter.

        // We assume uniform polarization of light,
        // and therefore average the contributions of s-polarized
        // and p-polarized light.
        const double Rs = PolarizedReflection(
            sourceRefractiveIndex,
            targetRefractiveIndex,
            cos_a1,
            cos_a2);

        const double Rp = PolarizedReflection(
            sourceRefractiveIndex,
            targetRefractiveIndex,
            cos_a2,
            cos_a1);

        outReflectionFactor = (Rs + Rp) / 2.0;

        // Whatever fraction of the light is NOT reflected
        // goes into refraction.  The incoming ray intensity
        // is thus diminished by this fraction.
        const Color nextRayIntensity =
            (1.0 - outReflectionFactor) * rayIntensity;

        // Follow the ray in the new direction from the intersection point.
        return TraceRay(
            intersection.point,
            refractDir,
            targetRefractiveIndex,
            nextRayIntensity,
            recursionDepth);
    }

    double Scene::PolarizedReflection(
        double n1,              // source material's index of refraction
        double n2,              // target material's index of refraction
        double cos_a1,          // incident or outgoing ray angle cosine
        double cos_a2) const    // outgoing or incident ray angle cosine
    {
        const double left  = n1 * cos_a1;
        const double right = n2 * cos_a2;
        double numer = left - right;
        double denom = left + right;
        denom *= denom;     // square the denominator
        if (denom < EPSILON)
        {
            // Assume complete reflection.
            return 1.0;
        }
        double reflection = (numer*numer) / denom;
        if (reflection > 1.0)
        {
            // Clamp to actual upper limit.
            return 1.0;
        }
        return reflection;
    }

    int PickClosestIntersection(
        const IntersectionList& list,
        Intersection& intersection)
    {
        // We pick the closest intersection, but we return
        // the number of intersections tied for first place
        // in that contest.  This allows the caller to
        // check for ambiguities in cases where that matters.

        const size_t count = list.size();
        switch (count)
        {
        case 0:
            // No intersection is available.
            // We leave 'intersection' unmodified.
            // The caller must check the return value
            // to know to avoid using 'intersection'.
            return 0;

        case 1:
            // There is exactly one intersection
            // in the given direction, so there is
            // no need to think very hard; just use it!
            intersection = list[0];
            return 1;

        default:
            // There are 2 or more intersections, so we need
            // to find the closest one, and look for ties.
            IntersectionList::const_iterator iter = list.begin();
            IntersectionList::const_iterator end  = list.end();
            IntersectionList::const_iterator closest = iter;
            int tieCount = 1;
            for (++iter; iter != end; ++iter)
            {
                const double diff = iter->distanceSquared - closest->distanceSquared;
                if (fabs(diff) < EPSILON)
                {
                    // Within tolerance of the closest so far,
                    // so consider this a tie.
                    ++tieCount;
                }
                else if (diff < 0.0)
                {
                    // This new intersection is definitely closer
                    // to the vantage point.
                    tieCount = 1;
                    closest = iter;
                }
            }
            intersection = *closest;

            // The caller may need to know if there was an ambiguity,
            // so report back the total number of closest intersections.
            return tieCount;
        }
    }

    // Searches for an intersections with any solid in the scene from the
    // vantage point in the given direction.  If none are found, the
    // function returns 0 and the 'intersection' parameter is left
    // unchanged.  Otherwise, returns the positive number of
    // intersections that lie at minimal distance from the vantage point
    // in that direction.  Usually this number will be 1 (a unique
    // intersection is closer than all the others) but it can be greater
    // if multiple intersections are equally close (e.g. the ray hitting
    // exactly at the corner of a cube could cause this function to
    // return 3).  If this function returns a value greater than zero,
    // it means the 'intersection' parameter has been filled in with the
    // closest intersection (or one of the equally closest intersections).
    int Scene::FindClosestIntersection(
        const Vector& vantage,
        const Vector& direction,
        Intersection& intersection) const
    {
        // Build a list of all intersections from all objects.
        cachedIntersectionList.clear();     // empty any previous contents
        SolidObjectList::const_iterator iter = solidObjectList.begin();
        SolidObjectList::const_iterator end  = solidObjectList.end();
        for (; iter != end; ++iter)
        {
            const SolidObject& solid = *(*iter);
            solid.AppendAllIntersections(
                vantage,
                direction,
                cachedIntersectionList);
        }
        return PickClosestIntersection(cachedIntersectionList, intersection);
    }


    // Returns true if nothing blocks a line drawn between point1 and point2.
    bool Scene::HasClearLineOfSight(
        const Vector& point1,
        const Vector& point2) const
    {
        // Subtract point2 from point1 to obtain the direction
        // from point1 to point2, along with the square of
        // the distance between the two points.
        const Vector dir = point2 - point1;
        const double gapDistanceSquared = dir.MagnitudeSquared();

        // Iterate through all the solid objects in this scene.
        SolidObjectList::const_iterator iter = solidObjectList.begin();
        SolidObjectList::const_iterator end  = solidObjectList.end();
        for (; iter != end; ++iter)
        {
            // If any object blocks the line of sight,
            // we can return false immediately.
            const SolidObject& solid = *(*iter);

            // Find the closest intersection from point1
            // in the direction toward point2.
            Intersection closest;
            if (0 != solid.FindClosestIntersection(point1, dir, closest))
            {
                // We found the closest intersection, but it is only
                // a blocker if it is closer to point1 than point2 is.
                // If the closest intersection is farther away than
                // point2, there is nothing on this object blocking
                // the line of sight.

                if (closest.distanceSquared < gapDistanceSquared)
                {
                    // We found a surface that is definitely blocking
                    // the line of sight.  No need to keep looking!
                    return false;
                }
            }
        }

        // We would not find any solid object that blocks the line of sight.
        return true;
    }

    // Generate an image of the scene and write it to the
    // specified output PNG file.
    // outPngFileName is the name of the PNG file to write the image to.
    // pixelsWide, pixelsHigh are the pixel dimensions of the output file.
    // The zoom is a positive number that controls the magnification of
    // the image: smaller values magnify the image more (zoom in),
    // and larger values shrink all the scenery to fit more objects
    // into the image (zoom out).
    // Adjust antiAliasFactor to increase the amount over oversampling
    // to make smoother (less jagged) looking images.
    // Generally, antiAliasFactor should be between 1 (fastest, but jagged)
    // and 4 (16 times slower, but very smooth looking).
    void Scene::SaveImage(
        const char *outPngFileName,
        size_t pixelsWide,
        size_t pixelsHigh,
        double zoom,
        size_t antiAliasFactor) const
    {
        // Oversample the image using the anti-aliasing factor.
        const size_t largePixelsWide = antiAliasFactor * pixelsWide;
        const size_t largePixelsHigh = antiAliasFactor * pixelsHigh;
        const size_t smallerDim =
            ((pixelsWide < pixelsHigh) ? pixelsWide : pixelsHigh);

        const double largeZoom  = antiAliasFactor * zoom * smallerDim;
        ImageBuffer buffer(largePixelsWide, largePixelsHigh, backgroundColor);

        // The camera is located at the origin.
        Vector camera(0.0, 0.0, 0.0);

        // The camera faces in the -z direction.
        // This allows the +x direction to be to the right,
        // and the +y direction to be upward.
        Vector direction(0.0, 0.0, -1.0);

        const Color fullIntensity(1.0, 1.0, 1.0);

        // We keep a list of (i,j) screen coordinates for pixels
        // we are not able to trace definitive rays for.
        // Later we will come back and fix these pixels.
        PixelList ambiguousPixelList;

        for (size_t i=0; i < largePixelsWide; ++i)
        {
            direction.x = (i - largePixelsWide/2.0) / largeZoom;
            for (size_t j=0; j < largePixelsHigh; ++j)
            {
                direction.y = (largePixelsHigh/2.0 - j) / largeZoom;

#if RAYTRACE_DEBUG_POINTS
                {
                    using namespace std;

                    // Assume no active debug point unless we find one below.
                    activeDebugPoint = NULL;

                    DebugPointList::const_iterator iter = debugPointList.begin();
                    DebugPointList::const_iterator end  = debugPointList.end();
                    for(; iter != end; ++iter)
                    {
                        if ((iter->iPixel == i) && (iter->jPixel == j))
                        {
                            cout << endl;
                            cout << "Hit breakpoint at (";
                            cout << i << ", " << j <<")" << endl;
                            activeDebugPoint = &(*iter);
                            break;
                        }
                    }
                }
#endif

                PixelData& pixel = buffer.Pixel(i,j);
                try
                {
                    Vector aim = (aimer != nullptr) ? aimer->Aim(direction) : direction;

                    // Trace a ray from the camera toward the given direction
                    // to figure out what color to assign to this pixel.
                    pixel.color = TraceRay(camera, aim, ambientRefraction, fullIntensity, 0);
                }
                catch (AmbiguousIntersectionException)
                {
                    // Getting here means that somewhere in the recursive
                    // code for tracing rays, there were multiple
                    // intersections that had minimum distance from a
                    // vantage point.  This can be really bad,
                    // for example causing a ray of light to reflect
                    // inward into a solid.

                    // Mark the pixel as ambiguous, so that any other
                    // ambiguous pixels nearby know not to use it.
                    pixel.isAmbiguous = true;

                    // Keep a list of all ambiguous pixel coordinates
                    // so that we can rapidly enumerate through them
                    // in the disambiguation pass.
                    ambiguousPixelList.push_back(PixelCoordinates(i, j));
                }
            }
        }

#if RAYTRACE_DEBUG_POINTS
        // Leave no chance of a dangling pointer into debug points.
        activeDebugPoint = NULL;
#endif

        // Go back and "heal" ambiguous pixels as best we can.
        PixelList::const_iterator iter = ambiguousPixelList.begin();
        PixelList::const_iterator end  = ambiguousPixelList.end();
        for (; iter != end; ++iter)
        {
            const PixelCoordinates& p = *iter;
            ResolveAmbiguousPixel(buffer, p.i, p.j);
        }

        // We want to scale the arbitrary range of
        // color component values to the range 0..255
        // allowed by PNG format.  We therefore find
        // the maximum red, green, or blue value anywhere
        // in the image.
        const double max = buffer.MaxColorValue();

        // Downsample the image buffer to an integer array of RGBA
        // values that LodePNG understands.
        const unsigned char OPAQUE_ALPHA_VALUE = 255;
        const unsigned BYTES_PER_PIXEL = 4;

        // The number of bytes in buffer to be passed to LodePNG.
        const unsigned RGBA_BUFFER_SIZE =
            pixelsWide * pixelsHigh * BYTES_PER_PIXEL;

        std::vector<unsigned char> rgbaBuffer(RGBA_BUFFER_SIZE);
        unsigned rgbaIndex = 0;
        const double patchSize = antiAliasFactor * antiAliasFactor;
        for (size_t j=0; j < pixelsHigh; ++j)
        {
            for (size_t i=0; i < pixelsWide; ++i)
            {
                Color sum(0.0, 0.0, 0.0);
                for (size_t di=0; di < antiAliasFactor; ++di)
                {
                    for (size_t dj=0; dj < antiAliasFactor; ++dj)
                    {
                        sum += buffer.Pixel(
                            antiAliasFactor*i + di,
                            antiAliasFactor*j + dj).color;
                    }
                }
                sum /= patchSize;

                // Convert to integer red, green, blue, alpha values,
                // all of which must be in the range 0..255.
                rgbaBuffer[rgbaIndex++] = ConvertPixelValue(sum.red,   max);
                rgbaBuffer[rgbaIndex++] = ConvertPixelValue(sum.green, max);
                rgbaBuffer[rgbaIndex++] = ConvertPixelValue(sum.blue,  max);
                rgbaBuffer[rgbaIndex++] = OPAQUE_ALPHA_VALUE;
            }
        }

        // Write the PNG file
        const unsigned error = lodepng::encode(
            outPngFileName,
            rgbaBuffer,
            pixelsWide,
            pixelsHigh);

        // If there was an encoding error, throw an exception.
        if (error != 0)
        {
            std::string message = "PNG encoder error: ";
            message += lodepng_error_text(error);
            throw ImagerException(message.c_str());
        }
    }

    // The following function searches through all solid objects
    // for the first solid (if any) that contains the given point.
    // In the case of ties, the solid that was inserted into the
    // scene first wins.  This arbitrary convention allows the
    // composer of a scene to decide which of multiple overlapping
    // objects should control the index of refraction for any
    // overlapping volumes of space.
    const SolidObject* Scene::PrimaryContainer(const Vector& point) const
    {
        SolidObjectList::const_iterator iter = solidObjectList.begin();
        SolidObjectList::const_iterator end  = solidObjectList.end();
        for (; iter != end; ++iter)
        {
            const SolidObject* solid = *iter;
            if (solid->Contains(point))
            {
                return solid;
            }
        }

        return NULL;
    }

    void Scene::ResolveAmbiguousPixel(
        ImageBuffer& buffer,
        size_t i,
        size_t j) const
    {
        // This function is called whenever SaveImage could not
        // figure out what color to assign to a pixel, because
        // multiple intersections were found that minimize the
        // distance to the vantage point.

        // Avoid going out of bounds with pixel coordinates.
        const size_t iMin = (i > 0) ? (i - 1) : i;
        const size_t iMax = (i < buffer.GetPixelsWide()-1) ? (i + 1) : i;
        const size_t jMin = (j > 0) ? (j - 1) : j;
        const size_t jMax = (j < buffer.GetPixelsHigh()-1) ? (j + 1) : j;

        // Look for surrounding unambiguous pixels.
        // Average their color values together.
        Color colorSum(0.0, 0.0, 0.0);
        int numFound = 0;
        for (size_t si = iMin; si <= iMax; ++si)
        {
            for (size_t sj = jMin; sj <= jMax; ++sj)
            {
                const PixelData& pixel = buffer.Pixel(si, sj);
                if (!pixel.isAmbiguous)
                {
                    ++numFound;
                    colorSum += pixel.color;
                }
            }
        }

        if (numFound > 0)   // avoid division by zero
        {
            colorSum /= numFound;
        }

        // "Airbrush" out the imperfection.
        // This is not perfect, but it looks a lot better
        // than leaving the pixel some arbitrary color,
        // and better than picking the wrong intersection
        // and following it into a crazy direction.
        buffer.Pixel(i, j).color = colorSum;
    }
}
