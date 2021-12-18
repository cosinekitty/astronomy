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

                colorSum = optics.GetMatteColor() * rayIntensity * CalculateMatte(intersection);
            }
        }

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
