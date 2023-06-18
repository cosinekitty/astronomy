/*
    altazsearch.cpp  -  Don Cross  -  2023-06-17

    https://github.com/cosinekitty/astronomy/discussions/308

    Problem: given a range of altitudes and azimuths that
    form a "window" on the sky, search for when the Moon
    enters/exits that window starting from a given search time.
*/

#include <cmath>
#include <cstdio>
#include <string>
#include <stdexcept>
#include "astronomy.h"


struct Event
{
    astro_time_t    time;
    double          azimuth = NAN;
    double          altitude = NAN;

    Event()
    {
        time.tt = time.ut = time.eps = time.psi = time.st = NAN;
    }

    void Print() const
    {
        char text[TIME_TEXT_BYTES];
        Astronomy_FormatTime(time, TIME_FORMAT_SECOND, text, sizeof(text));
        printf("%s az=%0.2lf alt=%0.2lf", text, azimuth, altitude);
    }
};


struct Solution
{
    bool valid = false;
    Event start;
    Event finish;

    void Print() const
    {
        if (valid)
        {
            printf("Start: ");
            start.Print();
            printf("; Finish: ");
            finish.Print();
            printf(".\n");
        }
        else
        {
            printf("No solution.\n");
        }
    }
};


void Verify(astro_status_t status, const char *message)
{
    if (status != ASTRO_SUCCESS)
        throw std::logic_error(std::string(message) + ": error " + std::to_string(static_cast<int>(status)));
}


class SearchProblem
{
private:
    astro_body_t body;
    astro_observer_t observer;
    double az1;
    double az2;
    double alt1;
    double alt2;
    astro_vector_t center;

    astro_horizon_t Position(astro_time_t time) const
    {
        // Get topocentric equatorial coordinates of body, using the Earth's equator of date.
        astro_equatorial_t equ = Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION);
        Verify(equ.status, "Equator");

        // Convert to observer's horizontal coordinates, correcting for atmospheric refraction.
        return Astronomy_Horizon(&time, observer, equ.ra, equ.dec, REFRACTION_NORMAL);
    }

    double AngularDistance(astro_time_t time) const
    {
        astro_horizon_t hor = Position(time);

        // Translate angular horizontal coordinates to a vector.
        // Do NOT remove the refraction from the vector (REFRACTION_NONE).
        astro_spherical_t sphere;
        sphere.status = ASTRO_SUCCESS;
        sphere.dist = 1.0;
        sphere.lat = hor.altitude;
        sphere.lon = hor.azimuth;
        astro_vector_t vec = Astronomy_VectorFromHorizon(sphere, time, REFRACTION_NONE);
        Verify(vec.status, "VectorFromHorizon");

        // Calculate the angle in degrees between the body and the center of the window.
        astro_angle_result_t result = Astronomy_AngleBetween(center, vec);
        Verify(result.status, "AngleBetween");

        return result.angle;
    }

    static astro_func_result_t DistanceSlopeCallback(void *context, astro_time_t time)
    {
        const SearchProblem& p = *static_cast<const SearchProblem *>(context);
        astro_func_result_t result;
        result.value = p.DistanceSlope(time);
        result.status = ASTRO_SUCCESS;
        return result;
    }

    double DistanceSlope(astro_time_t time) const
    {
        const double dt = 0.1 / 86400.0;
        astro_time_t t1 = Astronomy_AddDays(time, -dt);
        astro_time_t t2 = Astronomy_AddDays(time, +dt);
        double a1 = AngularDistance(t1);
        double a2 = AngularDistance(t2);
        return (a2 - a1) / (2 * dt);
    }

    bool IsInsideWindow(astro_time_t time) const
    {
        astro_horizon_t hor = Position(time);
        bool insideAzimuthLimits;
        if (az1 <= az2)
            insideAzimuthLimits = (az1 <= hor.azimuth && hor.azimuth <= az2);
        else
            insideAzimuthLimits = (az1 <= hor.azimuth || hor.azimuth <= az2);

        return insideAzimuthLimits && alt1 <= hor.altitude && hor.altitude <= alt2;
    }

    Solution FindBracket(astro_time_t closestTime) const
    {
        Solution solution;

        // If the closestTime is inside the window, we can find a bracket.
        // Otherwise, there is no bracket here.
        if (IsInsideWindow(closestTime))
        {
            const double dt = 10.0 / (24.0 * 60.0);     // 10 minutes, converted to days

            // Look backwards until we find a time before entering the window.
            // Do a binary search to find the moment when we enter the window.
            solution.start = FindTransition(closestTime, -dt);

            // Look forward until we find a time after leaving the window.
            // Do a binary search to find the moment we leave the window.
            solution.finish = FindTransition(closestTime, +dt);

            solution.valid = true;
        }

        return solution;
    }

    Event FindTransition(astro_time_t closestTime, double dt) const
    {
        // Find a bracket [t1, t2] that straddles being inside/outside the window.
        astro_time_t t1 = closestTime;
        astro_time_t t2 = Astronomy_AddDays(closestTime, dt);
        while (IsInsideWindow(t2))
        {
            t1 = t2;
            t2 = Astronomy_AddDays(t2, dt);
        }

        // Do a binary search to find the moment of transition, within tolerance.
        const double tolerance = 0.1 / (3600.0 * 24.0);     // one tenth of a second, expressed in days
        astro_time_t tm = Astronomy_TimeFromDays((t1.ut + t2.ut) / 2);
        while (fabs(t2.ut - t1.ut) > tolerance)
        {
            if (IsInsideWindow(tm))
                t1 = tm;
            else
                t2 = tm;
            tm = Astronomy_TimeFromDays((t1.ut + t2.ut) / 2);
        }

        astro_horizon_t hor = Position(tm);

        Event event;
        event.time = tm;
        event.altitude = hor.altitude;
        event.azimuth = hor.azimuth;
        return event;
    }

public:
    SearchProblem(astro_body_t _body, astro_observer_t _observer, double _az1, double _az2, double _alt1, double _alt2)
        : body(_body)
        , observer(_observer)
        , az1(_az1)
        , az2(_az2)
        , alt1(_alt1)
        , alt2(_alt2)
    {
        astro_time_t dummyTime = Astronomy_TimeFromDays(0.0);
        astro_spherical_t sphere;
        sphere.status = ASTRO_SUCCESS;
        sphere.dist = 1.0;
        sphere.lat = (_alt1 + _alt2) / 2;
        sphere.lon = (_az1 + _az2) / 2;
        center = Astronomy_VectorFromHorizon(sphere, dummyTime, REFRACTION_NONE);
    }

    Solution FindNext(astro_time_t startTime, double limitDays)
    {
        astro_time_t stopTime = Astronomy_AddDays(startTime, limitDays);
        const double stepDays = 1.0 / 24.0;     // one hour
        astro_time_t t1 = startTime;
        double m1 = DistanceSlope(t1);

        while (t1.ut < stopTime.ut)
        {
            astro_time_t t2 = Astronomy_AddDays(t1, stepDays);
            double m2 = DistanceSlope(t2);
            if (m1 <= 0.0 && m2 >= 0.0)
            {
                astro_search_result_t result = Astronomy_Search(DistanceSlopeCallback, this, t1, t2, 0.1);
                if (result.status == ASTRO_SUCCESS)
                {
                    // We found a time bracket [t1, t2] where the body passes closest to
                    // the center of the target window. Now search nearby for when the
                    // body enters and exits the window.
                    Solution solution = FindBracket(result.time);
                    if (solution.valid)
                        return solution;
                }
            }
            t1 = t2;
            m1 = m2;
        }

        // We could not find a solution.
        // Return a default-constructed Solution, which indicates failure.
        return Solution();
    }
};


int main(int argc, const char *argv[])
{
    // 2023-06-17T12:00:00Z  Moon    AZ = 74.03   ALT = 27.44
    // 2023-06-17T13:00:00Z  Moon    AZ = 79.22   ALT = 39.74
    astro_observer_t observer = Astronomy_MakeObserver(30.0, -80.0, 0.0);
    SearchProblem problem(BODY_MOON, observer, 74.0, 78.0, 25.0, 40.0);
    astro_time_t startTime = Astronomy_MakeTime(2023, 6, 17, 0, 0, 0.0);
    Solution solution = problem.FindNext(startTime, 2.0);
    solution.Print();
    return 0;
}

