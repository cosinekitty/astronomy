/**
    @preserve

    Astronomy library for JavaScript (browser and Node.js).
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2022 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

/**
 * @fileoverview Astronomy calculation library for browser scripting and Node.js.
 * @author Don Cross <cosinekitty@gmail.com>
 * @license MIT
 */
'use strict';

export type FlexibleDateTime = Date | number | AstroTime;

/**
 * @brief The speed of light in AU/day.
 */
export const C_AUDAY = 173.1446326846693;

/**
 * @brief The number of kilometers per astronomical unit.
 */
export const KM_PER_AU = 1.4959787069098932e+8;

/**
 * @brief The factor to convert degrees to radians = pi/180.
 */
export const DEG2RAD = 0.017453292519943296;

/**
 * @brief The factor to convert sidereal hours to radians = pi/12.
 */
export const HOUR2RAD = 0.2617993877991494365;

 /**
  * @brief The factor to convert radians to degrees = 180/pi.
  */
export const RAD2DEG = 57.295779513082321;

 /**
  * @brief The factor to convert radians to sidereal hours = 12/pi.
  */
export const RAD2HOUR = 3.819718634205488;


// Jupiter radius data are nominal values obtained from:
// https://www.iau.org/static/resolutions/IAU2015_English.pdf
// https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

/**
 * @brief The equatorial radius of Jupiter, expressed in kilometers.
 */
export const JUPITER_EQUATORIAL_RADIUS_KM = 71492.0;

 /**
  * @brief The polar radius of Jupiter, expressed in kilometers.
  */
export const JUPITER_POLAR_RADIUS_KM = 66854.0;

 /**
  * @brief The volumetric mean radius of Jupiter, expressed in kilometers.
  */
export const JUPITER_MEAN_RADIUS_KM = 69911.0;


// The radii of Jupiter's four major moons are obtained from:
// https://ssd.jpl.nasa.gov/?sat_phys_par

/**
 * @brief The mean radius of Jupiter's moon Io, expressed in kilometers.
 */
export const IO_RADIUS_KM = 1821.6;

/**
 * @brief The mean radius of Jupiter's moon Europa, expressed in kilometers.
 */
export const EUROPA_RADIUS_KM = 1560.8;

/**
 * @brief The mean radius of Jupiter's moon Ganymede, expressed in kilometers.
 */
export const GANYMEDE_RADIUS_KM = 2631.2;

/**
 * @brief The mean radius of Jupiter's moon Callisto, expressed in kilometers.
 */
export const CALLISTO_RADIUS_KM = 2410.3;


const DAYS_PER_TROPICAL_YEAR = 365.24217;
const J2000 = new Date('2000-01-01T12:00:00Z');
const PI2 = 2 * Math.PI;
const ARC = 3600 * (180 / Math.PI);     // arcseconds per radian
const ASEC2RAD = 4.848136811095359935899141e-6;
const ASEC180 = 180 * 60 * 60;              // arcseconds per 180 degrees (or pi radians)
const ASEC360 = 2 * ASEC180;                // arcseconds per 360 degrees (or 2*pi radians)
const ANGVEL = 7.2921150e-5;
const AU_PER_PARSEC = ASEC180 / Math.PI;    // exact definition of how many AU = one parsec
const SUN_MAG_1AU = -0.17 - 5*Math.log10(AU_PER_PARSEC);    // formula from JPL Horizons
const MEAN_SYNODIC_MONTH = 29.530588;       // average number of days for Moon to return to the same phase
const SECONDS_PER_DAY = 24 * 3600;
const MILLIS_PER_DAY = SECONDS_PER_DAY * 1000;
const SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592;

const SUN_RADIUS_KM = 695700.0;
const SUN_RADIUS_AU  = SUN_RADIUS_KM / KM_PER_AU;

const EARTH_FLATTENING = 0.996647180302104;
const EARTH_FLATTENING_SQUARED = EARTH_FLATTENING * EARTH_FLATTENING;
const EARTH_EQUATORIAL_RADIUS_KM = 6378.1366;
const EARTH_EQUATORIAL_RADIUS_AU = EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU;
const EARTH_POLAR_RADIUS_KM = EARTH_EQUATORIAL_RADIUS_KM * EARTH_FLATTENING;
const EARTH_MEAN_RADIUS_KM = 6371.0;    /* mean radius of the Earth's geoid, without atmosphere */
const EARTH_ATMOSPHERE_KM = 88.0;       /* effective atmosphere thickness for lunar eclipses */
const EARTH_ECLIPSE_RADIUS_KM = EARTH_MEAN_RADIUS_KM + EARTH_ATMOSPHERE_KM;

const MOON_EQUATORIAL_RADIUS_KM = 1738.1;
const MOON_MEAN_RADIUS_KM       = 1737.4;
const MOON_POLAR_RADIUS_KM      = 1736.0;
const MOON_EQUATORIAL_RADIUS_AU = (MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU);

const REFRACTION_NEAR_HORIZON = 34 / 60;        // degrees of refractive "lift" seen for objects near horizon
const EARTH_MOON_MASS_RATIO = 81.30056;

/*
    Masses of the Sun and outer planets, used for:
    (1) Calculating the Solar System Barycenter
    (2) Integrating the movement of Pluto

    https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf

    Page 10 in the above document describes the constants used in the DE405 ephemeris.
    The following are G*M values (gravity constant * mass) in [au^3 / day^2].
    This side-steps issues of not knowing the exact values of G and masses M[i];
    the products GM[i] are known extremely accurately.
*/
const SUN_GM     = 0.2959122082855911e-03;
const MERCURY_GM = 0.4912547451450812e-10;
const VENUS_GM   = 0.7243452486162703e-09;
const EARTH_GM   = 0.8887692390113509e-09;
const MARS_GM    = 0.9549535105779258e-10;
const JUPITER_GM = 0.2825345909524226e-06;
const SATURN_GM  = 0.8459715185680659e-07;
const URANUS_GM  = 0.1292024916781969e-07;
const NEPTUNE_GM = 0.1524358900784276e-07;
const PLUTO_GM   = 0.2188699765425970e-11;

const MOON_GM = EARTH_GM / EARTH_MOON_MASS_RATIO;

/**
 * @brief Returns the product of mass and universal gravitational constant of a Solar System body.
 *
 * For problems involving the gravitational interactions of Solar System bodies,
 * it is helpful to know the product GM, where G = the universal gravitational constant
 * and M = the mass of the body. In practice, GM is known to a higher precision than
 * either G or M alone, and thus using the product results in the most accurate results.
 * This function returns the product GM in the units au^3/day^2.
 * The values come from page 10 of a
 * [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).
 *
 * @param {Body} body
 *      The body for which to find the GM product.
 *      Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet.
 *      Any other value will cause an exception to be thrown.
 *
 * @returns {number}
 *      The mass product of the given body in au^3/day^2.
 */
export function MassProduct(body: Body): number {
    switch (body) {
        case Body.Sun:      return SUN_GM;
        case Body.Mercury:  return MERCURY_GM;
        case Body.Venus:    return VENUS_GM;
        case Body.Earth:    return EARTH_GM;
        case Body.Moon:     return MOON_GM;
        case Body.EMB:      return EARTH_GM + MOON_GM;
        case Body.Mars:     return MARS_GM;
        case Body.Jupiter:  return JUPITER_GM;
        case Body.Saturn:   return SATURN_GM;
        case Body.Uranus:   return URANUS_GM;
        case Body.Neptune:  return NEPTUNE_GM;
        case Body.Pluto:    return PLUTO_GM;
        default:
            throw `Do not know mass product for body: ${body}`;
    }
}


let ob2000: number;   // lazy-evaluated mean obliquity of the ecliptic at J2000, in radians
let cos_ob2000: number;
let sin_ob2000: number;

function VerifyBoolean(b: boolean): boolean {
    if (b !== true && b !== false) {
        console.trace();
        throw `Value is not boolean: ${b}`;
    }
    return b;
}

function VerifyNumber(x: number): number {
    if (!Number.isFinite(x)) {
        console.trace();
        throw `Value is not a finite number: ${x}`;
    }
    return x;
}

function Frac(x: number): number {
    return x - Math.floor(x);
}

/**
 * @brief Calculates the angle in degrees between two vectors.
 *
 * Given a pair of vectors, this function returns the angle in degrees
 * between the two vectors in 3D space.
 * The angle is measured in the plane that contains both vectors.
 *
 * @param {Vector} a
 *      The first of a pair of vectors between which to measure an angle.
 *
 * @param {Vector} b
 *      The second of a pair of vectors between which to measure an angle.
 *
 * @returns {number}
 *      The angle between the two vectors expressed in degrees.
 *      The value is in the range [0, 180].
 */
export function AngleBetween(a: Vector, b: Vector): number {
    const aa = (a.x*a.x + a.y*a.y + a.z*a.z);
    if (Math.abs(aa) < 1.0e-8)
        throw `AngleBetween: first vector is too short.`;

    const bb = (b.x*b.x + b.y*b.y + b.z*b.z);
    if (Math.abs(bb) < 1.0e-8)
        throw `AngleBetween: second vector is too short.`;

    const dot = (a.x*b.x + a.y*b.y + a.z*b.z) / Math.sqrt(aa * bb);

    if (dot <= -1.0)
        return 180;

    if (dot >= +1.0)
        return 0;

    return RAD2DEG * Math.acos(dot);
}

/**
 * @brief String constants that represent the solar system bodies supported by Astronomy Engine.
 *
 * The following strings represent solar system bodies supported by various Astronomy Engine functions.
 * Not every body is supported by every function; consult the documentation for each function
 * to find which bodies it supports.
 *
 * "Sun", "Moon", "Mercury", "Venus", "Earth", "Mars", "Jupiter",
 * "Saturn", "Uranus", "Neptune", "Pluto",
 * "SSB" (Solar System Barycenter),
 * "EMB" (Earth/Moon Barycenter)
 *
 * You can also use enumeration syntax for the bodies, like
 * `Astronomy.Body.Moon`, `Astronomy.Body.Jupiter`, etc.
 *
 * @enum {string}
 */
export enum Body {
    Sun     = 'Sun',
    Moon    = 'Moon',
    Mercury = 'Mercury',
    Venus   = 'Venus',
    Earth   = 'Earth',
    Mars    = 'Mars',
    Jupiter = 'Jupiter',
    Saturn  = 'Saturn',
    Uranus  = 'Uranus',
    Neptune = 'Neptune',
    Pluto   = 'Pluto',
    SSB     = 'SSB',          // Solar System Barycenter
    EMB     = 'EMB'           // Earth/Moon Barycenter
}


enum PrecessDirection {
    From2000,
    Into2000
}


interface PlanetInfo {
    OrbitalPeriod: number;
}

interface PlanetTable {
    [body: string]: PlanetInfo;
}

const Planet: PlanetTable = {
    Mercury: { OrbitalPeriod:    87.969 },
    Venus:   { OrbitalPeriod:   224.701 },
    Earth:   { OrbitalPeriod:   365.256 },
    Mars:    { OrbitalPeriod:   686.980 },
    Jupiter: { OrbitalPeriod:  4332.589 },
    Saturn:  { OrbitalPeriod: 10759.22  },
    Uranus:  { OrbitalPeriod: 30685.4   },
    Neptune: { OrbitalPeriod: 60189.0   },
    Pluto:   { OrbitalPeriod: 90560.0   }
};

type VsopModel = number[][][][];

interface VsopTable {
    [body: string]: VsopModel;
}

const vsop: VsopTable = {
    Mercury: $ASTRO_LIST_VSOP(Mercury),
    Venus:   $ASTRO_LIST_VSOP(Venus),
    Earth:   $ASTRO_LIST_VSOP(Earth),
    Mars:    $ASTRO_LIST_VSOP(Mars),
    Jupiter: $ASTRO_LIST_VSOP(Jupiter),
    Saturn:  $ASTRO_LIST_VSOP(Saturn),
    Uranus:  $ASTRO_LIST_VSOP(Uranus),
    Neptune: $ASTRO_LIST_VSOP(Neptune)
};

export function DeltaT_EspenakMeeus(ut: number): number {
    var u: number, u2: number, u3: number, u4: number, u5: number, u6: number, u7: number;

    /*
        Fred Espenak writes about Delta-T generically here:
        https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
        https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html

        He provides polynomial approximations for distant years here:
        https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html

        They start with a year value 'y' such that y=2000 corresponds
        to the UTC Date 15-January-2000. Convert difference in days
        to mean tropical years.
    */

    const y = 2000 + ((ut - 14) / DAYS_PER_TROPICAL_YEAR);

    if (y < -500) {
        u = (y - 1820) / 100;
        return -20 + (32 * u*u);
    }
    if (y < 500) {
        u = y / 100;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3;
        return 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3 - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6;
    }
    if (y < 1600) {
        u = (y - 1000) / 100;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3;
        return 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3 - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6;
    }
    if (y < 1700) {
        u = y - 1600;
        u2 = u*u; u3 = u*u2;
        return 120 - 0.9808*u - 0.01532*u2 + u3/7129.0;
    }
    if (y < 1800) {
        u = y - 1700;
        u2 = u*u; u3 = u*u2; u4 = u2*u2;
        return 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000;
    }
    if (y < 1860) {
        u = y - 1800;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3; u7 = u3*u4;
        return 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4 + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7;
    }
    if (y < 1900) {
        u = y - 1860;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3;
        return 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3 - 0.0004473624*u4 + u5/233174;
    }
    if (y < 1920) {
        u = y - 1900;
        u2 = u*u; u3 = u*u2; u4 = u2*u2;
        return -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4;
    }
    if (y < 1941) {
        u = y - 1920;
        u2 = u*u; u3 = u*u2;
        return 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3;
    }
    if (y < 1961) {
        u = y - 1950;
        u2 = u*u; u3 = u*u2;
        return 29.07 + 0.407*u - u2/233 + u3/2547;
    }
    if (y < 1986) {
        u = y - 1975;
        u2 = u*u; u3 = u*u2;
        return 45.45 + 1.067*u - u2/260 - u3/718;
    }
    if (y < 2005) {
        u = y - 2000;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3;
        return 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4 + 0.00002373599*u5;
    }
    if (y < 2050) {
        u = y - 2000;
        return 62.92 + 0.32217*u + 0.005589*u*u;
    }
    if (y < 2150) {
        u = (y-1820)/100;
        return -20 + 32*u*u - 0.5628*(2150 - y);
    }

    /* all years after 2150 */
    u = (y - 1820) / 100;
    return -20 + (32 * u*u);
}


export type DeltaTimeFunction = (ut: number) => number;


export function DeltaT_JplHorizons(ut: number): number {
    return DeltaT_EspenakMeeus(Math.min(ut, 17.0 * DAYS_PER_TROPICAL_YEAR));
}

let DeltaT: DeltaTimeFunction = DeltaT_EspenakMeeus;

export function SetDeltaTFunction(func: DeltaTimeFunction) {
    DeltaT = func;
}

/**
 * @ignore
 *
 * @brief Calculates Terrestrial Time (TT) from Universal Time (UT).
 *
 * @param {number} ut
 *      The Universal Time expressed as a floating point number of days since the 2000.0 epoch.
 *
 * @returns {number}
 *      A Terrestrial Time expressed as a floating point number of days since the 2000.0 epoch.
 */
function TerrestrialTime(ut: number): number {
    return ut + DeltaT(ut)/86400;
}

/**
 * @brief The date and time of an astronomical observation.
 *
 * Objects of type `AstroTime` are used throughout the internals
 * of the Astronomy library, and are included in certain return objects.
 * Use the constructor or the {@link MakeTime} function to create an `AstroTime` object.
 *
 * @property {Date} date
 *      The JavaScript Date object for the given date and time.
 *      This Date corresponds to the numeric day value stored in the `ut` property.
 *
 * @property {number} ut
 *      Universal Time (UT1/UTC) in fractional days since the J2000 epoch.
 *      Universal Time represents time measured with respect to the Earth's rotation,
 *      tracking mean solar days.
 *      The Astronomy library approximates UT1 and UTC as being the same thing.
 *      This gives sufficient accuracy for the precision requirements of this project.
 *
 * @property {number} tt
 *      Terrestrial Time in fractional days since the J2000 epoch.
 *      TT represents a continuously flowing ephemeris timescale independent of
 *      any variations of the Earth's rotation, and is adjusted from UT
 *      using a best-fit piecewise polynomial model devised by
 *      [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).
 */
export class AstroTime {
    date: Date;
    ut: number;
    tt: number;

    /**
     * @param {FlexibleDateTime} date
     *      A JavaScript Date object, a numeric UTC value expressed in J2000 days, or another AstroTime object.
     */
    constructor(date: FlexibleDateTime) {
        if (date instanceof AstroTime) {
            // Construct a clone of the AstroTime passed in.
            this.date = date.date;
            this.ut = date.ut;
            this.tt = date.tt;
            return;
        }

        const MillisPerDay = 1000 * 3600 * 24;

        if ((date instanceof Date) && Number.isFinite(date.getTime())) {
            this.date = date;
            this.ut = (date.getTime() - J2000.getTime()) / MillisPerDay;
            this.tt = TerrestrialTime(this.ut);
            return;
        }

        if (Number.isFinite(date)) {
            this.date = new Date(J2000.getTime() + <number>date*MillisPerDay);
            this.ut = <number>date;
            this.tt = TerrestrialTime(this.ut);
            return;
        }

        throw 'Argument must be a Date object, an AstroTime object, or a numeric UTC Julian date.';
    }

    /**
     * @brief Creates an `AstroTime` value from a Terrestrial Time (TT) day value.
     *
     * This function can be used in rare cases where a time must be based
     * on Terrestrial Time (TT) rather than Universal Time (UT).
     * Most developers will want to invoke `new AstroTime(ut)` with a universal time
     * instead of this function, because usually time is based on civil time adjusted
     * by leap seconds to match the Earth's rotation, rather than the uniformly
     * flowing TT used to calculate solar system dynamics. In rare cases
     * where the caller already knows TT, this function is provided to create
     * an `AstroTime` value that can be passed to Astronomy Engine functions.
     *
     * @param {number} tt
     *      The number of days since the J2000 epoch as expressed in Terrestrial Time.
     *
     * @returns {AstroTime}
     *      An `AstroTime` object for the specified terrestrial time.
     */
    static FromTerrestrialTime(tt: number): AstroTime {
        let time = new AstroTime(tt);
        for(;;) {
            const err = tt - time.tt;
            if (Math.abs(err) < 1.0e-12)
                return time;
            time = time.AddDays(err);
        }
    }

    /**
     * Formats an `AstroTime` object as an [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601)
     * date/time string in UTC, to millisecond resolution.
     * Example: `2018-08-17T17:22:04.050Z`
     * @returns {string}
     */
    toString(): string {
        return this.date.toISOString();
    }

    /**
     * Returns a new `AstroTime` object adjusted by the floating point number of days.
     * Does NOT modify the original `AstroTime` object.
     *
     * @param {number} days
     *      The floating point number of days by which to adjust the given date and time.
     *      Positive values adjust the date toward the future, and
     *      negative values adjust the date toward the past.
     *
     * @returns {AstroTime}
     */
    AddDays(days: number): AstroTime {
        // This is slightly wrong, but the error is tiny.
        // We really should be adding to TT, not to UT.
        // But using TT would require creating an inverse function for DeltaT,
        // which would be quite a bit of extra calculation.
        // I estimate the error is in practice on the order of 10^(-7)
        // times the value of 'days'.
        // This is based on a typical drift of 1 second per year between UT and TT.
        return new AstroTime(this.ut + days);
    }
}

function InterpolateTime(time1: AstroTime, time2: AstroTime, fraction: number): AstroTime {
    return new AstroTime(time1.ut + fraction*(time2.ut - time1.ut));
}

/**
 * @brief A `Date`, `number`, or `AstroTime` value that specifies the date and time of an astronomical event.
 *
 * `FlexibleDateTime` is a placeholder type that represents three different types
 * that may be passed to many Astronomy Engine functions: a JavaScript `Date` object,
 * a number representing the real-valued number of UT days since the J2000 epoch,
 * or an {@link AstroTime} object.
 *
 * This flexibility is for convenience of outside callers.
 * Internally, Astronomy Engine always converts a `FlexibleDateTime` parameter
 * to an `AstroTime` object by calling {@link MakeTime}.
 *
 * @typedef {Date | number | AstroTime} FlexibleDateTime
 */

/**
 * @brief Converts multiple date/time formats to `AstroTime` format.
 *
 * Given a Date object or a number days since noon (12:00) on January 1, 2000 (UTC),
 * this function creates an {@link AstroTime} object.
 *
 * Given an {@link AstroTime} object, returns the same object unmodified.
 * Use of this function is not required for any of the other exposed functions in this library,
 * because they all guarantee converting date/time parameters to `AstroTime`
 * as needed. However, it may be convenient for callers who need to understand
 * the difference between UTC and TT (Terrestrial Time). In some use cases,
 * converting once to `AstroTime` format and passing the result into multiple
 * function calls may be more efficient than passing in native JavaScript Date objects.
 *
 * @param {FlexibleDateTime} date
 *      A Date object, a number of UTC days since the J2000 epoch (noon on January 1, 2000),
 *      or an AstroTime object. See remarks above.
 *
 * @returns {AstroTime}
 */
export function MakeTime(date: FlexibleDateTime): AstroTime {
    if (date instanceof AstroTime) {
        return date;
    }
    return new AstroTime(date);
}

const iaudata = [
$ASTRO_IAU_DATA()
];

interface NutationAngles {
    dpsi: number;
    deps: number;
}

function iau2000b(time: AstroTime): NutationAngles {
    var i: number, t: number, el: number, elp: number, f: number, d: number, om: number, arg: number, dp: number, de: number, sarg: number, carg: number;
    var nals: number[], cls: number[];

    function mod(x: number): number {
        return (x % ASEC360) * ASEC2RAD;
    }

    t = time.tt / 36525;
    el  = mod(485868.249036 + t*1717915923.2178);
    elp = mod(1287104.79305 + t*129596581.0481);
    f   = mod(335779.526232 + t*1739527262.8478);
    d   = mod(1072260.70369 + t*1602961601.2090);
    om  = mod(450160.398036 - t*6962890.5431);
    dp = 0;
    de = 0;
    for (i=76; i >= 0; --i) {
        nals = iaudata[i][0];
        cls = iaudata[i][1];
        arg = (nals[0]*el + nals[1]*elp + nals[2]*f + nals[3]*d + nals[4]*om) % PI2;
        sarg = Math.sin(arg);
        carg = Math.cos(arg);
        dp += (cls[0] + cls[1]*t) * sarg + cls[2]*carg;
        de += (cls[3] + cls[4]*t) * carg + cls[5]*sarg;
    }
    return {
        dpsi: -0.000135 + (dp * 1.0e-7),
        deps: +0.000388 + (de * 1.0e-7)
    };
}

function mean_obliq(time: AstroTime): number {
    var t = time.tt / 36525;
    var asec = (
        (((( -  0.0000000434   * t
             -  0.000000576  ) * t
             +  0.00200340   ) * t
             -  0.0001831    ) * t
             - 46.836769     ) * t + 84381.406
    );
    return asec / 3600.0;
}

interface EarthTiltInfo {
    tt: number;
    dpsi: number;
    deps: number;
    ee: number;
    mobl: number;
    tobl: number;
}

var cache_e_tilt: EarthTiltInfo;

function e_tilt(time: AstroTime): EarthTiltInfo {
    if (!cache_e_tilt || Math.abs(cache_e_tilt.tt - time.tt) > 1.0e-6) {
        const nut = iau2000b(time);
        const mean_ob = mean_obliq(time);
        const true_ob = mean_ob + (nut.deps / 3600);
        cache_e_tilt = {
            tt: time.tt,
            dpsi: nut.dpsi,
            deps: nut.deps,
            ee: nut.dpsi * Math.cos(mean_ob * DEG2RAD) / 15,
            mobl: mean_ob,
            tobl: true_ob
        };
    }
    return cache_e_tilt;
}

function ecl2equ_vec(time: AstroTime, pos: ArrayVector): ArrayVector {
    var obl = mean_obliq(time) * DEG2RAD;
    var cos_obl = Math.cos(obl);
    var sin_obl = Math.sin(obl);
    return [
        pos[0],
        pos[1]*cos_obl - pos[2]*sin_obl,
        pos[1]*sin_obl + pos[2]*cos_obl
    ];
}

export let CalcMoonCount = 0;

function CalcMoon(time: AstroTime) {
    ++CalcMoonCount;

    const T = time.tt / 36525;

    interface PascalArray1 {
        min: number;
        array: number[];
    }

    interface PascalArray2 {
        min: number;
        array: PascalArray1[];
    }

    function DeclareArray1(xmin: number, xmax: number): PascalArray1 {
        const array = [];
        let i: number;
        for (i=0; i <= xmax-xmin; ++i) {
            array.push(0);
        }
        return {min:xmin, array:array};
    }

    function DeclareArray2(xmin: number, xmax: number, ymin: number, ymax: number): PascalArray2 {
        const array = [];
        for (let i=0; i <= xmax-xmin; ++i) {
            array.push(DeclareArray1(ymin, ymax));
        }
        return {min:xmin, array:array};
    }

    function ArrayGet2(a: PascalArray2, x: number, y: number) {
        const m = a.array[x - a.min];
        return m.array[y - m.min];
    }

    function ArraySet2(a: PascalArray2, x: number, y: number, v: number) {
        const m = a.array[x - a.min];
        m.array[y - m.min] = v;
    }

    let S: number, MAX: number, ARG: number, FAC: number, I: number, J: number, T2: number, DGAM: number, DLAM: number, N: number, GAM1C: number, SINPI: number, L0: number, L: number, LS: number, F: number, D: number, DL0: number, DL: number, DLS: number, DF: number, DD: number, DS: number;
    let coArray = DeclareArray2(-6, 6, 1, 4);
    let siArray = DeclareArray2(-6, 6, 1, 4);

    function CO(x: number, y: number) {
        return ArrayGet2(coArray, x, y);
    }

    function SI(x: number, y: number) {
        return ArrayGet2(siArray, x, y);
    }

    function SetCO(x: number, y: number, v: number) {
        return ArraySet2(coArray, x, y, v);
    }

    function SetSI(x: number, y: number, v: number) {
        return ArraySet2(siArray, x, y, v);
    }

    type ThetaFunc = (real:number, imag:number) => void;

    function AddThe(c1: number, s1: number, c2: number, s2: number, func:ThetaFunc): void {
        func(c1*c2 - s1*s2, s1*c2 + c1*s2);
    }

    function Sine(phi: number): number {
        return Math.sin(PI2 * phi);
    }

    T2 = T*T;
    DLAM = 0;
    DS = 0;
    GAM1C = 0;
    SINPI = 3422.7000;

    var S1 = Sine(0.19833+0.05611*T);
    var S2 = Sine(0.27869+0.04508*T);
    var S3 = Sine(0.16827-0.36903*T);
    var S4 = Sine(0.34734-5.37261*T);
    var S5 = Sine(0.10498-5.37899*T);
    var S6 = Sine(0.42681-0.41855*T);
    var S7 = Sine(0.14943-5.37511*T);
    DL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6;
    DL  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6;
    DLS =-6.40*S1                                   -1.89*S6;
    DF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7;
    DD  = DL0-DLS;
    DGAM  = (-3332E-9 * Sine(0.59734-5.37261*T)
              -539E-9 * Sine(0.35498-5.37899*T)
               -64E-9 * Sine(0.39943-5.37511*T));

    L0 = PI2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + DL0/ARC;
    L  = PI2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + DL /ARC;
    LS = PI2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + DLS/ARC;
    F  = PI2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + DF /ARC;
    D  = PI2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + DD /ARC;
    for (I=1; I<=4; ++I) {
        switch (I) {
            case 1: ARG=L;  MAX=4; FAC=1.000002208;               break;
            case 2: ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
            case 3: ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
            case 4: ARG=D;  MAX=6; FAC=1.0;                       break;
            default: throw `Internal error: I = ${I}`;      // persuade TypeScript that ARG, ... are all initialized before use.
        }
        SetCO(0, I, 1);
        SetCO(1, I, Math.cos(ARG) * FAC);
        SetSI(0, I, 0);
        SetSI(1, I, Math.sin(ARG) * FAC);
        for (J=2; J<=MAX; ++J) {
            AddThe(CO(J-1,I), SI(J-1,I), CO(1,I), SI(1,I), (c:number, s:number) => (SetCO(J,I,c), SetSI(J,I,s)));
        }
        for (J=1; J<=MAX; ++J) {
            SetCO(-J, I, CO(J, I));
            SetSI(-J, I, -SI(J, I));
        }
    }

    interface ComplexValue {
        x: number;
        y: number;
    }

    function Term(p: number, q: number, r: number, s: number): ComplexValue {
        var result = { x:1, y:0 };
        var I = [ 0, p, q, r, s ];      // I[0] is not used; it is a placeholder
        for (var k=1; k <= 4; ++k)
            if (I[k] !== 0)
                AddThe(result.x, result.y, CO(I[k], k), SI(I[k], k), (c:number, s:number) => (result.x=c, result.y=s));
        return result;
    }

    function AddSol(coeffl: number, coeffs: number, coeffg: number, coeffp: number, p: number, q: number, r: number, s: number): void {
        var result = Term(p, q, r, s);
        DLAM += coeffl * result.y;
        DS += coeffs * result.y;
        GAM1C += coeffg * result.x;
        SINPI += coeffp * result.x;
    }

$ASTRO_ADDSOL()

    function ADDN(coeffn: number, p: number, q: number, r: number, s: number) {
        return coeffn * Term(p, q, r, s).y;
    }

    N = 0;
    N += ADDN(-526.069, 0, 0,1,-2);
    N += ADDN(  -3.352, 0, 0,1,-4);
    N += ADDN( +44.297,+1, 0,1,-2);
    N += ADDN(  -6.000,+1, 0,1,-4);
    N += ADDN( +20.599,-1, 0,1, 0);
    N += ADDN( -30.598,-1, 0,1,-2);
    N += ADDN( -24.649,-2, 0,1, 0);
    N += ADDN(  -2.000,-2, 0,1,-2);
    N += ADDN( -22.571, 0,+1,1,-2);
    N += ADDN( +10.985, 0,-1,1,-2);

    DLAM += (
        +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
        +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
        +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
        +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
        +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
        +0.33*Sine(0.3132   +6.3368*T)
    );

    S = F + DS/ARC;

    let lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*Math.sin(S) - 6.24*Math.sin(3*S) + N;

    return {
        geo_eclip_lon: PI2 * Frac((L0+DLAM/ARC) / PI2),
        geo_eclip_lat: (Math.PI / (180 * 3600)) * lat_seconds,
        distance_au: (ARC * EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI)
    };
}

/**
 * @brief Lunar libration angles, returned by {@link Libration}.
 *
 * @property {number} elat
 *      Sub-Earth libration ecliptic latitude angle, in degrees.
 * @property {number} elon
 *      Sub-Earth libration ecliptic longitude angle, in degrees.
 * @property {number} mlat
 *      Moon's geocentric ecliptic latitude, in degrees.
 * @property {number} mlon
 *      Moon's geocentric ecliptic longitude, in degrees.
 * @property {number} dist_km
 *      Distance between the centers of the Earth and Moon in kilometers.
 * @property {number} diam_deg
 *      The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth.
 */
export class LibrationInfo {
    constructor(
        public elat: number,
        public elon: number,
        public mlat: number,
        public mlon: number,
        public dist_km: number,
        public diam_deg: number
    ) {}
}

/**
 * @brief Calculates the Moon's libration angles at a given moment in time.
 *
 * Libration is an observed back-and-forth wobble of the portion of the
 * Moon visible from the Earth. It is caused by the imperfect tidal locking
 * of the Moon's fixed rotation rate, compared to its variable angular speed
 * of orbit around the Earth.
 *
 * This function calculates a pair of perpendicular libration angles,
 * one representing rotation of the Moon in eclitpic longitude `elon`, the other
 * in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.
 *
 * This function also returns the geocentric position of the Moon
 * expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
 * distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
 * and the apparent angular diameter of the Moon `diam_deg`.
 *
 * @param {FlexibleDateTime} date
 *      A Date object, a number of UTC days since the J2000 epoch (noon on January 1, 2000),
 *      or an AstroTime object.
 *
 * @returns {LibrationInfo}
 */
export function Libration(date: FlexibleDateTime): LibrationInfo {
    const time = MakeTime(date);
    const t = time.tt / 36525.0;
    const t2 = t * t;
    const t3 = t2 * t;
    const t4 = t2 * t2;
    const moon = CalcMoon(time);
    const mlon = moon.geo_eclip_lon;
    const mlat = moon.geo_eclip_lat;
    const dist_km = moon.distance_au * KM_PER_AU;

    // Inclination angle
    const I = DEG2RAD * 1.543;

    // Moon's argument of latitude in radians.
    const f = DEG2RAD * NormalizeLongitude(93.2720950 + 483202.0175233*t - 0.0036539*t2 - t3/3526000 + t4/863310000);

    // Moon's ascending node's mean longitude in radians.
    const omega = DEG2RAD * NormalizeLongitude(125.0445479 - 1934.1362891*t + 0.0020754*t2 + t3/467441 - t4/60616000);

    // Sun's mean anomaly.
    const m = DEG2RAD * NormalizeLongitude(357.5291092 + 35999.0502909*t - 0.0001536*t2 + t3/24490000);

    // Moon's mean anomaly.
    const mdash = DEG2RAD * NormalizeLongitude(134.9633964 + 477198.8675055*t + 0.0087414*t2 + t3/69699 - t4/14712000);

    // Moon's mean elongation.
    const d = DEG2RAD * NormalizeLongitude(297.8501921 + 445267.1114034*t - 0.0018819*t2 + t3/545868 - t4/113065000);

    // Eccentricity of the Earth's orbit.
    const e = 1.0 - 0.002516*t - 0.0000074*t2;

    // Optical librations
    const w = mlon - omega;
    const a = Math.atan2(Math.sin(w)*Math.cos(mlat)*Math.cos(I) - Math.sin(mlat)*Math.sin(I), Math.cos(w)*Math.cos(mlat));
    const ldash = LongitudeOffset(RAD2DEG * (a - f));
    const bdash = Math.asin(-Math.sin(w)*Math.cos(mlat)*Math.sin(I) - Math.sin(mlat)*Math.cos(I));

    // Physical librations
    const k1 = DEG2RAD*(119.75 + 131.849*t);
    const k2 = DEG2RAD*(72.56 + 20.186*t);

    const rho = (
        -0.02752*Math.cos(mdash) +
        -0.02245*Math.sin(f) +
        +0.00684*Math.cos(mdash - 2*f) +
        -0.00293*Math.cos(2*f) +
        -0.00085*Math.cos(2*f - 2*d) +
        -0.00054*Math.cos(mdash - 2*d) +
        -0.00020*Math.sin(mdash + f) +
        -0.00020*Math.cos(mdash + 2*f) +
        -0.00020*Math.cos(mdash - f) +
        +0.00014*Math.cos(mdash + 2*f - 2*d)
    );

    const sigma = (
        -0.02816*Math.sin(mdash) +
        +0.02244*Math.cos(f) +
        -0.00682*Math.sin(mdash - 2*f) +
        -0.00279*Math.sin(2*f) +
        -0.00083*Math.sin(2*f - 2*d) +
        +0.00069*Math.sin(mdash - 2*d) +
        +0.00040*Math.cos(mdash + f) +
        -0.00025*Math.sin(2*mdash) +
        -0.00023*Math.sin(mdash + 2*f) +
        +0.00020*Math.cos(mdash - f) +
        +0.00019*Math.sin(mdash - f) +
        +0.00013*Math.sin(mdash + 2*f - 2*d) +
        -0.00010*Math.cos(mdash - 3*f)
    );

    const tau = (
        +0.02520*e*Math.sin(m) +
        +0.00473*Math.sin(2*mdash - 2*f) +
        -0.00467*Math.sin(mdash) +
        +0.00396*Math.sin(k1) +
        +0.00276*Math.sin(2*mdash - 2*d) +
        +0.00196*Math.sin(omega) +
        -0.00183*Math.cos(mdash - f) +
        +0.00115*Math.sin(mdash - 2*d) +
        -0.00096*Math.sin(mdash - d) +
        +0.00046*Math.sin(2*f - 2*d) +
        -0.00039*Math.sin(mdash - f) +
        -0.00032*Math.sin(mdash - m - d) +
        +0.00027*Math.sin(2*mdash - m - 2*d) +
        +0.00023*Math.sin(k2) +
        -0.00014*Math.sin(2*d) +
        +0.00014*Math.cos(2*mdash - 2*f) +
        -0.00012*Math.sin(mdash - 2*f) +
        -0.00012*Math.sin(2*mdash) +
        +0.00011*Math.sin(2*mdash - 2*m - 2*d)
    );

    const ldash2 = -tau + (rho*Math.cos(a) + sigma*Math.sin(a))*Math.tan(bdash);
    const bdash2 = sigma*Math.cos(a) - rho*Math.sin(a);
    const diam_deg = 2.0 * RAD2DEG * Math.atan(MOON_MEAN_RADIUS_KM / Math.sqrt(dist_km*dist_km - MOON_MEAN_RADIUS_KM*MOON_MEAN_RADIUS_KM));
    return new LibrationInfo(RAD2DEG*bdash + bdash2, ldash + ldash2, RAD2DEG*mlat, RAD2DEG*mlon, dist_km, diam_deg);
}

function rotate(rot: RotationMatrix, vec: ArrayVector): ArrayVector {
    return [
        rot.rot[0][0]*vec[0] + rot.rot[1][0]*vec[1] + rot.rot[2][0]*vec[2],
        rot.rot[0][1]*vec[0] + rot.rot[1][1]*vec[1] + rot.rot[2][1]*vec[2],
        rot.rot[0][2]*vec[0] + rot.rot[1][2]*vec[1] + rot.rot[2][2]*vec[2]
    ];
}

function precession(pos: ArrayVector, time: AstroTime, dir: PrecessDirection): ArrayVector {
    const r = precession_rot(time, dir);
    return rotate(r, pos);
}

function precession_posvel(state: StateVector, time: AstroTime, dir: PrecessDirection): StateVector {
    const r = precession_rot(time, dir);
    return RotateState(r, state);
}

function precession_rot(time: AstroTime, dir: PrecessDirection): RotationMatrix {
    const t = time.tt / 36525;

    let eps0 = 84381.406;

    let psia   = (((((-    0.0000000951  * t
                      +    0.000132851 ) * t
                      -    0.00114045  ) * t
                      -    1.0790069   ) * t
                      + 5038.481507    ) * t);

    let omegaa = (((((+    0.0000003337  * t
                      -    0.000000467 ) * t
                      -    0.00772503  ) * t
                      +    0.0512623   ) * t
                      -    0.025754    ) * t + eps0);

    let chia   = (((((-    0.0000000560  * t
                      +    0.000170663 ) * t
                      -    0.00121197  ) * t
                      -    2.3814292   ) * t
                      +   10.556403    ) * t);

    eps0   *= ASEC2RAD;
    psia   *= ASEC2RAD;
    omegaa *= ASEC2RAD;
    chia   *= ASEC2RAD;

    const sa = Math.sin(eps0);
    const ca = Math.cos(eps0);
    const sb = Math.sin(-psia);
    const cb = Math.cos(-psia);
    const sc = Math.sin(-omegaa);
    const cc = Math.cos(-omegaa);
    const sd = Math.sin(chia);
    const cd = Math.cos(chia);

    const xx =  cd*cb - sb*sd*cc;
    const yx =  cd*sb*ca + sd*cc*cb*ca - sa*sd*sc;
    const zx =  cd*sb*sa + sd*cc*cb*sa + ca*sd*sc;
    const xy = -sd*cb - sb*cd*cc;
    const yy = -sd*sb * ca + cd*cc*cb*ca - sa*cd*sc;
    const zy = -sd*sb * sa + cd*cc*cb*sa + ca*cd*sc;
    const xz =  sb*sc;
    const yz = -sc*cb * ca - sa*cc;
    const zz = -sc*cb * sa + cc*ca;

    if (dir === PrecessDirection.Into2000) {
        // Perform rotation from epoch to J2000.0.
        return new RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ]);
    }

    if (dir === PrecessDirection.From2000) {
        // Perform rotation from J2000.0 to epoch.
        return new RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ]);
    }

    throw 'Invalid precess direction';
}

function era(time: AstroTime): number {    // Earth Rotation Angle
    const thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut;
    const thet3 = time.ut % 1;
    let theta = 360 * ((thet1 + thet3) % 1);
    if (theta < 0) {
        theta += 360;
    }
    return theta;
}

interface SiderealTimeInfo {
    tt: number;
    st: number;
}

let sidereal_time_cache: SiderealTimeInfo;

function sidereal_time(time: AstroTime): number {          // calculates Greenwich Apparent Sidereal Time (GAST)
    if (!sidereal_time_cache || sidereal_time_cache.tt !== time.tt) {
        const t = time.tt / 36525;
        let eqeq = 15 * e_tilt(time).ee;    // Replace with eqeq=0 to get GMST instead of GAST (if we ever need it)
        const theta = era(time);
        const st = (eqeq + 0.014506 +
                 (((( -    0.0000000368   * t
                      -    0.000029956  ) * t
                      -    0.00000044   ) * t
                      +    1.3915817    ) * t
                      + 4612.156534     ) * t);

        let gst = ((st/3600 + theta) % 360) / 15;
        if (gst < 0) {
            gst += 24;
        }
        sidereal_time_cache = {
            tt: time.tt,
            st: gst
        };
    }
    return sidereal_time_cache.st;     // return sidereal hours in the half-open range [0, 24).
}

/**
 * @brief Calculates Greenwich Apparent Sidereal Time (GAST).
 *
 * Given a date and time, this function calculates the rotation of the
 * Earth, represented by the equatorial angle of the Greenwich prime meridian
 * with respect to distant stars (not the Sun, which moves relative to background
 * stars by almost one degree per day).
 * This angle is called Greenwich Apparent Sidereal Time (GAST).
 * GAST is measured in sidereal hours in the half-open range [0, 24).
 * When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
 * corrected at that time for precession and nutation of the Earth's axis.
 * In this context, the "equinox" is the direction in space where the Earth's
 * orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
 * at the location on the Earth's orbit of the (seasonal) March equinox.
 * As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
 * then starts over at 0.
 * To convert to degrees, multiply the return value by 15.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to find GAST.
 *
 * @returns {number}
 */
export function SiderealTime(date: FlexibleDateTime): number {
    const time = MakeTime(date);
    return sidereal_time(time);
}

function inverse_terra(ovec: ArrayVector, st: number): Observer {
    // Convert from AU to kilometers
    const x = ovec[0] * KM_PER_AU;
    const y = ovec[1] * KM_PER_AU;
    const z = ovec[2] * KM_PER_AU;
    const p = Math.hypot(x, y);
    let lon_deg:number, lat_deg:number, height_km:number;
    if (p < 1.0e-6) {
        // Special case: within 1 millimeter of a pole!
        // Use arbitrary longitude, and latitude determined by polarity of z.
        lon_deg = 0;
        lat_deg = (z > 0.0) ? +90 : -90;
        // Elevation is calculated directly from z.
        height_km = Math.abs(z) - EARTH_POLAR_RADIUS_KM;
    } else {
        const stlocl = Math.atan2(y, x);
        // Calculate exact longitude.
        lon_deg = (RAD2DEG * stlocl) - (15.0 * st);
        // Normalize longitude to the range (-180, +180].
        while (lon_deg <= -180)
            lon_deg += 360;
        while (lon_deg > +180)
            lon_deg -= 360;
        // Numerically solve for exact latitude, using Newton's Method.
        // Start with initial latitude estimate, based on a spherical Earth.
        let lat = Math.atan2(z, p);
        let cos:number, sin:number, denom:number;
        for(;;) {
            // Calculate the error function W(lat).
            // We try to find the root of W, meaning where the error is 0.
            cos = Math.cos(lat);
            sin = Math.sin(lat);
            const factor = (EARTH_FLATTENING_SQUARED-1)*EARTH_EQUATORIAL_RADIUS_KM;
            const cos2 = cos*cos;
            const sin2 = sin*sin;
            const radicand = cos2 + EARTH_FLATTENING_SQUARED*sin2;
            denom = Math.sqrt(radicand);
            const W = (factor*sin*cos)/denom - z*cos + p*sin;
            if (Math.abs(W) < 1.0e-12)
                break;  // The error is now negligible
            // Error is still too large. Find the next estimate.
            // Calculate D = the derivative of W with respect to lat.
            const D = factor*((cos2 - sin2)/denom - sin2*cos2*(EARTH_FLATTENING_SQUARED-1)/(factor*radicand)) + z*sin + p*cos;
            lat -= W/D;
        }
        // We now have a solution for the latitude in radians.
        lat_deg = RAD2DEG * lat;
        // Solve for exact height in meters.
        // There are two formulas I can use. Use whichever has the less risky denominator.
        const adjust = EARTH_EQUATORIAL_RADIUS_KM / denom;
        if (Math.abs(sin) > Math.abs(cos))
            height_km = z/sin - EARTH_FLATTENING_SQUARED*adjust;
        else
            height_km = p/cos - adjust;
    }
    return new Observer(lat_deg, lon_deg, 1000*height_km);
}

function terra(observer: Observer, st: number): TerraInfo {
    const phi = observer.latitude * DEG2RAD;
    const sinphi = Math.sin(phi);
    const cosphi = Math.cos(phi);
    const c = 1 / Math.hypot(cosphi, EARTH_FLATTENING*sinphi);
    const s = EARTH_FLATTENING_SQUARED * c;
    const ht_km = observer.height / 1000;
    const ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
    const ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
    const stlocl = (15*st + observer.longitude) * DEG2RAD;
    const sinst = Math.sin(stlocl);
    const cosst = Math.cos(stlocl);
    return {
        pos: [ach*cosphi*cosst/KM_PER_AU, ach*cosphi*sinst/KM_PER_AU, ash*sinphi/KM_PER_AU],
        vel: [-ANGVEL*ach*cosphi*sinst*86400/KM_PER_AU, ANGVEL*ach*cosphi*cosst*86400/KM_PER_AU, 0]
    };
}

function nutation(pos: ArrayVector, time: AstroTime, dir: PrecessDirection): ArrayVector {
    const r = nutation_rot(time, dir);
    return rotate(r, pos);
}

function nutation_posvel(state: StateVector, time: AstroTime, dir: PrecessDirection): StateVector {
    const r = nutation_rot(time, dir);
    return RotateState(r, state);
}

function nutation_rot(time: AstroTime, dir: PrecessDirection): RotationMatrix {
    const tilt = e_tilt(time);
    const oblm = tilt.mobl * DEG2RAD;
    const oblt = tilt.tobl * DEG2RAD;
    const psi  = tilt.dpsi * ASEC2RAD;
    const cobm = Math.cos(oblm);
    const sobm = Math.sin(oblm);
    const cobt = Math.cos(oblt);
    const sobt = Math.sin(oblt);
    const cpsi = Math.cos(psi);
    const spsi = Math.sin(psi);

    const xx =  cpsi;
    const yx = -spsi*cobm;
    const zx = -spsi*sobm;
    const xy =  spsi*cobt;
    const yy =  cpsi*cobm*cobt + sobm*sobt;
    const zy =  cpsi*sobm*cobt - cobm*sobt;
    const xz =  spsi*sobt;
    const yz =  cpsi*cobm*sobt - sobm*cobt;
    const zz =  cpsi*sobm*sobt + cobm*cobt;

    if (dir === PrecessDirection.From2000) {
        // convert J2000 to of-date
        return new RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ]);
    }

    if (dir === PrecessDirection.Into2000) {
        // convert of-date to J2000
        return new RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ]);
    }

    throw 'Invalid precess direction';
}

function gyration(pos: ArrayVector, time: AstroTime, dir: PrecessDirection) {
    // Combine nutation and precession into a single operation I call "gyration".
    // The order they are composed depends on the direction,
    // because both directions are mutual inverse functions.
    return (dir === PrecessDirection.Into2000) ?
        precession(nutation(pos, time, dir), time, dir) :
        nutation(precession(pos, time, dir), time, dir);
}

function gyration_posvel(state: StateVector, time: AstroTime, dir: PrecessDirection) {
    // Combine nutation and precession into a single operation I call "gyration".
    // The order they are composed depends on the direction,
    // because both directions are mutual inverse functions.
    return (dir === PrecessDirection.Into2000) ?
        precession_posvel(nutation_posvel(state, time, dir), time, dir) :
        nutation_posvel(precession_posvel(state, time, dir), time, dir);
}

function geo_pos(time: AstroTime, observer: Observer): ArrayVector {
    const gast = sidereal_time(time);
    const pos = terra(observer, gast).pos;
    return gyration(pos, time, PrecessDirection.Into2000);
}

/**
 * @brief A 3D Cartesian vector with a time attached to it.
 *
 * Holds the Cartesian coordinates of a vector in 3D space,
 * along with the time at which the vector is valid.
 *
 * @property {number} x        The x-coordinate expressed in astronomical units (AU).
 * @property {number} y        The y-coordinate expressed in astronomical units (AU).
 * @property {number} z        The z-coordinate expressed in astronomical units (AU).
 * @property {AstroTime} t     The time at which the vector is valid.
 */
export class Vector {
    constructor(
        public x: number,
        public y: number,
        public z: number,
        public t: AstroTime)
        {}

    /**
     * Returns the length of the vector in astronomical units (AU).
     * @returns {number}
     */
    Length(): number {
        return Math.hypot(this.x, this.y, this.z);
    }
}

/**
 * @brief A combination of a position vector, a velocity vector, and a time.
 *
 * Holds the state vector of a body at a given time, including its position,
 * velocity, and the time they are valid.
 *
 * @property {number} x        The position x-coordinate expressed in astronomical units (AU).
 * @property {number} y        The position y-coordinate expressed in astronomical units (AU).
 * @property {number} z        The position z-coordinate expressed in astronomical units (AU).
 * @property {number} vx       The velocity x-coordinate expressed in AU/day.
 * @property {number} vy       The velocity y-coordinate expressed in AU/day.
 * @property {number} vz       The velocity z-coordinate expressed in AU/day.
 * @property {AstroTime} t     The time at which the vector is valid.
 */
export class StateVector {
    constructor(
        public x: number,
        public y: number,
        public z: number,
        public vx: number,
        public vy: number,
        public vz: number,
        public t: AstroTime)
        {}
}

/**
 * @brief Holds spherical coordinates: latitude, longitude, distance.
 *
 * Spherical coordinates represent the location of
 * a point using two angles and a distance.
 *
 * @property {number} lat       The latitude angle: -90..+90 degrees.
 * @property {number} lon       The longitude angle: 0..360 degrees.
 * @property {number} dist      Distance in AU.
 */
export class Spherical {
    lat: number;
    lon: number;
    dist: number;

    constructor(lat: number, lon: number, dist: number) {
        this.lat  = VerifyNumber(lat);
        this.lon  = VerifyNumber(lon);
        this.dist = VerifyNumber(dist);
    }
}

/**
 * @brief Holds right ascension, declination, and distance of a celestial object.
 *
 * @property {number} ra
 *      Right ascension in sidereal hours: [0, 24).
 *
 * @property {number} dec
 *      Declination in degrees: [-90, +90].
 *
 * @property {number} dist
 *      Distance to the celestial object expressed in
 *      <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU).
 *
 * @property {Vector} vec
 *      The equatorial coordinates in cartesian form, using AU distance units.
 *      x = direction of the March equinox,
 *      y = direction of the June solstice,
 *      z = north.
 */
export class EquatorialCoordinates {
    ra: number;
    dec: number;
    dist: number;
    vec: Vector;

    constructor(ra: number, dec: number, dist: number, vec: Vector) {
        this.ra   = VerifyNumber(ra);
        this.dec  = VerifyNumber(dec);
        this.dist = VerifyNumber(dist);
        this.vec = vec;
    }
}

function IsValidRotationArray(rot: number[][]) {
    if (!(rot instanceof Array) || (rot.length !== 3))
        return false;

    for (let i=0; i < 3; ++i) {
        if (!(rot[i] instanceof Array) || (rot[i].length !== 3))
            return false;

        for (let j=0; j < 3; ++j)
            if (!Number.isFinite(rot[i][j]))
                return false;
    }

    return true;
}

type ArrayVector = [number, number, number];

interface TerraInfo {
    pos: ArrayVector;
    vel: ArrayVector;
}

/**
 * @brief Contains a rotation matrix that can be used to transform one coordinate system to another.
 *
 * @property {number[][]} rot
 *      A normalized 3x3 rotation matrix. For example, the identity matrix is represented
 *      as `[[1, 0, 0], [0, 1, 0], [0, 0, 1]]`.
 */
export class RotationMatrix {
    constructor(public rot: number[][]) {}
}

/**
 * @brief Creates a rotation matrix that can be used to transform one coordinate system to another.
 *
 * This function verifies that the `rot` parameter is of the correct format:
 * a number[3][3] array. It throws an exception if `rot` is not of that shape.
 * Otherwise it creates a new {@link RotationMatrix} object based on `rot`.
 *
 * @param {number[][]} rot
 *      An array [3][3] of numbers. Defines a rotation matrix used to premultiply
 *      a 3D vector to reorient it into another coordinate system.
 *
 * @returns {RotationMatrix}
 */
export function MakeRotation(rot: number[][]) {
    if (!IsValidRotationArray(rot))
        throw 'Argument must be a [3][3] array of numbers';

    return new RotationMatrix(rot);
}

/**
 * @brief Represents the location of an object seen by an observer on the Earth.
 *
 * Holds azimuth (compass direction) and altitude (angle above/below the horizon)
 * of a celestial object as seen by an observer at a particular location on the Earth's surface.
 * Also holds right ascension and declination of the same object.
 * All of these coordinates are optionally adjusted for atmospheric refraction;
 * therefore the right ascension and declination values may not exactly match
 * those found inside a corresponding {@link EquatorialCoordinates} object.
 *
 * @property {number} azimuth
 *      A horizontal compass direction angle in degrees measured starting at north
 *      and increasing positively toward the east.
 *      The value is in the range [0, 360).
 *      North = 0, east = 90, south = 180, west = 270.
 *
 * @property {number} altitude
 *      A vertical angle in degrees above (positive) or below (negative) the horizon.
 *      The value is in the range [-90, +90].
 *      The altitude angle is optionally adjusted upward due to atmospheric refraction.
 *
 * @property {number} ra
 *      The right ascension of the celestial body in sidereal hours.
 *      The value is in the reange [0, 24).
 *      If `altitude` was adjusted for atmospheric reaction, `ra`
 *      is likewise adjusted.
 *
 * @property {number} dec
 *      The declination of of the celestial body in degrees.
 *      The value in the range [-90, +90].
 *      If `altitude` was adjusted for atmospheric reaction, `dec`
 *      is likewise adjusted.
 */
export class HorizontalCoordinates {
    azimuth: number;
    altitude: number;
    ra: number;
    dec: number;

    constructor(azimuth: number, altitude: number, ra: number, dec: number) {
        this.azimuth  = VerifyNumber(azimuth);
        this.altitude = VerifyNumber(altitude);
        this.ra       = VerifyNumber(ra);
        this.dec      = VerifyNumber(dec);
    }
}

/**
 * @brief Ecliptic coordinates of a celestial body.
 *
 * The origin and date of the coordinate system may vary depending on the caller's usage.
 * In general, ecliptic coordinates are measured with respect to the mean plane of the Earth's
 * orbit around the Sun.
 * Includes Cartesian coordinates `(ex, ey, ez)` measured in
 * <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU)
 * and spherical coordinates `(elon, elat)` measured in degrees.
 *
 * @property {Vector} vec
 *      Ecliptic cartesian vector with components measured in astronomical units (AU).
 *      The x-axis is within the ecliptic plane and is oriented in the direction of the
 *      <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>.
 *      The y-axis is within the ecliptic plane and is oriented 90 degrees
 *      counterclockwise from the equinox, as seen from above the Sun's north pole.
 *      The z-axis is oriented perpendicular to the ecliptic plane,
 *      along the direction of the Sun's north pole.
 *
 * @property {number} elat
 *      The ecliptic latitude of the body in degrees.
 *      This is the angle north or south of the ecliptic plane.
 *      The value is in the range [-90, +90].
 *      Positive values are north and negative values are south.
 *
 * @property {number} elon
 *      The ecliptic longitude of the body in degrees.
 *      This is the angle measured counterclockwise around the ecliptic plane,
 *      as seen from above the Sun's north pole.
 *      This is the same direction that the Earth orbits around the Sun.
 *      The angle is measured starting at 0 from the equinox and increases
 *      up to 360 degrees.
 */
export class EclipticCoordinates {
    vec: Vector;
    elat: number;
    elon: number;

    constructor(vec: Vector, elat: number, elon: number) {
        this.vec = vec;
        this.elat = VerifyNumber(elat);
        this.elon = VerifyNumber(elon);
    }
}

function VectorFromArray(av: ArrayVector, time: AstroTime): Vector {
    return new Vector(av[0], av[1], av[2], time);
}

function vector2radec(pos: ArrayVector, time: AstroTime): EquatorialCoordinates {
    const vec = VectorFromArray(pos, time);
    const xyproj = vec.x*vec.x + vec.y*vec.y;
    const dist = Math.sqrt(xyproj + vec.z*vec.z);
    if (xyproj === 0) {
        if (vec.z === 0)
            throw 'Indeterminate sky coordinates';
        return new EquatorialCoordinates(0, (vec.z < 0) ? -90 : +90, dist, vec);
    }

    let ra = RAD2HOUR * Math.atan2(vec.y, vec.x);
    if (ra < 0)
        ra += 24;
    const dec = RAD2DEG * Math.atan2(pos[2], Math.sqrt(xyproj));
    return new EquatorialCoordinates(ra, dec, dist, vec);
}

function spin(angle: number, pos: ArrayVector): ArrayVector {
    const angr = angle * DEG2RAD;
    const c = Math.cos(angr);
    const s = Math.sin(angr);
    return [c*pos[0]+s*pos[1], c*pos[1]-s*pos[0], pos[2]];
}

/**
 * @brief Converts equatorial coordinates to horizontal coordinates.
 *
 * Given a date and time, a geographic location of an observer on the Earth, and
 * equatorial coordinates (right ascension and declination) of a celestial body,
 * returns horizontal coordinates (azimuth and altitude angles) for that body
 * as seen by that observer. Allows optional correction for atmospheric refraction.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to find horizontal coordinates.
 *
 * @param {Observer} observer
 *      The location of the observer for which to find horizontal coordinates.
 *
 * @param {number} ra
 *      Right ascension in sidereal hours of the celestial object,
 *      referred to the mean equinox of date for the J2000 epoch.
 *
 * @param {number} dec
 *      Declination in degrees of the celestial object,
 *      referred to the mean equator of date for the J2000 epoch.
 *      Positive values are north of the celestial equator and negative values are south.
 *
 * @param {string} refraction
 *      If omitted or has a false-like value (false, null, undefined, etc.)
 *      the calculations are performed without any correction for atmospheric
 *      refraction. If the value is the string `"normal"`,
 *      uses the recommended refraction correction based on Meeus "Astronomical Algorithms"
 *      with a linear taper more than 1 degree below the horizon. The linear
 *      taper causes the refraction to linearly approach 0 as the altitude of the
 *      body approaches the nadir (-90 degrees).
 *      If the value is the string `"jplhor"`, uses a JPL Horizons
 *      compatible formula. This is the same algorithm as `"normal"`,
 *      only without linear tapering; this can result in physically impossible
 *      altitudes of less than -90 degrees, which may cause problems for some applications.
 *      (The `"jplhor"` option was created for unit testing against data
 *      generated by JPL Horizons, and is otherwise not recommended for use.)
 *
 * @returns {HorizontalCoordinates}
 */
export function Horizon(date: FlexibleDateTime, observer: Observer, ra: number, dec: number, refraction?: string): HorizontalCoordinates {
    // based on NOVAS equ2hor()
    let time = MakeTime(date);
    VerifyObserver(observer);
    VerifyNumber(ra);
    VerifyNumber(dec);

    const sinlat = Math.sin(observer.latitude * DEG2RAD);
    const coslat = Math.cos(observer.latitude * DEG2RAD);
    const sinlon = Math.sin(observer.longitude * DEG2RAD);
    const coslon = Math.cos(observer.longitude * DEG2RAD);
    const sindc  = Math.sin(dec * DEG2RAD);
    const cosdc  = Math.cos(dec * DEG2RAD);
    const sinra  = Math.sin(ra * HOUR2RAD);
    const cosra  = Math.cos(ra * HOUR2RAD);

    // Calculate three mutually perpendicular unit vectors
    // in equatorial coordinates: uze, une, uwe.
    //
    // uze = The direction of the observer's local zenith (straight up).
    // une = The direction toward due north on the observer's horizon.
    // uwe = The direction toward due west on the observer's horizon.
    //
    // HOWEVER, these are uncorrected for the Earth's rotation due to the time of day.
    //
    // The components of these 3 vectors are as follows:
    // x = direction from center of Earth toward 0 degrees longitude (the prime meridian) on equator.
    // y = direction from center of Earth toward 90 degrees west longitude on equator.
    // z = direction from center of Earth toward the north pole.

    let uze: ArrayVector = [coslat*coslon, coslat*sinlon, sinlat];
    let une: ArrayVector = [-sinlat*coslon, -sinlat*sinlon, coslat];
    let uwe: ArrayVector = [sinlon, -coslon, 0];

    // Correct the vectors uze, une, uwe for the Earth's rotation by calculating
    // sideral time. Call spin() for each uncorrected vector to rotate about
    // the Earth's axis to yield corrected unit vectors uz, un, uw.
    // Multiply sidereal hours by -15 to convert to degrees and flip eastward
    // rotation of the Earth to westward apparent movement of objects with time.

    const spin_angle = -15 * sidereal_time(time);
    let uz = spin(spin_angle, uze);
    let un = spin(spin_angle, une);
    let uw = spin(spin_angle, uwe);

    // Convert angular equatorial coordinates (RA, DEC) to
    // cartesian equatorial coordinates in 'p', using the
    // same orientation system as uze, une, uwe.

    let p = [cosdc*cosra, cosdc*sinra, sindc];

    // Use dot products of p with the zenith, north, and west
    // vectors to obtain the cartesian coordinates of the body in
    // the observer's horizontal orientation system.
    // pz = zenith component [-1, +1]
    // pn = north  component [-1, +1]
    // pw = west   component [-1, +1]

    const pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2];
    const pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2];
    const pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2];

    // proj is the "shadow" of the body vector along the observer's flat ground.
    let proj = Math.hypot(pn, pw);

    // Calculate az = azimuth (compass direction clockwise from East.)
    let az: number;
    if (proj > 0) {
        // If the body is not exactly straight up/down, it has an azimuth.
        // Invert the angle to produce degrees eastward from north.
        az = -RAD2DEG * Math.atan2(pw, pn);
        if (az < 0) az += 360;
    } else {
        // The body is straight up/down, so it does not have an azimuth.
        // Report an arbitrary but reasonable value.
        az = 0;
    }

    // zd = the angle of the body away from the observer's zenith, in degrees.
    let zd = RAD2DEG * Math.atan2(proj, pz);
    let out_ra = ra;
    let out_dec = dec;

    if (refraction) {
        let zd0 = zd;
        let refr = Refraction(refraction, 90-zd);
        zd -= refr;
        if (refr > 0.0 && zd > 3.0e-4) {
            const sinzd = Math.sin(zd * DEG2RAD);
            const coszd = Math.cos(zd * DEG2RAD);
            const sinzd0 = Math.sin(zd0 * DEG2RAD);
            const coszd0 = Math.cos(zd0 * DEG2RAD);
            var pr = [];
            for (let j=0; j<3; ++j) {
                pr.push(((p[j] - coszd0 * uz[j]) / sinzd0)*sinzd + uz[j]*coszd);
            }
            proj = Math.hypot(pr[0], pr[1]);
            if (proj > 0) {
                out_ra = RAD2HOUR * Math.atan2(pr[1], pr[0]);
                if (out_ra < 0) {
                    out_ra += 24;
                }
            } else {
                out_ra = 0;
            }
            out_dec = RAD2DEG * Math.atan2(pr[2], proj);
        }
    }

    return new HorizontalCoordinates(az, 90-zd, out_ra, out_dec);
}


function VerifyObserver(observer: Observer): Observer {
    if (!(observer instanceof Observer)) {
        throw `Not an instance of the Observer class: ${observer}`;
    }
    VerifyNumber(observer.latitude);
    VerifyNumber(observer.longitude);
    VerifyNumber(observer.height);
    if (observer.latitude < -90 || observer.latitude > +90) {
        throw `Latitude ${observer.latitude} is out of range. Must be -90..+90.`;
    }
    return observer;
}


/**
 * @brief Represents the geographic location of an observer on the surface of the Earth.
 *
 * @property {number} latitude
 *      The observer's geographic latitude in degrees north of the Earth's equator.
 *      The value is negative for observers south of the equator.
 *      Must be in the range -90 to +90.
 *
 * @property {number} longitude
 *      The observer's geographic longitude in degrees east of the prime meridian
 *      passing through Greenwich, England.
 *      The value is negative for observers west of the prime meridian.
 *      The value should be kept in the range -180 to +180 to minimize floating point errors.
 *
 * @property {number} height
 *      The observer's elevation above mean sea level, expressed in meters.
 */
export class Observer {
    constructor(
        public latitude: number,
        public longitude: number,
        public height: number)
    {
        VerifyObserver(this);
    }
}

/**
 * @brief Returns apparent geocentric true ecliptic coordinates of date for the Sun.
 *
 * This function is used for calculating the times of equinoxes and solstices.
 *
 * <i>Geocentric</i> means coordinates as the Sun would appear to a hypothetical observer
 * at the center of the Earth.
 * <i>Ecliptic coordinates of date</i> are measured along the plane of the Earth's mean
 * orbit around the Sun, using the
 * <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>
 * of the Earth as adjusted for precession and nutation of the Earth's
 * axis of rotation on the given date.
 *
 * @param {FlexibleDateTime} date
 *      The date and time at which to calculate the Sun's apparent location as seen from
 *      the center of the Earth.
 *
 * @returns {EclipticCoordinates}
 */
export function SunPosition(date: FlexibleDateTime): EclipticCoordinates {
    // Correct for light travel time from the Sun.
    // This is really the same as correcting for aberration.
    // Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes!
    const time = MakeTime(date).AddDays(-1 / C_AUDAY);

    // Get heliocentric cartesian coordinates of Earth in J2000.
    const earth2000 = CalcVsop(vsop.Earth, time);

    // Convert to geocentric location of the Sun.
    const sun2000: ArrayVector = [-earth2000.x, -earth2000.y, -earth2000.z];

    // Convert to equator-of-date equatorial cartesian coordinates.
    const [gx, gy, gz] = gyration(sun2000, time, PrecessDirection.From2000);

    // Convert to ecliptic coordinates of date.
    const true_obliq = DEG2RAD * e_tilt(time).tobl;
    const cos_ob = Math.cos(true_obliq);
    const sin_ob = Math.sin(true_obliq);

    const vec = new Vector(gx, gy, gz, time);
    const sun_ecliptic = RotateEquatorialToEcliptic(vec, cos_ob, sin_ob);
    return sun_ecliptic;
}

/**
 * @brief Calculates equatorial coordinates of a Solar System body at a given time.
 *
 * Returns topocentric equatorial coordinates (right ascension and declination)
 * in one of two different systems: J2000 or true-equator-of-date.
 * Allows optional correction for aberration.
 * Always corrects for light travel time (represents the object as seen by the observer
 * with light traveling to the Earth at finite speed, not where the object is right now).
 * <i>Topocentric</i> refers to a position as seen by an observer on the surface of the Earth.
 * This function corrects for
 * <a href="https://en.wikipedia.org/wiki/Parallax">parallax</a>
 * of the object between a geocentric observer and a topocentric observer.
 * This is most significant for the Moon, because it is so close to the Earth.
 * However, it can have a small effect on the apparent positions of other bodies.
 *
 * @param {Body} body
 *      The body for which to find equatorial coordinates.
 *      Not allowed to be `Body.Earth`.
 *
 * @param {FlexibleDateTime} date
 *      Specifies the date and time at which the body is to be observed.
 *
 * @param {Observer} observer
 *      The location on the Earth of the observer.
 *
 * @param {bool} ofdate
 *      Pass `true` to return equatorial coordinates of date,
 *      i.e. corrected for precession and nutation at the given date.
 *      This is needed to get correct horizontal coordinates when you call
 *      {@link Horizon}.
 *      Pass `false` to return equatorial coordinates in the J2000 system.
 *
 * @param {bool} aberration
 *      Pass `true` to correct for
 *      <a href="https://en.wikipedia.org/wiki/Aberration_of_light">aberration</a>,
 *      or `false` to leave uncorrected.
 *
 * @returns {EquatorialCoordinates}
 *      The topocentric coordinates of the body as adjusted for the given observer.
 */
export function Equator(body: Body, date: FlexibleDateTime, observer: Observer, ofdate: boolean, aberration: boolean): EquatorialCoordinates {
    VerifyObserver(observer);
    VerifyBoolean(ofdate);
    VerifyBoolean(aberration);
    const time = MakeTime(date);
    const gc_observer = geo_pos(time, observer);
    const gc = GeoVector(body, time, aberration);
    const j2000: ArrayVector = [
        gc.x - gc_observer[0],
        gc.y - gc_observer[1],
        gc.z - gc_observer[2]
    ];

    if (!ofdate)
        return vector2radec(j2000, time);

    const datevect = gyration(j2000, time, PrecessDirection.From2000);
    return vector2radec(datevect, time);
}


/**
 * @brief Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.
 *
 * This function calculates a vector from the center of the Earth to
 * a point on or near the surface of the Earth, expressed in equatorial
 * coordinates. It takes into account the rotation of the Earth at the given
 * time, along with the given latitude, longitude, and elevation of the observer.
 *
 * The caller may pass `ofdate` as `true` to return coordinates relative to the Earth's
 * equator at the specified time, or `false` to use the J2000 equator.
 *
 * The returned vector has components expressed in astronomical units (AU).
 * To convert to kilometers, multiply the `x`, `y`, and `z` values by
 * the constant value {@link KM_PER_AU}.
 *
 * The inverse of this function is also available: {@link VectorObserver}.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the observer's position vector.
 *
 * @param {Observer} observer
 *      The geographic location of a point on or near the surface of the Earth.
 *
 * @param {boolean} ofdate
 *      Selects the date of the Earth's equator in which to express the equatorial coordinates.
 *      The caller may pass `false` to use the orientation of the Earth's equator
 *      at noon UTC on January 1, 2000, in which case this function corrects for precession
 *      and nutation of the Earth as it was at the moment specified by the `time` parameter.
 *      Or the caller may pass `true` to use the Earth's equator at `time`
 *      as the orientation.
 *
 * @returns {Vector}
 *      An equatorial vector from the center of the Earth to the specified location
 *      on (or near) the Earth's surface.
 */
export function ObserverVector(date: FlexibleDateTime, observer: Observer, ofdate: boolean): Vector {
    const time = MakeTime(date);
    const gast = sidereal_time(time);
    let ovec = terra(observer, gast).pos;
    if (!ofdate)
        ovec = gyration(ovec, time, PrecessDirection.Into2000);
    return VectorFromArray(ovec, time);
}

/**
 * @brief Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.
 *
 * This function calculates position and velocity vectors of an observer
 * on or near the surface of the Earth, expressed in equatorial
 * coordinates. It takes into account the rotation of the Earth at the given
 * time, along with the given latitude, longitude, and elevation of the observer.
 *
 * The caller may pass `ofdate` as `true` to return coordinates relative to the Earth's
 * equator at the specified time, or `false` to use the J2000 equator.
 *
 * The returned position vector has components expressed in astronomical units (AU).
 * To convert to kilometers, multiply the `x`, `y`, and `z` values by
 * the constant value {@link KM_PER_AU}.
 * The returned velocity vector has components expressed in AU/day.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the observer's position and velocity vectors.
 *
 * @param {Observer} observer
 *      The geographic location of a point on or near the surface of the Earth.
 *
 * @param {boolean} ofdate
 *      Selects the date of the Earth's equator in which to express the equatorial coordinates.
 *      The caller may pass `false` to use the orientation of the Earth's equator
 *      at noon UTC on January 1, 2000, in which case this function corrects for precession
 *      and nutation of the Earth as it was at the moment specified by the `time` parameter.
 *      Or the caller may pass `true` to use the Earth's equator at `time`
 *      as the orientation.
 *
 * @returns {StateVector}
 */
 export function ObserverState(date: FlexibleDateTime, observer: Observer, ofdate: boolean): StateVector {
    const time = MakeTime(date);
    const gast = sidereal_time(time);
    const svec = terra(observer, gast);
    const state = new StateVector(
        svec.pos[0], svec.pos[1], svec.pos[2],
        svec.vel[0], svec.vel[1], svec.vel[2],
        time
    );

    if (!ofdate)
        return gyration_posvel(state, time, PrecessDirection.Into2000);

    return state;
}

/**
 * @brief Calculates the geographic location corresponding to an equatorial vector.
 *
 * This is the inverse function of {@link ObserverVector}.
 * Given a geocentric equatorial vector, it returns the geographic
 * latitude, longitude, and elevation for that vector.
 *
 * @param {Vector} vector
 *      The geocentric equatorial position vector for which to find geographic coordinates.
 *      The components are expressed in Astronomical Units (AU).
 *      You can calculate AU by dividing kilometers by the constant {@link KM_PER_AU}.
 *      The time `vector.t` determines the Earth's rotation.
 *
 * @param {boolean} ofdate
 *      Selects the date of the Earth's equator in which `vector` is expressed.
 *      The caller may select `false` to use the orientation of the Earth's equator
 *      at noon UTC on January 1, 2000, in which case this function corrects for precession
 *      and nutation of the Earth as it was at the moment specified by `vector.t`.
 *      Or the caller may select `true` to use the Earth's equator at `vector.t`
 *      as the orientation.
 *
 * @returns {Observer}
 *      The geographic latitude, longitude, and elevation above sea level
 *      that corresponds to the given equatorial vector.
 */
 export function VectorObserver(vector: Vector, ofdate: boolean): Observer {
    const gast = sidereal_time(vector.t);
    let ovec:ArrayVector = [vector.x, vector.y, vector.z];
    if (!ofdate) {
        ovec = precession(ovec, vector.t, PrecessDirection.From2000);
        ovec = nutation(ovec, vector.t, PrecessDirection.From2000);
    }
    return inverse_terra(ovec, gast);
}

/**
 * @brief Calculates the gravitational acceleration experienced by an observer on the Earth.
 *
 * This function implements the WGS 84 Ellipsoidal Gravity Formula.
 * The result is a combination of inward gravitational acceleration
 * with outward centrifugal acceleration, as experienced by an observer
 * in the Earth's rotating frame of reference.
 * The resulting value increases toward the Earth's poles and decreases
 * toward the equator, consistent with changes of the weight measured
 * by a spring scale of a fixed mass moved to different latitudes and heights
 * on the Earth.
 *
 * @param {number} latitude
 *      The latitude of the observer in degrees north or south of the equator.
 *      By formula symmetry, positive latitudes give the same answer as negative
 *      latitudes, so the sign does not matter.
 *
 * @param {number} height
 *      The height above the sea level geoid in meters.
 *      No range checking is done; however, accuracy is only valid in the
 *      range 0 to 100000 meters.
 *
 * @returns {number}
 *      The effective gravitational acceleration expressed in meters per second squared [m/s^2].
 */
export function ObserverGravity(latitude: number, height: number): number {
    const s = Math.sin(latitude * DEG2RAD);
    const s2 = s*s;
    const g0 = 9.7803253359 * (1.0 + 0.00193185265241*s2) / Math.sqrt(1.0 - 0.00669437999013*s2);
    return g0 * (1.0 - (3.15704e-07 - 2.10269e-09*s2)*height + 7.37452e-14*height*height);
}

function RotateEquatorialToEcliptic(equ: Vector, cos_ob: number, sin_ob: number): EclipticCoordinates {
    // Rotate equatorial vector to obtain ecliptic vector.
    const ex =  equ.x;
    const ey =  equ.y*cos_ob + equ.z*sin_ob;
    const ez = -equ.y*sin_ob + equ.z*cos_ob;

    const xyproj = Math.hypot(ex, ey);
    let elon = 0;
    if (xyproj > 0) {
        elon = RAD2DEG * Math.atan2(ey, ex);
        if (elon < 0) elon += 360;
    }
    let elat = RAD2DEG * Math.atan2(ez, xyproj);
    let ecl = new Vector(ex, ey, ez, equ.t);
    return new EclipticCoordinates(ecl, elat, elon);
}

/**
 * @brief Converts equatorial Cartesian coordinates to ecliptic Cartesian and angular coordinates.
 *
 * Given J2000 equatorial Cartesian coordinates,
 * returns J2000 ecliptic latitude, longitude, and cartesian coordinates.
 * You can call {@link GeoVector} and pass the resulting vector to this function.
 *
 * @param {Vector} equ
 *      A vector in the J2000 equatorial coordinate system.
 *
 * @returns {EclipticCoordinates}
 */
export function Ecliptic(equ: Vector): EclipticCoordinates {
    // Based on NOVAS functions equ2ecl() and equ2ecl_vec().
    if (ob2000 === undefined) {
        // Lazy-evaluate and keep the mean obliquity of the ecliptic at J2000.
        // This way we don't need to crunch the numbers more than once.
        ob2000 = DEG2RAD * e_tilt(MakeTime(J2000)).mobl;
        cos_ob2000 = Math.cos(ob2000);
        sin_ob2000 = Math.sin(ob2000);
    }
    return RotateEquatorialToEcliptic(equ, cos_ob2000, sin_ob2000);
}

/**
 * @brief Calculates equatorial geocentric Cartesian coordinates for the Moon.
 *
 * Given a time of observation, calculates the Moon's position as a vector.
 * The vector gives the location of the Moon's center relative to the Earth's center
 * with x-, y-, and z-components measured in astronomical units.
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 * Based on the Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
 * which in turn derives from E. W. Brown's lunar theories.
 * Adapted from Turbo Pascal code from the book
 * <a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a>
 * by Montenbruck and Pfleger.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the Moon's geocentric position.
 *
 * @returns {Vector}
 */
export function GeoMoon(date: FlexibleDateTime): Vector {
    const time = MakeTime(date);
    const moon = CalcMoon(time);

    // Convert geocentric ecliptic spherical coords to cartesian coords.
    const dist_cos_lat = moon.distance_au * Math.cos(moon.geo_eclip_lat);
    const gepos: ArrayVector = [
        dist_cos_lat * Math.cos(moon.geo_eclip_lon),
        dist_cos_lat * Math.sin(moon.geo_eclip_lon),
        moon.distance_au * Math.sin(moon.geo_eclip_lat)
    ];

    // Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
    const mpos1 = ecl2equ_vec(time, gepos);

    // Convert from mean equinox of date to J2000...
    const mpos2 = precession(mpos1, time, PrecessDirection.Into2000);

    return new Vector(mpos2[0], mpos2[1], mpos2[2], time);
}

/**
 * @brief Calculates spherical ecliptic geocentric position of the Moon.
 *
 * Given a time of observation, calculates the Moon's geocentric position
 * in ecliptic spherical coordinates. Provides the ecliptic latitude and
 * longitude in degrees, and the geocentric distance in astronomical units (AU).
 * The ecliptic longitude is measured relative to the equinox of date.
 *
 * This algorithm is based on the Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
 * which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
 * It is adapted from Turbo Pascal code from the book
 * <a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a>
 * by Montenbruck and Pfleger.
 *
 * To calculate an equatorial J2000 vector instead, use {@link GeoMoon}.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the Moon's position.
 *
 * @returns {Spherical}
 */
export function EclipticGeoMoon(date: FlexibleDateTime): Spherical {
    const time = MakeTime(date);
    const moon = CalcMoon(time);
    return new Spherical(
        moon.geo_eclip_lat * RAD2DEG,
        moon.geo_eclip_lon * RAD2DEG,
        moon.distance_au
    );
}


/**
 * @brief Calculates equatorial geocentric position and velocity of the Moon at a given time.
 *
 * Given a time of observation, calculates the Moon's position and velocity vectors.
 * The position and velocity are of the Moon's center relative to the Earth's center.
 * The position (x, y, z) components are expressed in AU (astronomical units).
 * The velocity (vx, vy, vz) components are expressed in AU/day.
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 * If you need the Moon's position only, and not its velocity,
 * it is much more efficient to use {@link GeoMoon} instead.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the Moon's geocentric state.
 *
 * @returns {StateVector}
 */
export function GeoMoonState(date: FlexibleDateTime): StateVector {
    const time = MakeTime(date);

    // This is a hack, because trying to figure out how to derive a time
    // derivative for CalcMoon() would be extremely painful!
    // Calculate just before and just after the given time.
    // Average to find position, subtract to find velocity.
    const dt = 1.0e-5;   // 0.864 seconds

    const t1 = time.AddDays(-dt);
    const t2 = time.AddDays(+dt);
    const r1 = GeoMoon(t1);
    const r2 = GeoMoon(t2);

    return new StateVector(
        (r1.x + r2.x) / 2,
        (r1.y + r2.y) / 2,
        (r1.z + r2.z) / 2,
        (r2.x - r1.x) / (2 * dt),
        (r2.y - r1.y) / (2 * dt),
        (r2.z - r1.z) / (2 * dt),
        time
    );
}

/**
 * @brief Calculates the geocentric position and velocity of the Earth/Moon barycenter.
 *
 * Given a time of observation, calculates the geocentric position and velocity vectors
 * of the Earth/Moon barycenter (EMB).
 * The position (x, y, z) components are expressed in AU (astronomical units).
 * The velocity (vx, vy, vz) components are expressed in AU/day.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the EMB's geocentric state.
 *
 * @returns {StateVector}
 */
export function GeoEmbState(date: FlexibleDateTime): StateVector {
    const time = MakeTime(date);
    const s = GeoMoonState(time);
    const d = 1.0 + EARTH_MOON_MASS_RATIO;
    return new StateVector(s.x/d, s.y/d, s.z/d, s.vx/d, s.vy/d, s.vz/d, time);
}

function VsopFormula(formula: any, t: number, clamp_angle: boolean): number {
    let tpower = 1;
    let coord = 0;
    for (let series of formula) {
        let sum = 0;
        for (let [ampl, phas, freq] of series)
            sum += ampl * Math.cos(phas + (t * freq));
        let incr = tpower * sum;
        if (clamp_angle)
            incr %= PI2;    // improve precision for longitudes: they can be hundreds of radians
        coord += incr;
        tpower *= t;
    }
    return coord;
}

function VsopDeriv(formula: any, t: number) {
    let tpower = 1;   // t^s
    let dpower = 0;   // t^(s-1)
    let deriv = 0;
    let s = 0;
    for (let series of formula) {
        let sin_sum = 0;
        let cos_sum = 0;
        for (let [ampl, phas, freq] of series) {
            let angle = phas + (t * freq);
            sin_sum += ampl * freq * Math.sin(angle);
            if (s > 0)
                cos_sum += ampl * Math.cos(angle);
        }
        deriv += (s * dpower * cos_sum) - (tpower * sin_sum);
        dpower = tpower;
        tpower *= t;
        ++s;
    }
    return deriv;
}

const DAYS_PER_MILLENNIUM = 365250;
const LON_INDEX = 0;
const LAT_INDEX = 1;
const RAD_INDEX = 2;

function VsopRotate(eclip: ArrayVector): TerseVector {
    // Convert ecliptic cartesian coordinates to equatorial cartesian coordinates.
    return new TerseVector(
        eclip[0] + 0.000000440360*eclip[1] - 0.000000190919*eclip[2],
        -0.000000479966*eclip[0] + 0.917482137087*eclip[1] - 0.397776982902*eclip[2],
        0.397776982902*eclip[1] + 0.917482137087*eclip[2]
    );
}

function VsopSphereToRect(lon: number, lat: number, radius: number): ArrayVector {
    // Convert spherical coordinates to ecliptic cartesian coordinates.
    const r_coslat = radius * Math.cos(lat);
    const coslon = Math.cos(lon);
    const sinlon = Math.sin(lon);
    return [
        r_coslat * coslon,
        r_coslat * sinlon,
        radius * Math.sin(lat)
    ];
}

function CalcVsop(model: VsopModel, time: AstroTime): Vector {
    const t = time.tt / DAYS_PER_MILLENNIUM;   // millennia since 2000
    const lon = VsopFormula(model[LON_INDEX], t, true);
    const lat = VsopFormula(model[LAT_INDEX], t, false);
    const rad = VsopFormula(model[RAD_INDEX], t, false);
    const eclip = VsopSphereToRect(lon, lat, rad);
    return VsopRotate(eclip).ToAstroVector(time);
}

function CalcVsopPosVel(model: VsopModel, tt: number): body_state_t {
    const t = tt / DAYS_PER_MILLENNIUM;

    // Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
    const lon = VsopFormula(model[LON_INDEX], t, true);
    const lat = VsopFormula(model[LAT_INDEX], t, false);
    const rad = VsopFormula(model[RAD_INDEX], t, false);

    const dlon_dt = VsopDeriv(model[LON_INDEX], t);
    const dlat_dt = VsopDeriv(model[LAT_INDEX], t);
    const drad_dt = VsopDeriv(model[RAD_INDEX], t);

    // Use spherical coords and spherical derivatives to calculate
    // the velocity vector in rectangular coordinates.

    const coslon = Math.cos(lon);
    const sinlon = Math.sin(lon);
    const coslat = Math.cos(lat);
    const sinlat = Math.sin(lat);

    const vx = (
        + (drad_dt * coslat * coslon)
        - (rad * sinlat * coslon * dlat_dt)
        - (rad * coslat * sinlon * dlon_dt)
    );

    const vy = (
        + (drad_dt * coslat * sinlon)
        - (rad * sinlat * sinlon * dlat_dt)
        + (rad * coslat * coslon * dlon_dt)
    );

    const vz = (
        + (drad_dt * sinlat)
        + (rad * coslat * dlat_dt)
    );

    const eclip_pos = VsopSphereToRect(lon, lat, rad);

    // Convert speed units from [AU/millennium] to [AU/day].
    const eclip_vel: ArrayVector = [
        vx / DAYS_PER_MILLENNIUM,
        vy / DAYS_PER_MILLENNIUM,
        vz / DAYS_PER_MILLENNIUM
    ];

    // Rotate the vectors from ecliptic to equatorial coordinates.
    const equ_pos = VsopRotate(eclip_pos);
    const equ_vel = VsopRotate(eclip_vel);
    return new body_state_t(tt, equ_pos, equ_vel);
}

function AdjustBarycenter(ssb: Vector, time: AstroTime, body: Body, pmass: number): void {
    const shift = pmass / (pmass + SUN_GM);
    const planet = CalcVsop(vsop[body], time);
    ssb.x += shift * planet.x;
    ssb.y += shift * planet.y;
    ssb.z += shift * planet.z;
}

function CalcSolarSystemBarycenter(time: AstroTime): Vector {
    const ssb = new Vector(0.0, 0.0, 0.0, time);
    AdjustBarycenter(ssb, time, Body.Jupiter, JUPITER_GM);
    AdjustBarycenter(ssb, time, Body.Saturn,  SATURN_GM);
    AdjustBarycenter(ssb, time, Body.Uranus,  URANUS_GM);
    AdjustBarycenter(ssb, time, Body.Neptune, NEPTUNE_GM);
    return ssb;
}

// Pluto integrator begins ----------------------------------------------------

$ASTRO_PLUTO_TABLE()

class TerseVector {
    constructor(
        public x: number,
        public y: number,
        public z: number)
        {}

    clone(): TerseVector {
        return new TerseVector(this.x, this.y, this.z);
    }

    ToAstroVector(t: AstroTime) {
        return new Vector(this.x, this.y, this.z, t);
    }

    public static zero(): TerseVector {
        return new TerseVector(0, 0, 0);
    }

    quadrature() {
        return this.x*this.x + this.y*this.y + this.z*this.z;
    }

    add(other: TerseVector): TerseVector {
        return new TerseVector(this.x + other.x, this.y + other.y, this.z + other.z);
    }

    sub(other: TerseVector): TerseVector {
        return new TerseVector(this.x - other.x, this.y - other.y, this.z - other.z);
    }

    incr(other: TerseVector) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;
    }

    decr(other: TerseVector) {
        this.x -= other.x;
        this.y -= other.y;
        this.z -= other.z;
    }

    mul(scalar: number): TerseVector {
        return new TerseVector(scalar * this.x, scalar * this.y, scalar * this.z);
    }

    div(scalar: number): TerseVector {
        return new TerseVector(this.x / scalar, this.y / scalar, this.z / scalar);
    }

    mean(other: TerseVector): TerseVector {
        return new TerseVector(
            (this.x + other.x) / 2,
            (this.y + other.y) / 2,
            (this.z + other.z) / 2
        );
    }

    neg(): TerseVector {
        return new TerseVector(-this.x, -this.y, -this.z);
    }
}

class body_state_t {
    constructor(
        public tt: number,
        public r: TerseVector,
        public v: TerseVector)
        {}

    clone(): body_state_t {
        return new body_state_t(this.tt, this.r, this.v);
    }

    sub(other: body_state_t): body_state_t {
        return new body_state_t(this.tt, this.r.sub(other.r), this.v.sub(other.v));
    }
}

type BodyStateTableEntry = [number, ArrayVector, ArrayVector];

function BodyStateFromTable(entry: BodyStateTableEntry): body_state_t {
    let [ tt, [rx, ry, rz], [vx, vy, vz] ] = entry;
    return new body_state_t(tt, new TerseVector(rx, ry, rz), new TerseVector(vx, vy, vz));
}

function AdjustBarycenterPosVel(ssb: body_state_t, tt: number, body: Body, planet_gm: number): body_state_t {
    const shift = planet_gm / (planet_gm + SUN_GM);
    const planet = CalcVsopPosVel(vsop[body], tt);
    ssb.r.incr(planet.r.mul(shift));
    ssb.v.incr(planet.v.mul(shift));
    return planet;
}

function AccelerationIncrement(small_pos: TerseVector, gm: number, major_pos: TerseVector): TerseVector {
    const delta = major_pos.sub(small_pos);
    const r2 = delta.quadrature();
    return delta.mul(gm / (r2 * Math.sqrt(r2)));
}

class major_bodies_t {
    Jupiter: body_state_t;
    Saturn: body_state_t;
    Uranus: body_state_t;
    Neptune: body_state_t;
    Sun: body_state_t;

    constructor(tt: number) {
        // Accumulate the Solar System Barycenter position.

        let ssb = new body_state_t(tt, new TerseVector(0, 0, 0), new TerseVector(0, 0, 0));

        this.Jupiter = AdjustBarycenterPosVel(ssb, tt, Body.Jupiter, JUPITER_GM);
        this.Saturn  = AdjustBarycenterPosVel(ssb, tt, Body.Saturn,  SATURN_GM);
        this.Uranus  = AdjustBarycenterPosVel(ssb, tt, Body.Uranus,  URANUS_GM);
        this.Neptune = AdjustBarycenterPosVel(ssb, tt, Body.Neptune, NEPTUNE_GM);

        // Convert planets' [pos, vel] vectors from heliocentric to barycentric.

        this.Jupiter.r.decr(ssb.r);
        this.Jupiter.v.decr(ssb.v);

        this.Saturn.r.decr(ssb.r);
        this.Saturn.v.decr(ssb.v);

        this.Uranus.r.decr(ssb.r);
        this.Uranus.v.decr(ssb.v);

        this.Neptune.r.decr(ssb.r);
        this.Neptune.v.decr(ssb.v);

        // Convert heliocentric SSB to barycentric Sun.
        this.Sun = new body_state_t(tt, ssb.r.mul(-1), ssb.v.mul(-1));
    }

    Acceleration(pos: TerseVector): TerseVector {
        // Use barycentric coordinates of the Sun and major planets to calculate
        // the gravitational acceleration vector experienced at location 'pos'.
        let acc = AccelerationIncrement(pos, SUN_GM,     this.Sun.r);
        acc.incr(AccelerationIncrement (pos, JUPITER_GM, this.Jupiter.r));
        acc.incr(AccelerationIncrement (pos, SATURN_GM,  this.Saturn.r));
        acc.incr(AccelerationIncrement (pos, URANUS_GM,  this.Uranus.r));
        acc.incr(AccelerationIncrement (pos, NEPTUNE_GM, this.Neptune.r));
        return acc;
    }
}

/**
 * @ignore
 *
 * @brief The state of a body at an incremental step in a gravity simulation.
 *
 * This is an internal data structure used to represent the
 * position, velocity, and acceleration vectors of a body
 * in a gravity simulation at a given moment in time.
 *
 * @property tt
 *      The J2000 terrestrial time of the state [days].
 *
 * @property r
 *      The position vector [au].
 *
 * @property v
 *      The velocity vector [au/day].
 *
 * @property a
 *      The acceleration vector [au/day^2].
 */
class body_grav_calc_t {
    constructor(
        public tt: number,
        public r: TerseVector,
        public v: TerseVector,
        public a: TerseVector)
        {}

    clone(): body_grav_calc_t {
        return new body_grav_calc_t(this.tt, this.r.clone(), this.v.clone(), this.a.clone());
    }
}

class grav_sim_t {
    constructor(
        public bary: major_bodies_t,
        public grav: body_grav_calc_t)
        {}
}

function UpdatePosition(dt: number, r: TerseVector, v: TerseVector, a: TerseVector): TerseVector {
    return new TerseVector(
        r.x + dt*(v.x + dt*a.x/2),
        r.y + dt*(v.y + dt*a.y/2),
        r.z + dt*(v.z + dt*a.z/2)
    );
}

function UpdateVelocity(dt: number, v: TerseVector, a: TerseVector): TerseVector {
    return new TerseVector(
        v.x + dt*a.x,
        v.y + dt*a.y,
        v.z + dt*a.z,
    );
}

function GravSim(tt2: number, calc1: body_grav_calc_t): grav_sim_t {
    const dt = tt2 - calc1.tt;

    // Calculate where the major bodies (Sun, Jupiter...Neptune) will be at tt2.
    const bary2 = new major_bodies_t(tt2);

    // Estimate position of small body as if current acceleration applies across the whole time interval.
    const approx_pos = UpdatePosition(dt, calc1.r, calc1.v, calc1.a);

    // Calculate the average acceleration of the endpoints.
    // This becomes our estimate of the mean effective acceleration over the whole interval.
    const mean_acc = bary2.Acceleration(approx_pos).mean(calc1.a);

    // Refine the estimates of [pos, vel, acc] at tt2 using the mean acceleration.
    const pos = UpdatePosition(dt, calc1.r, calc1.v, mean_acc);
    const vel = calc1.v.add(mean_acc.mul(dt));
    const acc = bary2.Acceleration(pos);
    const grav = new body_grav_calc_t(tt2, pos, vel, acc);
    return new grav_sim_t(bary2, grav);
}

const pluto_cache: body_grav_calc_t[][] = [];

function ClampIndex(frac: number, nsteps: number) {
    const index = Math.floor(frac);
    if (index < 0)
        return 0;
    if (index >= nsteps)
        return nsteps-1;
    return index;
}

function GravFromState(entry: BodyStateTableEntry) {
    const state = BodyStateFromTable(entry);
    const bary = new major_bodies_t(state.tt);
    const r = state.r.add(bary.Sun.r);
    const v = state.v.add(bary.Sun.v);
    const a = bary.Acceleration(r);
    const grav = new body_grav_calc_t(state.tt, r, v, a);
    return new grav_sim_t(bary, grav);
}

function GetSegment(cache: body_grav_calc_t[][], tt: number): body_grav_calc_t[] | null {
    const t0: number = PlutoStateTable[0][0];

    if (tt < t0 || tt > PlutoStateTable[PLUTO_NUM_STATES-1][0]) {
        // Don't bother calculating a segment. Let the caller crawl backward/forward to this time.
        return null;
    }

    const seg_index = ClampIndex((tt - t0) / PLUTO_TIME_STEP, PLUTO_NUM_STATES-1);
    if (!cache[seg_index]) {
        const seg : body_grav_calc_t[] = cache[seg_index] = [];

        // Each endpoint is exact.
        seg[0] = GravFromState(PlutoStateTable[seg_index]).grav;
        seg[PLUTO_NSTEPS-1] = GravFromState(PlutoStateTable[seg_index + 1]).grav;

        // Simulate forwards from the lower time bound.
        let i: number;
        let step_tt = seg[0].tt;
        for (i=1; i < PLUTO_NSTEPS-1; ++i)
            seg[i] = GravSim(step_tt += PLUTO_DT, seg[i-1]).grav;

        // Simulate backwards from the upper time bound.
        step_tt = seg[PLUTO_NSTEPS-1].tt;
        var reverse = [];
        reverse[PLUTO_NSTEPS-1] = seg[PLUTO_NSTEPS-1];
        for (i=PLUTO_NSTEPS-2; i > 0; --i)
            reverse[i] = GravSim(step_tt -= PLUTO_DT, reverse[i+1]).grav;

        // Fade-mix the two series so that there are no discontinuities.
        for (i=PLUTO_NSTEPS-2; i > 0; --i) {
            const ramp = i / (PLUTO_NSTEPS-1);
            seg[i].r = seg[i].r.mul(1 - ramp).add(reverse[i].r.mul(ramp));
            seg[i].v = seg[i].v.mul(1 - ramp).add(reverse[i].v.mul(ramp));
            seg[i].a = seg[i].a.mul(1 - ramp).add(reverse[i].a.mul(ramp));
        }
    }

    return cache[seg_index];
}

function CalcPlutoOneWay(entry: BodyStateTableEntry, target_tt: number, dt: number): grav_sim_t {
    let sim = GravFromState(entry);
    const n = Math.ceil((target_tt - sim.grav.tt) / dt);
    for (let i=0; i < n; ++i)
        sim = GravSim((i+1 === n) ? target_tt : (sim.grav.tt + dt), sim.grav);
    return sim;
}

function CalcPluto(time: AstroTime, helio: boolean): StateVector {
    let r, v, bary;
    const seg = GetSegment(pluto_cache, time.tt);
    if (!seg) {
        // The target time is outside the year range 0000..4000.
        // Calculate it by crawling backward from 0000 or forward from 4000.
        // FIXFIXFIX - This is super slow. Could optimize this with extra caching if needed.
        let sim;
        if (time.tt < PlutoStateTable[0][0])
            sim = CalcPlutoOneWay(PlutoStateTable[0], time.tt, -PLUTO_DT);
        else
            sim = CalcPlutoOneWay(PlutoStateTable[PLUTO_NUM_STATES-1], time.tt, +PLUTO_DT);
        r = sim.grav.r;
        v = sim.grav.v;
        bary = sim.bary;
    } else {
        const left = ClampIndex((time.tt - seg[0].tt) / PLUTO_DT, PLUTO_NSTEPS-1);
        const s1 = seg[left];
        const s2 = seg[left+1];

        // Find mean acceleration vector over the interval.
        const acc = s1.a.mean(s2.a);

        // Use Newtonian mechanics to extrapolate away from t1 in the positive time direction.
        const ra = UpdatePosition(time.tt - s1.tt, s1.r, s1.v, acc);
        const va = UpdateVelocity(time.tt - s1.tt, s1.v, acc);

        // Use Newtonian mechanics to extrapolate away from t2 in the negative time direction.
        const rb = UpdatePosition(time.tt - s2.tt, s2.r, s2.v, acc);
        const vb = UpdateVelocity(time.tt - s2.tt, s2.v, acc);

        // Use fade in/out idea to blend the two position estimates.
        const ramp = (time.tt - s1.tt)/PLUTO_DT;
        r = ra.mul(1 - ramp).add(rb.mul(ramp));
        v = va.mul(1 - ramp).add(vb.mul(ramp));
    }

    if (helio) {
        // Convert barycentric vectors to heliocentric vectors.
        if (!bary)
            bary = new major_bodies_t(time.tt);
        r = r.sub(bary.Sun.r);
        v = v.sub(bary.Sun.v);
    }

    return new StateVector(r.x, r.y, r.z, v.x, v.y, v.z, time);
}

// Pluto integrator ends -----------------------------------------------------

// Jupiter Moons begins ------------------------------------------------------

type jm_term_t = [number, number, number];    // amplitude, phase, frequency
type jm_series_t = jm_term_t[];

interface jupiter_moon_t {
    mu: number;
    al: [number, number];
    a: jm_series_t;
    l: jm_series_t;
    z: jm_series_t;
    zeta: jm_series_t;
};

$ASTRO_JUPITER_MOONS();

/**
 * @brief Holds the positions and velocities of Jupiter's major 4 moons.
 *
 * The {@link JupiterMoons} function returns an object of this type
 * to report position and velocity vectors for Jupiter's largest 4 moons
 * Io, Europa, Ganymede, and Callisto. Each position vector is relative
 * to the center of Jupiter. Both position and velocity are oriented in
 * the EQJ system (that is, using Earth's equator at the J2000 epoch).
 * The positions are expressed in astronomical units (AU),
 * and the velocities in AU/day.
 *
 * @property {StateVector} io
 *      The position and velocity of Jupiter's moon Io.
 *
 * @property {StateVector} europa
 *      The position and velocity of Jupiter's moon Europa.
 *
 * @property {StateVector} ganymede
 *      The position and velocity of Jupiter's moon Ganymede.
 *
 * @property {StateVector} callisto
 *      The position and velocity of Jupiter's moon Callisto.
 */
export class JupiterMoonsInfo {
    constructor(
        public io: StateVector,
        public europa: StateVector,
        public ganymede: StateVector,
        public callisto: StateVector)
        {}
}

function JupiterMoon_elem2pv(
    time: AstroTime,
    mu: number,
    elem: [number, number, number, number, number, number]): StateVector {

    // Translation of FORTRAN subroutine ELEM2PV from:
    // https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

    const A  = elem[0];
    const AL = elem[1];
    const K  = elem[2];
    const H  = elem[3];
    const Q  = elem[4];
    const P  = elem[5];

    const AN = Math.sqrt(mu / (A*A*A));

    let CE: number, SE: number, DE: number;
    let EE = AL + K*Math.sin(AL) - H*Math.cos(AL);
    do {
        CE = Math.cos(EE);
        SE = Math.sin(EE);
        DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE);
        EE += DE;
    }
    while (Math.abs(DE) >= 1.0e-12);

    CE = Math.cos(EE);
    SE = Math.sin(EE);
    const DLE = H*CE - K*SE;
    const RSAM1 = -K*CE - H*SE;
    const ASR = 1.0/(1.0 + RSAM1);
    const PHI = Math.sqrt(1.0 - K*K - H*H);
    const PSI = 1.0/(1.0 + PHI);
    const X1 = A*(CE - K - PSI*H*DLE);
    const Y1 = A*(SE - H + PSI*K*DLE);
    const VX1 = AN*ASR*A*(-SE - PSI*H*RSAM1);
    const VY1 = AN*ASR*A*(+CE + PSI*K*RSAM1);
    const F2 = 2.0*Math.sqrt(1.0 - Q*Q - P*P);
    const P2 = 1.0 - 2.0*P*P;
    const Q2 = 1.0 - 2.0*Q*Q;
    const PQ = 2.0*P*Q;

    return new StateVector(
        X1*P2 + Y1*PQ,
        X1*PQ + Y1*Q2,
        (Q*Y1 - X1*P)*F2,
        VX1*P2 + VY1*PQ,
        VX1*PQ + VY1*Q2,
        (Q*VY1 - VX1*P)*F2,
        time
    );
}

function CalcJupiterMoon(time: AstroTime, m: jupiter_moon_t): StateVector {
    // This is a translation of FORTRAN code by Duriez, Lainey, and Vienne:
    // https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

    const t = time.tt + 18262.5;     // number of days since 1950-01-01T00:00:00Z

    // Calculate 6 orbital elements at the given time t
    const elem: [number, number, number, number, number, number] = [0, m.al[0] + (t * m.al[1]), 0, 0, 0, 0];

    for (let [amplitude, phase, frequency] of m.a)
        elem[0] += amplitude * Math.cos(phase + (t * frequency));

    for (let [amplitude, phase, frequency] of m.l)
        elem[1] += amplitude * Math.sin(phase + (t * frequency));

    elem[1] %= PI2;
    if (elem[1] < 0)
        elem[1] += PI2;

    for (let [amplitude, phase, frequency] of m.z) {
        const arg = phase + (t * frequency);
        elem[2] += amplitude * Math.cos(arg);
        elem[3] += amplitude * Math.sin(arg);
    }

    for (let [amplitude, phase, frequency] of m.zeta) {
        const arg = phase + (t * frequency);
        elem[4] += amplitude * Math.cos(arg);
        elem[5] += amplitude * Math.sin(arg);
    }

    // Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP).
    const state = JupiterMoon_elem2pv(time, m.mu, elem);

    // Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ).
    return RotateState(Rotation_JUP_EQJ, state);
}

/**
 * @brief Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.
 *
 * Calculates position and velocity vectors for Jupiter's moons
 * Io, Europa, Ganymede, and Callisto, at the given date and time.
 * The vectors are jovicentric (relative to the center of Jupiter).
 * Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
 * The position components are expressed in astronomical units (AU), and the
 * velocity components are in AU/day.
 *
 * To convert to heliocentric vectors, call {@link HelioVector}
 * with `Astronomy.Body.Jupiter` to get Jupiter's heliocentric position, then
 * add the jovicentric vectors. Likewise, you can call {@link GeoVector}
 * to convert to geocentric vectors.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate Jupiter's moons.
 *
 * @return {JupiterMoonsInfo}
 *      Position and velocity vectors of Jupiter's largest 4 moons.
 */
export function JupiterMoons(date: FlexibleDateTime): JupiterMoonsInfo {
    const time = new AstroTime(date);
    return new JupiterMoonsInfo(
        CalcJupiterMoon(time, JupiterMoonModel[0]),
        CalcJupiterMoon(time, JupiterMoonModel[1]),
        CalcJupiterMoon(time, JupiterMoonModel[2]),
        CalcJupiterMoon(time, JupiterMoonModel[3])
    );
}

// Jupiter Moons ends --------------------------------------------------------


/**
 * @brief Calculates a vector from the center of the Sun to the given body at the given time.
 *
 * Calculates heliocentric (i.e., with respect to the center of the Sun)
 * Cartesian coordinates in the J2000 equatorial system of a celestial
 * body at a specified time. The position is not corrected for light travel time or aberration.
 *
 * @param {Body} body
 *      One of the following values:
 *      `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`,
 *      `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`,
 *      `Body.Uranus`, `Body.Neptune`, `Body.Pluto`,
 *      `Body.SSB`, or `Body.EMB`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which the body's position is to be calculated.
 *
 * @returns {Vector}
 */
export function HelioVector(body: Body, date: FlexibleDateTime): Vector {
    var time = MakeTime(date);
    if (body in vsop)
        return CalcVsop(vsop[body], time);
    if (body === Body.Pluto) {
        const p = CalcPluto(time, true);
        return new Vector(p.x, p.y, p.z, time);
    }
    if (body === Body.Sun)
        return new Vector(0, 0, 0, time);
    if (body === Body.Moon) {
        var e = CalcVsop(vsop.Earth, time);
        var m = GeoMoon(time);
        return new Vector(e.x+m.x, e.y+m.y, e.z+m.z, time);
    }
    if (body === Body.EMB) {
        const e = CalcVsop(vsop.Earth, time);
        const m = GeoMoon(time);
        const denom = 1.0 + EARTH_MOON_MASS_RATIO;
        return new Vector(e.x+(m.x/denom), e.y+(m.y/denom), e.z+(m.z/denom), time);
    }
    if (body === Body.SSB)
        return CalcSolarSystemBarycenter(time);
    throw `HelioVector: Unknown body "${body}"`;
};

/**
 * @brief Calculates the distance between a body and the Sun at a given time.
 *
 * Given a date and time, this function calculates the distance between
 * the center of `body` and the center of the Sun.
 * For the planets Mercury through Neptune, this function is significantly
 * more efficient than calling {@link HelioVector} followed by taking the length
 * of the resulting vector.
 *
 * @param {Body} body
 *      A body for which to calculate a heliocentric distance:
 *      the Sun, Moon, or any of the planets.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the heliocentric distance.
 *
 * @returns {number}
 *      The heliocentric distance in AU.
 */
export function HelioDistance(body: Body, date: FlexibleDateTime): number {
    const time = MakeTime(date);
    if (body in vsop)
        return VsopFormula(vsop[body][RAD_INDEX], time.tt / DAYS_PER_MILLENNIUM, false);
    return HelioVector(body, time).Length();
}

/**
 * @brief Calculates a vector from the center of the Earth to the given body at the given time.
 *
 * Calculates geocentric (i.e., with respect to the center of the Earth)
 * Cartesian coordinates in the J2000 equatorial system of a celestial
 * body at a specified time. The position is always corrected for light travel time:
 * this means the position of the body is "back-dated" based on how long it
 * takes light to travel from the body to an observer on the Earth.
 * Also, the position can optionally be corrected for aberration, an effect
 * causing the apparent direction of the body to be shifted based on
 * transverse movement of the Earth with respect to the rays of light
 * coming from that body.
 *
 * @param {Body} body
 *      One of the following values:
 *      `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`,
 *      `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`,
 *      `Body.Uranus`, `Body.Neptune`, or `Body.Pluto`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which the body's position is to be calculated.
 *
 * @param {bool} aberration
 *      Pass `true` to correct for
 *      <a href="https://en.wikipedia.org/wiki/Aberration_of_light">aberration</a>,
 *      or `false` to leave uncorrected.
 *
 * @returns {Vector}
 */
export function GeoVector(body: Body, date: FlexibleDateTime, aberration: boolean): Vector {
    VerifyBoolean(aberration);
    const time = MakeTime(date);
    if (body === Body.Moon)
        return GeoMoon(time);
    if (body === Body.Earth)
        return new Vector(0, 0, 0, time);

    let earth: Vector | null = null;
    let h: Vector;
    let geo: Vector;
    let dt: number = 0;
    let ltime = time;

    // Correct for light-travel time, to get position of body as seen from Earth's center.
    for (let iter=0; iter < 10; ++iter) {
        h = HelioVector(body, ltime);

        if (aberration) {
            /*
                Include aberration, so make a good first-order approximation
                by backdating the Earth's position also.
                This is confusing, but it works for objects within the Solar System
                because the distance the Earth moves in that small amount of light
                travel time (a few minutes to a few hours) is well approximated
                by a line segment that substends the angle seen from the remote
                body viewing Earth. That angle is pretty close to the aberration
                angle of the moving Earth viewing the remote body.
                In other words, both of the following approximate the aberration angle:
                    (transverse distance Earth moves) / (distance to body)
                    (transverse speed of Earth) / (speed of light).
            */
            earth = CalcVsop(vsop.Earth, ltime);
        } else {
            if (!earth) {
                // No aberration, so calculate Earth's position once, at the time of observation.
                earth = CalcVsop(vsop.Earth, time);
            }
        }

        geo = new Vector(h.x-earth.x, h.y-earth.y, h.z-earth.z, time);
        let ltime2 = time.AddDays(-geo.Length() / C_AUDAY);
        dt = Math.abs(ltime2.tt - ltime.tt);
        if (dt < 1.0e-9)
            return geo;
        ltime = ltime2;
    }
    throw `Light-travel time solver did not converge: dt=${dt}`;
}

function ExportState(terse: body_state_t, time: AstroTime): StateVector {
    return new StateVector(terse.r.x, terse.r.y, terse.r.z, terse.v.x, terse.v.y, terse.v.z, time);
}

/**
 * @brief  Calculates barycentric position and velocity vectors for the given body.
 *
 * Given a body and a time, calculates the barycentric position and velocity
 * vectors for the center of that body at that time.
 * The vectors are expressed in equatorial J2000 coordinates (EQJ).
 *
 * @param {Body} body
 *      The celestial body whose barycentric state vector is to be calculated.
 *      Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets:
 *      `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
 *      `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate position and velocity.
 *
 *  @returns {StateVector}
 *      An object that contains barycentric position and velocity vectors.
 */
export function BaryState(body: Body, date: FlexibleDateTime): StateVector {
    const time = MakeTime(date);
    if (body === Body.SSB) {
        // Trivial case: the solar system barycenter itself.
        return new StateVector(0, 0, 0, 0, 0, 0, time);
    }

    if (body === Body.Pluto) {
        return CalcPluto(time, false);
    }

    // Find the barycentric positions and velocities for the 5 major bodies:
    // Sun, Jupiter, Saturn, Uranus, Neptune.
    const bary = new major_bodies_t(time.tt);

    switch (body) {
        case Body.Sun:      return ExportState(bary.Sun, time);
        case Body.Jupiter:  return ExportState(bary.Jupiter, time);
        case Body.Saturn:   return ExportState(bary.Saturn, time);
        case Body.Uranus:   return ExportState(bary.Uranus, time);
        case Body.Neptune:  return ExportState(bary.Neptune, time);

        case Body.Moon:
        case Body.EMB:
            const earth = CalcVsopPosVel(vsop[Body.Earth], time.tt);
            const state = (body === Body.Moon) ? GeoMoonState(time) : GeoEmbState(time);
            return new StateVector(
                state.x  + bary.Sun.r.x + earth.r.x,
                state.y  + bary.Sun.r.y + earth.r.y,
                state.z  + bary.Sun.r.z + earth.r.z,
                state.vx + bary.Sun.v.x + earth.v.x,
                state.vy + bary.Sun.v.y + earth.v.y,
                state.vz + bary.Sun.v.z + earth.v.z,
                time
            );
    }

    // Handle the remaining VSOP bodies: Mercury, Venus, Earth, Mars.
    if (body in vsop) {
        const planet = CalcVsopPosVel(vsop[body], time.tt);
        return new StateVector(
            bary.Sun.r.x + planet.r.x,
            bary.Sun.r.y + planet.r.y,
            bary.Sun.r.z + planet.r.z,
            bary.Sun.v.x + planet.v.x,
            bary.Sun.v.y + planet.v.y,
            bary.Sun.v.z + planet.v.z,
            time
        );
    }

    throw `BaryState: Unsupported body "${body}"`;
}


/**
 * @brief  Calculates heliocentric position and velocity vectors for the given body.
 *
 * Given a body and a time, calculates the position and velocity
 * vectors for the center of that body at that time, relative to the center of the Sun.
 * The vectors are expressed in equatorial J2000 coordinates (EQJ).
 * If you need the position vector only, it is more efficient to call {@link HelioVector}.
 * The Sun's center is a non-inertial frame of reference. In other words, the Sun
 * experiences acceleration due to gravitational forces, mostly from the larger
 * planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
 * kinetic energy, or other quantities that require a non-accelerating frame
 * of reference, consider using {@link BaryState} instead.
 *
 * @param {Body} body
 *      The celestial body whose heliocentric state vector is to be calculated.
 *      Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets:
 *      `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
 *      `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
 *
 *  @param {FlexibleDateTime} date
 *      The date and time for which to calculate position and velocity.
 *
 *  @returns {StateVector}
 *      An object that contains heliocentric position and velocity vectors.
 */
export function HelioState(body: Body, date: FlexibleDateTime): StateVector {
    const time = MakeTime(date);
    switch (body) {
        case Body.Sun:
            // Trivial case: the Sun is the origin of the heliocentric frame.
            return new StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time);

        case Body.SSB:
            // Calculate the barycentric Sun. Then the negative of that is the heliocentric SSB.
            const bary = new major_bodies_t(time.tt);
            return new StateVector(
                -bary.Sun.r.x,
                -bary.Sun.r.y,
                -bary.Sun.r.z,
                -bary.Sun.v.x,
                -bary.Sun.v.y,
                -bary.Sun.v.z,
                time
            );

        case Body.Mercury:
        case Body.Venus:
        case Body.Earth:
        case Body.Mars:
        case Body.Jupiter:
        case Body.Saturn:
        case Body.Uranus:
        case Body.Neptune:
            // Planets included in the VSOP87 model.
            const planet = CalcVsopPosVel(vsop[body], time.tt);
            return ExportState(planet, time);

        case Body.Pluto:
            return CalcPluto(time, true);

        case Body.Moon:
        case Body.EMB:
            const earth = CalcVsopPosVel(vsop.Earth, time.tt);
            const state = (body == Body.Moon) ? GeoMoonState(time) : GeoEmbState(time);
            return new StateVector(
                state.x  + earth.r.x,
                state.y  + earth.r.y,
                state.z  + earth.r.z,
                state.vx + earth.v.x,
                state.vy + earth.v.y,
                state.vz + earth.v.z,
                time
            );

        default:
            throw `HelioState: Unsupported body "${body}"`;
    }
}


function QuadInterp(tm: number, dt: number, fa: number, fm: number, fb: number) {
    let Q = (fb + fa)/2 - fm;
    let R = (fb - fa)/2;
    let S = fm;
    let x: number;

    if (Q == 0) {
        // This is a line, not a parabola.
        if (R == 0) {
            // This is a HORIZONTAL line... can't make progress!
            return null;
        }
        x = -S / R;
        if (x < -1 || x > +1) return null;  // out of bounds
    } else {
        // It really is a parabola. Find roots x1, x2.
        let u = R*R - 4*Q*S;
        if (u <= 0) return null;
        let ru = Math.sqrt(u);
        let x1 = (-R + ru) / (2 * Q);
        let x2 = (-R - ru) / (2 * Q);

        if (-1 <= x1 && x1 <= +1) {
            if (-1 <= x2 && x2 <= +1) return null;
            x = x1;
        } else if (-1 <= x2 && x2 <= +1) {
            x = x2;
        } else {
            return null;
        }
    }

    let t = tm + x*dt;
    let df_dt = (2*Q*x + R) / dt;
    return { x:x, t:t, df_dt:df_dt };
}

// Quirk: for some reason, I need to put the interface declaration *before* its
// documentation, or jsdoc2md will strip it out.

export interface SearchOptions {
    dt_tolerance_seconds? : number;
    init_f1? : number;
    init_f2? : number;
    iter_limit? : number;
}

/**
 * @brief Options for the {@link Search} function.
 *
 * @typedef {object} SearchOptions
 *
 * @property {number | undefined} dt_tolerance_seconds
 *      The number of seconds for a time window smaller than which the search
 *      is considered successful.  Using too large a tolerance can result in
 *      an inaccurate time estimate.  Using too small a tolerance can cause
 *      excessive computation, or can even cause the search to fail because of
 *      limited floating-point resolution.  Defaults to 1 second.
 *
 * @property {number | undefined} init_f1
 *      As an optimization, if the caller of {@link Search}
 *      has already calculated the value of the function being searched (the parameter `func`)
 *      at the time coordinate `t1`, it can pass in that value as `init_f1`.
 *      For very expensive calculations, this can measurably improve performance.
 *
 * @property {number | undefined} init_f2
 *      The same as `init_f1`, except this is the optional initial value of `func(t2)`
 *      instead of `func(t1)`.
 *
 * @property {number | undefined} iter_limit
 */

/**
 * @brief Finds the time when a function ascends through zero.
 *
 * Search for next time <i>t</i> (such that <i>t</i> is between `t1` and `t2`)
 * that `func(t)` crosses from a negative value to a non-negative value.
 * The given function must have "smooth" behavior over the entire inclusive range [`t1`, `t2`],
 * meaning that it behaves like a continuous differentiable function.
 * It is not required that `t1` &lt; `t2`; `t1` &gt; `t2`
 * allows searching backward in time.
 * Note: `t1` and `t2` must be chosen such that there is no possibility
 * of more than one zero-crossing (ascending or descending), or it is possible
 * that the "wrong" event will be found (i.e. not the first event after t1)
 * or even that the function will return `null`, indicating that no event was found.
 *
 * @param {function(AstroTime): number} func
 *      The function to find an ascending zero crossing for.
 *      The function must accept a single parameter of type {@link AstroTime}
 *      and return a numeric value.
 *
 * @param {AstroTime} t1
 *      The lower time bound of a search window.
 *
 * @param {AstroTime} t2
 *      The upper time bound of a search window.
 *
 * @param {SearchOptions | undefined} options
 *      Options that can tune the behavior of the search.
 *      Most callers can omit this argument.
 *
 * @returns {AstroTime | null}
 *      If the search is successful, returns the date and time of the solution.
 *      If the search fails, returns `null`.
 */
export function Search(
    f: (t: AstroTime) => number,
    t1: AstroTime,
    t2: AstroTime,
    options?: SearchOptions
): AstroTime | null {

    const dt_tolerance_seconds = VerifyNumber((options && options.dt_tolerance_seconds) || 1);

    const dt_days = Math.abs(dt_tolerance_seconds / SECONDS_PER_DAY);

    let f1 = (options && options.init_f1) || f(t1);
    let f2 = (options && options.init_f2) || f(t2);
    let fmid: number = NaN;

    let iter = 0;
    let iter_limit = (options && options.iter_limit) || 20;
    let calc_fmid = true;
    while (true) {
        if (++iter > iter_limit)
            throw `Excessive iteration in Search()`;

        let tmid = InterpolateTime(t1, t2, 0.5);
        let dt = tmid.ut - t1.ut;

        if (Math.abs(dt) < dt_days) {
            // We are close enough to the event to stop the search.
            return tmid;
        }

        if (calc_fmid)
            fmid = f(tmid);
        else
            calc_fmid = true;   // we already have the correct value of fmid from the previous loop

        // Quadratic interpolation:
        // Try to find a parabola that passes through the 3 points we have sampled:
        // (t1,f1), (tmid,fmid), (t2,f2).
        let q = QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2);

        // Did we find an approximate root-crossing?
        if (q) {
            // Evaluate the function at our candidate solution.
            let tq = MakeTime(q.t);
            let fq = f(tq);

            if (q.df_dt !== 0) {
                if (Math.abs(fq / q.df_dt) < dt_days) {
                    // The estimated time error is small enough that we can quit now.
                    return tq;
                }

                // Try guessing a tighter boundary with the interpolated root at the center.
                let dt_guess = 1.2 * Math.abs(fq / q.df_dt);
                if (dt_guess < dt/10) {
                    let tleft = tq.AddDays(-dt_guess);
                    let tright = tq.AddDays(+dt_guess);
                    if ((tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0) {
                        if ((tright.ut - t1.ut)*(tright.ut - t2.ut) < 0) {
                            let fleft = f(tleft);
                            let fright = f(tright);
                            if (fleft<0 && fright>=0) {
                                f1 = fleft;
                                f2 = fright;
                                t1 = tleft;
                                t2 = tright;
                                fmid = fq;
                                calc_fmid = false;
                                continue;
                            }
                        }
                    }
                }
            }
        }

        if (f1<0 && fmid>=0) {
            t2 = tmid;
            f2 = fmid;
            continue;
        }

        if (fmid<0 && f2>=0) {
            t1 = tmid;
            f1 = fmid;
            continue;
        }

        // Either there is no ascending zero-crossing in this range
        // or the search window is too wide.
        return null;
    }
}

function LongitudeOffset(diff: number): number {
    let offset = diff;
    while (offset <= -180) offset += 360;
    while (offset > 180) offset -= 360;
    return offset;
}

function NormalizeLongitude(lon: number): number {
    while (lon < 0) lon += 360;
    while (lon >= 360) lon -= 360;
    return lon;
}

/**
 * @brief Searches for when the Sun reaches a given ecliptic longitude.
 *
 * Searches for the moment in time when the center of the Sun reaches a given apparent
 * ecliptic longitude, as seen from the center of the Earth, within a given range of dates.
 * This function can be used to determine equinoxes and solstices.
 * However, it is usually more convenient and efficient to call {@link Seasons}
 * to calculate equinoxes and solstices for a given calendar year.
 * `SearchSunLongitude` is more general in that it allows searching for arbitrary longitude values.
 *
 * @param {number} targetLon
 *      The desired ecliptic longitude of date in degrees.
 *      This may be any value in the range [0, 360), although certain
 *      values have conventional meanings:
 *
 *      When `targetLon` is 0, finds the March equinox,
 *      which is the moment spring begins in the northern hemisphere
 *      and the beginning of autumn in the southern hemisphere.
 *
 *      When `targetLon` is 180, finds the September equinox,
 *      which is the moment autumn begins in the northern hemisphere and
 *      spring begins in the southern hemisphere.
 *
 *      When `targetLon` is 90, finds the northern solstice, which is the
 *      moment summer begins in the northern hemisphere and winter
 *      begins in the southern hemisphere.
 *
 *      When `targetLon` is 270, finds the southern solstice, which is the
 *      moment winter begins in the northern hemisphere and summer
 *      begins in the southern hemisphere.
 *
 * @param {FlexibleDateTime} dateStart
 *      A date and time known to be earlier than the desired longitude event.
 *
 * @param {number} limitDays
 *      A floating point number of days, which when added to `dateStart`,
 *      yields a date and time known to be after the desired longitude event.
 *
 * @returns {AstroTime | null}
 *      The date and time when the Sun reaches the apparent ecliptic longitude `targetLon`
 *      within the range of times specified by `dateStart` and `limitDays`.
 *      If the Sun does not reach the target longitude within the specified time range, or the
 *      time range is excessively wide, the return value is `null`.
 *      To avoid a `null` return value, the caller must pick a time window around
 *      the event that is within a few days but not so small that the event might fall outside the window.
 */
export function SearchSunLongitude(targetLon: number, dateStart: FlexibleDateTime, limitDays: number): AstroTime | null {
    function sun_offset(t: AstroTime): number {
        let pos = SunPosition(t);
        return LongitudeOffset(pos.elon - targetLon);
    }
    VerifyNumber(targetLon);
    VerifyNumber(limitDays);
    let t1 = MakeTime(dateStart);
    let t2 = t1.AddDays(limitDays);
    return Search(sun_offset, t1, t2, {dt_tolerance_seconds: 0.01});
}

/**
 * @brief Returns one body's ecliptic longitude with respect to another, as seen from the Earth.
 *
 * This function determines where one body appears around the ecliptic plane
 * (the plane of the Earth's orbit around the Sun) as seen from the Earth,
 * relative to the another body's apparent position.
 * The function returns an angle in the half-open range [0, 360) degrees.
 * The value is the ecliptic longitude of `body1` relative to the ecliptic
 * longitude of `body2`.
 *
 * The angle is 0 when the two bodies are at the same ecliptic longitude
 * as seen from the Earth. The angle increases in the prograde direction
 * (the direction that the planets orbit the Sun and the Moon orbits the Earth).
 *
 * When the angle is 180 degrees, it means the two bodies appear on opposite sides
 * of the sky for an Earthly observer.
 *
 * Neither `body1` nor `body2` is allowed to be `Body.Earth`.
 * If this happens, the function throws an exception.
 *
 * @param {Body} body1
 *      The first body, whose longitude is to be found relative to the second body.
 *
 * @param {Body} body2
 *      The second body, relative to which the longitude of the first body is to be found.
 *
 * @param {FlexibleDateTime} date
 *      The date and time of the observation.
 *
 * @returns {number}
 *      An angle in the range [0, 360), expressed in degrees.
 */
export function PairLongitude(body1: Body, body2: Body, date: FlexibleDateTime): number {
    if (body1 === Body.Earth || body2 === Body.Earth)
        throw 'The Earth does not have a longitude as seen from itself.';

    const time = MakeTime(date);

    const vector1 = GeoVector(body1, time, false);
    const eclip1 = Ecliptic(vector1);

    const vector2 = GeoVector(body2, time, false);
    const eclip2 = Ecliptic(vector2);

    return NormalizeLongitude(eclip1.elon - eclip2.elon);
}

/**
 * @brief Calculates the angular separation between the Sun and the given body.
 *
 * Returns the full angle seen from
 * the Earth, between the given body and the Sun.
 * Unlike {@link PairLongitude}, this function does not
 * project the body's "shadow" onto the ecliptic;
 * the angle is measured in 3D space around the plane that
 * contains the centers of the Earth, the Sun, and `body`.
 *
 * @param {Body} body
 *      The name of a supported celestial body other than the Earth.
 *
 * @param {FlexibleDateTime} date
 *      The time at which the angle from the Sun is to be found.
 *
 * @returns {number}
 *      An angle in degrees in the range [0, 180].
 */
export function AngleFromSun(body: Body, date: FlexibleDateTime): number {
    if (body == Body.Earth)
        throw 'The Earth does not have an angle as seen from itself.';

    const time = MakeTime(date);
    const sv = GeoVector(Body.Sun, time, true);
    const bv = GeoVector(body, time, true);
    const angle = AngleBetween(sv, bv);
    return angle;
}

/**
 * @brief Calculates heliocentric ecliptic longitude based on the J2000 equinox.
 *
 * @param {Body} body
 *      The name of a celestial body other than the Sun.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the ecliptic longitude.
 *
 * @returns {number}
 *      The ecliptic longitude angle of the body in degrees measured counterclockwise around the mean
 *      plane of the Earth's orbit, as seen from above the Sun's north pole.
 *      Ecliptic longitude starts at 0 at the J2000
 *      <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a> and
 *      increases in the same direction the Earth orbits the Sun.
 *      The returned value is always in the range [0, 360).
 */
export function EclipticLongitude(body: Body, date: FlexibleDateTime): number {
    if (body === Body.Sun)
        throw 'Cannot calculate heliocentric longitude of the Sun.';

    const hv = HelioVector(body, date);
    const eclip = Ecliptic(hv);
    return eclip.elon;
}

function VisualMagnitude(body: Body, phase: number, helio_dist: number, geo_dist: number): number {
    // For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212
    let c0: number, c1=0, c2=0, c3=0;
    switch (body) {
    case Body.Mercury:     c0 = -0.60; c1 = +4.98; c2 = -4.88; c3 = +3.02; break;
    case Body.Venus:
        if (phase < 163.6) {
            c0 = -4.47; c1 = +1.03; c2 = +0.57; c3 = +0.13;
        } else {
            c0 = 0.98; c1 = -1.02;
        }
        break;
    case Body.Mars:        c0 = -1.52; c1 = +1.60;                         break;
    case Body.Jupiter:     c0 = -9.40; c1 = +0.50;                         break;
    case Body.Uranus:      c0 = -7.19; c1 = +0.25;                         break;
    case Body.Neptune:     c0 = -6.87;                                     break;
    case Body.Pluto:       c0 = -1.00; c1 = +4.00;                         break;
    default: throw `VisualMagnitude: unsupported body ${body}`;
    }

    const x = phase / 100;
    let mag = c0 + x*(c1 + x*(c2 + x*c3));
    mag += 5*Math.log10(helio_dist * geo_dist);
    return mag;
}

function SaturnMagnitude(phase: number, helio_dist: number, geo_dist: number, gc: Vector, time: AstroTime) {
    // Based on formulas by Paul Schlyter found here:
    // http://www.stjarnhimlen.se/comp/ppcomp.html#15

    // We must handle Saturn's rings as a major component of its visual magnitude.
    // Find geocentric ecliptic coordinates of Saturn.
    const eclip = Ecliptic(gc);
    const ir = DEG2RAD * 28.06;   // tilt of Saturn's rings to the ecliptic, in radians
    const Nr = DEG2RAD * (169.51 + (3.82e-5 * time.tt));    // ascending node of Saturn's rings, in radians

    // Find tilt of Saturn's rings, as seen from Earth.
    const lat = DEG2RAD * eclip.elat;
    const lon = DEG2RAD * eclip.elon;
    const tilt = Math.asin(Math.sin(lat)*Math.cos(ir) - Math.cos(lat)*Math.sin(ir)*Math.sin(lon-Nr));
    const sin_tilt = Math.sin(Math.abs(tilt));

    let mag = -9.0 + 0.044*phase;
    mag += sin_tilt*(-2.6 + 1.2*sin_tilt);
    mag += 5*Math.log10(helio_dist * geo_dist);
    return { mag:mag, ring_tilt:RAD2DEG*tilt };
}

function MoonMagnitude(phase: number, helio_dist: number, geo_dist: number): number {
    // https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
    let rad = phase * DEG2RAD;
    let rad2 = rad * rad;
    let rad4 = rad2 * rad2;
    let mag = -12.717 + 1.49*Math.abs(rad) + 0.0431*rad4;

    const moon_mean_distance_au = 385000.6 / KM_PER_AU;
    let geo_au = geo_dist / moon_mean_distance_au;
    mag += 5*Math.log10(helio_dist * geo_au);
    return mag;
}

/**
 * @brief Information about the apparent brightness and sunlit phase of a celestial object.
 *
 * @property {AstroTime} time
 *      The date and time pertaining to the other calculated values in this object.
 *
 * @property {number} mag
 *      The <a href="https://en.wikipedia.org/wiki/Apparent_magnitude">apparent visual magnitude</a> of the celestial body.
 *
 * @property {number} phase_angle
 *      The angle in degrees as seen from the center of the celestial body between the Sun and the Earth.
 *      The value is always in the range 0 to 180.
 *      The phase angle provides a measure of what fraction of the body's face appears
 *      illuminated by the Sun as seen from the Earth.
 *      When the observed body is the Sun, the `phase` property is set to 0,
 *      although this has no physical meaning because the Sun emits, rather than reflects, light.
 *      When the phase is near 0 degrees, the body appears "full".
 *      When it is 90 degrees, the body appears "half full".
 *      And when it is 180 degrees, the body appears "new" and is very difficult to see
 *      because it is both dim and lost in the Sun's glare as seen from the Earth.
 *
 * @property {number} phase_fraction
 *      The fraction of the body's face that is illuminated by the Sun, as seen from the Earth.
 *      Calculated from `phase_angle` for convenience.
 *      This value ranges from 0 to 1.
 *
 * @property {number} helio_dist
 *      The distance between the center of the Sun and the center of the body in
 *      <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU).
 *
 * @property {number} geo_dist
 *      The distance between the center of the Earth and the center of the body in AU.
 *
 * @property {Vector} gc
 *      Geocentric coordinates: the 3D vector from the center of the Earth to the center of the body.
 *      The components are in expressed in AU and are oriented with respect to the J2000 equatorial plane.
 *
 * @property {Vector} hc
 *      Heliocentric coordinates: The 3D vector from the center of the Sun to the center of the body.
 *      Like `gc`, `hc` is expressed in AU and oriented with respect
 *      to the J2000 equatorial plane.
 *
 * @property {number | undefined} ring_tilt
 *      For Saturn, this is the angular tilt of the planet's rings in degrees away
 *      from the line of sight from the Earth. When the value is near 0, the rings
 *      appear edge-on from the Earth and are therefore difficult to see.
 *      When `ring_tilt` approaches its maximum value (about 27 degrees),
 *      the rings appear widest and brightest from the Earth.
 *      Unlike the <a href="https://ssd.jpl.nasa.gov/horizons.cgi">JPL Horizons</a> online tool,
 *      this library includes the effect of the ring tilt angle in the calculated value
 *      for Saturn's visual magnitude.
 *      For all bodies other than Saturn, the value of `ring_tilt` is `undefined`.
 */
export class IlluminationInfo {
    phase_fraction: number;

    constructor(
        public time: AstroTime,
        public mag: number,
        public phase_angle: number,
        public helio_dist: number,
        public geo_dist: number,
        public gc: Vector,
        public hc: Vector,
        public ring_tilt?: number)
    {
        this.phase_fraction = (1 + Math.cos(DEG2RAD * phase_angle)) / 2;
    }
}

/**
 * @brief Calculates visual magnitude and related information about a body.
 *
 * Calculates the phase angle, visual magnitude,
 * and other values relating to the body's illumination
 * at the given date and time, as seen from the Earth.
 *
 * @param {Body} body
 *      The name of the celestial body being observed.
 *      Not allowed to be `Body.Earth`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the illumination data for the given body.
 *
 * @returns {IlluminationInfo}
 */
export function Illumination(body: Body, date: FlexibleDateTime): IlluminationInfo {
    if (body === Body.Earth)
        throw `The illumination of the Earth is not defined.`;

    const time = MakeTime(date);
    const earth = CalcVsop(vsop.Earth, time);
    let phase: number;      // phase angle in degrees between Earth and Sun as seen from body
    let hc: Vector;         // vector from Sun to body
    let gc: Vector;         // vector from Earth to body
    let mag: number;        // visual magnitude

    if (body === Body.Sun) {
        gc = new Vector(-earth.x, -earth.y, -earth.z, time);
        hc = new Vector(0, 0, 0, time);
        phase = 0;      // a placeholder value; the Sun does not have an illumination phase because it emits, rather than reflects, light.
    } else {
        if (body === Body.Moon) {
            // For extra numeric precision, use geocentric moon formula directly.
            gc = GeoMoon(time);
            hc = new Vector(earth.x + gc.x, earth.y + gc.y, earth.z + gc.z, time);
        } else {
            // For planets, heliocentric vector is most direct to calculate.
            hc = HelioVector(body, date);
            gc = new Vector(hc.x - earth.x, hc.y - earth.y, hc.z - earth.z, time);
        }
        phase = AngleBetween(gc, hc);
    }

    let geo_dist = gc.Length();     // distance from body to center of Earth
    let helio_dist = hc.Length();   // distance from body to center of Sun
    let ring_tilt;   // only reported for Saturn

    if (body === Body.Sun) {
        mag = SUN_MAG_1AU + 5*Math.log10(geo_dist);
    } else if (body === Body.Moon) {
        mag = MoonMagnitude(phase, helio_dist, geo_dist);
    } else if (body === Body.Saturn) {
        const saturn = SaturnMagnitude(phase, helio_dist, geo_dist, gc, time);
        mag = saturn.mag;
        ring_tilt = saturn.ring_tilt;
    } else {
        mag = VisualMagnitude(body, phase, helio_dist, geo_dist);
    }

    return new IlluminationInfo(time, mag, phase, helio_dist, geo_dist, gc, hc, ring_tilt);
}

function SynodicPeriod(body: Body): number {
    if (body === Body.Earth)
        throw 'The Earth does not have a synodic period as seen from itself.';

    if (body === Body.Moon)
        return MEAN_SYNODIC_MONTH;

    // Calculate the synodic period of the planet from its and the Earth's sidereal periods.
    // The sidereal period of a planet is how long it takes to go around the Sun in days, on average.
    // The synodic period of a planet is how long it takes between consecutive oppositions
    // or conjunctions, on average.

    let planet = Planet[body];
    if (!planet)
        throw `Not a valid planet name: ${body}`;

    // See here for explanation of the formula:
    // https://en.wikipedia.org/wiki/Elongation_(astronomy)#Elongation_period

    const Te = Planet.Earth.OrbitalPeriod;
    const Tp = planet.OrbitalPeriod;
    const synodicPeriod = Math.abs(Te / (Te/Tp - 1));

    return synodicPeriod;
}

/**
 * @brief Searches for when the Earth and a given body reach a relative ecliptic longitude separation.
 *
 * Searches for the date and time the relative ecliptic longitudes of
 * the specified body and the Earth, as seen from the Sun, reach a certain
 * difference. This function is useful for finding conjunctions and oppositions
 * of the planets. For the opposition of a superior planet (Mars, Jupiter, ..., Pluto),
 * or the inferior conjunction of an inferior planet (Mercury, Venus),
 * call with `targetRelLon` = 0. The 0 value indicates that both
 * planets are on the same ecliptic longitude line, ignoring the other planet's
 * distance above or below the plane of the Earth's orbit.
 * For superior conjunctions, call with `targetRelLon` = 180.
 * This means the Earth and the other planet are on opposite sides of the Sun.
 *
 * @param {Body} body
 *      Any planet other than the Earth.
 *
 * @param {number} targetRelLon
 *      The desired angular difference in degrees between the ecliptic longitudes
 *      of `body` and the Earth. Must be in the range (-180, +180].
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to find the next occurrence of the
 *      body and the Earth reaching the desired relative longitude.
 *
 * @returns {AstroTime}
 *      The time when the Earth and the body next reach the specified relative longitudes.
 */
export function SearchRelativeLongitude(body: Body, targetRelLon: number, startDate: FlexibleDateTime): AstroTime {
    VerifyNumber(targetRelLon);
    const planet = Planet[body];
    if (!planet)
        throw `Cannot search relative longitude because body is not a planet: ${body}`;

    if (body === Body.Earth)
        throw 'Cannot search relative longitude for the Earth (it is always 0)';

    // Determine whether the Earth "gains" (+1) on the planet or "loses" (-1)
    // as both race around the Sun.
    const direction = (planet.OrbitalPeriod > Planet.Earth.OrbitalPeriod) ? +1 : -1;

    function offset(t: AstroTime): number {
        const plon = EclipticLongitude(body, t);
        const elon = EclipticLongitude(Body.Earth, t);
        const diff = direction * (elon - plon);
        return LongitudeOffset(diff - targetRelLon);
    }

    let syn = SynodicPeriod(body);
    let time = MakeTime(startDate);

    // Iterate until we converge on the desired event.
    // Calculate the error angle, which will be a negative number of degrees,
    // meaning we are "behind" the target relative longitude.
    let error_angle = offset(time);
    if (error_angle > 0) error_angle -= 360;    // force searching forward in time

    for (let iter=0; iter < 100; ++iter) {
        // Estimate how many days in the future (positive) or past (negative)
        // we have to go to get closer to the target relative longitude.
        let day_adjust = (-error_angle/360) * syn;
        time = time.AddDays(day_adjust);
        if (Math.abs(day_adjust) * SECONDS_PER_DAY < 1)
            return time;

        let prev_angle = error_angle;
        error_angle = offset(time);

        if (Math.abs(prev_angle) < 30) {
            // Improve convergence for Mercury/Mars (eccentric orbits)
            // by adjusting the synodic period to more closely match the
            // variable speed of both planets in this part of their respective orbits.
            if (prev_angle !== error_angle) {
                let ratio = prev_angle / (prev_angle - error_angle);
                if (ratio > 0.5 && ratio < 2.0)
                    syn *= ratio;
            }
        }
    }

    throw `Relative longitude search failed to converge for ${body} near ${time.toString()} (error_angle = ${error_angle}).`;
}

/**
 * @brief Determines the moon's phase expressed as an ecliptic longitude.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the moon's phase.
 *
 * @returns {number}
 *      A value in the range [0, 360) indicating the difference
 *      in ecliptic longitude between the center of the Sun and the
 *      center of the Moon, as seen from the center of the Earth.
 *      Certain longitude values have conventional meanings:
 *
 * * 0 = new moon
 * * 90 = first quarter
 * * 180 = full moon
 * * 270 = third quarter
 */
export function MoonPhase(date: FlexibleDateTime): number {
    return PairLongitude(Body.Moon, Body.Sun, date);
}

/**
 * @brief Searches for the date and time that the Moon reaches a specified phase.
 *
 * Lunar phases are defined in terms of geocentric ecliptic longitudes
 * with respect to the Sun.  When the Moon and the Sun have the same ecliptic
 * longitude, that is defined as a new moon. When the two ecliptic longitudes
 * are 180 degrees apart, that is defined as a full moon.
 * To enumerate quarter lunar phases, it is simpler to call
 * {@link SearchMoonQuarter} once, followed by repeatedly calling
 * {@link NextMoonQuarter}. `SearchMoonPhase` is only
 * necessary for finding other lunar phases than the usual quarter phases.
 *
 * @param {number} targetLon
 *      The difference in geocentric ecliptic longitude between the Sun and Moon
 *      that specifies the lunar phase being sought. This can be any value
 *      in the range [0, 360). Here are some helpful examples:
 *      0 = new moon,
 *      90 = first quarter,
 *      180 = full moon,
 *      270 = third quarter.
 *
 * @param {FlexibleDateTime} dateStart
 *      The beginning of the window of time in which to search.
 *
 * @param {number} limitDays
 *      The floating point number of days after `dateStart`
 *      that limits the window of time in which to search.
 *
 * @returns {AstroTime | null}
 *      If the specified lunar phase occurs after `dateStart`
 *      and before `limitDays` days after `dateStart`,
 *      this function returns the date and time of the first such occurrence.
 *      Otherwise, it returns `null`.
 */
export function SearchMoonPhase(targetLon: number, dateStart: FlexibleDateTime, limitDays: number): AstroTime | null {
    function moon_offset(t: AstroTime): number {
        let mlon: number = MoonPhase(t);
        return LongitudeOffset(mlon - targetLon);
    }

    VerifyNumber(targetLon);
    VerifyNumber(limitDays);

    // To avoid discontinuities in the moon_offset function causing problems,
    // we need to approximate when that function will next return 0.
    // We probe it with the start time and take advantage of the fact
    // that every lunar phase repeats roughly every 29.5 days.
    // There is a surprising uncertainty in the quarter timing,
    // due to the eccentricity of the moon's orbit.
    // I have seen more than 0.9 days away from the simple prediction.
    // To be safe, we take the predicted time of the event and search
    // +/-1.5 days around it (a 3.0-day wide window).
    // But we must return null if the final result goes beyond limitDays after dateStart.
    const uncertainty = 1.5;

    let ta = MakeTime(dateStart);
    let ya = moon_offset(ta);
    if (ya > 0) ya -= 360;  // force searching forward in time, not backward
    let est_dt = -(MEAN_SYNODIC_MONTH*ya)/360;
    let dt1 = est_dt - uncertainty;
    if (dt1 > limitDays)
        return null;   // not possible for moon phase to occur within the specified window
    let dt2 = Math.min(limitDays, est_dt + uncertainty);
    let t1 = ta.AddDays(dt1);
    let t2 = ta.AddDays(dt2);
    return Search(moon_offset, t1, t2);
}

/**
 * @brief A quarter lunar phase, along with when it occurs.
 *
 * @property {number} quarter
 *      An integer as follows:
 *      0 = new moon,
 *      1 = first quarter,
 *      2 = full moon,
 *      3 = third quarter.
 *
 * @property {AstroTime} time
 *      The date and time of the quarter lunar phase.
 */
export class MoonQuarter {
    constructor(
        public quarter: number,
        public time: AstroTime)
        {}
}

/**
 * @brief Finds the first quarter lunar phase after the specified date and time.
 *
 * The quarter lunar phases are: new moon, first quarter, full moon, and third quarter.
 * To enumerate quarter lunar phases, call `SearchMoonQuarter` once,
 * then pass its return value to {@link NextMoonQuarter} to find the next
 * `MoonQuarter`. Keep calling `NextMoonQuarter` in a loop,
 * passing the previous return value as the argument to the next call.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which to find the first quarter lunar phase.
 *
 * @returns {MoonQuarter}
 */
export function SearchMoonQuarter(dateStart: FlexibleDateTime): MoonQuarter {
    // Determine what the next quarter phase will be.
    let phaseStart = MoonPhase(dateStart);
    let quarterStart = Math.floor(phaseStart / 90);
    let quarter = (quarterStart + 1) % 4;
    let time = SearchMoonPhase(90 * quarter, dateStart, 10);
    if (!time)
        throw 'Cannot find moon quarter';
    return new MoonQuarter(quarter, time);
}

/**
 * @brief Finds the next quarter lunar phase in a series.
 *
 * Given a {@link MoonQuarter} object, finds the next consecutive
 * quarter lunar phase. See remarks in {@link SearchMoonQuarter}
 * for explanation of usage.
 *
 * @param {MoonQuarter} mq
 *      The return value of a prior call to {@link MoonQuarter} or `NextMoonQuarter`.
 *
 * @returns {MoonQuarter}
 */
export function NextMoonQuarter(mq: MoonQuarter): MoonQuarter {
    // Skip 6 days past the previous found moon quarter to find the next one.
    // This is less than the minimum possible increment.
    // So far I have seen the interval well contained by the range (6.5, 8.3) days.
    let date = new Date(mq.time.date.getTime() + 6*MILLIS_PER_DAY);
    return SearchMoonQuarter(date);
}

function BodyRadiusAu(body: Body): number {
    // For the purposes of calculating rise/set times,
    // only the Sun and Moon appear large enough to an observer
    // on the Earth for their radius to matter.
    // All other bodies are treated as points.
    switch (body) {
        case Body.Sun:   return SUN_RADIUS_AU;
        case Body.Moon:  return MOON_EQUATORIAL_RADIUS_AU;
        default:      return 0;
    }
}

/**
 * @brief Finds the next rise or set time for a body.
 *
 * Finds a rise or set time for the given body as
 * seen by an observer at the specified location on the Earth.
 * Rise time is defined as the moment when the top of the body
 * is observed to first appear above the horizon in the east.
 * Set time is defined as the moment the top of the body
 * is observed to sink below the horizon in the west.
 * The times are adjusted for typical atmospheric refraction conditions.
 *
 * @param {Body} body
 *      The body to find the rise or set time for.
 *
 * @param {Observer} observer
 *      Specifies the geographic coordinates and elevation above sea level of the observer.
 *
 * @param {number} direction
 *      Either +1 to find rise time or -1 to find set time.
 *      Any other value will cause an exception to be thrown.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which the specified rise or set time is to be found.
 *
 * @param {number} limitDays
 *      The fractional number of days after `dateStart` that limits
 *      when the rise or set time is to be found.
 *
 * @returns {AstroTime | null}
 *      The date and time of the rise or set event, or null if no such event
 *      occurs within the specified time window.
 */
export function SearchRiseSet(
    body: Body,
    observer: Observer,
    direction: number,
    dateStart: FlexibleDateTime,
    limitDays: number): AstroTime | null
{
    let body_radius_au:number = BodyRadiusAu(body);

    function peak_altitude(t: AstroTime): number {
        // Return the angular altitude above or below the horizon
        // of the highest part (the peak) of the given object.
        // This is defined as the apparent altitude of the center of the body plus
        // the body's angular radius.
        // The 'direction' variable in the enclosing function controls
        // whether the angle is measured positive above the horizon or
        // positive below the horizon, depending on whether the caller
        // wants rise times or set times, respectively.

        const ofdate = Equator(body, t, observer, true, true);
        const hor = Horizon(t, observer, ofdate.ra, ofdate.dec);
        const alt = hor.altitude + RAD2DEG*(body_radius_au / ofdate.dist) + REFRACTION_NEAR_HORIZON;
        return direction * alt;
    }

    return InternalSearchAltitude(body, observer, direction, dateStart, limitDays, peak_altitude);
}

/**
 * @brief Finds the next time a body reaches a given altitude.
 *
 * Finds when the given body ascends or descends through a given
 * altitude angle, as seen by an observer at the specified location on the Earth.
 * By using the appropriate combination of `direction` and `altitude` parameters,
 * this function can be used to find when civil, nautical, or astronomical twilight
 * begins (dawn) or ends (dusk).
 *
 * Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
 * the horizon. To find civil dawn, pass +1 for `direction` and -6 for `altitude`.
 *
 * Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
 * To find civil dusk, pass -1 for `direction` and -6 for `altitude`.
 *
 * Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.
 *
 * Astronomical twilight uses -18 degrees as the `altitude` value.
 *
 * @param {Body} body
 *      The body for which to find the altitude event.
 *      Can be the Sun, Moon, or any planet other than the Earth.
 *
 * @param {Observer} observer
 *      Specifies the geographic coordinates and elevation above sea level of the observer.
 *
 * @param {number} direction
 *      Either +1 to find when the body ascends through the altitude,
 *      or -1 for when the body descends through the altitude.
 *      Any other value will cause an exception to be thrown.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which the specified altitude event is to be found.
 *
 * @param {number} limitDays
 *      The fractional number of days after `dateStart` that limits
 *      when the altitude event is to be found. Must be a positive number.
 *
 * @param {number} altitude
 *      The desired altitude angle of the body's center above (positive)
 *      or below (negative) the observer's local horizon, expressed in degrees.
 *      Must be in the range [-90, +90].
 *
 * @returns {AstroTime | null}
 *      The date and time of the altitude event, or null if no such event
 *      occurs within the specified time window.
 */
 export function SearchAltitude(
    body: Body,
    observer: Observer,
    direction: number,
    dateStart: FlexibleDateTime,
    limitDays: number,
    altitude: number): AstroTime | null
{
    if (!Number.isFinite(altitude) || altitude < -90 || altitude > +90)
        throw `Invalid altitude angle: ${altitude}`;

    function altitude_error(t: AstroTime): number {
        const ofdate = Equator(body, t, observer, true, true);
        const hor = Horizon(t, observer, ofdate.ra, ofdate.dec);
        return direction * (hor.altitude - altitude);
    }

    return InternalSearchAltitude(body, observer, direction, dateStart, limitDays, altitude_error);
}

function InternalSearchAltitude(
    body: Body,
    observer: Observer,
    direction: number,
    dateStart: FlexibleDateTime,
    limitDays: number,
    altitude_error: (t: AstroTime) => number):    AstroTime | null
{
    VerifyObserver(observer);
    VerifyNumber(limitDays);
    if (limitDays <= 0)
        throw `Invalid value for limitDays: ${limitDays}`;

    if (body === Body.Earth)
        throw 'Cannot find altitude event for the Earth.';

    // See if the body is currently above/below the horizon.
    // If we are looking for next rise time and the body is below the horizon,
    // we use the current time as the lower time bound and the next culmination
    // as the upper bound.
    // If the body is above the horizon, we search for the next bottom and use it
    // as the lower bound and the next culmination after that bottom as the upper bound.
    // The same logic applies for finding set times, only we swap the hour angles.
    // The peak_altitude() function already considers the 'direction' parameter.

    let ha_before: number, ha_after: number;
    if (direction === +1) {
        ha_before = 12;     // minimum altitude (bottom) happens BEFORE the body rises.
        ha_after = 0;       // maximum altitude (culmination) happens AFTER the body rises.
    } else if (direction === -1) {
        ha_before = 0;      // culmination happens BEFORE the body sets.
        ha_after = 12;      // bottom happens AFTER the body sets.
    } else {
        throw `Invalid direction parameter ${direction} -- must be +1 or -1`;
    }

    let time_start = MakeTime(dateStart);
    let time_before: AstroTime;
    let evt_before: HourAngleEvent;
    let evt_after: HourAngleEvent;
    let error_before = altitude_error(time_start);
    let error_after: number;
    if (error_before > 0) {
        // We are past the sought event, so we have to wait for the next "before" event (culm/bottom).
        evt_before = SearchHourAngle(body, observer, ha_before, time_start);
        time_before = evt_before.time;
        error_before = altitude_error(time_before);
    } else {
        // We are before or at the sought event, so we find the next "after" event (bottom/culm),
        // and use the current time as the "before" event.
        time_before = time_start;
    }
    evt_after = SearchHourAngle(body, observer, ha_after, time_before);
    error_after = altitude_error(evt_after.time);

    while (true) {
        if (error_before <= 0 && error_after > 0) {
            // Search between evt_before and evt_after for the desired event.
            let tx = Search(altitude_error, time_before, evt_after.time, {init_f1:error_before, init_f2:error_after});
            if (tx)
                return tx;
        }

        // If we didn't find the desired event, use time_after to find the next before-event.
        evt_before = SearchHourAngle(body, observer, ha_before, evt_after.time);
        evt_after = SearchHourAngle(body, observer, ha_after, evt_before.time);
        if (evt_before.time.ut >= time_start.ut + limitDays)
            return null;

        time_before = evt_before.time;
        error_before = altitude_error(evt_before.time);
        error_after = altitude_error(evt_after.time);
    }
}


/**
 * @brief Horizontal position of a body upon reaching an hour angle.
 *
 * Returns information about an occurrence of a celestial body
 * reaching a given hour angle as seen by an observer at a given
 * location on the surface of the Earth.
 *
 * @property {AstroTime} time
 *      The date and time of the celestial body reaching the hour angle.
 *
 * @property {HorizontalCoordinates} hor
 *      Topocentric horizontal coordinates for the body
 *      at the time indicated by the `time` property.
 */
export class HourAngleEvent {
    constructor(
        public time: AstroTime,
        public hor: HorizontalCoordinates)
        {}
}

/**
 * @brief Finds when a body will reach a given hour angle.
 *
 * Finds the next time the given body is seen to reach the specified
 * <a href="https://en.wikipedia.org/wiki/Hour_angle">hour angle</a>
 * by the given observer.
 * Providing `hourAngle` = 0 finds the next maximum altitude event (culmination).
 * Providing `hourAngle` = 12 finds the next minimum altitude event.
 * Note that, especially close to the Earth's poles, a body as seen on a given day
 * may always be above the horizon or always below the horizon, so the caller cannot
 * assume that a culminating object is visible nor that an object is below the horizon
 * at its minimum altitude.
 *
 * @param {Body} body
 *      A celestial body other than the Earth.
 *
 * @param {Observer} observer
 *      Specifies the geographic coordinates and elevation above sea level of the observer.
 *
 * @param {number} hourAngle
 *      The hour angle expressed in
 *      <a href="https://en.wikipedia.org/wiki/Sidereal_time">sidereal</a>
 *      hours for which the caller seeks to find the body attain.
 *      The value must be in the range [0, 24).
 *      The hour angle represents the number of sidereal hours that have
 *      elapsed since the most recent time the body crossed the observer's local
 *      <a href="https://en.wikipedia.org/wiki/Meridian_(astronomy)">meridian</a>.
 *      This specifying `hourAngle` = 0 finds the moment in time
 *      the body reaches the highest angular altitude in a given sidereal day.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which the desired hour angle crossing event
 *      is to be found.
 *
 * @returns {HourAngleEvent}
 */
export function SearchHourAngle(body: Body, observer: Observer, hourAngle: number, dateStart: FlexibleDateTime): HourAngleEvent {
    VerifyObserver(observer);
    let time = MakeTime(dateStart);
    let iter = 0;

    if (body === Body.Earth)
        throw 'Cannot search for hour angle of the Earth.';

    VerifyNumber(hourAngle);
    if (hourAngle < 0.0 || hourAngle >= 24.0)
        throw `Invalid hour angle ${hourAngle}`;

    while (true) {
        ++iter;

        // Calculate Greenwich Apparent Sidereal Time (GAST) at the given time.
        let gast = sidereal_time(time);

        let ofdate = Equator(body, time, observer, true, true);

        // Calculate the adjustment needed in sidereal time to bring
        // the hour angle to the desired value.
        let delta_sidereal_hours = ((hourAngle + ofdate.ra - observer.longitude/15) - gast) % 24;
        if (iter === 1) {
            // On the first iteration, always search forward in time.
            if (delta_sidereal_hours < 0)
                delta_sidereal_hours += 24;
        } else {
            // On subsequent iterations, we make the smallest possible adjustment,
            // either forward or backward in time.
            if (delta_sidereal_hours < -12)
                delta_sidereal_hours += 24;
            else if (delta_sidereal_hours > +12)
                delta_sidereal_hours -= 24;
        }

        // If the error is tolerable (less than 0.1 seconds), stop searching.
        if (Math.abs(delta_sidereal_hours) * 3600 < 0.1) {
            const hor = Horizon(time, observer, ofdate.ra, ofdate.dec, 'normal');
            return new HourAngleEvent(time, hor);
        }

        // We need to loop another time to get more accuracy.
        // Update the terrestrial time adjusting by sidereal time.
        let delta_days = (delta_sidereal_hours / 24) * SOLAR_DAYS_PER_SIDEREAL_DAY;
        time = time.AddDays(delta_days);
    }
}

/**
 * @brief When the seasons change for a given calendar year.
 *
 * Represents the dates and times of the two solstices
 * and the two equinoxes in a given calendar year.
 * These four events define the changing of the seasons on the Earth.
 *
 * @property {AstroTime} mar_equinox
 *      The date and time of the March equinox in the given calendar year.
 *      This is the moment in March that the plane of the Earth's equator passes
 *      through the center of the Sun; thus the Sun's declination
 *      changes from a negative number to a positive number.
 *      The March equinox defines
 *      the beginning of spring in the northern hemisphere and
 *      the beginning of autumn in the southern hemisphere.
 *
 * @property {AstroTime} jun_solstice
 *      The date and time of the June solstice in the given calendar year.
 *      This is the moment in June that the Sun reaches its most positive
 *      declination value.
 *      At this moment the Earth's north pole is most tilted most toward the Sun.
 *      The June solstice defines
 *      the beginning of summer in the northern hemisphere and
 *      the beginning of winter in the southern hemisphere.
 *
 * @property {AstroTime} sep_equinox
 *      The date and time of the September equinox in the given calendar year.
 *      This is the moment in September that the plane of the Earth's equator passes
 *      through the center of the Sun; thus the Sun's declination
 *      changes from a positive number to a negative number.
 *      The September equinox defines
 *      the beginning of autumn in the northern hemisphere and
 *      the beginning of spring in the southern hemisphere.
 *
 * @property {AstroTime} dec_solstice
 *      The date and time of the December solstice in the given calendar year.
 *      This is the moment in December that the Sun reaches its most negative
 *      declination value.
 *      At this moment the Earth's south pole is tilted most toward the Sun.
 *      The December solstice defines
 *      the beginning of winter in the northern hemisphere and
 *      the beginning of summer in the southern hemisphere.
 */
export class SeasonInfo {
    constructor(
        public mar_equinox: AstroTime,
        public jun_solstice: AstroTime,
        public sep_equinox: AstroTime,
        public dec_solstice: AstroTime)
        {}
}

/**
 * @brief Finds the equinoxes and solstices for a given calendar year.
 *
 * @param {number | AstroTime} year
 *      The integer value or `AstroTime` object that specifies
 *      the UTC calendar year for which to find equinoxes and solstices.
 *
 * @returns {SeasonInfo}
 */
export function Seasons(year: (number | AstroTime)): SeasonInfo {
    function find(targetLon: number, month: number, day: number): AstroTime {
        let startDate = new Date(Date.UTC(<number>year, month-1, day));
        let time = SearchSunLongitude(targetLon, startDate, 20);
        if (!time)
            throw `Cannot find season change near ${startDate.toISOString()}`;
        return time;
    }

    if ((year instanceof Date) && Number.isFinite(year.getTime()))
        year = year.getUTCFullYear();

    if (!Number.isSafeInteger(year))
        throw `Cannot calculate seasons because year argument ${year} is neither a Date nor a safe integer.`;

    let mar_equinox  = find(  0,  3, 10);
    let jun_solstice = find( 90,  6, 10);
    let sep_equinox  = find(180,  9, 10);
    let dec_solstice = find(270, 12, 10);

    return new SeasonInfo(mar_equinox, jun_solstice, sep_equinox, dec_solstice);
}

/**
 * @brief The viewing conditions of a body relative to the Sun.
 *
 * Represents the angular separation of a body from the Sun as seen from the Earth
 * and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.
 *
 * @property {AstroTime} time
 *      The date and time of the observation.
 *
 * @property {string}  visibility
 *      Either `"morning"` or `"evening"`,
 *      indicating when the body is most easily seen.
 *
 * @property {number}  elongation
 *      The angle in degrees, as seen from the center of the Earth,
 *      of the apparent separation between the body and the Sun.
 *      This angle is measured in 3D space and is not projected onto the ecliptic plane.
 *      When `elongation` is less than a few degrees, the body is very
 *      difficult to see from the Earth because it is lost in the Sun's glare.
 *      The elongation is always in the range [0, 180].
 *
 * @property {number}  ecliptic_separation
 *      The absolute value of the difference between the body's ecliptic longitude
 *      and the Sun's ecliptic longitude, both as seen from the center of the Earth.
 *      This angle measures around the plane of the Earth's orbit (the ecliptic),
 *      and ignores how far above or below that plane the body is.
 *      The ecliptic separation is measured in degrees and is always in the range [0, 180].
 *
 * @see {@link Elongation}
 */
export class ElongationEvent {
    constructor(
        public time: AstroTime,
        public visibility: string,
        public elongation: number,
        public ecliptic_separation: number)
        {}
}

/**
 * @brief Calculates the viewing conditions of a body relative to the Sun.
 *
 * Calculates angular separation of a body from the Sun as seen from the Earth
 * and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.
 * See the return type {@link ElongationEvent} for details.
 *
 * This function is helpful for determining how easy
 * it is to view a planet away from the Sun's glare on a given date.
 * It also determines whether the object is visible in the morning or evening;
 * this is more important the smaller the elongation is.
 * It is also used to determine how far a planet is from opposition, conjunction, or quadrature.
 *
 * @param {Body} body
 *      The name of the observed body. Not allowed to be `Body.Earth`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time of the observation.
 *
 * @returns {ElongationEvent}
 */
export function Elongation(body: Body, date: FlexibleDateTime): ElongationEvent {
    let time = MakeTime(date);

    let lon = PairLongitude(body, Body.Sun, time);
    let vis: string;
    if (lon > 180) {
        vis = 'morning';
        lon = 360 - lon;
    } else {
        vis = 'evening';
    }
    let angle = AngleFromSun(body, time);
    return new ElongationEvent(time, vis, angle, lon);
}

/**
 * @brief Finds the next time Mercury or Venus reaches maximum elongation.
 *
 * Searches for the next maximum elongation event for Mercury or Venus
 * that occurs after the given start date. Calling with other values
 * of `body` will result in an exception.
 * Maximum elongation occurs when the body has the greatest
 * angular separation from the Sun, as seen from the Earth.
 * Returns an `ElongationEvent` object containing the date and time of the next
 * maximum elongation, the elongation in degrees, and whether
 * the body is visible in the morning or evening.
 *
 * @param {Body} body
 *      Either `Body.Mercury` or `Body.Venus`.
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to search for the next maximum elongation event.
 *
 * @returns {ElongationEvent}
 */
export function SearchMaxElongation(body: Body, startDate: FlexibleDateTime): ElongationEvent {
    const dt = 0.01;

    function neg_slope(t: AstroTime): number {
        // The slope de/dt goes from positive to negative at the maximum elongation event.
        // But Search() is designed for functions that ascend through zero.
        // So this function returns the negative slope.
        const t1 = t.AddDays(-dt/2);
        const t2 = t.AddDays(+dt/2);
        let e1 = AngleFromSun(body, t1);
        let e2 = AngleFromSun(body, t2);
        let m = (e1-e2)/dt;
        return m;
    }

    let startTime = MakeTime(startDate);

    interface InferiorPlanetEntry {
        s1: number;
        s2: number;
    }

    interface InferiorPlanetTable {
        [body: string]: InferiorPlanetEntry;
    }

    const table: InferiorPlanetTable = {
        Mercury : { s1:50.0, s2:85.0 },
        Venus :   { s1:40.0, s2:50.0 }
    };

    const planet:InferiorPlanetEntry = table[body];
    if (!planet)
        throw 'SearchMaxElongation works for Mercury and Venus only.';

    let iter = 0;
    while (++iter <= 2) {
        // Find current heliocentric relative longitude between the
        // inferior planet and the Earth.
        let plon = EclipticLongitude(body, startTime);
        let elon = EclipticLongitude(Body.Earth, startTime);
        let rlon = LongitudeOffset(plon - elon);    // clamp to (-180, +180]

        // The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        // because there is a cusp there that causes a discontinuity in the derivative.
        // So we need to guard against searching near such times.

        let rlon_lo: number, rlon_hi: number, adjust_days: number;
        if (rlon >= -planet.s1 && rlon < +planet.s1 ) {
            // Seek to the window [+s1, +s2].
            adjust_days = 0;
            // Search forward for the time t1 when rel lon = +s1.
            rlon_lo = +planet.s1;
            // Search forward for the time t2 when rel lon = +s2.
            rlon_hi = +planet.s2;
        } else if (rlon >= +planet.s2 || rlon < -planet.s2 ) {
            // Seek to the next search window at [-s2, -s1].
            adjust_days = 0;
            // Search forward for the time t1 when rel lon = -s2.
            rlon_lo = -planet.s2;
            // Search forward for the time t2 when rel lon = -s1.
            rlon_hi = -planet.s1;
        } else if (rlon >= 0) {
            // rlon must be in the middle of the window [+s1, +s2].
            // Search BACKWARD for the time t1 when rel lon = +s1.
            adjust_days = -SynodicPeriod(body) / 4;
            rlon_lo = +planet.s1;
            rlon_hi = +planet.s2;
            // Search forward from t1 to find t2 such that rel lon = +s2.
        } else {
            // rlon must be in the middle of the window [-s2, -s1].
            // Search BACKWARD for the time t1 when rel lon = -s2.
            adjust_days = -SynodicPeriod(body) / 4;
            rlon_lo = -planet.s2;
            // Search forward from t1 to find t2 such that rel lon = -s1.
            rlon_hi = -planet.s1;
        }

        let t_start = startTime.AddDays(adjust_days);
        let t1 = SearchRelativeLongitude(body, rlon_lo, t_start);
        let t2 = SearchRelativeLongitude(body, rlon_hi, t1);

        // Now we have a time range [t1,t2] that brackets a maximum elongation event.
        // Confirm the bracketing.
        let m1 = neg_slope(t1);
        if (m1 >= 0)
            throw `SearchMaxElongation: internal error: m1 = ${m1}`;

        let m2 = neg_slope(t2);
        if (m2 <= 0)
            throw `SearchMaxElongation: internal error: m2 = ${m2}`;

        // Use the generic search algorithm to home in on where the slope crosses from negative to positive.
        let tx = Search(neg_slope, t1, t2, {init_f1:m1, init_f2:m2, dt_tolerance_seconds:10});
        if (!tx)
            throw `SearchMaxElongation: failed search iter ${iter} (t1=${t1.toString()}, t2=${t2.toString()})`;

        if (tx.tt >= startTime.tt)
            return Elongation(body, tx);

        // This event is in the past (earlier than startDate).
        // We need to search forward from t2 to find the next possible window.
        // We never need to search more than twice.
        startTime = t2.AddDays(1);
    }

    throw `SearchMaxElongation: failed to find event after 2 tries.`;
}

/**
 * @brief Searches for the date and time Venus will next appear brightest as seen from the Earth.
 *
 * @param {Body} body
 *      Currently only `Body.Venus` is supported.
 *      Mercury's peak magnitude occurs at superior conjunction, when it is impossible to see from Earth,
 *      so peak magnitude events have little practical value for that planet.
 *      The Moon reaches peak magnitude very close to full moon, which can be found using
 *      {@link SearchMoonQuarter} or {@link SearchMoonPhase}.
 *      The other planets reach peak magnitude very close to opposition,
 *      which can be found using {@link SearchRelativeLongitude}.
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to find the next peak magnitude event.
 *
 * @returns {IlluminationInfo}
 */
export function SearchPeakMagnitude(body: Body, startDate: FlexibleDateTime): IlluminationInfo {
    if (body !== Body.Venus)
        throw 'SearchPeakMagnitude currently works for Venus only.';

    const dt = 0.01;

    function slope(t: AstroTime): number {
        // The Search() function finds a transition from negative to positive values.
        // The derivative of magnitude y with respect to time t (dy/dt)
        // is negative as an object gets brighter, because the magnitude numbers
        // get smaller. At peak magnitude dy/dt = 0, then as the object gets dimmer,
        // dy/dt > 0.
        const t1 = t.AddDays(-dt/2);
        const t2 = t.AddDays(+dt/2);
        const y1 = Illumination(body, t1).mag;
        const y2 = Illumination(body, t2).mag;
        const m = (y2-y1) / dt;
        return m;
    }

    let startTime = MakeTime(startDate);

    // s1 and s2 are relative longitudes within which peak magnitude of Venus can occur.
    const s1 = 10.0;
    const s2 = 30.0;

    let iter = 0;
    while (++iter <= 2) {
        // Find current heliocentric relative longitude between the
        // inferior planet and the Earth.
        let plon = EclipticLongitude(body, startTime);
        let elon = EclipticLongitude(Body.Earth, startTime);
        let rlon = LongitudeOffset(plon - elon);    // clamp to (-180, +180]

        // The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        // because there is a cusp there that causes a discontinuity in the derivative.
        // So we need to guard against searching near such times.

        let rlon_lo: number, rlon_hi: number, adjust_days: number;
        if (rlon >= -s1 && rlon < +s1) {
            // Seek to the window [+s1, +s2].
            adjust_days = 0;
            // Search forward for the time t1 when rel lon = +s1.
            rlon_lo = +s1;
            // Search forward for the time t2 when rel lon = +s2.
            rlon_hi = +s2;
        } else if (rlon >= +s2 || rlon < -s2 ) {
            // Seek to the next search window at [-s2, -s1].
            adjust_days = 0;
            // Search forward for the time t1 when rel lon = -s2.
            rlon_lo = -s2;
            // Search forward for the time t2 when rel lon = -s1.
            rlon_hi = -s1;
        } else if (rlon >= 0) {
            // rlon must be in the middle of the window [+s1, +s2].
            // Search BACKWARD for the time t1 when rel lon = +s1.
            adjust_days = -SynodicPeriod(body) / 4;
            rlon_lo = +s1;
            // Search forward from t1 to find t2 such that rel lon = +s2.
            rlon_hi = +s2;
        } else {
            // rlon must be in the middle of the window [-s2, -s1].
            // Search BACKWARD for the time t1 when rel lon = -s2.
            adjust_days = -SynodicPeriod(body) / 4;
            rlon_lo = -s2;
            // Search forward from t1 to find t2 such that rel lon = -s1.
            rlon_hi = -s1;
        }

        let t_start = startTime.AddDays(adjust_days);
        let t1 = SearchRelativeLongitude(body, rlon_lo, t_start);
        let t2 = SearchRelativeLongitude(body, rlon_hi, t1);

        // Now we have a time range [t1,t2] that brackets a maximum magnitude event.
        // Confirm the bracketing.
        let m1 = slope(t1);
        if (m1 >= 0)
            throw `SearchPeakMagnitude: internal error: m1 = ${m1}`;

        let m2 = slope(t2);
        if (m2 <= 0)
            throw `SearchPeakMagnitude: internal error: m2 = ${m2}`;

        // Use the generic search algorithm to home in on where the slope crosses from negative to positive.
        let tx = Search(slope, t1, t2, {init_f1:m1, init_f2:m2, dt_tolerance_seconds:10});
        if (!tx)
            throw `SearchPeakMagnitude: failed search iter ${iter} (t1=${t1.toString()}, t2=${t2.toString()})`;

        if (tx.tt >= startTime.tt)
            return Illumination(body, tx);

        // This event is in the past (earlier than startDate).
        // We need to search forward from t2 to find the next possible window.
        // We never need to search more than twice.
        startTime = t2.AddDays(1);
    }

    throw `SearchPeakMagnitude: failed to find event after 2 tries.`;
}

/**
 * @brief The two kinds of apsis: pericenter (closest) and apocenter (farthest).
 *
 * `Pericenter`: The body is at its closest distance to the object it orbits.
 * `Apocenter`:  The body is at its farthest distance from the object it orbits.
 *
 * @enum {number}
 */
export enum ApsisKind {
    Pericenter = 0,
    Apocenter = 1
}

/**
 * @brief A closest or farthest point in a body's orbit around its primary.
 *
 * For a planet orbiting the Sun, apsis is a perihelion or aphelion, respectively.
 * For the Moon orbiting the Earth, apsis is a perigee or apogee, respectively.
 *
 * @property {AstroTime} time
 *      The date and time of the apsis.
 *
 * @property {ApsisKind} kind
 *      For a closest approach (perigee or perihelion), `kind` is `ApsisKind.Pericenter`.
 *      For a farthest distance event (apogee or aphelion), `kind` is `ApsisKind.Apocenter`.
 *
 * @property {number} dist_au
 *      The distance between the centers of the two bodies in astronomical units (AU).
 *
 * @property {number} dist_km
 *      The distance between the centers of the two bodies in kilometers.
 *
 * @see {@link SearchLunarApsis}
 * @see {@link NextLunarApsis}
 * @see {@link SearchPlanetApsis}
 * @see {@link NextPlanetApsis}
 */
export class Apsis {
    dist_km: number;

    constructor(
        public time: AstroTime,
        public kind: ApsisKind,
        public dist_au: number)
    {
        this.dist_km = dist_au * KM_PER_AU;
    }
}

/**
 * @brief Finds the next perigee or apogee of the Moon.
 *
 * Finds the next perigee (closest approach) or apogee (farthest remove) of the Moon
 * that occurs after the specified date and time.
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to find the next perigee or apogee.
 *
 * @returns {Apsis}
 */
export function SearchLunarApsis(startDate: FlexibleDateTime): Apsis {
    const dt = 0.001;

    function distance_slope(t: AstroTime): number {
        let t1 = t.AddDays(-dt/2);
        let t2 = t.AddDays(+dt/2);

        let r1 = CalcMoon(t1).distance_au;
        let r2 = CalcMoon(t2).distance_au;

        let m = (r2-r1) / dt;
        return m;
    }

    function negative_distance_slope(t: AstroTime): number {
        return -distance_slope(t);
    }

    // Check the rate of change of the distance dr/dt at the start time.
    // If it is positive, the Moon is currently getting farther away,
    // so start looking for apogee.
    // Conversely, if dr/dt < 0, start looking for perigee.
    // Either way, the polarity of the slope will change, so the product will be negative.
    // Handle the crazy corner case of exactly touching zero by checking for m1*m2 <= 0.

    let t1 = MakeTime(startDate);
    let m1 = distance_slope(t1);
    const increment = 5;      // number of days to skip in each iteration

    for (var iter = 0; iter * increment < 2 * MEAN_SYNODIC_MONTH; ++iter) {
        let t2 = t1.AddDays(increment);
        let m2 = distance_slope(t2);

        if (m1 * m2 <= 0) {
            // The time range [t1, t2] contains an apsis.
            // Figure out whether it is perigee or apogee.

            if (m1 < 0 || m2 > 0) {
                // We found a minimum distance event: perigee.
                // Search the time range [t1, t2] for the time when the slope goes
                // from negative to positive.
                let tx = Search(distance_slope, t1, t2, {init_f1:m1, init_f2:m2});
                if (!tx)
                    throw 'SearchLunarApsis INTERNAL ERROR: perigee search failed!';

                let dist = CalcMoon(tx).distance_au;
                return new Apsis(tx, 0, dist);
            }

            if (m1 > 0 || m2 < 0) {
                // We found a maximum distance event: apogee.
                // Search the time range [t1, t2] for the time when the slope goes
                // from positive to negative.
                let tx = Search(negative_distance_slope, t1, t2, {init_f1:-m1, init_f2:-m2});
                if (!tx)
                    throw 'SearchLunarApsis INTERNAL ERROR: apogee search failed!';

                let dist = CalcMoon(tx).distance_au;
                return new Apsis(tx, 1, dist);
            }

            // This should never happen; it should not be possible for consecutive
            // times t1 and t2 to both have zero slope.
            throw 'SearchLunarApsis INTERNAL ERROR: cannot classify apsis event!';
        }

        t1 = t2;
        m1 = m2;
    }

    // It should not be possible to fail to find an apsis within 2 synodic months.
    throw 'SearchLunarApsis INTERNAL ERROR: could not find apsis within 2 synodic months of start date.';
}

/**
 * @brief Finds the next lunar apsis (perigee or apogee) in a series.
 *
 * Given a lunar apsis returned by an initial call to {@link SearchLunarApsis},
 * or a previous call to `NextLunarApsis`, finds the next lunar apsis.
 * If the given apsis is a perigee, this function finds the next apogee, and vice versa.
 *
 * @param {Apsis} apsis
 *      A lunar perigee or apogee event.
 *
 * @returns {Apsis}
 *      The successor apogee for the given perigee, or the successor perigee for the given apogee.
 */
export function NextLunarApsis(apsis: Apsis): Apsis {
    const skip = 11;    // number of days to skip to start looking for next apsis event
    let next = SearchLunarApsis(apsis.time.AddDays(skip));
    if (next.kind + apsis.kind !== 1)
        throw `NextLunarApsis INTERNAL ERROR: did not find alternating apogee/perigee: prev=${apsis.kind} @ ${apsis.time.toString()}, next=${next.kind} @ ${next.time.toString()}`;
    return next;
}

function PlanetExtreme(body: Body, kind: ApsisKind, start_time: AstroTime, dayspan: number): Apsis {
    const direction = (kind === ApsisKind.Apocenter) ? +1.0 : -1.0;
    const npoints = 10;

    for(;;) {
        const interval = dayspan / (npoints - 1);

        // iterate until uncertainty is less than one minute
        if (interval < 1.0 / 1440.0) {
            const apsis_time = start_time.AddDays(interval / 2.0);
            const dist_au = HelioDistance(body, apsis_time);
            return new Apsis(apsis_time, kind, dist_au);
        }

        let best_i = -1;
        let best_dist = 0.0;
        for (let i=0; i < npoints; ++i) {
            const time = start_time.AddDays(i * interval);
            const dist = direction * HelioDistance(body, time);
            if (i==0 || dist > best_dist) {
                best_i = i;
                best_dist = dist;
            }
        }

        /* Narrow in on the extreme point. */
        start_time = start_time.AddDays((best_i - 1) * interval);
        dayspan = 2.0 * interval;
    }
}

function BruteSearchPlanetApsis(body: Body, startTime: AstroTime): Apsis {
    /*
        Neptune is a special case for two reasons:
        1. Its orbit is nearly circular (low orbital eccentricity).
        2. It is so distant from the Sun that the orbital period is very long.
        Put together, this causes wobbling of the Sun around the Solar System Barycenter (SSB)
        to be so significant that there are 3 local minima in the distance-vs-time curve
        near each apsis. Therefore, unlike for other planets, we can't use an optimized
        algorithm for finding dr/dt = 0.
        Instead, we use a dumb, brute-force algorithm of sampling and finding min/max
        heliocentric distance.

        There is a similar problem in the TOP2013 model for Pluto:
        Its position vector has high-frequency oscillations that confuse the
        slope-based determination of apsides.
    */

    /*
        Rewind approximately 30 degrees in the orbit,
        then search forward for 270 degrees.
        This is a very cautious way to prevent missing an apsis.
        Typically we will find two apsides, and we pick whichever
        apsis is ealier, but after startTime.
        Sample points around this orbital arc and find when the distance
        is greatest and smallest.
    */
    const npoints = 100;
    const t1 = startTime.AddDays(Planet[body].OrbitalPeriod * ( -30 / 360));
    const t2 = startTime.AddDays(Planet[body].OrbitalPeriod * (+270 / 360));
    let t_min = t1;
    let t_max = t1;
    let min_dist = -1.0;
    let max_dist = -1.0;
    const interval = (t2.ut - t1.ut) / (npoints - 1);

    for (let i=0; i < npoints; ++i) {
        const time = t1.AddDays(i * interval);
        const dist = HelioDistance(body, time);
        if (i === 0) {
            max_dist = min_dist = dist;
        } else {
            if (dist > max_dist) {
                max_dist = dist;
                t_max = time;
            }
            if (dist < min_dist) {
                min_dist = dist;
                t_min = time;
            }
        }
    }

    const perihelion = PlanetExtreme(body, 0, t_min.AddDays(-2*interval), 4*interval);
    const aphelion   = PlanetExtreme(body, 1, t_max.AddDays(-2*interval), 4*interval);
    if (perihelion.time.tt >= startTime.tt) {
        if (aphelion.time.tt >= startTime.tt && aphelion.time.tt < perihelion.time.tt)
            return aphelion;
        return perihelion;
    }
    if (aphelion.time.tt >= startTime.tt)
        return aphelion;
    throw 'Internal error: failed to find Neptune apsis.';
}

/**
 * @brief Finds the next perihelion or aphelion of a planet.
 *
 * Finds the date and time of a planet's perihelion (closest approach to the Sun)
 * or aphelion (farthest distance from the Sun) after a given time.
 *
 * Given a date and time to start the search in `startTime`, this function finds the
 * next date and time that the center of the specified planet reaches the closest or farthest point
 * in its orbit with respect to the center of the Sun, whichever comes first
 * after `startTime`.
 *
 * The closest point is called *perihelion* and the farthest point is called *aphelion*.
 * The word *apsis* refers to either event.
 *
 * To iterate through consecutive alternating perihelion and aphelion events,
 * call `SearchPlanetApsis` once, then use the return value to call
 * {@link NextPlanetApsis}. After that, keep feeding the previous return value
 * from `NextPlanetApsis` into another call of `NextPlanetApsis`
 * as many times as desired.
 *
 * @param {Body} body
 *      The planet for which to find the next perihelion/aphelion event.
 *      Not allowed to be `Body.Sun` or `Body.Moon`.
 *
 * @param {FlexibleDateTime} startTime
 *      The date and time at which to start searching for the next perihelion or aphelion.
 *
 * @returns {Apsis}
 *      The next perihelion or aphelion that occurs after `startTime`.
 */
export function SearchPlanetApsis(body: Body, startTime: FlexibleDateTime): Apsis {
    startTime = MakeTime(startTime);
    if (body === Body.Neptune || body === Body.Pluto)
        return BruteSearchPlanetApsis(body, startTime);

    function positive_slope(t: AstroTime): number {
        const dt = 0.001;
        let t1 = t.AddDays(-dt/2);
        let t2 = t.AddDays(+dt/2);
        let r1 = HelioDistance(body, t1);
        let r2 = HelioDistance(body, t2);
        let m = (r2-r1) / dt;
        return m;
    }

    function negative_slope(t: AstroTime): number {
        return -positive_slope(t);
    }

    const orbit_period_days = Planet[body].OrbitalPeriod;
    const increment = orbit_period_days / 6.0;
    let t1 = startTime;
    let m1 = positive_slope(t1);
    for (let iter = 0; iter * increment < 2.0 * orbit_period_days; ++iter) {
        const t2 = t1.AddDays(increment);
        const m2 = positive_slope(t2);
        if (m1 * m2 <= 0.0) {
            /* There is a change of slope polarity within the time range [t1, t2]. */
            /* Therefore this time range contains an apsis. */
            /* Figure out whether it is perihelion or aphelion. */

            let slope_func: (t: AstroTime) => number;
            let kind: ApsisKind;
            if (m1 < 0.0 || m2 > 0.0) {
                /* We found a minimum-distance event: perihelion. */
                /* Search the time range for the time when the slope goes from negative to positive. */
                slope_func = positive_slope;
                kind = ApsisKind.Pericenter;
            } else if (m1 > 0.0 || m2 < 0.0) {
                /* We found a maximum-distance event: aphelion. */
                /* Search the time range for the time when the slope goes from positive to negative. */
                slope_func = negative_slope;
                kind = ApsisKind.Apocenter;
            } else {
                /* This should never happen. It should not be possible for both slopes to be zero. */
                throw "Internal error with slopes in SearchPlanetApsis";
            }

            const search = Search(slope_func, t1, t2);
            if (!search)
                throw "Failed to find slope transition in planetary apsis search.";

            const dist = HelioDistance(body, search);
            return new Apsis(search, kind, dist);
        }
        /* We have not yet found a slope polarity change. Keep searching. */
        t1 = t2;
        m1 = m2;
    }
    throw "Internal error: should have found planetary apsis within 2 orbital periods.";
}

/**
 * @brief Finds the next planetary perihelion or aphelion event in a series.
 *
 * This function requires an {@link Apsis} value obtained from a call
 * to {@link SearchPlanetApsis} or `NextPlanetApsis`.
 * Given an aphelion event, this function finds the next perihelion event, and vice versa.
 * See {@link SearchPlanetApsis} for more details.
 *
 * @param {Body} body
 *      The planet for which to find the next perihelion/aphelion event.
 *      Not allowed to be `Body.Sun` or `Body.Moon`.
 *      Must match the body passed into the call that produced the `apsis` parameter.
 *
 * @param {Apsis} apsis
 *      An apsis event obtained from a call to {@link SearchPlanetApsis} or `NextPlanetApsis`.
 *
 * @returns {Apsis}
 *      Same as the return value for {@link SearchPlanetApsis}.
 */
export function NextPlanetApsis(body: Body, apsis: Apsis): Apsis {
    if (apsis.kind !== ApsisKind.Pericenter && apsis.kind !== ApsisKind.Apocenter)
        throw `Invalid apsis kind: ${apsis.kind}`;

    /* skip 1/4 of an orbit before starting search again */
    const skip = 0.25 * Planet[body].OrbitalPeriod;
    const time = apsis.time.AddDays(skip);
    const next = SearchPlanetApsis(body, time);

    /* Verify that we found the opposite apsis from the previous one. */
    if (next.kind + apsis.kind !== 1)
        throw `Internal error: previous apsis was ${apsis.kind}, but found ${next.kind} for next apsis.`;

    return next;
}

/**
 * @brief Calculates the inverse of a rotation matrix.
 *
 * Given a rotation matrix that performs some coordinate transform,
 * this function returns the matrix that reverses that transform.
 *
 * @param {RotationMatrix} rotation
 *      The rotation matrix to be inverted.
 *
 * @returns {RotationMatrix}
 *      The inverse rotation matrix.
 */
export function InverseRotation(rotation: RotationMatrix): RotationMatrix {
    return new RotationMatrix([
        [rotation.rot[0][0], rotation.rot[1][0], rotation.rot[2][0]],
        [rotation.rot[0][1], rotation.rot[1][1], rotation.rot[2][1]],
        [rotation.rot[0][2], rotation.rot[1][2], rotation.rot[2][2]]
    ]);
}

/**
 * @brief Creates a rotation based on applying one rotation followed by another.
 *
 * Given two rotation matrices, returns a combined rotation matrix that is
 * equivalent to rotating based on the first matrix, followed by the second.
 *
 * @param {RotationMatrix} a
 *      The first rotation to apply.
 *
 * @param {RotationMatrix} b
 *      The second rotation to apply.
 *
 * @returns {RotationMatrix}
 *      The combined rotation matrix.
 */
export function CombineRotation(a: RotationMatrix, b: RotationMatrix): RotationMatrix {
    /*
        Use matrix multiplication: c = b*a.
        We put 'b' on the left and 'a' on the right because,
        just like when you use a matrix M to rotate a vector V,
        you put the M on the left in the product M*V.
        We can think of this as 'b' rotating all the 3 column vectors in 'a'.
    */

    return new RotationMatrix([
        [
            b.rot[0][0]*a.rot[0][0] + b.rot[1][0]*a.rot[0][1] + b.rot[2][0]*a.rot[0][2],
            b.rot[0][1]*a.rot[0][0] + b.rot[1][1]*a.rot[0][1] + b.rot[2][1]*a.rot[0][2],
            b.rot[0][2]*a.rot[0][0] + b.rot[1][2]*a.rot[0][1] + b.rot[2][2]*a.rot[0][2]
        ],
        [
            b.rot[0][0]*a.rot[1][0] + b.rot[1][0]*a.rot[1][1] + b.rot[2][0]*a.rot[1][2],
            b.rot[0][1]*a.rot[1][0] + b.rot[1][1]*a.rot[1][1] + b.rot[2][1]*a.rot[1][2],
            b.rot[0][2]*a.rot[1][0] + b.rot[1][2]*a.rot[1][1] + b.rot[2][2]*a.rot[1][2]
        ],
        [
            b.rot[0][0]*a.rot[2][0] + b.rot[1][0]*a.rot[2][1] + b.rot[2][0]*a.rot[2][2],
            b.rot[0][1]*a.rot[2][0] + b.rot[1][1]*a.rot[2][1] + b.rot[2][1]*a.rot[2][2],
            b.rot[0][2]*a.rot[2][0] + b.rot[1][2]*a.rot[2][1] + b.rot[2][2]*a.rot[2][2]
        ]
    ]);
}

/**
 * @brief Creates an identity rotation matrix.
 *
 * Returns a rotation matrix that has no effect on orientation.
 * This matrix can be the starting point for other operations,
 * such as using a series of calls to {@link Pivot} to
 * create a custom rotation matrix.
 *
 * @returns {RotationMatrix}
 *      The identity matrix.
 */
 export function IdentityMatrix(): RotationMatrix {
     return new RotationMatrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
     ]);
 }

 /**
 * @brief Re-orients a rotation matrix by pivoting it by an angle around one of its axes.
 *
 * Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
 * this function pivots the rotation matrix by that angle around that coordinate axis.
 *
 * For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
 * to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
 * of a telescope camera pointed at a given body, you can use `Astronomy_Pivot` twice:
 * (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
 * western axis by the body's altitude angle. The resulting rotation matrix will then
 * reorient ECL coordinates to the orientation of your telescope camera.
 *
 * @param {RotationMatrix} rotation
 *      The input rotation matrix.
 *
 * @param {number} axis
 *      An integer that selects which coordinate axis to rotate around:
 *      0 = x, 1 = y, 2 = z. Any other value will cause an exception.
 *
 * @param {number} angle
 *      An angle in degrees indicating the amount of rotation around the specified axis.
 *      Positive angles indicate rotation counterclockwise as seen from the positive
 *      direction along that axis, looking towards the origin point of the orientation system.
 *      Any finite number of degrees is allowed, but best precision will result from
 *      keeping `angle` in the range [-360, +360].
 *
 * @returns {RotationMatrix}
 *      A pivoted matrix object.
 */
export function Pivot(rotation: RotationMatrix, axis: 0 | 1 | 2, angle: number): RotationMatrix {
    // Check for an invalid coordinate axis.
    if (axis !== 0 && axis !== 1 && axis !== 2)
        throw `Invalid axis ${axis}. Must be [0, 1, 2].`;

    const radians = VerifyNumber(angle) * DEG2RAD;
    const c = Math.cos(radians);
    const s = Math.sin(radians);

    /*
        We need to maintain the "right-hand" rule, no matter which
        axis was selected. That means we pick (i, j, k) axis order
        such that the following vector cross product is satisfied:
        i x j = k
    */
    const i = (axis + 1) % 3;
    const j = (axis + 2) % 3;
    const k = axis;

    let rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

    rot[i][i] = c*rotation.rot[i][i] - s*rotation.rot[i][j];
    rot[i][j] = s*rotation.rot[i][i] + c*rotation.rot[i][j];
    rot[i][k] = rotation.rot[i][k];

    rot[j][i] = c*rotation.rot[j][i] - s*rotation.rot[j][j];
    rot[j][j] = s*rotation.rot[j][i] + c*rotation.rot[j][j];
    rot[j][k] = rotation.rot[j][k];

    rot[k][i] = c*rotation.rot[k][i] - s*rotation.rot[k][j];
    rot[k][j] = s*rotation.rot[k][i] + c*rotation.rot[k][j];
    rot[k][k] = rotation.rot[k][k];

    return new RotationMatrix(rot);
}

/**
 * @brief Converts spherical coordinates to Cartesian coordinates.
 *
 * Given spherical coordinates and a time at which they are valid,
 * returns a vector of Cartesian coordinates. The returned value
 * includes the time, as required by `AstroTime`.
 *
 * @param {Spherical} sphere
 *      Spherical coordinates to be converted.
 *
 * @param {FlexibleDateTime} time
 *      The time that should be included in the returned vector.
 *
 * @returns {Vector}
 *      The vector form of the supplied spherical coordinates.
 */
export function VectorFromSphere(sphere: Spherical, time: FlexibleDateTime): Vector {
    time = MakeTime(time);
    const radlat = sphere.lat * DEG2RAD;
    const radlon = sphere.lon * DEG2RAD;
    const rcoslat = sphere.dist * Math.cos(radlat);
    return new Vector(
        rcoslat * Math.cos(radlon),
        rcoslat * Math.sin(radlon),
        sphere.dist * Math.sin(radlat),
        time
    );
}

/**
 * @brief Given an equatorial vector, calculates equatorial angular coordinates.
 *
 * @param {Vector} vec
 *      A vector in an equatorial coordinate system.
 *
 * @returns {EquatorialCoordinates}
 *      Angular coordinates expressed in the same equatorial system as `vec`.
 */
export function EquatorFromVector(vec: Vector): EquatorialCoordinates {
    const sphere = SphereFromVector(vec);
    return new EquatorialCoordinates(sphere.lon / 15, sphere.lat, sphere.dist, vec);
}

/**
 * @brief Converts Cartesian coordinates to spherical coordinates.
 *
 * Given a Cartesian vector, returns latitude, longitude, and distance.
 *
 * @param {Vector} vector
 *      Cartesian vector to be converted to spherical coordinates.
 *
 * @returns {Spherical}
 *      Spherical coordinates that are equivalent to the given vector.
 */
export function SphereFromVector(vector: Vector): Spherical {
    const xyproj = vector.x*vector.x + vector.y*vector.y;
    const dist = Math.sqrt(xyproj + vector.z*vector.z);
    let lat: number, lon: number;
    if (xyproj === 0.0) {
        if (vector.z === 0.0)
            throw 'Zero-length vector not allowed.';
        lon = 0.0;
        lat = (vector.z < 0.0) ? -90.0 : +90.0;
    } else {
        lon = RAD2DEG * Math.atan2(vector.y, vector.x);
        if (lon < 0.0)
            lon += 360.0;
        lat = RAD2DEG * Math.atan2(vector.z, Math.sqrt(xyproj));
    }
    return new Spherical(lat, lon, dist);
}

function ToggleAzimuthDirection(az: number): number {
    az = 360.0 - az;
    if (az >= 360.0)
        az -= 360.0;
    else if (az < 0.0)
        az += 360.0;
    return az;
}

/**
 * @brief Converts Cartesian coordinates to horizontal coordinates.
 *
 * Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
 *
 * *IMPORTANT:* This function differs from {@link SphereFromVector} in two ways:
 * - `SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
 *   from north (e.g., west = +90), but this function represents a clockwise rotation
 *   (e.g., east = +90). The difference is because `SphereFromVector` is intended
 *   to preserve the vector "right-hand rule", while this function defines azimuth in a more
 *   traditional way as used in navigation and cartography.
 * - This function optionally corrects for atmospheric refraction, while `SphereFromVector` does not.
 *
 * The returned object contains the azimuth in `lon`.
 * It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
 *
 * The altitude is stored in `lat`.
 *
 * The distance to the observed object is stored in `dist`,
 * and is expressed in astronomical units (AU).
 *
 * @param {Vector} vector
 *      Cartesian vector to be converted to horizontal coordinates.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @returns {Spherical}
 */
export function HorizonFromVector(vector: Vector, refraction: string): Spherical {
    const sphere = SphereFromVector(vector);
    sphere.lon = ToggleAzimuthDirection(sphere.lon);
    sphere.lat += Refraction(refraction, sphere.lat);
    return sphere;
}


/**
 * @brief Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.
 *
 * @param {Spherical} sphere
 *      A structure that contains apparent horizontal coordinates:
 *      `lat` holds the refracted azimuth angle,
 *      `lon` holds the azimuth in degrees clockwise from north,
 *      and `dist` holds the distance from the observer to the object in AU.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the observation. This is needed because the returned
 *      vector object requires a valid time value when passed to certain other functions.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @returns {Vector}
 *      A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
 */
export function VectorFromHorizon(sphere: Spherical, time: FlexibleDateTime, refraction: string): Vector {
    time = MakeTime(time);

    /* Convert azimuth from clockwise-from-north to counterclockwise-from-north. */
    const lon = ToggleAzimuthDirection(sphere.lon);

    /* Reverse any applied refraction. */
    const lat = sphere.lat + InverseRefraction(refraction, sphere.lat);

    const xsphere = new Spherical(lat, lon, sphere.dist);
    return VectorFromSphere(xsphere, time);
}


/**
 * @brief Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
 *
 * Given an altitude angle and a refraction option, calculates
 * the amount of "lift" caused by atmospheric refraction.
 * This is the number of degrees higher in the sky an object appears
 * due to the lensing of the Earth's atmosphere.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @param {number} altitude
 *      An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90.
 *
 * @returns {number}
 *      The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing.
 */
export function Refraction(refraction: string, altitude: number): number {
    let refr: number;

    VerifyNumber(altitude);

    if (altitude < -90.0 || altitude > +90.0)
        return 0.0;     /* no attempt to correct an invalid altitude */

    if (refraction === 'normal' || refraction === 'jplhor') {
        // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        // JPL Horizons says it uses refraction algorithm from
        // Meeus "Astronomical Algorithms", 1991, p. 101-102.
        // I found the following Go implementation:
        // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        // This is a translation from the function "Saemundsson" there.
        // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        // This is important because the 'refr' formula below goes crazy near hd = -5.11.
        let hd = altitude;
        if (hd < -1.0)
            hd = -1.0;

        refr = (1.02 / Math.tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60.0;

        if (refraction === 'normal' && altitude < -1.0) {
            // In "normal" mode we gradually reduce refraction toward the nadir
            // so that we never get an altitude angle less than -90 degrees.
            // When horizon angle is -1 degrees, the factor is exactly 1.
            // As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
            refr *= (altitude + 90.0) / 89.0;
        }
    } else if (!refraction) {
        // The caller does not want refraction correction.
        refr = 0.0;
    } else {
        throw `Invalid refraction option: ${refraction}`;
    }

    return refr;
}

/**
 * @brief Calculates the inverse of an atmospheric refraction angle.
 *
 * Given an observed altitude angle that includes atmospheric refraction,
 * calculates the negative angular correction to obtain the unrefracted
 * altitude. This is useful for cases where observed horizontal
 * coordinates are to be converted to another orientation system,
 * but refraction first must be removed from the observed position.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @param {number} bent_altitude
 *      The apparent altitude that includes atmospheric refraction.
 *
 * @returns {number}
 *      The angular adjustment in degrees to be added to the
 *      altitude angle to correct for atmospheric lensing.
 *      This will be less than or equal to zero.
 */
export function InverseRefraction(refraction: string, bent_altitude: number): number {
    if (bent_altitude < -90.0 || bent_altitude > +90.0)
        return 0.0;     /* no attempt to correct an invalid altitude */

    /* Find the pre-adjusted altitude whose refraction correction leads to 'altitude'. */
    let altitude = bent_altitude - Refraction(refraction, bent_altitude);

    for(;;) {
        /* See how close we got. */
        let diff = (altitude + Refraction(refraction, altitude)) - bent_altitude;
        if (Math.abs(diff) < 1.0e-14)
            return altitude - bent_altitude;

        altitude -= diff;
    }
}

/**
 * @brief Applies a rotation to a vector, yielding a rotated vector.
 *
 * This function transforms a vector in one orientation to a vector
 * in another orientation.
 *
 * @param {RotationMatrix} rotation
 *      A rotation matrix that specifies how the orientation of the vector is to be changed.
 *
 * @param {Vector} vector
 *      The vector whose orientation is to be changed.
 *
 * @returns {Vector}
 *      A vector in the orientation specified by `rotation`.
 */
export function RotateVector(rotation: RotationMatrix, vector: Vector): Vector
{
    return new Vector(
        rotation.rot[0][0]*vector.x + rotation.rot[1][0]*vector.y + rotation.rot[2][0]*vector.z,
        rotation.rot[0][1]*vector.x + rotation.rot[1][1]*vector.y + rotation.rot[2][1]*vector.z,
        rotation.rot[0][2]*vector.x + rotation.rot[1][2]*vector.y + rotation.rot[2][2]*vector.z,
        vector.t
    );
}


/**
 * @brief Applies a rotation to a state vector, yielding a rotated vector.
 *
 * This function transforms a state vector in one orientation to a vector
 * in another orientation.
 *
 * @param {RotationMatrix} rotation
 *      A rotation matrix that specifies how the orientation of the state vector is to be changed.
 *
 * @param {StateVector} state
 *      The state vector whose orientation is to be changed.
 *      Both the position and velocity components are transformed.
 *
 * @return {StateVector}
 *      A state vector in the orientation specified by `rotation`.
 */
export function RotateState(rotation: RotationMatrix, state: StateVector): StateVector
{
    return new StateVector(
        rotation.rot[0][0]*state.x + rotation.rot[1][0]*state.y + rotation.rot[2][0]*state.z,
        rotation.rot[0][1]*state.x + rotation.rot[1][1]*state.y + rotation.rot[2][1]*state.z,
        rotation.rot[0][2]*state.x + rotation.rot[1][2]*state.y + rotation.rot[2][2]*state.z,

        rotation.rot[0][0]*state.vx + rotation.rot[1][0]*state.vy + rotation.rot[2][0]*state.vz,
        rotation.rot[0][1]*state.vx + rotation.rot[1][1]*state.vy + rotation.rot[2][1]*state.vz,
        rotation.rot[0][2]*state.vx + rotation.rot[1][2]*state.vy + rotation.rot[2][2]*state.vz,

        state.t
    );
}



/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQJ to ECL.
 */
export function Rotation_EQJ_ECL(): RotationMatrix {
    /* ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians. */
    const c = 0.9174821430670688;    /* cos(ob) */
    const s = 0.3977769691083922;    /* sin(ob) */
    return new RotationMatrix([
        [1,  0,  0],
        [0, +c, -s],
        [0, +s, +c]
    ]);
}


/**
 * @brief Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts ECL to EQJ.
 */
export function Rotation_ECL_EQJ(): RotationMatrix {
    /* ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians. */
    const c = 0.9174821430670688;    /* cos(ob) */
    const s = 0.3977769691083922;    /* sin(ob) */
    return new RotationMatrix([
        [ 1,  0,  0],
        [ 0, +c, +s],
        [ 0, -s, +c]
    ]);
}


/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param {FlexibleDateTime} time
 *      The date and time at which the Earth's equator defines the target orientation.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQJ to EQD at `time`.
 */
export function Rotation_EQJ_EQD(time: FlexibleDateTime): RotationMatrix {
    time = MakeTime(time);
    const prec = precession_rot(time, PrecessDirection.From2000);
    const nut = nutation_rot(time, PrecessDirection.From2000);
    return CombineRotation(prec, nut);
}


/**
 * @brief Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @param {FlexibleDateTime} time
 *      The date and time at which the Earth's equator defines the source orientation.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQD at `time` to EQJ.
 */
export function Rotation_EQD_EQJ(time: FlexibleDateTime): RotationMatrix {
    time = MakeTime(time);
    const nut = nutation_rot(time, PrecessDirection.Into2000);
    const prec = precession_rot(time, PrecessDirection.Into2000);
    return CombineRotation(nut, prec);
}


/**
 * @brief Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: HOR = horizontal system.
 *
 * Use `HorizonFromVector` to convert the return value
 * to a traditional altitude/azimuth pair.
 *
 * @param {FlexibleDateTime} time
 *      The date and time at which the Earth's equator applies.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQD to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
export function Rotation_EQD_HOR(time: FlexibleDateTime, observer: Observer): RotationMatrix {
    time = MakeTime(time);
    const sinlat = Math.sin(observer.latitude * DEG2RAD);
    const coslat = Math.cos(observer.latitude * DEG2RAD);
    const sinlon = Math.sin(observer.longitude * DEG2RAD);
    const coslon = Math.cos(observer.longitude * DEG2RAD);

    const uze: ArrayVector = [coslat * coslon, coslat * sinlon, sinlat];
    const une: ArrayVector = [-sinlat * coslon, -sinlat * sinlon, coslat];
    const uwe: ArrayVector = [sinlon, -coslon, 0];

    const spin_angle = -15 * sidereal_time(time);
    const uz = spin(spin_angle, uze);
    const un = spin(spin_angle, une);
    const uw = spin(spin_angle, uwe);

    return new RotationMatrix([
        [un[0], uw[0], uz[0]],
        [un[1], uw[1], uz[1]],
        [un[2], uw[2], uz[2]],
    ]);
}


/**
 * @brief Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param {FlexibleDateTime} time
 *      The date and time at which the Earth's equator applies.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts HOR to EQD at `time` and for `observer`.
 */
export function Rotation_HOR_EQD(time: FlexibleDateTime, observer: Observer): RotationMatrix {
    const rot = Rotation_EQD_HOR(time, observer);
    return InverseRotation(rot);
}


/**
 * @brief Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQJ = equatorial system, using equator at the J2000 epoch.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the observation.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts HOR to EQJ at `time` and for `observer`.
 */
export function Rotation_HOR_EQJ(time: FlexibleDateTime, observer: Observer): RotationMatrix {
    time = MakeTime(time);
    const hor_eqd = Rotation_HOR_EQD(time, observer);
    const eqd_eqj = Rotation_EQD_EQJ(time);
    return CombineRotation(hor_eqd, eqd_eqj);
}


/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using the equator at the J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * Use {@link HorizonFromVector} to convert the return value
 * to a traditional altitude/azimuth pair.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the desired horizontal orientation.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
export function Rotation_EQJ_HOR(time: FlexibleDateTime, observer: Observer): RotationMatrix {
    const rot = Rotation_HOR_EQJ(time, observer);
    return InverseRotation(rot);
}


/**
 * @brief Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of date.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the source equator.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQD to ECL.
 */
export function Rotation_EQD_ECL(time: FlexibleDateTime): RotationMatrix {
    const eqd_eqj = Rotation_EQD_EQJ(time);
    const eqj_ecl = Rotation_EQJ_ECL();
    return CombineRotation(eqd_eqj, eqj_ecl);
}


/**
 * @brief Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of date.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the desired equator.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts ECL to EQD.
 */
export function Rotation_ECL_EQD(time: FlexibleDateTime): RotationMatrix {
    const rot = Rotation_EQD_ECL(time);
    return InverseRotation(rot);
}


/**
 * @brief Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * Use {@link HorizonFromVector} to convert the return value
 * to a traditional altitude/azimuth pair.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the desired horizontal orientation.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts ECL to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
export function Rotation_ECL_HOR(time: FlexibleDateTime, observer: Observer): RotationMatrix {
    time = MakeTime(time);
    const ecl_eqd = Rotation_ECL_EQD(time);
    const eqd_hor = Rotation_EQD_HOR(time, observer);
    return CombineRotation(ecl_eqd, eqd_hor);
}


/**
 * @brief Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param {FlexibleDateTime} time
 *      The date and time of the horizontal observation.
 *
 * @param {Observer} observer
 *      The location of the horizontal observer.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts HOR to ECL.
 */
export function Rotation_HOR_ECL(time: FlexibleDateTime, observer: Observer): RotationMatrix {
    const rot = Rotation_ECL_HOR(time, observer);
    return InverseRotation(rot);
}


/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to galactic (GAL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using the equator at the J2000 epoch.
 * Target: GAL = galactic system (IAU 1958 definition).
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQJ to GAL.
 */
export function Rotation_EQJ_GAL(): RotationMatrix {
    // This rotation matrix was calculated by the following script
    // in this same source code repository:
    // demo/python/galeqj_matrix.py
    return new RotationMatrix([
        [-0.0548624779711344, +0.4941095946388765, -0.8676668813529025],
        [-0.8734572784246782, -0.4447938112296831, -0.1980677870294097],
        [-0.4838000529948520, +0.7470034631630423, +0.4559861124470794]
    ]);
}


/**
 * @brief Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: GAL = galactic system (IAU 1958 definition).
 * Target: EQJ = equatorial system, using the equator at the J2000 epoch.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts GAL to EQJ.
 */
export function Rotation_GAL_EQJ(): RotationMatrix {
    // This rotation matrix was calculated by the following script
    // in this same source code repository:
    // demo/python/galeqj_matrix.py
    return new RotationMatrix([
        [-0.0548624779711344, -0.8734572784246782, -0.4838000529948520],
        [+0.4941095946388765, -0.4447938112296831, +0.7470034631630423],
        [-0.8676668813529025, -0.1980677870294097, +0.4559861124470794]
    ]);
}


$ASTRO_CONSTEL()

let ConstelRot: RotationMatrix;
let Epoch2000 : AstroTime;

/**
 * @brief Reports the constellation that a given celestial point lies within.
 *
 * @property {string} symbol
 *      3-character mnemonic symbol for the constellation, e.g. "Ori".
 *
 * @property {string} name
 *      Full name of constellation, e.g. "Orion".
 *
 * @property {number} ra1875
 *      Right ascension expressed in B1875 coordinates.
 *
 * @property {number} dec1875
 *      Declination expressed in B1875 coordinates.
 */
export class ConstellationInfo {
    constructor(
        public symbol: string,
        public name: string,
        public ra1875: number,
        public dec1875: number)
        {}
}


/**
 * @brief Determines the constellation that contains the given point in the sky.
 *
 * Given J2000 equatorial (EQJ) coordinates of a point in the sky,
 * determines the constellation that contains that point.
 *
 * @param {number} ra
 *      The right ascension (RA) of a point in the sky, using the J2000 equatorial system.
 *
 * @param {number} dec
 *      The declination (DEC) of a point in the sky, using the J2000 equatorial system.
 *
 * @returns {ConstellationInfo}
 *      An object that contains the 3-letter abbreviation and full name
 *      of the constellation that contains the given (ra,dec), along with
 *      the converted B1875 (ra,dec) for that point.
 */
export function Constellation(ra: number, dec: number): ConstellationInfo {
    VerifyNumber(ra);
    VerifyNumber(dec);
    if (dec < -90 || dec > +90)
        throw 'Invalid declination angle. Must be -90..+90.';

    // Clamp right ascension to [0, 24) sidereal hours.
    ra %= 24.0;
    if (ra < 0.0)
        ra += 24.0;

    // Lazy-initialize rotation matrix.
    if (!ConstelRot) {
        // Need to calculate the B1875 epoch. Based on this:
        // https://en.wikipedia.org/wiki/Epoch_(astronomy)#Besselian_years
        // B = 1900 + (JD - 2415020.31352) / 365.242198781
        // I'm interested in using TT instead of JD, giving:
        // B = 1900 + ((TT+2451545) - 2415020.31352) / 365.242198781
        // B = 1900 + (TT + 36524.68648) / 365.242198781
        // TT = 365.242198781*(B - 1900) - 36524.68648 = -45655.741449525
        // But the AstroTime constructor wants UT, not TT.
        // Near that date, I get a historical correction of ut-tt = 3.2 seconds.
        // That gives UT = -45655.74141261017 for the B1875 epoch,
        // or 1874-12-31T18:12:21.950Z.
        ConstelRot = Rotation_EQJ_EQD(new AstroTime(-45655.74141261017));
        Epoch2000 = new AstroTime(0);
    }

    // Convert coordinates from J2000 to B1875.
    const sph2000 = new Spherical(dec, 15.0 * ra, 1.0);
    const vec2000 = VectorFromSphere(sph2000, Epoch2000);
    const vec1875 = RotateVector(ConstelRot, vec2000);
    const equ1875 = EquatorFromVector(vec1875);

    // Search for the constellation using the B1875 coordinates.
    const fd = 10 / (4 * 60);   // conversion factor from compact units to DEC degrees
    const fr = fd / 15;         // conversion factor from compact units to RA  sidereal hours
    for (let b of ConstelBounds) {
        // Convert compact angular units to RA in hours, DEC in degrees.
        const dec = b[3] * fd;
        const ra_lo = b[1] * fr;
        const ra_hi = b[2] * fr;
        if (dec <= equ1875.dec && ra_lo <= equ1875.ra && equ1875.ra < ra_hi) {
            const c = ConstelNames[b[0]];
            return new ConstellationInfo(c[0], c[1], equ1875.ra, equ1875.dec);
        }
    }

    // This should never happen!
    throw 'Unable to find constellation for given coordinates.';
}


/**
 * @brief The different kinds of lunar/solar eclipses..
 *
 * `Penumbral`: A lunar eclipse in which only the Earth's penumbra falls on the Moon. (Never used for a solar eclipse.)
 * `Partial`: A partial lunar/solar eclipse.
 * `Annular`: A solar eclipse in which the entire Moon is visible against the Sun, but the Sun appears as a ring around the Moon. (Never used for a lunar eclipse.)
 * `Total`: A total lunar/solar eclipse.
 *
 * @enum {string}
 */
export enum EclipseKind {
    Penumbral = 'penumbral',
    Partial = 'partial',
    Annular = 'annular',
    Total = 'total'
}


/**
 * @brief Returns information about a lunar eclipse.
 *
 * Returned by {@link SearchLunarEclipse} or {@link NextLunarEclipse}
 * to report information about a lunar eclipse event.
 * When a lunar eclipse is found, it is classified as penumbral, partial, or total.
 * Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
 * by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
 * Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
 * Total eclipses occur when the entire Moon passes into the Earth's umbra.
 *
 * The `kind` field thus holds one of the enum values `EclipseKind.Penumbral`, `EclipseKind.Partial`,
 * or `EclipseKind.Total`, depending on the kind of lunar eclipse found.
 *
 * Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.
 *
 * Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
 * of the eclipse, which is half of the amount of time the eclipse spends in each
 * phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
 * By converting from minutes to days, and subtracting/adding with `peak`, the caller
 * may determine the date and time of the beginning/end of each eclipse phase.
 *
 * @property {EclipseKind} kind
 *      The type of lunar eclipse found.
 *
 * @property {AstroTime} peak
 *      The time of the eclipse at its peak.
 *
 * @property {number} sd_penum
 *      The semi-duration of the penumbral phase in minutes.
 *
 * @property {number} sd_partial
 *      The semi-duration of the penumbral phase in minutes, or 0.0 if none.
 *
 * @property {number} sd_total
 *      The semi-duration of the penumbral phase in minutes, or 0.0 if none.
 *
 */
export class LunarEclipseInfo {
    constructor(
        public kind: EclipseKind,
        public peak: AstroTime,
        public sd_penum: number,
        public sd_partial: number,
        public sd_total: number)
        {}
}

/**
 * @ignore
 *
 * @brief Represents the relative alignment of the Earth and another body, and their respective shadows.
 *
 * This is an internal data structure used to assist calculation of
 * lunar eclipses, solar eclipses, and transits of Mercury and Venus.
 *
 * Definitions:
 *
 * casting body = A body that casts a shadow of interest, possibly striking another body.
 *
 * receiving body = A body on which the shadow of another body might land.
 *
 * shadow axis = The line passing through the center of the Sun and the center of the casting body.
 *
 * shadow plane = The plane passing through the center of a receiving body,
 * and perpendicular to the shadow axis.
 *
 * @property {AstroTime} time
 *      The time associated with the shadow calculation.
 *
 * @property {number} u
 *      The distance [au] between the center of the casting body and the shadow plane.
 *
 * @property {number} r
 *      The distance [km] between center of receiving body and the shadow axis.
 *
 * @property {number} k
 *      The umbra radius [km] at the shadow plane.
 *
 * @property {number} p
 *      The penumbra radius [km] at the shadow plane.
 *
 * @property {Vector} target
 *      The location in space where we are interested in determining how close a shadow falls.
 *      For example, when calculating lunar eclipses, `target` would be the center of the Moon
 *      expressed in geocentric coordinates. Then we can evaluate how far the center of the Earth's
 *      shadow cone approaches the center of the Moon.
 *      The vector components are expressed in [au].
 *
 * @property {Vector} dir
 *      The direction in space that the shadow points away from the center of a shadow-casting body.
 *      This vector lies on the shadow axis and points away from the Sun.
 *      In other words: the direction light from the Sun would be traveling,
 *      except that the center of a body (Earth, Moon, Mercury, or Venus) is blocking it.
 *      The distance units do not matter, because the vector will be normalized.
 */
class ShadowInfo {
    constructor(
        public time: AstroTime,
        public u: number,
        public r: number,
        public k: number,
        public p: number,
        public target: Vector,
        public dir: Vector)
        {}
}


function CalcShadow(body_radius_km: number, time: AstroTime, target: Vector, dir: Vector): ShadowInfo {
    const u = (dir.x*target.x + dir.y*target.y + dir.z*target.z) / (dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    const dx = (u * dir.x) - target.x;
    const dy = (u * dir.y) - target.y;
    const dz = (u * dir.z) - target.z;
    const r = KM_PER_AU * Math.hypot(dx, dy, dz);
    const k = +SUN_RADIUS_KM - (1.0 + u)*(SUN_RADIUS_KM - body_radius_km);
    const p = -SUN_RADIUS_KM + (1.0 + u)*(SUN_RADIUS_KM + body_radius_km);
    return new ShadowInfo(time, u, r, k, p, target, dir);
}


function EarthShadow(time: AstroTime): ShadowInfo {
    const e = CalcVsop(vsop.Earth, time);
    const m = GeoMoon(time);
    return CalcShadow(EARTH_ECLIPSE_RADIUS_KM, time, m, e);
}


function MoonShadow(time: AstroTime): ShadowInfo {
    // This is a variation on the logic in _EarthShadow().
    // Instead of a heliocentric Earth and a geocentric Moon,
    // we want a heliocentric Moon and a lunacentric Earth.
    const h = CalcVsop(vsop.Earth, time);    // heliocentric Earth
    const m = GeoMoon(time);       // geocentric Moon
    // Calculate lunacentric Earth.
    const e = new Vector(-m.x, -m.y, -m.z, m.t);
    // Convert geocentric moon to heliocentric Moon.
    m.x += h.x;
    m.y += h.y;
    m.z += h.z;
    return CalcShadow(MOON_MEAN_RADIUS_KM, time, e, m);
}


function LocalMoonShadow(time: AstroTime, observer: Observer): ShadowInfo {
    // Calculate observer's geocentric position.
    // For efficiency, do this first, to populate the earth rotation parameters in 'time'.
    // That way they can be recycled instead of recalculated.
    const pos = geo_pos(time, observer);
    const h = CalcVsop(vsop.Earth, time);     // heliocentric Earth
    const m = GeoMoon(time);        // geocentric Moon

    // Calculate lunacentric location of an observer on the Earth's surface.
    const o = new Vector(pos[0] - m.x, pos[1] - m.y, pos[2] - m.z, time);

    // Convert geocentric moon to heliocentric Moon.
    m.x += h.x;
    m.y += h.y;
    m.z += h.z;

    return CalcShadow(MOON_MEAN_RADIUS_KM, time, o, m);
}


function PlanetShadow(body: Body, planet_radius_km: number, time: AstroTime): ShadowInfo {
    // Calculate light-travel-corrected vector from Earth to planet.
    const g = GeoVector(body, time, false);

    // Calculate light-travel-corrected vector from Earth to Sun.
    const e = GeoVector(Body.Sun, time, false);

    // Deduce light-travel-corrected vector from Sun to planet.
    const p = new Vector(g.x - e.x, g.y - e.y, g.z - e.z, time);

    // Calcluate Earth's position from the planet's point of view.
    e.x = -g.x;
    e.y = -g.y;
    e.z = -g.z;

    return CalcShadow(planet_radius_km, time, e, p);
}


function ShadowDistanceSlope(shadowfunc: (t:AstroTime) => ShadowInfo, time: AstroTime): number {
    const dt = 1.0 / 86400.0;
    const t1 = time.AddDays(-dt);
    const t2 = time.AddDays(+dt);
    const shadow1 = shadowfunc(t1);
    const shadow2 = shadowfunc(t2);
    return (shadow2.r - shadow1.r) / dt;
}


function PlanetShadowSlope(body: Body, planet_radius_km: number, time: AstroTime): number {
    const dt = 1.0 / 86400.0;
    const shadow1 = PlanetShadow(body, planet_radius_km, time.AddDays(-dt));
    const shadow2 = PlanetShadow(body, planet_radius_km, time.AddDays(+dt));
    return (shadow2.r - shadow1.r) / dt;
}


function PeakEarthShadow(search_center_time: AstroTime): ShadowInfo {
    const window = 0.03;        /* initial search window, in days, before/after given time */
    const t1 = search_center_time.AddDays(-window);
    const t2 = search_center_time.AddDays(+window);
    const tx = Search((time: AstroTime) => ShadowDistanceSlope(EarthShadow, time), t1, t2);
    if (!tx)
        throw 'Failed to find peak Earth shadow time.';
    return EarthShadow(tx);
}


function PeakMoonShadow(search_center_time: AstroTime): ShadowInfo {
    const window = 0.03;        /* initial search window, in days, before/after given time */
    const t1 = search_center_time.AddDays(-window);
    const t2 = search_center_time.AddDays(+window);
    const tx = Search((time: AstroTime) => ShadowDistanceSlope(MoonShadow, time), t1, t2);
    if (!tx)
        throw 'Failed to find peak Moon shadow time.';
    return MoonShadow(tx);
}


function PeakPlanetShadow(body: Body, planet_radius_km: number, search_center_time: AstroTime): ShadowInfo {
    // Search for when the body's shadow is closest to the center of the Earth.
    const window = 1.0;     // days before/after inferior conjunction to search for minimum shadow distance.
    const t1 = search_center_time.AddDays(-window);
    const t2 = search_center_time.AddDays(+window);
    const tx = Search((time: AstroTime) => PlanetShadowSlope(body, planet_radius_km, time), t1, t2);
    if (!tx)
        throw 'Failed to find peak planet shadow time.';
    return PlanetShadow(body, planet_radius_km, tx);
}


function PeakLocalMoonShadow(search_center_time: AstroTime, observer: Observer): ShadowInfo {
    // Search for the time near search_center_time that the Moon's shadow comes
    // closest to the given observer.
    const window = 0.2;
    const t1 = search_center_time.AddDays(-window);
    const t2 = search_center_time.AddDays(+window);
    function shadowfunc(time: AstroTime): ShadowInfo {
        return LocalMoonShadow(time, observer);
    }
    const time = Search((time: AstroTime) => ShadowDistanceSlope(shadowfunc, time), t1, t2);
    if (!time)
        throw `PeakLocalMoonShadow: search failure for search_center_time = ${search_center_time}`;
    return LocalMoonShadow(time, observer);
}


function ShadowSemiDurationMinutes(center_time: AstroTime, radius_limit: number, window_minutes: number): number {
    // Search backwards and forwards from the center time until shadow axis distance crosses radius limit.
    const window = window_minutes / (24.0 * 60.0);
    const before = center_time.AddDays(-window);
    const after  = center_time.AddDays(+window);
    const t1 = Search((time: AstroTime) => -(EarthShadow(time).r - radius_limit), before, center_time);
    const t2 = Search((time: AstroTime) => +(EarthShadow(time).r - radius_limit), center_time, after);
    if (!t1 || !t2)
        throw 'Failed to find shadow semiduration';
    return (t2.ut - t1.ut) * ((24.0 * 60.0) / 2.0);    // convert days to minutes and average the semi-durations.
}


function MoonEclipticLatitudeDegrees(time: AstroTime): number {
    const moon = CalcMoon(time);
    return RAD2DEG * moon.geo_eclip_lat;
}


/**
 * @brief Searches for a lunar eclipse.
 *
 * This function finds the first lunar eclipse that occurs after `startTime`.
 * A lunar eclipse may be penumbral, partial, or total.
 * See {@link LunarEclipseInfo} for more information.
 * To find a series of lunar eclipses, call this function once,
 * then keep calling {@link NextLunarEclipse} as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for starting the search for a lunar eclipse.
 *
 * @returns {LunarEclipseInfo}
 */
export function SearchLunarEclipse(date: FlexibleDateTime): LunarEclipseInfo {
    const PruneLatitude = 1.8;   /* full Moon's ecliptic latitude above which eclipse is impossible */
    let fmtime = MakeTime(date);
    for (let fmcount = 0; fmcount < 12; ++fmcount) {
        /* Search for the next full moon. Any eclipse will be near it. */
        const fullmoon = SearchMoonPhase(180, fmtime, 40);
        if (!fullmoon)
            throw 'Cannot find full moon.';

        /*
            Pruning: if the full Moon's ecliptic latitude is too large,
            a lunar eclipse is not possible. Avoid needless work searching for
            the minimum moon distance.
        */
       const eclip_lat = MoonEclipticLatitudeDegrees(fullmoon);
       if (Math.abs(eclip_lat) < PruneLatitude) {
           /* Search near the full moon for the time when the center of the Moon */
           /* is closest to the line passing through the centers of the Sun and Earth. */
           const shadow = PeakEarthShadow(fullmoon);
           if (shadow.r < shadow.p + MOON_MEAN_RADIUS_KM) {
               /* This is at least a penumbral eclipse. We will return a result. */
               let kind = EclipseKind.Penumbral;
               let sd_total = 0.0;
               let sd_partial = 0.0;
               let sd_penum = ShadowSemiDurationMinutes(shadow.time, shadow.p + MOON_MEAN_RADIUS_KM, 200.0);

               if (shadow.r < shadow.k + MOON_MEAN_RADIUS_KM) {
                   /* This is at least a partial eclipse. */
                   kind = EclipseKind.Partial;
                   sd_partial = ShadowSemiDurationMinutes(shadow.time, shadow.k + MOON_MEAN_RADIUS_KM, sd_penum);

                   if (shadow.r + MOON_MEAN_RADIUS_KM < shadow.k) {
                       /* This is a total eclipse. */
                       kind = EclipseKind.Total;
                       sd_total = ShadowSemiDurationMinutes(shadow.time, shadow.k - MOON_MEAN_RADIUS_KM, sd_partial);
                   }
               }
               return new LunarEclipseInfo(kind, shadow.time, sd_penum, sd_partial, sd_total);
           }
       }

       /* We didn't find an eclipse on this full moon, so search for the next one. */
       fmtime = fullmoon.AddDays(10);
    }

    /* This should never happen because there are always at least 2 full moons per year. */
    throw 'Failed to find lunar eclipse within 12 full moons.';
}


/**
 * @brief Reports the time and geographic location of the peak of a solar eclipse.
 *
 * Returned by {@link SearchGlobalSolarEclipse} or {@link NextGlobalSolarEclipse}
 * to report information about a solar eclipse event.
 *
 * The eclipse is classified as partial, annular, or total, depending on the
 * maximum amount of the Sun's disc obscured, as seen at the peak location
 * on the surface of the Earth.
 *
 * The `kind` field thus holds one of the values `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
 * A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
 * An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
 * to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
 * A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
 * observes either a total or annular eclipse.
 *
 * If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
 * fields give the geographic coordinates of the center of the Moon's shadow projected
 * onto the daytime side of the Earth at the instant of the eclipse's peak.
 * If `kind` has any other value, `latitude` and `longitude` are undefined and should
 * not be used.
 *
 * @property {EclipseKind} kind
 *     One of the following enumeration values: `EclipseKind.Partial`, `EclipseKind.Annular`, `EclipseKind.Total`.
 *
 * @property {AstroTime} peak
 *     The date and time when the solar eclipse is darkest.
 *     This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.
 *
 * @property {number} distance
 *     The distance in kilometers between the axis of the Moon's shadow cone
 *     and the center of the Earth at the time indicated by `peak`.
 *
 * @property {number | undefined} latitude
 *     If `kind` holds `EclipseKind.Total`, the geographic latitude in degrees
 *     where the center of the Moon's shadow falls on the Earth at the
 *     time indicated by `peak`; otherwise, `latitude` holds `undefined`.
 *
 * @property {number | undefined} longitude
 *     If `kind` holds `EclipseKind.Total`, the geographic longitude in degrees
 *     where the center of the Moon's shadow falls on the Earth at the
 *     time indicated by `peak`; otherwise, `longitude` holds `undefined`.
 */
export class GlobalSolarEclipseInfo {
    constructor(
        public kind: EclipseKind,
        public peak: AstroTime,
        public distance: number,
        public latitude?: number,
        public longitude?: number)
        {}
}


function EclipseKindFromUmbra(k: number): EclipseKind {
    // The umbra radius tells us what kind of eclipse the observer sees.
    // If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular.
    // HACK: I added a tiny bias (14 meters) to match Espenak test data.
    return (k > 0.014) ? EclipseKind.Total : EclipseKind.Annular;
}


function GeoidIntersect(shadow: ShadowInfo): GlobalSolarEclipseInfo {
    let kind = EclipseKind.Partial;
    let peak = shadow.time;
    let distance = shadow.r;
    let latitude: number | undefined;       // left undefined for partial eclipses
    let longitude: number | undefined;      // left undefined for partial eclipses

    // We want to calculate the intersection of the shadow axis with the Earth's geoid.
    // First we must convert EQJ (equator of J2000) coordinates to EQD (equator of date)
    // coordinates that are perfectly aligned with the Earth's equator at this
    // moment in time.
    const rot = Rotation_EQJ_EQD(shadow.time);
    const v = RotateVector(rot, shadow.dir);       // shadow-axis vector in equator-of-date coordinates
    const e = RotateVector(rot, shadow.target);    // lunacentric Earth in equator-of-date coordinates

    // Convert all distances from AU to km.
    // But dilate the z-coordinates so that the Earth becomes a perfect sphere.
    // Then find the intersection of the vector with the sphere.
    // See p 184 in Montenbruck & Pfleger's "Astronomy on the Personal Computer", second edition.
    v.x *= KM_PER_AU;
    v.y *= KM_PER_AU;
    v.z *= KM_PER_AU / EARTH_FLATTENING;
    e.x *= KM_PER_AU;
    e.y *= KM_PER_AU;
    e.z *= KM_PER_AU / EARTH_FLATTENING;

    // Solve the quadratic equation that finds whether and where
    // the shadow axis intersects with the Earth in the dilated coordinate system.
    const R = EARTH_EQUATORIAL_RADIUS_KM;
    const A = v.x*v.x + v.y*v.y + v.z*v.z;
    const B = -2.0 * (v.x*e.x + v.y*e.y + v.z*e.z);
    const C = (e.x*e.x + e.y*e.y + e.z*e.z) - R*R;
    const radic = B*B - 4*A*C;

    if (radic > 0.0) {
        // Calculate the closer of the two intersection points.
        // This will be on the day side of the Earth.
        const u = (-B - Math.sqrt(radic)) / (2 * A);

        // Convert lunacentric dilated coordinates to geocentric coordinates.
        const px = u*v.x - e.x;
        const py = u*v.y - e.y;
        const pz = (u*v.z - e.z) * EARTH_FLATTENING;

        // Convert cartesian coordinates into geodetic latitude/longitude.
        const proj = Math.hypot(px, py) * EARTH_FLATTENING_SQUARED;
        if (proj == 0.0)
            latitude = (pz > 0.0) ? +90.0 : -90.0;
        else
            latitude = RAD2DEG * Math.atan(pz / proj);

        // Adjust longitude for Earth's rotation at the given UT.
        const gast = sidereal_time(peak);
        longitude = (RAD2DEG*Math.atan2(py, px) - (15*gast)) % 360.0;
        if (longitude <= -180.0)
            longitude += 360.0;
        else if (longitude > +180.0)
            longitude -= 360.0;

        // We want to determine whether the observer sees a total eclipse or an annular eclipse.
        // We need to perform a series of vector calculations...
        // Calculate the inverse rotation matrix, so we can convert EQD to EQJ.
        const inv = InverseRotation(rot);

        // Put the EQD geocentric coordinates of the observer into the vector 'o'.
        // Also convert back from kilometers to astronomical units.
        let o = new Vector(px / KM_PER_AU, py / KM_PER_AU, pz / KM_PER_AU, shadow.time);

        // Rotate the observer's geocentric EQD back to the EQJ system.
        o = RotateVector(inv, o);

        // Convert geocentric vector to lunacentric vector.
        o.x += shadow.target.x;
        o.y += shadow.target.y;
        o.z += shadow.target.z;

        // Recalculate the shadow using a vector from the Moon's center toward the observer.
        const surface = CalcShadow(MOON_POLAR_RADIUS_KM, shadow.time, o, shadow.dir);

        // If we did everything right, the shadow distance should be very close to zero.
        // That's because we already determined the observer 'o' is on the shadow axis!
        if (surface.r > 1.0e-9 || surface.r < 0.0)
            throw `Unexpected shadow distance from geoid intersection = ${surface.r}`;

        kind = EclipseKindFromUmbra(surface.k);
    }

    return new GlobalSolarEclipseInfo(kind, peak, distance, latitude, longitude);
}


/**
 * @brief Searches for the next lunar eclipse in a series.
 *
 * After using {@link SearchLunarEclipse} to find the first lunar eclipse
 * in a series, you can call this function to find the next consecutive lunar eclipse.
 * Pass in the `peak` value from the {@link LunarEclipseInfo} returned by the
 * previous call to `SearchLunarEclipse` or `NextLunarEclipse`
 * to find the next lunar eclipse.
 *
 * @param {FlexibleDateTime} prevEclipseTime
 *      A date and time near a full moon. Lunar eclipse search will start at the next full moon.
 *
 * @returns {LunarEclipseInfo}
 */
export function NextLunarEclipse(prevEclipseTime: FlexibleDateTime): LunarEclipseInfo {
    prevEclipseTime = MakeTime(prevEclipseTime);
    const startTime = prevEclipseTime.AddDays(10);
    return SearchLunarEclipse(startTime);
}

/**
 * @brief Searches for a solar eclipse visible anywhere on the Earth's surface.
 *
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See {@link GlobalSolarEclipseInfo} for more information.
 * To find a series of solar eclipses, call this function once,
 * then keep calling {@link NextGlobalSolarEclipse} as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * @param {FlexibleDateTime} startTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @returns {GlobalSolarEclipseInfo}
 */
export function SearchGlobalSolarEclipse(startTime: FlexibleDateTime): GlobalSolarEclipseInfo {
    startTime = MakeTime(startTime);
    const PruneLatitude = 1.8;      // Moon's ecliptic latitude beyond which eclipse is impossible
    // Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
    let nmtime = startTime;
    let nmcount: number;
    for (nmcount=0; nmcount < 12; ++nmcount) {
        // Search for the next new moon. Any eclipse will be near it.
        const newmoon = SearchMoonPhase(0.0, nmtime, 40.0);
        if (!newmoon)
            throw 'Cannot find new moon';

        // Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
        const eclip_lat = MoonEclipticLatitudeDegrees(newmoon);
        if (Math.abs(eclip_lat) < PruneLatitude) {
            // Search near the new moon for the time when the center of the Earth
            // is closest to the line passing through the centers of the Sun and Moon.
            const shadow = PeakMoonShadow(newmoon);
            if (shadow.r < shadow.p + EARTH_MEAN_RADIUS_KM) {
                // This is at least a partial solar eclipse visible somewhere on Earth.
                // Try to find an intersection between the shadow axis and the Earth's oblate geoid.
                return GeoidIntersect(shadow);
            }
        }

        // We didn't find an eclipse on this new moon, so search for the next one.
        nmtime = newmoon.AddDays(10.0);
    }

    // Safety valve to prevent infinite loop.
    // This should never happen, because at least 2 solar eclipses happen per year.
    throw 'Failed to find solar eclipse within 12 full moons.';
}


/**
 * @brief Searches for the next global solar eclipse in a series.
 *
 * After using {@link SearchGlobalSolarEclipse} to find the first solar eclipse
 * in a series, you can call this function to find the next consecutive solar eclipse.
 * Pass in the `peak` value from the {@link GlobalSolarEclipseInfo} returned by the
 * previous call to `SearchGlobalSolarEclipse` or `NextGlobalSolarEclipse`
 * to find the next solar eclipse.
 *
 * @param {FlexibleDateTime} prevEclipseTime
 *      A date and time near a new moon. Solar eclipse search will start at the next new moon.
 *
 * @returns {GlobalSolarEclipseInfo}
 */
export function NextGlobalSolarEclipse(prevEclipseTime: FlexibleDateTime): GlobalSolarEclipseInfo {
    prevEclipseTime = MakeTime(prevEclipseTime);
    const startTime = prevEclipseTime.AddDays(10.0);
    return SearchGlobalSolarEclipse(startTime);
}


/**
 * @brief Holds a time and the observed altitude of the Sun at that time.
 *
 * When reporting a solar eclipse observed at a specific location on the Earth
 * (a "local" solar eclipse), a series of events occur. In addition
 * to the time of each event, it is important to know the altitude of the Sun,
 * because each event may be invisible to the observer if the Sun is below
 * the horizon.
 *
 * If `altitude` is negative, the event is theoretical only; it would be
 * visible if the Earth were transparent, but the observer cannot actually see it.
 * If `altitude` is positive but less than a few degrees, visibility will be impaired by
 * atmospheric interference (sunrise or sunset conditions).
 *
 * @property {AstroTime} time
 *      The date and time of the event.
 *
 * @property {number} altitude
 *      The angular altitude of the center of the Sun above/below the horizon, at `time`,
 *      corrected for atmospheric refraction and expressed in degrees.
 */
export class EclipseEvent {
    constructor(
        public time: AstroTime,
        public altitude: number)
        {}
}


/**
 * @brief Information about a solar eclipse as seen by an observer at a given time and geographic location.
 *
 * Returned by {@link SearchLocalSolarEclipse} or {@link NextLocalSolarEclipse}
 * to report information about a solar eclipse as seen at a given geographic location.
 *
 * When a solar eclipse is found, it is classified by setting `kind`
 * to `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
 * A partial solar eclipse is when the Moon does not line up directly enough with the Sun
 * to completely block the Sun's light from reaching the observer.
 * An annular eclipse occurs when the Moon's disc is completely visible against the Sun
 * but the Moon is too far away to completely block the Sun's light; this leaves the
 * Sun with a ring-like appearance.
 * A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
 * Sun just right to completely block all sunlight from reaching the observer.
 *
 * There are 5 "event" fields, each of which contains a time and a solar altitude.
 * Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
 * The fields `partial_begin` and `partial_end` are always set, and indicate when
 * the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
 * `total_begin` and `total_end` indicate when the total/annular phase begins/ends.
 * When an event field is valid, the caller must also check its `altitude` field to
 * see whether the Sun is above the horizon at the time indicated by the `time` field.
 * See {@link EclipseEvent} for more information.
 *
 * @property {EclipseKind} kind
 *      The type of solar eclipse found: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
 *
 * @property {EclipseEvent} partial_begin
 *      The time and Sun altitude at the beginning of the eclipse.
 *
 * @property {EclipseEvent | undefined} total_begin
 *      If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise undefined.
 *
 * @property {EclipseEvent} peak
 *      The time and Sun altitude when the eclipse reaches its peak.
 *
 * @property {EclipseEvent | undefined} total_end
 *      If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise undefined.
 *
 * @property {EclipseEvent} partial_end
 *      The time and Sun altitude at the end of the eclipse.
 */
export class LocalSolarEclipseInfo {
    constructor(
        public kind: EclipseKind,
        public partial_begin: EclipseEvent,
        public total_begin: EclipseEvent | undefined,
        public peak : EclipseEvent,
        public total_end: EclipseEvent | undefined,
        public partial_end: EclipseEvent)
        {}
}


function local_partial_distance(shadow: ShadowInfo): number {
    return shadow.p - shadow.r;
}

function local_total_distance(shadow: ShadowInfo): number {
    // Must take the absolute value of the umbra radius 'k'
    // because it can be negative for an annular eclipse.
    return Math.abs(shadow.k) - shadow.r;
}


function LocalEclipse(shadow: ShadowInfo, observer: Observer): LocalSolarEclipseInfo {
    const PARTIAL_WINDOW = 0.2;
    const TOTAL_WINDOW = 0.01;
    const peak = CalcEvent(observer, shadow.time);
    let t1 = shadow.time.AddDays(-PARTIAL_WINDOW);
    let t2 = shadow.time.AddDays(+PARTIAL_WINDOW);
    const partial_begin = LocalEclipseTransition(observer, +1.0, local_partial_distance, t1, shadow.time);
    const partial_end   = LocalEclipseTransition(observer, -1.0, local_partial_distance, shadow.time, t2);
    let total_begin: EclipseEvent | undefined;
    let total_end: EclipseEvent | undefined;
    let kind: EclipseKind;

    if (shadow.r < Math.abs(shadow.k)) {     // take absolute value of 'k' to handle annular eclipses too.
        t1 = shadow.time.AddDays(-TOTAL_WINDOW);
        t2 = shadow.time.AddDays(+TOTAL_WINDOW);
        total_begin = LocalEclipseTransition(observer, +1.0, local_total_distance, t1, shadow.time);
        total_end = LocalEclipseTransition(observer, -1.0, local_total_distance, shadow.time, t2);
        kind = EclipseKindFromUmbra(shadow.k);
    } else {
        kind = EclipseKind.Partial;
    }

    return new LocalSolarEclipseInfo(kind, partial_begin, total_begin, peak, total_end, partial_end);
}

type ShadowFunc = (shadow: ShadowInfo) => number;

function LocalEclipseTransition(observer: Observer, direction: number, func: ShadowFunc, t1: AstroTime, t2: AstroTime): EclipseEvent {
    function evaluate(time: AstroTime): number {
        const shadow = LocalMoonShadow(time, observer);
        return direction * func(shadow);
    }
    const search = Search(evaluate, t1, t2);
    if (!search)
        throw "Local eclipse transition search failed.";
    return CalcEvent(observer, search);
}

function CalcEvent(observer: Observer, time: AstroTime): EclipseEvent {
    const altitude = SunAltitude(time, observer);
    return new EclipseEvent(time, altitude);
}

function SunAltitude(time: AstroTime, observer: Observer): number {
    const equ = Equator(Body.Sun, time, observer, true, true);
    const hor = Horizon(time, observer, equ.ra, equ.dec, 'normal');
    return hor.altitude;
}


/**
 * @brief Searches for a solar eclipse visible at a specific location on the Earth's surface.
 *
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See {@link LocalSolarEclipseInfo} for more information.
 *
 * To find a series of solar eclipses, call this function once,
 * then keep calling {@link NextLocalSolarEclipse} as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * IMPORTANT: An eclipse reported by this function might be partly or
 * completely invisible to the observer due to the time of day.
 * See {@link LocalSolarEclipseInfo} for more information about this topic.
 *
 * @param {FlexibleDateTime} startTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @param {Observer} observer
 *      The geographic location of the observer.
 *
 * @returns {LocalSolarEclipseInfo}
 */
export function SearchLocalSolarEclipse(startTime: FlexibleDateTime, observer: Observer): LocalSolarEclipseInfo {
    startTime = MakeTime(startTime);
    VerifyObserver(observer);
    const PruneLatitude = 1.8;   /* Moon's ecliptic latitude beyond which eclipse is impossible */

    /* Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth. */
    let nmtime = startTime;
    for(;;) {
        /* Search for the next new moon. Any eclipse will be near it. */
        const newmoon = SearchMoonPhase(0.0, nmtime, 40.0);
        if (!newmoon)
            throw 'Cannot find next new moon';

        /* Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible. */
        const eclip_lat = MoonEclipticLatitudeDegrees(newmoon);
        if (Math.abs(eclip_lat) < PruneLatitude) {
            /* Search near the new moon for the time when the observer */
            /* is closest to the line passing through the centers of the Sun and Moon. */
            const shadow = PeakLocalMoonShadow(newmoon, observer);
            if (shadow.r < shadow.p) {
                /* This is at least a partial solar eclipse for the observer. */
                const eclipse = LocalEclipse(shadow, observer);

                /* Ignore any eclipse that happens completely at night. */
                /* More precisely, the center of the Sun must be above the horizon */
                /* at the beginning or the end of the eclipse, or we skip the event. */
                if (eclipse.partial_begin.altitude > 0.0 || eclipse.partial_end.altitude > 0.0)
                    return eclipse;
            }
        }

        /* We didn't find an eclipse on this new moon, so search for the next one. */
        nmtime = newmoon.AddDays(10.0);
    }
}


/**
 * @brief Searches for the next local solar eclipse in a series.
 *
 * After using {@link SearchLocalSolarEclipse} to find the first solar eclipse
 * in a series, you can call this function to find the next consecutive solar eclipse.
 * Pass in the `peak` value from the {@link LocalSolarEclipseInfo} returned by the
 * previous call to `SearchLocalSolarEclipse` or `NextLocalSolarEclipse`
 * to find the next solar eclipse.
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See {@link LocalSolarEclipseInfo} for more information.
 *
 * @param {FlexibleDateTime} prevEclipseTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @param {Observer} observer
 *      The geographic location of the observer.
 *
 * @returns {LocalSolarEclipseInfo}
 */
export function NextLocalSolarEclipse(prevEclipseTime: FlexibleDateTime, observer: Observer): LocalSolarEclipseInfo {
    prevEclipseTime = MakeTime(prevEclipseTime);
    const startTime = prevEclipseTime.AddDays(10.0);
    return SearchLocalSolarEclipse(startTime, observer);
}


/**
 * @brief Information about a transit of Mercury or Venus, as seen from the Earth.
 *
 * Returned by {@link SearchTransit} or {@link NextTransit} to report
 * information about a transit of Mercury or Venus.
 * A transit is when Mercury or Venus passes between the Sun and Earth so that
 * the other planet is seen in silhouette against the Sun.
 *
 * The calculations are performed from the point of view of a geocentric observer.
 *
 * @property {AstroTime} start
 *      The date and time at the beginning of the transit.
 *      This is the moment the planet first becomes visible against the Sun in its background.
 *
 * @property {AstroTime} peak
 *      When the planet is most aligned with the Sun, as seen from the Earth.
 *
 * @property {AstroTime} finish
 *      The date and time at the end of the transit.
 *      This is the moment the planet is last seen against the Sun in its background.
 *
 * @property {number} separation
 *      The minimum angular separation, in arcminutes, between the centers of the Sun and the planet.
 *      This angle pertains to the time stored in `peak`.
 */
export class TransitInfo {
    constructor(
        public start: AstroTime,
        public peak: AstroTime,
        public finish: AstroTime,
        public separation: number)
        {}
}


function PlanetShadowBoundary(time: AstroTime, body: Body, planet_radius_km: number, direction: number): number {
    const shadow = PlanetShadow(body, planet_radius_km, time);
    return direction * (shadow.r - shadow.p);
}


function PlanetTransitBoundary(body: Body, planet_radius_km: number, t1: AstroTime, t2: AstroTime, direction: number): AstroTime {
    // Search for the time the planet's penumbra begins/ends making contact with the center of the Earth.
    const tx = Search((time: AstroTime) => PlanetShadowBoundary(time, body, planet_radius_km, direction), t1, t2);
    if (!tx)
        throw 'Planet transit boundary search failed';

    return tx;
}


/**
 * @brief Searches for the first transit of Mercury or Venus after a given date.
 *
 * Finds the first transit of Mercury or Venus after a specified date.
 * A transit is when an inferior planet passes between the Sun and the Earth
 * so that the silhouette of the planet is visible against the Sun in the background.
 * To continue the search, pass the `finish` time in the returned structure to
 * {@link NextTransit}.
 *
 * @param {Body} body
 *      The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
 *
 * @param {FlexibleDateTime} startTime
 *      The date and time for starting the search for a transit.
 *
 * @returns {TransitInfo}
 */
export function SearchTransit(body: Body, startTime: FlexibleDateTime): TransitInfo {
    startTime = MakeTime(startTime);
    const threshold_angle = 0.4;     // maximum angular separation to attempt transit calculation
    const dt_days = 1.0;

    // Validate the planet and find its mean radius.
    let planet_radius_km: number;
    switch (body)
    {
        case Body.Mercury:
            planet_radius_km = 2439.7;
            break;

        case Body.Venus:
            planet_radius_km = 6051.8;
            break;

        default:
            throw `Invalid body: ${body}`;
    }

    let search_time = startTime;
    for(;;) {
        // Search for the next inferior conjunction of the given planet.
        // This is the next time the Earth and the other planet have the same
        // ecliptic longitude as seen from the Sun.
        const conj = SearchRelativeLongitude(body, 0.0, search_time);

        // Calculate the angular separation between the body and the Sun at this time.
        const conj_separation = AngleFromSun(body, conj);

        if (conj_separation < threshold_angle) {
            // The planet's angular separation from the Sun is small enough
            // to consider it a transit candidate.
            // Search for the moment when the line passing through the Sun
            // and planet are closest to the Earth's center.
            const shadow = PeakPlanetShadow(body, planet_radius_km, conj);

            if (shadow.r < shadow.p) {      // does the planet's penumbra touch the Earth's center?
                // Find the beginning and end of the penumbral contact.
                const time_before = shadow.time.AddDays(-dt_days);
                const start = PlanetTransitBoundary(body, planet_radius_km, time_before, shadow.time, -1.0);
                const time_after = shadow.time.AddDays(+dt_days);
                const finish = PlanetTransitBoundary(body, planet_radius_km, shadow.time, time_after, +1.0);
                const min_separation = 60.0 * AngleFromSun(body, shadow.time);
                return new TransitInfo(start, shadow.time, finish, min_separation);
            }
        }

        // This inferior conjunction was not a transit. Try the next inferior conjunction.
        search_time = conj.AddDays(10.0);
    }
}


/**
 * @brief Searches for the next transit of Mercury or Venus in a series.
 *
 * After calling {@link SearchTransit} to find a transit of Mercury or Venus,
 * this function finds the next transit after that.
 * Keep calling this function as many times as you want to keep finding more transits.
 *
 * @param {Body} body
 *      The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
 *
 * @param {FlexibleDateTime} prevTransitTime
 *      A date and time near the previous transit.
 *
 * @returns {TransitInfo}
 */
export function NextTransit(body: Body, prevTransitTime: FlexibleDateTime): TransitInfo {
    prevTransitTime = MakeTime(prevTransitTime);
    const startTime = prevTransitTime.AddDays(100.0);
    return SearchTransit(body, startTime);
}


/**
 * @brief Indicates whether a crossing through the ecliptic plane is ascending or descending.
 *
 * `Invalid` is a placeholder for an unknown or missing node.
 * `Ascending` indicates a body passing through the ecliptic plane from south to north.
 * `Descending` indicates a body passing through the ecliptic plane from north to south.
 *
 * @enum {number}
 */
export enum NodeEventKind {
    Invalid = 0,
    Ascending = +1,
    Descending = -1
}

/**
 * @brief Information about an ascending or descending node of a body.
 *
 * This object is returned by {@link SearchMoonNode} and {@link NextMoonNode}
 * to report information about the center of the Moon passing through the ecliptic plane.
 *
 * @property {NodeEventKind} kind   Whether the node is ascending (south to north) or descending (north to south).
 * @property {AstroTime} time       The time when the body passes through the ecliptic plane.
 */
export class NodeEventInfo {
    constructor(
        public kind: NodeEventKind,
        public time: AstroTime
    ) {}
}

const MoonNodeStepDays = +10.0;   // a safe number of days to step without missing a Moon node

/**
 * @brief Searches for a time when the Moon's center crosses through the ecliptic plane.
 *
 * Searches for the first ascending or descending node of the Moon after `startTime`.
 * An ascending node is when the Moon's center passes through the ecliptic plane
 * (the plane of the Earth's orbit around the Sun) from south to north.
 * A descending node is when the Moon's center passes through the ecliptic plane
 * from north to south. Nodes indicate possible times of solar or lunar eclipses,
 * if the Moon also happens to be in the correct phase (new or full, respectively).
 * Call `SearchMoonNode` to find the first of a series of nodes.
 * Then call {@link NextMoonNode} to find as many more consecutive nodes as desired.
 *
 * @param {FlexibleDateTime} startTime
 *      The date and time for starting the search for an ascending or descending node of the Moon.
 *
 * @returns {NodeEventInfo}
 */
export function SearchMoonNode(startTime: FlexibleDateTime): NodeEventInfo {
    // Start at the given moment in time and sample the Moon's ecliptic latitude.
    // Step 10 days at a time, searching for an interval where that latitude crosses zero.

    let time1: AstroTime = MakeTime(startTime);
    let eclip1: Spherical = EclipticGeoMoon(time1);

    for(;;) {
        const time2: AstroTime = time1.AddDays(MoonNodeStepDays);
        const eclip2: Spherical = EclipticGeoMoon(time2);
        if (eclip1.lat * eclip2.lat <= 0.0) {
            // There is a node somewhere inside this closed time interval.
            // Figure out whether it is an ascending node or a descending node.
            const kind = (eclip2.lat > eclip1.lat) ? NodeEventKind.Ascending : NodeEventKind.Descending;
            const result = Search(t => kind * EclipticGeoMoon(t).lat, time1, time2);
            if (!result)
                throw `Could not find moon node.`;   // should never happen
            return new NodeEventInfo(kind, result);
        }
        time1 = time2;
        eclip1 = eclip2;
    }
}

/**
 * @brief Searches for the next time when the Moon's center crosses through the ecliptic plane.
 *
 * Call {@link SearchMoonNode} to find the first of a series of nodes.
 * Then call `NextMoonNode` to find as many more consecutive nodes as desired.
 *
 * @param {NodeEventInfo} prevNode
 *      The previous node found from calling {@link SearchMoonNode} or `NextMoonNode`.
 *
 * @returns {NodeEventInfo}
 */
export function NextMoonNode(prevNode: NodeEventInfo): NodeEventInfo {
    const time = prevNode.time.AddDays(MoonNodeStepDays);
    const node = SearchMoonNode(time);
    switch (prevNode.kind) {
        case NodeEventKind.Ascending:
            if (node.kind !== NodeEventKind.Descending)
                throw `Internal error: previous node was ascending, but this node was: ${node.kind}`;
            break;

        case NodeEventKind.Descending:
            if (node.kind !== NodeEventKind.Ascending)
                throw `Internal error: previous node was descending, but this node was: ${node.kind}`;
            break;

        default:
            throw `Previous node has an invalid node kind: ${prevNode.kind}`;
    }
    return node;
}


/**
 * @brief Information about a body's rotation axis at a given time.
 *
 * This structure is returned by {@link RotationAxis} to report
 * the orientation of a body's rotation axis at a given moment in time.
 * The axis is specified by the direction in space that the body's north pole
 * points, using angular equatorial coordinates in the J2000 system (EQJ).
 *
 * Thus `ra` is the right ascension, and `dec` is the declination, of the
 * body's north pole vector at the given moment in time. The north pole
 * of a body is defined as the pole that lies on the north side of the
 * [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
 * regardless of the body's direction of rotation.
 *
 * The `spin` field indicates the angular position of a prime meridian
 * arbitrarily recommended for the body by the International Astronomical
 * Union (IAU).
 *
 * The fields `ra`, `dec`, and `spin` correspond to the variables
 * α0, δ0, and W, respectively, from
 * [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
 * The field `north` is a unit vector pointing in the direction of the body's north pole.
 * It is expressed in the equatorial J2000 system (EQJ).
 *
 * @property {number} ra
 *      The J2000 right ascension of the body's north pole direction, in sidereal hours.
 *
 * @property {number} dec
 *      The J2000 declination of the body's north pole direction, in degrees.
 *
 * @property {number} spin
 *      Rotation angle of the body's prime meridian, in degrees.
 *
 * @property {Vector} north
 *      A J2000 dimensionless unit vector pointing in the direction of the body's north pole.
 */
export class AxisInfo {
    constructor(
        public ra: number,
        public dec: number,
        public spin: number,
        public north: Vector)
        {}
}

function EarthRotationAxis(time: AstroTime): AxisInfo {
    // Unlike the other planets, we have a model of precession and nutation
    // for the Earth's axis that provides a north pole vector.
    // So calculate the vector first, then derive the (RA,DEC) angles from the vector.

    // Start with a north pole vector in equator-of-date coordinates: (0,0,1).
    // Convert the vector into J2000 coordinates.
    const pos2 = nutation([0, 0, 1], time, PrecessDirection.Into2000);
    const nvec = precession(pos2, time, PrecessDirection.Into2000);
    const north = new Vector(nvec[0], nvec[1], nvec[2], time);

    // Derive angular values: right ascension and declination.
    const equ = EquatorFromVector(north);

    // Use a modified version of the era() function that does not trim to 0..360 degrees.
    // This expression is also corrected to give the correct angle at the J2000 epoch.
    const spin = 190.41375788700253 + (360.9856122880876 * time.ut);

    return new AxisInfo(equ.ra, equ.dec, spin, north);
}

/**
 * @brief Calculates information about a body's rotation axis at a given time.
 * Calculates the orientation of a body's rotation axis, along with
 * the rotation angle of its prime meridian, at a given moment in time.
 *
 * This function uses formulas standardized by the IAU Working Group
 * on Cartographics and Rotational Elements 2015 report, as described
 * in the following document:
 *
 * https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf
 *
 * See {@link AxisInfo} for more detailed information.
 *
 * @param {Body} body
 *      One of the following values:
 *      `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`,
 *      `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
 *
 * @param {FlexibleDateTime} date
 *      The time at which to calculate the body's rotation axis.
 *
 * @returns {AxisInfo}
 */
export function RotationAxis(body: Body, date: FlexibleDateTime): AxisInfo {
    const time = MakeTime(date);
    const d = time.tt;
    const T = d / 36525.0;
    let ra:number, dec:number, w:number;

    switch (body)
    {
    case Body.Sun:
        ra = 286.13;
        dec = 63.87;
        w = 84.176 + (14.1844 * d);
        break;

    case Body.Mercury:
        ra = 281.0103 - (0.0328 * T);
        dec = 61.4155 - (0.0049 * T);
        w = (
            329.5988
            + (6.1385108 * d)
            + (0.01067257 * Math.sin(DEG2RAD*(174.7910857 + 4.092335*d)))
            - (0.00112309 * Math.sin(DEG2RAD*(349.5821714 + 8.184670*d)))
            - (0.00011040 * Math.sin(DEG2RAD*(164.3732571 + 12.277005*d)))
            - (0.00002539 * Math.sin(DEG2RAD*(339.1643429 + 16.369340*d)))
            - (0.00000571 * Math.sin(DEG2RAD*(153.9554286 + 20.461675*d)))
        );
        break;

    case Body.Venus:
        ra = 272.76;
        dec = 67.16;
        w = 160.20 - (1.4813688 * d);
        break;

    case Body.Earth:
        return EarthRotationAxis(time);

    case Body.Moon:
        // See page 8, Table 2 in:
        // https://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf
        const E1  = DEG2RAD * (125.045 -  0.0529921*d);
        const E2  = DEG2RAD * (250.089 -  0.1059842*d);
        const E3  = DEG2RAD * (260.008 + 13.0120009*d);
        const E4  = DEG2RAD * (176.625 + 13.3407154*d);
        const E5  = DEG2RAD * (357.529 +  0.9856003*d);
        const E6  = DEG2RAD * (311.589 + 26.4057084*d);
        const E7  = DEG2RAD * (134.963 + 13.0649930*d);
        const E8  = DEG2RAD * (276.617 +  0.3287146*d);
        const E9  = DEG2RAD * (34.226  +  1.7484877*d);
        const E10 = DEG2RAD * (15.134  -  0.1589763*d);
        const E11 = DEG2RAD * (119.743 +  0.0036096*d);
        const E12 = DEG2RAD * (239.961 +  0.1643573*d);
        const E13 = DEG2RAD * (25.053  + 12.9590088*d);

        ra = (
            269.9949 + 0.0031*T
            - 3.8787*Math.sin(E1)
            - 0.1204*Math.sin(E2)
            + 0.0700*Math.sin(E3)
            - 0.0172*Math.sin(E4)
            + 0.0072*Math.sin(E6)
            - 0.0052*Math.sin(E10)
            + 0.0043*Math.sin(E13)
        );

        dec = (
            66.5392 + 0.0130*T
            + 1.5419*Math.cos(E1)
            + 0.0239*Math.cos(E2)
            - 0.0278*Math.cos(E3)
            + 0.0068*Math.cos(E4)
            - 0.0029*Math.cos(E6)
            + 0.0009*Math.cos(E7)
            + 0.0008*Math.cos(E10)
            - 0.0009*Math.cos(E13)
        );

        w = (
            38.3213 + (13.17635815 - 1.4e-12*d)*d
            + 3.5610*Math.sin(E1)
            + 0.1208*Math.sin(E2)
            - 0.0642*Math.sin(E3)
            + 0.0158*Math.sin(E4)
            + 0.0252*Math.sin(E5)
            - 0.0066*Math.sin(E6)
            - 0.0047*Math.sin(E7)
            - 0.0046*Math.sin(E8)
            + 0.0028*Math.sin(E9)
            + 0.0052*Math.sin(E10)
            + 0.0040*Math.sin(E11)
            + 0.0019*Math.sin(E12)
            - 0.0044*Math.sin(E13)
        );
        break;

    case Body.Mars:
        ra = (
            317.269202 - 0.10927547*T
            + 0.000068 * Math.sin(DEG2RAD*(198.991226 + 19139.4819985*T))
            + 0.000238 * Math.sin(DEG2RAD*(226.292679 + 38280.8511281*T))
            + 0.000052 * Math.sin(DEG2RAD*(249.663391 + 57420.7251593*T))
            + 0.000009 * Math.sin(DEG2RAD*(266.183510 + 76560.6367950*T))
            + 0.419057 * Math.sin(DEG2RAD*(79.398797 + 0.5042615*T))
        );

        dec = (
            54.432516 - 0.05827105*T
            + 0.000051*Math.cos(DEG2RAD*(122.433576 + 19139.9407476*T))
            + 0.000141*Math.cos(DEG2RAD*(43.058401 + 38280.8753272*T))
            + 0.000031*Math.cos(DEG2RAD*(57.663379 + 57420.7517205*T))
            + 0.000005*Math.cos(DEG2RAD*(79.476401 + 76560.6495004*T))
            + 1.591274*Math.cos(DEG2RAD*(166.325722 + 0.5042615*T))
        );

        w = (
            176.049863 + 350.891982443297*d
            + 0.000145*Math.sin(DEG2RAD*(129.071773 + 19140.0328244*T))
            + 0.000157*Math.sin(DEG2RAD*(36.352167 + 38281.0473591*T))
            + 0.000040*Math.sin(DEG2RAD*(56.668646 + 57420.9295360*T))
            + 0.000001*Math.sin(DEG2RAD*(67.364003 + 76560.2552215*T))
            + 0.000001*Math.sin(DEG2RAD*(104.792680 + 95700.4387578*T))
            + 0.584542*Math.sin(DEG2RAD*(95.391654 + 0.5042615*T))
        );
        break;

    case Body.Jupiter:
        const Ja = DEG2RAD*(99.360714 + 4850.4046*T);
        const Jb = DEG2RAD*(175.895369 + 1191.9605*T);
        const Jc = DEG2RAD*(300.323162 + 262.5475*T);
        const Jd = DEG2RAD*(114.012305 + 6070.2476*T);
        const Je = DEG2RAD*(49.511251 + 64.3000*T);

        ra = (
            268.056595 - 0.006499*T
            + 0.000117*Math.sin(Ja)
            + 0.000938*Math.sin(Jb)
            + 0.001432*Math.sin(Jc)
            + 0.000030*Math.sin(Jd)
            + 0.002150*Math.sin(Je)
        );

        dec = (
            64.495303 + 0.002413*T
            + 0.000050*Math.cos(Ja)
            + 0.000404*Math.cos(Jb)
            + 0.000617*Math.cos(Jc)
            - 0.000013*Math.cos(Jd)
            + 0.000926*Math.cos(Je)
        );

        w = 284.95 + 870.536*d;
        break;

    case Body.Saturn:
        ra = 40.589 - 0.036*T;
        dec = 83.537 - 0.004*T;
        w = 38.90 + 810.7939024*d;
        break;

    case Body.Uranus:
        ra = 257.311;
        dec = -15.175;
        w = 203.81 - 501.1600928*d;
        break;

    case Body.Neptune:
        const N = DEG2RAD*(357.85 + 52.316*T);
        ra = 299.36 + 0.70*Math.sin(N);
        dec = 43.46 - 0.51*Math.cos(N);
        w = 249.978 + 541.1397757*d - 0.48*Math.sin(N);
        break;

    case Body.Pluto:
        ra = 132.993;
        dec = -6.163;
        w = 302.695 + 56.3625225*d;
        break;

    default:
        throw `Invalid body: ${body}`;
    }

    // Calculate the north pole vector using the given angles.
    const radlat = dec * DEG2RAD;
    const radlon = ra * DEG2RAD;
    const rcoslat = Math.cos(radlat);
    const north = new Vector(
        rcoslat * Math.cos(radlon),
        rcoslat * Math.sin(radlon),
        Math.sin(radlat),
        time
    );

    return new AxisInfo(ra/15, dec, w, north);
}


/**
 * @brief Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.
 *
 * Given a more massive "major" body and a much less massive "minor" body,
 * calculates one of the five Lagrange points in relation to the minor body's
 * orbit around the major body. The parameter `point` is an integer that
 * selects the Lagrange point as follows:
 *
 * 1 = the Lagrange point between the major body and minor body.
 * 2 = the Lagrange point on the far side of the minor body.
 * 3 = the Lagrange point on the far side of the major body.
 * 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
 * 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
 *
 * The function returns the state vector for the selected Lagrange point
 * in equatorial J2000 coordinates (EQJ), with respect to the center of the
 * major body.
 *
 * To calculate Sun/Earth Lagrange points, pass in `Body.Sun` for `major_body`
 * and `Body.EMB` (Earth/Moon barycenter) for `minor_body`.
 * For Lagrange points of the Sun and any other planet, pass in that planet
 * (e.g. `Body.Jupiter`) for `minor_body`.
 * To calculate Earth/Moon Lagrange points, pass in `Body.Earth` and `Body.Moon`
 * for the major and minor bodies respectively.
 *
 * In some cases, it may be more efficient to call {@link LagrangePointFast},
 * especially when the state vectors have already been calculated, or are needed
 * for some other purpose.
 *
 * @param {number} point
 *      An integer 1..5 that selects which of the Lagrange points to calculate.
 *
 * @param {FlexibleDateTime} date
 *      The time at which the Lagrange point is to be calculated.
 *
 * @param {Body} major_body
 *      The more massive of the co-orbiting bodies: `Body.Sun` or `Body.Earth`.
 *
 * @param {Body} minor_body
 *      The less massive of the co-orbiting bodies. See main remarks.
 *
 * @returns {StateVector}
 *      The position and velocity of the selected Lagrange point with respect to the major body's center.
 */
export function LagrangePoint(
    point: number,
    date: FlexibleDateTime,
    major_body: Body,
    minor_body: Body
): StateVector {
    const time = MakeTime(date);
    const major_mass = MassProduct(major_body);
    const minor_mass = MassProduct(minor_body);

    let major_state: StateVector;
    let minor_state: StateVector;

    // Calculate the state vectors for the major and minor bodies.
    if (major_body === Body.Earth && minor_body === Body.Moon) {
        // Use geocentric calculations for more precision.
        // The Earth's geocentric state is trivial.
        major_state = new StateVector(0, 0, 0, 0, 0, 0, time);
        minor_state = GeoMoonState(time);
    } else {
        major_state = HelioState(major_body, time);
        minor_state = HelioState(minor_body, time);
    }

    return LagrangePointFast(
        point,
        major_state,
        major_mass,
        minor_state,
        minor_mass
    );
}


/**
 * @brief Calculates one of the 5 Lagrange points from body masses and state vectors.
 *
 * Given a more massive "major" body and a much less massive "minor" body,
 * calculates one of the five Lagrange points in relation to the minor body's
 * orbit around the major body. The parameter `point` is an integer that
 * selects the Lagrange point as follows:
 *
 * 1 = the Lagrange point between the major body and minor body.
 * 2 = the Lagrange point on the far side of the minor body.
 * 3 = the Lagrange point on the far side of the major body.
 * 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
 * 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
 *
 * The caller passes in the state vector and mass for both bodies.
 * The state vectors can be in any orientation and frame of reference.
 * The body masses are expressed as GM products, where G = the universal
 * gravitation constant and M = the body's mass. Thus the units for
 * `major_mass` and `minor_mass` must be au^3/day^2.
 * Use {@link MassProduct} to obtain GM values for various solar system bodies.
 *
 * The function returns the state vector for the selected Lagrange point
 * using the same orientation as the state vector parameters `major_state` and `minor_state`,
 * and the position and velocity components are with respect to the major body's center.
 *
 * Consider calling {@link LagrangePoint}, instead of this function, for simpler usage in most cases.
 *
 * @param {number} point
 *      A value 1..5 that selects which of the Lagrange points to calculate.
 *
 * @param {StateVector} major_state
 *      The state vector of the major (more massive) of the pair of bodies.
 *
 * @param {number} major_mass
 *      The mass product GM of the major body.
 *
 * @param {StateVector} minor_state
 *      The state vector of the minor (less massive) of the pair of bodies.
 *
 * @param {number} minor_mass
 *      The mass product GM of the minor body.
 *
 * @returns {StateVector}
 *      The position and velocity of the selected Lagrange point with respect to the major body's center.
 */
export function LagrangePointFast(
    point: number,
    major_state: StateVector,
    major_mass: number,
    minor_state: StateVector,
    minor_mass: number
): StateVector {
    const cos_60 = 0.5;
    const sin_60 = 0.8660254037844386;   /* sqrt(3) / 2 */

    if (point < 1 || point > 5)
        throw `Invalid lagrange point ${point}`;

    if (!Number.isFinite(major_mass) || major_mass <= 0.0)
        throw 'Major mass must be a positive number.';

    if (!Number.isFinite(minor_mass) || minor_mass <= 0.0)
        throw 'Minor mass must be a negative number.';

    // Find the relative position vector <dx, dy, dz>.
    let dx = minor_state.x - major_state.x;
    let dy = minor_state.y - major_state.y;
    let dz = minor_state.z - major_state.z;
    const R2 = (dx*dx + dy*dy + dz*dz);

    // R = Total distance between the bodies.
    const R = Math.sqrt(R2);

    // Find the relative velocity vector <vx, vy, vz>.
    const vx = minor_state.vx - major_state.vx;
    const vy = minor_state.vy - major_state.vy;
    const vz = minor_state.vz - major_state.vz;

    let p: StateVector;
    if (point === 4 || point === 5) {
        // For L4 and L5, we need to find points 60 degrees away from the
        // line connecting the two bodies and in the instantaneous orbital plane.
        // Define the instantaneous orbital plane as the unique plane that contains
        // both the relative position vector and the relative velocity vector.

        // Take the cross product of position and velocity to find a normal vector <nx, ny, nz>.
        const nx = dy*vz - dz*vy;
        const ny = dz*vx - dx*vz;
        const nz = dx*vy - dy*vx;

        // Take the cross product normal*position to get a tangential vector <ux, uy, uz>.
        let ux = ny*dz - nz*dy;
        let uy = nz*dx - nx*dz;
        let uz = nx*dy - ny*dx;

        // Convert the tangential direction vector to a unit vector.
        const U = Math.sqrt(ux*ux + uy*uy + uz*uz);
        ux /= U;
        uy /= U;
        uz /= U;

        // Convert the relative position vector into a unit vector.
        dx /= R;
        dy /= R;
        dz /= R;

        // Now we have two perpendicular unit vectors in the orbital plane: 'd' and 'u'.

        // Create new unit vectors rotated (+/-)60 degrees from the radius/tangent directions.
        const vert = (point == 4) ? +sin_60 : -sin_60;

        // Rotated radial vector
        const Dx = cos_60*dx + vert*ux;
        const Dy = cos_60*dy + vert*uy;
        const Dz = cos_60*dz + vert*uz;

        // Rotated tangent vector
        const Ux = cos_60*ux - vert*dx;
        const Uy = cos_60*uy - vert*dy;
        const Uz = cos_60*uz - vert*dz;

        // Calculate L4/L5 positions relative to the major body.
        const px = R * Dx;
        const py = R * Dy;
        const pz = R * Dz;

        // Use dot products to find radial and tangential components of the relative velocity.
        const vrad = vx*dx + vy*dy + vz*dz;
        const vtan = vx*ux + vy*uy + vz*uz;

        // Calculate L4/L5 velocities.
        const pvx = vrad*Dx + vtan*Ux;
        const pvy = vrad*Dy + vtan*Uy;
        const pvz = vrad*Dz + vtan*Uz;

        p = new StateVector(px, py, pz, pvx, pvy, pvz, major_state.t);
    } else {
        // Calculate the distances of each body from their mutual barycenter.
        // r1 = negative distance of major mass from barycenter (e.g. Sun to the left of barycenter)
        // r2 = positive distance of minor mass from barycenter (e.g. Earth to the right of barycenter)
        const r1 = -R * (minor_mass / (major_mass + minor_mass));
        const r2 = +R * (major_mass / (major_mass + minor_mass));

        // Calculate the square of the angular orbital speed in [rad^2 / day^2].
        const omega2 = (major_mass + minor_mass) / (R2*R);

        // Use Newton's Method to numerically solve for the location where
        // outward centrifugal acceleration in the rotating frame of reference
        // is equal to net inward gravitational acceleration.
        // First derive a good initial guess based on approximate analysis.
        let scale:number, numer1:number, numer2:number;
        if (point === 1 || point === 2) {
            scale = (major_mass / (major_mass + minor_mass)) * Math.cbrt(minor_mass / (3.0 * major_mass));
            numer1 = -major_mass;    // The major mass is to the left of L1 and L2.
            if (point == 1) {
                scale = 1.0 - scale;
                numer2 = +minor_mass;    // The minor mass is to the right of L1.
            } else {
                scale = 1.0 + scale;
                numer2 = -minor_mass;    // The minor mass is to the left of L2.
            }
        } else if (point === 3) {
            scale = ((7.0/12.0)*minor_mass - major_mass) / (minor_mass + major_mass);
            numer1 = +major_mass;    // major mass is to the right of L3.
            numer2 = +minor_mass;    // minor mass is to the right of L3.
        } else {
            throw `Invalid Langrage point ${point}. Must be an integer 1..5.`;
        }

        // Iterate Newton's Method until it converges.
        let x = R*scale - r1;
        let deltax: number;
        do {
            const dr1 = x - r1;
            const dr2 = x - r2;
            const accel = omega2*x + numer1/(dr1*dr1) + numer2/(dr2*dr2);
            const deriv = omega2 - 2*numer1/(dr1*dr1*dr1) - 2*numer2/(dr2*dr2*dr2);
            deltax = accel/deriv;
            x -= deltax;
        }
        while (Math.abs(deltax/R) > 1.0e-14);
        scale = (x - r1) / R;
        p = new StateVector(scale*dx, scale*dy, scale*dz, scale*vx, scale*vy, scale*vz, major_state.t);
    }
    return p;
}


/**
 * @brief A simulation of zero or more small bodies moving through the Solar System.
 *
 * This class calculates the movement of arbitrary small bodies,
 * such as asteroids or comets, that move through the Solar System.
 * It does so by calculating the gravitational forces on the bodies
 * from the Sun and planets. The user of this class supplies a
 * list of initial positions and velocities for the bodies.
 * Then the class can update the positions and velocities over small
 * time steps.
 *
 */
export class GravitySimulator {
    private originBody: Body;
    private prev: GravSimEndpoint;
    private curr: GravSimEndpoint;

    /**
     * @brief Creates a gravity simulation object.
     *
     * @param {Body} originBody
     *      Specifies the origin of the reference frame.
     *      All position vectors and velocity vectors will use `originBody`
     *      as the origin of the coordinate system.
     *      This origin applies to all the input vectors provided in the
     *      `bodyStates` parameter of this function, along with all
     *      output vectors returned by {@link GravitySimulator#Update}.
     *      Most callers will want to provide one of the following:
     *      `Body.Sun` for heliocentric coordinates,
     *      `Body.SSB` for solar system barycentric coordinates,
     *      or `Body.Earth` for geocentric coordinates. Note that the
     *      gravity simulator does not correct for light travel time;
     *      all state vectors are tied to a Newtonian "instantaneous" time.
     *
     * @param {FlexibleDateTime} date
     *      The initial time at which to start the simulation.
     *
     * @param {StateVector[]} bodyStates
     *      An array of zero or more initial state vectors (positions and velocities)
     *      of the small bodies to be simulated.
     *      The caller must know the positions and velocities of the small bodies at an initial moment in time.
     *      Their positions and velocities are expressed with respect to `originBody`, using equatorial
     *      J2000 orientation (EQJ).
     *      Positions are expressed in astronomical units (AU).
     *      Velocities are expressed in AU/day.
     *      All the times embedded within the state vectors must exactly match `date`,
     *      or this constructor will throw an exception.
     */
    constructor(
        originBody: Body,
        date: FlexibleDateTime,
        bodyStates: StateVector[])
    {
        const time = MakeTime(date);
        this.originBody = originBody;

        // Verify that the state vectors have matching times.
        for (let b of bodyStates) {
            if (b.t.tt !== time.tt) {
                throw 'Inconsistent times in bodyStates';
            }
        }

        // Create a stub list that we append to later.
        // We just need the stub to put into `this.curr`.
        const smallBodyList: body_grav_calc_t[] = [];

        // Calculate the states of the Sun and planets.
        const largeBodyDict = this.CalcSolarSystem(time);
        this.curr = new GravSimEndpoint(time, largeBodyDict, smallBodyList);

        // Convert origin-centric bodyStates vectors into barycentric body_grav_calc_t array.
        const o = this.InternalBodyState(originBody);
        for (let b of bodyStates) {
            const r = new TerseVector(b.x + o.r.x, b.y + o.r.y, b.z + o.r.z);
            const v = new TerseVector(b.vx + o.v.x, b.vy + o.v.y, b.vz + o.v.z);
            const a = TerseVector.zero();
            smallBodyList.push(new body_grav_calc_t(time.tt, r, v, a));
        }

        // Calculate the net acceleration experienced by the small bodies.
        this.CalcBodyAccelerations();

        // To prepare for a possible swap operation, duplicate the current state into the previous state.
        this.prev = this.Duplicate();
    }

    /**
     * @brief The body that was selected as the coordinate origin when this simulator was created.
     */
    public get OriginBody(): Body {
        return this.originBody;
    }

    /**
     * @brief The time represented by the current step of the gravity simulation.
     */
    public get Time(): AstroTime {
        return this.curr.time;
    }

    /**
     * Advances a gravity simulation by a small time step.
     *
     * Updates the simulation of the user-supplied small bodies
     * to the time indicated by the `date` parameter.
     * Returns an array of state vectors for the simulated bodies.
     * The array is in the same order as the original array that
     * was used to construct this simulator object.
     * The positions and velocities in the returned array are
     * referenced to the `originBody` that was used to construct
     * this simulator.
     *
     * @param {FlexibleDateTime} date
     *      A time that is a small increment away from the current simulation time.
     *      It is up to the developer to figure out an appropriate time increment.
     *      Depending on the trajectories, a smaller or larger increment
     *      may be needed for the desired accuracy. Some experimentation may be needed.
     *      Generally, bodies that stay in the outer Solar System and move slowly can
     *      use larger time steps. Bodies that pass into the inner Solar System and
     *      move faster will need a smaller time step to maintain accuracy.
     *      Some experimentation may be necessary to find a good value.
     *      The `date` value may be after or before the current simulation time
     *      to move forward or backward in time.
     */
    public Update(date: FlexibleDateTime): StateVector[] {
        const time = MakeTime(date);
        const dt = time.tt - this.curr.time.tt;
        if (dt === 0.0) {
            // Special case: the time has not changed, so skip the usual physics calculations.
            // This allows another way for the caller to query the current body states.
            // It is also necessary to avoid dividing by `dt` if `dt` is zero.
            // To prepare for a possible swap operation, duplicate the current state into the previous state.
            this.prev = this.Duplicate();
        } else {
            // Exchange the current state with the previous state. Then calculate the new current state.
            this.Swap();

            // Update the current time.
            this.curr.time = time;

            // Now that the time is set, it is safe to update the Solar System.
            this.curr.gravitators = this.CalcSolarSystem(time);

            // Estimate the position of each small body as if their existing
            // accelerations apply across the whole time interval.
            for (let i = 0; i < this.curr.bodies.length; ++i) {
                const p = this.prev.bodies[i];
                this.curr.bodies[i].r = UpdatePosition(dt, p.r, p.v, p.a);
            }

            // Calculate the acceleration experienced by the small bodies at
            // their respective approximate next locations.
            this.CalcBodyAccelerations();

            for (let i = 0; i < this.curr.bodies.length; ++i) {
                // Calculate the average of the acceleration vectors
                // experienced by the previous body positions and
                // their estimated next positions.
                // These become estimates of the mean effective accelerations
                // over the whole interval.
                const p = this.prev.bodies[i];
                const c = this.curr.bodies[i];
                const acc = p.a.mean(c.a);

                // Refine the estimates of position and velocity at the next time step,
                // using the mean acceleration as a better approximation of the
                // continuously changing acceleration acting on each body.
                c.tt = time.tt;
                c.r = UpdatePosition(dt, p.r, p.v, acc);
                c.v = UpdateVelocity(dt, p.v, acc);
            }

            // Re-calculate accelerations experienced by each body.
            // These will be needed for the next simulation step (if any).
            // Also, they will be potentially useful if some day we add
            // a function to query the acceleration vectors for the bodies.
            this.CalcBodyAccelerations();
        }

        // Translate our internal calculations of body positions and velocities
        // into state vectors that the caller can understand.
        // We have to convert the internal type body_grav_calc_t to the public type StateVector.
        // Also convert from barycentric coordinates to coordinates based on the selected origin body.
        const bodyStates: StateVector[] = [];
        const ostate = this.InternalBodyState(this.originBody);
        for (let bcalc of this.curr.bodies) {
            bodyStates.push(new StateVector(
                bcalc.r.x - ostate.r.x,
                bcalc.r.y - ostate.r.y,
                bcalc.r.z - ostate.r.z,
                bcalc.v.x - ostate.v.x,
                bcalc.v.y - ostate.v.y,
                bcalc.v.z - ostate.v.z,
                time
            ));
        }
        return bodyStates;
    }

    /**
     * Exchange the current time step with the previous time step.
     *
     * Sometimes it is helpful to "explore" various times near a given
     * simulation time step, while repeatedly returning to the original
     * time step. For example, when backdating a position for light travel
     * time, the caller may wish to repeatedly try different amounts of
     * backdating. When the backdating solver has converged, the caller
     * wants to leave the simulation in its original state.
     *
     * This function allows a single "undo" of a simulation, and does so
     * very efficiently.
     *
     * Usually this function will be called immediately after a matching
     * call to {@link GravitySimulator#Update}. It has the effect of rolling
     * back the most recent update. If called twice in a row, it reverts
     * the swap and thus has no net effect.
     *
     * The constructor initializes the current state and previous
     * state to be identical. Both states represent the `time` parameter that was
     * passed into the constructor. Therefore, `Swap` will
     * have no effect from the caller's point of view when passed a simulator
     * that has not yet been updated by a call to {@link GravitySimulator#Update}.
     */
    public Swap(): void {
        const swap = this.curr;
        this.curr = this.prev;
        this.prev = swap;
    }

    /**
     * Get the position and velocity of a Solar System body included in the simulation.
     *
     * In order to simulate the movement of small bodies through the Solar System,
     * the simulator needs to calculate the state vectors for the Sun and planets.
     *
     * If an application wants to know the positions of one or more of the planets
     * in addition to the small bodies, this function provides a way to obtain
     * their state vectors. This is provided for the sake of efficiency, to avoid
     * redundant calculations.
     *
     * The state vector is returned relative to the position and velocity
     * of the `originBody` parameter that was passed to this object's constructor.
     *
     * @param {Body} body
     *      The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune.
     */
    public SolarSystemBodyState(body: Body): StateVector {
        const bstate = this.InternalBodyState(body);
        const ostate = this.InternalBodyState(this.originBody);
        return ExportState(bstate.sub(ostate), this.curr.time);
    }

    private InternalBodyState(body: Body): body_state_t {
        if (body === Body.SSB)
            return new body_state_t(this.curr.time.tt, TerseVector.zero(), TerseVector.zero());

        const bstate = this.curr.gravitators[body];
        if (bstate)
            return bstate;

        throw `Invalid body: ${body}`;
    }

    private CalcSolarSystem(time: AstroTime): BodyStateTable {
        const dict: BodyStateTable = {};

        // Start with the SSB at zero position and velocity.
        const ssb = new body_state_t(time.tt, TerseVector.zero(), TerseVector.zero());

        // Calculate the heliocentric position of each planet, and adjust the SSB
        // based each planet's pull on the Sun.
        dict[Body.Mercury] = AdjustBarycenterPosVel(ssb, time.tt, Body.Mercury, MERCURY_GM);
        dict[Body.Venus  ] = AdjustBarycenterPosVel(ssb, time.tt, Body.Venus,   VENUS_GM);
        dict[Body.Earth  ] = AdjustBarycenterPosVel(ssb, time.tt, Body.Earth,   EARTH_GM + MOON_GM);
        dict[Body.Mars   ] = AdjustBarycenterPosVel(ssb, time.tt, Body.Mars,    MARS_GM);
        dict[Body.Jupiter] = AdjustBarycenterPosVel(ssb, time.tt, Body.Jupiter, JUPITER_GM);
        dict[Body.Saturn ] = AdjustBarycenterPosVel(ssb, time.tt, Body.Saturn,  SATURN_GM);
        dict[Body.Uranus ] = AdjustBarycenterPosVel(ssb, time.tt, Body.Uranus,  URANUS_GM);
        dict[Body.Neptune] = AdjustBarycenterPosVel(ssb, time.tt, Body.Neptune, NEPTUNE_GM);

        // Convert planet states from heliocentric to barycentric.
        for (let body in dict) {
            dict[body].r.decr(ssb.r);
            dict[body].v.decr(ssb.v);
        }

        // Convert heliocentric SSB to barycentric Sun.
        dict[Body.Sun] = new body_state_t(time.tt, ssb.r.neg(), ssb.v.neg());

        return dict;
    }

    private CalcBodyAccelerations(): void {
        // Calculate the gravitational acceleration experienced by the simulated small bodies.
        for (let b of this.curr.bodies) {
            b.a = TerseVector.zero();
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Sun    ].r,  SUN_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Mercury].r,  MERCURY_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Venus  ].r,  VENUS_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Earth  ].r,  EARTH_GM + MOON_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Mars   ].r,  MARS_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Jupiter].r,  JUPITER_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Saturn ].r,  SATURN_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Uranus ].r,  URANUS_GM);
            this.AddAcceleration(b.a, b.r, this.curr.gravitators[Body.Neptune].r,  NEPTUNE_GM);
        }
    }

    private AddAcceleration(acc: TerseVector, smallPos: TerseVector, majorPos: TerseVector, gm: number): void {
        const dx = majorPos.x - smallPos.x;
        const dy = majorPos.y - smallPos.y;
        const dz = majorPos.z - smallPos.z;
        const r2 = dx*dx + dy*dy + dz*dz;
        const pull = gm / (r2 * Math.sqrt(r2));
        acc.x += dx * pull;
        acc.y += dy * pull;
        acc.z += dz * pull;
    }

    private Duplicate(): GravSimEndpoint {
        // Copy the current state into the previous state, so that both become the same moment in time.
        const gravitators: BodyStateTable = {};
        for (let body in this.curr.gravitators) {
            gravitators[body] = this.curr.gravitators[body].clone();
        }

        const bodies: body_grav_calc_t[] = [];
        for (let b of this.curr.bodies) {
            bodies.push(b.clone());
        }

        return new GravSimEndpoint(this.curr.time, gravitators, bodies);
    }
}

interface BodyStateTable {
    [body: string]: body_state_t;
}

class GravSimEndpoint {
    constructor(
        public time: AstroTime,
        public gravitators: BodyStateTable,
        public bodies: body_grav_calc_t[]
    ) {}
}
