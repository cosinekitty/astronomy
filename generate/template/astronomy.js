/* astronomy.js - Don Cross - cosinekitty@gmail.com */
'use strict';

(function(Astronomy){
'use strict';
const j2000 = new Date('2000-01-01T12:00:00Z');
const T0 = 2451545.0;
const MJD_BASIS = 2400000.5;            // mjd + MJD_BASIS = jd
const Y2000_IN_MJD = T0 - MJD_BASIS;    // the 2000.0 epoch expressed in MJD
const PI2 = 2 * Math.PI;
const ARC = 3600 * (180 / Math.PI);     // arcseconds per radian
const ERAD = 6378136.6;                 // mean earth radius
const AU = 1.4959787069098932e+11;      // astronomical unit in meters
const C_AUDAY = 173.1446326846693;      // speed of light in AU/day
const ASEC2RAD = 4.848136811095359935899141e-6;
const DEG2RAD = 0.017453292519943296;
const RAD2DEG = 57.295779513082321;
const ASEC360 = 1296000;
const ANGVEL = 7.2921150e-5;
const AU_KM = 1.4959787069098932e+8;

function Frac(x) {
    return x - Math.floor(x);
}

Astronomy.Bodies = [
    'Sun',
    'Moon',
    'Mercury',
    'Venus',
    'Earth',
    'Mars',
    'Jupiter',
    'Saturn',
    'Uranus',
    'Neptune',
    'Pluto'
];

const vsop = {
    Mercury: $ASTRO_JS_VSOP(Mercury),
    Venus:   $ASTRO_JS_VSOP(Venus),
    Earth:   $ASTRO_JS_VSOP(Earth),
    Mars:    $ASTRO_JS_VSOP(Mars),
    Jupiter: $ASTRO_JS_VSOP(Jupiter),
    Saturn:  $ASTRO_JS_VSOP(Saturn),
    Uranus:  $ASTRO_JS_VSOP(Uranus),
    Neptune: $ASTRO_JS_VSOP(Neptune)
};

const cheb = {
    Pluto:  $ASTRO_JS_CHEBYSHEV(8)
};

const DT = $ASTRO_JS_DELTA_T();

function DeltaT(mjd) {
    // Calculate the difference TT-UT for the given date/time, expressed
    // as a Modified Julian Date.
    // DT[i] = { mjd: 58484.0, dt: 69.34 }
    // Check end ranges. If outside the known bounds, clamp to the closest known value.

    if (mjd <= DT[0].mjd) {
        return DT[0].dt;
    }

    if (mjd >= DT[DT.length-1].mjd) {
        return DT[DT.length-1].dt;
    }

    // Do a binary search to find the pair of indexes this mjd lies between.
    
    let lo = 0;
    let hi = DT.length-2;   // make sure there is always an array element after the one we are looking at
    while (lo <= hi) {
        let c = (lo + hi) >> 1;
        if (mjd < DT[c].mjd) {
            hi = c-1;
        } else if (mjd > DT[c+1].mjd) {
            lo = c+1;
        } else {
            let frac = (mjd - DT[c].mjd) / (DT[c+1].mjd - DT[c].mjd);
            return DT[c].dt + frac*(DT[c+1].dt - DT[c].dt);
        }
    }

    // This should never happen if the binary search algorithm is correct.
    throw `Could not find Delta-T value for MJD=${mjd}`;
}

function Time(date) {
    const MillisPerDay = 1000 * 3600 * 24;

    if (date instanceof Date) {
        this.date = date;
        this.ut = (date - j2000) / MillisPerDay;        // ut = days since J2000.0
        this.jd_utc = this.ut + T0;

        let dt = DeltaT(this.ut + Y2000_IN_MJD) / 86400;
        this.tt = this.ut + dt;
        this.jd_tt = this.jd_utc + dt;
        return;
    }

    if (typeof date === 'number') {
        this.date = new Date(j2000 - (-date)*MillisPerDay);
        this.ut = date;
        this.jd_utc = this.ut + T0;

        let dt = DeltaT(this.ut + Y2000_IN_MJD) / 86400.0;
        this.tt = this.ut + dt;
        this.jd_tt = this.jd_utc + dt;
        return;
    }

    throw 'AstroTime() argument must be a Date object, a Time object, or a numeric UTC Julian date.';
}

Time.prototype.SubtractDays = function(days) {
    // FIXFIXFIX: rework to subtract TT instead of UT. This requires an inverse function for DeltaT.
    return new Time(this.ut - days);
}

function AstroTime(date) {
    if (date instanceof Time) {
        return date;
    }
    return new Time(date);
}

var nals_t = [
    [ 0,    0,    0,    0,    1],
    [ 0,    0,    2,   -2,    2],
    [ 0,    0,    2,    0,    2],
    [ 0,    0,    0,    0,    2],
    [ 0,    1,    0,    0,    0],
    [ 0,    1,    2,   -2,    2],
    [ 1,    0,    0,    0,    0],
    [ 0,    0,    2,    0,    1],
    [ 1,    0,    2,    0,    2],
    [ 0,   -1,    2,   -2,    2],
    [ 0,    0,    2,   -2,    1],
    [-1,    0,    2,    0,    2],
    [-1,    0,    0,    2,    0],
    [ 1,    0,    0,    0,    1],
    [-1,    0,    0,    0,    1],
    [-1,    0,    2,    2,    2],
    [ 1,    0,    2,    0,    1],
    [-2,    0,    2,    0,    1],
    [ 0,    0,    0,    2,    0],
    [ 0,    0,    2,    2,    2],
    [ 0,   -2,    2,   -2,    2],
    [-2,    0,    0,    2,    0],
    [ 2,    0,    2,    0,    2],
    [ 1,    0,    2,   -2,    2],
    [-1,    0,    2,    0,    1],
    [ 2,    0,    0,    0,    0],
    [ 0,    0,    2,    0,    0],
    [ 0,    1,    0,    0,    1],
    [-1,    0,    0,    2,    1],
    [ 0,    2,    2,   -2,    2],
    [ 0,    0,   -2,    2,    0],
    [ 1,    0,    0,   -2,    1],
    [ 0,   -1,    0,    0,    1],
    [-1,    0,    2,    2,    1],
    [ 0,    2,    0,    0,    0],
    [ 1,    0,    2,    2,    2],
    [-2,    0,    2,    0,    0],
    [ 0,    1,    2,    0,    2],
    [ 0,    0,    2,    2,    1],
    [ 0,   -1,    2,    0,    2],
    [ 0,    0,    0,    2,    1],
    [ 1,    0,    2,   -2,    1],
    [ 2,    0,    2,   -2,    2],
    [-2,    0,    0,    2,    1],
    [ 2,    0,    2,    0,    1],
    [ 0,   -1,    2,   -2,    1],
    [ 0,    0,    0,   -2,    1],
    [-1,   -1,    0,    2,    0],
    [ 2,    0,    0,   -2,    1],
    [ 1,    0,    0,    2,    0],
    [ 0,    1,    2,   -2,    1],
    [ 1,   -1,    0,    0,    0],
    [-2,    0,    2,    0,    2],
    [ 3,    0,    2,    0,    2],
    [ 0,   -1,    0,    2,    0],
    [ 1,   -1,    2,    0,    2],
    [ 0,    0,    0,    1,    0],
    [-1,   -1,    2,    2,    2],
    [-1,    0,    2,    0,    0],
    [ 0,   -1,    2,    2,    2],
    [-2,    0,    0,    0,    1],
    [ 1,    1,    2,    0,    2],
    [ 2,    0,    0,    0,    1],
    [-1,    1,    0,    1,    0],
    [ 1,    1,    0,    0,    0],
    [ 1,    0,    2,    0,    0],
    [-1,    0,    2,   -2,    1],
    [ 1,    0,    0,    0,    2],
    [-1,    0,    0,    1,    0],
    [ 0,    0,    2,    1,    2],
    [-1,    0,    2,    4,    2],
    [-1,    1,    0,    1,    1],
    [ 0,   -2,    2,   -2,    1],
    [ 1,    0,    2,    2,    1],
    [-2,    0,    2,    2,    2],
    [-1,    0,    0,    0,    2],
    [ 1,    1,    2,   -2,    2]
];

var cls_t = [
    [-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0],
    [ -13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0],
    [  -2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0],
    [   2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0],
    [   1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0],
    [   -516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0],
    [    711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0],
    [   -387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0],
    [   -301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0],
    [    215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0],
    [    128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0],
    [    123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0],
    [    156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0],
    [     63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0],
    [    -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0],
    [    -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0],
    [    -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0],
    [     45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0],
    [     63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0],
    [    -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0],
    [     32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0],
    [    -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0],
    [    -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0],
    [     28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0],
    [     20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0],
    [     29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0],
    [     25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0],
    [    -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0],
    [     15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0],
    [    -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0],
    [     21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0],
    [    -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0],
    [    -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0],
    [    -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0],
    [     16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0],
    [     -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0],
    [    -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0],
    [      7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0],
    [     -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0],
    [     -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0],
    [     -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0],
    [      5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0],
    [      6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0],
    [     -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0],
    [     -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0],
    [     -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0],
    [     -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0],
    [      7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0],
    [      4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0],
    [      6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0],
    [      3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0],
    [      4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0],
    [     -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0],
    [     -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0],
    [      4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0],
    [     -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0],
    [     -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0],
    [     -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0],
    [     -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0],
    [     -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0],
    [     -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0],
    [      2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0],
    [      2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0],
    [      3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0],
    [     -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0],
    [      3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0],
    [     -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0],
    [     -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0],
    [      4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0],
    [      1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0],
    [     -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0],
    [      1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0],
    [     -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0],
    [     -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0],
    [      1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0],
    [      1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0],
    [      1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0]
];

function iau2000b(jd_high, jd_low) {
    var i, t, el, elp, f, d, om, arg, dp, de, sarg, carg;

    function mod(x) {
        return (x % ASEC360) * ASEC2RAD;
    }

    t = ((jd_high - T0) + jd_low) / 36525.0;
    el  = mod(485868.249036 + t * 1717915923.2178);
    elp = mod(1287104.79305 + t * 129596581.0481);
    f   = mod(335779.526232 + t * 1739527262.8478);
    d   = mod(1072260.70369 + t * 1602961601.2090);
    om  = mod(450160.398036 - t * 6962890.5431);
    dp = 0;
    de = 0;
    for (i=76; i >= 0; --i) {
        arg = (nals_t[i][0]*el + nals_t[i][1]*elp + nals_t[i][2]*f + nals_t[i][3]*d + nals_t[i][4]*om) % PI2;
        sarg = Math.sin(arg);
        carg = Math.cos(arg);
        dp += (cls_t[i][0] + cls_t[i][1] * t) * sarg + cls_t[i][2] * carg;
        de += (cls_t[i][3] + cls_t[i][4] * t) * carg + cls_t[i][5] * sarg;
    }
    return {
        dpsi: (-0.000135 * ASEC2RAD) + (dp * 1.0e-7 * ASEC2RAD),
        deps: (+0.000388 * ASEC2RAD) + (de * 1.0e-7 * ASEC2RAD)
    };
 }

function nutation_angles(t) {
    var t1 = t * 36525;
    var nut = iau2000b(T0, t1);
    return { dpsi: nut.dpsi/ASEC2RAD, deps: nut.deps/ASEC2RAD };
}

function mean_obliq(jd) {
    var t = (jd - T0) / 36525;
    return (
        (((( -  0.0000000434   * t
             -  0.000000576  ) * t
             +  0.00200340   ) * t
             -  0.0001831    ) * t
             - 46.836769     ) * t + 84381.406
    );
}

var cache_e_tilt;

function e_tilt(jd) {
    if (!cache_e_tilt || Math.abs(cache_e_tilt.jd - jd) > 1.0e-6) {
        const t = (jd - T0) / 36525;
        const nut = nutation_angles(t);
        const mean_obl_seconds = mean_obliq(jd);
        const mean_ob = mean_obl_seconds / 3600;
        const true_ob = (mean_obl_seconds + nut.deps) / 3600;
        cache_e_tilt = {
            jd: jd,
            dpsi: nut.dpsi,
            deps: nut.deps,
            ee: nut.dpsi * Math.cos(mean_ob * DEG2RAD) / 15,
            mobl: mean_ob,
            tobl: true_ob
        };
    }
    return cache_e_tilt;
}

function ecl2equ_vec(jd, pos) {
    var obl = e_tilt(jd).mobl * DEG2RAD;
    var cos_obl = Math.cos(obl);
    var sin_obl = Math.sin(obl);
    return [
        pos[0],
        pos[1]*cos_obl - pos[2]*sin_obl,
        pos[1]*sin_obl + pos[2]*cos_obl
    ];
}

function CalcMoon(T) {
    function DeclareArray1(xmin, xmax) {
        var array = [];
        var i;
        for (i=0; i <= xmax-xmin; ++i) {
            array.push(0);
        }
        return {min:xmin, array:array};
    }

    function DeclareArray2(xmin, xmax, ymin, ymax) {
        var array = [];
        var i;
        for (i=0; i <= xmax-xmin; ++i) {
            array.push(DeclareArray1(ymin, ymax));
        }
        return {min:xmin, array:array};
    }

    function ArrayGet2(a, x, y) {
        var m = a.array[x - a.min];
        return m.array[y - m.min];
    }

    function ArraySet2(a, x, y, v) {
        var m = a.array[x - a.min];
        m.array[y - m.min] = v;
    }

    var S, MAX, ARG, FAC, I, J, T2, DGAM, DLAM, N, GAM1C, SINPI, L0, L, LS, F, D, DL0, DL, DLS, DF, DD, DS;
    var coArray = DeclareArray2(-6, 6, 1, 4);
    var siArray = DeclareArray2(-6, 6, 1, 4);

    function CO(x, y) {
        return ArrayGet2(coArray, x, y);
    }

    function SI(x, y) {
        return ArrayGet2(siArray, x, y);
    }

    function SetCO(x, y, v) {
        return ArraySet2(coArray, x, y, v);
    }

    function SetSI(x, y, v) {
        return ArraySet2(siArray, x, y, v);
    }

    function AddThe(c1, s1, c2, s2, func) {
        return func(c1*c2 - s1*s2, s1*c2 + c1*s2);
    }

    function Sine(phi) {
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
    for (I=1; I<=4; ++I)
    {
        switch (I)
        {
            case 1: ARG=L;  MAX=4; FAC=1.000002208;               break;
            case 2: ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
            case 3: ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
            case 4: ARG=D;  MAX=6; FAC=1.0;                       break;
        }
        SetCO(0, 1, 1);
        SetCO(1, I, Math.cos(ARG) * FAC);
        SetSI(0, I, 0);
        SetSI(1, I, Math.sin(ARG) * FAC);
        for (J=2; J<=MAX; ++J) {
            AddThe(CO(J-1,I), SI(J-1,I), CO(1,I), SI(1,I), (c, s) => (SetCO(J,I,c), SetSI(J,I,s)));
        }
        for (J=1; J<=MAX; ++J) {
            SetCO(-J, I, CO(J, I));
            SetSI(-J, I, -SI(J, I));
        }
    }

    function Term(p, q, r, s) {
        var result = { x:1, y:0 };
        var I = [ null, p, q, r, s ];
        for (var k=1; k <= 4; ++k)
            if (I[k] !== 0) 
                AddThe(result.x, result.y, CO(I[k], k), SI(I[k], k), (c, s) => (result.x=c, result.y=s));
        return result;
    }

    function AddSol(coeffl, coeffs, coeffg, coeffp, p, q, r, s) {
        var result = Term(p, q, r, s);
        DLAM += coeffl * result.y;
        DS += coeffs * result.y;
        GAM1C += coeffg * result.x;
        SINPI += coeffp * result.x;
    }

    AddSol(    13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4);
    AddSol(     0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3);
    AddSol(  2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2);
    AddSol(  -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1);
    AddSol(     1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4);
    AddSol(   191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2);
    AddSol(    -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1);
    AddSol( 22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0);
    AddSol(    18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1);
    AddSol( -4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2);
    AddSol(    +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3);
    AddSol(   -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4);
    AddSol(    -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6);
    AddSol(    -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4);
    AddSol(   -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2);
    AddSol(    18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1);
    AddSol(  -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0);
    AddSol(     0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1);
    AddSol(  -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2);
    AddSol(    -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4);
    AddSol(     0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4);
    AddSol(    14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2);
    AddSol(    -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1);
    AddSol(   769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0);
    AddSol(    +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1);
    AddSol(  -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2);
    AddSol(    +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3);
    AddSol(   -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4);
    AddSol(    -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6);
    AddSol(    -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2);
    AddSol(    +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1);
    AddSol(  -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0);
    AddSol(  -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2);
    AddSol(     0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3);
    AddSol(    -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4);

    AddSol(     0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4);
    AddSol(    14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2);
    AddSol(   147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0);
    AddSol(    -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1);
    AddSol(    28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2);
    AddSol(    -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3);
    AddSol(     0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4);
    AddSol(    -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2);
    AddSol(    -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0);
    AddSol(    -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2);
    AddSol(    -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2);
    AddSol(     0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1);
    AddSol(  -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0);
    AddSol(     0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1);
    AddSol(   -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2);
    AddSol(     0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3);
    AddSol(    +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4);
    AddSol(     1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2);
    AddSol(    36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0);
    AddSol(   -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2);
    AddSol(    -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4);
    AddSol(    -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6);
    AddSol(    -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2);
    AddSol(    -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0);
    AddSol(    -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2);
    AddSol(    -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4);
    AddSol(     1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2);
    AddSol(     9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0);
    AddSol(    -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1);
    AddSol(    -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2);
    AddSol(     0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4);
    AddSol(    -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0);
    AddSol(    -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2);
    AddSol(    -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4);
    AddSol(    +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2);
    AddSol(    +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0);
    AddSol(    +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2);
    AddSol(    -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2);
    AddSol(    -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2);
    AddSol(   -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0);
    AddSol(    -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2);
    AddSol(    -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4);
    AddSol(    -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2);
    AddSol(    39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0);
    AddSol(     9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2);
    AddSol(     0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4);

    AddSol(     0.415,    0.10, 0.0  ,  0.0013,0, 1, 2, 0);
    AddSol(    -2.152,   -2.26, 0.0  , -0.0066,0, 1, 2,-2);
    AddSol(    -1.440,   -1.30, 0.0  , +0.0014,0, 1,-2, 2);
    AddSol(     0.384,   -0.04, 0.0  ,  0.0   ,0, 1,-2,-2);
    AddSol(    +1.938,   +3.60,-0.145, +0.0401,4, 0, 0, 0);
    AddSol(    -0.952,   -1.58,+0.052, -0.0130,4, 0, 0,-2);
    AddSol(    -0.551,   -0.94,+0.032, -0.0097,3, 1, 0, 0);
    AddSol(    -0.482,   -0.57,+0.005, -0.0045,3, 1, 0,-2);
    AddSol(     0.681,    0.96,-0.026,  0.0115,3,-1, 0, 0);
    AddSol(    -0.297,   -0.27, 0.002, -0.0009,2, 2, 0,-2);
    AddSol(     0.254,   +0.21,-0.003,  0.0   ,2,-2, 0,-2);
    AddSol(    -0.250,   -0.22, 0.004,  0.0014,1, 3, 0,-2);
    AddSol(    -3.996,    0.0 , 0.0  , +0.0004,2, 0, 2, 0);
    AddSol(     0.557,   -0.75, 0.0  , -0.0090,2, 0, 2,-2);
    AddSol(    -0.459,   -0.38, 0.0  , -0.0053,2, 0,-2, 2);
    AddSol(    -1.298,    0.74, 0.0  , +0.0004,2, 0,-2, 0);
    AddSol(     0.538,    1.14, 0.0  , -0.0141,2, 0,-2,-2);
    AddSol(     0.263,    0.02, 0.0  ,  0.0   ,1, 1, 2, 0);
    AddSol(     0.426,   +0.07, 0.0  , -0.0006,1, 1,-2,-2);
    AddSol(    -0.304,   +0.03, 0.0  , +0.0003,1,-1, 2, 0);
    AddSol(    -0.372,   -0.19, 0.0  , -0.0027,1,-1,-2, 2);
    AddSol(    +0.418,    0.0 , 0.0  ,  0.0   ,0, 0, 4, 0);
    AddSol(    -0.330,   -0.04, 0.0  ,  0.0   ,3, 0, 2, 0);

    function ADDN(coeffn, p, q, r, s) {
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

    var lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*Math.sin(S) - 6.24*Math.sin(3*S) + N;

    return {
        geo_eclip_lon: PI2 * Frac((L0+DLAM/ARC) / PI2),
        geo_eclip_lat: (Math.PI / (180 * 3600)) * lat_seconds,
        distance_au: (ARC * (ERAD / AU)) / (0.999953253 * SINPI)
    };
}

function precession(jd_tdb1, pos1, jd_tdb2) {
    var xx, yx, zx, xy, yy, zy, xz, yz, zz;
    var eps0 = 84381.406;
    var t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;

    if ((jd_tdb1 != T0) && (jd_tdb2 != T0))
        throw 'One of the two epochs must be J2000.0.';

    t = (jd_tdb2 - jd_tdb1) / 36525;
    if (jd_tdb2 === T0)
        t = -t;

    psia   = (((((-    0.0000000951  * t
                 +    0.000132851 ) * t
                 -    0.00114045  ) * t
                 -    1.0790069   ) * t
                 + 5038.481507    ) * t);

    omegaa = (((((+    0.0000003337  * t
                 -    0.000000467 ) * t
                 -    0.00772503  ) * t
                 +    0.0512623   ) * t
                 -    0.025754    ) * t + eps0);

    chia   = (((((-    0.0000000560  * t
                 +    0.000170663 ) * t
                 -    0.00121197  ) * t
                 -    2.3814292   ) * t
                 +   10.556403    ) * t);

    eps0 = eps0 * ASEC2RAD;
    psia = psia * ASEC2RAD;
    omegaa = omegaa * ASEC2RAD;
    chia = chia * ASEC2RAD;

    sa = Math.sin(eps0);
    ca = Math.cos(eps0);
    sb = Math.sin(-psia);
    cb = Math.cos(-psia);
    sc = Math.sin(-omegaa);
    cc = Math.cos(-omegaa);
    sd = Math.sin(chia);
    cd = Math.cos(chia);

    xx =  cd * cb - sb * sd * cc;
    yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;
    zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
    xy = -sd * cb - sb * cd * cc;
    yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;
    zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
    xz =  sb * sc;
    yz = -sc * cb * ca - sa * cc;
    zz = -sc * cb * sa + cc * ca;

    if (jd_tdb2 == T0)
    { 
        // Perform rotation from epoch to J2000.0.
        return [
            xx * pos1[0] + xy * pos1[1] + xz * pos1[2],
            yx * pos1[0] + yy * pos1[1] + yz * pos1[2],
            zx * pos1[0] + zy * pos1[1] + zz * pos1[2]
        ];
    }

    // Perform rotation from J2000.0 to epoch.
    return [
        xx * pos1[0] + yx * pos1[1] + zx * pos1[2],
        xy * pos1[0] + yy * pos1[1] + zy * pos1[2],
        xz * pos1[0] + yz * pos1[1] + zz * pos1[2]
    ]; 
}

function era(tt_ut1) {    // Earth Rotation Angle
    const thet1 = 0.7790572732640 + 0.00273781191135448 * (tt_ut1 - T0);
    const thet3 = tt_ut1 % 1;
    let theta = 360 * ((thet1 + thet3) % 1);
    if (theta < 0) {
        theta += 360;
    }
    return theta;
}

function sidereal_time(time, gst_type) {
    const t = (time.jd_tt - T0) / 36525;
    let eqeq;

    switch (gst_type) {
    case 'gmst':
        eqeq = 0.0;
        break;

    case 'gast':
        eqeq = 15 * e_tilt(time.jd_tt).ee;
        break;

    default:
        throw 'gst_type must be either "gmst" or "gast"';
    }
    
    const theta = era(time.jd_utc);

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
    return gst;
}

function terra(observer, st) {
    const erad_km = ERAD / 1000;
    const df = 1 - 0.003352819697896;    // flattening of the Earth
    const df2 = df * df;
    const phi = observer.latitude * DEG2RAD;
    const sinphi = Math.sin(phi);
    const cosphi = Math.cos(phi);
    const c = 1 / Math.sqrt(cosphi*cosphi + df2*sinphi*sinphi);
    const s = df2 * c;
    const ht_km = observer.height / 1000;
    const ach = erad_km*c + ht_km;
    const ash = erad_km*s + ht_km;
    const stlocl = (15*st + observer.longitude) * DEG2RAD;
    const sinst = Math.sin(stlocl);
    const cosst = Math.cos(stlocl);
    return {
        pos: [ach*cosphi*cosst/AU_KM, ach*cosphi*sinst/AU_KM, ash*sinphi/AU_KM],
        vel: [-ANGVEL*ach*cosphi*sinst*86400, ANGVEL*ach*cosphi*cosst*86400, 0]
    };
}

function nutation(jd_tdb, direction, pos) {
    const tilt = e_tilt(jd_tdb);
    const oblm = tilt.mobl * DEG2RAD;
    const oblt = tilt.tobl * DEG2RAD;
    const psi = tilt.dpsi * ASEC2RAD;
    const cobm = Math.cos(oblm);
    const sobm = Math.sin(oblm);
    const cobt = Math.cos(oblt);
    const sobt = Math.sin(oblt);
    const cpsi = Math.cos(psi);
    const spsi = Math.sin(psi);

    const xx = cpsi;
    const yx = -spsi * cobm;
    const zx = -spsi * sobm;
    const xy = spsi * cobt;
    const yy = cpsi * cobm * cobt + sobm * sobt;
    const zy = cpsi * sobm * cobt - cobm * sobt;
    const xz = spsi * sobt;
    const yz = cpsi * cobm * sobt - sobm * cobt;
    const zz = cpsi * sobm * sobt + cobm * cobt; 

    if (direction === 0) {
        // forward rotation
        return [
            xx * pos[0] + yx * pos[1] + zx * pos[2],
            xy * pos[0] + yy * pos[1] + zy * pos[2],
            xz * pos[0] + yz * pos[1] + zz * pos[2]
        ];
    }

    // inverse rotation
    return [
        xx * pos[0] + xy * pos[1] + xz * pos[2],
        yx * pos[0] + yy * pos[1] + yz * pos[2],
        zx * pos[0] + zy * pos[1] + zz * pos[2]
    ];
}

function geo_pos(time, observer) {
    const gmst = sidereal_time(time, 'gmst');
    const tilt = e_tilt(time.jd_tt);
    const gast = gmst + tilt.ee/3600;
    const pos1 = terra(observer, gast).pos;
    const pos2 = nutation(time.jd_tt, -1, pos1);
    const pos3 = precession(time.jd_tt, pos2, T0);
    return pos3;
}

function vector2radec(pos)
{
    const xyproj = pos[0]*pos[0] + pos[1]*pos[1];
    const dist = Math.sqrt(xyproj + pos[2]*pos[2]);
    if (xyproj === 0)
    {
        if (pos[2] === 0)
            throw 'Indeterminate sky coordinates';

        if (pos[2] < 0)
            return { ra:0, dec:-90, dist:dist };

        return { ra:0, dec:+90, dist:dist };
    }

    let ra = Math.atan2(pos[1], pos[0]) / (DEG2RAD * 15);
    if (ra < 0) {
        ra += 24;
    }
    let dec = Math.atan2(pos[2], Math.sqrt(xyproj)) / DEG2RAD;
    return { ra:ra, dec:dec, dist:dist };
}

function spin(angle, pos1) {
    const angr = angle * DEG2RAD;
    const cosang = Math.cos(angr);
    const sinang = Math.sin(angr);
    const xx = cosang;
    const yx = sinang;
    const zx = 0;
    const xy = -sinang;
    const yy = cosang;
    const zy = 0;
    const xz = 0;
    const yz = 0;
    const zz = 1;
    let pos2 = [
        xx*pos1[0] + yx*pos1[1] + zx*pos1[2],
        xy*pos1[0] + yy*pos1[1] + zy*pos1[2],
        xz*pos1[0] + yz*pos1[1] + zz*pos1[2]
    ];
    return pos2;
}

function ter2cel(time, vec1) {
    const gast = sidereal_time(time, 'gast');
    let vec2 = spin(-15 * gast, vec1);
    return vec2;
}

function refract(zd_obs, location) {
    if (zd_obs < 0.1 || zd_obs > 91.0) {
        return 0.0;
    }
    let pressure = 1010 * Math.exp(-location.height / 9100);
    let celsius = 10;
    let kelvin = 273 + celsius;
    let h = 90 - zd_obs;
    let angle = (h + 7.31/(h + 4.4));
    let r = 0.016667 / Math.tan(angle * DEG2RAD);
    let refr = r * (0.28 * pressure / kelvin);
    return refr;
}

Astronomy.Horizon = function(date, location, ra, dec, refraction) {     // based on NOVAS equ2hor()
    let time = AstroTime(date);

    const sinlat = Math.sin(location.latitude * DEG2RAD);
    const coslat = Math.cos(location.latitude * DEG2RAD);
    const sinlon = Math.sin(location.longitude * DEG2RAD);
    const coslon = Math.cos(location.longitude * DEG2RAD);
    const sindc = Math.sin(dec * DEG2RAD);
    const cosdc = Math.cos(dec * DEG2RAD);
    const sinra = Math.sin(ra * 15 * DEG2RAD);
    const cosra = Math.cos(ra * 15 * DEG2RAD);
    let uze = [coslat*coslon, coslat*sinlon, sinlat];
    let une = [-sinlat*coslon, -sinlat*sinlon, coslat];
    let uwe = [sinlon, -coslon, 0];

    let uz = ter2cel(time, uze);
    let un = ter2cel(time, une);
    let uw = ter2cel(time, uwe);

    let p = [cosdc*cosra, cosdc*sinra, sindc];

    const pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2];
    const pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2];
    const pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2];

    let proj = Math.sqrt(pn*pn + pw*pw);
    let az = 0;
    if (proj > 0) {
        az = -Math.atan2(pw, pn) * RAD2DEG;
        if (az < 0) az += 360;
        if (az >= 360) az -= 360;
    }
    let zd = Math.atan2(proj, pz) * RAD2DEG;
    let out_ra = ra;
    let out_dec = dec;

    if (refraction) {
        let zd1, refr, j;
        let zd0 = zd;

        if (refraction === 'novas') {
            do {
                zd1 = zd;
                refr = refract(zd, location);
                zd = zd0 - refr;
            } while (Math.abs(zd - zd1) > 3.0e-5);
        } else if (refraction === 'sae' || refraction === 'jplhor') {
            // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
            // JPL Horizons says it uses refraction algorithm from 
            // Meeus "Astronomical Algorithms", 1991, p. 101-102.
            // I found the following Go implementation:
            // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
            // This is a translation from the function "Saemundsson" there.
            let hd = 90 - zd;
            if (refraction === 'jplhor') {
                // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
                // The only reason I'm including this is for a unit test that compares against
                // JPL Horizons data. I recommend most users use the 'sae' option if they want refraction.
                hd = Math.max(-1, hd);
            }
            refr = (1.02 / Math.tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60;
            zd -= refr;
        } else {
            throw 'If specified, refraction must be one of: "novas", "jplhor", "sae".';
        }

        if (refr > 0.0 && zd > 3.0e-4) {
            const sinzd = Math.sin(zd * DEG2RAD);
            const coszd = Math.cos(zd * DEG2RAD);
            const sinzd0 = Math.sin(zd0 * DEG2RAD);
            const coszd0 = Math.cos(zd0 * DEG2RAD);
            var pr = [];
            for (j=0; j<3; ++j) {
                pr.push(((p[j] - coszd0 * uz[j]) / sinzd0)*sinzd + uz[j]*coszd);
            }
            proj = Math.sqrt(pr[0]*pr[0] + pr[1]*pr[1]);
            if (proj > 0) {
                out_ra = Math.atan2(pr[1], pr[0]) * RAD2DEG / 15;
                if (out_ra < 0) {
                    out_ra += 24;
                }
                if (out_ra >= 24) {
                    out_ra -= 24;
                }
            } else {
                out_ra = 0;
            }
            out_dec = Math.atan2(pr[2], proj) * RAD2DEG;
        }
    }

    return { azimuth:az, altitude:90-zd, ra:out_ra, dec:out_dec };
}

function Observer(latitude_degrees, longitude_degrees, height_in_meters) {
    this.latitude = latitude_degrees;
    this.longitude = longitude_degrees;
    this.height = height_in_meters;
}

Astronomy.MakeObserver = function(latitude_degrees, longitude_degrees, height_in_meters) {
    return new Observer(latitude_degrees, longitude_degrees, height_in_meters);
}

Astronomy.SkyPos = function(gc_vector, observer) {     // based on NOVAS place()
    const gc_observer = geo_pos(gc_vector.t, observer);        // vector from geocenter to observer
    const j2000_vector = [
        gc_vector.x - gc_observer[0], 
        gc_vector.y - gc_observer[1], 
        gc_vector.z - gc_observer[2] ];

    let j2000_radec = vector2radec(j2000_vector);

    let pos7 = precession(T0, j2000_vector, gc_vector.t.jd_tt);
    let ofdate_vector = nutation(gc_vector.t.jd_tt, 0, pos7);
    let ofdate_radec = vector2radec(ofdate_vector);

    let sky = {
        t: gc_vector.t,
        j2000: j2000_radec,
        ofdate: ofdate_radec
    }
    return sky;
}

Astronomy.GeoMoon = function(date) {
    var time = AstroTime(date);
    var t = time.tt / 36525;   // t = centuries since J2000.0 epoch
    var moon = CalcMoon(t);

    // Convert geocentric ecliptic spherical coords to cartesian coords.
    var dist_cos_lat = moon.distance_au * Math.cos(moon.geo_eclip_lat);    
    var gepos = [
        dist_cos_lat * Math.cos(moon.geo_eclip_lon),
        dist_cos_lat * Math.sin(moon.geo_eclip_lon),
        moon.distance_au * Math.sin(moon.geo_eclip_lat)
    ];

    // Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
    var mpos1 = ecl2equ_vec(time.jd_tt, gepos);

    // Convert from mean equinox of date to J2000...
    var mpos2 = precession(time.jd_tt, mpos1, T0);

    return { t:time, x:mpos2[0], y:mpos2[1], z:mpos2[2] };
}

function CalcVsop(model, time) {
    var spher = [], eclip, r_coslat;
    var t = time.tt / 365250;   // millennia since 2000
    var formula, series, term, tpower, sum, coord;
    for (formula of model) {
        tpower = 1;
        coord = 0;
        for (series of formula) {
            sum = 0;
            for (term of series) {
                sum += term[0] * Math.cos(term[1] + (t * term[2]));
            }
            coord += tpower * sum;
            tpower *= t;
        }
        spher.push(coord);
    }

    // Convert spherical coordinates to ecliptic cartesian coordinates.
    r_coslat = spher[2] * Math.cos(spher[1]);
    eclip = [
        r_coslat * Math.cos(spher[0]),
        r_coslat * Math.sin(spher[0]),
        spher[2] * Math.sin(spher[1])
    ];

    // Convert ecliptic cartesian coordinates to equatorial cartesian coordinates.
    return {
        t: time,
        x: eclip[0] + 0.000000440360*eclip[1] - 0.000000190919*eclip[2],
        y: -0.000000479966*eclip[0] + 0.917482137087*eclip[1] - 0.397776982902*eclip[2],
        z: 0.397776982902*eclip[1] + 0.917482137087*eclip[2]
    };
}

function ChebScale(t_min, t_max, t) {
    return (2*t - (t_max + t_min)) / (t_max - t_min);
}

function CalcChebyshev(model, time) {
    var record, x, k, d, sum, p0, p1, p2, pos;

    // Search for a record that overlaps the Julian Date 'jd'.
    for (record of model) {
        x = ChebScale(record.jd, record.jd + record.ndays, time.jd_tt);
        if (-1 <= x && x <= +1) {
            pos = [];
            for (d=0; d < 3; ++d) {
                p0 = 1;
                sum = record.coeff[0][d];
                p1 = x;
                sum += record.coeff[1][d] * p1;
                for (k=2; k < record.coeff.length; ++k) {
                    p2 = (2 * x * p1) - p0;
                    sum += record.coeff[k][d] * p2;
                    p0 = p1;
                    p1 = p2;
                }
                pos.push(sum - record.coeff[0][d]/2);
            }
            return { t:time, x:pos[0], y:pos[1], z:pos[2] };
        }
    }
    throw `Cannot extrapolate Chebyshev model for given Terrestrial Time: ${time.jd_tt}`;
}

Astronomy.HelioVector = function(body, date) {
    var time = AstroTime(date);
    if (body in vsop) {
        return CalcVsop(vsop[body], time);
    }
    if (body in cheb) {
        return CalcChebyshev(cheb[body], time);
    }
    if (body === 'Sun') {
        return { t:time, x:0, y:0, z:0 };
    }
    if (body === 'Moon') {
        var e = CalcVsop(vsop.Earth, time);
        var m = Astronomy.GeoMoon(time);
        return { t:time, x:e.x+m.x, y:e.y+m.y, z:e.z+m.z };
    }
    throw "Unknown body";
};

Astronomy.GeoVector = function(body, date) {
    const time = AstroTime(date);
    if (body === 'Moon') {
        return Astronomy.GeoMoon(time);
    }
    if (body === 'Earth') {
        return { t:time, x:0, y:0, z:0 };
    }

    const e = CalcVsop(vsop.Earth, time);

    // Correct for light-travel time, to get position of body as seen from Earth's center.
    let h, geo, ltravel, dt;
    let ltime = time;
    for (let iter=0; iter < 10; ++iter) {
        h = Astronomy.HelioVector(body, ltime);
        geo = { t:time, x:h.x-e.x, y:h.y-e.y, z:h.z-e.z, iter:iter };
        ltravel = Math.sqrt(geo.x*geo.x + geo.y*geo.y + geo.z*geo.z) / C_AUDAY;
        let ltime2 = time.SubtractDays(ltravel);
        dt = Math.abs(ltime2.tt - ltime.tt);
        if (dt < 1.0e-9) {
            return geo;
        }
        ltime = ltime2;
    }
    throw `Light-travel time solver did not converge: dt=${dt}`;
}

})(typeof exports==='undefined' ? (this.Astronomy={}) : exports);
