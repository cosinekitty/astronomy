/*
    vsop.h  -  Don Cross  -  2019-03-24

    Loads and calculates planetary positions using VSOP87 analytic models.
    See:  https://en.wikipedia.org/wiki/VSOP_(planets)
*/
#ifndef __DDC_VSOP_H
#define __DDC_VSOP_H

typedef enum
{
    VSOP_INVALID_VERSION = -1,
    VSOP_ELLIPTIC_J2000 = 0,        /* VSOP87.*  files */
    VSOP_HELIO_RECT_J2000 = 1,      /* VSOP87A.* files */
    VSOP_HELIO_SPHER_J2000 = 2,     /* VSOP87B.* files */
    VSOP_HELIO_RECT_DATE = 3,       /* VSOP87C.* files */
    VSOP_HELIO_SPHER_DATE = 4,      /* VSOP87D.* files */
    VSOP_BARY_RECT_J2000 = 5,       /* VSOP87E.* files */
}
vsop_version_t;

typedef struct 
{
    double amplitude;
    double phase;
    double frequency;
}
vsop_term_t;

typedef struct 
{
    int nterms_total;       /* total number of terms in the VSOP87 model */
    int nterms_calc;        /* number of terms at the front to actually calculate (allows faster, lower-res calc) */
    vsop_term_t *term; 
}
vsop_series_t;

#define VSOP_MAX_SERIES 10

typedef struct
{
    int nseries_total;      /* total number of series in the VSOP87 model */
    int nseries_calc;       /* number of series at the front to actually calculate (allows faster, lower-rest calc) */
    vsop_series_t series[VSOP_MAX_SERIES];
}
vsop_formula_t;

#define VSOP_MIN_COORDS 3
#define VSOP_MAX_COORDS 6

typedef enum    /* values created for compatibility with NOVAS; these are *NOT* the body codes used by VSOP! */
{
    VSOP_INVALID_BODY = -1,
    VSOP_MERCURY =  0,
    VSOP_VENUS   =  1,
    VSOP_EMB     =  2,       /* Earth/Moon barycenter, to match NOVAS */
    VSOP_MARS    =  3,
    VSOP_JUPITER =  4,
    VSOP_SATURN  =  5,
    VSOP_URANUS  =  6,
    VSOP_NEPTUNE =  7,
    // NOTE: VSOP does not provide Pluto or Moon.
    VSOP_SUN     = 10,
    VSOP_EARTH   = 11,      /* weird value so as not to conflict with NOVAS body values */
}
vsop_body_t;

#define VSOP_BODY_LIMIT 12

typedef struct 
{
    vsop_version_t version;
    vsop_body_t body;
    int ncoords;
    vsop_formula_t formula[VSOP_MAX_COORDS];
}
vsop_model_t;

void VsopInit(vsop_model_t *model);    /* optional usage: create invalid null model that can be safely freed */
void VsopFreeModel(vsop_model_t *model);

int VsopLoadModel(vsop_model_t *model, const char *inFileName);
int VsopCalc(const vsop_model_t *model, double jd, double rect_j2000_equatorial_pos[3]);
int VsopTruncate(vsop_model_t *model, double jd1, double jd2, double amplitudeThreshold);
int VsopTermCount(const vsop_model_t *model);
int VsopWriteTrunc(const vsop_model_t *model, const char *outFileName);
int VsopReadTrunc(vsop_model_t *model, const char *inFileName);

#endif /* __DDC_VSOP_H */
