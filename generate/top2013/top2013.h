/*
    top2013.h  -  Don Cross  -  2020-06-30

    Loads and calculates planetary positions using TOP2013 analytic models.
    See:  https://en.wikipedia.org/wiki/VSOP_(planets)
*/
#ifndef __DDC_TOP2013_H
#define __DDC_TOP2013_H

typedef struct
{
    double  k;      /* mu coefficient */
    double  c;      /* cosine coefficient */
    double  s;      /* sine coefficient */
    double  p;      /* period expressed in years */
    int     rc;     /* hack: final digit rounding adjustment for 'c' */
    int     rs;     /* hack: final digit rounding adjustment for 's' */
}
top_term_t;


typedef struct
{
    int nterms_total;
    int nterms_calc;
    top_term_t *terms;
}
top_series_t;

#define TOP_NCOORDS  6
#define TOP_NSERIES 13

typedef struct
{
    int planet;     /* 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune, 9=Pluto */
    top_series_t series[TOP_NCOORDS][TOP_NSERIES];
}
top_model_t;


typedef struct
{
    double a;           /* AU */
    double lambda;      /* rad */
    double k;           /* 1 */
    double h;           /* 1 */
    double q;           /* 1 */
    double p;           /* 1 */
}
top_elliptical_t;


typedef struct
{
    /* Position in AU */
    double  x;
    double  y;
    double  z;

    /* velocity in AU/day */
    double  vx;
    double  vy;
    double  vz;
}
top_rectangular_t;


typedef struct
{
    /* a sortable item that estimates how much a given term contributes to an elliptical element */
    double  magnitude;
    int     s;
    int     t;
}
top_contrib_t;


typedef struct
{
    int             nterms;
    top_contrib_t  *array;
}
top_contrib_list_t;


typedef struct
{
    /* an independent sortable list for each of the elliptical formulas */
    top_contrib_list_t list[TOP_NCOORDS];
}
top_contrib_map_t;


void TopInitModel(top_model_t *model);
void TopFreeModel(top_model_t *model);
int  TopCloneModel(top_model_t *copy, const top_model_t* original);
int  TopLoadModel(top_model_t *model, const char *filename, int planet);
int  TopSaveModel(const top_model_t *model, const char *filename);
int  TopWriteModel(const top_model_t *model, FILE *outfile);
int  TopCalcElliptical(const top_model_t *model, double tt, top_elliptical_t *ellip);
int  TopEcliptic(int planet, const top_elliptical_t *ellip, top_rectangular_t *ecl);
int  TopEquatorial(const top_rectangular_t *ecl, top_rectangular_t *equ);

void TopInitContribMap(top_contrib_map_t *map);
int  TopMakeContribMap(const top_model_t *model, top_contrib_map_t *map, double millennia);
void TopFreeContribMap(top_contrib_map_t *map);

void TopSquash(
    top_model_t *copy,
    const top_model_t *original,
    const top_contrib_map_t *map,
    double amplitude[TOP_NCOORDS]);

#endif /* __DDC_TOP2013_H */
