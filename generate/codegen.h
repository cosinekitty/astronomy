/*
    codegen.h  -  Don Cross  -  2019-03-25

    Generate astronomy source code optimized for size
    and within 1-arcminute precision over the years 1900..2100.
*/
#ifndef __DDC_ASTRO_CODEGEN_H
#define __DDC_ASTRO_CODEGEN_H

int GenerateCode(
    const char *outCodeFileName,
    const char *inTemplateFileName,
    const char *dataPath);

#endif /* __DDC_ASTRO_CODEGEN_H */
