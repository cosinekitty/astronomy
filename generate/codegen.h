/*
    MIT License

    Copyright (c) 2019-2021 Don Cross <cosinekitty@gmail.com>

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
#ifndef __DDC_ASTRO_CODEGEN_H
#define __DDC_ASTRO_CODEGEN_H

#define CHECK(x)    do{if(0 != (error = (x))) goto fail;}while(0)
#define FAIL(...)   do{fprintf(stderr, __VA_ARGS__); error = 1; goto fail;}while(0)

typedef enum
{
    CODEGEN_LANGUAGE_UNKNOWN,
    CODEGEN_LANGUAGE_JS,
    CODEGEN_LANGUAGE_C,
    CODEGEN_LANGUAGE_CSHARP,
    CODEGEN_LANGUAGE_PYTHON
}
cg_language_t;

int GenerateCode(
    cg_language_t language,
    const char *outCodeFileName,
    const char *inTemplateFileName,
    const char *dataPath);

double ExtrapolatedDeltaT(int year);

#endif /* __DDC_ASTRO_CODEGEN_H */
