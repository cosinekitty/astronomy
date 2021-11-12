/*
    pluto_gravsim.h  -  Don Cross <cosinekitty@gmail.com>
    https://github.com/cosinekitty/astronomy

    Constants that tune the gravitational simulation of Pluto's orbit.
*/

#ifndef __ASTRONOMY_PLUTO_GRAVSIM_H
#define __ASTRONOMY_PLUTO_GRAVSIM_H

#define PLUTO_NUM_STATES    51
#define PLUTO_TT1           (-730000)     /* 0001-04-30T12:00:00.000Z */
#define PLUTO_TT2           (+730000)     /* 3998-09-03T12:00:00.000Z */

#if ((PLUTO_TT2 - PLUTO_TT1) % (PLUTO_NUM_STATES - 1)) != 0
    #error PLUTO_TIME_STEP ratio must be an integer.
#endif

#define PLUTO_TIME_STEP     ((PLUTO_TT2 - PLUTO_TT1) / (PLUTO_NUM_STATES - 1))
#define PLUTO_DT            146

#if PLUTO_TIME_STEP % PLUTO_DT != 0
    #error Invalid combination of Pluto time step, time increment.
#endif

#define PLUTO_NSTEPS    ((PLUTO_TIME_STEP / PLUTO_DT) + 1)

#endif
