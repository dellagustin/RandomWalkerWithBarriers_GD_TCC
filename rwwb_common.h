#ifndef RWWB_COMMON_H_INCLUDED
#define RWWB_COMMON_H_INCLUDED

#include <stdio.h>
#include <math.h>

#ifndef BOOL
typedef unsigned int BOOL;
#endif

#ifndef TRUE
#define TRUE (-1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

#ifndef __min
#define __min(arg1, arg2) ((arg1)<(arg2)?arg1:arg2)
#endif

// returns a random number between 0.0 and fMaxValue
double randomNumber(double fMaxValue);

// returns a random angle between 0.0 and 2*pi
double randomAngle();

// set the random number seed, if pnSeed != NULL set *pnSeed to the seed used
// return TRUE at success, otherwise return FALSE
BOOL randomSeed(unsigned int *pnSeed = NULL);

// returns the square of fValue
double square(double fValue);


// definitions

inline double randomNumber(double fMaxValue)
{
	return fMaxValue*((double)rand()/((double)RAND_MAX));
}

inline double randomAngle()
{
    return randomNumber(2.0*M_PI);
}

inline BOOL randomSeed(unsigned int *pnSeed)
{
    FILE *pRandomFile;
    unsigned int nSeed;

    pRandomFile = fopen("/dev/urandom", "rb");

    if(pRandomFile)
    {
		fread(&nSeed, sizeof(nSeed), 1, pRandomFile);
        fclose(pRandomFile);

        if(pnSeed)
        {
            *pnSeed = nSeed;
        }

        srand(nSeed);

        return TRUE;
    }

    return FALSE;
}

inline double square(double fValue)
{
	return fValue * fValue;
}

#endif // RWWB_COMMON_H_INCLUDED
