#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// general helper function declaration {{

// returns a random number between 0.0 and fMaxValue
double randomNumber(double fMaxValue);

// returns a random angle between 0.0 and 2*pi
double randomAngle();

// set the random number seed, if pnSeed != NULL set *pnSeed to the seed used
// return TRUE at success, otherwise return FALSE
BOOL randomSeed(unsigned int *pnSeed = NULL);

// returns the square of fValue
double square(double fValue);

// }} general helper function declaration

// general helper function definition {{

double randomNumber(double fMaxValue)
{
	return fMaxValue*((double)rand()/((double)RAND_MAX));
}

double randomAngle()
{
    return randomNumber(2.0*M_PI);
}

BOOL randomSeed(unsigned int *pnSeed)
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

double square(double fValue)
{
	return fValue * fValue;
}

// }} general helper function definition

// CWalkerCellBounds {{

// class that represents the bounds of a cell
class CWalkerCellBounds
{
	// constructor
public:
	CWalkerCellBounds();

	// methods
	double height() const;
	double width() const;
	double area() const;

	// data
public:
	double m_fxmin;
	double m_fxmax;
	double m_fymin;
	double m_fymax;
};

// constructor
CWalkerCellBounds::CWalkerCellBounds()
{
	m_fxmin = 0.0;
	m_fxmax = 0.0;
	m_fymin = 0.0;
	m_fymax = 0.0;
}

double CWalkerCellBounds::height() const
{
	return m_fymax - m_fymin;
}

double CWalkerCellBounds::width() const
{
	return m_fxmax - m_fxmin;
}

double CWalkerCellBounds::area() const
{
    return height()*width();
}

// }} CWalkerCellBounds

// CWalkerCell {{

// class that represent the cell position of a walker
class CWalkerCell
{
	// constructor
public:
	CWalkerCell();

public:
	int m_nx;
	int m_ny;
};

// constructor
CWalkerCell::CWalkerCell()
{
	m_nx = 0;
	m_ny = 0;
}

// }} CWalkerCell

// class that represents a point in two dimensions
class CWalkerPoint
{
	// constructor;
public:
	CWalkerPoint();

	double squareDistanceTo(const CWalkerPoint&) const;

	// data
public:
	double m_fx;
	double m_fy;
};

double CWalkerPoint::squareDistanceTo(const CWalkerPoint& otherPoint) const
{
	return (square(otherPoint.m_fx - m_fx) + square(otherPoint.m_fy - m_fy));
}

// constructor
CWalkerPoint::CWalkerPoint()
{
	m_fx = 0.0;
	m_fx = 0.0;
}

// CWalkerCircularBarrier implementation {{

class CWalkerCircularBarrier
{
public:
	// constructor
	CWalkerCircularBarrier();

	// this method iterate with a moving point, it will set pointPosEnd acording to the iteration, or return
	// FALSE if it is an invalid point.
	// BOOL iterateWithMovingPoint(const CWalkerPoint& pointPosIni, CWalkerPoint& pointPosEnd) const;

	BOOL isInside(const CWalkerPoint&) const;

public:
	CWalkerPoint m_position;
	double m_radius;
};

CWalkerCircularBarrier::CWalkerCircularBarrier()
{
	m_radius = 1.0;
}

BOOL CWalkerCircularBarrier::isInside(const CWalkerPoint& point) const
{
	return (point.squareDistanceTo(m_position) < square(m_radius));
}

// }} end of CWalkerCircularBarrier implementation

// CWalker implementation {{

// class that represent each walker
class CWalker
{
public:
	CWalker();

	// return the square distance from the origin for this walker
	double squareDistanceFromOrigin(const CWalkerCellBounds* pcCellBounds = NULL) const;

	// iterate the walker with nSteps, and calls manageCellPosition if pcCellBounds != NULL (for each step
	void iterate(int nSteps = 1, const CWalkerCellBounds* pcCellBounds = NULL);

	// manage the cell position if the walker point is outside bounds
	void manageCellPosition(const CWalkerCellBounds& cellBounds);

public:
	CWalkerPoint m_origin;
	CWalkerPoint m_position;
	CWalkerCell m_cell;
};

CWalker::CWalker()
{
}

double CWalker::squareDistanceFromOrigin(const CWalkerCellBounds* pcCellBounds) const
{
	// TODO still not usin the cell

	double fDeltaX = m_position.m_fx - m_origin.m_fx;
	double fDeltaY = m_position.m_fy - m_origin.m_fy;

	if(pcCellBounds)
	{
		fDeltaX += pcCellBounds->width()*((double)m_cell.m_nx);
		fDeltaY += pcCellBounds->height()*((double)m_cell.m_ny);
	}

	return square(fDeltaX) + square(fDeltaY);
}

void CWalker::iterate(int nSteps, const CWalkerCellBounds* pcCellBounds)
{
	int i;

	for(i = 0; i < nSteps; i++)
	{
		double fRandomAngle = randomAngle();
		m_position.m_fx += cos(fRandomAngle);
		m_position.m_fy += sin(fRandomAngle);

		if(pcCellBounds)
		{
			manageCellPosition(*pcCellBounds);
		}
	}
}

void CWalker::manageCellPosition(const CWalkerCellBounds& cellBounds)
{
	if(m_position.m_fx < cellBounds.m_fxmin)
	{
		m_position.m_fx += cellBounds.width();
		m_cell.m_nx--;
	}

	if(m_position.m_fx > cellBounds.m_fxmax)
	{
		m_position.m_fx -= cellBounds.width();
		m_cell.m_nx++;
	}

	if(m_position.m_fy < cellBounds.m_fymin)
	{
		m_position.m_fy += cellBounds.width();
		m_cell.m_ny--;
	}

	if(m_position.m_fy > cellBounds.m_fymax)
	{
		m_position.m_fy -= cellBounds.width();
		m_cell.m_ny++;
	}
}

// }} end of CWalker implementation

// verify if point lies inside one of the barriers of the array pcBarriers
BOOL isInsideBarrier(const CWalkerPoint& point, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers);

// do the walker iteration
// returns TRUE if the walker interacted with the barriers
// {{

// simple iteration, ignore barriers
BOOL iterateWalkerWithBarriers0(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds*);

// iterate until the walker is in a valid position (reset the position at every attempt).
BOOL iterateWalkerWithBarriers1(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds*);

// try to move the walker, if theres a berrir blocking resets the position and returns.
BOOL iterateWalkerWithBarriers2(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds*);

// reflects the walker at the barriers
BOOL iterateWalkerWithBarriers3(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds*);

// }}

BOOL isInsideBarrier(const CWalkerPoint& point, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers)
{
    unsigned int i = 0;
    for(i = 0; i < nBarriers; i++)
    {
        if(pcBarrier[i].isInside(point))
        {
            return TRUE;
        }
    }

    return FALSE;
}

BOOL iterateWalkerWithBarriers0(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds* pCellBounds)
{
    walker.iterate(1, pCellBounds);
    return FALSE;
}

BOOL iterateWalkerWithBarriers1(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds* pCellBounds)
{
    CWalker oldWalker = walker;
    BOOL bOutOfBarriers;
    BOOL bBouncedOnBarrier = FALSE;

    // itetates until it falls out of a barrier
    do
    {
        walker.iterate(1, pCellBounds);
        // fprintf(stdout, "# AfterIterate. j:%d, i:%d, wp:%f, %f\n", j, i, pWalkerArray[i].m_position.m_fx, pWalkerArray[i].m_position.m_fy);

        bOutOfBarriers = TRUE;

        if(isInsideBarrier(walker.m_position, pcBarrier, nBarriers))
        {
            bOutOfBarriers = FALSE;
            bBouncedOnBarrier = TRUE;
            walker = oldWalker;
        }
    }
    while(!bOutOfBarriers);

    return bBouncedOnBarrier;
}

BOOL iterateWalkerWithBarriers2(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds* pCellBounds)
{
    CWalker oldWalker = walker;

    walker.iterate(1, pCellBounds);

    if(isInsideBarrier(walker.m_position, pcBarrier, nBarriers))
    {
        walker = oldWalker;

        return TRUE;
    }

    return FALSE;
}

BOOL iterateWalkerWithBarriers3(CWalker &walker, const CWalkerCircularBarrier* pcBarrier, unsigned int nBarriers, const CWalkerCellBounds* pCellBounds)
{
    // these will need some math...


    return FALSE;
}

int main(int argc, const char* argv[])
{
	CWalker *pWalkerArray = NULL;
	CWalkerCircularBarrier *pWalkerBarrierArray = NULL;
	CWalkerCellBounds cellBounds;
	CWalkerPoint pntOrigin;
    unsigned int nWalkers;
    unsigned int nIterations;
    unsigned int nBarriers;
    unsigned int nInterferenceCounter;
    unsigned int nBarrierMaxPlacementAttempts;
    unsigned int nBarrierPlacementAttempts;
    unsigned int nBarrierActionMethod;
    unsigned int nSeed;
    unsigned int i, j; //, k; // iteration steps
    double fBarrierRadius;
    double fOccupiedArea;
    double fOccupiedAreaRatio;
	FILE *pStatisticsOutStream = stdout;
	FILE *pWalkersOutStream = NULL;
	FILE *pBarriersOutStream = NULL;

	// default parameters
	nBarrierActionMethod = 1;
	nBarriers = 200;
    nWalkers = 1000;
    nIterations = 1000;
	cellBounds.m_fxmax = 100.0;
	cellBounds.m_fymax= 100.0;
	fBarrierRadius = 1.0;
	fOccupiedArea = 0.0;
	nSeed = RAND_MAX;
	// these is a guess, needs a little calculation to be more accurate, currently based on attempts with nBarriers = 400 and fBarrierRadius = 2.0
	nBarrierMaxPlacementAttempts = 100000;

	// read the command line arguments
	for(i = 0; i < (unsigned int)argc; i++)
	{
		if(argv[i][0] != '-')
			continue;

		if(!strcmp(argv[i]+1, "w")) // number of walkers
		{
			nWalkers = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "i")) // number of iterations
		{
			nIterations = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "b")) // number of barriers
		{
			nBarriers = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "br")) // barrier radius
		{
			fBarrierRadius = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "m")) // barrier action method
		{
			nBarrierActionMethod = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "s")) // new seed
		{
		    // here the whole range can be requested
		    char *pDummy;
			nSeed = strtoul(argv[i+1], &pDummy, 10);
		}
		else if(!strcmp(argv[i]+1, "sof")) // out file for statistics data such as <r2>, <r>, <x>, <y>, etc...
		{
            pStatisticsOutStream = fopen(argv[i+1], "wt");

            if(!pStatisticsOutStream)
            {
                fprintf(stderr, "Failed to open file for statistics data storage.\n");
                return 0;
            }
		}
		else if(!strcmp(argv[i]+1, "bof")) // out file for barrier position...
		{
            pBarriersOutStream = fopen(argv[i+1], "wt");

            if(!pBarriersOutStream)
            {
                fprintf(stderr, "Failed to open file for barriers data storage.\n");
                return 0;
            }
		}
    }

	// out vars
	nInterferenceCounter = 0;

    if(nSeed == RAND_MAX)
    {
        // randomize random seed
        randomSeed(&nSeed);
    }

	// print some information about the simulation
	fprintf(pStatisticsOutStream, "# data generated by random walker with barriers\n");
	fprintf(pStatisticsOutStream, "# by guilherme dellagustin\n");
	fprintf(pStatisticsOutStream, "# random seed: %u\n", nSeed);
	fprintf(pStatisticsOutStream, "# walkers: %d\n", nWalkers);
	fprintf(pStatisticsOutStream, "# iterations: %d\n", nIterations);
	fprintf(pStatisticsOutStream, "# barriers: %d\n", nBarriers);
	fprintf(pStatisticsOutStream, "# barriers radius: %f\n", fBarrierRadius);
	fprintf(pStatisticsOutStream, "# min bounds: %f, %f\n", cellBounds.m_fxmin, cellBounds.m_fymin);
	fprintf(pStatisticsOutStream, "# max bounds: %f, %f\n", cellBounds.m_fxmax, cellBounds.m_fymax);
	fprintf(pStatisticsOutStream, "# barrier action method: %d\n", nBarrierActionMethod);

	// allocate memory
	if(nWalkers)
		pWalkerArray = new CWalker[nWalkers];

	if(nBarriers)
		pWalkerBarrierArray = new CWalkerCircularBarrier[nBarriers];

	// initiate data

	// place the barriers {{

	nBarrierPlacementAttempts = 0;

	for(i = 0; i < nBarriers; i++)
	{
	    BOOL bHasInterference;

        j = 0;

        do
	    {
	        // just some extra info...
	        nBarrierPlacementAttempts ++;

            bHasInterference = FALSE;

	        pWalkerBarrierArray[i].m_position.m_fx = cellBounds.m_fxmin + randomNumber(cellBounds.width());
            pWalkerBarrierArray[i].m_position.m_fy = cellBounds.m_fymin + randomNumber(cellBounds.height());
            pWalkerBarrierArray[i].m_radius = fBarrierRadius;

            // avoid barriers touching the cell walls
            if(
            pWalkerBarrierArray[i].m_position.m_fx < cellBounds.m_fxmin + pWalkerBarrierArray[i].m_radius ||
            pWalkerBarrierArray[i].m_position.m_fx > cellBounds.m_fxmax - pWalkerBarrierArray[i].m_radius ||
            pWalkerBarrierArray[i].m_position.m_fy < cellBounds.m_fymin + pWalkerBarrierArray[i].m_radius ||
            pWalkerBarrierArray[i].m_position.m_fy > cellBounds.m_fymax - pWalkerBarrierArray[i].m_radius
            )
            {
                bHasInterference = TRUE;
            }
            else
            {
                // test if the barriers are overllaped
                for(j = 0; j < i; j++)
                {
                    if(pWalkerBarrierArray[j].m_position.squareDistanceTo(pWalkerBarrierArray[i].m_position) < square(pWalkerBarrierArray[j].m_radius + pWalkerBarrierArray[i].m_radius))
                    {
                        bHasInterference = TRUE;
                        break;
                    }
                }
            }
	    }
	    while(bHasInterference && j < nBarrierMaxPlacementAttempts);
		// fprintf(stdout, "%f\t%f\n", pWalkerBarrierArray[i].m_position.m_fx, pWalkerBarrierArray[i].m_position.m_fy);

		fOccupiedArea += M_PI*square(pWalkerBarrierArray[i].m_radius);

		if(pBarriersOutStream)
		{
		    fprintf(pBarriersOutStream, "%f\t%f\t%f\n", pWalkerBarrierArray[j].m_position.m_fx, pWalkerBarrierArray[j].m_position.m_fy, pWalkerBarrierArray[j].m_radius);
		}
	}

    fOccupiedAreaRatio = fOccupiedArea/(cellBounds.area());
	fprintf(pStatisticsOutStream, "# attempts to place barriers: %d\n", nBarrierPlacementAttempts);
	fprintf(pStatisticsOutStream, "# ocupied area ratio: %f\n", fOccupiedAreaRatio);

	// }} place the barriers

	do
	{
		pntOrigin.m_fx = cellBounds.m_fxmin + randomNumber(cellBounds.width());
		pntOrigin.m_fy = cellBounds.m_fymin + randomNumber(cellBounds.height());
	}
	while(isInsideBarrier(pntOrigin, pWalkerBarrierArray, nBarriers));

	for(i = 0; i < nWalkers; i++)
	{
		pWalkerArray[i].m_origin = pntOrigin;
		pWalkerArray[i].m_position.m_fx = pWalkerArray[i].m_origin.m_fx;
		pWalkerArray[i].m_position.m_fy = pWalkerArray[i].m_origin.m_fy;
    }

    // iteration in 'time'
	for(j = 0; j < nIterations; j++)
    {
		double fAverageSquareDistance = 0.0;
		double fAverageX = 0.0;
		double fAverageY = 0.0;
		unsigned int nOldInterferenceCounter = nInterferenceCounter;

		// iteration in walkers
		for(i = 0 ; i < nWalkers; i++)
		{
			BOOL bBouncedOnBarrier = FALSE;

			switch(nBarrierActionMethod)
			{
            default:
            case 0:
                bBouncedOnBarrier = iterateWalkerWithBarriers0(pWalkerArray[i], pWalkerBarrierArray, nBarriers, &cellBounds);
                break;
            case 1:
                bBouncedOnBarrier = iterateWalkerWithBarriers1(pWalkerArray[i], pWalkerBarrierArray, nBarriers, &cellBounds);
                break;
            case 2:
                bBouncedOnBarrier = iterateWalkerWithBarriers2(pWalkerArray[i], pWalkerBarrierArray, nBarriers, &cellBounds);
                break;
            case 3:
                bBouncedOnBarrier = iterateWalkerWithBarriers3(pWalkerArray[i], pWalkerBarrierArray, nBarriers, &cellBounds);
                break;
			}

			if(bBouncedOnBarrier)
				nInterferenceCounter++;

			fAverageSquareDistance += pWalkerArray[i].squareDistanceFromOrigin(&cellBounds);
			fAverageX += pWalkerArray[i].m_position.m_fx;
			fAverageY += pWalkerArray[i].m_position.m_fy;
		}

		fAverageSquareDistance /= (double)nWalkers;
		fAverageX /= (double)nWalkers;
		fAverageY /= (double)nWalkers;
		fprintf(pStatisticsOutStream, "%d\t%f\t%f\t%f\t\%u\t%u\t%f\n", j, fAverageSquareDistance, fAverageX, fAverageY, nInterferenceCounter, nInterferenceCounter - nOldInterferenceCounter, fAverageSquareDistance*fOccupiedAreaRatio);
    }

    if(pStatisticsOutStream && pStatisticsOutStream != stdout)
    {
        fclose(pStatisticsOutStream);
        pStatisticsOutStream = NULL;
    }

    if(pBarriersOutStream && pBarriersOutStream != stdout)
    {
        fclose(pBarriersOutStream);
        pBarriersOutStream = NULL;
    }

    if(pWalkersOutStream && pWalkersOutStream != stdout)
    {
        fclose(pWalkersOutStream);
        pWalkersOutStream = NULL;
    }

	// deallocate memory
	delete [] pWalkerArray ;

	return 1;
}
