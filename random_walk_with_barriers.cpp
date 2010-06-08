#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "rwwb_common.h"

// general helper function declaration {{

void close_stream(FILE*&);

// }} general helper function declaration

// general helper function definition {{

void close_stream(FILE*& pStream)
{
	if(pStream && pStream != stdout && pStream != stderr)
    {
        fclose(pStream);
    }

    pStream = NULL;
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
	m_fy = 0.0;
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

	// return the square distance from the origin for this walker
	double distanceFromOrigin(const CWalkerCellBounds* pcCellBounds = NULL) const;

	// iterate the walker with nSteps, and calls manageCellPosition if pcCellBounds != NULL (for each step
	void iterate(int nSteps = 1, const CWalkerCellBounds* pcCellBounds = NULL);

	// manage the cell position if the walker point is outside bounds
	void manageCellPosition(const CWalkerCellBounds& cellBounds);

	CWalkerPoint realPosition(const CWalkerCellBounds& cellBounds) const;

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
	CWalkerPoint position;

	position = pcCellBounds ? realPosition(*pcCellBounds) : m_position;

	double fDeltaX = position.m_fx - m_origin.m_fx;
	double fDeltaY = position.m_fy - m_origin.m_fy;

	return square(fDeltaX) + square(fDeltaY);
}

double CWalker::distanceFromOrigin(const CWalkerCellBounds* pcCellBounds) const
{
	return sqrt(squareDistanceFromOrigin(pcCellBounds));
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

CWalkerPoint CWalker::realPosition(const CWalkerCellBounds& cellBounds) const
{
	CWalkerPoint retPoint;

	retPoint = m_position;

	retPoint.m_fx += ((double)m_cell.m_nx)*cellBounds.width();
	retPoint.m_fy += ((double)m_cell.m_ny)*cellBounds.height();

	return retPoint;
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

void report_walker(unsigned int nTime, unsigned int nWalker, const CWalker* pWalkerArray, const CWalkerCellBounds& cellBounds, FILE *pStream)
{
	CWalkerPoint realPosition;
	realPosition = pWalkerArray[nWalker].realPosition(cellBounds);

	fprintf(pStream, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
		nTime, nWalker,
		realPosition.m_fx, realPosition.m_fy,
		pWalkerArray[nWalker].m_position.m_fx, pWalkerArray[nWalker].m_position.m_fy,
		pWalkerArray[nWalker].m_origin.m_fx, pWalkerArray[nWalker].m_origin.m_fy,
		pWalkerArray[nWalker].squareDistanceFromOrigin(&cellBounds),
		pWalkerArray[nWalker].distanceFromOrigin(&cellBounds));
}

int main(int argc, const char* argv[])
{
	CWalker *pWalkerArray = NULL;
	CWalkerCircularBarrier *pWalkerBarrierArray = NULL;
	CWalkerCellBounds cellBounds;
	CWalkerPoint pntOrigin;
    unsigned int nWalkers;
    unsigned int nIterations;
    unsigned int nFirstIterationToPrint;
    unsigned int nBarriers;
    unsigned int nInterferenceCounter;
    unsigned int nBarrierMaxPlacementAttempts;
    unsigned int nBarrierPlacementAttempts;
    unsigned int nBarrierActionMethod;
    unsigned int nSeed;
    unsigned int nStepsToShowErr;
    unsigned int i, j, k; // iteration steps
    double fBarrierRadius;
    double fOccupiedArea;
    double fOccupiedAreaRatio;
    double fTotalBarrierPerimeter;
    double fDesiredOcupiedArea;
    FILE *pStatisticsOutStream = stdout;
	FILE *pWalkersOutStream = NULL;
	FILE *pWalkersStartOutStream = NULL;
	FILE *pWalkersEndOutStream = NULL;
	FILE *pBarriersOutStream = NULL;
	struct timeval time_data_start;
	struct timeval time_data_current;
	struct timeval time_data_old;

	// read the start time
	gettimeofday(&time_data_start, NULL);

	// default parameters
	nBarrierActionMethod = 1;
	nBarriers = 200;
    nWalkers = 1000;
    nIterations = 1000;
    nFirstIterationToPrint = 0;
	cellBounds.m_fxmax = 100.0;
	cellBounds.m_fymax= 100.0;
	fBarrierRadius = 1.0;
	fOccupiedArea = 0.0;
	fDesiredOcupiedArea = -1.0;
	fTotalBarrierPerimeter = 0.0;
	nSeed = RAND_MAX;
	// these is a guess, needs a little calculation to be more accurate, currently based on attempts with nBarriers = 400 and fBarrierRadius = 2.0
	nBarrierMaxPlacementAttempts = 1000000;

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
                fprintf(stderr, "# Failed to open file for statistics data storage.\n");
                return 0;
            }
		}
		else if(!strcmp(argv[i]+1, "bof")) // out file for barrier position...
		{
            pBarriersOutStream = fopen(argv[i+1], "wt");

            if(!pBarriersOutStream)
            {
                fprintf(stderr, "# Failed to open file for barriers data storage.\n");
                return 0;
            }
		}
		else if(!strcmp(argv[i]+1, "wsof")) // out file for barrier position...
		{
            pWalkersStartOutStream = fopen(argv[i+1], "wt");

            if(!pWalkersStartOutStream)
            {
                fprintf(stderr, "# Failed to open file for walkers starting position data storage.\n");
                return 0;
            }
		}
		else if(!strcmp(argv[i]+1, "weof")) // out file for barrier position...
		{
            pWalkersEndOutStream = fopen(argv[i+1], "wt");

            if(!pWalkersEndOutStream)
            {
                fprintf(stderr, "# Failed to open file for walkers starting position data storage.\n");
                return 0;
            }
		}
		else if(!strcmp(argv[i]+1, "wof")) // out file for barrier position...
		{
            pWalkersOutStream = fopen(argv[i+1], "wt");

            if(!pWalkersOutStream)
            {
                fprintf(stderr, "# Failed to open file for walkers position data storage.\n");
                return 0;
            }
		}
    }

    nStepsToShowErr = nIterations/10;

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

	gettimeofday(&time_data_old, NULL);
	fprintf(stderr, "# Placing %u Barriers...\n", nBarriers);

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
	    while(bHasInterference && nBarrierPlacementAttempts < nBarrierMaxPlacementAttempts);

		if(nBarrierPlacementAttempts >= nBarrierMaxPlacementAttempts)
		{
            fprintf(stderr, "# Maximum attemps to place barriers reached: %u.\n", nBarrierPlacementAttempts);
            // break;
            goto clean_up;
		}
		else
		{
			fOccupiedArea += M_PI*square(pWalkerBarrierArray[i].m_radius);
			fTotalBarrierPerimeter += 2.0*M_PI*pWalkerBarrierArray[i].m_radius;

			if(pBarriersOutStream)
			{
				fprintf(pBarriersOutStream, "%f\t%f\t%f\n", pWalkerBarrierArray[j].m_position.m_fx, pWalkerBarrierArray[j].m_position.m_fy, pWalkerBarrierArray[j].m_radius);

				// desperate measures
				for(k = 0; k < 100; k++)
				{
					fprintf(pBarriersOutStream, "%f\t%f\t%f\n",
					pWalkerBarrierArray[j].m_position.m_fx + pWalkerBarrierArray[j].m_radius * cos(2.0*M_PI*((double)k)/100.0),
					pWalkerBarrierArray[j].m_position.m_fy + pWalkerBarrierArray[j].m_radius * sin(2.0*M_PI*((double)k)/100.0),
					pWalkerBarrierArray[j].m_radius);
				}
			}
		}
	}

	close_stream(pBarriersOutStream);

	gettimeofday(&time_data_current, NULL);
	fprintf(stderr, "# %u Barriers Placed in %u seconds.\n", nBarriers, (unsigned int)(time_data_current.tv_sec - time_data_old.tv_sec));

    fOccupiedAreaRatio = fOccupiedArea/(cellBounds.area());
	fprintf(pStatisticsOutStream, "# attempts to place barriers: %d\n", nBarrierPlacementAttempts);
	fprintf(pStatisticsOutStream, "# ocupied area ratio: %f\n", fOccupiedAreaRatio);
	fprintf(pStatisticsOutStream, "# ocupied area: %f\n", fOccupiedArea);
	fprintf(pStatisticsOutStream, "# total barrier perimeter: %f\n", fTotalBarrierPerimeter);
	fprintf(pStatisticsOutStream, "# average barrier perimeter: %f\n", fTotalBarrierPerimeter/((double)nBarriers));

	// }} place the barriers

	// place the walkers {{

	gettimeofday(&time_data_old, NULL);
    fprintf(stderr, "# Placing %u Walkers.\n", nWalkers);

    if(pWalkersStartOutStream)
	{
		fprintf(pWalkersStartOutStream, "# Walkers Starting Position File\n");
		fprintf(pWalkersStartOutStream, "# Walker\tx\ty\n");
	}

	if(pWalkersOutStream)
	{
		fprintf(pWalkersOutStream, "# Walkers Position File\n");
		fprintf(pWalkersOutStream, "# t\tWalker\tx\ty\txo\tyo\tr2\n");
	}

	for(i = 0; i < nWalkers; i++)
	{
		do
		{
			pWalkerArray[i].m_origin.m_fx = cellBounds.m_fxmin + randomNumber(cellBounds.width());
			pWalkerArray[i].m_origin.m_fy = cellBounds.m_fymin + randomNumber(cellBounds.height());
		}
		while(isInsideBarrier(pWalkerArray[i].m_origin, pWalkerBarrierArray, nBarriers));

		pWalkerArray[i].m_position = pWalkerArray[i].m_origin;

		if(pWalkersStartOutStream)
		{
			report_walker(0, i, pWalkerArray, cellBounds, pWalkersStartOutStream);
		}

		if(pWalkersOutStream)
		{
			report_walker(0, i, pWalkerArray, cellBounds, pWalkersOutStream);
		}
	}

    close_stream(pWalkersStartOutStream);

    gettimeofday(&time_data_current, NULL);
    fprintf(stderr, "# %u Walkers Placed in %u seconds.\n", nWalkers, (unsigned int)(time_data_current.tv_sec - time_data_old.tv_sec));

    // }} place the walkers

	gettimeofday(&time_data_old, NULL);
    fprintf(stderr, "# Starting simulation.\n");

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

			// write walker information to site
			if(pWalkersOutStream)
			{
				report_walker(j+1, i, pWalkerArray, cellBounds, pWalkersOutStream);
			}

			if(bBouncedOnBarrier)
				nInterferenceCounter++;

			fAverageSquareDistance += pWalkerArray[i].squareDistanceFromOrigin(&cellBounds);
			fAverageX += pWalkerArray[i].m_position.m_fx - pWalkerArray[i].m_origin.m_fx;
			fAverageY += pWalkerArray[i].m_position.m_fy - pWalkerArray[i].m_origin.m_fx;
		}

		fAverageSquareDistance /= (double)nWalkers;
		fAverageX /= (double)nWalkers;
		fAverageY /= (double)nWalkers;

		// write statistical information to file
		if(j > nFirstIterationToPrint)
		{
            fprintf(pStatisticsOutStream, "%d\t%f\t%u\t%u\t%f\t%f\t%u\n",
            j,
            fAverageSquareDistance,
            nInterferenceCounter,
            nInterferenceCounter - nOldInterferenceCounter,
            fOccupiedAreaRatio,
            fTotalBarrierPerimeter,
            nBarriers
            );
		}

		// print some information to stderr every nStepsToShowErr
		if(!(j % nStepsToShowErr))
		{
			gettimeofday(&time_data_current, NULL);
			fprintf(stderr, "%.2f %% Complete in %d seconds.\n", 100.0*((double)j)/((double)nIterations), (unsigned int)(time_data_current.tv_sec - time_data_old.tv_sec));
			fprintf(stderr, "%d\t%f\t%u\t%u\t%f\t%f\t%u\n",
            j,
            fAverageSquareDistance,
            nInterferenceCounter,
            nInterferenceCounter - nOldInterferenceCounter,
            fOccupiedAreaRatio,
            fTotalBarrierPerimeter,
            nBarriers
            );
		}
    }

	gettimeofday(&time_data_current, NULL);
    fprintf(stderr, "# Simulation complete in %u segundos.\n", (unsigned int)(time_data_current.tv_sec - time_data_old.tv_sec));

    if(pWalkersEndOutStream)
	{
		fprintf(stderr, "# Writing final state of walkers.\n");

		for(i = 0; i < nWalkers; i++)
		{
			report_walker(nIterations, i, pWalkerArray, cellBounds, pWalkersEndOutStream);
		}
	}

clean_up:

	close_stream(pStatisticsOutStream);
	close_stream(pWalkersOutStream);
	close_stream(pWalkersEndOutStream);

	// deallocate memory
	delete [] pWalkerArray;
	delete [] pWalkerBarrierArray;

	return 1;
}
