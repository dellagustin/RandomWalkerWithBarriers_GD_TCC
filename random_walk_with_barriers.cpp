#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef BOOL
typedef unsigned int BOOL;
#endif

// returns a random number between 0.0 and fMaxValue 
double randomNumber(double fMaxValue);

// returns a random angle between 0.0 and 2*pi
double randomAngle();

// returns the square of fValue 
double square(double fValue);

double randomNumber(double fMaxValue)
{
	return fMaxValue*((double)rand()/((double)RAND_MAX)); 
}

double randomAngle()
{
    return randomNumber(2.0*M_PI);
}

double square(double fValue)
{
	return fValue * fValue;
}

// 
class CWalkerCellBounds
{
	// constructor
public:
	CWalkerCellBounds();
	
	// methods
	double height() const;
	double width() const;

	// data
public:
	double m_fxmin;
	double m_fxmax;
	double m_fymin;
	double m_fymax;
};

double CWalkerCellBounds::height() const
{
	return m_fymax - m_fymin;
}

double CWalkerCellBounds::width() const
{
	return m_fxmax - m_fxmin;
}

// constructor
CWalkerCellBounds::CWalkerCellBounds()
{
	m_fxmin = 0.0;
	m_fxmax = 0.0;
	m_fymin = 0.0;
	m_fymax = 0.0;
}

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

/*
BOOL iterateWithMovingPoint(const CWalkerPoint& pointPosIni, CWalkerPoint& pointPosEnd) const
{
	if(pointPosEnd.squareDistanceTo(m_position) < square(m_radius))
	{
		return FALSE;
	}
	
	return TRUE;
}
*/

BOOL CWalkerCircularBarrier::isInside(const CWalkerPoint& point) const
{
	return point.squareDistanceTo(m_position) < square(m_radius);
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

int main()
{
	CWalker *pWalkerArray;
	CWalkerCircularBarrier *pWalkerBarrierArray;
	CWalkerCellBounds cellBounds;
    int nWalkers, nIterations, nBarriers, i, j, k;
    double fBarrierRadius;
	FILE *pOutStream = stdout;
	
	nBarriers = 200;
    nWalkers = 2000;
    nIterations = 16000;
	cellBounds.m_fxmax = 100.0;
	cellBounds.m_fymax= 100.0;
	fBarrierRadius = 1.0;
	
	// print some information about the simulation
	fprintf(pOutStream, "# data generated by random walker with barriers\n"); 
	fprintf(pOutStream, "# by guilherme dellagustin\n");
	fprintf(pOutStream, "# walkers: %d\n", nWalkers);
	fprintf(pOutStream, "# iterations: %d\n", nIterations); 
	fprintf(pOutStream, "# min bounds: %f, %f\n", cellBounds.m_fxmin, cellBounds.m_fymin); 
	fprintf(pOutStream, "# max bounds: %f, %f\n", cellBounds.m_fxmax, cellBounds.m_fymax); 
	    
    // allocate memory
	pWalkerArray = new CWalker[nWalkers];
	pWalkerBarrierArray = new CWalkerCircularBarrier[nBarriers];
    	
	// initiate data
	for(i = 0; i < nWalkers; i++)
	{
		pWalkerArray[i].m_origin.m_fx = cellBounds.m_fxmin + 0.5 * cellBounds.width();
		pWalkerArray[i].m_origin.m_fy = cellBounds.m_fymin + 0.5 * cellBounds.height();
		pWalkerArray[i].m_position.m_fx = pWalkerArray[i].m_origin.m_fx;
		pWalkerArray[i].m_position.m_fy = pWalkerArray[i].m_origin.m_fy;
	}

	for(i = 0; i < nBarriers; i++)
	{
		pWalkerBarrierArray[i].m_position.m_fx = cellBounds.m_fxmin + randomNumber(cellBounds.width());
		pWalkerBarrierArray[i].m_position.m_fy = cellBounds.m_fymin + randomNumber(cellBounds.height());
		pWalkerBarrierArray[i].m_radius = 1.0;
	}
	
    // iteration in 'time'
	for(j = 0; j < nIterations; j++)
    {
		double fAverageSquareDistance = 0.0;
		double fAverageX = 0.0;
		double fAverageY = 0.0;
		
		// iteration in walkers
		for(i = 0 ; i < nWalkers; i++)
		{
			CWalkerPoint oldPoint = pWalkerArray[i].m_position;
			pWalkerArray[i].iterate(1, &cellBounds);
			
			BOOL bOutOfBarriers;
			
			do
			{
				bOutOfBarriers = TRUE;
				
				for(k = 0; k < nBarriers; k++)
				{
					// if(pWalkerBarrier
				}
			}
			while(bOutOfBarriers)
			
			// 
			fAverageSquareDistance += pWalkerArray[i].squareDistanceFromOrigin(&cellBounds);
			fAverageX += pWalkerArray[i].m_position.m_fx;
			fAverageY += pWalkerArray[i].m_position.m_fy;
		}
		
		fAverageSquareDistance /= (double)nWalkers;
		fAverageX /= (double)nWalkers;
		fAverageY /= (double)nWalkers;
		
		fprintf(pOutStream, "%d\t%f\t%f\t%f\n", j, fAverageSquareDistance, fAverageX, fAverageY);
    }    
    
	// deallocate memory
	delete [] pWalkerArray ;
}