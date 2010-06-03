#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rwwb_common.h"

BOOL get_field(unsigned int uiField, const char *pcszInStr, char *pszField, int nMaxFieldSize)
{
	/*
	int i, j, nLength;
	nLength = strlen(pcszInsStr);

	for(i = 0, j = 0; i < nLenght; i++)
	{
		if(pcszInStr[i] == '\t')
		{
			j++;
		}

		if(j > nFields)
			break;
	}
	*/

	unsigned int i;
	char *pszStrText = (char*)pcszInStr;
	const char *pcszFieldSeparator = "\t";

	for(i = 0; i < uiField; i++)
	{
		pszStrText = strstr(pszStrText, pcszFieldSeparator);

		if(!pszStrText)
			return FALSE;

		pszStrText++;
	}

	char *pszStrEndText = strstr(pszStrText, pcszFieldSeparator);

	if(!pszStrEndText)
	{
		strcpy(pszField, pszStrText);
	}
	else
	{
		int nCharsToCopy = pszStrEndText - pszStrText;
		strncpy(pszField, pszStrText, nCharsToCopy);
		pszField[nCharsToCopy] = 0;
	}

	return TRUE;
}

int main(int argc, char *argv[])
{
	char *szInFile;
	// BOOL bUseStdin;
	// BOOL bNormalize;
	FILE* pInStream;
	unsigned int *puiCount;
	unsigned int uiField;
	unsigned int uiIntervals;
	unsigned int i;
	double fMaxValue;
	double fLinearStep;

	szInFile = NULL;
	pInStream = stdin;
	uiField = 6;
	uiIntervals = 100;
	fMaxValue = 0.0;

	// read the command line arguments
	for(i = 0; i < (unsigned int)argc; i++)
	{
		if(argv[i][0] != '-')
			continue;

		if(!strcmp(argv[i]+1, "if")) // number of iterations
		{
			if(i+1<(unsigned int)argc)
			szInFile = strdup(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "i")) // number of iterations
		{
			if(i+1<(unsigned int)argc)
			uiIntervals = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "f")) // number of iterations
		{
			if(i+1<(unsigned int)argc)
			uiField = atoi(argv[i+1]);
		}
	}

	// temporary
	// szInFile = strdup("walker_data_final.dat");

	// memory allocation
	puiCount = new unsigned int[uiIntervals];
	memset(puiCount, 0, uiIntervals*sizeof(unsigned int));

	if(szInFile)
	{
		pInStream = fopen(szInFile, "rt");

		if(!pInStream)
		{
			fprintf(stderr, "# can not open file: %s", szInFile);
			goto clean_up;
		}
	}

	char szLineRead[1024];
	while(fgets(szLineRead, sizeof(szLineRead)-1, pInStream))
	{
		if(szLineRead[0] == '#')
			continue;

		char szField[128];

		if(!get_field(uiField, szLineRead, szField, sizeof(szField)-1))
			continue;

		double fValue = atof(szField);

		if(fValue > fMaxValue)
			fMaxValue = fValue;
	}

	fseek(pInStream, 0, SEEK_SET);

	while(fgets(szLineRead, sizeof(szLineRead)-1, pInStream))
	{
		if(szLineRead[0] == '#')
			continue;

		char szField[128];

		if(!get_field(uiField, szLineRead, szField, sizeof(szField)-1))
			continue;

		double fValue = atof(szField);

		unsigned int uiIndex = ((double)uiIntervals)*fValue/fMaxValue;
		uiIndex = __min(uiIndex, (uiIntervals-1));
		puiCount[uiIndex]++;
	}

	fLinearStep = fMaxValue*((double)uiIntervals);

	for(i = 0; i < uiIntervals; i++)
	{
		fprintf(stdout, "%f\t%f\t%f\t%d\n", ((double)i)*fLinearStep, (((double)i)+0.5)*fLinearStep, (((double)i)+1.0)*fLinearStep, puiCount[i]);
	}

	clean_up:

	if(pInStream && pInStream != stdin)
	{
		fclose(pInStream);
	}

	if(puiCount)
	{
		delete puiCount;
	}

	if(szInFile)
		delete szInFile;
}
