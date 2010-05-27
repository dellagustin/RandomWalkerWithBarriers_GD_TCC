// #include <iostream>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef int BOOL;

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

struct main_data
{
	int nWalkers;
    int nIterations;
    int nBarriers;
    int nBarriersMin;
    int nBarriersMax;
    int nScripts;
    double fBarrierRadius;
    double fBarrierRadiusMin;
    double fBarrierRadiusMax;
    double fDesiredOccupiedArea;
    BOOL bByArea;
    char *pszBarrierOutFileName;
};

void parse_cmd_line(int argc, char* argv[], main_data& data)
{
	int i;
	for(i = 0; i < argc; i++)
	{
		if(argv[i][0]!='-')
			continue;

		if(!strcmp(argv[i]+1, "ua")) // use fixed area
		{
			data.bByArea = TRUE;
		}
		else if(!strcmp(argv[i]+1, "da")) // desired area
		{
			if(i+1 >= argc) continue;
			data.fDesiredOccupiedArea = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i]+1, "ns")) // desired area
		{
			if(i+1 >= argc) continue;
			data.nScripts = atoi(argv[i+1]);
		}
	}
}

int main(int argc, char* argv[])
{
    struct main_data data;
    char **ppszFileNames;
    memset(&data, 0, sizeof(data));

    data.nWalkers = 8000;
    data.nIterations = 8000;
	data.bByArea = TRUE;
    data.fDesiredOccupiedArea = 5000;
    data.nScripts = 20;
    data.nBarriersMax = 1600;
    data.nBarriersMin = 1;
	data.fBarrierRadiusMax = 50.0;
    data.fBarrierRadiusMax = 1.0;

    parse_cmd_line(argc, argv, data);

    ppszFileNames= new char*[data.nScripts];
    memset(ppszFileNames, 0, data.nScripts*sizeof(char*));

	if(data.bByArea)
	{
		// needs some error verification here in the limit values
		int i;
		for(i = 0; i < data.nScripts; i++)
		{
			do
			{
				data.nBarriers = data.nBarriersMin + rand()%(data.nBarriersMax-data.nBarriersMin);
				data.fBarrierRadius = sqrt((data.fDesiredOccupiedArea/((double)data.nBarriers))/M_PI);
			}
			while(data.fBarrierRadius >= data.fBarrierRadiusMin && data.fBarrierRadius <= data.fBarrierRadiusMax);

			char szFileName[256];
			sprintf(szFileName, "w=%d_i=%d_b=%d_br=%.5f.qss", data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius);
			FILE *pOutFile = fopen(szFileName, "wt");

			if(!pOutFile)
				continue;

			fprintf(pOutFile, "#!/bin/bash\n\n");
			fprintf(pOutFile, "#PBS -j oe\n");
			fprintf(pOutFile, "#PBS -l ncpus=1\n");
			fprintf(pOutFile, "#PBS -q q_um_dia\n");
			fprintf(pOutFile, "#PBS -N rwwb_w%di%db%dbr%.5f\n\n", data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius);
			fprintf(pOutFile, "work_dir=$PBS_O_WORKDIR\n");
			fprintf(pOutFile, "echo $PBS_O_WORKDIR\n");
			fprintf(pOutFile, "cd $work_dir\n\n");
			fprintf(pOutFile, "./rwwbcompile\n");
			fprintf(pOutFile, "./randomwwb -w %d -i %d -b %d -br %f > test_file.dat\n", data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius);

			fclose(pOutFile);
		}
	}

    return 0;
}
