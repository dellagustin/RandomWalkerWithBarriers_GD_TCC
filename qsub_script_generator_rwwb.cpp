#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "rwwb_common.h"

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
    // char **ppszFileNames;
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

    if(data.bByArea)
	{
		// ppszFileNames= new char*[data.nScripts];
		// memset(ppszFileNames, 0, data.nScripts*sizeof(char*));
		FILE *pQSubCallerOutStream = fopen("qsub_caller", "wt");

		if(!pQSubCallerOutStream)
		{
			fprintf(stderr, "Impossible to create qsub_caller file.\n");
			return 0;
		}

		fprintf(pQSubCallerOutStream, "#!/bin/bash\n");

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

			int j;
			for(j = 0; j < 2; j++)
			{
				char szFileName[256];
				sprintf(szFileName, "s_%d_w_%d_i_%d_b_%d_br_%.5f_m%d.qss", i, data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius, j+1);
				fprintf(pQSubCallerOutStream, "qsub %s\n", szFileName);
				fprintf(pQSubCallerOutStream, "sleep 1\n");
				FILE *pOutFile = fopen(szFileName, "wt");

				if(!pOutFile)
					continue;

				fprintf(pOutFile, "#!/bin/bash\n\n");
				fprintf(pOutFile, "#PBS -j oe\n");
				fprintf(pOutFile, "#PBS -l ncpus=1\n");
				fprintf(pOutFile, "#PBS -q q_um_dia\n");
				fprintf(pOutFile, "#PBS -N rwwb_s%d_w%d_i%d_b%d_br%.5f_m%d\n\n", i, data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius, j+1);
				fprintf(pOutFile, "work_dir=$PBS_O_WORKDIR\n");
				fprintf(pOutFile, "echo $PBS_O_WORKDIR\n");
				fprintf(pOutFile, "cd $work_dir\n\n");
				fprintf(pOutFile, "./rwwbcompile\n");
				fprintf(pOutFile, "./randomwwb -w %d -i %d -b %d -br %f -m %d", data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius, j+1);
				fprintf(pOutFile, " > res_s%d_w%d_i%d_b%d_br%.5f_m%d.dat", i, data.nWalkers, data.nIterations, data.nBarriers, data.fBarrierRadius, j+1);

				fclose(pOutFile);
			}
		}

		fclose(pQSubCallerOutStream);
	}

    return 0;
}
