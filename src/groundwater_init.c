/* 
groundwater_init.c
read ground water depth information if it is available

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v6.1.
Copyright 2020, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"

int groundwater_init(siteconst_struct *sitec, control_struct *ctrl)
{
	int errorCode = 0;
	int i;
	int n_data = 0;
	file gwd_file;

	/* ------------------------------------------------------ */
	/* using varying whole plant mortality values */

	/* gwd flag: constans or varying gwd from file */
	ctrl->GWD_flag = 0;
	if (!errorCode)
	{
		if (ctrl->spinup == 0) /* normal run */
		{
			strcpy(gwd_file.name, "groundwater_normal.txt");
			if (!file_open(&gwd_file, 'j', 1))
			{
				ctrl->GWD_flag = 1;
			}
			n_data = ctrl->simyears * nDAYS_OF_YEAR;
		}
		else /* spinup and transient run */
		{

			strcpy(gwd_file.name, "groundwater_spinup.txt");

			if (!file_open(&gwd_file, 'j', 1))
			{
				ctrl->GWD_flag = 1;
			}
			n_data = nDAYS_OF_YEAR;
		}
	}

	if (!errorCode && ctrl->GWD_flag)
	{
		/* allocate space for the annual gwd array */
		sitec->gwd_array = (double *)malloc(ctrl->simyears * nDAYS_OF_YEAR * sizeof(double));
#ifdef ERROR_PRINT
		if (!sitec->gwd_array)
		{
			printf("ERROR allocating for gwd array, groundwater_init()\n");
			errorCode = 219;
		}
#endif
		/* read year and gwd for each simday in each simyear */
		for (i = 0; !errorCode && i < n_data; i++)
		{
			sitec->gwd_array[i] = DATA_GAP;
			if (fscanf(gwd_file.ptr, "%*i%lf", &(sitec->gwd_array[i])) == EOF)
			{
				printf("ERROR reading annual gwd array, groundwater_init()\n");
				printf("Note: file must contain a pair of values for each simday\n");
				printf("simyear: year and gwd.\n");
				errorCode = 21901;
			}
			if (sitec->gwd_array[i] < 0.0)
			{
				printf("ERROR in groundwater_init(): negative groundwater data OR too short GW-file\n");
				errorCode = 21902;
			}
		}
		fclose(gwd_file.ptr);
	}
	else /* if no changing data constant sitec parameter are used */
	{
		sitec->gwd_array = 0;
	}

	/* -------------------------------------------*/

	return (errorCode);
}
