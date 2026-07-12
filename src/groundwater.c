/* 
groundwater.c
calculate the effect of gound water depth in the soil

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

int groundwater(const siteconst_struct *sitec, const soilprop_struct *sprop, control_struct *ctrl, epvar_struct *epv, wstate_struct *ws, wflux_struct *wf)
{
	int groundwater = 0;
	int layerSAT;
	int errorCode = 0;
	int layer = 0;
	double weight;
	double soilw_sat, soilw_befGW, soilw_aftGW;

	/* if gound water depth information is available using in multilayer calculation */
	if (ctrl->GWD_flag)
	{
		/* in spinup and transient run only 365 values are used (daily averages) */
		if (ctrl->spinup == 0)
			ctrl->gwd_act = sitec->gwd_array[(ctrl->simyr * nDAYS_OF_YEAR) + ctrl->yday];
		else
			ctrl->gwd_act = sitec->gwd_array[ctrl->yday];

		/* calculate processes layer to layer (i0: actual layer, i1:deeper layer)  - EXCEPT OF THE BOTTOM LAYER */
		while (groundwater == 0 && layer < N_SOILLAYERS)
		{
			/* actual groundwater level is above the lower boundary of bottom layer (3m) */
			if (ctrl->gwd_act < sitec->soillayer_depth[layer])
			{

				soilw_sat = sprop->vwc_sat[layer] * sitec->soillayer_thickness[layer] * water_density;
				soilw_befGW = ws->soilw[layer];

				weight = (sitec->soillayer_depth[layer] - ctrl->gwd_act) / sitec->soillayer_thickness[layer];

				soilw_aftGW = soilw_sat * weight + soilw_befGW * (1 - weight);
				wf->soilw_from_GW[layer] = soilw_aftGW - soilw_befGW;

				ws->soilw[layer] = soilw_aftGW;
				epv->vwc[layer] = ws->soilw[layer] / water_density / sitec->soillayer_thickness[layer];

				groundwater = 1;

				/* soil layers below the groundwater level are saturated - net water gain from soil system */
				for (layerSAT = 1; layer + layerSAT < N_SOILLAYERS; layerSAT++)
				{
					soilw_sat = sprop->vwc_sat[layer + layerSAT] * sitec->soillayer_thickness[layer + layerSAT] * water_density;
					wf->soilw_from_GW[layer + layerSAT] = soilw_sat - ws->soilw[layer + layerSAT];
					ws->soilw[layer + layerSAT] = soilw_sat;
					epv->vwc[layer + layerSAT] = sprop->vwc_sat[layer + layerSAT];
				}
			}
			layer += 1;
		}
	}
	else
		ctrl->gwd_act = DATA_GAP;

	return (errorCode);
}
