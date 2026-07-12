/* 
ploughing.c
do ploughing  - decrease the plant material (leafc, leafn, canopy water)

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
#include "pointbgc_struct.h"
#include "bgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"

int mulch(siteconst_struct *sitec, soilprop_struct *sprop, metvar_struct *metv, epvar_struct *epv, cstate_struct *cs, wstate_struct *ws)
{

	/* mulching parameters */
	int layer_mulcheff, layer, errorCode;
	double soilw_SUM, tsoil_SUM;
	double paramCOV1, paramCOV2, paramCOV3, paramRED1;

	errorCode = 0;
	tsoil_SUM = soilw_SUM = 0;

	layer_mulcheff = 2;
	paramCOV1 = 100;
	paramCOV2 = 0.75;
	paramCOV3 = 0.75;
	paramRED1 = 0.5;

	/**********************************************************************************************/
	/* I. CALCULATING the effect of mulch on the top soil layers: UNIFORM DISTRIBUTION OF VWC and tsoil */

	epv->mulch_coverage = paramCOV1 * pow((cs->litr1c[0] + cs->cwdc[0]) / paramCOV2, paramCOV3);
	epv->evapREDmulch = pow(paramRED1, epv->mulch_coverage / 100);

	if (cs->litr1c[0] > 0)
	{
		for (layer = 0; layer < layer_mulcheff; layer++)
		{
			tsoil_SUM += metv->tsoil[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[layer_mulcheff - 1];
		}

		for (layer = 0; layer < layer_mulcheff; layer++)
		{
			metv->tsoil[layer] = tsoil_SUM;
		}

		metv->tsoil_surface_pre = metv->tsoil[0];
		metv->tsoil_surface = metv->tsoil[0];
	}

	return (errorCode);
}
