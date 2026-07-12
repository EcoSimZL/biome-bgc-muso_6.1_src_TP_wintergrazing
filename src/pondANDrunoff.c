/* 
richards.c
calculation of soil water content layer by layer taking into account soil hydrological processes 
(percolation, diffusion)  using inamic discretization level in order to taking account the change of the diffusivity and conductivity depending on soil moisture 
(after rain event the soil moisture rises suddenly causing high diffusivity and conductivity. If we assumed that these diffusivity 
 and conductivity values were valid all the day we would overestimate the velocity diffusion and percolation process) 

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v6.1.
Copyright 2020, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_constants.h"
#include "bgc_func.h"

int pondANDrunoff(const siteconst_struct *sitec, const soilprop_struct *sprop, epvar_struct *epv, wstate_struct *ws, wflux_struct *wf)
{

	/* internal variables */
	int errorCode = 0;
	double RCN, coeff_soiltype, coeff_soilmoist, coeff_runoff, INFILT, soilw_SAT;

	/* sum of potential infitration */
	INFILT = (wf->prcp_to_soilw + wf->snoww_to_soilw + wf->canopyw_to_soilw + wf->IRG_to_prcp + ws->pond_water);

	/* --------------------------------------------------------------------------------------------------------*/
	/* 1. POND WATER */

	/* if there is pond water in the area */
	if (ws->pond_water)
	{
		soilw_SAT = sprop->vwc_sat[0] * sitec->soillayer_thickness[0] * water_density;

		/* if potential infiltration plus actual soil water content is greater than saturation -> pond water */
		if (ws->soilw[0] + INFILT > soilw_SAT)
		{
			wf->pondw_to_soilw = soilw_SAT - ws->soilw[0];
			ws->soilw[0] = soilw_SAT;
			if (INFILT - wf->pondw_to_soilw > sprop->pondmax)
			{
				wf->prcp_to_runoff = INFILT - wf->pondw_to_soilw - sprop->pondmax;
				ws->pond_water = sprop->pondmax;
			}
			else
			{
				ws->pond_water = INFILT - wf->pondw_to_soilw;
				wf->prcp_to_runoff = 0;
			}
			wf->infiltPOT = wf->pondw_to_soilw;
		}
		else
		{
			wf->infiltPOT = INFILT;
			wf->pondw_to_soilw = ws->pond_water;
			ws->pond_water = 0;
		}
	}
	else
	{
		RCN = sprop->RCN;
		coeff_soiltype = 254 * (100 / RCN - 1);

		coeff_soilmoist = 0.15 * ((sprop->vwc_sat[0] - epv->vwc[0]) / (sprop->vwc_sat[0] - sprop->vwc_hw[0]));

		coeff_runoff = coeff_soiltype * coeff_soilmoist;

		if (INFILT > coeff_runoff)
		{
			wf->prcp_to_runoff = pow(INFILT - coeff_runoff, 2) / (INFILT + (1 - coeff_soilmoist) * coeff_soiltype);
		}
		else
		{
			wf->prcp_to_runoff = 0;
		}

		wf->infiltPOT = (wf->prcp_to_soilw + wf->snoww_to_soilw + wf->canopyw_to_soilw + wf->IRG_to_prcp + wf->pondw_to_soilw - wf->prcp_to_runoff);
	}

	return (errorCode);
}
