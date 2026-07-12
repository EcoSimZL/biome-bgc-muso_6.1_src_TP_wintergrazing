/* 
potEVAP_to_actEVAP.c
calculation of actual soil evaporation from potential soil evaporation regarding to dry day limit (Joe Rictchie, 1970)

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

int potEVAP_to_actEVAP(control_struct *ctrl, const siteconst_struct *sitec, const soilprop_struct *sprop,
					   epvar_struct *epv, wstate_struct *ws, wflux_struct *wf)
{

	/* internal variables */
	int errorCode = 0;
	double infiltPOT, soilw_hw0, evap_lack;
	infiltPOT = wf->infiltPOT;

	/* 1. first evaporation phase - no limit */
	if (ws->soil_evapCUM1 < sprop->soilevapLIM)
	{
		epv->dsr = 0;

		if (infiltPOT >= ws->soil_evapCUM1)
		{
			ws->soil_evapCUM1 = 0;
		}
		else
		{
			ws->soil_evapCUM1 -= infiltPOT;
		}

		if (!errorCode && evapPHASE1toPHASE2(sprop, epv, ws, wf))
		{
			printf("\n");
			printf("ERROR in evapPHASE1toPHASE2() in POTevap_to_actEVAP.c\n");
			errorCode = 1;
		}
	}
	/* 2. second evaporation phase - limit of dsr: ws->soil_evapCUM1 > soilevapLIM */
	else
	{
		if (infiltPOT >= ws->soil_evapCUM2)
		{
			epv->dsr = 0;

			/* berszivárgó csapadék feltölt CUM2-t, a maradék meg a CUM1-bõl amennyit csak tud */
			infiltPOT -= ws->soil_evapCUM2;
			ws->soil_evapCUM2 = 0;

			if (infiltPOT > sprop->soilevapLIM)
				ws->soil_evapCUM1 = 0;
			else
				ws->soil_evapCUM1 -= infiltPOT;

			if (!errorCode && evapPHASE1toPHASE2(sprop, epv, ws, wf))
			{
				printf("\n");
				printf("ERROR in evapPHASE1toPHASE2() in POTevap_to_actEVAP.c\n");
				errorCode = 1;
			}
		}
		else
		{
			epv->dsr += 1;
			wf->soilw_evap = 3.5 * pow(epv->dsr, 0.5) - ws->soil_evapCUM2;

			if (wf->soilw_evap > wf->evapPOT)
				wf->soilw_evap = wf->evapPOT;

			/* beszivárgásnak köszönhetõen visszamegy */
			ws->soil_evapCUM2 += wf->soilw_evap - infiltPOT;
			epv->dsr = pow((ws->soil_evapCUM2 / 3.5), 2);
		}
	}

	/* control */
	if (wf->soilw_evap < 0 || epv->dsr < 0)
	{
		printf("ERROR in soil evaporation calculation (potEVAP_to_actEVAP)\n");
		errorCode = 1;
	}

	/* 3. Soil evaporation can not be larger than the extractable soil water in the top layer */

	/* actual soil water content at theoretical lower limit of water content: hygroscopic water content */
	soilw_hw0 = sprop->vwc_hw[0] * sitec->soillayer_thickness[0] * water_density;

	/* evap_lack: control parameter to avoid negative soil water content (due to overestimated evaporation + dry soil) */
	evap_lack = wf->soilw_evap - (ws->soilw[0] - soilw_hw0);

	/* theoretical lower limit of water content: hygroscopic water content. */
	if (evap_lack > 0)
	{
		if (ws->soil_evapCUM2 >= sprop->soilevapLIM)
		{
			if (ws->soil_evapCUM2 > wf->soilw_evap)
			{
				ws->soil_evapCUM2 = ws->soil_evapCUM2 - wf->soilw_evap + (ws->soilw[0] - soilw_hw0);
				epv->dsr = pow((ws->soil_evapCUM2 / 3.5), 2);
				wf->soilw_evap = (ws->soilw[0] - soilw_hw0);
			}
			else
			{
				ws->soil_evapCUM1 = ws->soil_evapCUM1 - (wf->soilw_evap - ws->soil_evapCUM2);

				ws->soil_evapCUM2 = (ws->soilw[0] - soilw_hw0) - sprop->soilevapLIM;
				if (ws->soil_evapCUM2 < 0)
					ws->soil_evapCUM2 = 0;

				ws->soil_evapCUM1 = ws->soil_evapCUM1 + (ws->soilw[0] - soilw_hw0);
				if (ws->soil_evapCUM1 < sprop->soilevapLIM)
					ws->soil_evapCUM2 = sprop->soilevapLIM;

				epv->dsr = pow((ws->soil_evapCUM2 / 3.5), 2);

				wf->soilw_evap = ws->soilw[0] - soilw_hw0;
			}
		}
		else
		{
			ws->soil_evapCUM1 = ws->soil_evapCUM1 - wf->soilw_evap + (ws->soilw[0] - soilw_hw0);
			wf->soilw_evap = ws->soilw[0] - soilw_hw0;
		}

		/* limitevap_flag: flag of WARNING writing in log file (only at first time) */
		if (fabs(evap_lack) > CRIT_PREC && !ctrl->limitevap_flag)
			ctrl->limitevap_flag = 1;
	}

	ws->soilw[0] -= wf->soilw_evap;
	epv->vwc[0] = ws->soilw[0] / water_density / sitec->soillayer_thickness[0];

	return (errorCode);
}

int evapPHASE1toPHASE2(const soilprop_struct *sprop, epvar_struct *epv, wstate_struct *ws, wflux_struct *wf)
{
	/* internal variables */
	int errorCode = 0;

	ws->soil_evapCUM1 += wf->evapPOT;
	if (ws->soil_evapCUM1 > sprop->soilevapLIM)
	{
		wf->soilw_evap = wf->evapPOT - 0.4 * (ws->soil_evapCUM1 - sprop->soilevapLIM);
		ws->soil_evapCUM2 = 0.6 * (ws->soil_evapCUM1 - sprop->soilevapLIM);
		epv->dsr = pow(ws->soil_evapCUM2 / 3.5, 2);
		ws->soil_evapCUM1 = sprop->soilevapLIM;
	}
	else
	{
		wf->soilw_evap = wf->evapPOT;
	}

	wf->soilw_evap *= epv->evapREDmulch;

	return (errorCode);
}