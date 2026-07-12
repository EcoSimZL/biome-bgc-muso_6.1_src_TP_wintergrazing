/* 
soilstress_calculation.c
calculation of drought and saturation soil water content stress layer by layers 

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

// 消除两个if -ZL
// #define soilstress_debug 1

int soilstress_calculation(const siteconst_struct *sitec, const soilprop_struct *sprop, const epconst_struct *epc,
						   epvar_struct *epv, wstate_struct *ws, wflux_struct *wf)
{
	int layer;
	double m_vwcR_layer, m_SWCstress_avg, m_Nstress_avg;

	int errorCode = 0;
	m_SWCstress_avg = m_Nstress_avg = 0;

	/* 1. Calculation SWC-stress (can occur only after emergence  - if root is presence) */
	if (epv->rootdepth)
	{
		/*--------------------------------------------*/
		/* 1.1 SWC-stress calculation based on VWC   */
		if (epc->soilstress_flag == 0)
		{
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				m_vwcR_layer = 1;

				if (epv->vwc[layer] > sprop->vwc_wp[layer] && (sprop->vwc_sat[layer] - epv->vwc[layer]) > CRIT_PREC)
				{
					/* DROUGHT STRESS */
					if (epv->vwc[layer] <= epv->vwc_crit1[layer])
					{
						m_vwcR_layer = pow((epv->vwc[layer] - sprop->vwc_wp[layer]) / (epv->vwc_crit1[layer] - sprop->vwc_wp[layer]), sprop->q_soilstress);
					}

					/* ANOXIC CONDITION STRESS */
					if (epv->vwc[layer] >= epv->vwc_crit2[layer])
					{
						m_vwcR_layer = (sprop->vwc_sat[layer] - epv->vwc[layer]) / (sprop->vwc_sat[layer] - epv->vwc_crit2[layer]);
					}
				}
				else
					m_vwcR_layer = 0;

				/* lower limit for saturation: m_fullstress2 */
				if (epv->vwc[layer] >= epv->vwc_crit1[layer] && m_vwcR_layer < epc->m_fullstress2)
					m_vwcR_layer = epc->m_fullstress2;

#ifdef soilstress_debug
				/* if VWCratio_crit is set to -9999: no soilstress calculation -> m_SWCstress = 1 */
				if (epc->VWCratio_crit1 == DATA_GAP && epv->vwc[layer] <= sprop->vwc_fc[layer])
					m_vwcR_layer = 1;
				if (epc->VWCratio_crit2 == DATA_GAP && epv->vwc[layer] >= sprop->vwc_fc[layer])
					m_vwcR_layer = 1;
#endif

				epv->m_SWCstress_layer[layer] = m_vwcR_layer;
				m_SWCstress_avg += epv->m_SWCstress_layer[layer] * epv->rootlength_prop[layer];
			}

			epv->m_SWCstress = m_SWCstress_avg;
		}
		/* 1.2 SWC-stress calculation based on transpiration demand-possibitiy  */
		else
		{
			if (epv->n_rootlayers)
			{
				for (layer = 0; layer < epv->n_rootlayers; layer++)
				{
					if (wf->soilw_transDEMAND[layer])
						epv->m_SWCstress_layer[layer] = wf->soilw_trans[layer] / wf->soilw_transDEMAND[layer];
					else
						epv->m_SWCstress_layer[layer] = 0;

					/* control */
					if (epv->m_SWCstress_layer[layer] < 0 || epv->m_SWCstress_layer[layer] > 1)
					{
						printf("\n");
						printf("ERROR: soilstress calculation error in multilayer_hydrolpocess.c:\n");
						errorCode = 1;
					}

					m_SWCstress_avg += epv->m_SWCstress_layer[layer] * epv->rootlength_prop[layer];
				}
			}
			else
			{
				for (layer = 0; layer < N_SOILLAYERS; layer++)
					epv->m_SWCstress_layer[layer] = 1;
				m_SWCstress_avg = 1;
			}
		}
	}
	else
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			epv->m_SWCstress_layer[layer] = 1;
			ws->soilw_avail[layer] = 0;
		}

		m_SWCstress_avg = 1;
	}

	epv->m_SWCstress = m_SWCstress_avg;

	/****************************************************************************************/
	/* 2. calculating cumulative SWC stress and extreme temperature effect */

	/* 2.1 calculating WATER STRESS DAYS regarding to the average soil moisture conditions */

	// 因为如果有任何一层亏缺，那么就会开始计算累积，这样的话，表层是非常容易亏缺，很容易造成，水分亏缺，因此总和不以等于1作为判断，而是有大于0.85-ZL-20211011
	// 如果m_SWCstress没有压力，那么length为0
	// 原文:if (epv->m_SWCstress == 1) epv->SWCstressLENGTH = 0;
	// 因为现在根据参数(generic_EPC_set_6.1)，将水分亏缺定义为原来的一半，因此将这里的>=0.85改回原来的==1-ZL-211012
	// 还是修改为0.85比较好，因为如果是1，需要所有各层同时满足，这个对于多层来说是比较难-ZL-211012
	if (epv->m_SWCstress >= 0.85)
		epv->SWCstressLENGTH = 0;
	if (epv->n_actphen > epc->n_emerg_phenophase)
	{
		epv->SWCstressLENGTH += (1 - epv->m_SWCstress);
		epv->cumSWCstress += (1 - epv->m_SWCstress);
	}
	else
	{
		epv->SWCstressLENGTH = 0;
		epv->cumSWCstress = 0;
	}

	// 在手册里写的一个天数的限制，即crit是一个天数，但从这里来看，并不是如此
	// 但这样看起来更加合理，因此length不仅仅是长度的单位，更是强度的积累
	if (epv->SWCstressLENGTH < epc->SWCstressLENGTH_crit)
		epv->m_SWCstressLENGTH = 1;
	else
		epv->m_SWCstressLENGTH = epc->SWCstressLENGTH_crit / epv->cumSWCstress;

	/****************************************************************************************/
	/* 3. N-stress based on immobilization ratio */

	if (epv->rootdepth)
		for (layer = 0; layer < N_SOILLAYERS; layer++)
			m_Nstress_avg += epv->IMMOBratio[layer] * epv->rootlength_prop[layer];
	else
		m_Nstress_avg = 1;

	if (1 - m_Nstress_avg < CRIT_PREC)
		m_Nstress_avg = 1;

	if (epv->n_actphen > epc->n_emerg_phenophase)
		epv->cumNstress += (1 - m_Nstress_avg);
	else
		epv->cumNstress = 0;

	return (errorCode);
}