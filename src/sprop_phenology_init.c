/* 
sprop_phenology_init.c
read sprop file for pointbgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v6.1.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2020, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
20.03.01 Galina Churkina added variable "sum" substituting  t1+t2+t3 in IF statement,
which gave an error.
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
#include "bgc_func.h"

int sprop_phenology_init(file init, soilprop_struct *sprop, control_struct *ctrl, epconst_struct *epc)
{
	int errorCode = 0;
	int layer, scanflag;
	int dofilecloseSOILPROP = 1;
	file sprop_file;
	char key[] = "SOIL_PHENOLOGY_FILE";
	char keyword[STRINGSIZE];
	char header[STRINGSIZE];

	/* original soil C:N values */
	double soilCN12_orig = 12;
	double soilCN34_orig = 10;

	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** SOIL_FILE                                                       ** 
	**                                                                 **
	********************************************************************/

	/* scan for the SOIL file keyword, exit if not next */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword for control data\n");
		errorCode = 208;
	}
	if (!errorCode && strcmp(keyword, key))
	{
		printf("Expecting keyword --> %s in file %s\n", key, init.name);
		errorCode = 208;
	}
	/* open simple SOIL file  */
	if (!errorCode && scan_open(init, &sprop_file, 'r', 1))
	{
		printf("ERROR opening sproponst file, sprop_phenology_init()\n");
		dofilecloseSOILPROP = 0;
		errorCode = 20800;
	}

	/* first scan sprop keyword to ensure proper *.sprop format */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading header, sprop_phenology_init()\n");
		errorCode = 20800;
	}

	/****************************************************************************************************************/
	/* dividing line from file */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 1. dividing line, sprop_phenology_init()\n");
		errorCode = 20801;
	}
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 1. dividing line, sprop_phenology_init()\n");
		errorCode = 20801;
	}
	/****************************************************************************************************************/
	/* NITROGEN AND DECOMPOSITION PARAMETERS  */
	/****************************************************************************************************************/
	/*  denetirification proportion and mobilen proportion
		original: BBGC constant - new version: can be set in SOIL file*/
	if (!errorCode && scan_value(sprop_file, &sprop->denitr_coeff, 'd'))
	{
		printf("ERROR reading denitr_coeff: sprop_phenology_init()\n");
		errorCode = 20802;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->netMiner_to_nitrif, 'd'))
	{
		printf("ERROR reading netMiner_to_nitrif: sprop_phenology_init()\n");
		errorCode = 20803;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->maxNitrif_rate, 'd'))
	{
		printf("ERROR reading maxNitrif_rate: sprop_phenology_init()\n");
		errorCode = 20803;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->N2Ocoeff_nitrif, 'd'))
	{
		printf("ERROR reading N2Ocoeff_nitrif: sprop_phenology_init()\n");
		errorCode = 20804;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->NH4_mobilen_prop, 'd'))
	{
		printf("ERROR reading NH4_mobilen_prop: sprop_phenology_init()\n");
		errorCode = 20805;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->texDEP_N2Oratio, 'd'))
	{
		printf("ERROR reading texDEP_N2Oratio: sprop_phenology_init()\n");
		errorCode = 20806;
	}

	sprop->critWFPS_denitr = 0.55;

	if (!errorCode && scan_value(sprop_file, &sprop->efolding_depth, 'd'))
	{
		printf("ERROR reading e-folding depth: sprop_phenology_init()\n");
		errorCode = 20807;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->SOIL1_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL1_dissolv_prop: sprop_phenology_init()\n");
		errorCode = 20808;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->SOIL2_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL2_dissolv_prop: sprop_phenology_init()\n");
		errorCode = 20808;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->SOIL3_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL3_dissolv_prop: sprop_phenology_init()\n");
		errorCode = 20808;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->SOIL4_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL4_dissolv_prop: sprop_phenology_init()\n");
		errorCode = 20808;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->minWFPS, 'd'))
	{
		printf("ERROR reading WFPS_min: sprop_phenology_init()\n");
		errorCode = 20809;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->opt1WFPS, 'd'))
	{
		printf("ERROR reading WFPS_opt1: sprop_phenology_init()\n");
		errorCode = 20809;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->opt2WFPS, 'd'))
	{
		printf("ERROR reading WFPS_opt2: sprop_phenology_init()\n");
		errorCode = 20809;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->WFPSscalar_min, 'd'))
	{
		printf("ERROR reading WFPSscalar_min: sprop_phenology_init()\n");
		errorCode = 20809;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->soil4_CN, 'd'))
	{
		printf("ERROR reading SOIL4_CN: sprop_phenology_init()\n");
		errorCode = 20809;
	}
	sprop->soil1_CN = sprop->soil4_CN * (soilCN12_orig / soilCN34_orig);
	sprop->soil2_CN = sprop->soil4_CN * (soilCN12_orig / soilCN34_orig);
	sprop->soil3_CN = sprop->soil4_CN;

	/*********************************************************************************************************************/
	/* dividing line from file */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 2. dividing line, sprop_phenology_init()\n");
		errorCode = 20810;
	}
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 2. dividing line, sprop_phenology_init()\n");
		errorCode = 20810;
	}
	/*********************************************************************************************************************/
	/* RATE SCALARS */
	/*********************************************************************************************************************/

	/* respiration fractions for fluxes between compartments  */
	if (!errorCode && scan_value(sprop_file, &sprop->rfl1s1, 'd'))
	{
		printf("ERROR reading rfl1s1: sprop_phenology_init()\n");
		errorCode = 20811;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->rfl2s2, 'd'))
	{
		printf("ERROR reading rfl2s2: sprop_phenology_init()\n");
		errorCode = 20811;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->rfl4s3, 'd'))
	{
		printf("ERROR reading rfl4s3: sprop_phenology_init()\n");
		errorCode = 20811;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->rfs1s2, 'd'))
	{
		printf("ERROR reading rfs1s2: sprop_phenology_init()\n");
		errorCode = 20811;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->rfs2s3, 'd'))
	{
		printf("ERROR reading rfs2s3: sprop_phenology_init()\n");
		errorCode = 20811;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->rfs3s4, 'd'))
	{
		printf("ERROR reading rfs3s4: sprop_phenology_init()\n");
		errorCode = 20811;
	}

	/* 	base values of rate constants are (1/day)   */
	if (!errorCode && scan_value(sprop_file, &sprop->kl1_base, 'd'))
	{
		printf("ERROR reading kl1_base: sprop_phenology_init()\n");
		errorCode = 20812;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->kl2_base, 'd'))
	{
		printf("ERROR reading kl2_base: sprop_phenology_init()\n");
		errorCode = 20812;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->kl4_base, 'd'))
	{
		printf("ERROR reading kl4_base: sprop_phenology_initology_init()\n");
		errorCode = 20812;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->ks1_base, 'd'))
	{
		printf("ERROR reading ks1_base: sprop_phenology_initology_initology_initology_init()\n");
		errorCode = 20812;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->ks2_base, 'd'))
	{
		printf("ERROR reading ks2_base: sprop_phenology_init()\n");
		errorCode = 20812;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->ks3_base, 'd'))
	{
		printf("ERROR reading ks3_base: sprop_phenology_init()\n");
		errorCode = 20811;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->ks4_base, 'd'))
	{
		printf("ERROR reading ks4_base: sprop_phenology_init()\n");
		errorCode = 20812;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->kfrag_base, 'd'))
	{
		printf("ERROR reading kfrag_base: sprop_phenology_init()\n");
		errorCode = 20812;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 3. dividing line, sprop_phenology_init()\n");
		errorCode = 20813;
	}
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 3. dividing line, sprop_phenology_initology_init()\n");
		errorCode = 20813;
	}
	/*********************************************************************************************************************/
	/* CH4 PARAMETERS */
	/*********************************************************************************************************************/

	if (!errorCode && scan_value(sprop_file, &sprop->C_pBD1, 'd'))
	{
		printf("ERROR reading C_pBD1: sprop_phenology_initology_initology_init()\n");
		errorCode = 20814;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->C_pBD2, 'd'))
	{
		printf("ERROR reading C_pBD2: sprop_phenology_initology_initology_init()\n");
		errorCode = 20814;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->C_pVWC1, 'd'))
	{
		printf("ERROR reading C_pVWC1: sprop_phenology_init()\n");
		errorCode = 20814;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->C_pVWC2, 'd'))
	{
		printf("ERROR reading C_pVWC2: sprop_phenology_init()\n");
		errorCode = 20814;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->C_pVWC3, 'd'))
	{
		printf("ERROR reading C_pVWC3: sprop_phenology_initology_init()\n");
		errorCode = 20814;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->C_pVWC4, 'd'))
	{
		printf("ERROR reading C_pVWC4: sprop_phenology_init()\n");
		errorCode = 20814;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->C_pTS, 'd'))
	{
		printf("ERROR reading C_pTS: sprop_phenology_init()\n");
		errorCode = 20814;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 4. dividing line, sprop_phenology_init()\n");
		errorCode = 20815;
	}
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 4. dividing line, sprop_phenology_init()\n");
		errorCode = 20815;
	}
	/*********************************************************************************************************************/
	/* SOIL PARAMETERS */

	/*  soil depth */
	if (!errorCode && scan_value(sprop_file, &sprop->soildepth, 'd'))
	{
		printf("ERROR reading soildepth: sprop_phenology_init()\n");
		errorCode = 20816;
	}

	if (!errorCode && sprop->soildepth <= 0)
	{
		printf("ERROR in sprop data in SOIL file: negative or zero soildepth, sprop_phenology_init()\n");
		errorCode = 2081601;
	}

	/*  ratio of bare soil evaporation and pot.evaporation */
	if (!errorCode && scan_value(sprop_file, &sprop->soilevapLIM, 'd'))
	{
		printf("ERROR reading soilevapLIM: sprop_phenology_init()\n");
		errorCode = 20817;
	}

	/*  maximum height of pond water */
	if (!errorCode && scan_value(sprop_file, &sprop->pondmax, 'd'))
	{
		printf("ERROR reading maximum height of pond water: sprop_phenology_init()\n");
		errorCode = 20818;
	}

	/*  maximum height of pond water */
	if (!errorCode && scan_value(sprop_file, &sprop->q_soilstress, 'd'))
	{
		printf("ERROR reading q_soilstress: sprop_phenology_init()\n");
		errorCode = 20818;
	}

	/* runoff parameter (Campbell and Diaz)  */
	if (!errorCode && scan_value(sprop_file, &sprop->RCN_mes, 'd'))
	{
		printf("ERROR reading measured runoff curve number: sprop_phenology_init()\n");
		errorCode = 20819;
	}

	/* aerodynamic paramter (Wallace and Holwill, 1997) */
	if (!errorCode && scan_value(sprop_file, &sprop->aerodyn_resist, 'd'))
	{
		printf("ERROR reading aerodyn_resist: sprop_phenology_init()\n");
		errorCode = 20819;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 5. dividing line, sprop_phenology_init()\n");
		errorCode = 20820;
	}
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 5. dividing line, sprop_phenology_init()\n");
		errorCode = 20820;
	}
	/*********************************************************************************************************************/
	/* SOIL COMPOSITION AND CHARACTERISTIC VALUES (-9999: no measured data) */

	/* SAND array - mulilayer soil  */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->sand[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading percent sand in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20821;
		}
	}

	/* SILT array - mulilayer soil   */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->silt[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading percent silt in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20822;
		}
	}

	/* pH array - mulilayer soil */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->pH[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soil pH in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20823;
		}
	}

	/* measured bulk density    */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->BD_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading BD_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20824;
		}
	}

	/* measured critical VWC values - saturation    */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->vwc_sat_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading vwc_sat_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20825;
		}
	}

	/* measured critical VWC values - field capacity     */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->vwc_fc_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading vwc_fc_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20826;
		}
	}

	/* measured critical VWC values - wilting point    */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->vwc_wp_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading vwc_wp_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20827;
		}
	}

	/* measured critical VWC values - hygr. water    */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->vwc_hw_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading vwc_hw_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20828;
		}
	}

	/* measured drainage coeff */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->drain_coeff_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading drain_coeff_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20829;
		}
	}

	/* measured hydraulic conductivity  */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1)
			scanflag = 1;
		if (!errorCode && scan_array(sprop_file, &(sprop->conduct_sat_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading conduct_sat_mes in layer %i, sprop_phenology_init()\n", layer);
			errorCode = 20830;
		}
	}

	/* CONTROL - measured VWC values (all or none should to be set) */
	if (sprop->BD_mes[0] == DATA_GAP || sprop->vwc_sat_mes[0] == DATA_GAP ||
		sprop->vwc_fc_mes[0] == DATA_GAP || sprop->vwc_wp_mes[0] == DATA_GAP || sprop->vwc_hw_mes[0] == DATA_GAP)
	{

		if (sprop->BD_mes[0] == DATA_GAP && sprop->vwc_sat_mes[0] == DATA_GAP &&
			sprop->vwc_fc_mes[0] == DATA_GAP && sprop->vwc_wp_mes[0] == DATA_GAP && sprop->vwc_hw_mes[0] == DATA_GAP)
		{
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				sprop->BD_mes[layer] = DATA_GAP;
				sprop->vwc_sat_mes[layer] = DATA_GAP;
				sprop->vwc_fc_mes[layer] = DATA_GAP;
				sprop->vwc_wp_mes[layer] = DATA_GAP;
				sprop->vwc_hw_mes[layer] = DATA_GAP;
			}
		}
		else
		{
			printf("ERROR in measured soil data in SOI file\n");
			printf("All or none measured bulk density and critical VWC data should to be set by the user\n");
			errorCode = 20831;
		}
	}

	/* MULTILAYER SOIL CALCULATION  */
	if (!errorCode && multilayer_soilcalc(ctrl, sprop))
	{
		printf("\n");
		printf("ERROR in multilayer_soilcalc() in sprop_phenology_init.c\n");
		errorCode = 20832;
	}

	/* -------------------------------------------*/
	if (dofilecloseSOILPROP)
		fclose(sprop_file.ptr);

	// ?txt??SOS?EOS??-ZL
	if (!errorCode && scan_value(init, epc->SOS_EOS_name, 's'))
	{
		printf("ERROR reading outfile prefix: output_init(), output_init.c\n");
		errorCode = 21601;
	}

	return (errorCode);
}

int soilb_estimation(double sand, double silt, double *soil_b, double *vwc_sat, double *vwc_fc, double *vwc_wp, double *BD, double *RCN, int *soiltype)
{

	int st = -1;
	int errorCode = 0;

	double clay = 100 - sand - silt;

	double soilb_array[12] = {3.45, 4.11, 5.26, 6.12, 5.39, 4.02, 7.63, 7.71, 8.56, 9.22, 10.45, 12.46};
	double VWCsat_array[12] = {0.4, 0.42, 0.44, 0.46, 0.48, 0.49, 0.5, 0.505, 0.51, 0.515, 0.52, 0.525};
	double VWCfc_array[12] = {0.155, 0.190, 0.250, 0.310, 0.360, 0.380, 0.390, 0.405, 0.420, 0.435, 0.445, 0.460};
	double VWCwp_array[12] = {0.030, 0.050, 0.090, 0.130, 0.170, 0.190, 0.205, 0.220, 0.240, 0.260, 0.275, 0.290};
	double BD_array[12] = {1.6, 1.58, 1.56, 1.54, 1.52, 1.5, 1.48, 1.46, 1.44, 1.42, 1.4, 1.38};
	double RCN_array[12] = {50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72};

	if (silt + 1.5 * clay < 15) //sand
		st = 0;

	if ((silt + 1.5 * clay >= 15) && (silt + 2 * clay < 30)) //loamy_sand
		st = 1;

	if ((clay >= 7 && clay < 20) && (sand > 52) && (silt + 2 * clay >= 30)) //sandy_loam
		st = 2;

	if (clay < 7 && silt < 50 && silt + 2 * clay >= 30) //sandy_loam
		st = 2;

	if (clay >= 7 && clay < 27 && silt >= 28 && silt < 50 && sand <= 52) //loam
		st = 3;

	if (silt >= 50 && clay >= 12 && clay < 27) //silt loam
		st = 4;

	if (silt >= 50 && silt < 80 && clay < 12) //silt loam
		st = 4;

	if (silt >= 80 && clay < 12) //silt
		st = 5;

	if (clay >= 20 && clay < 35 && silt < 28 && sand > 45) //sandy_clay_loam
		st = 6;

	if (clay >= 27 && clay < 40 && sand > 20 && sand <= 45) //clay_loam
		st = 7;

	if (clay >= 27 && clay < 40 && sand <= 20) //silty clay_loam
		st = 8;

	if (clay >= 35 && sand > 45) //sandy_clay
		st = 9;

	if (clay >= 40 && silt >= 40) //silty_clay
		st = 10;

	if (clay >= 40 && sand <= 45 && silt < 40) //clay
		st = 11;

	*soiltype = st;

	if (st >= 0)
	{
		*soil_b = soilb_array[st];
		*vwc_sat = VWCsat_array[st];
		*vwc_fc = VWCfc_array[st];
		*vwc_wp = VWCwp_array[st];
		*BD = BD_array[st];
		*RCN = RCN_array[st];
	}
	else
	{
		printf("ERROR in soiltype founding\n");
		errorCode = 1;
	}

	return (errorCode);
}

int multilayer_soilcalc(control_struct *ctrl, soilprop_struct *sprop)
{
	int layer;

	double sand, silt, clay, psi_fc, psi_wp, psi_sat;
	double soil_b, BD, RCN, vwc_sat, vwc_fc, vwc_wp, vwc_hw, hydr_conduct_sat, hydr_diffus_sat, hydr_conduct_fc, hydr_diffus_fc;
	int errorCode = 0;

	soil_b = BD = RCN = psi_sat = vwc_sat = vwc_wp = vwc_hw = hydr_conduct_sat = hydr_diffus_sat = hydr_conduct_fc = hydr_diffus_fc = 0;

	/* -------------------------------------------------------------------------------------------------------------------------------*/

	/* 1. estimated soil water potential at hygroscopic water in MPa (1MPa = 10000cm)  (fc: pF = 2.5; wp: pF = 4.2) */
	sprop->psi_hw = pow(10, pF_hygroscopw) / (-10000);

	/* 2. Calculate the soil pressure-volume coefficients from texture data (soil water content, soil water potential and Clapp-Hornberger parameter) 
		    - different estimation methods 4: modell estimtaion, talajharomszog, measured data )*/

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{

		sand = sprop->sand[layer];
		silt = sprop->silt[layer];
		clay = 100 - sand - silt;

		sprop->clay[layer] = clay;

		/*  2.1 CONTROL to avoid negative data  */
		if (sand < 0 || silt < 0 || clay < 0)
		{
			if (!errorCode)
				printf("ERROR in site data in SOI file: negative sand/silt/clay data, sprop_phenology_init()\n");
			if (!errorCode)
				printf("Check values in initialization file.\n");
			errorCode = 1;
		}

		if (sand == 0 && silt == 0 && sprop->pH[layer] == 0)
		{
			if (!errorCode)
				printf("ERROR in site data in SOI file: user-defined sand and silt and pH data is zero in layer %i\n", layer);
			if (!errorCode)
				printf("Check values in initialization file.\n");
			errorCode = 1;
		}

		/* 2.2 saturation value of soil water potential */
		psi_sat = -(exp((1.54 - 0.0095 * sand + 0.0063 * silt) * log(10.0)) * 9.8e-5);

		if (soilb_estimation(sand, silt, &soil_b, &vwc_sat, &vwc_fc, &vwc_wp, &BD, &RCN, &ctrl->soiltype))
		{
			if (!errorCode)
			{
				printf("\n");
				printf("ERROR: soilb_estimation() for sprop_phenology_init \n");
			}
			errorCode = 1;
		}

		/* using the 4M model soil properties */
		if (sprop->vwc_sat_mes[layer] == (double)DATA_GAP)
		{
			/* VWC_fc, VWC_wp VWC_HW are the functions of VWC_sat and soil_b (new HW function from Fodor Nándor */
			vwc_hw = (0.01 * (100 - sand) + 0.06 * clay) / 100.;
			//vwc_hw = vwc_sat * (log(soil_b) / log(sprop->psi_hw/psi_sat));
			psi_fc = exp(vwc_sat / vwc_fc * log(soil_b)) * psi_sat;
			psi_wp = exp(vwc_sat / vwc_wp * log(soil_b)) * psi_sat;

			/* control for soil type with high clay content  */
			if (vwc_sat - vwc_fc < 0.001)
				vwc_fc = vwc_sat - 0.001;
			if (vwc_fc - vwc_wp < 0.001)
				vwc_wp = vwc_fc - 0.001;
			if (vwc_wp - vwc_hw < 0.001)
				vwc_hw = vwc_wp - 0.001;
		}
		/* measured soil water content at wilting point, field capacity and saturation in m3/m3 */
		else
		{

			vwc_sat = sprop->vwc_sat_mes[layer];
			vwc_fc = sprop->vwc_fc_mes[layer];
			vwc_wp = sprop->vwc_wp_mes[layer];
			vwc_hw = sprop->vwc_hw_mes[layer];
			BD = sprop->BD_mes[layer];

			psi_fc = exp(vwc_sat / vwc_fc * log(soil_b)) * psi_sat;
			psi_wp = exp(vwc_sat / vwc_wp * log(soil_b)) * psi_sat;
		}

		/* 2.4 CONTROL - measured VWC values: SAT>FC>WP>HW */
		if ((vwc_sat - vwc_fc) < 0.001 || (vwc_fc - vwc_wp) < 0.001 || (vwc_wp - vwc_hw) < 0.001 || vwc_hw < 0.001 || vwc_sat > 1.0)
		{
			if (!errorCode)
				printf("ERROR in measured VWC data in SOI file\n");
			if (!errorCode)
				printf("rules: VWC_sat > VWC_fc; VWC_fc > VWC_wp; VWC_wp > VWC_hw; VWC_sat <1.0, VWC_hw>0.01\n");
			errorCode = 1;
		}

		/* 2.5 hydr. conduct and diffusivity at saturation and field capacity(Cosby et al. 1984)*/

		/* if measured data is availabe (in cm/day) */
		if (sprop->conduct_sat_mes[layer] != DATA_GAP)
			hydr_conduct_sat = sprop->conduct_sat_mes[layer] / 100 / NSEC_IN_DAY;
		else
		{
			hydr_conduct_sat = 7.05556 * 1e-6 * pow(10, (-0.6 + 0.0126 * sand - 0.0064 * clay));
			hydr_conduct_sat = (50 * exp(-0.075 * clay) + 200 * exp(-0.075 * (100 - sand))) / 100 / NSEC_IN_DAY;
		}

		hydr_diffus_sat = (soil_b * hydr_conduct_sat * (-100 * psi_sat)) / vwc_sat;
		hydr_conduct_fc = hydr_conduct_sat * pow(vwc_fc / vwc_sat, 2 * soil_b + 3);
		hydr_diffus_fc = (((soil_b * hydr_conduct_sat * (-100 * psi_sat))) / vwc_sat) *
						 pow(vwc_fc / vwc_sat, soil_b + 2);

		sprop->soil_b[layer] = soil_b;
		sprop->BD[layer] = BD;
		sprop->psi_sat[layer] = psi_sat;
		sprop->psi_fc[layer] = psi_fc;
		sprop->psi_wp[layer] = psi_wp;
		sprop->vwc_sat[layer] = vwc_sat;
		sprop->vwc_fc[layer] = vwc_fc;
		sprop->vwc_wp[layer] = vwc_wp;
		sprop->vwc_hw[layer] = vwc_hw;
		sprop->hydr_conduct_sat[layer] = hydr_conduct_sat;
		sprop->hydr_diffus_sat[layer] = hydr_diffus_sat;
		sprop->hydr_conduct_fc[layer] = hydr_conduct_fc;
		sprop->hydr_diffus_fc[layer] = hydr_diffus_fc;

		if (layer == 0)
		{
			if (sprop->RCN_mes == DATA_GAP)
				sprop->RCN = RCN;
			else
				sprop->RCN = sprop->RCN_mes;
		}
	}
	return (errorCode);
}