/* 
GSI_calculation.c
Calculation of pheonological state of vegetation (onset and offset day) from Tmin, VPD, daylength and cumulative Tavg 
based on literure (Jolly et al, 2005) and own method. The goal is to replace prephenology.c subroutines becuse the calulation
of the model-defined onset and offset day does not work correctly

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
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"
#include "misc_func.h"

// 函数声明ivalue and 寻找第二大的数字-20201204-ZL
double ivalue(double value, double limit1, double limit2);
int search_second_max(int **array, int n);
int search_min(int **array, int n);

// Chen的秋季物候，直接添加到GSI_calculation中，参考文献：A new process-based model for predicting autumn phenology: How is leaf senescence controlled by photoperiod and temperature coupling?

// 控制是否调试-ZL
// #define DEBUG_GSI 1

int GSI_calculation(const metarr_struct *metarr, const siteconst_struct *sitec, epconst_struct *epc,
					phenarray_struct *phenarr, control_struct *ctrl)

{
	int errorCode = 0;
	int ny, yday, back;

	int firstdayLP = 240; /* theoretically first day of litterfall */

	/* flags and counters */
	int onday_flag = 0;
	int offday_flag = 0;
	int onday = 0;
	int offday = 0;
	int nyears = ctrl->simyears;
	int n_yday = nDAYS_OF_YEAR;

	/*  enviromental conditions taken account to calculate onset and offset days */
	double tmax_act, tmin_act, tavg_act, vpd_act, dayl_act, heatsum_act;

	/* threshold limits for each variable, between assuming that phenological activity varied linearly from inactive to unconstrained */
	double base_temp = epc->base_temp;
	double heatsum_limit1 = epc->heatsum_limit1;
	double heatsum_limit2 = epc->heatsum_limit2;
	double tmin_limit1 = epc->tmin_limit1;
	double tmin_limit2 = epc->tmin_limit2;
	double vpd_limit1 = epc->vpd_limit1;
	double vpd_limit2 = epc->vpd_limit2;
	double dayl_limit1 = epc->dayl_limit1;
	double dayl_limit2 = epc->dayl_limit2;
	int n_moving_avg = epc->n_moving_avg - 1;  /* moving averages are calculated from indicatiors 
											      to avoid the effects of single extreme events; -1: number to index */
	double epc_limit_SGS = epc->GSI_limit_SGS; /* when GSI index fisrt time greater that limit -> start of the growing season */
	double epc_limit_EGS = epc->GSI_limit_EGS; /* when GSI index fisrt time less that limit -> end of the growing season */

	/* indexes for the different variables and total index (multiplication of  partial indexes)*/
	double tmin_index = 0;
	double vpd_index = 0;
	double dayl_index = 0;
	double epc_index = 0;

	/* at the presence of snow cover no vegetation period (calculating snow cover from precipitation, Tavg and srad) */
	//DONE 增加了prcp_sum/mean_30 + index, prcp/sw-limits1/2-20201203(后续通过EPC进行读取)-ZL prcp_sum删除了[后期]-20210512
	double prcp_limit1 = epc->prcp_limit1 / 10.0;
	double prcp_limit2 = epc->prcp_limit2 / 10.0;
	double sw_limit1 = epc->sw_limit1;
	double sw_limit2 = epc->sw_limit2;
	// double prcp_limit1 = 0.2 / 10.0, prcp_limit2 = 2.0 / 10.0, sw_limit1 = 170, sw_limit2 = 320;
	int second_max, min;

	// 增加秋季物候的控制-并读取参数-ZL
	int start_day_find;
	int start_day;
	// 1 for T, 2 for P
	int start_option = epc->start_option;
	double start = epc->start;
	double a = epc->a;
	double b = epc->b;
	double S_threshold = epc->S_threshold;
	int TPM_use = epc->TPM_use_id;
	double Rsen = 0, Rsen_tmp = 0;
	// TPM在一年内是否完成
	int TPM_done = 0;

	//增加了prcp_mean_30,prcp_index,sw_index-0309-ZL
	double snowcover, snowcover_limit, prcp_act, srad_act, prcp_mean_30, prcp_index, sw_index;
	double rn, tmelt, rmelt, melt, snow_loss, snow_plus, rsub;
	double albedo_sw = sitec->albedo_sw;

	/* to calculate moving average from total index values */
	double epc_index_SUM = 0;
	double epc_index_avg = 0;
	double epc_index_total = 0;

	onday_flag = 0;
	offday_flag = 1;

	/* limitation of snow cover */
	snowcover = 0;

	snowcover_limit = epc->snowcover_limit;

	if (!errorCode)
	{
		/* allocate space for the onday_arr and offday_arr: first column - year, second column: day*/
		phenarr->onday_arr = (int **)malloc(nyears * sizeof(int *));
		phenarr->offday_arr = (int **)malloc(nyears * sizeof(int *));

		for (ny = 0; ny < nyears; ny++)
		{
			phenarr->onday_arr[ny] = (int *)malloc(2 * sizeof(int));
			phenarr->offday_arr[ny] = (int *)malloc(2 * sizeof(int));
			phenarr->onday_arr[ny] = (int *)malloc(2 * sizeof(int));
			phenarr->offday_arr[ny] = (int *)malloc(2 * sizeof(int));
		}
#ifdef ERROR_PRINT
		if (!phenarr->onday_arr)
		{
			printf("\n");
			printf("ERROR allocating for onday_arr, prephenology()\n");
			errorCode = 1;
		}

		if (!phenarr->offday_arr)
		{
			printf("\n");
			printf("ERROR allocating for offday_arr, prephenology()\n");
			errorCode = 1;
		}
#endif
	}

	for (ny = 0; ny < nyears; ny++)
	{
		onday = 0;
		offday = 0;

		// 20210518-ZL
		//start_day = 0;
		start_day_find = 0;
		TPM_done = 0;
		Rsen_tmp = 0;
		Rsen = 0;

#ifdef DEBUG_GSI
		// 清空prcp_sum_30-20210106-ZL
		// for(int day_prcp = 0; day_prcp <= 364; day_prcp++){

		// }

		//测试-为什么GSI一直输出失败在大部分年份-20210104-ZL
		// if (ny == 37)
		// {
		FILE *fp;
		fp = fopen("E:/GMUSO-DESKTOP/6.1/biome-bgc-muso_6.1_Linux_src/biome-bgc-muso_6.1_Linux_src/OUTPUT/GSI-Sun.txt", "w");
		// }
#endif

		for (yday = 0; yday < nDAYS_OF_YEAR; yday++)
		{
			/* ******************************************************************* */
			/* 1. calculation of snow loss and plus*/

			tmax_act = metarr->tmax[ny * n_yday + yday];
			tmin_act = metarr->tmin[ny * n_yday + yday];
			tavg_act = (tmax_act + tmin_act) / 2.;
			prcp_act = metarr->prcp[ny * n_yday + yday];
			srad_act = metarr->swavgfd[ny * n_yday + yday];
			dayl_act = metarr->dayl[ny * n_yday + yday];

			/* canopy transmitted radiaiton: convert from W/m2 --> KJ/m2/d */
			rn = srad_act * (1.0 - albedo_sw) * dayl_act * sn_abs * 0.001;

			/* 1.a. snow loss: melting and sublim. */
			if (tavg_act > 0.0)
			{ /* temperature melt from snowpack */
				tmelt = tcoef * tavg_act;
				/* radiaiton melt from snowpack */
				rmelt = rn / lh_fus;
				melt = tmelt + rmelt;
				/* snowmelt decreases actual snow pool */
				snow_loss = melt;
			}
			else
			{ /* sublimation from snowpack */
				rsub = rn / lh_sub;
				/* sublimation decreases actual snow pool */
				snow_loss = rsub;
			}
			if (snow_loss > snowcover)
				snow_loss = snowcover;

			/* 1.b. snow plus: melting and sublim. */
			if (tavg_act < 0.0)
			{
				snow_plus = prcp_act * 10;
			}
			else
				snow_plus = 0;

			snowcover -= snow_loss;
			snowcover += snow_plus;

			/* ******************************************************************* */
			/* 2. calculation of indexes (based on moving averages of evironmental parameters) */

			if (yday < n_moving_avg)
			{
				heatsum_act = 0;
				tmin_index = 0;
				// 增加计算30天降水滑动平均-20201203-ZL
				prcp_index = 0;
				epc_index_avg = 0;
				epc_index_total = 0;
			}
			else
			{
				epc_index_SUM = 0;
				heatsum_act = 0;
				for (back = 0; back <= n_moving_avg; back++)
				{
					/* search actual values of variables */
					tmax_act = metarr->tmax[ny * n_yday + yday - (n_moving_avg - back)];
					tmin_act = metarr->tmin[ny * n_yday + yday - (n_moving_avg - back)];
					tavg_act = (tmax_act + tmin_act) / 2.;
					vpd_act = metarr->vpd[ny * n_yday + yday - (n_moving_avg - back)];
					dayl_act = metarr->dayl[ny * n_yday + yday - (n_moving_avg - back)];
					// 考虑srad-2021204-ZL
					srad_act = metarr->swavgfd[ny * n_yday + yday - (n_moving_avg - back)];
					prcp_mean_30 = metarr->prcp30_ra[ny * n_yday + yday - (n_moving_avg - back)];

					/* ******************************************************************* */
					/* 1.1 calculation of heatsum regarding to the basic temperature and n_moving_avg long period */
					if (tavg_act > base_temp)
					{
						heatsum_act += (tavg_act - base_temp);
					}

					/* ******************************************************************* */
					/* 1.2. calculation of indexes regarding to the different variables */
					/* !!!!!!!!!!!!!!!!!!!!!!!!!!  A: tmin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
					tmin_index = ivalue(tmin_act, tmin_limit1, tmin_limit2);

					/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!  B: vpd !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
					vpd_index = ivalue(vpd_act, vpd_limit1, vpd_limit2);

					/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  C: dayl !!!!!!!!!!!!!!!!!!!!!!!!!!! */
					dayl_index = ivalue(dayl_act, dayl_limit1, dayl_limit2);

					// 计算Swad-20201204-ZL
					sw_index = ivalue(srad_act, sw_limit1, sw_limit2);

					// 计算降水的限制-20210106-ZL
					prcp_index = ivalue(prcp_mean_30, prcp_limit1, prcp_limit2);

					// prcp代替VPD，20201204-ZL
					// 原语句 epc_index = tmin_index * vpd_index * dayl_index;
					epc_index = tmin_index * prcp_index * dayl_index * sw_index;
					epc_index_SUM += epc_index;

				} /* endfor - calculating indexes for the n_moving_average long period  */

				epc_index_avg = epc_index_SUM / (n_moving_avg + 1);
				// 忽略heatsum的影响，20201204-ZL
				// 原语句e pc_index_total = epc_index_avg * heatsum_index;
				epc_index_total = epc_index_avg;

			} /* endelse - calculating indexes */

			/* 3. calculation of autumn phenology-ZL */
			if (TPM_use == 1 && TPM_done == 0)
			{
				if (start_day_find == 0)
				{
					// option选1，为T，2，为P
					if (start_option == 1)
					{
						if (yday >= 200 && tmin_act <= start)
						{
							//start_day = yday;
							// 已经找到
							start_day_find = 1;
						}
					}
					else
					{
						if (yday >= 173 && (dayl_act / 3600) <= start)
						{
							//start_day = yday;
							start_day_find = 1;
						}
					}
				}
				else if (start_day_find == 1)
				{
					Rsen_tmp = 1 / (1 + exp(a * ((tmin_act * (dayl_act / 3600.0)) - b)));
					Rsen += Rsen_tmp;
				}
			}

			//监测GSI为什么失败-20210104-ZL
			// fprintf(fp, "%f ", tmin_act);
#ifdef DEBUG_GSI
			if (ny == 37)
			{
				// FILE *fp;
				// fp = fopen("E:/GMUSO-DESKTOP/6.1/biome-bgc-muso_6.1_Linux_src/biome-bgc-muso_6.1_Linux_src/OUTPUT/GSI-Sun.txt", "w");

				fprintf(fp, "%f ", epc_index);
				fprintf(fp, "%f ", epc_index_avg);
				fprintf(fp, "%f ", prcp_index);
				// fprintf(fp, "%f ", tmin_index);
				// fprintf(fp, "%f ", epc_index_total);
				fprintf(fp, "%f ", prcp_act);
				// fprintf(fp, "%f ", prcp_sum_30);
				fprintf(fp, "%f ", prcp_mean_30);
				fputc('\r', fp);
				fputc('\n', fp);

				// fclose(fp);
			}
#endif
			//测试-为什么GSI一直输出失败在大部分年份-20210104-ZL

			// 判断春季物候
			if (onday_flag == 0 && offday_flag == 1 && epc_index_total > epc_limit_SGS && yday < firstdayLP && snowcover <= snowcover_limit)
			{
				onday_flag = 1;
				offday_flag = 0;
				onday = yday;

				if (ctrl->onscreen && ctrl->spinup == 0)
					printf("Year: %i\n", ctrl->simstartyear + ny);
				if (ctrl->onscreen && ctrl->spinup == 0)
					printf("onday  - %i\n", onday);
			}

			// 判断秋季物候，如果使用TPM_use，那么就不计算秋季物候了(GSI)，依旧使用理论的起始点240-ZL
			if (TPM_use == 1)
			{
				if (yday > firstdayLP && onday_flag == 1 && offday_flag == 0 && Rsen >= S_threshold)
				{
					onday_flag = 0;
					offday_flag = 1;
					offday = yday;
					// 将TPM_done改为1，认为完成了工作-ZL
					TPM_done = 1;

					if (ctrl->onscreen && ctrl->spinup == 0)
						printf("offday - %i\n", offday);
				}
			}
			else
			{
				if (onday_flag == 1 && offday_flag == 0 && epc_index_total < epc_limit_EGS && yday > firstdayLP)
				{
					onday_flag = 0;
					offday_flag = 1;
					offday = yday;

					if (ctrl->onscreen && ctrl->spinup == 0)
						printf("offday - %i\n", offday);
				}
			}

			/* if vegetation period has not ended until the last day of year, the offday is equal to the last day of year */
			if (yday == nDAYS_OF_YEAR - 1 && offday == 0)
			{ /* if vegetation period has not began */
				if (onday_flag == 0)
				{
					onday = yday - 2;
					ctrl->vegper_flag = 1;
				}
				onday_flag = 0;
				offday_flag = 1;
				// 如果没有找到秋季物候的时间，则TPMT和TPMP都改为0，这样是为了在调参的过程中不影响最优值-ZL-20210916
				// offday = yday;
				offday = 0;
			}

		} /* endfor - simdays */

#ifdef DEBUG_GSI
		//监测GSI为什么失败-20210104-ZL
		// if (ny == 37)
		// {
		fclose(fp);
		// }
		//监测GSI为什么失败-20210104-ZL
#endif
		// 依照Sun的方法，将返青期向前5天，将枯黄期向后推迟10天-ZL-20210916
		// 如果是使用TPM方法的话，那么生长季结束日期就不减了-ZL-20210916
		if (TPM_use == 1)
		{
			phenarr->onday_arr[ny][0] = ny + ctrl->simstartyear;
			phenarr->onday_arr[ny][1] = onday - 5;

			phenarr->offday_arr[ny][0] = ny + ctrl->simstartyear;
			phenarr->offday_arr[ny][1] = offday;
		}
		else
		{
			phenarr->onday_arr[ny][0] = ny + ctrl->simstartyear;
			phenarr->onday_arr[ny][1] = onday - 5;

			phenarr->offday_arr[ny][0] = ny + ctrl->simstartyear;
			phenarr->offday_arr[ny][1] = offday + 10;
		}

	} /* endfor - simyears */

	//替换phenarr中异常值-ZL
	//在进行物候调参的过程中，这一输出会影响调参结果，因此注释掉，在正常的运行过程中需要打开-ZL-20210719-[Phenology]
	/* 	second_max = search_second_max(phenarr->onday_arr, nyears);
	min = search_min(phenarr->offday_arr, nyears);

	for (int i = 0; i < nyears; i++)
	{
		if (phenarr->offday_arr[i][1] == 364)
		{
			phenarr->offday_arr[i][1] = min + 15;
		}

		if (phenarr->onday_arr[i][1] == 362)
		{
			phenarr->onday_arr[i][1] = second_max - 15;
		}
	} */

	return (errorCode);

} /* end - subroutine */

// 函数定义ivalue-20201204-ZL
double ivalue(double value, double limit1, double limit2)
{
	double index = 0;

	if (value < limit1)
	{
		index = 0;
	}
	else
	{
		if (value < limit2)
		{
			index = (value - limit1) / (limit2 - limit1);
		}
		else
		{
			index = 1;
		}
	}

	return (index);
}

// 函数定义search_second_max_min-20210107-ZL
int search_second_max(int **array, int n)
{
	int max1, max2, i;
	max1 = array[0][1];
	for (i = 1; i < n; i++)
	{
		if (array[i][1] > max1)
		{
			max2 = max1;
			max1 = array[i][1];
		}
		else if (array[i][1] == max1)
		{
			/* do-nothing */
		}
		else
		{
			if (i == 1)
				max2 = array[i][1];
			else if (array[i][1] > max2)
				max2 = array[i][1];
		}
	}
	return max2;
}

// 函数定义search_min-20210107-ZL
int search_min(int **array, int n)
{
	int min, i;

	min = array[0][1];
	for (i = 1; i < n; i++)
	{
		if (min > array[i][1])
			min = array[i][1];
	}

	return min;
}
