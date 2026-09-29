/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang 2020-01-14 creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    sum.cpp
 * @brief		Encode the residual data
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		    2020-01-14
 *
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>

#include "gcoders/sum.h"
#include "gcoders/gcoder.h"
#include "gutils/gleo.h"

using namespace std;
namespace great
{

	/** @brief default constructor. */
	t_sum::t_sum()
	{
	}

	/**
	* @brief default constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_sum::t_sum(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_sum::constructor");

	}

	/**
	* @brief default constructor.
	*/
	t_sum::~t_sum()
	{

	}


	int t_sum::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		try
		{
			int i = 0;
			vector<string> sat;
			vector<string> site;

			vector<string> sat_slr;
			vector<string> site_slr;
			/*_ss.seekg(0, _ss.end);
			long len = (long)_ss.tellg();
			len = len - _ss_position;
			_ss.seekg(_ss_position, ios_base::beg);*/

			if (_ss_position == 0)
			{
				for (auto itt = _data.begin(); itt != _data.end(); itt++)
				{

					if (itt->second->id_type() == t_gdata::ALLSUM)
					{
						vector<string> sat_init = ((t_gallsum*)itt->second)->getSat();
						vector<string> site_init = ((t_gallsum*)itt->second)->getSite();
						map<t_gobscombtype, map<int, map<int, double>>> rmstable = ((t_gallsum*)itt->second)->getRMS();
						map<int, map<int, int>> nepotable = ((t_gallsum*)itt->second)->getNEPO();

						map<string, map<int, map<int, double>>> rmstable_slr = ((t_gallsum*)itt->second)->getRMS_SLR();
						map<int, map<int, int>> nepotable_slr = ((t_gallsum*)itt->second)->getNEPO_SLR();

						//add for slr
						for (int m = 0; m < sat_init.size(); m++)
						{
							if (!IS_SUPPORTED_LEO(strs2strl(sat_init[m])))
							{
								sat.push_back(sat_init[m]);
							}
							else
							{
								sat_slr.push_back(sat_init[m]);
							}
						}
						for (int m = 0; m < site_init.size(); m++)
						{
							if (IS_SUPPORTED_LEO(strs2strl(site_init[m])))
							{
								site.push_back(site_init[m]);
							}
							else
							{
								site_slr.push_back(site_init[m]);
							}
						}

						for (i = 0; i < 2; i++)
						{
							if (i == 0)
							{
								_ss << "RMS of each station-satellite pair, phase in mm" << endl;
							}
							else
							{
								_ss << "RMS of each station-satellite pair, range in dm" << endl;
							}

							_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

							for (int isat = 0; isat < sat.size(); isat++)
							{

								_ss << setw(7) << right << sat[isat];


							}
							_ss << endl;

							for (int ista = 0; ista < site.size(); ista++)
							{

								_ss << " " << setw(4) << right << site[ista];
								_ss << " ";
								for (int isat = 0; isat <= sat.size(); isat++)
								{

									if (i == 0)
									{
										//mm
										_ss << setw(7) << right << (int)(rmstable[t_gobscombtype("LC")][ista+1][isat] * 1000);
									}
									else
									{
										//dm
										_ss << setw(7) << right << (int)(rmstable[t_gobscombtype("PC")][ista+1][isat] * 10);
									}



								}
								_ss << endl;


							}

							_ss << " " << setw(4) << right << "SUMM";
							_ss << " ";
							for (int isat = 0; isat <= sat.size(); isat++)
							{

								if (i == 0)
								{
									//mm
									_ss << setw(7) << right << (int)(rmstable[t_gobscombtype("LC")][0][isat] * 1000);
								}
								else
								{
									//dm
									_ss << setw(7) << right << (int)(rmstable[t_gobscombtype("PC")][0][isat] * 10);
								}



							}
							_ss << endl;

							_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

							for (int isat = 0; isat < sat.size(); isat++)
							{

								_ss << setw(7) << right << sat[isat];


							}
							_ss << endl;
							_ss << endl;


						}

						_ss << "# of Epochs for each station-satellite pair" << endl;

						_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

						for (int isat = 0; isat < sat.size(); isat++)
						{

							_ss << setw(7) << right << sat[isat];


						}
						_ss << endl;

						for (int ista = 0; ista < site.size(); ista++)
						{

							_ss << " " << setw(4) << right << site[ista];
							_ss << " ";
							for (int isat = 0; isat <= sat.size(); isat++)
							{
								
							     _ss << setw(7) << right << (int)nepotable[ista+1][isat];
								

							}
							_ss << endl;


						}

						_ss << " " << setw(4) << right << "SUMM";
						_ss << " ";
						for (int isat = 0; isat <= sat.size(); isat++)
						{

							_ss << setw(7) << right << (int)nepotable[0][isat];



						}
						_ss << endl;

						_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

						for (int isat = 0; isat < sat.size(); isat++)
						{

							_ss << setw(7) << right << sat[isat];


						}
						_ss << endl;
						
						//add for SLR
						if (sat_slr.size() != 0)
						{
							_ss << endl;
							_ss << "SLR RMS of each station-satellite pair, in mm" << endl;

							_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

							for (int isat = 0; isat < sat_slr.size(); isat++)
							{

								_ss << setw(7) << right << sat_slr[isat];


							}
							_ss << endl;

							for (int ista = 0; ista < site_slr.size(); ista++)
							{

								_ss << " " << setw(4) << right << site_slr[ista];
								_ss << " ";
								for (int isat = 0; isat <= sat_slr.size(); isat++)
								{

									//mm
									_ss << setw(7) << right << (int)(rmstable_slr["SLR"][ista + 1][isat] * 1000);

								}
								_ss << endl;


							}

							_ss << " " << setw(4) << right << "SUMM";
							_ss << " ";
							for (int isat = 0; isat <= sat_slr.size(); isat++)
							{

								_ss << setw(7) << right << (int)(rmstable_slr["SLR"][0][isat] * 1000);

							}
							_ss << endl;

							_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

							for (int isat = 0; isat < sat_slr.size(); isat++)
							{

								_ss << setw(7) << right << sat_slr[isat];


							}
							_ss << endl;
							_ss << endl;

							_ss << "# of Epochs for each station-satellite pair" << endl;

							_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

							for (int isat = 0; isat < sat_slr.size(); isat++)
							{

								_ss << setw(7) << right << sat_slr[isat];


							}
							_ss << endl;

							for (int ista = 0; ista < site_slr.size(); ista++)
							{

								_ss << " " << setw(4) << right << site_slr[ista];
								_ss << " ";
								for (int isat = 0; isat <= sat_slr.size(); isat++)
								{

									_ss << setw(7) << right << (int)nepotable_slr[ista + 1][isat];


								}
								_ss << endl;


							}

							_ss << " " << setw(4) << right << "SUMM";
							_ss << " ";
							for (int isat = 0; isat <= sat_slr.size(); isat++)
							{

								_ss << setw(7) << right << (int)nepotable_slr[0][isat];

							}
							_ss << endl;

							_ss << " " << setw(4) << right << "NAME" << " " << setw(7) << right << "SUMM";

							for (int isat = 0; isat < sat_slr.size(); isat++)
							{

								_ss << setw(7) << right << sat_slr[isat];


							}
							_ss << endl;
						}
						

					}
				}
			}


		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_result::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_result::encode_data throw exception" << endl;
			}
			return -1;
		}

		int size = _fill_buffer(buff, sz);

		_mutex.unlock();  return size;

	}

}