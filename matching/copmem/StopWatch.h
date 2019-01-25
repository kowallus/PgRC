/* ================================================================= *
 *  StopWatch.h : Header file with supporting class definitions      *
 *                                                                   *
 *  StopWatch: A tool to measure program performance
 *                                                                   *
 *  Copyright (c) 2018, Szymon Grabowski and Wojciech Bieniecki      *
 *  All rights reserved                                              *
 *                                                                   * 
 *  This program is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public License as   *
 *  published by the Free Software Foundation, either version 3 of   *
 *  the License, or (at your option) any later version.              *
 *                                                                   *
 *  This program is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 *  GNU General Public License for more details.                     *
 *                                                                   *
 *  You should have received a copy of the GNU General Public        *
 *  License along with this program.                                 *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'license', which is part of this source code package.       *
 * ================================================================= */

#pragma once
#include <chrono>
#include <vector>
class CStopWatch
{
private:
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::steady_clock::time_point t_temp;
	int running;
	std::vector<double> elapsed;

public:
	CStopWatch();
	~CStopWatch();
	void start();
	double stop();
	void resume();
	double totalTime();
};

