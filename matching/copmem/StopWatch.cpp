/* ================================================================= *
 *  StopWatch.cpp : Main class                                       *
 *                                                                   *
 *  StopWatch: A tool to measure program performance                 *
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

#include "StopWatch.h"
CStopWatch::CStopWatch(){
	running = 0;
}


CStopWatch::~CStopWatch(){
}


void CStopWatch::start(){
	running = 1;
	elapsed.clear();
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::stop(){
	t2  = std::chrono::steady_clock::now();
	running = 0;
	double el = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
	elapsed.push_back(el);
	return el;
}


void CStopWatch::resume(){
	running = 1;
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::totalTime(){
	double sum = 0;
	for (std::vector<double>::iterator it = elapsed.begin(); it != elapsed.end(); ++it) {
		sum += *it;
	}
	return sum;
}
