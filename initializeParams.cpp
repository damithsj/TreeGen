/* ----------------------------------------------------------------------------
libconfig - A library for processing structured configuration files
Copyright (C) 2005-2010  Mark A Lindner

This file is part of libconfig.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 2.1 of
the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Library General Public
License along with this library; if not, see
<http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <cmath>
#include "initializeParams.h"



// definitionas of extern variables. otherwise they cannot use here.
#pragma region Externs defined in initializeParams.h
int constants::VERT_LENGTH;
int constants::HORZ_LENGTH;
int constants::DEPTH_LENGTH;

double constants::CUTOFF;
double constants::ORP;

double constants::MIN_X;
double constants::MAX_X;
double constants::MIN_Y;
double constants::MAX_Y;
double constants::MIN_Z;
double constants::MAX_Z;

double constants::NETA;

int constants::MAX_LEVELS;

double constants::THETA;
double constants::BETA;
double constants::DELTA;

double constants::BR_SCALE_FACTOR;
double constants::BUD_RADIUS;
double constants::BUD_LENGTH;

// Random factors added to growth parameters
double constants::R_THETA;
double constants::R_BETA;
double constants::R_DELTA;
double constants::R_SF; // scale factor

bool constants::SHOW_TREE_GROWTH;

double constants::SOURCE_X_START;
double constants::SOURCE_X_END;
double constants::SOURCE_Y_START;
double constants::SOURCE_Y_END;
double constants::SOURCE_Z_START;
double constants::SOURCE_Z_END;

std::vector<obstacle> constants::OBSTACLES;

#pragma endregion


using namespace std;
using namespace libconfig;
using namespace constants;

int initializeParameters()
{


#pragma region Read the file. If there is an error, report it and exit
	Config cfg;

	// Read the file. If there is an error, report it and exit.
	try
	{
		cfg.readFile("../conf/initialize-params.cfg");
	}
	catch(const FileIOException &fioex)
	{
		std::cerr << "I/O error while reading file." << std::endl;
		return(EXIT_FAILURE);
	}
	catch(const ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
			<< " - " << pex.getError() << std::endl;
		return(EXIT_FAILURE);
	}
#pragma endregion
	// assign the root
	const Setting& root = cfg.getRoot();
	
	// start lookup the parameters

#pragma region Definitions of the Grid
	// VERT_LENGTH
	try
	{
		constants::VERT_LENGTH = cfg.lookup("VERT_LENGTH");
		//cout << "VERT_LENGTH: " << constants::VERT_LENGTH << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'VERT_LENGTH' setting in configuration file." << endl;
	}

	// HORZ_LENGTH
	try
	{
		constants::HORZ_LENGTH = cfg.lookup("HORZ_LENGTH");
		//cout << "HORZ_LENGTH: " << constants::HORZ_LENGTH << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'HORZ_LENGTH' setting in configuration file." << endl;
	}

	// DEPTH_LENGTH
	try
	{
		constants::DEPTH_LENGTH = cfg.lookup("DEPTH_LENGTH");
		//cout << "DEPTH_LENGTH: " << constants::DEPTH_LENGTH << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'DEPTH_LENGTH' setting in configuration file." << endl;
	}

#pragma endregion

#pragma region SOR parameters
	// CUTOFF
	try
	{
		constants::CUTOFF = cfg.lookup("CUTOFF");
		//cout << "CUTOFF: " << constants::CUTOFF << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'CUTOFF' setting in configuration file." << endl;
	}

	// ORP
	try
	{
		constants::ORP = cfg.lookup("ORP");
		//cout << "ORP: " << constants::ORP << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'ORP' setting in configuration file." << endl;
	}

#pragma endregion

#pragma region Tree growth boundary
	// MIN_X
	try
	{
		constants::MIN_X = cfg.lookup("MIN_X");
		//cout << "MIN_X: " << constants::MIN_X << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MIN_X' setting in configuration file." << endl;
	}

	// MAX_X
	try
	{
		constants::MAX_X = cfg.lookup("MAX_X");
		//cout << "MAX_X: " << constants::MAX_X << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MAX_X' setting in configuration file." << endl;
	}

	// MIN_Y
	try
	{
		constants::MIN_Y = cfg.lookup("MIN_Y");
		//cout << "MIN_Y: " << constants::MIN_Y << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MIN_Y' setting in configuration file." << endl;
	}

	// MAX_Y
	try
	{
		constants::MAX_Y = cfg.lookup("MAX_Y");
		//cout << "MAX_Y: " << constants::MAX_Y << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MAX_Y' setting in configuration file." << endl;
	}

	// MIN_Z
	try
	{
		constants::MIN_Z = cfg.lookup("MIN_Z");
		//cout << "MIN_Z: " << constants::MIN_Z << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MIN_Z' setting in configuration file." << endl;
	}

	// MAX_Z
	try
	{
		constants::MAX_Z = cfg.lookup("MAX_Z");
		//cout << "MAX_Z: " << constants::MAX_Z << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MAX_Z' setting in configuration file." << endl;
	}
#pragma endregion

#pragma region Neta; Max levels
	// NETA
	try
	{
		constants::NETA = cfg.lookup("NETA");
		//cout << "NETA: " << constants::NETA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'NETA' setting in configuration file." << endl;
	}

	// MAX_LEVELS
	try
	{
		constants::MAX_LEVELS = cfg.lookup("MAX_LEVELS");
		//cout << "MAX_LEVELS: " << constants::MAX_LEVELS << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'MAX_LEVELS' setting in configuration file." << endl;
	}
#pragma endregion

#pragma region Branch angles
	//THETA
	try
	{
		double __theta;
		__theta = cfg.lookup("THETA");
		constants::THETA  = __theta*PI/180;
		//cout << "THETA: " << constants::THETA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'THETA' setting in configuration file." << endl;
	}

	//BETA
	try
	{
		double __beta;
		__beta = cfg.lookup("BETA");
		constants::BETA  = __beta*PI/180;
		//cout << "BETA: " << constants::BETA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'BETA' setting in configuration file." << endl;
	}

	//DELTA
	try
	{
		double __delta;
		__delta = cfg.lookup("DELTA");
		constants::DELTA  = __delta*PI/180;
		//cout << "DELTA: " << constants::DELTA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'DELTA' setting in configuration file." << endl;
	}
#pragma endregion

#pragma region Random variations
	//R_THETA
	try
	{
		constants::R_THETA = cfg.lookup("R_THETA");
		//cout << "R_THETA: " << constants::R_THETA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'R_THETA' setting in configuration file." << endl;
	}

	//R_BETA
	try
	{
		constants::R_BETA = cfg.lookup("R_BETA");
		//cout << "R_BETA: " << constants::R_BETA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'R_BETA' setting in configuration file." << endl;
	}

	//R_DELTA
	try
	{
		constants::R_DELTA = cfg.lookup("R_DELTA");
		//cout << "R_DELTA: " << constants::R_DELTA << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'R_DELTA' setting in configuration file." << endl;
	}

	//R_SF
	try
	{
		constants::R_SF = cfg.lookup("R_SF");
		//cout << "R_SF: " << constants::R_SF << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'R_SF' setting in configuration file." << endl;
	}
#pragma endregion

#pragma region Scale factor; Show Tree Growth 
	// BR_SCALE_FACTOR
	try
	{
		constants::BR_SCALE_FACTOR = cfg.lookup("BR_SCALE_FACTOR");
		//cout << "BR_SCALE_FACTOR: " << constants::BR_SCALE_FACTOR << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'BR_SCALE_FACTOR' setting in configuration file." << endl;
	}

	// BUD_RADIUS
	// commented this and added along with BUD_LENGTH
	//try
	//{
	//	constants::BUD_RADIUS = cfg.lookup("BUD_RADIUS");
	//	//cout << "BUD_RADIUS: " << constants::BUD_RADIUS << endl << endl;
	//}
	//catch(const SettingNotFoundException &nfex)
	//{
	//	cerr << "No 'BUD_RADIUS' setting in configuration file." << endl;
	//}

	// SHOW_TREE_GROWTH
	try
	{
		constants::SHOW_TREE_GROWTH = cfg.lookup("SHOW_TREE_GROWTH");
		//cout << "SHOW_TREE_GROWTH: " << constants::SHOW_TREE_GROWTH << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'SHOW_TREE_GROWTH' setting in configuration file." << endl;
	}
#pragma endregion
	// BUD_LENGTH - calculate using above values
	constants::BUD_LENGTH = sqrt(((MAX_Y - MIN_Y)*(MAX_Y - MIN_Y))/(VERT_LENGTH*VERT_LENGTH) + ((MAX_X - MIN_X)*(MAX_X - MIN_X))/(HORZ_LENGTH*HORZ_LENGTH) + ((MAX_Z - MIN_Z)*(MAX_Z - MIN_Z))/(DEPTH_LENGTH*DEPTH_LENGTH))*1.5;
	constants::BUD_RADIUS = BUD_LENGTH/8;
#pragma region Source configurations
	// Source configurations	
	try
	{
		const Setting &XSOURCE = root["LIGHT_SOURCE"]["X"];	
		XSOURCE.lookupValue("START", constants::SOURCE_X_START);
		XSOURCE.lookupValue("END", constants::SOURCE_X_END);

		const Setting &YSOURCE = root["LIGHT_SOURCE"]["Y"];	
		YSOURCE.lookupValue("START", constants::SOURCE_Y_START);
		YSOURCE.lookupValue("END", constants::SOURCE_Y_END);

		const Setting &ZSOURCE = root["LIGHT_SOURCE"]["Z"];	
		ZSOURCE.lookupValue("START", constants::SOURCE_Z_START);
		ZSOURCE.lookupValue("END", constants::SOURCE_Z_END);

		//return(EXIT_SUCCESS);
	}
	catch(const SettingNotFoundException &nfex)
	{
		// Ignore.
	}
#pragma endregion

//#pragma region Obstacle configurations
//	// Obstacle configurations
//	//const Setting& root = cfg.getRoot();
//	try
//	{
//		const Setting &XOBS = root["OBSTACLE"]["X"];	
//		XOBS.lookupValue("START", constants::OBS_X_START);
//		XOBS.lookupValue("END", constants::OBS_X_END);
//
//		const Setting &YOBS = root["OBSTACLE"]["Y"];	
//		YOBS.lookupValue("START", constants::OBS_Y_START);
//		YOBS.lookupValue("END", constants::OBS_Y_END);
//
//		const Setting &ZOBS = root["OBSTACLE"]["Z"];	
//		ZOBS.lookupValue("START", constants::OBS_Z_START);
//		ZOBS.lookupValue("END", constants::OBS_Z_END);
//
//		//return(EXIT_SUCCESS);
//	}
//	catch(const SettingNotFoundException &nfex)
//	{
//		// Ignore.
//	}
//#pragma endregion

#pragma region Obstacles
	try
	{
		const Setting &OBSTACLES = root["OBSTACLES"];
		int count = OBSTACLES.getLength();

		for(int i = 0; i < count; ++i)
		{
			const Setting &OBSTACLE = OBSTACLES[i];
			double XStart, XEnd, YStart, YEnd, ZStart, ZEnd;
			obstacle obs;
			
			OBSTACLE.lookupValue("XSTART", obs.obsXStart);
			OBSTACLE.lookupValue("XEND", obs.obsXEnd);

			OBSTACLE.lookupValue("YSTART", obs.obsYStart);
			OBSTACLE.lookupValue("YEND", obs.obsYEnd);

			OBSTACLE.lookupValue("ZSTART", obs.obsZStart);
			OBSTACLE.lookupValue("ZEND", obs.obsZEnd);
			constants::OBSTACLES.push_back(obs);	
		}
	}
	catch(const SettingNotFoundException &nfex)
	{
		// Ignore.
	}
#pragma endregion


	return(EXIT_SUCCESS);


}