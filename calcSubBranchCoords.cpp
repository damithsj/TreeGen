//calcSubBranchCoords.cpp
#include "calcSubBranchCoords.h"
#include <iostream>
#include <cmath>

using namespace std;
//this function generates the end coordinates of the child branch using the mother branch
//Source for the equation:
//Control of Development in the Bifurcating Branch System of Tabebuia rosea: A Computer Simulation
//Rolf Borchert; Hisao Honda
//Botanical Gazette, Vol. 145, No. 2. (Jun., 1984), pp. 184-195
void calcSubBranchCoords(double xA, double yA, double zA, double xB, double yB, double zB, double& subX, double& subY, double& subZ, double angle, double incAngle, double scaleFactor)
{
// calculate u, v
double u = xB - xA;
double v = yB - yA;
double w = zB - zA;

double L = sqrt(u*u + v*v + w*w);
double M = sqrt(u*u + v*v);

//calculate sub branch coords using Honda I model
subX = xB + scaleFactor * (u*cos(angle) - (L*v*cos(incAngle) + u*w*sin(incAngle))*sin(angle)/M);
subY = yB + scaleFactor * (v*cos(angle) + (L*u*cos(incAngle) - v*w*sin(incAngle))*sin(angle)/M);
subZ = zB + scaleFactor * (w*cos(angle) + M*sin(incAngle)*sin(angle));
}

