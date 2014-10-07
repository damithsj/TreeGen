//calcSubBranchCoords.cpp
#include "calcSubBranchCoords.h"
#include <iostream>
#include <cmath>

using namespace std;

void calcSubBranchCoords(double xA, double yA, double zA, double xB, double yB, double zB, double& subX, double& subY, double& subZ, double angle, double incAngle, double scaleFactor)
{
// calculate u, v
double u = xB - xA;
double v = yB - yA;
double w = zB - zA;

double L = sqrt(u*u + v*v + w*w);
double M = sqrt(u*u + v*v);

//calculate sub branch coords using Honda model
subX = xB + scaleFactor * (u*cos(angle) - (L*v*cos(incAngle) + u*w*sin(incAngle))*sin(angle)/M);
subY = yB + scaleFactor * (v*cos(angle) + (L*u*cos(incAngle) - v*w*sin(incAngle))*sin(angle)/M);
subZ = zB + scaleFactor * (w*cos(angle) + M*sin(incAngle)*sin(angle));
}

void calcSubBranchCoords1(double xA, double yA, double xB, double yB, double& subX, double& subY, double angle, double branchLength)
{
	double u = xB - xA;
	double v = yB - yA;

	double L1 = abs(sqrt(u*u + v*v));
	subX = xB + (branchLength/L1)*(xB*cos(angle) + yB*sin(angle));
	subY = yB + (branchLength/L1)*(yB*cos(angle) - xB*sin(angle));
}
