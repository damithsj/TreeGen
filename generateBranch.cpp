#include "initializeParams.h"
#include "calcSubBranchCoords.h"
#include "mathUtil.h"

#include <cmath>
#include <iostream>

using namespace std;
using namespace constants;

#define PI 3.14159265

//outout the new branch point giving the mother branch points (start and end)
point generateBranch(point A, point B, char branch_type)
{
	point sub_;
	double subX;
	double subY;
	double subZ;
	double scaleFactor = 1;

	double theta_ = THETA + THETA*(-1 + 2*unifRand())*R_THETA; // variation +/- R_THETA
	double beta_ = BETA + BETA*(-1 + 2*unifRand())*R_BETA; // variation +/- R_BETA
	double delta_ = DELTA + DELTA*(-1 + 2*unifRand())*R_DELTA; // variation +/- R_DELTA

	switch ( branch_type )
	{
	case 'T': // Terminal branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, theta_, 0.0, scaleFactor);
		break;
	case 't': // negative Terminal branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, -theta_, 0.0, scaleFactor);
		break;
	case 'L': // Lateral branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, beta_, delta_, scaleFactor);
		break;
	case 'l': // negative Lateral branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, -beta_, -delta_, scaleFactor);
		break;
	default:
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, theta_, 0.0, scaleFactor);
	}

	sub_.x = subX;
	sub_.y = subY;
	sub_.z = subZ;
	return sub_;
}

point generateBranch(point A, point B, char branch_type, double scaleFactor)
{
	point sub_;
	double subX;
	double subY;
	double subZ;

	double theta_ = THETA + THETA*(-1 + 2*unifRand())*R_THETA; // variation +/- R_THETA
	double beta_ = BETA + BETA*(-1 + 2*unifRand())*R_BETA; // variation +/- R_BETA
	double delta_ = DELTA + DELTA*(-1 + 2*unifRand())*R_DELTA; // variation +/- R_DELTA
	double scaleFactor_ = scaleFactor + scaleFactor*(-1 + 2*unifRand())*R_SF; // variation +/- R_SF

	switch ( branch_type )
	{
	case 'T': // Terminal branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, theta_, 0.0, scaleFactor_);
		break;
	case 't': // negative Terminal branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, -theta_, 0.0, scaleFactor_);
		break;
	case 'L': // Lateral branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, beta_, delta_, scaleFactor_);
		break;
	case 'l': // negative Lateral branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, -beta_, -delta_, scaleFactor_);
		break;
	default:
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, theta_, 0.0, scaleFactor_);
	}

	sub_.x = subX;
	sub_.y = subY;
	sub_.z = subZ;
	return sub_;
}

point generateBranch(point A, point B, node N)
{
	// get sth nr and branch type from node;
	point sub_;
	double subX;
	double subY;
	double subZ;

	double theta_ = THETA + THETA*(-1 + 2*unifRand())*R_THETA; // variation +/- R_THETA
	double beta_ = BETA + BETA*(-1 + 2*unifRand())*R_BETA; // variation +/- R_BETA
	double delta_ = DELTA + DELTA*(-1 + 2*unifRand())*R_DELTA; // variation +/- R_DELTA
	double scaleFactor_ = BR_SCALE_FACTOR + BR_SCALE_FACTOR*(-1 + 2*unifRand())*R_SF; // variation +/- R_SF

	// changed according to the formula given by Viennot. et.al (ACM 1989)
	//double branchLength = BUD_LENGTH*N.strahler_number*scaleFactor_;
	double c1 = 1.6;
	double branchLength = BUD_LENGTH*c1*N.strahler_number*scaleFactor_; // for linear growth (added BUD_LENGTH and scale factor by my own) . Make sure to change main.() if you doing a change here
	//double branchLength = BUD_LENGTH*c1*N.strahler_number*N.strahler_number*scaleFactor_; // for quadratic growth
	
	double u = B.x - A.x;
	double v = B.y - A.y;
	double w = B.z - A.z;

	double L1 = abs(sqrt(u*u + v*v + w*w));

	double scaleFactor = branchLength/ L1;


	switch ( N.branch_type )
	{
	case 'T': // Terminal branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, theta_, 0.0, scaleFactor);
		break;
	case 't': // negative Terminal branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, -theta_, 0.0, scaleFactor);
		break;
	case 'L': // Lateral branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, beta_, delta_, scaleFactor);
		break;
	case 'l': // negative Lateral branch
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, -beta_, -delta_, scaleFactor);
		break;
	default:
		calcSubBranchCoords(A.x, A.y, A.z, B.x, B.y, B.z, subX, subY, subZ, theta_, 0.0, scaleFactor);
	}

	sub_.x = subX;
	sub_.y = subY;
	sub_.z = subZ;
	return sub_;
}
