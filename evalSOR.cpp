#include <iostream>
#include <cmath>
#include "evalSOR.h"
#include "initializeParams.h"
#include "definitions.h"

using namespace std;
using namespace constants;

// This function will evaluate the potential of each cell in the grid using Laplace equation
// Successive over-relaxation (SOR) method was used to calculate the solution
// More Info regarding Laplace and SOR
// http://www.damtp.cam.ac.uk/lab/people/sd/lectures/nummeth98/pdes.htm

//grid - actual 3D grid which used for the DBM
//source_distrib - grid which has same dims as grid and has vaule 1 for cells which occupied (by source, obstacle or branch)
void evalSOR (dVec3D (&grid), dVec3D source_distrib)
{
	bool terminate = false;
	double residual;

	while (terminate == false)
	{
		terminate = true;
		for (int i = 1; i < HORZ_LENGTH - 1; i++)
		{
			for (int j = 1; j < VERT_LENGTH - 1; j++)
			{
				for (int k = 1; k < DEPTH_LENGTH - 1; k++)
				{
					if (source_distrib[i][j][k] != 1) // check whether the current cell is occupied
					{
						// here we go baby, relaxxxx...
						residual = grid[i][j][k-1] + grid[i][j][k+1] + grid[i][j-1][k] + grid[i][j+1][k] + grid[i-1][j][k] + grid[i+1][j][k] - 6*grid[i][j][k];
						grid[i][j][k] = grid[i][j][k] + ORP/6*residual;

						if (terminate == true)
						{
							if (abs(residual) > CUTOFF)
							{
								//cout << residual << endl;
								terminate = false;
							}
						}
					}
				}
			}
		}
	}
}