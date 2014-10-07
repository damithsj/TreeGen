#include <cmath>
#include <vector>
#include <ctime>

#include "mathUtil.h"

using namespace std;

//linspace function in MatLab. hacked the code ;)
// ...\toolbox\matlab\elmat\linspace.m
//generates N points between min and max.
vector<double> linspace(double min, double max, int n)
{
	vector<double> result;
	// vector iterator
	int iterator = 0;

	for (int i = 0; i <= n-2; i++)	
	{
		double temp = min + i*(max-min)/(floor((double)n) - 1);
		result.insert(result.begin() + iterator, temp);
		iterator += 1;
	}

	//iterator += 1;

	result.insert(result.begin() + iterator, max);
	return result;
}

// return a uniform number in [0,1].
double unifRand()
{
	//seed(); // seed is not called since we dont want to intialize generator each time
    return rand() / double(RAND_MAX);
}
// Reset the random number generator with the system clock.
void seed()
{
    srand(time(0));
}

bool evalProbability(double p)
{
	if (unifRand() < p)
		return true;
	else
		return false;

}