#include <iostream>

using namespace std;

class Optimization_Toolbox
{
public:
	static bool calculate_parretto_frontier(double * x, double * y, bool * paretto_set, unsigned int number_of_values, bool maximize_x, bool maximize_y, ostream &err_msg = cout);
	
	//the function will calculate the parretto frontier of a dataset. It will take an array of x and y coordinates. Maximize_x and Maximize_y are to denote whether the values need to maximized or minimized
	//Use quicksort on them and pick the best point.
	//From there it will go left and right to determine all of the points which belong in the parretto frontier. The function will return the paretto_set marking positions in which points are part of the 
	//paretto frontier as "true".
};

bool Optimization_Toolbox::calculate_parretto_frontier(double * x, double * y, bool * paretto_set, unsigned int number_of_values, bool maximize_x, bool maximize_y, ostream &err_msg)
{
	return true;
}