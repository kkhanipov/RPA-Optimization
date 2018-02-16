#include <iostream>
#include "Array_Sequences.h"
#include "Primer_Set.h"
#include "Optimization_Toolbox.h"
#include <time.h>
using namespace std;

void test_1()
{
	double * x = new double[100];
	double * y = new double[100];
	bool * pareto = new bool[100];
	srand(time(0));
	for (int i = 0; i < 100; i++)
	{
		x[i] = (double)rand();
		y[i] = (double)rand();
		pareto[i] = 0;
	}
	
	Optimization_Toolbox::calculate_pareto_frontier(x, y, pareto, 100, false, true);
}

int main()
{

	Array_Sequences * as;
	as = new Array_Sequences("test.fasta");
	as->show_Statistics();
	as->show_All();

	test_1();
	system("PAUSE");
	return 1;
}