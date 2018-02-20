#include <iostream>
#include "Array_Sequences.h"
#include "Primer_Set.h"
#include "Optimization_Toolbox.h"
#include <time.h>
#include "PCR_Profile.h"
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
	Primer_Set * primers = new Primer_Set("primers.fasta", 2048);

	Primer_Set ** individual_primers;
	individual_primers = new Primer_Set *[2048];
	PCR_Profile ** individual_PCR_profiles;
	individual_PCR_profiles = new PCR_Profile *[2048];
	for (int i = 0; i < 2048; i++)
	{
		individual_primers[i] = new Primer_Set(1,6);
		individual_primers[i]->add_primer(primers->get_primer_as_value(i));
		individual_PCR_profiles[i] = new PCR_Profile(individual_primers[i], as->get_pointer_to_sequence_object(0));
		fstream out;
		out.open("dump.txt");
		individual_PCR_profiles[i]->show_statistics(out);
		out.close();
	}

	
	
	system("PAUSE");
	return 1;
}