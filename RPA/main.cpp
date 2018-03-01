#include <iostream>
#include "Array_Sequences.h"
#include "Primer_Set.h"
#include "Optimization_Toolbox.h"
#include <time.h>
#include "PCR_Profile.h"
#include <fstream>
#include <cstdlib>
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
bool prepare_pareto(PCR_Profile ** pcr_profiles_4_comparison, unsigned int num_pcr_profiles, unsigned int & optimal_solution, ostream &err_msg=cout)
{
	double *genome_coverage, *genome_overrepresentation;
	bool *pareto_set;

	genome_coverage = new double[num_pcr_profiles];
	genome_overrepresentation = new double[num_pcr_profiles];
	pareto_set = new bool[num_pcr_profiles];

	for (int i = 0; i < num_pcr_profiles; i++)
	{
		genome_coverage[i] = pcr_profiles_4_comparison[i]->get_total_lenght_long_amplicons();
		genome_overrepresentation[i] = pcr_profiles_4_comparison[i]->get_total_lenght_short_amplicons();
		pareto_set[i] = false;
	}

	Optimization_Toolbox::calculate_pareto_frontier(genome_coverage, genome_overrepresentation, pareto_set, num_pcr_profiles, true, false, err_msg);
	
	int position = -1;
	
	int *tmp_pareto = new int[num_pcr_profiles];
	memset(tmp_pareto, 0, sizeof(int)*num_pcr_profiles);
	int n_pareto = 0;
	for (int i = 0; i < num_pcr_profiles; i++) if (pareto_set[i] == true) tmp_pareto[n_pareto++] = i;
	assert(n_pareto);
	position = tmp_pareto[rand() % n_pareto];
	delete[]tmp_pareto;


	optimal_solution = position;
	delete[] genome_coverage;
	delete[] genome_overrepresentation;
	delete[] pareto_set;
	return true;
}
int main()
{
	srand(time(0));

	Array_Sequences * as;
	//as = new Array_Sequences("sequence.fasta");
	as = new Array_Sequences("synthetic.fasta"); assert(as);

	as->show_Statistics();
	//as->show_All();
	Primer_Set * primers = new Primer_Set("primers.fasta"); assert(primers);
	unsigned int number_of_individual_primers = 2048;
	Primer_Set ** individual_primers;
	individual_primers = new Primer_Set *[number_of_individual_primers]; assert(individual_primers);
	PCR_Profile ** individual_PCR_profiles;
	individual_PCR_profiles = new PCR_Profile *[number_of_individual_primers]; assert(individual_PCR_profiles);
	ofstream log_out;
	log_out.open("log_file.txt");

	if (!log_out.is_open()) cout << "couldnt open log file" << endl;
	//pareralize this loop
	
	individual_primers[0] = new Primer_Set(1, 6); assert(individual_primers[0]);
	individual_primers[0]->add_primer(primers->get_primer_as_value(0));

	individual_PCR_profiles[0] = new PCR_Profile(individual_primers[0], as->get_pointer_to_sequence_object(0)); assert(individual_PCR_profiles[0]);
	unsigned int count = 0;
	//#pragma omp parallel for
	for (int i = 1; i < number_of_individual_primers; i++)
	{
		individual_primers[i] = new Primer_Set(1, 6); assert(individual_primers[i]);
		individual_primers[i]->add_primer(primers->get_primer_as_value(i));
		individual_PCR_profiles[i] = new PCR_Profile(individual_primers[i], as->get_pointer_to_sequence_object(0), individual_PCR_profiles[0]->get_pointer_to_pos_strand_sequence_int_profile());
		assert(individual_PCR_profiles[i]);
		
		//#pragma omp critical
		{
			count++;
			//if (count % 300 == 0)cout << "Processed " << count << " primers" << endl; 
		}
				
	}
	//for (int i = 0; i < number_of_individual_primers; i++)individual_PCR_profiles[i]->show_All();
	
	unsigned int pareto_index;
	prepare_pareto(individual_PCR_profiles, number_of_individual_primers, pareto_index);
	PCR_Profile * pareto_PCR_profile = new PCR_Profile(individual_PCR_profiles[pareto_index]); assert(pareto_PCR_profile);
	
	pareto_PCR_profile->show_All();

	while (true)
	{
		count = 0;
		PCR_Profile ** temp_pareto_PCR_profile;
		temp_pareto_PCR_profile = new PCR_Profile *[number_of_individual_primers]; assert(temp_pareto_PCR_profile);
		//#pragma omp parallel for
		for (int i = 0; i < number_of_individual_primers; i++)
		{
			temp_pareto_PCR_profile[i] = new PCR_Profile(pareto_PCR_profile, individual_PCR_profiles[i]); assert(temp_pareto_PCR_profile[i]);
			//for (int j = 0; j < temp_pareto_PCR_profile[i]->get_profile_length(); j++)
			//{
			//	if (temp_pareto_PCR_profile[i]->get_pointer_to_primer_locations()[j] != 0) cout << temp_pareto_PCR_profile[i]->get_pointer_to_primer_locations()[j];
			//}
			//cout << endl;
			
			//#pragma omp critical
			{
				count++;
				//if (count % 300 == 0)cout << "Processed " << count << " primers" << endl; 
			}


		}
		prepare_pareto(temp_pareto_PCR_profile, number_of_individual_primers, pareto_index);
		cout << "Pareto finished" << endl;
		delete pareto_PCR_profile;
		pareto_PCR_profile = new PCR_Profile(temp_pareto_PCR_profile[pareto_index]); assert(pareto_PCR_profile);
		//for (int i = 0; i < pareto_PCR_profile->get_profile_length(); i++)
		//{
		//	if (pareto_PCR_profile->get_pointer_to_primer_locations()[i] != 0) cout << pareto_PCR_profile->get_pointer_to_primer_locations()[i];
		//}
		//cout << endl;
		pareto_PCR_profile->show_All();

		for (int i = 0; i<number_of_individual_primers; i++) delete temp_pareto_PCR_profile[i];
		delete[]temp_pareto_PCR_profile;

	}


	log_out.close();
	system("PAUSE");
	return 1;
}