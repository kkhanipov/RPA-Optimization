#include "Primer_Set.h"
#include "Sequence.h"
#include <iostream>

using namespace std;
// Object PCR profile will be used to store the primer sets and information about how they are mapping to the sequence of interest

class PCR_Profile
{
	Primer_Set * p_set; // this will be the primer set used to create the primer locations profile

	unsigned int profile_length; //length of the primer_locations array, should be the size of the sequence of interest
	int * primer_locations; 
	// an array the size of the seqeunce of interest to denote where the primers from the p_set are mapping. +1 will denote forward mapping primer location, 
	//0 will denote a non mapped location, -1 will denote a reverse mapping primer location, -5 will denote positions which are not available due to unknown nucletodies 
	//being present there, 5 will denote if there are primers in both directions at the same location.

	Sequence * sequence;
	//sequence to use for the primer_locations_profile
	

	unsigned int number_forward_primers; // calculated by calculate_statistics(), contains the number of forward primers in the primer locations profile
	unsigned int number_reverse_primers;// calculated by calculate_statistics(), contains the number of reverse primers in the primer locations profile
	unsigned int number_short_amplicons; // calculated by calculate_statistics(), contains the number of amplicons shorter than designed by a global variable in the primer locations profile
	unsigned int number_long_amplicons;// calculated by calculate_statistics(), contains the number of amplicons of a good size than designed by a global variable in the primer locations profile
	unsigned int total_lenght_short_amplicons; // calculated by calculate_statistics(),  total length of the primer locations profile covered by short amplicons
	unsigned int total_lenght_long_amplicons;// calculated by calculate_statistics(),  total length of the primer locations profile covered by long amplicons

public:
	PCR_Profile(PCR_Profile * _pcr_profile, ostream &err_msg = cout);
	PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg = cout);


	bool calculate_statistics(ostream & out = cout, ostream &err_msg = cout); //calculates statistics on the primer_locations array

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); // gives the values of the calculate_statistics();
	bool show_All(ostream & out = cout, ostream &err_msg = cout); // first shows the show_statistics() and then prints the primer_location profile

	unsigned int get_number_forward_primers() { return number_forward_primers;}
	unsigned int get_number_reverse_primers() { return number_reverse_primers; }
	unsigned int get_number_short_amplicons() { return number_short_amplicons; }
	unsigned int get_number_long_amplicons() { return number_long_amplicons; }
	unsigned int get_total_lenght_short_amplicons() { return total_lenght_short_amplicons; }
	unsigned int get_total_lenght_long_amplicons() { return total_lenght_long_amplicons; }

};

PCR_Profile::PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg)
{
	number_forward_primers = 0;
	number_reverse_primers = 0;
	number_short_amplicons = 0;
	number_long_amplicons = 0;
	total_lenght_short_amplicons = 0;
	total_lenght_long_amplicons = 0;
	p_set = new Primer_Set(_primer_set, err_msg);
	sequence = new Sequence(_sequence, err_msg);

	profile_length = sequence->get_sequence_length()-p_set->get_primer_length();
	primer_locations = new int[profile_length];
}