#ifndef PCR_Profile_H
#define PCR_Profile_H


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
	unsigned int * pos_strand_sequence_int_profile;
	//unsigned int * neg_strand_sequence_int_profile;

	Sequence * sequence;
	//sequence to use for the primer_locations_profile
	

	unsigned int number_forward_primers; // calculated by calculate_statistics(), contains the number of forward primers in the primer locations profile
	unsigned int number_reverse_primers;// calculated by calculate_statistics(), contains the number of reverse primers in the primer locations profile
	unsigned int number_short_amplicons; // calculated by calculate_statistics(), contains the number of amplicons shorter than designed by a global variable in the primer locations profile
	unsigned int number_long_amplicons;// calculated by calculate_statistics(), contains the number of amplicons of a good size than designed by a global variable in the primer locations profile
	unsigned int total_lenght_short_amplicons; // calculated by calculate_statistics(),  total length of the primer locations profile covered by short amplicons
	unsigned int total_lenght_long_amplicons;// calculated by calculate_statistics(),  total length of the primer locations profile covered by long amplicons

	unsigned int min_primer_distance;
	unsigned int max_primer_distance;

public:
	PCR_Profile(PCR_Profile * _pcr_profile_a, PCR_Profile * _pcr_profile_b, ostream &err_msg = cout);
	PCR_Profile(PCR_Profile * _pcr_profile, ostream &err_msg = cout);
	PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg = cout);

	bool PCR_profile_calculation(ostream &err_msg = cout);
	bool calculate_statistics(ostream & out = cout, ostream &err_msg = cout); //calculates statistics on the primer_locations array

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); // gives the values of the calculate_statistics();
	bool show_All(ostream & out = cout, ostream &err_msg = cout); // first shows the show_statistics() and then prints the primer_location profile

	unsigned int get_number_forward_primers() { return number_forward_primers;}
	unsigned int get_number_reverse_primers() { return number_reverse_primers; }
	unsigned int get_number_short_amplicons() { return number_short_amplicons; }
	unsigned int get_number_long_amplicons() { return number_long_amplicons; }
	unsigned int get_total_lenght_short_amplicons() { return total_lenght_short_amplicons; }
	unsigned int get_total_lenght_long_amplicons() { return total_lenght_long_amplicons; }
	unsigned int get_profile_length() { return profile_length; }
	Primer_Set * get_pointer_to_primer_set() { return p_set; }
	int * get_pointer_to_primer_locations() { return primer_locations; }
};
PCR_Profile::PCR_Profile(PCR_Profile * _pcr_profile_a, PCR_Profile * _pcr_profile_b, ostream &err_msg)
{
	number_forward_primers = 0;
	number_reverse_primers = 0;
	number_short_amplicons = 0;
	number_long_amplicons = 0;
	total_lenght_short_amplicons = 0;
	total_lenght_long_amplicons = 0;
	min_primer_distance = 50;
	max_primer_distance = 1500;
	if (_pcr_profile_a->get_profile_length() != _pcr_profile_b->get_profile_length())
	{
		err_msg << "ERROR: PCR_Profile::PCR_Profile(...) Profile Lenght of Profiles A and B is not the same" << endl;
		throw("PCR Profiles must have same profile lenght");
	}
	profile_length = _pcr_profile_a->get_profile_length();
	unsigned int max_number_of_primers = _pcr_profile_a->get_pointer_to_primer_set()->get_number_of_primers() + _pcr_profile_b->get_pointer_to_primer_set()->get_number_of_primers();
	p_set = new Primer_Set(max_number_of_primers);

	for (int i = 0; i < _pcr_profile_a->get_pointer_to_primer_set()->get_number_of_primers(); i++)
	{
		p_set->add_primer(_pcr_profile_a->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	for (int i = 0; i < _pcr_profile_b->get_pointer_to_primer_set()->get_number_of_primers(); i++)
	{
		p_set->add_primer(_pcr_profile_b->get_pointer_to_primer_set()->get_primer_as_value(i));
	}
	primer_locations = new int[profile_length];
	for (int i = 0; i < profile_length; i++)
	{
		primer_locations[i] = _pcr_profile_a->get_pointer_to_primer_locations()[i];
	}
	for (int i = 0; i < profile_length; i++)
	{
		if (primer_locations[i] == -5) continue;
		if (primer_locations[i] == 5) continue;
		if (primer_locations[i] == 1 && _pcr_profile_b->get_pointer_to_primer_locations()[i] == 1) continue;
		if (primer_locations[i] == -1 && _pcr_profile_b->get_pointer_to_primer_locations()[i] == -1) continue;
		if (primer_locations[i] + _pcr_profile_b->get_pointer_to_primer_locations()[i] == 0)primer_locations[i] = 5;
		else _pcr_profile_b->get_pointer_to_primer_locations()[i];
	}


}

PCR_Profile::PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg)
{
	number_forward_primers = 0;
	number_reverse_primers = 0;
	number_short_amplicons = 0;
	number_long_amplicons = 0;
	total_lenght_short_amplicons = 0;
	total_lenght_long_amplicons = 0;
	min_primer_distance = 50;
	max_primer_distance = 1500;

	p_set = new Primer_Set(_primer_set, err_msg);
	sequence = new Sequence(_sequence, err_msg);

	profile_length = sequence->get_sequence_length()-p_set->get_primer_length();
	primer_locations = new int[profile_length];
	pos_strand_sequence_int_profile = new unsigned int [profile_length];
	
}

bool PCR_Profile::show_statistics(ostream & out, ostream &err_msg)
{
	out << "number_forward_primers " << get_number_forward_primers() << endl;
	out << "number_reverse_primers " << get_number_reverse_primers() << endl;
	out << "number_short_amplicons " << get_number_short_amplicons() << endl;
	out << "number_long_amplicons " << get_number_long_amplicons() << endl;
	out << "total_lenght_short_amplicons " << get_total_lenght_short_amplicons() << endl;
	out << "total_lenght_long_amplicons " << get_total_lenght_long_amplicons() << endl;
}
bool PCR_Profile::show_All(ostream & out, ostream &err_msg)
{
	show_statistics();
}

bool PCR_Profile::PCR_profile_calculation(ostream &err_msg)
{
	for (int j = 0; j < profile_length; j++)
	{
		primer_locations[j] = 0;
	}

	for (int i = 0; i < p_set->get_number_of_primers; i++)
	{
		for (int j = 0; j < profile_length; j++)
		{
			if (pos_strand_sequence_int_profile[j] == 999999999)
			{
				primer_locations[j] = -5;
				continue;
			}

			if (pos_strand_sequence_int_profile[j] == p_set->get_pointer_to_primer_array[i])
			{
				if (primer_locations[j] == -1 || primer_locations[j] == 5) primer_locations[j] = 5;
				else primer_locations[j] = 1;
			}
			if (pos_strand_sequence_int_profile[j] == p_set->convert_primer_to_reverse_complement(i))
			{
				if (primer_locations[j] == 1 || primer_locations[j] == 5) primer_locations[j] = 5;
				else primer_locations[j] = -1;
			}

		}

	}
}
bool PCR_Profile::calculate_statistics(ostream & out, ostream &err_msg)
{
	unsigned int start_position = 0;
	unsigned int end_position = 0;
	for (int j = 0; j < profile_length; j++)
	{
		if (start_position == 0)
		{
			for (int i = j; i < profile_length; i++)
			{
				if (primer_locations[i] == 1 || primer_locations[i] == 5)
				{
					start_position = i;
				}
			}
		}
		for (int i = start_position + 1; i < profile_length; i++)
		{
			if (primer_locations[i] == -1 || primer_locations[i] == 5)
			{
				end_position = i;
				unsigned int distance = end_position - start_position;
				if (distance >= min_primer_distance && distance <= max_primer_distance)
				{
					total_lenght_long_amplicons += distance;
					number_long_amplicons++;
				}
				if (distance < min_primer_distance)
				{
					total_lenght_short_amplicons += distance;
					number_short_amplicons++;
				}
				if(primer_locations[i] == 5)start_position = end_position;
				else start_position = 0;
				j = end_position;
			}
		}
	}

}


#endif