#ifndef PCR_Profile_H
#define PCR_Profile_H


#include "Primer_Set.h"
#include "Sequence.h"
#include <iostream>

using namespace std;
// Object PCR profile will be used to store the primer sets and information about how they are mapping to the sequence of interest

class PCR_Profile
{
	static unsigned int min_primer_distance;
	static unsigned int max_primer_distance;
	Primer_Set * p_set; // this will be the primer set used to create the primer locations profile

	unsigned int profile_length; //length of the primer_locations array, should be the size of the sequence of interest
	int * primer_locations; 
	// an array the size of the seqeunce of interest to denote where the primers from the p_set are mapping. +1 will denote forward mapping primer location, 
	//0 will denote a non mapped location, -1 will denote a reverse mapping primer location, -5 will denote positions which are not available due to unknown nucletodies 
	//being present there, 5 will denote if there are primers in both directions at the same location.
	
	unsigned int * pos_strand_sequence_int_profile;

	Sequence * sequence;
	//sequence to use for the primer_locations_profile
	struct Stats
	{
		unsigned int number_forward_primers; // calculated by calculate_statistics(), contains the number of forward primers in the primer locations profile
		unsigned int number_reverse_primers;// calculated by calculate_statistics(), contains the number of reverse primers in the primer locations profile
		unsigned int number_short_amplicons; // calculated by calculate_statistics(), contains the number of amplicons shorter than designed by a global variable in the primer locations profile
		unsigned int number_long_amplicons;// calculated by calculate_statistics(), contains the number of amplicons of a good size than designed by a global variable in the primer locations profile
		unsigned int total_lenght_short_amplicons; // calculated by calculate_statistics(),  total length of the primer locations profile covered by short amplicons
		unsigned int total_lenght_long_amplicons;// calculated by calculate_statistics(),  total length of the primer locations profile covered by long amplicons
		unsigned int total_lenght_too_long_amplicons;
	} stats;

public:
	PCR_Profile(PCR_Profile * _pcr_profile_a, PCR_Profile * _pcr_profile_b, ostream &err_msg = cout);
	PCR_Profile(PCR_Profile * _pcr_profile, ostream &err_msg = cout);
	PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg = cout);
	PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, unsigned int * _pos_strand_sequence_int_profile, ostream &err_msg = cout);
	~PCR_Profile();
	bool PCR_profile_calculation(ostream &err_msg = cout);
	bool calculate_statistics(ostream & out = cout, ostream &err_msg = cout); //calculates statistics on the primer_locations array

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); // gives the values of the calculate_statistics();
	bool show_All(ostream & out = cout, ostream &err_msg = cout); // first shows the show_statistics() and then prints the primer_location profile

	unsigned int get_number_forward_primers() { return stats.number_forward_primers;}
	unsigned int get_number_reverse_primers() { return stats.number_reverse_primers; }
	unsigned int get_number_short_amplicons() { return stats.number_short_amplicons; }
	unsigned int get_number_long_amplicons() { return stats.number_long_amplicons; }
	unsigned int get_total_lenght_short_amplicons() { return stats.total_lenght_short_amplicons; }
	unsigned int get_total_lenght_long_amplicons() { return stats.total_lenght_long_amplicons; }
	unsigned int get_profile_length() { return profile_length; }
	unsigned int get_total_lenght_too_long_amplicons() { return stats.total_lenght_too_long_amplicons; }
	Primer_Set * get_pointer_to_primer_set() { return p_set; }
	unsigned int get_number_of_primers() { return p_set->get_number_of_primers(); }
	int * get_pointer_to_primer_locations() { return primer_locations; }
	unsigned int * get_pointer_to_pos_strand_sequence_int_profile() { return pos_strand_sequence_int_profile; }

	Stats get_Stats() { return stats; }
};

unsigned int PCR_Profile::min_primer_distance = 100;
unsigned int PCR_Profile::max_primer_distance = 500;

PCR_Profile::~PCR_Profile()
{
	if (sequence != NULL) delete sequence;
	
	if (p_set != NULL) delete p_set;

	if (pos_strand_sequence_int_profile != NULL) delete[] pos_strand_sequence_int_profile;

	if (primer_locations != NULL) delete[] primer_locations;


}
PCR_Profile::PCR_Profile(PCR_Profile * _pcr_profile_a, PCR_Profile * _pcr_profile_b, ostream &err_msg)
{
	stats = {};
	
	if (_pcr_profile_a->get_profile_length() != _pcr_profile_b->get_profile_length())
	{
		err_msg << "ERROR: PCR_Profile::PCR_Profile(...) Profile Lenght of Profiles A and B is not the same" << endl;
		err_msg << "PCR Profiles must have same profile lenght" << endl;
		assert(NULL);
	}
	profile_length = _pcr_profile_a->get_profile_length();
	unsigned int max_number_of_primers = _pcr_profile_a->get_number_of_primers() + _pcr_profile_b->get_number_of_primers();
	sequence = NULL;
	p_set = new Primer_Set(max_number_of_primers, _pcr_profile_a->get_pointer_to_primer_set()->get_primer_length());
	assert(p_set);
	pos_strand_sequence_int_profile = NULL;

	for (int i = 0; i < _pcr_profile_a->get_number_of_primers(); i++)
	{
		p_set->add_primer(_pcr_profile_a->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	for (int i = 0; i < _pcr_profile_b->get_number_of_primers(); i++)
	{
		p_set->add_primer(_pcr_profile_b->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	primer_locations = new int[profile_length]; assert(primer_locations);

	for (int i = 0; i < profile_length; i++)
	{
		primer_locations[i] = _pcr_profile_a->get_pointer_to_primer_locations()[i];
		//if (primer_locations[i] != 0)cout << i << "\t" << primer_locations[i]<<"\t";
	}
	//cout << endl;

	for (int i = 0; i < profile_length; i++)
	{
		if (primer_locations[i] == -5) continue;
		if (primer_locations[i] == 5) continue;
		if (primer_locations[i] == 1 && _pcr_profile_b->get_pointer_to_primer_locations()[i] == 1) continue;
		if (primer_locations[i] == -1 && _pcr_profile_b->get_pointer_to_primer_locations()[i] == -1) continue;
		if (primer_locations[i] + _pcr_profile_b->get_pointer_to_primer_locations()[i] == 0 && abs(primer_locations[i])==1 )primer_locations[i] = 5;
		if (primer_locations[i]==0) primer_locations[i]=_pcr_profile_b->get_pointer_to_primer_locations()[i];
	}
	//for (int i = 0; i < profile_length; i++)
	//{
	//	if (primer_locations[i] != 0)cout << i << "\t" << primer_locations[i] << "\t";
	//}
	//cout << endl;
	calculate_statistics();

}

PCR_Profile::PCR_Profile(PCR_Profile * _pcr_profile, ostream &err_msg)
{
	stats = _pcr_profile->get_Stats();

	profile_length = _pcr_profile->get_profile_length();
	unsigned int max_number_of_primers = _pcr_profile->get_number_of_primers();
	sequence = NULL;
	pos_strand_sequence_int_profile = NULL;
	p_set = new Primer_Set(max_number_of_primers, _pcr_profile->get_pointer_to_primer_set()->get_primer_length());

	for (int i = 0; i < _pcr_profile->get_number_of_primers(); i++)
	{
		p_set->add_primer(_pcr_profile->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	primer_locations = new int[profile_length]; assert(primer_locations);
	for (int i = 0; i < profile_length; i++)
	{
		primer_locations[i] = _pcr_profile->get_pointer_to_primer_locations()[i];
	}
}
PCR_Profile::PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, unsigned int * _pos_strand_sequence_int_profile, ostream &err_msg)
{
	stats = {};

	sequence = NULL;
	p_set = new Primer_Set(_primer_set, err_msg); assert(p_set);
	
	profile_length = _sequence->get_sequence_length() - p_set->get_primer_length();
	primer_locations = new int[profile_length]; assert(primer_locations);
	pos_strand_sequence_int_profile = _pos_strand_sequence_int_profile;
	PCR_profile_calculation();
	pos_strand_sequence_int_profile = NULL;
	calculate_statistics();
}
PCR_Profile::PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg)
{
	stats = {};

	sequence = NULL;
	p_set = new Primer_Set(_primer_set, err_msg); assert(p_set);
	
	assert(_sequence->get_sequence_length() >= p_set->get_primer_length());
	profile_length = _sequence->get_sequence_length()-p_set->get_primer_length();
	primer_locations = new int[profile_length]; assert(primer_locations);
	pos_strand_sequence_int_profile = new unsigned int[profile_length]; assert(pos_strand_sequence_int_profile);
	unsigned int primer_val = 0;
	for (int i = 0; i < profile_length; i++)
	{
		primer_val = 0;
		p_set->convert_primer_txt_to_int(&_sequence->get_pointer_to_sequence()[i], primer_val);
		pos_strand_sequence_int_profile[i] = primer_val;
	}

	PCR_profile_calculation();
	
	calculate_statistics();
}

bool PCR_Profile::show_statistics(ostream & out, ostream &err_msg)
{
	out << "number_forward_primers " << get_number_forward_primers() << endl;
	out << "number_reverse_primers " << get_number_reverse_primers() << endl;
	out << "number_short_amplicons " << get_number_short_amplicons() << endl;
	out << "number_long_amplicons " << get_number_long_amplicons() << endl;
	out << "total_lenght_short_amplicons " << get_total_lenght_short_amplicons() << endl;
	out << "total_lenght_long_amplicons " << get_total_lenght_long_amplicons() << endl;
	out << "total_lenght_too_long_amplicons " << get_total_lenght_too_long_amplicons() << endl;
	return true;
}
bool PCR_Profile::show_All(ostream & out, ostream &err_msg)
{
	out << "Number of Primers " << p_set->get_number_of_primers() << endl;
	for (int i = 0; i < p_set->get_number_of_primers(); i++)
	{
		out << "Primers " << p_set->get_pointer_to_primer_array()[i] << endl;
	}
	show_statistics(out, err_msg);
	return true;
}

bool PCR_Profile::PCR_profile_calculation(ostream &err_msg)
{
	for (int j = 0; j < profile_length; j++)
	{
		primer_locations[j] = 0;
	}

	for (int i = 0; i < p_set->get_number_of_primers(); i++)
	{
		//err_msg << "num of primers" << p_set->get_number_of_primers() << endl;
		//err_msg << "Primer " << p_set->get_pointer_to_primer_array()[i]<<endl;
		//err_msg << "Reverse complement " << p_set->get_pointer_to_reverse_complement_primer_array()[i] << endl;


		for (int j = 0; j < profile_length; j++)
		{
			if (pos_strand_sequence_int_profile[j] == 555555)
			{
				primer_locations[j] = -5;
				continue;
			}

			if (pos_strand_sequence_int_profile[j] == p_set->get_pointer_to_primer_array()[i])
			{
				if (primer_locations[j] == -1 || primer_locations[j] == 5) primer_locations[j] = 5;
				else primer_locations[j] = 1;
			}
			if (pos_strand_sequence_int_profile[j] == p_set->get_pointer_to_reverse_complement_primer_array()[i])
			{
				if (primer_locations[j] == 1 || primer_locations[j] == 5) primer_locations[j] = 5;
				else primer_locations[j] = -1;
			}

		}
	}
	return true;
}

//#define DEBUG
bool PCR_Profile::calculate_statistics(ostream & out, ostream &err_msg)
{
#ifdef DEBUG
	char msg[64];
	static int call_count = 0;
	sprintf(msg,"island_%d.txt", call_count++);
	ofstream out_f(msg);
	out_f << "id\tbegin\tend\tdist" << endl;
#endif
	int start_position = 0;
	int end_position = 0;
	bool start_found;
	for (int j = 0; j < profile_length; j++)
	{
		start_found = false;
		if (start_position == 0)
		{
			for (int i = j; i < profile_length; i++)
			{
				if (primer_locations[i] == 1 || primer_locations[i] == 5)
				{
					start_found = true;
					start_position = i;
					stats.number_forward_primers++;
					break;
				}
			}
			if (!start_found) return true;
		}
		for (int i = start_position + 1; i < profile_length; i++)
		{
			if (primer_locations[i] == -1 || primer_locations[i] == 5)
			{
				end_position = i;
				stats.number_reverse_primers++;
				unsigned int distance = end_position - start_position+p_set->get_primer_length();
#ifdef DEBUG
				out_f << number_reverse_primers << '\t' << start_position << '\t' << end_position << '\t' << distance << endl;
#endif
			//	cout << "Distance " << distance << endl;
				if (distance >= min_primer_distance && distance <= max_primer_distance)
				{
					stats.total_lenght_long_amplicons += distance;
					stats.number_long_amplicons++;
				}
				if (distance < min_primer_distance)
				{
					stats.total_lenght_short_amplicons += distance;
					stats.number_short_amplicons++;
				}
				if(distance>max_primer_distance) stats.total_lenght_too_long_amplicons += distance;

				if (primer_locations[i] == 5)
				{
					stats.number_forward_primers++;
					start_position = end_position;
				}
				else start_position = 0;
				break;
			}
			j = i;
		}
	}

#ifdef DEBUG
	out_f.close();
#endif

	return true;
}
#ifdef DEBUG
#undef DEBUG
#endif

#endif