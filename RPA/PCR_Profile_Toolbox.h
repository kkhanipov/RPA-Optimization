#ifndef PCR_Profile_Toolbox_H
#define PCR_Profile_Toolbox_H

#include "PCR_Profile.h"
#include <iostream>

using namespace std;
class PCR_Profile_Toolbox
{
public:
	static bool merge_pcr_profiles(PCR_Profile * & resulting, PCR_Profile * input_a, PCR_Profile * input_b, ostream &err_msg = cout);
	//takes two PCR profiles and merges them into a new third one. As part of the process, it recalculates the statistics in the new PCR profile.
	static bool convert_primer_txt_to_int(char * _primer, unsigned int primer_length, unsigned int & primer_value, ostream & err_msg=cout);
	
};

bool PCR_Profile_Toolbox::merge_pcr_profiles(PCR_Profile * & resulting, PCR_Profile * input_a, PCR_Profile * input_b, ostream &err_msg)
{
	return true;
}

bool PCR_Profile_Toolbox::convert_primer_txt_to_int(char * _primer, unsigned int primer_length, unsigned int & primer_value, ostream & err_msg)
{
	primer_value = 0;
	char * _primer_value = new char[primer_length]; assert(_primer_value);
	for (int i = 0; i < primer_length; i++)
	{
		switch (_primer[i])
		{
		case 'A':_primer_value[i] = '1'; break;
		case 'T':_primer_value[i] = '2'; break;
		case 'C':_primer_value[i] = '3'; break;
		case 'G':_primer_value[i] = '4'; break;
		default: strcpy(_primer_value, "555555"); primer_value = atoi(_primer_value);  return true;
		}
	}
	primer_value = atoi(_primer_value);
	return true;
}


#endif