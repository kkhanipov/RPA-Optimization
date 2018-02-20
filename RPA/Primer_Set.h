#ifndef Primer_Set_H
#define Primer_Set_H



#include <iostream>
#include "Array_Sequences.h"


using namespace std;
//this object will hold a list of primers to be used in the computation. It will be part of the PCR profile object. 

class Primer_Set
{
	unsigned int primer_length; // the length of the primers

	unsigned int * primer; // contains a list of primers in their integer value
	unsigned int max_number_of_primers; // the array size of the list

	unsigned int number_of_primers; // the number of primers in the list

public:

	Primer_Set(unsigned int _max_number_of_primers, unsigned int _primer_length, ostream & err_msg = cout);
	Primer_Set(char * filename, unsigned int _max_number_of_primers, ostream & err_msg = cout); // ability to read a list or primers from a FASTA file
	Primer_Set(Primer_Set * primer_set, ostream & err_msg = cout);
	Primer_Set(unsigned int * primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout); //constructor to make a primer set from a array of primers converted to integers
	Primer_Set(char ** primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout);//constructor to make a primer set from a array of primers in char array

	bool write_to_file(char * filename, ostream & err_msg = cout); //writes the list of primers in FASTA format to file

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); //prints out the nucleotide contribution of the sequence

	bool show_All(ostream & out = cout, ostream &err_msg = cout); //prints out the show_statistics() and then the actual sequences

	unsigned int get_primer_as_value(unsigned int position) { return primer[position]; }

	//	bool get_primer_as_txt(unsigned int position, char *& primer, ostream & err_msg = cout); //copies the primers to the char ** primers object as text

	bool delete_primer(int position, ostream & err_msg = cout); //ability to mark a primer for deletion at a specific position in the array. 
//	bool delete_primer(int primer_value, ostream & err_msg = cout); //mark a primer for deletion based on its value
//	bool delete_primer(char * primer, ostream & err_msg = cout); //mark a primer for deletion based on its sequence

	bool add_primer(unsigned int primer_value, ostream & err_msg = cout); // add a primer to the array based on its integer value
	bool add_primer(char * primer, ostream & err_msg = cout); // add a primer to the array based on its sequence

	bool convert_primer_int_to_txt(unsigned int primer_value, char *& primer, ostream & err_msg = cout);
	bool convert_primer_txt_to_int(char * primer, unsigned int & primer_value, ostream & err_msg = cout);
	unsigned int convert_primer_to_reverse_complement(int position, ostream & err_msg = cout);
	unsigned int * get_pointer_to_primer_array() { return primer; };

	unsigned int get_number_of_primers() { return number_of_primers; };
	unsigned int get_primer_length() { return primer_length; };

};
Primer_Set::Primer_Set(unsigned int _max_number_of_primers, unsigned int _primer_length, ostream & err_msg)
{
	primer_length = _primer_length;
	number_of_primers = 0;
	max_number_of_primers = _max_number_of_primers;
	primer = new unsigned int[_max_number_of_primers];
}
Primer_Set::Primer_Set(Primer_Set * primer_set, ostream & err_msg)
{
	primer_length = primer_set->get_primer_length();
	number_of_primers = 0;
	max_number_of_primers = 4096;
	primer = new unsigned int[max_number_of_primers];
	for(int i=0;i<primer_set->get_number_of_primers();i++)add_primer(primer_set->get_primer_as_value(i));
}
bool Primer_Set::convert_primer_txt_to_int(char * _primer, unsigned int & primer_value, ostream & err_msg)
{
	primer_value = 0;
	char * _primer_value;
	_primer_value = new char[primer_length];
	for (int i = 0; i < 6; i++)
	{
		switch (_primer[i])
		{
		case 'A':_primer_value[i] = '1'; break;
		case 'T':_primer_value[i] = '2'; break;
		case 'C':_primer_value[i] = '3'; break;
		case 'G':_primer_value[i]=	'4'; break;
		default: strcpy(_primer_value, "555555"); primer_value = atoi(_primer_value);  return true;
		}
	}
	primer_value = atoi(_primer_value);
	return true;
}
Primer_Set::Primer_Set(char * filename, unsigned int _max_number_of_primers, ostream & err_msg)
{
	Array_Sequences * as = new Array_Sequences(filename, err_msg);
	number_of_primers = as->get_number_of_sequences();
	if (number_of_primers == 0)
	{
		err_msg << "Primer_Set::Primer_Set(...) 0 primers were read from the file";
		throw("could not read file with primers");
	}
	primer_length = as->get_pointer_to_sequence_object(0)->get_sequence_length();
	max_number_of_primers = _max_number_of_primers;
	primer = new unsigned int[_max_number_of_primers];
	unsigned int primer_val = 555555;
	for (int i = 0; i < number_of_primers; i++)
	{
		primer_val = 555555;
		convert_primer_txt_to_int(as->get_pointer_to_sequence_object(i)->get_pointer_to_sequence(), primer_val);
		primer[i] = primer_val;
	}
	as->~Array_Sequences();
}

bool Primer_Set::write_to_file(char * filename, ostream & err_msg)
{
	ofstream out;
	out.open(filename);
	if (!out.is_open())
	{
		err_msg << "ERROR: bool Primer_Set::write_to_file(), could not open " << filename << " file" << endl;
		return false;
	}
	if (number_of_primers == 0)
	{
		err_msg << "ERROR: bool Primer_Set::write_to_file() there are no primers to write" << endl;
		return false;
	}
	for (int i = 0; i < number_of_primers; i++)
	{
		out << ">" << i << endl;
		char * text_primer;
		if (!convert_primer_int_to_txt(primer[i], text_primer))
		{
			err_msg << "ERROR: bool Primer_Set::write_to_file() could not convert primer " << i << " to text" << endl;
		}
		out << text_primer << endl;
	}
	out.close();
	return true;
}

bool Primer_Set::show_statistics(ostream & out, ostream &err_msg)
{
	int count_A = 0;
	int count_C = 0;
	int count_T = 0;
	int count_G = 0;
	int count_N = 0;

	out << "Primer length is " << primer_length << endl;
	out << "Number of primers is " << number_of_primers << endl;
	for (int i = 0; i < number_of_primers; i++)
	{
		char * text = NULL;
		
		if (!convert_primer_int_to_txt(primer[i], text))
		{
			err_msg << "ERROR: bool Primer_Set::show_statistics() could not convert primer " << i << " to text" << endl;
		}
		for (int i = 0; i < primer_length; i++)
		{
			switch (text[i])
			{
			case 'A':case 'a': count_A++; break;
			case 'T':case 't': count_T++; break;
			case 'C':case 'c': count_C++; break;
			case 'G':case 'g': count_G++; break;
			case 'N':case 'n': count_N++; break;
			default:
				err_msg << "ERROR: Sequence::Sequence ==> unexpected character: sequence[" << i << "]=" << text[i] << endl;
				return false;
			}
		}
	}
	out << "Nucleotide Contribution is:" << endl;
	out << "A :" << count_A << endl;
	out << "C :" << count_C << endl;
	out << "T :" << count_T << endl;
	out << "G :" << count_G << endl;
	out << "N :" << count_N << endl;

	return true;
}

bool Primer_Set::convert_primer_int_to_txt(unsigned int primer_value, char *& _primer, ostream & err_msg)
{
	_primer = new char[primer_length+1];
	char * conversion_primer = new char[primer_length];
	_itoa((int)primer_value, conversion_primer, 10);
	for (int i = 0; i < 6; i++)
	{
		switch (conversion_primer[i])
		{
		case '1':_primer[i] = 'A'; break;
		case '2':_primer[i] = 'T'; break;
		case '3':_primer[i] = 'C'; break;
		case '4':_primer[i] = 'G'; break;
		default: _primer[i] = 'N';
		}
	}
	_primer[primer_length] = '\0';
	
	return true;
}
unsigned int Primer_Set::convert_primer_to_reverse_complement(int position, ostream & err_msg)
{
	char * conversion_primer = new char[primer_length];
	_itoa((int)primer[position], conversion_primer, 10);
	int ii = 0;
	for (int i = 5; i >= 0; i--)
	{
		
		switch (conversion_primer[i])
		{
		case '1':conversion_primer[ii] = 'T'; break;
		case '2':conversion_primer[ii] = 'A'; break;
		case '3':conversion_primer[ii] = 'G'; break;
		case '4':conversion_primer[ii] = 'C'; break;
		default: conversion_primer[ii] = 'N';
		}
		ii++;
	}
	return atoi(conversion_primer);
}
bool Primer_Set::show_All(ostream & out, ostream &err_msg)
{
	show_statistics();

	for (int i = 0; i < number_of_primers; i++)
	{
		char * text = NULL;
		if (!convert_primer_int_to_txt(primer[i], text))
		{
			err_msg << "ERROR: bool Primer_Set::show_statistics() could not convert primer " << i << " to text" << endl;
		}
		out << "Primer id: " << i << " Primer Value :" << primer[i] << " text: ";
		for (int j = 0; j < primer_length; j++)out << text[j];
		out << endl;
	}
	return true;
}
bool Primer_Set::add_primer(unsigned int primer_value, ostream & err_msg)
{
	if (primer_value == NULL)
	{
		err_msg << "ERROR: bool Primer_Set::add_primer() no primer_value==NULL" << endl;
		return false;
	}
	if (number_of_primers > max_number_of_primers - 1)
	{
		err_msg << "ERROR: bool Primer_Set::add_primer() number_of_primers > max_number_of_primers - 1" << endl;
		return false;
	}
	for (int i = 0; i < number_of_primers; i++)
	{
		if (primer[i] == primer_value)
		{
			err_msg << "ERROR: bool Primer_Set::add_primer() primer already included in the list" << endl;
			return true;
		}
	}
	primer[number_of_primers] = primer_value;
	number_of_primers++;
}
bool Primer_Set::add_primer(char * _primer, ostream & err_msg)
{
	if (_primer == NULL)
	{
		err_msg << "ERROR: bool Primer_Set::add_primer() no primer_value==NULL" << endl;
		return false;
	}
	if (number_of_primers > max_number_of_primers - 1)
	{
		err_msg << "ERROR: bool Primer_Set::add_primer() number_of_primers > max_number_of_primers - 1" << endl;
		return false;
	}
	unsigned int primer_value = 0;
	convert_primer_txt_to_int(_primer, primer_value);

	for (int i = 0; i < number_of_primers; i++)
	{
		if (primer[i] == primer_value)
		{
			err_msg << "ERROR: bool Primer_Set::add_primer() primer already included in the list" << endl;
			return true;
		}
	}
	primer[number_of_primers] = primer_value;
	number_of_primers++;
}




#endif // !Primer_Set_H