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

	Primer_Set(unsigned int _max_number_of_primers, ostream & err_msg = cout);
	Primer_Set(char * filename, unsigned int _max_number_of_primers,ostream & err_msg = cout); // ability to read a list or primers from a FASTA file
	Primer_Set(Primer_Set * primer_set, ostream & err_msg = cout);
	Primer_Set(unsigned int * primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout); //constructor to make a primer set from a array of primers converted to integers
	Primer_Set(char ** primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout);//constructor to make a primer set from a array of primers in char array

	bool write_to_file(char * filename, ostream & err_msg = cout); //writes the list of primers in FASTA format to file

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); //prints out the nucleotide contribution of the sequence
	 
	bool show_All(ostream & out = cout, ostream &err_msg = cout); //prints out the show_statistics() and then the actual sequences

	unsigned int get_primer_as_value(unsigned int position, ostream & err_msg = cout); 

	bool get_primer_as_txt(char * primers, ostream & err_msg = cout); //copies the array of primers to the char ** primers object as text

	bool delete_primer(int position, ostream & err_msg = cout); //ability to mark a primer for deletion at a specific position in the array. 
//	bool delete_primer(int primer_value, ostream & err_msg = cout); //mark a primer for deletion based on its value
//	bool delete_primer(char * primer, ostream & err_msg = cout); //mark a primer for deletion based on its sequence

	bool add_primer(unsigned int primer_value, ostream & err_msg = cout); // add a primer to the array based on its integer value
	bool add_primer(char * primer, ostream & err_msg = cout); // add a primer to the array based on its sequence

	bool convert_primer_int_to_txt(unsigned int primer_value, char * primer, ostream & err_msg = cout);
	bool convert_primer_txt_to_int(char * primer, unsigned int & primer_value,  ostream & err_msg = cout);

	unsigned int * get_pointer_to_primer_array() { return primer };

	unsigned int get_number_of_primers() { return number_of_primers; };
	unsigned int get_primer_length() { return primer_length };

};
Primer_Set::Primer_Set(unsigned int _max_number_of_primers, ostream & err_msg)
{
	primer_length = 0;
	number_of_primers = 0;
	max_number_of_primers = _max_number_of_primers;
	primer = new unsigned int[_max_number_of_primers];
}

Primer_Set::Primer_Set(char * filename, unsigned int _max_number_of_primers, ostream & err_msg)
{
	Array_Sequences * as = new Array_Sequences(filename, err_msg);
	number_of_primers=as->get_number_of_sequences();
	if (number_of_primers == 0)
	{
		err_msg << "Primer_Set::Primer_Set(...) 0 primers were read from the file";
		throw("could not read file with primers");
	}
	primer_length = as->get_pointer_to_sequence_object(0)->get_sequence_length();
	max_number_of_primers = _max_number_of_primers;
	primer = new unsigned int[_max_number_of_primers];
	unsigned int primer_val = 999999999;
	for (int i = 0; i < number_of_primers; i++)
	{
		primer_val = 999999999;
		convert_primer_txt_to_int(as->get_pointer_to_sequence_object(i)->get_pointer_to_sequence(), primer_val);
		primer[i] = primer_val;
	}
	as->~Array_Sequences();
}