#pragma once
//this object will hold a list of primers to be used in the computation. It will be part of the PCR profile object. 

class Primer_Set
{
	unsigned int primer_length; // the length of the primers

	unsigned int * primer; // contains a list of primers in their integer value
	unsigned int max_number_of_primers; // the array size of the list

	unsigned int number_of_primers; // the number of primers in the list
	
public:

	Primer_Set(unsigned int _max_number_of_primers; ostream & err_msg = cout);
	Primer_Set(char * filename, ostream & err_msg = cout); // ability to read a list or primers from a FASTA file
	Primer_Set(Primer_Set * primer_set, ostream & err_msg = cout);
	Primer_Set(unsigned int * primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout); //constructor to make a primer set from a array of primers converted to integers
	Primer_Set(char ** primers, unsigned int number_primers, ostream & err_msg = cout);//constructor to make a primer set from a array of primers in char array

	bool write_to_file(char * filename, ostream & err_msg = cout); //writes the list of primers in FASTA format to file

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); //prints out the nucleotide contribution of the sequence
	 
	bool show_All(ostream & out = cout, ostream &err_msg = cout); //prints out the show_statistics() and then the actual sequences

	unsigned int get_primer_as_value(unsigned int position, ostream & err_msg = cout); //copies the array of primers to the int * primers object

	bool get_primer_as_txt(char * primers, ostream & err_msg = cout); //copies the array of primers to the char ** primers object as text

	bool delete_primer(int position, ostream & err_msg = cout); //ability to mark a primer for deletion at a specific position in the array. 
//	bool delete_primer(int primer_value, ostream & err_msg = cout); //mark a primer for deletion based on its value
//	bool delete_primer(char * primer, ostream & err_msg = cout); //mark a primer for deletion based on its sequence

	bool add_primer(unsigned int primer_value, ostream & err_msg = cout); // add a primer to the array based on its integer value
	bool add_primer(char * primer, ostream & err_msg = cout); // add a primer to the array based on its sequence

	bool convert_primer_int_to_txt(iunsigned nt primer_value, char * primer, ostream & err_msg = cout);
	bool convert_primer_txt_to_int(char * primer, unsignedin & primer_value,  ostream & err_msg = cout);

	unsigned int * get_pointer_to_primer_array(ostream & err_msg = cout);

	unsigned int get_number_of_primers(ostream & err_msg = cout);
	unsigned int get_primer_length(ostream & err_msg = cout);




};