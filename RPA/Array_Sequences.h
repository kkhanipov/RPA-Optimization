#ifndef Array_Sequences_H
#define Array_Sequences_H

#include "Sequence.h"
#include <iostream>
#include <fstream>

using namespace std;

//This class is meant to store a list of sequences read in from FASTA files. 
//The class uses class Sequences for each individual sequence, and this is an array that ties them all together.

class Array_Sequences
{
	Sequence ** sequence; //an array of pointers to sequence objects
	unsigned int max_number_of_sequences; //this is the array size
	unsigned int number_of_sequences; // number of sequences located inside the array
	unsigned int max_sequence_length;

public:
	Array_Sequences(char * filename, ostream &err_msg = cout); //constructor to be able to read a FASTA file and create the individual sequences and array of sequences
	Array_Sequences(unsigned int _max_number_of_sequences, ostream &err_msg = cout); //constructor to be able to read a FASTA file and create the individual sequences and array of sequences

	Array_Sequences(Array_Sequences * _arr_seq, ostream &err_msg = cout);
	~Array_Sequences();
	bool show_Statistics(ostream & out = cout, ostream &err_msg = cout); //show basic statistics about the array of sequences, number of sequences present,  and their nucleotide contribution

	bool show_All(ostream & out = cout, ostream &err_msg = cout); // shows the basic statistics about the array of sequences and then prints to screen the actual sequences

	bool add_sequence(Sequence * _sequence, ostream &err_msg = cout); //ability to add a sequence to the array of sequences by sequence object
	bool add_sequence(char * _sequence, unsigned int _sequence_length, ostream &err_msg = cout); //ability to add a sequence to the array of sequences by giving the actual sequence
	
	unsigned int get_number_of_sequences() { return number_of_sequences; }
	Sequence * get_pointer_to_sequence_object(unsigned int position) { return sequence[position]; }
};
Array_Sequences::~Array_Sequences()
{
	if (sequence != NULL)
	{

		for (int i = 0; i < number_of_sequences; i++) sequence[i]->~Sequence();
		delete[] sequence;
	}
}
Array_Sequences::Array_Sequences(unsigned int _max_number_of_sequences, ostream &err_msg)
{
	number_of_sequences = 0;
	max_number_of_sequences = _max_number_of_sequences;
	max_sequence_length = 4000000000;
	sequence = new Sequence*[max_number_of_sequences];
}

Array_Sequences::Array_Sequences(char * filename, ostream &err_msg)
{
	number_of_sequences = 0;
	max_number_of_sequences = 5000;
	max_sequence_length = 4000000000;
	sequence = new Sequence*[max_number_of_sequences];
	// Check if input file is opened successfully
	if (filename == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> filename == NULL" << endl; throw ("No File Name Provided");
	}

	ifstream in;
	in.open(filename);

	if (!in.is_open())
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> Cannot open FASTA file" << endl; throw ("Failed to open file");
	}

	char * _header_text = new char[10000 + 1];
	if (_header_text == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> Not enough memory (_header_text == NULL)" << endl; in.close();
		throw ("Failed to read header");
	}

	char * _sequence = new char[max_sequence_length + 1];
	if (_sequence == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==>  Not enough memory (_sequence == NULL)" << endl;
		in.close(); delete[] _header_text; 	throw ("Failed to read sequence");
	}


	unsigned int i;
	unsigned int ii;

	// Read sequences in file
	while (!in.eof())
	{
		in.getline(_header_text, max_sequence_length);
		if (_header_text[0] == '\0' || _header_text[0] == EOF) break;

		in.getline(_sequence, max_sequence_length + 1);

		for (i = 0; i<max_sequence_length + 1; i++) if (_sequence[i] == '\0') break;
				
		sequence[number_of_sequences] = new Sequence(_sequence, i);
		

		if (sequence[number_of_sequences] == NULL)
		{
			err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> Not enough memory (sequence[" << i << "] == NULL)" << endl;
			in.close(); delete[] _header_text; delete[] _sequence; 
			for (ii = 0; ii < i; ii++) if (sequence[i] != NULL) {} ; // Delete all allocated sequences
			delete[] sequence; throw ("Not enough memory to allocate for sequence");
		}

		
		number_of_sequences++;
		// Increase array size if needed
		if (number_of_sequences == max_number_of_sequences - 1)
		{
			//realloc stuff
		}
	}


	delete[] _header_text;
	delete[] _sequence;
}

bool Array_Sequences::add_sequence(char * _sequence, unsigned int _sequence_length, ostream &err_msg)
{
	if (number_of_sequences > max_number_of_sequences - 1)
	{
		err_msg << "ERROR: Array_Sequences::add_sequence(...) ==> Not array size needs to be increased == NULL)" << endl;
		return false;
	}
	sequence[number_of_sequences] = new Sequence(_sequence, _sequence_length);
	number_of_sequences++;
	return true;
}

bool Array_Sequences::add_sequence(Sequence * _sequence, ostream &err_msg)
{
	if (number_of_sequences > max_number_of_sequences - 1)
	{
		err_msg << "ERROR: Array_Sequences::add_sequence(...) ==> Not array size needs to be increased == NULL)" << endl;
		return false;
	}
	sequence[number_of_sequences] = new Sequence(_sequence->get_pointer_to_sequence(), _sequence->get_sequence_length());
	number_of_sequences++;
	return true;
}
bool Array_Sequences::show_Statistics(ostream & out, ostream &err_msg)
{
	for (int i = 0; i < number_of_sequences; i++)
	{
		sequence[i]->show_statistics();
	}
	return true;
}
bool Array_Sequences::show_All(ostream & out, ostream &err_msg)
{
	for (int i = 0; i < number_of_sequences; i++)
	{
		sequence[i]->show_All();
	}
	return true;
}

#endif