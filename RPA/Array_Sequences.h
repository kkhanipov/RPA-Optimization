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


public:
	Array_Sequences(char * filename, ostream &err_msg = cout); //constructor to be able to read a FASTA file and create the individual sequences and array of sequences
	Array_Sequences(unsigned int _max_number_of_sequences, ostream &err_msg = cout); //constructor to be able to read a FASTA file and create the individual sequences and array of sequences

	Array_Sequences(Array_Sequences * _arr_seq, ostream &err_msg = cout);

	bool show_Statistics(ostream & out = cout, ostream &err_msg = cout); //show basic statistics about the array of sequences, number of sequences present,  and their nucleotide contribution

	bool show_All(ostream & out = cout, ostream &err_msg = cout); // shows the basic statistics about the array of sequences and then prints to screen the actual sequences

	bool add_sequence(Sequence * sequence, ostream &err_msg = cout); //ability to add a sequence to the array of sequences by sequence object
	bool add_sequence(char * sequence, ostream &err_msg = cout); //ability to add a sequence to the array of sequences by giving the actual sequence

	unsigned int get_number_of_sequences(ostream &err_msg = cout);
	Sequence * get_pointer_to_sequence_object(unsigned int position, ostream &err_msg = cout);
};

Array_Sequences::Array_Sequences(char * filename, ostream &err_msg)
{
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

	unsigned int max_sequence_length = 4000000000;

	char * _sequence = new char[max_sequence_length + 1];
	if (_sequence == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==>  Not enough memory (_sequence == NULL)" << endl;
		in.close(); delete[] _header_text; 	throw ("Failed to read sequence");
	}


	unsigned int i;

	
	// Create sequence length array
	sequence_length = new unsigned int[array_size];
	if (sequence_length == NULL)
	{
		err_msg_out << "ERROR: Array_of_Sequences::read_Array_of_Sequences_From_FASTA_File(...) ==> Not enough memory (sequence_length == NULL)" << endl;
		in.close(); delete[] _header_text; delete[] _sequence; delete[] copy_number; return false;
	}


	// Create sequences array
	sequence = new char *[array_size];
	if (sequence == NULL)
	{
		err_msg_out << "ERROR: Array_of_Sequences::read_Array_of_Sequences_From_FASTA_File(...) ==> Not enough memory (sequence == NULL)" << endl;
		in.close(); delete[] _header_text; delete[] _sequence; delete[] copy_number; delete[] sequence_length; return false;
	}

	// Initialize 
	for (i = 0; i<array_size; i++) { copy_number[i] = 0; sequence_length[i] = 0; sequence[i] = NULL; }

	unsigned int ii;

	unsigned int position = 0;

	// Read sequences in file
	while (!in.eof())
	{
		in.getline(_header_text, max_text_length);
		if (_header_text[0] == '\0' || _header_text[0] == EOF) break;

		in.getline(_sequence, max_sequence_length + 1);

		for (i = 0; i<max_sequence_length + 1; i++) if (_sequence[i] == '\0') break;
		sequence_length[position] = i;

		sequence[position] = new char[i + 1];
		if (sequence[position] == NULL)
		{
			err_msg_out << "ERROR: Array_of_Sequences::read_Array_of_Sequences_From_FASTA_File(...) ==> Not enough memory (sequence[" << i << "] == NULL)" << endl;
			in.close(); delete[] _header_text; delete[] _sequence;  delete[] copy_number; delete[] sequence_length;
			for (ii = 0; ii<i; ii++) if (sequence[i] != NULL) delete[] sequence[ii]; // Delete all allocated sequences
			delete[] sequence; return false;
		}

		sequence[position][i] = '\0';

		// Copy new sequence in the object
		for (i = 0; i< sequence_length[position]; i++) sequence[position][i] = _sequence[i];
		sequence[position][sequence_length[position]] = '\0';
		copy_number[position] = 1;
		position++;

		// Increase array size if needed
		if (position == array_size - 1)
		{
			if (!z_increase_Array_Size())
			{
				err_msg_out << "ERROR: Array_of_Sequences::read_Array_of_Sequences_From_FASTA_File(...) ==> z_increase_Array_Size() == false" << endl;
				in.close(); delete[] _header_text; delete[] _sequence; delete[] copy_number; delete[] sequence_length;
				for (i = 0; i<array_size; i++) if (sequence[i] != NULL) delete[] sequence[i];	// Delete all allocated sequences
				delete[] sequence; return false;
			}
		}
	}


	delete[] _header_text;
	delete[] _sequence;

	return true;
}