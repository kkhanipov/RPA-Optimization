#include <iostream>

using namespace std;

// this class is meant to store sequences read from a FASTA file
class Sequence
{
	char * sequence;
	unsigned int seq_length;

public:
	
	Sequence(char * sequence, unsigned int sequence_length, ostream & err_msg=cout);

	bool show_statistics(ostream & out=cout, ostream &err_msg=cout); //prints out the nucleotide contribution of the sequence
	bool show_All(ostream & out = cout, ostream &err_msg = cout); //prints out the show_statistics() and then the actual sequence

	char * get_pointer_to_sequence() {return sequence;};

};

Sequence::Sequence(char * _sequence, unsigned int _sequence_length, ostream & err_msg)
{
	if (_sequence_length == 0)
	{
		err_msg << "ERROR: Sequence::Sequence ==> _sequence_length == 0" << endl; throw ("Sequence could not be allocated b/c length is 0");
	}
	sequence = new char[_sequence_length];
	seq_length = _sequence_length;
	if (sequence == NULL)
	{
		err_msg << "ERROR: Sequence::Sequence ==> sequence == NULL" << endl; throw ("Sequence could not be allocated");
	}
	for (int i = 0; i < _sequence_length; i++)
	{
		sequence[i] = _sequence[i];
	}
}

bool Sequence::show_statistics(ostream & out, ostream &err_msg)
{
	int count_A = 0;
	int count_C = 0;
	int count_T = 0;
	int count_G = 0;
	int count_N = 0;
	for (int i = 0; i < seq_length; i++)
	{
		switch (sequence[i])
		{
		case 'A':case 'a': count_A++; break;
		case 'T':case 't': count_T++; break;
		case 'C':case 'c': count_C++; break;
		case 'G':case 'g': count_G++; break;
		case 'N':case 'n': count_N++; break;
		default:
			err_msg << "ERROR: Sequence::Sequence ==> unexpected character: sequence[" << i << "]=" << sequence[i] << endl;
			return false;
		}
	}

	out << "Nucleotide Contribution is:" << endl;
	out << "A :"<< count_A << endl;
	out << "C :" << count_C << endl;
	out << "T :" << count_T << endl;
	out << "G :" << count_G << endl;
	out << "N :" << count_N << endl;

	return true;
}

bool Sequence::show_All(ostream & out, ostream &err_msg)
{
	show_statistics();

	for (int i = 0; i < seq_length; i++)
	{
		out << sequence[i];
	}
	out << endl;

	return true;
}