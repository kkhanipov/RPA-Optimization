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

	char * get_pointer_to_sequence() {return sequence; };

};

Sequence::Sequence(char * sequence, unsigned int sequence_length, ostream & err_msg);
{

}
