#include <iostream>
#include "Array_Sequences.h"

using namespace std;

int main()
{
	Array_Sequences * as;
	as = new Array_Sequences("test.fasta");
	as->show_Statistics();
	as->show_All();

	system("PAUSE");
	return 1;
}