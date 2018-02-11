#include "PCR_Profile.h"
#include <iostream>

using namespace std;
class PCR_Profile_Toolbox
{
public:
	static bool merge_pcr_profiles(PCR_Profile * & resulting, PCR_Profile * input_a, PCR_Profile * input_b, ostream &err_msg = cout);
	//takes two PCR profiles and merges them into a new third one. As part of the process, it recalculates the statistics in the new PCR profile.

};