

#ifndef ICAT_HPP
#define ICAT_HPP

class iCaT
{
public:
	iCaT();
	~iCaT(){};



	double compute_ICaT(double V, double dt);
	double ICaT;
	double gCaT, v_cat;

	double d_CaT, f_CaT; //gates for ICaT
};



#endif