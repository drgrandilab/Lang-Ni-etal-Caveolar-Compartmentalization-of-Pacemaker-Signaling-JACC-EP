

#ifndef RyR_HPP

#define RyR_HPP

#include "xor_rand.hpp"



class RyR
{
public:
	RyR(int RYR_Num = 100, int ID = 0, int rand_seed = 0,  double Jmax_set = 0.0147 * 18);
	~RyR();

	int RyR_1, RyR_2, RyR_3;
	xor_rand random_num_RyR;
	double Po;

	int N_RyR;

	double Jmax;

	double Update_RyR_stochastic(double dt, double Caj, double CaSRj);

	int Update_RyR_rates_stochastic(double num, double p);

	void set_Jmax(double Jmax_set) {
		Jmax = Jmax_set;
	}

	double MaxSR = 15;
	double MinSR = 1;
	double ec50SR = 450;
	double hkosrca = 2.5;


	double Ku = 5.0;
	double Kb = 0.005;
	// double tauu = 1250.0;  // could change this number to reduce fractoriness  // 14:21:56, Wed, 27-November-2019, By Haibo
	double tauu = 1250.0/3.0; // // 17:50:56, Thu, 09-January-2020, By Haibo
	// double tauu = 1250.0/4.0; /// 21:14:23, Mon, 22-June-2020, By Haibo
	double taub = 0.5;
	double tauc1 = 2.0;
	double tauc2 = 0.3;
	double hh = 10.0;
	double KK = 1400;


	double taup = 0.022;
	double BCSQN = 400;
	// double tautr = 5.0;
	double Kc = 600.0;
	double nM = 15;
	double nD = 35;
	double rhoinf = 5000;
	double BCSQN0 = 400;
	// double Kcp = 100;
	// double Kcp = 10;//100;
	// double Kcp = 10;//100;
	double Kcp = 5;//100;
	// double Kcp = 20;//100;

	double pedk12 = 0.000001;
	double pedk43 = 0.000001;
	// double gcabk = 0.0002513;
	// double qslcap = 2.35 * 4;
};



#endif