#include "RyR.hpp"


RyR::RyR(int RYR_Num, int ID,  int rand_seed, double Jmax_set)
	: random_num_RyR(rand_seed, ID), Jmax(Jmax_set)
{

	N_RyR = RYR_Num;
	RyR_1 = 0 + int (5.0 * random_num_RyR.gen_rand_uint()/ (double)(UINT_MAX));
	RyR_2 = 0;
	RyR_3 = 0;

}



RyR::~RyR() {}


double RyR::Update_RyR_stochastic(double dt, double Caj, double CaSRj) {


	// KK = 2400;  // was KK = 1400; KK = 2400 cause two peaks in INCX
	// KK = 2000;  // was KK = 1400; KK = 2400 cause two peaks in INCX
	//RyR
	double rho = rhoinf * pow(double(CaSRj / 1000.0), hh) / (pow(KK / 1000, hh) + pow(double(CaSRj / 1000), hh));
	if (rho < 0.0000001)rho = 0.0000001; /// to avoid div small (for cjsr<200)
	double MM = (sqrt(1 + 8 * rho * BCSQN) - 1) / (4 * rho * BCSQN);
	// double ncjsr = MM * nM + (1 - MM) * nD;

	     // double rhopri=(hh*rhoinf*pow(double(CaSRj),double(hh-1))*(pow(double(KK),double(hh))+pow(double(CaSRj),double(hh)))-rhoinf*pow(double(CaSRj),double(hh))*hh*pow(double(CaSRj),double(hh-1)))/pow((pow(double(KK),double(hh))+pow(double(CaSRj),double(hh))),2);

	     // double rhopri=(hh*rhoinf*pow(double(CaSRj/1000.0),double(hh-1))*(pow(double(KK/1000.0),double(hh))+pow(double(CaSRj/1000.0),double(hh)))-rhoinf*pow(double(CaSRj/1000.0),double(hh))*hh*pow(double(CaSRj/1000.0),double(hh-1)))/pow((pow(double(KK/1000.0),double(hh))+pow(double(CaSRj/1000.0),(hh))),2)/1000.0;

	//gpu code
	// double rhopri = (hh * rhoinf * pow(double(CaSRj / 1000.0), double(hh - 1)) * (pow(double(KK / 1000.0), double(hh)) + pow(double(CaSRj / 1000.0), double(hh))) - rhoinf * pow(double(CaSRj / 1000.0), double(hh)) * hh * pow(double(CaSRj / 1000.0), double(hh - 1))) / pow((pow(double(KK / 1000.0), double(hh)) + pow(double(CaSRj / 1000.0), double(hh))), 2);
	// rhopri *= 0.001;


	// double dMdc = (((1.0 / 2.0) * pow(double(1 + 8 * rho * BCSQN), double(-1.0 / 2.0)) * 8 * rhopri * BCSQN * 4 * rho * BCSQN) - (sqrt(1 + 8 * rho * BCSQN) - 1) * 4 * rhopri * BCSQN) / pow(4 * rho * BCSQN, 2);
	// double dndc = dMdc * (nM - nD);

	// double Betajsr = 1 / (1 + (Kc * BCSQN * ncjsr + dndc * (CaSRj * Kc + CaSRj * CaSRj)) / ((Kc + CaSRj) * (Kc + CaSRj)));


	Po = (RyR_2 + RyR_3) / 100.0;

	double cp2 = Caj * Caj;

	#ifdef ___USE_ORG_PARAM
		double k12 = Ku * cp2;
		double k43 = Kb * cp2;
	#else

	#ifdef ___KOSRCA



	// 		double MaxSR = 15;
	// double MinSR = 1;
	// double ec50SR = 450;
	// double hkosrca = 2.5;

		double kCaSR = MaxSR - (MaxSR - MinSR) / (1 + pow((ec50SR / CaSRj), hkosrca));
		double koSRCa = 1 / kCaSR;
	#else
		/*const*/ double koSRCa = 1;
	#endif
		// ec50SR=900;
		// ec50SR = 1400;
		double kCaSR = MaxSR - (MaxSR - MinSR) / (1 + pow((ec50SR / CaSRj), hkosrca));
		/*double*/ koSRCa = 1 / kCaSR;

		// koSRCa=1;
		double sgmd = cp2 / (Kcp * Kcp + cp2);
		double k12 = /*1.5**/koSRCa * Ku * sgmd + pedk12*10;//+ 0.2*3.8e-4* cp2;  // 18:39:11, Thu, 22-August-2019, By Haibo enhanced leak
		// double k12 = 2.5*koSRCa * Ku * sgmd + pedk12*10;//+ 0.2*3.8e-4* cp2;  // 18:39:11, Thu, 22-August-2019, By Haibo enhanced leak
		double k43 = koSRCa * Kb * sgmd + pedk43*10;//+ 0.2*1e-5 * cp2 ;  // 18:39:11, Thu, 22-August-2019, By Haibo enhanced leak


		// double k12 = 3.8e-4* cp2 + pedk12;
		// double k43 = 1e-5 * cp2 + pedk43;

	#endif





	double k14 = MM / taub * BCSQN / BCSQN0;
	double k21 = 1 / tauc1/koSRCa;
	double k23 = MM / taub * BCSQN / BCSQN0;
	double k41 = 1 / tauu*koSRCa;
	double k34 = 1 / tauc2/koSRCa;
	double k32 = k41 * k12 / k43;

	int RyR_4 = N_RyR - RyR_1 - RyR_2 - RyR_3;
	/*#ifndef __RYR_UNIFORM
		RyR_4 = int(nryr * sigma) - RyR_1 - RyR_2 - RyR_3;
	#endif*/
	int ryr12 = Update_RyR_rates_stochastic(RyR_1, k12 * dt);
	int ryr14 = Update_RyR_rates_stochastic(RyR_1, k14 * dt);
	int ryr21 = Update_RyR_rates_stochastic(RyR_2, k21 * dt);
	int ryr23 = Update_RyR_rates_stochastic(RyR_2, k23 * dt);
	int ryr43 = Update_RyR_rates_stochastic(RyR_4, k43 * dt);
	int ryr41 = Update_RyR_rates_stochastic(RyR_4, k41 * dt);
	int ryr34 = Update_RyR_rates_stochastic(RyR_3, k34 * dt);
	int ryr32 = Update_RyR_rates_stochastic(RyR_3, k32 * dt);
	RyR_1 = RyR_1 - (ryr12 + ryr14) + (ryr21 + ryr41);
	RyR_2 = RyR_2 - (ryr21 + ryr23) + (ryr12 + ryr32);
	RyR_3 = RyR_3 - (ryr34 + ryr32) + (ryr43 + ryr23);

	// check to make sure greater than 0 // 16:23:39, Thu, 06-December-2018, By Haibo
	if (RyR_1 < 0 || RyR_2 < 0 || RyR_3 < 0 || RyR_1 + RyR_2 + RyR_3 > N_RyR)
	{
		//          cout<<"RyR is negative "<<RyR_1<<"\t"<<RyR_2<<"\t"<<RyR_3<<endl;
		if (RyR_1 < 0)
		{
			if (random_num_RyR.gen_rand_uint() % 2)
				RyR_2 += RyR_1;
			RyR_1 = 0;
		}
		if (RyR_2 < 0)
		{
			if (random_num_RyR.gen_rand_uint() % 2)
			{
				RyR_1 += RyR_2;
				if (RyR_1 < 0)RyR_1 = 0;
			}
			else
				RyR_3 += RyR_2;
			RyR_2 = 0;
		}
		if (RyR_3 < 0)
		{
			if (random_num_RyR.gen_rand_uint() % 2)
			{
				RyR_2 += RyR_3;
				if (RyR_2 < 0)RyR_2 = 0;
			}
			RyR_3 = 0;
		}
		if (RyR_1 + RyR_2 + RyR_3 > N_RyR)
		{
			RyR_4 = N_RyR - (RyR_1 + RyR_2 + RyR_3);
			if (random_num_RyR.gen_rand_uint() % 2)
			{
				RyR_3 += RyR_4;
				if (RyR_3 < 0)
				{
					RyR_3 -= RyR_4;
					RyR_1 += RyR_4;
					if (RyR_1 < 0)
					{
						RyR_1 -= RyR_4;
						RyR_2 += RyR_4;
					}
				}
			}
			else
			{
				RyR_1 += RyR_4;
				if (RyR_1 < 0)
				{
					RyR_1 -= RyR_4;
					RyR_3 += RyR_4;
					if (RyR_3 < 0)
					{
						RyR_3 -= RyR_4;
						RyR_2 += RyR_4;
					}
				}
			}
		}
	}

	return Po;
}






// this function simulates the transition of RyR in the model/ stochastically .
int RyR::Update_RyR_rates_stochastic(double num, double p)
{
	int res;
	double lambda = num * p;


	if (lambda > 12)
	{
		//Gaussian
		double x1, x2, w;
		do
		{
			x1 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			x2 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);
		w = sqrt((-2.0 * log(w)) / w);
		double y1 = x1 * w;
		//double y2=x2*w;
		res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
		res = int(res + 0.5); //round
	}
	else if (100 * p < 6.6 + 52 * pow(num, double(-0.5)))
	{
		//Poisson
		double L = exp(-lambda);
		double k = 0;
		double pp = 1;
		do
		{
			k++;
			double u = random_num_RyR.gen_rand();
			pp *= u;
		} while (pp >= L);
		res = k - 1;
	}
	else
	{
		//Gaussian
		double x1, x2, w;
		do
		{
			x1 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			x2 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);
		w = sqrt((-2.0 * log(w)) / w);
		double y1 = x1 * w;
		//double y2=x2*w;
		res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
		res = int(res + 0.5); //round
	}

	if (res < 0)
		res = 0;

	return res;
}
