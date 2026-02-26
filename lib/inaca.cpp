

#include "inaca.hpp"


inaca::inaca(double nai, double nao, double V, double cao) {

	const double z = V * F / (R * T);
	t1 = (Kmcai * 0.001) * nao * nao * nao * (1 + (nai * nai * nai / (Kmnai * Kmnai * Kmnai)));
	x1a = exp(eta * z) * nai * nai * nai * cao;
	x1b = exp((eta - 1) * z) * nao * nao * nao;
	x2 = (1 + ksat * exp((eta - 1) * z));
}


double inaca::compute_NCX(double nai,  double nao, double Caj, double cao) {

	const double Kmcao = 1.3;
	const double Kmnao = 87.5;
	const double Kda = 0.11;


	double x3 =  vnaca; //x3; //x3 is vnaca, x3_new=x3/svr, svr is the rate #surface CRU / #all CRU
	double Ka = 1 / (1 + (Kda * Kda * Kda / (Caj * Caj * Caj)));
	double t2 = Kmnao * Kmnao * Kmnao * (Caj * 0.001) * (1 + Caj / Kmcai);
	double t3 = Kmcao * nai * nai * nai + nai * nai * nai * cao + nao * nao * nao * (Caj * 0.001);
	return  (Ka * x3 * (x1a - x1b * (Caj * 0.001)) / ((t1 + t2 + t3) * x2));

}


double inaca::compute_NCX_version_2(double v, double nai,  double nao, double Caj, double cao) {

	// Caj = Caj/3000.0;  // convert mM to uM  // 18:18:41, Wed, 19-December-2018, By Haibo  scale Caj by dividing 3000
	Caj = Caj/1000.0;  // convert mM to uM  // 18:18:41, Wed, 19-December-2018, By Haibo

	double R = 8.314472; //  [J/mol*K]
	double T = 310.5; //  [K] 
	double F = 96.4845; //  [C/kmol]
	double FRT = F/(R*T);
	// print FRT
	// inakmax = 1 * 1.85 * 0.077;
	double kmnap = 14.0;
	double kmkp = 1.4;
	double K1ni = 395.3;
	double K1no = 1628.0;
	double K2ni = 2.289;
	double K2no = 561.4;
	double K3ni = 26.44;
	double K3no = 4.663;
	// double Kci = 0.0207;
	double Kci = 0.0207;// *3;  // 18:58:05, Wed, 27-November-2019, By Haibo
	// double Kci = 0.0207*3;  // 18:58:05, Wed, 27-November-2019, By Haibo
	double Kco = 3.663;
	double Kcni = 26.44;
	double Qci = 0.1369;
	double Qco = 0.0;
	double Qn = 0.4315;
	double tdifca = 0.04;
	double Ttr = 40.0;
	double kNaCa = 5.5;

	double di = 1 + (Caj / Kci) * (1 + exp(-Qci * v * FRT) + nai / Kcni) + (nai / K1ni) * (1 + (nai / K2ni) * (1 + nai / K3ni));

	double doo =  1 + (cao / Kco) * (1 + exp(Qco * v * FRT)) +(nao / K1no) * (1 + (nao / K2no) * (1 + nao / K3no))  ;
	// print doo
	double k43 = nai / (K3ni + nai);
	double k12 = (Caj / Kci) * exp(-Qci * v * FRT) / di;
	double k14 = (nai / K1ni) * (nai / K2ni) * (1 + nai / K3ni) * exp(Qn * v * FRT / 2.0) / di;
	double k41 = exp(-Qn * v * FRT / 2.0);
	double k34 = nao / (K3no + nao);
	double k21 = (cao / Kco) * exp(Qco * v * FRT) / doo;
	double k23 = (nao / K1no) * (nao / K2no) * (1 + nao / K3no) * exp(-Qn * v * FRT / 2.0) / doo;
	double k32 = exp(Qn * v * FRT / 2.0);
	double x1 = k34 * k41 * (k23 + k21) + k21 * k32 * (k43 + k41);
	double x2 = k43 * k32 * (k14 + k12) + k41 * k12 * (k34 + k32);
	double x3 = k43 * k14 * (k23 + k21) + k12 * k23 * (k43 + k41);
	double x4 = k34 * k23 * (k14 + k12) + k21 * k14 * (k34 + k32);

	// return /*1.5*2**/40.0 * kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);  // 16:08:26, Wed, 27-November-2019, By Haibo
	// return 1.5*2*40.0 * kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);
	return 1.5*/*2**/40.0 * kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);

}



double inaca::compute_NCX_version_2_backup(double v, double nai,  double nao, double Caj, double cao) {

	// Caj = Caj/3000.0;  // convert mM to uM  // 18:18:41, Wed, 19-December-2018, By Haibo  scale Caj by dividing 3000
	Caj = Caj/1000.0;  // convert mM to uM  // 18:18:41, Wed, 19-December-2018, By Haibo  scale Caj by dividing 3000

	double R = 8.314472; //  [J/mol*K]
	double T = 310.5; //  [K] 
	double F = 96.4845; //  [C/kmol]
	double FRT = F/(R*T);
	// print FRT
	// inakmax = 1 * 1.85 * 0.077;
	double kmnap = 14.0;
	double kmkp = 1.4;
	double K1ni = 395.3;
	double K1no = 1628.0;
	double K2ni = 2.289;
	double K2no = 561.4;
	double K3ni = 26.44;
	double K3no = 4.663;
	// double Kci = 0.0207;
	double Kci = 0.0207*3;  // 18:58:05, Wed, 27-November-2019, By Haibo
	double Kco = 3.663;
	double Kcni = 26.44;
	double Qci = 0.1369;
	double Qco = 0.0;
	double Qn = 0.4315;
	double tdifca = 0.04;
	double Ttr = 40.0;
	double kNaCa = 5.5;

	double di = 1 + (Caj / Kci) * (1 + exp(-Qci * v * FRT) + nai / Kcni) + (nai / K1ni) * (1 + (nai / K2ni) * (1 + nai / K3ni));

	double doo =  1 + (cao / Kco) * (1 + exp(Qco * v * FRT)) +(nao / K1no) * (1 + (nao / K2no) * (1 + nao / K3no))  ;
	// print doo
	double k43 = nai / (K3ni + nai);
	double k12 = (Caj / Kci) * exp(-Qci * v * FRT) / di;
	double k14 = (nai / K1ni) * (nai / K2ni) * (1 + nai / K3ni) * exp(Qn * v * FRT / 2.0) / di;
	double k41 = exp(-Qn * v * FRT / 2.0);
	double k34 = nao / (K3no + nao);
	double k21 = (cao / Kco) * exp(Qco * v * FRT) / doo;
	double k23 = (nao / K1no) * (nao / K2no) * (1 + nao / K3no) * exp(-Qn * v * FRT / 2.0) / doo;
	double k32 = exp(Qn * v * FRT / 2.0);
	double x1 = k34 * k41 * (k23 + k21) + k21 * k32 * (k43 + k41);
	double x2 = k43 * k32 * (k14 + k12) + k41 * k12 * (k34 + k32);
	double x3 = k43 * k14 * (k23 + k21) + k12 * k23 * (k43 + k41);
	double x4 = k34 * k23 * (k14 + k12) + k21 * k14 * (k34 + k32);

	// return /*1.5*2**/40.0 * kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);  // 16:08:26, Wed, 27-November-2019, By Haibo
	// return 1.5*2*40.0 * kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);
	return 1.5*2*40.0 * kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);

}