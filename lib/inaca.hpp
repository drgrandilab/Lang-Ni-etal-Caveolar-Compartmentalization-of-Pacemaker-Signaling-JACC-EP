#ifndef INACA_HPP
#define INACA_HPP

#include <math.h>

class inaca
{
public:
	inaca(double nai, double nao, double V, double cao);

	double compute_NCX(double nai,  double nao, double Caj, double cao);
	double compute_NCX_version_2(double v, double nai,  double nao, double Caj, double cao);
	double compute_NCX_version_2_backup(double v, double nai,  double nao, double Caj, double cao);
	~inaca() {};
	//NCX (ca independent part)
	const double Kmnai = 12.3; // [mM]
	// const double nao = 136; //  [mM]
	const double Kmcai = 3.59; //  [uM]  // 18:39:00, Wed, 19-December-2018, By Haibo
	const double ksat = 0.27;
	const double eta = 0.35;

	constexpr static double F = 96.5; // [C/mmol]
	constexpr static double R = 8.314; // [J/mol*K]
	constexpr static double T = 308; // [K]
	constexpr static double rtf = R * T / F; //~26.5
	constexpr static double rtf2 = R * T / (2 * F);
	constexpr static  double vnaca = 7.452 *4/** 4*50*/;

	double t1 ;
	double x1a;
	double x1b;
	double x2;
};


#endif