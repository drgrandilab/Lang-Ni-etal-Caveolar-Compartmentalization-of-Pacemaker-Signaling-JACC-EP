#ifndef ICAL13_HPP
#define ICAL13_HPP

#include <vector>
#include "xor_rand.hpp"


class ical13
{
public:
	ical13(int NCaL_set = 8, int ID = 0, int rand_seed = 0);
	~ical13();

	int NCaL;
	int open_NCaL_num;

	double gcal13;

	std::vector<int> fl13;
	std::vector<int> dl13;
	std::vector<int> fca;

	xor_rand random_num_LTCC;
	double ecal;
	double alpha_dl, beta_dl, alpha_fl, beta_fl, alpha_fca, beta_fca, alpha_dl_12, beta_dl_12;
	double v_cal, VDB;
	// double tau_dl, tau_fl, tau_fca;
	double kmfca;
	double tauS_DB_ca, tauF_DB_ca;
	double vDB;
	int update_ical13_gating_rates(double V, double Caj);
	int update_ical13_stochastic(std::vector<double> rand, double dt);
	int update_ical13_stochastic_version_2(std::vector<double> rand, double dt);
	int update_ical13_stochastic_version_2(double dt, double V, double Caj);
	void reset_random_generator(unsigned int seed = 0, unsigned int id = 0);
};


#endif