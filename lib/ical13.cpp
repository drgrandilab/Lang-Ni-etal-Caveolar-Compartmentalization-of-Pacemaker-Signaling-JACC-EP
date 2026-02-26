
#include "ical13.hpp"
#include <math.h>
#include <stdio.h>


ical13::ical13(int NCaL_set, int ID, int rand_seed)
	: random_num_LTCC(rand_seed, ID)  // 1 for seed of IcaL generator, different from RyR
{
	v_cal = 0;
	vDB   = 0;
	ecal  = 47.0;
	kmfca = 0.8;// 0.35;
	gcal13 = /*1.5* */ 0.0030 * 4.0 * 1.5 * 2.0 / 3.0 * 1.261; // from Stefano's code
	tauS_DB_ca = 1.0;
	tauF_DB_ca = 1.0;
	alpha_fca = 0.021;

	open_NCaL_num = 0;

	NCaL = NCaL_set;

	for (int i = 0; i < NCaL; ++i) {
		fl13.push_back(1);
		dl13.push_back(0);
		fca .push_back(1);
	}
}

ical13::~ical13() {

}


int ical13::update_ical13_gating_rates(double V, double Caj) {

// ICaL - L-type Ca channel isoforms Cav1.2/Cav1.3 currents ***************
	if (fabs(V + v_cal + vDB) <= 0.001) {
		alpha_dl = -28.39 * (V + 35.0 + v_cal + vDB) / (exp(-(V + 35.0 + v_cal + vDB) / 2.5) - 1.0) + 408.173;
	}
	else if (fabs(V + 35.0 + v_cal + vDB) <= 0.001) {
		alpha_dl = 70.975 - 84.9 * (V + v_cal + vDB) / (exp(-0.208 * (V + v_cal + vDB)) - 1.0);
	}
	else {// if (fabs(V)>0.001 && fabs(V+35.0)>0.001),
		alpha_dl = -28.39 * (V + 35.0 + v_cal + vDB) / (exp(-(V + 35.0 + v_cal + vDB) / 2.5) - 1.0) - 84.9 * (V + v_cal + vDB) / (exp(-0.208 * (V + v_cal + vDB)) - 1.0);
	}
	if (fabs(V - 5.0 + v_cal + vDB) <= 0.001) {
		beta_dl = 28.575;
	} else {
		beta_dl = 11.43 * (V - 5.0 + v_cal + vDB) / (exp(0.4 * (V - 5.0 + v_cal + vDB)) - 1.0);
	}
	double tau = 2000.0 / (alpha_dl + beta_dl);

	// printf("%f\n", tau);

	double dl13_inf = 1.0 / (1 + exp(-(V + 13.5 + v_cal + vDB) / 6.0));

	alpha_dl = dl13_inf / tau;
	beta_dl = (1 - dl13_inf) / tau;

	// tau_fl = 7.4 + 45.77 * exp(-0.5 * (V + 28.1 + v_bas + vDB) * (V + 28.1 + v_bas + vDB) / (11 * 11));
	tau = tauS_DB_ca * (7.4 + 45.77 * exp(-0.5 * (V + 28.1 + v_cal + vDB) * (V + 28.1 + v_cal + vDB) / (11 * 11)));
	// dl13_inf = 1.0 / (1 + exp(-(V + 13.5) / 6.0));
	double fl13_inf = 1.0 / (1 + exp((V + 35.0 + v_cal + vDB) / 7.3));

	alpha_fl = fl13_inf / tau;
	beta_fl = (1 - fl13_inf) / tau;

	// dl12_inf = 1.0 / (1 + exp(-(V + 3.0 + v_cal + vDB) / 5.0));
	// fl12_inf = 1.0 / (1 + exp((V + 36.0 + v_cal + vDB) / 4.6));
	// dl12_dot = (dl12_inf - dl12) / tau_dl;
	// fl12_dot = (fl12_inf - fl12) / tau_fl;
	// double fca_inf = kmfca / (kmfca + Caj);
	// double Caj_bar = 1.0;// was 6.0; // From Song et al. 2015 Biophysical J
	double Caj_bar = 3.0;// was 6.0; // From Song et al. 2015 Biophysical J
	// double Caj_bar = 1.0;// // 20:08:26, Mon, 22-June-2020, By Haibo was 6.0; // From Song et al. 2015 Biophysical J
	double fca_inf = 1.0 / (1.0  +  Caj / Caj_bar) / (1.0 + Caj / Caj_bar); 


	const double cat = 0.5;
	// double fca = 1.0 / (1.0 + pow(double(cat / Caj), 3));

	double fca = 1.0 / (1.0 + pow(double(cat / Caj), 3));

	//      double s1=0.0182688*fca;
	double s1 = 0.02 * fca; //Juan code param
	// fca_inf = fca;

	// taufca = fca_inf / alpha_fca;
	tau = 15.0;// 15 ms from Song et al. 2015 BJ; tauF_DB_ca * fca_inf / 0.021;

	alpha_fca = fca_inf / tau ;
	beta_fca = (1 - fca_inf) / tau;

	alpha_dl_12 = 0.3;  // from Colman et al. 2017 Plos comp and Song et al. Biophysical J 2015.
	beta_dl_12 = 6.0;
	return 1;
}



int ical13::update_ical13_stochastic(std::vector<double> rand, double dt)  {

	int open_NCaL_num = 0;

	for (int i = 0; i < NCaL; ++i) {
		double rand_value = rand[i]; // assign rand to  randomly generated number
		if (rand_value > 1.0) printf("RAND GREATER THAN 1");
		else if (rand_value < 0.0) printf("RAND LESS THAN 0");


		if (dl13[i] == 0) { // closed state 0 voltage activation
			if (rand_value <= alpha_dl * dt) dl13[i] = 1; // if rand if less than alpha rate va state 0-1, transition to state 1
			else {
				if (fl13[i] == 0) // if voltage inactivation state is in 0
				{
					if (rand_value <= (alpha_fl + alpha_dl) * dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + alpha_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + alpha_dl) * dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + alpha_dl) * dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + alpha_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + alpha_dl) * dt) fca[i] = 0;
						}
					}
				}
			}
		} else if (dl13[i] == 1)
		{
			if (rand_value <= beta_dl * dt) dl13[i] = 0;
			else
			{
				if (fl13[i] == 0)
				{
					if (rand_value <= (alpha_fl + beta_dl)*dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + beta_dl)*dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + beta_dl)*dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + beta_dl)*dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + beta_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + beta_dl ) * dt) fca[i] = 0;
						}
					}
				}
			}
		}

		if (dl13[i] == 1 and fl13[i] == 1 and fca[i] == 1) open_NCaL_num ++;
	}

	return open_NCaL_num;

}




// LTCC stochastic
int ical13::update_ical13_stochastic_version_2(std::vector<double> rand, double dt)
{
	int open_NCaL_num = 0;

	for (int i = 0; i < NCaL; i++)
	{
		double rand_value = rand[i]; // assign rand to  randomly generated number
		if (rand_value > 1.0) printf("RAND GREATER THAN 1");
		else if (rand_value < 0.0) printf("RAND LESS THAN 0");

		// printf("alpha_dl = %f\n", alpha_dl*dt);
		// printf("alpha_dl = %f\n", beta_dl*dt);
		// This works by cyling through states, if no transition occurs in first state, check for transition in second state, then third
		if (dl13[i] == 0) // closed state 0 voltage activation
		{


			if (rand_value <= alpha_dl * dt) dl13[i] = 1; // if rand_value if less than alpha rate va state 0-1, transition to state 1
			else
			{
				if (fl13[i] == 0) // if voltage inactivation state is in 0
				{
					if (rand_value <= (alpha_fl + alpha_dl) * dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + alpha_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + alpha_dl) * dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + alpha_dl) * dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + alpha_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + alpha_dl) * dt) fca[i] = 0;
						}
					}
				}
			}
		}
		else if (dl13[i] == 1)
		{


			if (rand_value <= beta_dl * dt) {dl13[i] = 0;}
			else if (rand_value <= (alpha_dl_12 + beta_dl) * dt) {dl13[i] = 2;}
			else
			{
				if (fl13[i] == 0)
				{
					if (rand_value <= (alpha_fl + alpha_dl_12 + beta_dl) * dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + alpha_dl_12 + beta_dl) * dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 0;
						}
					}
				}
			}
		}
		else if (dl13[i] == 2)
		{
			if (rand_value <= beta_dl_12 * dt) dl13[i] = 1;
			else
			{
				if (fl13[i] == 0)
				{
					if (rand_value <= (alpha_fl + beta_dl_12)*dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + beta_dl_12)*dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + beta_dl_12)*dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + beta_dl_12)*dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + beta_dl_12) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + beta_dl_12 ) * dt) fca[i] = 0;
						}
					}
				}
			}
		}

		if (dl13[i] == 2 and fl13[i] == 1 and fca[i] == 1) {open_NCaL_num ++; } // if open state 2, inactivation state 1 for v and Ca, then channel is open

	} // end for
	return open_NCaL_num;
} // End LTCC stochastic




// LTCC stochastic
int ical13::update_ical13_stochastic_version_2(double dt, double V, double Caj)
{
	update_ical13_gating_rates(V, Caj);


	int open_NCaL_num = 0;

	for (int i = 0; i < NCaL; i++)
	{
		double rand_value = random_num_LTCC.gen_rand(); // assign rand to  randomly generated number
		if (rand_value > 1.0) printf("RAND GREATER THAN 1");
		else if (rand_value < 0.0) printf("RAND LESS THAN 0");

		// printf("alpha_dl = %f\n", alpha_dl*dt);
		// printf("alpha_dl = %f\n", beta_dl*dt);
		// This works by cyling through states, if no transition occurs in first state, check for transition in second state, then third
		if (dl13[i] == 0) // closed state 0 voltage activation
		{


			if (rand_value <= alpha_dl * dt) dl13[i] = 1; // if rand_value if less than alpha rate va state 0-1, transition to state 1
			else
			{
				if (fl13[i] == 0) // if voltage inactivation state is in 0
				{
					if (rand_value <= (alpha_fl + alpha_dl) * dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + alpha_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + alpha_dl) * dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + alpha_dl) * dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + alpha_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + alpha_dl) * dt) fca[i] = 0;
						}
					}
				}
			}
		}
		else if (dl13[i] == 1)
		{


			if (rand_value <= beta_dl * dt) {dl13[i] = 0;}
			else if (rand_value <= (alpha_dl_12 + beta_dl) * dt) {dl13[i] = 2;}
			else
			{
				if (fl13[i] == 0)
				{
					if (rand_value <= (alpha_fl + alpha_dl_12 + beta_dl) * dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + alpha_dl_12 + beta_dl) * dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + alpha_dl_12 + beta_dl) * dt) fca[i] = 0;
						}
					}
				}
			}
		}
		else if (dl13[i] == 2)
		{
			if (rand_value <= beta_dl_12 * dt) dl13[i] = 1;
			else
			{
				if (fl13[i] == 0)
				{
					if (rand_value <= (alpha_fl + beta_dl_12)*dt) fl13[i] = 1;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + alpha_fl + beta_dl_12)*dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + alpha_fl + beta_dl_12)*dt) fca[i] = 0;
						}
					}
				}
				else if (fl13[i] == 1)
				{
					if (rand_value <= (beta_fl + beta_dl_12)*dt) fl13[i] = 0;
					else
					{
						if (fca[i] == 0)
						{
							if (rand_value <= (alpha_fca + beta_fl + beta_dl_12) * dt) fca[i] = 1;
						}
						else if (fca[i] == 1)
						{
							if (rand_value <= (beta_fca + beta_fl + beta_dl_12 ) * dt) fca[i] = 0;
						}
					}
				}
			}
		}

		if (dl13[i] == 2 and fl13[i] == 1 and fca[i] == 1) {open_NCaL_num ++; } // if open state 2, inactivation state 1 for v and Ca, then channel is open

	} // end for
	return open_NCaL_num;
} // End LTCC stochastic



void ical13::reset_random_generator(unsigned int seed, unsigned int id) {
	random_num_LTCC.reset(seed, id);
}