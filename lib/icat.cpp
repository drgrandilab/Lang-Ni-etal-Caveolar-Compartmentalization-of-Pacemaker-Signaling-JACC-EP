/*
* Implementation of ICaT

*/


#include "icat.hpp"
#include "ExplicitSolver.hpp"
#include <math.h>


iCaT::iCaT() {
	// ICaT parameters
	gCaT =  1.890 * 0.75 * 0.01862 / 0.025;  // nS/pF
	v_cat = 0.0; //

	d_CaT = 0.0017968124;
	f_CaT = 0.3660355707;
	// 	y[2] = 0.0017968124;
	// y[3] = 0.3660355707;
}


// this function computes ICaT and update gating variables for iCaT; // 15:16:15, Thu, 13-December-2018, By Haibo
double iCaT::compute_ICaT(double V, double dt) {

	const double ecat = 45.0;
	double tau_dt = 1.0 / (1.068 * exp((V + 26.3 + v_cat) / 30.0) + 1.068 * exp(-(V + 26.3 + v_cat) / 30.0));
	double dt_inf = 1.0 / (1.0 + exp(-(V + 26.0 + v_cat) / 6.0));
	// double dt_dot = (dt_inf - d_CaT) / tau_dt;

	d_CaT = Rush_Larsen(dt, d_CaT, dt_inf, tau_dt);
	double tau_ft = 1.0 / (0.0153 * exp(-(V + 61.7 + v_cat) / 83.3) + 0.015 * exp((V + 61.7 + v_cat) / 15.38));
	double ft_inf = 1.0 / (1.0 + exp((V + 61.7 + v_cat) / 5.6));
	// double ft_dot = (ft_inf - f_CaT) / tau_ft;
	f_CaT = Rush_Larsen(dt, f_CaT, ft_inf, tau_ft);

	ICaT = gCaT * f_CaT * d_CaT * (V - ecat);
	return ICaT;
}



