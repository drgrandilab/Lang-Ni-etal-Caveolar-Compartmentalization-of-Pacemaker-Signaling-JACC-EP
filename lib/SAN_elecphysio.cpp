
#include "SAN_elecphysio.hpp"

#include "ExplicitSolver.hpp"

SAN_elecphysio::SAN_elecphysio() {
	y[0] = 0.6530744974;
	y[1] = 0.5555484667;
	y[2] = 0.0017968124;
	y[3] = 0.3660355707;
	y[4] = 0.2532371183;
	y[5] = 0.9227779217;
	y[6] = 0.0067122334;
	y[7] = 0.0011943240;
	y[8] = 0.8922702053;
	y[9] = 0.0000381319;
	y[10] = 0.9419590387;
	y[11] = 0.7021037357;
	y[12] = 0.1123190249;
	y[13] = 0.3689412442;
	y[14] = 0.0206150195;
	y[15] = 0.4119481339;
	y[16] = 0.2219871208;
	y[17] = 0.0196667018;
	y[18] = 0.0108493621;
	y[19] = 0.5909220268;
	y[20] = 0.0061225524;
	y[21] = 0.7619838052;
	y[22] = 0.0000002414;
	y[23] = 0.0000000754;
	y[24] = 0.0068997684;
	y[25] = 0.1393274499;
	y[26] = 0.7603298944;
	y[27] = 0.0320868804;
	y[28] = 0.0150587275;
	y[29] = 0.0566119105;
	y[30] = 0.0000747469;
	y[31] = 0.0000347190;
	y[32] = 0.0533031555;
	y[33] = 1.6270005811;
	y[34] = 8.0372421686;
	y[35] = 139.8854603066;
	y[36] = -63.9188580935;


	for (int i = 0; i < 40; ++i)
	{
		ydot[i] = 0.0;  // set to zero;
	}

	double tmp [] = {0.951, 1.000, 1.260, 1.890, 1.000, 1.261, 1.023, 0.985, 1.463, 0.972, 1.244, 1.246, 1.082, 1.060, 1.223, 0.952, 1.265, 0.893}; // for sensitivity analysis
	par_SA.assign(tmp, tmp + 18);

	/*for (std::vector<double>::iterator i = par_SA.begin(); i != par_SA.end(); ++i)
	{
		std::cout<< *i << std::endl;
	}
	*/
	int static  p[] = {3, 0, 0, 0, 0, -40};
	//// Input parameters
	m_ind = p[1 - 1];
	// Model index: 0 for Kharche model, 1 for updated model (BPS 2017)
	Na_clamp =0;// p[2 - 1]; // 0 for free Na, 1 for Na clamp
	DB = p[3 - 1]; // 0 for control, 1 for DB
	OLD = p[4 - 1]; // 0 for control, 1 for OLD
	prot_index = p[5 - 1]; // 0 for no stimulation, 1 for voltage-clamp protocol, 2 for AP-clamp
	// 3 for no stimulation & NKA block at 10 s
	prot_input = p[6 - 1]; // for protocol parameter (V-clamp/NKA block)

	buffer_shannon = 0;

	if (buffer_shannon == 1) {
		ConcTC  = 0.070;
		kfTC    = 32.7;
		kbTC    = 0.0196;
		ConcTMC = 0.140;
		kfTMC   = 2.37;
		kbTMC   = 0.000032;
		kfTMM   = 0.003;
		kbTMM   = 0.00333;
		ConcCM  = 0.024;
		kfCM    = 34;
		kbCM    = 0.238;
		ConcCQ  = 2.7;
		kfCQ    = 100;
		kbCQ    = 65;
	}

	if (ISO == 1)  {//  see description in the paper (results are not exactly the same)
		vhalf_gh = 92;


		gcal12 = (0.0010 * 4.0 * 1.5) * 1.45;
		gcal13 = (0.0030 * 4.0 * 1.5) * 1.45;
		v_cal = 5;
		kmnap = 11;
		gkr = (0.8 * 0.002955) * 1.12;
		v_kr = 5;
		koca = 20;
		pumpkmf = pumpkmf / 2.0; // paper
	}



	if (m_ind == 3) {// -> Updated (CINC) with new If & OPTIMIZATION // BS 2017 poster
		gto = gto * 2.5;
		gsus = gsus * 2.5;
		gcal12 = gcal12 * 0;
		gcal13 = gcal13 * 2.0 / 3.0;
		v_cal = 10.0;
		slope_gh = 13.3;
		vhalf_gh = 114.4;
		// Optimization
		// opt_factors = [0.951,1.000,1.260,1.890,1.000,1.261,1.023,0.985,1.463,0.972,1.244,1.246,1.082,1.060,1.223,0.952,1.265,0.893]; // BS 2017
		// //opt_factors = [0.7185,0.6660,0.5137,5.2825,1.0000,2.2112,1.1820,1.0436,1.2067,0.8858,2.2686,1.6898,0.5639,0.8407,1.5052,0.3388,2.2326,0.7365]; // target v2
		// par_SA = par_SA.*opt_factors;
	}

	ih = 0.0;
	ina_ttxr = 0.0;
	ina_ttxs = 0.0;
	ical12 = 0.0;
	ical13 = 0.0;
	iks = 0.0;
	ikr = 0.0;
	ik1 = 0.0;
	ist = 0.0;
	ib = 0.0;
	icat = 0.0;
	inak = 0.0;
	isus = 0.0;
	inaca = 0.0;
	ito = 0.0;
	icap = 0.0;
}
SAN_elecphysio::~SAN_elecphysio() {
	;
}



int SAN_elecphysio::update_state_FE(double dt) {
	for (int i = 0; i < 40; ++i) {
		y[i] += ydot[i] * dt;
	}
	return 1;
}


int SAN_elecphysio::update_Na_and_K_currents(double t) {
	// State variables
	double dst       = y[1 - 1];
	double fst       = y[2 - 1];
	double dt        = y[3 - 1];
	double ft        = y[4 - 1];
	double ikr_act   = y[5 - 1];
	double ikr_inact = y[6 - 1];
	double iks_act   = y[7 - 1];
	double dl13      = y[8 - 1];
	double fl13      = y[9 - 1];
	double dl12      = y[10 - 1];
	double fl12      = y[11 - 1];
	double fca       = y[12 - 1];
	double m_ttxs = y[13 - 1];
	double h_ttxs = y[14 - 1];
	double j_ttxs = y[15 - 1];
	double m_ttxr = y[16 - 1];
	double h_ttxr = y[17 - 1];
	double j_ttxr = y[18 - 1];
	double y_1_2  = y[19 - 1];
	double q           = y[20 - 1];
	double r           = y[21 - 1];
	double resting     = y[22 - 1];
	double open        = y[23 - 1];
	double inactivated = y[24 - 1];
	double Ftc         = y[25 - 1];
	double Ftmc        = y[26 - 1];
	double Ftmm        = y[27 - 1];
	double Fcms        = y[28 - 1];
	double Fcmi        = y[29 - 1];
	double Fcq         = y[30 - 1];
	double casub       = y[31 - 1];
	double cai         = y[32 - 1];
	double carel       = y[33 - 1];
	double caup = y[34 - 1]; //+(1.3-y[34))*(t>55000-1];
	double nai = y[35 - 1];
	double ki = y[36 - 1];
	double v  = y[37 - 1];


	double ena = (R * T / F) * log(nao / nai);
	double ek  = (R * T / F) * log(ko / ki);
	double eks = (R * T / F) * log((ko + 0.12 * nao) / (ki + 0.12 * nai));


	double qa = 1.0 / (1.0 + exp(-(v + 67.0) / 5.0));
	double alphaqa = 1.0 / (0.15 * exp(-(v) / 11.0) + 0.2 * exp(-(v) / 700.0));
	double betaqa = 1.0 / (16.0 * exp((v) / 8.0) + 15.0 * exp((v) / 50.0));
	double tauqa = 1.0 / (alphaqa + betaqa);
	double alphaqi = 0.15 * 1.0 / (3100.0 * exp((v + 10.0) / 13.0) + 700.3 * exp((v + 10.0) / 70.0));
	double betaqi =  0.15 * 1.0 / (95.7 * exp(-(v + 10.0) / 10.0) + 50.0 * exp(-(v + 10.0) / 700.0)) + 0.000229 / (1 + exp(-(v + 10.0) / 5.0));
	double qi = alphaqi / (alphaqi + betaqi);
	double tauqi = 1.0 / (alphaqi + betaqi);
	double dst_dot = (qa - dst) / tauqa;
	double fst_dot = (qi - fst) / tauqi;
	ist = par_SA[1 - 1] * gst * dst * fst * (v - eist);

	//// INa - Na channel isoforms Nav1.1/Nav1.5 currents ***********************
	double fna = (9.52e-02 * exp(-6.3e-2 * (v + 34.4)) / (1 + 1.66 * exp(-0.225 * (v + 63.7)))) + 8.69e-2;
	double m3_inf_ttxr = 1.0 / (1.0 + exp(-(v + 45.213705) / 7.219547));
	double h_inf_ttxr = 1.0 / (1.0 + exp(-(v + 62.578120 ) / (-6.084036)));
	double m3_inf_ttxs = 1.0 / (1.0 + exp(-(v + 36.097331 - 5.0) / 5.0));
	double h_inf_ttxs = 1.0 / (1.0 + exp((v + 56.0) / 3.0));
	double m_inf_ttxr = pow(m3_inf_ttxr, 0.333);
	double m_inf_ttxs = pow(m3_inf_ttxs, 0.333);
	double tau_m = 1000.0 * ((0.6247e-03 / (0.832 * exp(-0.335 * (v + 56.7)) + 0.627 * exp(0.082 * (v + 65.01)))) + 0.0000492);
	double tau_h = 1000.0 * (((3.717e-06 * exp(-0.2815 * (v + 17.11))) / (1 + 0.003732 * exp(-0.3426 * (v + 37.76)))) + 0.0005977);
	double tau_j = 1000.0 * (((0.00000003186 * exp(-0.6219 * (v + 18.8))) / (1 + 0.00007189 * exp(-0.6683 * (v + 34.07)))) + 0.003556);
	double m_ttxs_dot = (m_inf_ttxs - m_ttxs) / tau_m;
	double h_ttxs_dot = (h_inf_ttxs - h_ttxs) / tau_h;
	double j_ttxs_dot = (h_inf_ttxs - j_ttxs) / tau_j;
	double hs = (1.0 - fna) * h_ttxs + fna * j_ttxs;
	double tau_mr = 1000.0 * ((0.6247e-03 / (0.832 * exp(-0.335 * (v + 56.7)) + 0.627 * exp(0.082 * (v + 65.01)))) + 0.0000492);
	double tau_hr = 1000.0 * (((3.717e-06 * exp(-0.2815 * (v + 17.11))) / (1 + 0.003732 * exp(-0.3426 * (v + 37.76)))) + 0.0005977);
	double tau_jr = 1000.0 * (((0.00000003186 * exp(-0.6219 * (v + 18.8))) / (1 + 0.00007189 * exp(-0.6683 * (v + 34.07)))) + 0.003556);
	double m_ttxr_dot = (m_inf_ttxr - m_ttxr) / tau_mr;
	double h_ttxr_dot = (h_inf_ttxr - h_ttxr) / tau_hr;
	double j_ttxr_dot = (h_inf_ttxr - j_ttxr) / tau_jr;
	double hsr = (1.0 - fna) * h_ttxr + fna * j_ttxr;
	double f_ina_ttxs = 0;
	double f_ina_ttxr = 0;
	if (fabs(v) > 0.005) {
		f_ina_ttxs = m_ttxs * m_ttxs * m_ttxs * hs * nao * (F * F / (R * T)) * ((exp((v - ena) * F / (R * T)) - 1.0) / (exp(v * F / (R * T)) - 1.0)) * v;
		f_ina_ttxr = m_ttxr * m_ttxr * m_ttxr * hsr * nao * (F * F / (R * T)) * ((exp((v - enattxr) * F / (R * T)) - 1.0) / (exp(v * F / (R * T)) - 1.0)) * v;
	} else {
		f_ina_ttxs = m_ttxs * m_ttxs * m_ttxs * hs * nao * F * ((exp((v - ena) * F / (R * T)) - 1.0));
		f_ina_ttxr = m_ttxr * m_ttxr * m_ttxr * hsr * nao * F * ((exp((v - enattxr) * F / (R * T)) - 1.0));
	}
	ina_ttxs = par_SA[2 - 1] * gna_ttxs * f_ina_ttxs;
	ina_ttxr = par_SA[3 - 1] * gna_ttxr * f_ina_ttxr;

	// //  If - Hyperpolarization-activated current *******************************
	double ih_pK = 0.6167; // ih_PNa = 0.3833;
	//y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/16.3));
	double y_inf = 1.0 / (1.0 + exp((v + vhalf_gh) / slope_gh));
	double tau_y_1_2 = 1 * (1.5049 / (exp(-(v + 590.3) * 0.01094) + exp((v - 85.1) / 17.2)));

	double y_1_2_dot = (y_inf - y_1_2) / tau_y_1_2;
	ihk  = If_scale * ih_pK * par_SA[7 - 1] * gh * y_1_2 * (v - ek);
	ihna = If_scale * (1 - ih_pK) * par_SA[7 - 1] * gh * y_1_2 * (v - ena);
	ih = (ihk + ihna);


	//  IK1 - Time-independent K current ***************************************
	double xk1inf = 1.0 / (1.0 + exp(0.070727 * (v - ek)));
	ik1 = par_SA[-1 + 8] * gk1 * xk1inf * (ko / (ko + 0.228880)) * (v - ek);

	//  IKr - Rapid delayed rectifying K current *******************************
	double ikr_act_inf = 1.0 / (1.0 + exp(-(v + 21.173694 + v_kr) / 9.757086));
	double tau_ikr_act = 0.699821 / (0.003596 * exp((v) / 15.339290) + 0.000177 * exp(-(v) / 25.868423));
	double ikr_act_dot = (ikr_act_inf - ikr_act) / tau_ikr_act;
	double ikr_inact_inf = 1.0 / (1.0 + exp((v + 20.758474 - 4.0) / (19.0)));
	double tau_ikr_inact = 0.2 + 0.9 * 1.0 / (0.1 * exp(v / 54.645) + 0.656 * exp(v / 106.157));
	double ikr_inact_dot = (ikr_inact_inf - ikr_inact) / tau_ikr_inact;
	ikr = par_SA[-1 + 9] * gkr * ikr_act * ikr_inact * (v - ek);

	//  IKs - Slow delayed rectifying K current ********************************
	double iks_act_inf = 1.0 / (1.0 + exp(-(v - 20.876040) / 11.852723));
	double tau_iks_act =  1000.0 / (13.097938 / (1.0 + exp(-(v - 48.910584) / 10.630272)) + exp(-(v) / 35.316539));
	double iks_act_dot = (iks_act_inf - iks_act) / tau_iks_act;
	iks = par_SA[-1 + 10] * gks * iks_act * iks_act * (v - eks);

	//  Ito - Transient component of 4-AP-sensitive currents *******************
	double q_inf = 1.0 / (1.0 + exp((v + 49.0) / 13.0));
	double tau_q = (6.06 + 39.102 / (0.57 * exp(-0.08 * (v + 44.0)) + 0.065 * exp(0.1 * (v + 45.93)))) / 0.67;
	double q_dot = (q_inf - q) / tau_q;
	double r_inf = 1.0 / (1.0 + exp(-(v - 19.3) / 15.0));
	double tau_r = (2.75 + 14.40516 / (1.037 * exp(0.09 * (v + 30.61)) + 0.369 * exp(-0.12 * (v + 23.84)))) / 0.303;
	double r_dot = (r_inf - r) / tau_r;
	ito = par_SA[-1 + 11] * gto * q * r * (v - ek);

	//  Isus - Sustained component of 4-AP-sensitive currents ******************
	isus = par_SA[-1 + 12] * gsus * r * (v - ek);

	//  Ib - Background Na, Ca and K currents **********************************
	ibna = par_SA[-1 + 13] * gbna * (v - ena);
	ibk = gbk * (v - ek);

	//    INaK - Na-K pump current ***********************************************
	inak = par_SA[15 - 1] * inakmax * ((pow(ko, 1.2)) / (pow(kmkp, 1.2) + pow(ko, 1.2))) * (pow(nai, 1.3) / (pow(kmnap, 1.3) + pow(nai, 1.3))) / (1.0 + exp(-(v - ena + 120.0) / 30.0));


	//     Jrel = par_SA(17)*ks*open*(caup - casub);
	//     Fcq_dot = kfCQ*caup*(1.0-Fcq)-kbCQ*Fcq; // SRcomp
	//     casub_dot = ((-ca_flux+Jrel*vup)/vsub-Jcadif - ConcCM*Fcms_dot);  // SRcomp
	//     cai_dot = ((Jcadif*vsub-Jup*vup)/vi - (ConcCM*Fcmi_dot + ConcTC*Ftc_dot + ConcTMC*Ftmc_dot)); // SRcomp
	//     carel_dot = 0; // SRcomp
	//     caup_dot = (Jup  - Jrel - ConcCQ*Fcq_dot); // SRcomp


	// double I_app = 0.0;

	// v_dot = - total_current/capacitance;

	//  Output
	ydot[-1 + 1] = dst_dot;
	ydot[-1 + 2] = fst_dot;
	// ydot[-1 + 3] = dt_dot;
	// ydot[-1 + 4] = ft_dot;
	ydot[-1 + 5] = ikr_act_dot;
	ydot[-1 + 6] = ikr_inact_dot;
	ydot[-1 + 7] = iks_act_dot;
	// ydot[-1 + 8] = dl13_dot;
	// ydot[-1 + 9] = fl13_dot;
	// ydot[-1 + 10] = dl12_dot;
	// ydot[-1 + 11] = fl12_dot;
	// ydot[-1 + 12] = fca_dot;
	ydot[-1 + 13] = m_ttxs_dot;
	ydot[-1 + 14] = h_ttxs_dot;
	ydot[-1 + 15] = j_ttxs_dot;
	ydot[-1 + 16] = m_ttxr_dot;
	ydot[-1 + 17] = h_ttxr_dot;
	ydot[-1 + 18] = j_ttxr_dot;
	ydot[-1 + 19] = y_1_2_dot;
	ydot[-1 + 20] = q_dot;
	ydot[-1 + 21] = r_dot;
	// ydot[-1 + 22] = resting_dot;
	// ydot[-1 + 23] = open_dot;
	// ydot[-1 + 24] = inactivated_dot;
	// ydot[-1 + 25] = Ftc_dot;
	// ydot[-1 + 26] = Ftmc_dot;
	// ydot[-1 + 27] = Ftmm_dot;
	// ydot[-1 + 28] = Fcms_dot;
	// ydot[-1 + 29] = Fcmi_dot;
	// ydot[-1 + 30] = Fcq_dot;

	// // ydot[-1 + 31] = casub_dot; //*(t<4860); // MODIFIED
	// ydot[-1 + 32] = cai_dot * (1 - ((prot_index == 3) | (prot_index == 4))); //*(t<5000);//(t<20000);  // MODIFIED
	// ydot[-1 + 33] = carel_dot; // MODIFIED
	// ydot[-1 + 34] = caup_dot; //*(t<55000);  // MODIFIED
	// ydot[-1 + 35] = (1 - Na_clamp) * nai_dot;
	// ydot[-1 + 36] = 0 * ki_dot;
	// ydot[-1 + 37] = v_dot;
	return 1;
}



int SAN_elecphysio::update_Na_K_concentration(double t) {

	double nai_tot = ihna + ina_ttxr + ina_ttxs + 3.0 * inak + 3.0 * inaca + ist + ibna;
	double nai_dot = (-nai_tot) / (F * vi);

	double ki_tot = ihk + iks + ikr + ik1 + ibk - 2.0 * inak + isus + ito;
	double ki_dot = (-ki_tot) / (F * vi);


	// // ydot[-1 + 31] = casub_dot; //*(t<4860); // MODIFIED
	// ydot[-1 + 32] = cai_dot * (1 - ((prot_index == 3) | (prot_index == 4))); //*(t<5000);//(t<20000);  // MODIFIED
	// ydot[-1 + 33] = carel_dot; // MODIFIED
	// ydot[-1 + 34] = caup_dot; //*(t<55000);  // MODIFIED
	ydot[-1 + 35] = (1 - Na_clamp) * nai_dot;
	ydot[-1 + 36] = 0 * ki_dot;  // K+ clamped. 



	// ydot[-1 + 37] = v_dot;
	return 1;
}


int SAN_elecphysio::com_total_current(double t) {

	double casub       = y[31 - 1];
	double v = y[37-1];

	ib = (ibna + ibca + ibk);
	double total_current = ih + ina_ttxr + ina_ttxs + ical12 + ical13 + iks + ikr + ik1 + ist + ib + icat + inak + isus + inaca + ito + icap;
	double v_dot = - (total_current) / capacitance;
	ydot[-1 + 37] = v_dot;
	return 1;
}




