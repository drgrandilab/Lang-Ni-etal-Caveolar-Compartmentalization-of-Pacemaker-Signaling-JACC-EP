


#include "subcell.hpp"


double CSubcell::get_jSR_inst_buffering(double CaSR) {
	double rho = rhoinf * pow(double(CaSR / 1000.0), hh) / (pow(KK / 1000, hh) + pow(double(CaSR / 1000), hh));
	if (rho < 0.0000001)rho = 0.0000001; /// to avoid div small (for cjsr<200)
	double MM = (sqrt(1 + 8 * rho * BCSQN) - 1) / (4 * rho * BCSQN);
	double ncjsr = MM * nM + (1 - MM) * nD;

	//      double rhopri=(hh*rhoinf*pow(double(CaSR),double(hh-1))*(pow(double(KK),double(hh))+pow(double(CaSR),double(hh)))-rhoinf*pow(double(CaSR),double(hh))*hh*pow(double(CaSR),double(hh-1)))/pow((pow(double(KK),double(hh))+pow(double(CaSR),double(hh))),2);

	//      double rhopri=(hh*rhoinf*pow(double(CaSR/1000.0),double(hh-1))*(pow(double(KK/1000.0),double(hh))+pow(double(CaSR/1000.0),double(hh)))-rhoinf*pow(double(CaSR/1000.0),double(hh))*hh*pow(double(CaSR/1000.0),double(hh-1)))/pow((pow(double(KK/1000.0),double(hh))+pow(double(CaSR/1000.0),(hh))),2)/1000.0;

	//gpu code
	double rhopri = (hh * rhoinf * pow(double(CaSR / 1000.0), double(hh - 1)) * (pow(double(KK / 1000.0), double(hh)) + pow(double(CaSR / 1000.0), double(hh))) - rhoinf * pow(double(CaSR / 1000.0), double(hh)) * hh * pow(double(CaSR / 1000.0), double(hh - 1))) / pow((pow(double(KK / 1000.0), double(hh)) + pow(double(CaSR / 1000.0), double(hh))), 2);
	rhopri *= 0.001;


	double dMdc = (((1.0 / 2.0) * pow(double(1 + 8 * rho * BCSQN), double(-1.0 / 2.0)) * 8 * rhopri * BCSQN * 4 * rho * BCSQN) - (sqrt(1 + 8 * rho * BCSQN) - 1) * 4 * rhopri * BCSQN) / pow(4 * rho * BCSQN, 2);
	double dndc = dMdc * (nM - nD);

	double Betajsr = 1 / (1 + (Kc * BCSQN * ncjsr + dndc * (CaSR * Kc + CaSR * CaSR)) / ((Kc + CaSR) * (Kc + CaSR)));

	return Betajsr;
}



double CSubcell::get_cleft_caj_inst_buffering(double Caj) {

	//Instantaneous buffering functions
	const double KCAM = 7.0;
	const double BCAM = 24.0;
	const double KSR = 0.6;
	const double BSR = 47.0;
	const double KMCa = 0.033;
	const double BMCa = 140.0;
	const double KMMg = 3.64;
	const double BMMg = 140.0;
	const double KSLH = 0.3;
	const double BSLH = 13.4;

	double CAM = BCAM * KCAM / ((Caj + KCAM) * (Caj + KCAM));
	double SR = BSR * KSR / ((Caj + KSR) * (Caj + KSR));
	double MCa = BMCa * KMCa / ((Caj + KMCa) * (Caj + KMCa));
	double MMg = BMMg * KMMg / ((Caj + KMMg) * (Caj + KMMg));
	double SLH = BSLH * KSLH / ((Caj + KSLH) * (Caj + KSLH)); // only for cs
	// double Betap = 1 / (1 + CAM + SR + MCa + MMg + SLH);
	double Betap = 1 / (1 + CAM /*+ SR + MCa + MMg*/ + SLH);

	return Betap;

}


double CSubcell::get_submem_casl_inst_buffering(double Casl) {

	//Instantaneous buffering functions
	const double KCAM = 7.0;
	const double BCAM = 24.0;
	const double KSR = 0.6;
	const double BSR = 47.0;
	const double KMCa = 0.033;
	const double BMCa = 140.0;
	const double KMMg = 3.64;
	const double BMMg = 140.0;
	const double KSLH = 0.3;
	const double BSLH = 13.4;

	double   CAM = BCAM * KCAM / ((Casl + KCAM) * (Casl + KCAM));
	double  SR   = BSR * KSR / ((Casl + KSR) * (Casl + KSR));
	double   MCa = BMCa * KMCa / ((Casl + KMCa) * (Casl + KMCa));
	double   MMg = BMMg * KMMg / ((Casl + KMMg) * (Casl + KMMg));
	double SLH   = BSLH * KSLH / ((Casl + KSLH) * (Casl + KSLH)); // only for cs
	// double Betas = 1 / (1 + CAM + SR + MCa + MMg + SLH);
	double Betas = 1 / (1 + CAM /*+ SR + MCa + MMg */+ SLH);

	return Betas;

}




double CSubcell::get_cytosol_cai_inst_buffering(double Cai) {
	//Instantaneous buffering functions
	const double KCAM = 7.0;
	const double BCAM = 24.0;
	const double KSR = 0.6;
	const double BSR = 47.0;
	const double KMCa = 0.033;
	const double BMCa = 140.0;
	const double KMMg = 3.64;
	const double BMMg = 140.0;
	const double KSLH = 0.3;
	const double BSLH = 13.4;

	//Troponin C dynamic buffering current ITCi and ITCs


	double CAM = BCAM * KCAM / ((Cai + KCAM) * (Cai + KCAM));
	double SR = BSR * KSR / ((Cai + KSR) * (Cai + KSR));
	double MCa = BMCa * KMCa / ((Cai + KMCa) * (Cai + KMCa));
	double MMg = BMMg * KMMg / ((Cai + KMMg) * (Cai + KMMg));
	// double Betai = 1 / (1 + CAM + SR + MCa /*+ MMg*/);  // MMg for Mgi only  removing MMg would increase Cai
	double Betai = 1 / (1 + CAM + SR + MCa + MMg);

	return Betai;
}


double CSubcell::calculate_dynamic_buffer_cytosol(int id, double Cai, double dt) {
	double Bmax_TnClow = 70.0e-3;
	int ISO=0;
	double koff_tncl =  (1.0 + 0.5 * ISO) * 19.6e-3;  // according to the matlab code
	double kon_tncl = 32.7;
	double kon_tnchca  = 2.37;
	double koff_tnchca = 0.032e-3;
	double kon_tnchmg  = 3.0e-3;
	double koff_tnchmg = 3.33e-3;
	double kon_cam     = 34.0;
	double Bmax_CaM    = 24.0e-3;
	double koff_cam    = 238.0e-3;
	double kon_myoca   = 13.8;
	double Bmax_myosin = 140.0e-3;
	double koff_myoca  = 0.46e-3;
	double kon_myomg   = 0.0157;
	double koff_myomg  = 0.057e-3;
	double kon_sr      = 100.0;
	double Bmax_SR     = 19.0 * 0.9e-3;
	double koff_sr     = 60.0e-3;
	const double Mgi = 1.;    //  Intracellular Mg  [mM]


	//  Cytosolic Ca Buffers
	double df_TC  = kon_tncl * Cai* (Bmax_TnClow - Tropc_vec[id] ) - koff_tncl * Tropc_vec[id] ;      //  TnCL      [mM/ms]
	// ydot[19]  = kon_tnchca * Cai* (Bmax_TnChigh - y[19] - y[20] ) - koff_tnchca * y[19] ; //  TnCHc     [mM/ms]
	// ydot[20]  = kon_tnchmg * Mgi * (Bmax_TnChigh - y[19] - y[20] ) - koff_tnchmg * y[20] ; //  TnCHm     [mM/ms]
	double df_CMi  = kon_cam * Cai* (Bmax_CaM - CaM_vec[id]) - koff_cam * CaM_vec[id];           //  CaM       [mM/ms]
	double df_TMC  = kon_myoca * Cai* (Bmax_myosin - Myosin_Ca_vec[id] - Myosin_Mg_vec[id] ) - koff_myoca * Myosin_Ca_vec[id] ; //  Myosin_ca [mM/ms]
	double df_TMM  = kon_myomg * Mgi * (Bmax_myosin - Myosin_Ca_vec[id] - Myosin_Mg_vec[id] ) - koff_myomg * Myosin_Mg_vec[id] ; //  Myosin_mg [mM/ms]
	double df_SRb  = kon_sr * Cai* (Bmax_SR - SRB_vec[id]) - koff_sr * SRB_vec[id];              //  SRB       [mM/ms]

	Tropc_vec[id]     += df_TC*dt;
	CaM_vec[id]       += df_CMi*dt;
	Myosin_Ca_vec[id] += df_TMC*dt;
	Myosin_Mg_vec[id] += df_TMM*dt;
	SRB_vec[id]       += df_SRb*dt;

	double J_CaB_cytosol = df_TC + df_CMi + df_TMC + /*df_TMM + */df_SRb;//ydot[18] + ydot[19] + ydot[20] + ydot[21] + ydot[22] + ydot[23] + ydot[24]; // sum(ydot(19:25));

	return J_CaB_cytosol;
}

// double CMtot=0.045e3;
// double TCtot=0.031e3;
// double TMCtot=0.062e3;


// df_TC = kf_TC*Cai *(1-f_TC) - kb_TC*f_TC;
// df_TMC = kf_TMC*Cai *(1-(f_TMC+f_TMM)) - kb_TMC*f_TMC;
// df_TMM = kf_TMM*Mgi *(1-(f_TMC+f_TMM)) - kb_TMM*f_TMM;
// df_CMi = kf_CM*Cai *(1-f_CMi) - kb_CM*f_CMi;


// (CMtot*df_CMi+ TCtot*df_TC + TMCtot*df_TMC)


