

#ifndef SAN_ELECPHYSIO_HPP
#define SAN_ELECPHYSIO_HPP

using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
// %% Model Parameters - Kharche et al model

// R = 8.314472; % [J/mol*K]
// T = 310.5; % [K] 
// F = 96.4845; % [C/kmol]
// FRT = F/(R*T);

// capacitance = 0.025; % [nF] 0.032 in kurata model
// %vcell=3.0; % [pL]
// %l_cell=66.3767257; % [um]
// %r_cell=3.792956; % [um]
// vsub = 0.03328117;
// vi = 1.34671883;
// vrel = 0.0036;
// vup = 0.0348;

// % vcell = 3;
// % vsub = 0.01*vcell; % 0.03
// % vi = 0.46*vcell-vsub; % 1.35
// % vrel = 0.0012*vcell; % 0.0036
// % vup = 0.0116*vcell; % 0.0348

// Mgi = 2.5; % [mM]
// nao = 140; % [mM]
// cao = 1.8; % [mM]
// ko = 5.4; % [mM]



class SAN_elecphysio
{
public:
	SAN_elecphysio();
	~SAN_elecphysio();
	void pacex(double stim = 0);
	constexpr static int N = 101;
	 double vc;
	 double stim;
	 double stimduration;

	int update_single_time_step(double);
	int update_state_FE(double);
	int update_Na_and_K_currents(double);
	int update_Ca_currents(double);
	int com_INaCa(double);
	int com_SR_flux(double);
	int com_buffer(double);
	int update_Ca_concentration(double);
	int update_Na_K_concentration(double);
	int com_total_current(double);
	int update(double);


	double y[40];
	double ydot[40];


	//conditions


	////Model Flags
	// Model for INa
	double flagMina; // w/ 1 Markov INa model, w/ 0 H&H formulation for fast and late INa

	double drug_index;
	double drug;

	// AF
	double AF;
	// Right ATRIUM
	double RA;
	// other
	double epi; // EPI or ENDO?
	double bGal;
	double EAD;

	double Ach;//Acetylcholine concentration (uM)
 	// ih + ina_ttxr + ina_ttxs + ical12 + ical13 + iks + ikr + ik1 + ist + ib + icat + inak + isus + inaca + ito

	double ito;
	double ikr;
	double iks;
	double ik1;
	double isus, ibna, ibca, ibk, ib;
	double icap;
	double ihk, ihna, ih;  //funny current
	double ist; // I've no idea what this current is... at the moment;// 12:16:54, Wed, 05-December-2018, By Haibo
	double ina_ttxs, ina_ttxr;
	double icat;
	double ical12, ical13;
	double inak, inaca;
	double Jrel, Jup, Jtr;




	// %% Model Parameters - Kharche et al model

	double constexpr static R = 8.314472; //  [J/mol*K]
	double constexpr static T = 310.5; //  [K] 
	double constexpr static F = 96.4845; //  [C/kmol]
	double constexpr static FRT = F/(R*T);

	double constexpr static capacitance = 0.025; //  [nF] 0.032 in kurata model
	// vcell=3.0; //  [pL]
	// l_cell=66.3767257; //  [um]
	// r_cell=3.792956; //  [um]
	double constexpr static vsub = 0.03328117;
	double constexpr static vi = 1.34671883;
	double constexpr static vrel = 0.0036;
	double constexpr static vup = 0.0348;

	//  vcell = 3;
	//  vsub = 0.01*vcell; //  0.03
	//  vi = 0.46*vcell-vsub; //  1.35
	//  vrel = 0.0012*vcell; //  0.0036
	//  vup = 0.0116*vcell; //  0.0348

	double constexpr static Mgi = 2.5; //  [mM]
	double constexpr static nao = 140; //  [mM]
	double constexpr static cao = 1.8; //  [mM]
	double constexpr static ko = 5.4; //  [mM]


	std::vector<double> par_SA;

	std::vector<int> p;


	double gst = 0.00006;
	double eist = 17.0;
	double gbna = 0.0001215;
	double gbca = 0.000015;
	double gbk = 0.0000025;
	double gk1 = 0.229 * 0.0039228 * 0.9;
	double gks = 0.000299;
	double ecal = 47.0;
	double kmfca = 0.00035;
	double alpha_fca = 0.021;
	// all_ica_multiplier=1.0;
	double ecat = 45.0;
	double enattxr = 41.5761;
	// multiplier2=1.0;
	double gsus = 0.00039060;
	double inakmax_multiplier = 1;
	// if prot_index == 9
	//     NKA_block = prot_input; // 0.50 or 0.75//
	//     inakmax_multiplier = (1-(t>10e3)*NKA_block); // NKA block after 10 s
	//     //inakmax_multiplier = (1-(t>10e3)*(t<40e3)*0.50); // transient NKA block
	// end
	double inakmax = inakmax_multiplier * 1.85 * 0.077;
	double kmnap = 14.0;
	double kmkp = 1.4;
	double K1ni = 395.3;
	double K1no = 1628;
	double K2ni = 2.289;
	double K2no = 561.4;
	double K3ni = 26.44;
	double K3no = 4.663;
	double Kci = 0.0207;
	double Kco = 3.663;
	double Kcni = 26.44;
	double Qci = 0.1369;
	double Qco = 0.0;
	double Qn = 0.4315;
	double tdifca = 0.04;
	double Ttr = 40.0;

	//  Buffer
	double ConcTC = 0.031;
	double kfTC = 88.8;
	double kbTC = 0.446;
	double ConcTMC = 0.062;
	double kfTMC = 237.7;
	double kbTMC = 0.00751;
	double kfTMM = 2.277;
	double kbTMM = 0.751;
	double ConcCM = 0.045;
	double kfCM = 237.7;
	double kbCM = 0.542;
	double ConcCQ = 10.0;
	double kfCQ = 0.534;
	double kbCQ = 0.445;



	double koca = 10.0;
	double kom = 0.06;
	double kica = 0.5;
	double kim = 0.005;
	double eca50sr = 0.45;
	double maxsr = 15.0;
	double minsr = 1.0;
	double hsrr = 2.5;
	double pumphill = 2.0;

	double gna_ttxs = 0.1 * 5.925e-05;
	double gna_ttxr = 0.1 * 5.925e-05;
	double gcal12 = 0.0010 * 4.0 * 1.5;
	double gcal13 = 0.0030 * 4.0 * 1.5;
	double gcat = 0.75 * 0.01862;
	double gh = 0.0057;
	double vhalf_gh = 106.8;
	double gkr = 1 * (0.8 * 0.002955);
	double gto = 0.00492;
	double kNaCa = 5.5; //  *(t < 4860);
	double Pup = 0.04;
	double ks = 1300000; //  *(t < 4860);
	//  pumpkmf = 0.000246; //  paper
	//  pumpkmr = 3.29; //  paper
	double pumpkmf = 0.00008; //  c code
	double pumpkmr = 4.5; //  c code

	double slope_gh = 16.3;
	double tauS_DB_ca = 1;
	double tauF_DB_ca = 1;
	double v_cal = 0;
	double v_cat = 0;
	double vDB = 0;
	double v_kr = 0;

	int ISO = 0;

		//// Input parameters
	int m_ind = 3;
	// Model index: 0 for Kharche model, 1 for updated model (BPS 2017)
	int Na_clamp = 0; // 0 for free Na, 1 for Na clamp
	int DB = 0; // 0 for control, 1 for DB
	int OLD = 0; // 0 for control, 1 for OLD
	int prot_index = 0; // 0 for no stimulation, 1 for voltage-clamp protocol, 2 for AP-clamp
	// 3 for no stimulation & NKA block at 10 s
	int prot_input = 0; // for protocol parameter (V-clamp/NKA block)

	int buffer_shannon = 0;


	double If_scale = 1.0;

};



#endif



	// double _I_kp;
	// double _I_kur;
	// double _I_KAch;
	// double _I_CaK;

	// double _I_Catot;
	// double _I_ncx;


	// double gnabar;//Na Current conductance
	// double gnabbar;//Na Background Current conductance
	// double gnakbar;//Na/K Pump Current conductance
	// double gkrbar;//Rapidly Activating K Current conductance
	// double gksbar;//Slowly Activating K Current conductance
	// double gkpbar;//Plateau K Current conductance
	// double gkachbar;//Muscarinic-Receptor-Activated K Current conductance
	// double gtobar;//Transient Outward K Current conductance
	// double gkurbar;//Ultra rapid delayed rectifier Outward K Current conductance
	// double gkibar;//Time-Independent K Current conductance
	// double vupbar;//SERCA strength
	// double gcabbar;//Ca Background conductance
	// double gpcabar;//Sarcolemmal Ca Pump Current conductance
	// double gncxbar;//Na/Ca Exchanger flux conductance
	// double gcabar;//L-type Calcium Current conductance

	// double tauff;//L-type Ca channel recovery time constant
	// double taujj;//Na channel recovery time constant
	// // int update__currents(double);



	// void comp_na_buffer(void);
	// double comp_ca_buffer_J_CaB_cytosol(void);
	// double comp_ca_buffer_J_CaB_junction(void);
	// double comp_ca_buffer_J_CaB_sl(void);

	// double comp_I_Ca_junc(void);
	// double comp_I_Ca_sl(void);
	// double comp_I_CaK(void);
	// double comp_I_CaNa_junc(void);
	// double comp_I_CaNa_sl(void);


	// double comp_sr_release(void);
	// double comp_serca(void);
	// double comp_sr_leak(void);


	// double comp_ina(double *I_Na_junc, double *I_Na_sl);
	// double comp_I_nabk_junc(void);
	// double comp_I_nabk_sl(void);
	// double comp_inak_junc(void);
	// double comp_inak_sl(void);


	// double comp_ito(void);
	// double comp_ikr(void);
	// double comp_iks(void);
	// double comp_ikp(void);
	// double comp_ikach(void);
	// double comp_ik1(void);
	// double comp_ikur(void);


	// double comp_I_pca_junc(void);
	// double comp_I_pca_sl(void);
	// double comp_I_cabk_junc(void);
	// double comp_I_cabk_sl(void);

	// double comp_incx_junc(void);
	// double comp_incx_sl(void);


	// double comp_I_ClCa_junc(void);
	// double comp_I_ClCa_sl(void);
	// double comp_I_Clbk(void);
	// double comp_I_ClCFTR(void);

	// // Nernst Potentials
	// double ena_junc;     // [mV]
	// double ena_sl;       // [mV]
	// double ek;         // [mV]
	// double eca_junc;   // [mV]
	// double eca_sl;     // [mV]
	// double ecl;            // [mV]

	// double dt, ddt;

	// double vold;

	// void pace(double stim = 0);
	// double setdt(double dtt) {dt = dtt; return dt;}
	// double getdt(void) {return dt;}
	// int getdim(void) {return N;}
	// double getvc(void) {return vc;}
	// double getstim(void) {return stim;}
	// double getstimduration(void) {return stimduration;}
	// // double *y;
	// // double &v, &ci, &cnsr;
	// double getvmax(void) {return 10;}
	// double getvmin(void) {return -80;}

	// void prepare(double pcl = 500, int iter = 0);