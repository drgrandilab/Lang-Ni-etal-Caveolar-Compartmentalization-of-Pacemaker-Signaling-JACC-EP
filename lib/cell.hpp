#ifndef ___CELL_H
#define ___CELL_H

using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

class CCell{
private:
  void pacex(double stim=0);
  static const int N=101;
  static const double vc;
  static const double stim;
  static const double stimduration;
  
  
  
  void comp_na_buffer(void);
  double comp_ca_buffer_J_CaB_cytosol(void);
  double comp_ca_buffer_J_CaB_junction(void);
  double comp_ca_buffer_J_CaB_sl(void);

  double comp_I_Ca_junc(void);
  double comp_I_Ca_sl(void);
  double comp_I_CaK(void);
  double comp_I_CaNa_junc(void);
  double comp_I_CaNa_sl(void);

  
  double comp_sr_release(void);
  double comp_serca(void);
  double comp_sr_leak(void);


  double comp_ina(double *I_Na_junc, double *I_Na_sl);
  double comp_I_nabk_junc(void);
  double comp_I_nabk_sl(void);
  double comp_inak_junc(void);
  double comp_inak_sl(void);


  double comp_ito(void);
  double comp_ikr(void);
  double comp_iks(void);
  double comp_ikp(void);
  double comp_ikach(void);
  double comp_ik1(void);
  double comp_ikur(void);
  
  
  double comp_I_pca_junc(void);
  double comp_I_pca_sl(void);
  double comp_I_cabk_junc(void);
  double comp_I_cabk_sl(void);
  
  double comp_incx_junc(void);
  double comp_incx_sl(void);


  double comp_I_ClCa_junc(void);
  double comp_I_ClCa_sl(void);
  double comp_I_Clbk(void);
  double comp_I_ClCFTR(void);

  // Nernst Potentials
  double ena_junc;     // [mV]
  double ena_sl;       // [mV]
  double ek;         // [mV]
  double eca_junc;   // [mV]
  double eca_sl;     // [mV]
  double ecl;            // [mV]
  
  double dt,ddt;
  
  double vold;

public:
  void pace(double stim=0);
  double setdt(double dtt){dt=dtt;return dt;}
  double getdt(void){return dt;}
  int getdim(void){return N;}
  double getvc(void){return vc;}
  double getstim(void){return stim;}
  double getstimduration(void){return stimduration;}
  CCell(void);
  virtual ~CCell();
  CCell& operator=(const CCell& cell);
  double *y;
  double *ydot;
  double &v,&ci,&cnsr;

  double getvmax(void){return 10;}
  double getvmin(void){return -80;}

  void prepare(double pcl=500, int iter=0);

  //conditions

  // ISO administration (0 or 1)
  double ISO;

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


  double _I_to;
  double _I_kr;
  double _I_ks;
  double _I_ki;
  double _I_kp;
  double _I_kur;
  double _I_KAch;
  double _I_CaK;

  double _I_Catot;
  double _I_ncx;
  
  
  double gnabar;//Na Current conductance
  double gnabbar;//Na Background Current conductance
  double gnakbar;//Na/K Pump Current conductance
  double gkrbar;//Rapidly Activating K Current conductance
  double gksbar;//Slowly Activating K Current conductance
  double gkpbar;//Plateau K Current conductance
  double gkachbar;//Muscarinic-Receptor-Activated K Current conductance
  double gtobar;//Transient Outward K Current conductance
  double gkurbar;//Ultra rapid delayed rectifier Outward K Current conductance
  double gkibar;//Time-Independent K Current conductance
  double vupbar;//SERCA strength
  double gcabbar;//Ca Background conductance
  double gpcabar;//Sarcolemmal Ca Pump Current conductance
  double gncxbar;//Na/Ca Exchanger flux conductance
  double gcabar;//L-type Calcium Current conductance

  double tauff;//L-type Ca channel recovery time constant
  double taujj;//Na channel recovery time constant
};

#endif /* ___CELL_H */
