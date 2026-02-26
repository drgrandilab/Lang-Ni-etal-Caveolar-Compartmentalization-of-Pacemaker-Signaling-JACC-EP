#include "cell.hpp"
#include "const.hpp"


const double CCell::vc=-80;
const double CCell::stim=12.5;
const double CCell::stimduration=5;



CCell::CCell(void) : y(new double[N]), ydot(new double[N]),v(y[39]), ci(y[38]), cnsr(y[31]) {
// initial conditions
  dt=0.1;
  vold = -80;
  double tmp[] = {0,
0.0202919327386704,0.936490564532038,0.936304480468966,5.42243465425516e-06,0.999659741420995,0.0387439718487863,0.0286342994782812,0.00405157400000000,0.994551100000000,0.000648056012992726,0.973061283727168,0.00193864192924699,0.00419293036123511,0.809901538287505,1.32068131238543e-06,3.09984986234518e-07,3.76351356810290,0.821291689576353,0.0163997256637698,0.124384775412331,0.00737380343482592,0.000611244439499118,0.00302488782583208,0.136472047162854,0.00399721388615559,0.0117402146171874,0.0188546757489768,0.0962772087252904,0.176838896529319,1.05807454671198,0.446032394313757,9.90313914499844,9.90459022726698,9.90483249546865,120,0.000281071550870317,0.000204821658962907,0.000182893525702204,-81.7495598043969,0.994600000000000,1,0.00150000000000000,0.0244000000000000,0.149400000000000,0.407100000000000,0.416100000000000,0,0.000100000000000000,0.000600000000000000,0.000800000000000000,0,0,0,0,0,0,0,0.000149496595416145,0.988960482230953,0.00264830219853136,0.173614838476469,-2164.01487897279,0.802586369629388,0.0286375258151476,0.000416030212867675,1.16241361165146e-06,8.11698968573295e-07,2.89626479052790e-08,4.20742890084225e-10,1.17510183512479e-12,0.105921364476028,0.00378635125822485,6.03834795404533e-05,0.0585360319738570,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (int i=0;i<N;i++) {
    y[i]=tmp[i];
  }

//  gna=


  // ISO administration (0 or 1)
  ISO = 0;

  ////Model Flags
  // Model for INa
  flagMina = 1; // w/ 1 Markov INa model, w/ 0 H&H formulation for fast and late INa
  // AF
  AF = 0;
  // Right ATRIUM
  RA = 0;
  // other
  epi = 1; // EPI or ENDO?
  bGal = 1;
  EAD = 0;

  drug_index=0;
  drug=0;

  Ach=0.1;//Acetylcholine concentration (uM)



  gnabar=1;//Na Current conductance
  gnabbar=1;//Na Background Current conductance
  gnakbar=1;//Na/K Pump Current conductance
  gkrbar=1;//Rapidly Activating K Current conductance
  gksbar=1;//Slowly Activating K Current conductance
  gkpbar=1;//Plateau K Current conductance
  gkachbar=1;//Muscarinic-Receptor-Activated K Current conductance
  gtobar=1;//Transient Outward K Current conductance
  gkurbar=1;//Ultra rapid delayed rectifier Outward K Current conductance
  gkibar=1;//Time-Independent K Current conductance
  vupbar=1;//SERCA strength
  gcabbar=1;//Ca Background conductance
  gpcabar=1;//Sarcolemmal Ca Pump Current conductance
  gncxbar=1;//Na/Ca Exchanger flux conductance
  gcabar=1;//L-type Calcium Current conductance
  
  tauff=1;
  taujj=1;


}
CCell::~CCell() {
    delete[] y;
    delete[] ydot;
}
void CCell::prepare(double pcl, int iter) {
  int Tn=pcl*iter/dt, bcln=pcl/dt, durn=stimduration/dt;
  for (int tn=0;tn<Tn;tn++)
  {
    if (tn%bcln < durn)
      pace(stim);
    else
      pace();
  }
}
CCell& CCell::operator=(const CCell& cell) {
  if (&cell!=this)
  {
    for (int i=0;i<N;i++)
    {
      y[i]=cell.y[i];
    }
    vold=cell.vold;
    dt=cell.dt;
    gnabar=cell.gnabar;//Na Current conductance
    gnabbar=cell.gnabbar;//Na Background Current conductance
    gnakbar=cell.gnakbar;//Na/K Pump Current conductance
    gkrbar=cell.gkrbar;//Rapidly Activating K Current conductance
    gksbar=cell.gksbar;//Slowly Activating K Current conductance
    gkpbar=cell.gkpbar;//Plateau K Current conductance
    gkachbar=cell.gkachbar;//Muscarinic-Receptor-Activated K Current conductance
    gtobar=cell.gtobar;//Transient Outward K Current conductance
    gkurbar=cell.gkurbar;//Ultra rapid delayed rectifier Outward K Current conductance
    gkibar=cell.gkibar;//Time-Independent K Current conductance
    vupbar=cell.vupbar;//SERCA strength
    gcabbar=cell.gcabbar;//Ca Background conductance
    gpcabar=cell.gpcabar;//Sarcolemmal Ca Pump Current conductance
    gncxbar=cell.gncxbar;//Na/Ca Exchanger flux conductance
    gcabar=cell.gcabar;//L-type Calcium Current conductance
  }
  return(*this);
}
void CCell::pace(double st) {
// -------------time step adjustment ------------------------
  for (int loop=0;loop<20;loop++) {
    double v=y[39];
    double dv=(vold-v)/dt;
    vold=v;
    if(fabs(dv)>25.0)// then finer time step when dv/dt large
    {
      ddt=dt/(10*20);
      for (int iii=0;iii<10;iii++) {
          pacex(st);
      }
    } else
    {
      ddt=dt/20;
      pacex(st);
    }
  }
  return;
}

void CCell::pacex(double st) {
  // Nernst Potentials
  ena_junc = (1/FoRT)*log(Nao/y[32]);     // [mV]
  ena_sl = (1/FoRT)*log(Nao/y[33]);       // [mV]
  ek = (1/FoRT)*log(Ko/y[35]);         // [mV]
  eca_junc = (1/FoRT/2)*log(Cao/y[36]);   // [mV]
  eca_sl = (1/FoRT/2)*log(Cao/y[37]);     // [mV]
  ecl = (1/FoRT)*log(Cli/Clo);            // [mV]



  comp_na_buffer();
  double J_CaB_cytosol=comp_ca_buffer_J_CaB_cytosol();
  double J_CaB_junction=comp_ca_buffer_J_CaB_junction();
  double J_CaB_sl=comp_ca_buffer_J_CaB_sl();


  //// Ion concentrations
  // SR Ca Concentrations
  double J_SRCarel=comp_sr_release();
  double J_serca=comp_serca();
  double J_SRleak=comp_sr_leak();
  ydot[30] = kon_csqn*y[31]*(Bmax_Csqn-y[30])-koff_csqn*y[30];       // Csqn      [mM/ms]
  ydot[31] = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot[30];         // Ca_sr     [mM/ms] %Ratio 3 leak current
  // ydot(31) = 0;

  // Na Concentrations
  double I_Na_junc, I_Na_sl;
  double I_Na=comp_ina(&I_Na_junc, &I_Na_sl);

  double I_nabk_junc=comp_I_nabk_junc();
  double I_nabk_sl=comp_I_nabk_sl();
  double I_nabk = I_nabk_junc+I_nabk_sl;

  double I_ncx_junc=comp_incx_junc();
  double I_ncx_sl=comp_incx_sl();
  double I_ncx = I_ncx_junc+I_ncx_sl;

  double I_nak_junc=comp_inak_junc();
  double I_nak_sl=comp_inak_sl();
  double I_nak = I_nak_junc+I_nak_sl;

  double I_CaNa_junc=comp_I_CaNa_junc();
  double I_CaNa_sl=comp_I_CaNa_sl();
  
  double I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   // [uA/uF]
  double I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   // [uA/uF]

  ydot[32] = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y[33]-y[32])-ydot[17];
  ydot[33] = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y[32]-y[33])+J_na_slmyo/Vsl*(y[34]-y[33])-ydot[18];
  // ydot[32] = 0;
  // ydot[33] = 0;
  ydot[34] = J_na_slmyo/Vmyo*(y[33]-y[34]);             // [mM/msec]
  // ydot[34] = 0;

  // K Concentration
  double I_to=comp_ito();
  double I_kr=comp_ikr();
  double I_ks=comp_iks();
  double I_ki=comp_ik1();
  double I_kp=comp_ikp();
  double I_kur=comp_ikur();
  double I_KAch=comp_ikach();
  double I_CaK=comp_I_CaK();

  
  double I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur+1*I_KAch;     // [uA/uF] %SVP: added IKur
  ydot[35] = 0; // -I_K_tot*Cmem/(Vmyo*Frdy);            // [mM/msec]

  // Ca Concentrations
  double I_cabk_junc=comp_I_cabk_junc();
  double I_cabk_sl=comp_I_cabk_sl();
  double I_cabk = I_cabk_junc+I_cabk_sl;

  double I_pca_junc=comp_I_pca_junc();
  double I_pca_sl=comp_I_pca_sl();
  double I_pca = I_pca_junc+I_pca_sl;

  double I_Ca_junc=comp_I_Ca_junc();
  double I_Ca_sl=comp_I_Ca_sl();
  
  double I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   // [uA/uF]
  double I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            // [uA/uF]
  ydot[36] = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y[37]-y[36])-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  // Ca_j
  ydot[37] = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y[36]-y[37]) + J_ca_slmyo/Vsl* (y[38]-y[37])-J_CaB_sl;   // Ca_sl
  // ydot[36]=0;
  // ydot[37]=0;
  ydot[38] = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y[37]-y[38]);  // [mM/msec]
  // ydot[38]=0;

  //Cl currents
  double I_ClCa_junc=comp_I_ClCa_junc();
  double I_ClCa_sl=comp_I_ClCa_sl();
  double I_ClCa = I_ClCa_junc+I_ClCa_sl;
  double I_Clbk=comp_I_Clbk();
  double I_ClCFTR=comp_I_ClCFTR();

  //// Membrane Potential
  double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
  double I_Cl_tot = I_ClCa+I_Clbk+I_ClCFTR;
  double I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
  double I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot; // [pA/pF]
  ydot[39] = -(I_tot-st);
  double vmax = ydot[39];

  //update variables
  for (int i=0;i<N;i++) {
    y[i]+=ydot[i]*ddt;
    //cout<<i<<"\t"<<y[i]<<"\t"<<ydot[i]<<endl;
  }

  _I_to=I_to;
  _I_kr=I_kr;
  _I_ks=I_ks;
  _I_ki=I_ki;
  _I_kp=I_kp;
  _I_kur=I_kur;
  _I_KAch=I_KAch;
  _I_CaK=I_CaK;
  _I_ncx=I_ncx;
  double I_Ca = I_Ca_junc+I_Ca_sl;
  double I_CaNa = I_CaNa_junc+I_CaNa_sl;
  double I_Catot = I_Ca+I_CaK+I_CaNa;
  _I_Catot=I_Catot;

}
double CCell::comp_ina(double *I_Na_junc, double *I_Na_sl) {
  double GNa_hh = 23*(1-0.1*AF);//  % [mS/uF]
  double GNa = 10*(1-0.1*AF);//  % [mS/uF]
  //// Membrane Currents
  // Fast I_Na
  double mss = 1 / (1 + pow(exp( -(56.86 + y[39]) / 9.03), 2));
  double taum = 0.1292 * exp(-(pow((y[39]+45.79) / 15.54, 2))) + 0.06487 * exp(-(pow((y[39]-4.823)/51.12, 2)));
  double ah = (y[39] >= -40) * (0) + (y[39] < -40) * (0.057 * exp( -(y[39] + 80) / 6.8 ));
  double bh = (y[39] >= -40) * (0.77 / (0.13*(1 + exp( -(y[39] + 10.66) / 11.1 )))) + (y[39] < -40) * ((2.7 * exp (0.079*y[39])) + 3.1*pow(10,5) * exp(0.3485 * y[39]));
  double tauh = 1 / (ah + bh);
  double hss= 1 / (pow(1 + exp((y[39] + 71.55)/7.43), 2));
  double aj = (y[39] >= -40) * (0)+(y[39] < -40) * (((-2.5428 * pow(10, 4) *exp(0.2444*y[39]) - 6.948* pow(10,-6)  * exp(-0.04391*y[39])) * (y[39] + 37.78)) /(1 + exp( 0.311 * (y[39] + 79.23) )));
  double bj = (y[39] >= -40) * ((0.6 * exp( 0.057 * y[39])) / (1 + exp( -0.1 * (y[39] + 32) ))) + (y[39]<-40) * ((0.02424 * exp( -0.01052 * y[39] )) / (1 + exp( -0.1378 * (y[39] + 40.14) )));
  double tauj = 1 / (aj + bj);
  tauj*=taujj;
  double jss = 1 /(pow (1 + exp(y[39] +71.55)/7.43,2));
  ydot[1] = (mss - y[1]) / taum;
  ydot[2] = (hss - y[2]) / tauh;
  ydot[3] = (jss - y[3]) / tauj;

  double I_Na_junc1 = Fjunc*GNa_hh*pow((y[1]),3)*y[2]*y[3]*(y[39]-ena_junc);
  double I_Na_sl1 = Fsl*GNa_hh*pow((y[1]),3)*y[2]*y[3]*(y[39]-ena_sl);
  //I_Na1 = I_Na_junc1+I_Na_sl1;
    
  //Late I_Na
  double GNaL = 0.0025*AF;
  double aml = 0.32*(y[39]+47.13)/(1-exp(-0.1*(y[39]+47.13)));
  double bml = 0.08*exp(-y[39]/11);
  double hlinf = 1/(1+exp((y[39]+91)/6.1));
  double tauhl = 600;
  ydot[60] = aml*(1-y[60])-bml*y[60];
  ydot[61] = (hlinf-y[61])/tauhl;
  double I_NaL_junc = Fjunc*GNaL*pow((y[60]),3)*y[61]*(y[39]-ena_junc);
  double I_NaL_sl = Fsl*GNaL*pow((y[60]),3)*y[61]*(y[39]-ena_sl);
  double I_NaL = I_NaL_junc + I_NaL_sl;

//  if (t<9050)
//      ydot[62] = 0;
//  else
      ydot[62] = I_NaL;


  // drug
  double portion = 1/(1+pow(10,(pH-pKa)));
  double drug_charged = drug * portion;
  double drug_neutral = drug * (1-portion);
  double kd_open = kd0 * exp(dd*y[39]*Frdy/(R*Temp));
  double kd_open_b = kd0_b * exp(dd*y[39]*Frdy/(R*Temp));

  // charged drug
  double kon = drug_charged * diffusion;
  double koff = kd_open * diffusion;
  double kcon = kon;
  double kcoff = koff;
  double kbon = kon; // bursting
  double kboff = kd_open_b * diffusion;
  double kcbon = kbon;
  double kcboff = kboff;

  // neutral drug
  double k_on = drug_neutral * diffusion;
  double k_off = k_off_0 * diffusion;
  double ki_on = k_on/1;
  double ki_off = ki_off_0 * diffusion;
  double kc_on = k_on/1;
  double kc_off = kc_off_0 * diffusion;
  // kb_on = k_on; // bursting
  // kb_off = k_off;
  // kbc_on = kc_on;
  // kbc_off = kc_off;

  // Drug Free
  double P1a1=3.802;
  double alphaNa1 = Tfactor_INa * P1a1/(P2a1*exp(-(y[39]+P3a1)/P4a1)+P5a1*exp(-(y[39]+P3a1)/P6a1));
  double alphaNa2 = Tfactor_INa * P1a1/(P2a1*exp(-(y[39]+P3a1)/P4a2)+P5a2*exp(-(y[39]+P3a1)/P6a1));
  double alphaNa3 = Tfactor_INa * P1a1/(P2a1*exp(-(y[39]+P3a1)/P4a3)+P5a3*exp(-(y[39]+P3a1)/P6a1));
  double betaNa1 = Tfactor_INa * P1b1*exp(-(y[39]+P3a1)/P2b1); // shift
  double betaNa2 = Tfactor_INa * P1b2*exp(-(y[39]-P2b2)/P2b1);
  double betaNa3 = Tfactor_INa * P1b3*exp(-(y[39]-P2b3)/P2b1);
  double alphaNa4 = Tfactor_INa * 1/(P1a4*exp(-(y[39]+P4a4)/P2a4)+P3a4);
  double alphaNa5 = Tfactor_INa * P1a5*exp(-(y[39]+P4a4)/P2a5);
  double betaNa5 = Tfactor_INa * (P1b5+P2b5*(y[39]+P4a4));
  double betaNa6 = Tfactor_INa * P1b6*exp(-(y[39])/P2b6);
  double alphaNa7 = Tfactor_INa * P1a7*exp((y[39])/P2a7);
  double betaNa7 = Tfactor_INa * P1b7*exp(-(y[39])/P2b7);
  double alphaNa8 = Tfactor_INa * P1a8;
  double betaNa8 = Tfactor_INa * P1b8;
  double alphaNa6 = alphaNa4/P1a6;
  double betaNa4 = (alphaNa3*alphaNa4*alphaNa5)/(betaNa3*betaNa5); // REV

  // Charged Drug
  double alphaNa1_c = alphaNa1; // constrained
  double alphaNa2_c = alphaNa2; // constrained
  double alphaNa3_c = Pa3_c * alphaNa3;                      // can be changed
  double betaNa1_c = betaNa1; // constrained
  double betaNa2_c = betaNa2; // constrained
  //betaNa3_c = betaNa3; // constrained (REV)
  double alphaNa4_c = Pa4_c * alphaNa4;                      // can be changed
  double alphaNa5_c = Pa5_c * alphaNa5;                      // can be changed
  double betaNa5_c = Pb5_c * betaNa5;                        // can be changed
  double alphaNa6_c = Pa6_c * alphaNa6;                      // can be changed
  double betaNa6_c = Pb6_c * betaNa6;                        // can be changed
  double alphaNa7_c = Pa7_c * alphaNa7;                      // can be changed
  double betaNa7_c = Pb7_c * betaNa7;                        // can be changed
  double alphaNa8_c = alphaNa8; // constrained
  double betaNa8_c = betaNa8; // constrained
      //betaNa4_c = betaNa4; // constrained (REV)

  // Neutral Drug
  double alphaNa1_n = alphaNa1; // constrained
  double alphaNa2_n = alphaNa2; // constrained
  double alphaNa3_n = Pa3_n * alphaNa3;                      // can be changed
  double betaNa1_n = betaNa1; // constrained
  double betaNa2_n = betaNa2; // constrained
  //betaNa3_n = betaNa3; // constrained (REV)
  double alphaNa4_n = Pa4_n * alphaNa4;                      // can be changed
  //alphaNa5_n = alphaNa5; // constrained (REV)
  double betaNa5_n = Pb5_n * betaNa5;                        // can be changed
  double alphaNa6_n = Pa6_n * alphaNa6;                      // can be changed
  double alphaNa7_n = Pa7_n * alphaNa7;                      // can be changed
  double betaNa6_n = alphaNa6_n*betaNa6/alphaNa6; // constrained (REV)
  double betaNa7_n = alphaNa7_n*betaNa7/alphaNa7; // constrained (REV)
  double alphaNa8_n = alphaNa8; // constrained
  double betaNa8_n = betaNa8; // constrained
  //betaNa4_n = betaNa4; // constrained (REV)

 // Microscopic reversibility (REV)
 double betaNa3_c;
 if (drug == 0 || drug_charged == 0 ) {
      betaNa3_c = 0;
  }
  else {
    betaNa3_c = ( betaNa3 * kcon * koff * alphaNa3_c ) / ( kon * kcoff * alphaNa3);
  }
  double betaNa4_c;
  if ( betaNa3_c == 0 ) {
    betaNa4_c = 0;
  }
  else {
    betaNa4_c = ( alphaNa3_c * alphaNa4_c * alphaNa5_c ) / ( betaNa3_c * betaNa5_c );
  }
  
  double alphaNa5_n;
  if ( drug == 0 || drug_neutral == 0 )
    alphaNa5_n = 0;
  else
    alphaNa5_n = ( ki_off * alphaNa5 * kc_on * betaNa5_n ) / ( ki_on * kc_off * betaNa5 );
  double betaNa3_n;
  if ( drug == 0 || drug_neutral == 0 )
    betaNa3_n = 0;
  else
      betaNa3_n = ( betaNa3 * kc_on * alphaNa3_n * k_off ) / ( kc_off * alphaNa3 * k_on );
 double betaNa4_n;
  if ( betaNa3_n == 0 )
    betaNa4_n = 0;
  else
    betaNa4_n = ( alphaNa5_n * alphaNa3_n * alphaNa4_n ) / ( betaNa5_n * betaNa3_n );

  // INa Markov Model State Variables /////////////////////////////////
  double CNa3 = y[63];
  double CNa2 = y[64];
  double CNa1 = y[65];
  double ONa = y[66];
  double LCNa3 = y[67];
  double LCNa2 = y[68];
  double LCNa1 = y[69];
  double LONa = y[70];
  double ICNa3 = y[71];
  double ICNa2 = y[72];
  double IFNa = y[73];
  double I1Na = y[74];
  //I2Na = (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3));

  double CNa3_c = y[75];
  double CNa2_c = y[76];
  double CNa1_c = y[77];
  double ONa_c = y[78];
  double LCNa3_c = y[79];
  double LCNa2_c = y[80];
  double LCNa1_c = y[81];
  double LONa_c = y[82];
  double ICNa3_c = y[83];
  double ICNa2_c = y[84];
  double IFNa_c = y[85];
  double I1Na_c = y[86];
  double I2Na_c = y[87];

  double CNa3_n = y[88];
  double CNa2_n = y[89];
  double CNa1_n = y[90];
  double ONa_n = y[91];
  double LCNa3_n = y[92];
  double LCNa2_n = y[93];
  double LCNa1_n = y[94];
  double LONa_n = y[95];
  double ICNa3_n = y[96];
  double ICNa2_n = y[97];
  double IFNa_n = y[98];
  double I1Na_n = y[99];
  double I2Na_n = y[100];

  double I2Na = (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3+ONa_c+CNa1_c+CNa2_c+CNa3_c+IFNa_c+I1Na_c+I2Na_c+ICNa2_c+ICNa3_c+LONa_c+LCNa1_c+LCNa2_c+LCNa3_c+ONa_n+CNa1_n+CNa2_n+CNa3_n+IFNa_n+I1Na_n+I2Na_n+ICNa2_n+ICNa3_n+LONa_n+LCNa1_n+LCNa2_n+LCNa3_n));

  // INa Markov Model ODEs //////////////////////////////
  // Drug Free
  double coeff_CNa2  = (betaNa1+alphaNa2+betaNa5+alphaNa8 +kcon+kc_on);
  double coeff_CNa1  = (betaNa2+alphaNa3+betaNa5+alphaNa8 +kcon+kc_on);
  double coeff_ONa   = (betaNa3+alphaNa4+alphaNa8 +kon+k_on);
  double coeff_IFNa  = (betaNa4+alphaNa5+alphaNa6+betaNa2 +ki_on);
  double coeff_I1Na  = (betaNa6+alphaNa7 +ki_on);
  double coeff_CNa3  = (alphaNa1+betaNa5+alphaNa8 +kcon+kc_on);
  double coeff_ICNa2 = (betaNa1+alphaNa2+alphaNa5 +ki_on);
  double coeff_ICNa3 = (alphaNa1+alphaNa5 +ki_on);
  double coeff_LONa  = (betaNa8+betaNa3 +kbon+k_on);
  double coeff_LCNa1 = (betaNa8+betaNa2+alphaNa3 +kcbon+kc_on);
  double coeff_LCNa2 = (betaNa8+betaNa1+alphaNa2 +kcbon+kc_on);
  double coeff_LCNa3 = (betaNa8+alphaNa1 +kcbon+kc_on);
  //coeff_I2Na  = (betaNa7 +ki_on);

  double dCNa2  = kcoff*CNa2_c+kc_off*CNa2_n+ betaNa8*LCNa2+alphaNa1*CNa3+betaNa2*CNa1+alphaNa5*ICNa2-(coeff_CNa2)*CNa2;
  double dCNa1  = kcoff*CNa1_c+kc_off*CNa1_n+ betaNa8*LCNa1+alphaNa2*CNa2+betaNa3*ONa+alphaNa5*IFNa-(coeff_CNa1)*CNa1;
  double dONa   = koff*ONa_c+k_off*ONa_n+ betaNa8*LONa+alphaNa3*CNa1+betaNa4*IFNa-(coeff_ONa)*ONa;
  double dIFNa  = ki_off*IFNa_n+ alphaNa4*ONa+betaNa5*CNa1+betaNa6*I1Na+alphaNa2*ICNa2-(coeff_IFNa)*IFNa;
  double dI1Na  = ki_off*I1Na_n+ alphaNa6*IFNa+betaNa7*I2Na-(coeff_I1Na)*I1Na;
  double dCNa3  = kcoff*CNa3_c+kc_off*CNa3_n+ betaNa8*LCNa3+betaNa1*CNa2+alphaNa5*ICNa3-(coeff_CNa3)*CNa3;
  double dICNa2 = ki_off*ICNa2_n+ alphaNa1*ICNa3+betaNa2*IFNa+betaNa5*CNa2-(coeff_ICNa2)*ICNa2;
  double dICNa3 = ki_off*ICNa3_n+ betaNa1*ICNa2+betaNa5*CNa3-(coeff_ICNa3)*ICNa3;
  double dLONa  = kboff*LONa_c+k_off*LONa_n+ alphaNa3*LCNa1+alphaNa8*ONa-(coeff_LONa)*LONa;
  double dLCNa1 = kcboff*LCNa1_c+kc_off*LCNa1_n+ alphaNa8*CNa1+alphaNa2*LCNa2+betaNa3*LONa-(coeff_LCNa1)*LCNa1;
  double dLCNa2 = kcboff*LCNa2_c+kc_off*LCNa2_n+ betaNa2*LCNa1+alphaNa8*CNa2+alphaNa1*LCNa3-(coeff_LCNa2)*LCNa2;
  double dLCNa3 = kcboff*LCNa3_c+kc_off*LCNa3_n+ alphaNa8*CNa3+betaNa1*LCNa2-(coeff_LCNa3)*LCNa3;
  //dI2Na  = ki_off*I2Na_n+ alphaNa7*I1Na-(coeff_I2Na)*I2Na;

  // Charged Drug
  double coeff_CNa2_c  = (betaNa1_c+alphaNa2_c+betaNa5_c+alphaNa8_c +kcoff);
  double coeff_CNa1_c  = (betaNa2_c+alphaNa3_c+betaNa5_c+alphaNa8_c +kcoff);
  double coeff_ONa_c   = (betaNa3_c+alphaNa4_c+alphaNa8_c +koff);
  double coeff_IFNa_c  = (betaNa4_c+alphaNa5_c+alphaNa6_c+betaNa2_c);
  double coeff_I1Na_c  = (betaNa6_c+alphaNa7_c);
  double coeff_CNa3_c  = (alphaNa1_c+betaNa5_c+alphaNa8_c +kcoff);
  double coeff_ICNa2_c = (betaNa1_c+alphaNa2_c+alphaNa5_c);
  double coeff_ICNa3_c = (alphaNa1_c+alphaNa5_c);
  double coeff_LONa_c  = (betaNa8_c+betaNa3_c +kboff);
  double coeff_LCNa1_c = (betaNa8_c+betaNa2_c+alphaNa3_c +kcboff);
  double coeff_LCNa2_c = (betaNa8_c+betaNa1_c+alphaNa2_c +kcboff);
  double coeff_LCNa3_c = (betaNa8_c+alphaNa1_c +kcboff);
  double coeff_I2Na_c  = (betaNa7_c);

  double dCNa2_c  = kcon*CNa2+ betaNa8_c*LCNa2_c+alphaNa1_c*CNa3_c+betaNa2_c*CNa1_c+alphaNa5_c*ICNa2_c-(coeff_CNa2_c)*CNa2_c;
  double dCNa1_c  = kcon*CNa1+ betaNa8_c*LCNa1_c+alphaNa2_c*CNa2_c+betaNa3_c*ONa_c+alphaNa5_c*IFNa_c-(coeff_CNa1_c)*CNa1_c;
  double dONa_c   = kon*ONa+ betaNa8_c*LONa_c+alphaNa3_c*CNa1_c+betaNa4_c*IFNa_c-(coeff_ONa_c)*ONa_c;
  double dIFNa_c  = alphaNa4_c*ONa_c+betaNa5_c*CNa1_c+betaNa6_c*I1Na_c+alphaNa2_c*ICNa2_c-(coeff_IFNa_c)*IFNa_c;
  double dI1Na_c  = alphaNa6_c*IFNa_c+betaNa7_c*I2Na_c-(coeff_I1Na_c)*I1Na_c;
  double dCNa3_c  = kcon*CNa3+ betaNa8_c*LCNa3_c+betaNa1_c*CNa2_c+alphaNa5_c*ICNa3_c-(coeff_CNa3_c)*CNa3_c;
  double dICNa2_c = alphaNa1_c*ICNa3_c+betaNa2_c*IFNa_c+betaNa5_c*CNa2_c-(coeff_ICNa2_c)*ICNa2_c;
  double dICNa3_c = betaNa1_c*ICNa2_c+betaNa5_c*CNa3_c-(coeff_ICNa3_c)*ICNa3_c;
  double dLONa_c  = kbon*LONa+ alphaNa3_c*LCNa1_c+alphaNa8_c*ONa_c-(coeff_LONa_c)*LONa_c;
  double dLCNa1_c = kcbon*LCNa1+ alphaNa8_c*CNa1_c+alphaNa2_c*LCNa2_c+betaNa3_c*LONa_c-(coeff_LCNa1_c)*LCNa1_c;
  double dLCNa2_c = kcbon*LCNa2+ betaNa2_c*LCNa1_c+alphaNa8_c*CNa2_c+alphaNa1_c*LCNa3_c-(coeff_LCNa2_c)*LCNa2_c;
  double dLCNa3_c = kcbon*LCNa3+ alphaNa8_c*CNa3_c+betaNa1_c*LCNa2_c-(coeff_LCNa3_c)*LCNa3_c;
  double dI2Na_c  = alphaNa7_c*I1Na_c-(coeff_I2Na_c)*I2Na_c;

  // Neutral Drug
  double coeff_CNa2_n  = (betaNa1_n+alphaNa2_n+betaNa5_n+alphaNa8_n +kc_off);
  double coeff_CNa1_n  = (betaNa2_n+alphaNa3_n+betaNa5_n+alphaNa8_n +kc_off);
  double coeff_ONa_n   = (betaNa3_n+alphaNa4_n+alphaNa8_n +k_off);
  double coeff_IFNa_n  = (betaNa4_n+alphaNa5_n+alphaNa6_n+betaNa2_n +ki_off);
  double coeff_I1Na_n  = (betaNa6_n+alphaNa7_n +ki_off);
  double coeff_CNa3_n  = (alphaNa1_n+betaNa5_n+alphaNa8_n +kc_off);
  double coeff_ICNa2_n = (betaNa1_n+alphaNa2_n+alphaNa5_n +ki_off);
  double coeff_ICNa3_n = (alphaNa1_n+alphaNa5_n +ki_off);
  double coeff_LONa_n  = (betaNa8_n+betaNa3_n +k_off);
  double coeff_LCNa1_n = (betaNa8_n+betaNa2_n+alphaNa3_n +kc_off);
  double coeff_LCNa2_n = (betaNa8_n+betaNa1_n+alphaNa2_n +kc_off);
  double coeff_LCNa3_n = (betaNa8_n+alphaNa1_n +kc_off);
  double coeff_I2Na_n  = (betaNa7_n +ki_off);

  double dCNa2_n  = kc_on*CNa2+ betaNa8_n*LCNa2_n+alphaNa1_n*CNa3_n+betaNa2_n*CNa1_n+alphaNa5_n*ICNa2_n-(coeff_CNa2_n)*CNa2_n;
  double dCNa1_n  = kc_on*CNa1+ betaNa8_n*LCNa1_n+alphaNa2_n*CNa2_n+betaNa3_n*ONa_n+alphaNa5_n*IFNa_n-(coeff_CNa1_n)*CNa1_n;
  double dONa_n   = k_on*ONa+ betaNa8_n*LONa_n+alphaNa3_n*CNa1_n+betaNa4_n*IFNa_n-(coeff_ONa_n)*ONa_n;
  double dIFNa_n  = ki_on*IFNa+ alphaNa4_n*ONa_n+betaNa5_n*CNa1_n+betaNa6_n*I1Na_n+alphaNa2_n*ICNa2_n-(coeff_IFNa_n)*IFNa_n;
  double dI1Na_n  = ki_on*I1Na+ alphaNa6_n*IFNa_n+betaNa7_n*I2Na_n-(coeff_I1Na_n)*I1Na_n;
  double dCNa3_n  = kc_on*CNa3+ betaNa8_n*LCNa3_n+betaNa1_n*CNa2_n+alphaNa5_n*ICNa3_n-(coeff_CNa3_n)*CNa3_n;
  double dICNa2_n = ki_on*ICNa2+ alphaNa1_n*ICNa3_n+betaNa2_n*IFNa_n+betaNa5_n*CNa2_n-(coeff_ICNa2_n)*ICNa2_n;
  double dICNa3_n = ki_on*ICNa3+ betaNa1_n*ICNa2_n+betaNa5_n*CNa3_n-(coeff_ICNa3_n)*ICNa3_n;
  double dLONa_n  = k_on*LONa+ alphaNa3_n*LCNa1_n+alphaNa8_n*ONa_n-(coeff_LONa_n)*LONa_n;
  double dLCNa1_n = kc_on*LCNa1+ alphaNa8_n*CNa1_n+alphaNa2_n*LCNa2_n+betaNa3_n*LONa_n-(coeff_LCNa1_n)*LCNa1_n;
  double dLCNa2_n = kc_on*LCNa2+ betaNa2_n*LCNa1_n+alphaNa8_n*CNa2_n+alphaNa1_n*LCNa3_n-(coeff_LCNa2_n)*LCNa2_n;
  double dLCNa3_n = kc_on*LCNa3+ alphaNa8_n*CNa3_n+betaNa1_n*LCNa2_n-(coeff_LCNa3_n)*LCNa3_n;
  double dI2Na_n  = ki_on*I2Na+ alphaNa7_n*I1Na_n-(coeff_I2Na_n)*I2Na_n;

  // INa Markov Model dydt /////////////////////////////
  ydot[63]=dCNa3;
  ydot[64]=dCNa2;
  ydot[65]=dCNa1;
  ydot[66]=dONa;
  ydot[67]=dLCNa3;
  ydot[68]=dLCNa2;
  ydot[69]=dLCNa1;
  ydot[70]=dLONa;
  ydot[71]=dICNa3;
  ydot[72]=dICNa2;
  ydot[73]=dIFNa;
  ydot[74]=dI1Na;
  ydot[75]=dCNa3_c;
  ydot[76]=dCNa2_c;
  ydot[77]=dCNa1_c;
  ydot[78]=dONa_c;
  ydot[79]=dLCNa3_c;
  ydot[80]=dLCNa2_c;
  ydot[81]=dLCNa1_c;
  ydot[82]=dLONa_c;
  ydot[83]=dICNa3_c;
  ydot[84]=dICNa2_c;
  ydot[85]=dIFNa_c;
  ydot[86]=dI1Na_c;
  ydot[87]=dI2Na_c;
  ydot[88]=dCNa3_n;
  ydot[89]=dCNa2_n;
  ydot[90]=dCNa1_n;
  ydot[91]=dONa_n;
  ydot[92]=dLCNa3_n;
  ydot[93]=dLCNa2_n;
  ydot[94]=dLCNa1_n;
  ydot[95]=dLONa_n;
  ydot[96]=dICNa3_n;
  ydot[97]=dICNa2_n;
  ydot[98]=dIFNa_n;
  ydot[99]=dI1Na_n;
  ydot[100]=dI2Na_n;
  

  // INa Markov Model Output Current ///////////////////////
  double I_Na_junc2 = Fjunc*GNa*(ONa+LONa)*(y[39]-ena_junc);
  double I_Na_sl2 = Fsl*GNa*(ONa+LONa)*(y[39]-ena_sl);
  // I_Na2 = I_Na_junc2+I_Na_sl2;

  // Compute Total INa (HH or Markov) //////////////////////
  *I_Na_junc = gnabar*( (I_Na_junc1+I_NaL_junc)*(1-flagMina)+I_Na_junc2*flagMina);
  *I_Na_sl = gnabar*((I_Na_sl1+I_NaL_sl)*(1-flagMina)+I_Na_sl2*flagMina);
  double I_Na = *I_Na_junc+*I_Na_sl;
  return I_Na;
}

  //// I_nabk: Na Background Current
double CCell::comp_I_nabk_junc(void) {
  double I_nabk_junc = gnabbar*(Fjunc*GNaB*(y[39]-ena_junc));
  return I_nabk_junc;
}
double CCell::comp_I_nabk_sl(void) {
  double I_nabk_sl = gnabbar*(Fsl*GNaB*(y[39]-ena_sl));
  return I_nabk_sl;
}

  //// I_nak: Na/K Pump Current
double CCell::comp_inak_junc(void) {
  double KmNaip = 11*(1-0.25*ISO);         // [mM]11
  double sigma = (exp(Nao/67.3)-1)/7;
  double fnak = 1/(1+0.1245*exp(-0.1*y[39]*FoRT)+0.0365*sigma*exp(-y[39]*FoRT));
  double I_nak_junc = gnakbar*(1*Fjunc*IbarNaK*fnak*Ko /(1+pow(KmNaip/(y[32]),4)) /(Ko+KmKo));
  return I_nak_junc;
}
double CCell::comp_inak_sl(void) {
  double KmNaip = 11*(1-0.25*ISO);         // [mM]11
  double sigma = (exp(Nao/67.3)-1)/7;
  double fnak = 1/(1+0.1245*exp(-0.1*y[39]*FoRT)+0.0365*sigma*exp(-y[39]*FoRT));
  double I_nak_sl = gnakbar*(1*Fsl*IbarNaK*fnak*Ko /(1+pow((KmNaip/y[33]),4)) /(Ko+KmKo));
  return I_nak_sl;
}


  //// I_kr: Rapidly Activating K Current
double CCell::comp_ikr(void) {
  double factor_rano_kr;
  double IC50_kr;
  if ((flagMina == 1) && (drug_index == 1) ){// ranolazine effect
    IC50_kr = 35*(1e-6);
    factor_rano_kr = 1/(1+(drug/IC50_kr));
  }
  else{
    factor_rano_kr = 1;
  }
  double gkr = 0.035*sqrt(Ko/5.4)*factor_rano_kr;

  // gkr = 0.035*sqrt(Ko/5.4);
  double xrss = 1/(1+exp(-(y[39]+10)/5));
  double tauxr = 550/(1+exp((-22-y[39])/9))*6/(1+exp((y[39]-(-11))/9))+230/(1+exp((y[39]-(-40))/20));
  ydot[12] = (xrss-y[12])/tauxr;
  double rkr = 1/(1+exp((y[39]+74)/24));
  double I_kr = gkrbar*(gkr*y[12]*rkr*(y[39]-ek));

  return I_kr;
}

  //// I_ks: Slowly Activating K Current
double CCell::comp_iks(void) {
  double eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y[35]+pNaK*y[34]));
  double gks_junc = (1+1*AF+2*ISO)*0.0035*1;
  double gks_sl = (1+1*AF+2*ISO)*0.0035*1; //FRA
  double  xsss = 1 / (1+exp(-(y[39]+40*ISO + 3.8)/14.25)); // fitting Fra
  double tauxs = 990.1/(1+exp(-(y[39]+40*ISO+2.436)/14.12));
  ydot[13] = (xsss-y[13])/tauxs;
  double I_ks_junc = gksbar*(Fjunc*gks_junc*pow((y[13]),2)*(y[39]-eks));
  double I_ks_sl = gksbar*(Fsl*gks_sl*pow((y[13]),2)*(y[39]-eks));
  double I_ks = I_ks_junc+I_ks_sl;
  return I_ks;
}

  //// I_kp: Plateau K Current
double CCell::comp_ikp(void) {
  double kp_kp = 1/(1+exp(7.488-y[39]/5.98));
  double I_kp_junc = gkpbar*(Fjunc*gkp*kp_kp*(y[39]-ek));
  double I_kp_sl = gkpbar*(Fsl*gkp*kp_kp*(y[39]-ek));
  double I_kp = I_kp_junc+I_kp_sl;
  return I_kp;
}

  //// I_k,ach: Muscarinic-Receptor-Activated K Current
double CCell::comp_ikach(void) {
  double I_KAch = gkachbar*(1/(1+pow((0.03/Ach),2.1))*(0.08+0.4/(1+exp((y[39]+91)/12)))*(y[39]-ek));
  return I_KAch;
}

  //// I_to: Transient Outward K Current (slow and fast components)
double CCell::comp_ito(void) {
  // modified for human myocytes
  double GtoFast = (1.0-0.7*AF)*0.165*1.0; //nS/pF maleckar; //human atrium

  // 11/12/09; changed Itof to that from maleckar/giles/2009;
  // removed I_tos
  // atrium

  // equations for activation;
  double xtoss = ( (1)/ ( 1 + exp( -(y[39]+1.0)/11.0 ) ) );
  double tauxtof = 3.5*exp(-(pow((y[39]/30.0),2.0)))+1.5;

  // equations for inactivation;
  double ytoss = ( (1.0)/ ( 1 + exp( (y[39]+40.5)/11.5) ) ) ;
  double tauytof = 25.635*exp(-(pow(((y[39]+52.45)/15.8827),2.0)))+24.14; //14.14

  ydot[10] = (xtoss-y[10])/tauxtof;
  ydot[11] = (ytoss-y[11])/tauytof;
  double I_tof = gtobar*(1.0*GtoFast*y[10]*y[11]*(y[39]-ek));

  double I_to = I_tof;
  return I_to;
}
  //// I_kur: Ultra rapid delayed rectifier Outward K Current
  // Equation for IKur; from Maleckar et al. 2009 - EG
  // atrium
double CCell::comp_ikur(void) {
  // equations for activation;
  double Gkur = (1.0-0.5*AF)*(1+2*ISO)* 0.045*(1+0.2*RA); //nS/pF maleckar 0.045
  double xkurss = ( (1)/ ( 1 + exp( (y[39]+6)/-8.6 ) ) );
  double tauxkur = 9/(1+exp((y[39]+5)/12.0))+0.5;

  // equations for inactivation;
  double ykurss = ( (1)/ ( 1 + exp( (y[39]+7.5)/10 ) ) );
  double tauykur = 590/(1+exp((y[39]+60)/10.0))+3050;

  ydot[58] = (xkurss-y[58])/tauxkur;
  ydot[59] = (ykurss-y[59])/tauykur;
  double I_kur = gkurbar*(1*Gkur*y[58]*y[59]*(y[39]-ek));
  return I_kur;
}
  //// I_ki: Time-Independent K Current
double CCell::comp_ik1(void) {
  double aki = 1.02/(1+exp(0.2385*(y[39]-ek-59.215)));
  double bki =(0.49124*exp(0.08032*(y[39]+5.476-ek)) + exp(0.06175*(y[39]-ek-594.31))) /(1 + exp(-0.5143*(y[39]-ek+4.753)));
  double kiss = aki/(aki+bki);

  // I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(y[39]-ek);
  // SVP 11/11/09
  // multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
  // resting potential
  double I_ki = gkibar*((1+1*AF)*0.0525*sqrt(Ko/5.4)*kiss*(y[39]-ek));
  return I_ki;
}

  //// I_ClCa & I_Clbk: Ca-activated Cl Current and Background Cl Current
double CCell::comp_I_ClCa_junc(void) {
  double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y[36])*(y[39]-ecl);
  return I_ClCa_junc;
}
double CCell::comp_I_ClCa_sl(void) {
  double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y[37])*(y[39]-ecl);
  return I_ClCa_sl;
}
double CCell::comp_I_Clbk(void) {
  double I_Clbk = GClB*(y[39]-ecl);
  return I_Clbk;
}
double CCell::comp_I_ClCFTR(void) {
  double I_ClCFTR = GClCFTR*(y[39]-ecl);
  return I_ClCFTR;
}

  //// Na and Ca Buffering
void CCell::comp_na_buffer(void) {
  // Na Buffers
  ydot[17] = kon_na*y[32]*(Bmax_Naj-y[17])-koff_na*y[17];        // NaBj      [mM/ms]
  ydot[18] = kon_na*y[33]*(Bmax_Nasl-y[18])-koff_na*y[18];       // NaBsl     [mM/ms]
}
double CCell::comp_ca_buffer_J_CaB_cytosol(void) {
  double koff_tncl = (1+0.5*ISO)*19.6e-3; // [1/ms]
  // Cytosolic Ca Buffers
  ydot[19] = kon_tncl*y[38]*(Bmax_TnClow-y[19])-koff_tncl*y[19];            // TnCL      [mM/ms]
  ydot[20] = kon_tnchca*y[38]*(Bmax_TnChigh-y[20]-y[21])-koff_tnchca*y[20]; // TnCHc     [mM/ms]
  ydot[21] = kon_tnchmg*Mgi*(Bmax_TnChigh-y[20]-y[21])-koff_tnchmg*y[21];   // TnCHm     [mM/ms]
  ydot[22] = kon_cam*y[38]*(Bmax_CaM-y[22])-koff_cam*y[22];                 // CaM       [mM/ms]
  ydot[23] = kon_myoca*y[38]*(Bmax_myosin-y[23]-y[24])-koff_myoca*y[23];    // Myosin_ca [mM/ms]
  ydot[24] = kon_myomg*Mgi*(Bmax_myosin-y[23]-y[24])-koff_myomg*y[24];      // Myosin_mg [mM/ms]
  ydot[25] = kon_sr*y[38]*(Bmax_SR-y[25])-koff_sr*y[25];                    // SRB       [mM/ms]
  //J_CaB_cytosol = sum(ydot(19:25)); % wrong formulation
  double J_CaB_cytosol = ydot[19]+ydot[20]+ydot[22]+ydot[23]+ydot[25];
  return J_CaB_cytosol;
}

double CCell::comp_ca_buffer_J_CaB_junction(void) {
  // Junctional Ca Buffers
  ydot[26] = kon_sll*y[36]*(Bmax_SLlowj-y[26])-koff_sll*y[26];       // SLLj      [mM/ms]
  ydot[28] = kon_slh*y[36]*(Bmax_SLhighj-y[28])-koff_slh*y[28];      // SLHj      [mM/ms]
  double J_CaB_junction = ydot[26]+ydot[28];
  return J_CaB_junction;
}

double CCell::comp_ca_buffer_J_CaB_sl(void) {
  // SL Ca Buffers
  ydot[27] = kon_sll*y[37]*(Bmax_SLlowsl-y[27])-koff_sll*y[27];      // SLLsl     [mM/ms]
  ydot[29] = kon_slh*y[37]*(Bmax_SLhighsl-y[29])-koff_slh*y[29];     // SLHsl     [mM/ms]
  double J_CaB_sl = ydot[27]+ydot[29];
  return J_CaB_sl;
}

double CCell::comp_sr_release(void) {
  //// SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
  double koCa = 10+20*AF+10*ISO*(1-AF); // [mM^-2 1/ms]   %default 10   modified 20
  double MaxSR = 15;
  double MinSR = 1;
  double kCaSR = MaxSR - (MaxSR-MinSR)/(1+pow((ec50SR/y[31]),2.5));
  double koSRCa = (1)*koCa/kCaSR;//
  double kiSRCa = kiCa*kCaSR;
  double RI = 1-y[14]-y[15]-y[16];
  ydot[14] = (kim*RI-kiSRCa*y[36]*y[14])-(koSRCa*pow(y[36],2)*y[14]-kom*y[15]);   // R
  ydot[15] = (koSRCa*pow(y[36],2)*y[14]-kom*y[15])-(kiSRCa*y[36]*y[15]-kim*y[16]);// O
  ydot[16] = (kiSRCa*y[36]*y[15]-kim*y[16])-(kom*y[16]-koSRCa*pow(y[36],2)*RI);   // I
  double J_SRCarel = ks*y[15]*(y[31]-y[36]);          // [mM/ms]
  return J_SRCarel;
}
double CCell::comp_serca(void) {
  double Kmf = (2.5-1.5*ISO)*0.246e-3; // [mM] default 2.5-1.25*ISO
  double J_serca = vupbar*(pow(Q10SRCaP,Qpow)*Vmax_SRCaP* (pow((y[38]/Kmf),hillSRCaP)-pow((y[31]/Kmr),hillSRCaP))/(pow(1+(y[38]/Kmf),hillSRCaP)+pow((y[31]/Kmr),hillSRCaP)));
  return J_serca;
}

double CCell::comp_sr_leak(void) {
  double J_SRleak = (1.0+0.25*AF)*5.348e-6*(y[31]-y[36]);           //   [mM/ms]}
  return J_SRleak;
}

  //// I_cabk: Ca Background Current
double CCell::comp_I_cabk_junc(void) {
  double I_cabk_junc = gcabbar*(Fjunc*GCaB*(y[39]-eca_junc));
  return I_cabk_junc;
}
double CCell::comp_I_cabk_sl(void) {
  double I_cabk_sl = gcabbar*(Fsl*GCaB*(y[39]-eca_sl));
  return I_cabk_sl;
}


  //// I_pca: Sarcolemmal Ca Pump Current
double CCell::comp_I_pca_junc(void) {
  double I_pca_junc = gpcabar*(Fjunc*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(y[36],1.6)/(pow(KmPCa,1.6)+pow(y[36],1.6)));
  return I_pca_junc;
}
double CCell::comp_I_pca_sl(void) {
  double I_pca_sl = gpcabar*(Fsl*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(y[37],1.6)/(pow((KmPCa),1.6)+pow((y[37]),1.6)));
  return I_pca_sl;
}

  //// I_ncx: Na/Ca Exchanger flux
double CCell::comp_incx_junc(void) {
  double IbarNCX = (1+0.4*AF)*3.15;      // [uA/uF] 5.5 before - 9 in rabbit
  double Ka_junc = 1/(1+pow((Kdact/y[36]),2));
  double s1_junc = exp(nu*y[39]*FoRT)*pow((y[32]),3)*Cao;
  double s2_junc = exp((nu-1)*y[39]*FoRT)*pow((Nao),3)*y[36];
  double s3_junc = KmCai*pow(Nao,3)*(1+pow((y[32]/KmNai),3)) + pow(KmNao,3)*y[36]*(1+y[36]/KmCai)+KmCao*pow(y[32],3)+pow(y[32],3)*Cao+pow(Nao,3)*y[36];
  double I_ncx_junc = gncxbar*(Fjunc*IbarNCX*pow(Q10NCX,Qpow)*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y[39]*FoRT)));
  return I_ncx_junc;
}

double CCell::comp_incx_sl(void) {
  double IbarNCX = (1+0.4*AF)*3.15;      // [uA/uF] 5.5 before - 9 in rabbit
  double Ka_sl = 1/(1+pow((Kdact/y[37]),2));
  double s1_sl = exp(nu*y[39]*FoRT)*pow((y[33]),3)*Cao;
  double s2_sl = exp((nu-1)*y[39]*FoRT)*pow(Nao,3)*y[37];
  double s3_sl = KmCai*pow(Nao,3)*(1+pow((y[33]/KmNai),3)) + pow(KmNao,3)*y[37]*(1+y[37]/KmCai)+KmCao*pow(y[33],3)+pow(y[33],3)*Cao+pow(Nao,3)*y[37];
  double I_ncx_sl = gncxbar*(Fsl*IbarNCX*pow(Q10NCX,Qpow)*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y[39]*FoRT)));
  return I_ncx_sl;
}




double CCell::comp_I_Ca_junc(void) {
  //// I_Ca: L-type Calcium Current
  double dss = 1/(1+exp(-(y[39]+3*ISO+9)/6)); //in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
  double taud = 1*dss*(1-exp(-(y[39]+3*ISO+9)/6))/(0.035*(y[39]+3*ISO+9));
  double fss = 1/(1+exp((y[39]+3*ISO+30)/7))+0.2/(1+exp((50-y[39]-3*ISO)/20)); //in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
  double tauf = 1/(0.0197*exp(-pow((0.0337*(y[39]+3*ISO+25)),2 ))+0.02);
  tauf*=tauff;
  ydot[4] = (dss-y[4])/taud;
  ydot[5] = (fss-y[5])/tauf;
  ydot[6] = 1.7*y[36]*(1-y[6])-1*11.9e-3*y[6]; // fCa_junc   koff!!!!!!!!
  ydot[7] = 1.7*y[37]*(1-y[7])-1*11.9e-3*y[7]; // fCa_sl
  double fcaCaj= 0.1/(1+(0.01/y[36]));
  fcaCaj= 0;
  double pCa = (1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       // [cm/sec]
  double ibarca_j = pCa*4*(y[39]*Frdy*FoRT) * (0.341*y[36]*exp(2*y[39]*FoRT)-0.341*Cao) /(exp(2*y[39]*FoRT)-1);
  double I_Ca_junc = gcabar*((Fjunc_CaL*ibarca_j*y[4]*y[5]*((1-y[6])+fcaCaj)*pow(Q10CaL,Qpow))*0.45);
  return I_Ca_junc;
}
double CCell::comp_I_Ca_sl(void) {
  double fcaCaMSL= 0.1/(1+(0.01/y[37]));
  fcaCaMSL=0;
  double pCa = (1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       // [cm/sec]
  double ibarca_sl = pCa*4*(y[39]*Frdy*FoRT) * (0.341*y[37]*exp(2*y[39]*FoRT)-0.341*Cao) /(exp(2*y[39]*FoRT)-1);
  double I_Ca_sl = gcabar*((Fsl_CaL*ibarca_sl*y[4]*y[5]*((1-y[7])+fcaCaMSL)*pow(Q10CaL,Qpow))*0.45);
  return I_Ca_sl;
}
double CCell::comp_I_CaK(void) {
  double pK = (1+0.5*ISO)*(1-0.5*AF)*1.35e-7;        // [cm/sec]
  double fcaCaj= 0.1/(1+(0.01/y[36]));
  fcaCaj= 0;
  double fcaCaMSL= 0.1/(1+(0.01/y[37]));
  fcaCaMSL=0;
  double ibark = pK*(y[39]*Frdy*FoRT)*(0.75*y[35]*exp(y[39]*FoRT)-0.75*Ko) /(exp(y[39]*FoRT)-1);
  double I_CaK = gcabar*((ibark*y[4]*y[5]*(Fjunc_CaL*(fcaCaj+(1-y[6]))+Fsl_CaL*(fcaCaMSL+(1-y[7])))*pow(Q10CaL,Qpow))*0.45);
  return I_CaK;
}
double CCell::comp_I_CaNa_junc(void) {
  double pNa = (1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       // [cm/sec]
  double fcaCaj= 0.1/(1+(0.01/y[36]));
  fcaCaj= 0;
  double ibarna_j = pNa*(y[39]*Frdy*FoRT) *(0.75*y[32]*exp(y[39]*FoRT)-0.75*Nao)  /(exp(y[39]*FoRT)-1);
  double I_CaNa_junc = gcabar*((Fjunc_CaL*ibarna_j*y[4]*y[5]*((1-y[6])+fcaCaj)*pow(Q10CaL,Qpow))*0.45);
  return I_CaNa_junc;
}
double CCell::comp_I_CaNa_sl(void) {
  double pNa = (1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       // [cm/sec]
  double fcaCaMSL= 0.1/(1+(0.01/y[37]));
  fcaCaMSL=0;
  double ibarna_sl = pNa*(y[39]*Frdy*FoRT) *(0.75*y[33]*exp(y[39]*FoRT)-0.75*Nao)  /(exp(y[39]*FoRT)-1);
  double I_CaNa_sl = gcabar*((Fsl_CaL*ibarna_sl*y[4]*y[5]*((1-y[7])+fcaCaMSL)*pow(Q10CaL,Qpow))*0.45);
  return I_CaNa_sl;
}


