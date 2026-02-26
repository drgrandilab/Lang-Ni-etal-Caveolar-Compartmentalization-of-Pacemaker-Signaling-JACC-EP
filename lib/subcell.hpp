//Juan G Restrepo Model

// standard parameter is Ca leak param; if you don't want ___USE_ORG_PARAM
// no TC in Cs
// Cs has buffer if you don't want, use ___NO_CS_BUFFER
// fine mesh is optional// to use fmesh = 3 or 5 ; for 5 dt=0.005


//compile Options
//___USE_ORG_PARAM
//___NO_CS_BUFFER
#define ___NCX
#define __RYR_UNIFORM
// #define ___DEBUG
//___CPDIFF:use small dt; don't use with ___NO_CS_BUFFER

// *) xw: build a flag for removing the axial tubule
// *) xw: add the flag for indentify the location tubule
// #define __UNIFORM_TUBULE // *) close all the change of tubule distribution
#define ___CPDIFF
// #define __UNIFORM_TUBULE
#define ___UNIFORM

#ifndef ___SUBCELL_H
#define ___SUBCELL_H
//#endif
#define _USE_MATH_DEFINES

#include <iostream>
using namespace std;
#include <cmath>
#include <fstream>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include "ical13.hpp"
#include "icat.hpp"
#include "RyR.hpp"
#include "xor_rand.hpp"
#include <vector>
#include "inaca.hpp"
#include "LTCC_unitary.hpp"



class CSubcell {
public:

  //cell parameters
  int nxny;

  int finemesh3;
  int nnx;//the number of finemesh (x)
  int nny;//the number of finemesh (y)
  int nnz;//the number of finemesh (z)
  int nnxnny;

  double svr;// *) xw: set a rate of # CRU segments/ # all segments: surface to volume ratio

  double vi;
  double vs;
  static const  double vp_ave;


  double tausT;
  double tausL;
  double taumninv;
  double tauiT;//D/dx/dx*dt
  double tauiL;
  double taunsrT;
  double taunsrL;

  void computeIci(void);
  void computeIcnsr(void);
  double *Ici, *Icnsr;
#ifdef ___NO_CS_BUFFER
  double *csmn;
  void computecsmn(void);
#else
  void computeIcs(void);
  double *Ics;
  double *Idps;
#endif

  int bino(double num, double p, int ii);
  double calcvp(double mean, double std, double lim1, double lim2, int ii);

  double dt;
  double vup, kup, KNSR;
  double vnaca;
  double Jmax;
  double taup;
  double tausi;
  double tautr;
  double gca;//scale 1

  double cao;  // extracellular concentration of Na/Ca;
  double nao;  // extracellular concentration of Na/Ca;

  double gleak;
  double gcabk;
  double qslcap;
  double BCSQN;

  //  RyR parameters
  double Kc; // Dissociation constant of CSQN
  double nM; // Buffering capacity of CSQN monomers
  double nD; // Buffering capacity of CSQN dimers
  double hh; // Dimerization exponent
  double KK; // Dimerization constant
  double rhoinf; //Asymptotic ratio of dimers to monomers

  double Ku; // CSQN-unbound opening rate
  double Kb; // CSQN-bound opening rate
  double tauu; // CSQN unbinding timescale
  double taub; // CSQN binding timescale
  double tauc1;
  double tauc2;
  double BCSQN0; // CSQN concentration

  int seed;
  bool initialized;
  int bc;


  //#ifdef ___SIGMOID
  double Kcp;
  double pedk12;
  double pedk43;
  //#endif
  int NoCaL;//# of CaL channel per CaRU



  double MaxSR;
  double MinSR;
  double ec50SR;
  double hkosrca;







  void delarray(void);

#ifdef ___EGTA
  double BEGTA;
#endif

// public:
// *) xw: add the layer paramters
  int layer;//the numbers of layers with surface channels on each side
// *)xw: send the new RyR expression outside to the pace.cc
#ifndef __RYR_UNIFORM
  double *sigma;
#endif

  unsigned int *xsx, *xsy, *xsz, *xsw;
  int nx;//the number of CRU (x)
  int ny;//the number of CRU (y)
  int nz;//the number of CRU (z)
  int finemesh;// X times fine mesh (must be odd number , must be >3)
  double xi;

  int n;//the number of CRU x*y*z
  int nn;//the number of fine mesh

  static const double vjsr;
  double vnsr;

  static const double Cmem; //  *) xw: add the C_membrane

  std::vector<ical13> ical13_vec;
  std::vector<ical13> ical13_vec_sl;
  std::vector<ical13> ical13_vec_ci;

  int LTCC_alpha = 8; // num of LTCC in Cleft space; 
  int LTCC_gamma = 0; // num of LTCC in cytosolic space; 
  int LTCC_beta = 0;
  // int LTCC_alpha = 8;


  std::vector<xor_rand> random_num;

  std::vector<int>tubule_ID_vec;
  std::vector<RyR> RyR_vec;

  CSubcell(int sizex = 34, int sizey = 8, int sizez = 8, int fmesh = 1, double xii = 1);
  void init(double initci = 0.1, double initcj = 750,int LTCC_alpha_in = 8, int LTCC_gamma_in = 0, int rand_seed = 0); // sc.init(0.1, 1000)
  //  allocate array ci[nn], cs[nn], cp[nn], cjsr[n], cnsr[nn], cati[nn]
  //  no ___DETERMINISTIC ryr parameters ryr1[n], ryr2[n], ryr3[n], nryr[n]
  //  vp[n], Jmaxx[n], cscp1[nn], cscp2[nn], Ici[nn], Icnsr[nn]
  //  no ___NO_CS_BUFFER Ics[nn], Idps[nn], Ics[id]=cscp1[id]=0
  //  ___EGTA ifdef -> caEGTAi[nn], caEGTAs[nn]
  //  Itr[nn], crupos[n]
  //  random number xsx[n], xsy[n], xsz[n], xsw[n]
  //  initial ci=cs=initci, cnsr=initcj, Icnsr=0.
  //  Ici[id]= cscp2[id]= Itr[id]= 0
  virtual ~CSubcell();

  //  double *ci,*cs,*cp,*cjsr,*cnsr,*cati,*cats;
  double *ci, *cs, *cp, *cjsr, *cnsr, *cati;
  double *cscp1, *cscp2, *Itr;
#ifdef ___DETERMINISTIC
  double *c1, *c2, *i1ca, *i1ba, *i2ca, *i2ba, *fryr1, *fryr2, *fryr3;
#else
  int *y, *ryr1, *ryr2, *ryr3, *nryr;
#endif

  double * Tropc_vec, * CaM_vec, * Myosin_Ca_vec, * Myosin_Mg_vec, * SRB_vec; // buffers

  // *) xw: add the SERCA and Ir matrix for output
  double *j_serca, *j_ryr;
  // *) xw: add output of ica and NCX (default cleft area only) of each CRU
  double *ica_array, *ncx_array, *ica_array_sl, *ica_array_ci;
  // *) xw: set the array for outputting icabk and ipca
  //  *) xw: if add the ncx in the sl, output it for each CRU and each milisecond
#ifdef ___NCX
  double *ncx_array_sl;
#endif
  double *icabk_array, *jpca_array;
  // *) xw: declare a flag for open and close currents in cleft area and sacrolemmal area
  // *) xw: 1 means CRU is coupled with tubule, 0 means CRU is uncupled with tubule

  double *vp;
  double *Jmaxx;

// #ifndef __UNIFORM_TUBULE
  double *tubule_flag;
  double hp_ryr;
// #endif

  int *crupos;


  int *CRU_type;

  void pace(double v, double nai);
  void pace_new(double v, double nai);
  //set Nerst parameters

  double iupave, icaave, incxave, irave, ileakave, icabkave, islcapave;
  double ica_stan, incx_stan; //  *) xw: converted I [pA/pF] from sumica and sumncx [uM/ms]
  double icabk_stan, ipca_stan; //  *) xw: converted I [pA/pF] from sumjcabk and sumjslcap [uM/ms]
  double ir_ss, ir_ct;
  double compute_avg_ci(void);
  double compute_avg_cs(void);
  double compute_avg_cp(void);
  double compute_avg_cnsr(void);
  double compute_avg_cjsr(void);



  double ncx_scale = 1.0;
  double taup_scale = 1.0;
  double ICaL_scale = 1.0;
  double ICaT_scale = 1.0;

  double RyR_release_scale = 1.0;
  // double BkCa_alpha, Cap_alpha, alpha_CaT;
  double BkCa_alpha = 1 * 0.11;
  double Cap_alpha = 1 * 0.11;
  double alpha_CaT = 1 * 0.5;
  double gamma_CaT = 0;  // gamma in cytosol; SL - beta = 1 - alpha - gamma;  alpha in cleft space

  double NCXalpha= 0.5;
  double NCX_gamma = 0.0;
  double SERCA_scale = 1.0;

  int num_open_ICaL;


  bool CLAMP_Cai = false;


  CSubcell& operator=(const CSubcell& sc);//  the type of operator function and sc are both CSubcell

#ifndef __UNIFORM_TUBULE
  void sethyryr(int newhp_ryr) {hp_ryr = newhp_ryr;};
  double gethyryr(void) {return hp_ryr;};
#endif

  void setlayer(int newlayer) {layer = newlayer;};
  int getlayer(void) {return layer;}; // *) xw: add the layer setting functon
  double getsvr(void) {return svr;}; // *) xw: return the surface volume ratio value
  void setdt(double newdt) {dt = newdt;}; //set the time step as 0.01, in pace.cc
  double convertica(void) {return ica_stan;} // *) xw: return the standard Ica
  double convertincx(void) {return incx_stan;} // *) xw: return the standard NCX flux
  double converticabk(void) {return icabk_stan;} // *) xw: return the standard icabk
  double convertipca(void) {return ipca_stan;} // *) xw: return the standard icabk
  void set_Ttubule(int id, double tubule) {
    tubule_flag[id] = tubule;
  }; //*)xw: set the tubule value, 1 -> tubule, 0 -> no tubule
  double getdt(void) {return dt;};
  void setgca(double newgca) {gca = newgca;};
  double getgca(void) {return gca;};
  void setgncx(double newgncx) {vnaca = newgncx;};
  double getgncx(void) {return vnaca;};
  void setJmax(double newJmax);
  double getJmax(void) {return Jmax;};
  void setvup(double newvup) {vup *= newvup;};
  // void setvup(double newvup){vup=newvup;};
  double getvup(void) {return vup;};
  void setkup(double newkup) {kup = newkup;};
  double getkup(void) {return kup;};
  void setKNSR(double newKNSR) {KNSR = newKNSR;};
  double getKNSR(void) {return KNSR;};

  void setgleak(double newgleak) {gleak = newgleak;};
  double getgleak(void) {return gleak;};
  void setcao(double newcao) {cao = newcao;};
  double getcao(void) {return cao;};
  void setBCSQN(double newBCSQN) {BCSQN = newBCSQN;};
  double getBCSQN(void) {return BCSQN;};

  void setgcabk(double newgcabk) {gcabk = newgcabk;};
  double getgcabk(void) {return gcabk;};
  void setqslcap(double newqslcap) {qslcap = newqslcap;};
  double getqslcap(void) {return qslcap;};

  void setKc(double newKc) {Kc = newKc;};
  double getKc(void) {return Kc;};
  void setnM(double newnM) {nM = newnM;};
  double getnM(void) {return nM;};
  void setnD(double newnD) {nD = newnD;};
  double getnD(void) {return nD;};
  void sethh(double newhh) {hh = newhh;};
  double gethh(void) {return hh;};
  void setKK(double newKK) {KK = newKK;};
  double getKK(void) {return KK;};
  void setrhoinf(double newrhoinf) {rhoinf = newrhoinf;};
  double getrhoinf(void) {return rhoinf;};
  void setKu(double newKu) {Ku = newKu;};
  double getKu(void) {return Ku;};
  void setKb(double newKb) {Kb = newKb;};
  double getKb(void) {return Kb;};
  void settauu(double newtauu) {tauu = newtauu;};
  double gettauu(void) {return tauu;};
  void settaub(double newtaub) {taub = newtaub;};
  double gettaub(void) {return taub;};
  void settauc1(double newtauc1) {tauc1 = newtauc1;};
  double gettauc1(void) {return tauc1;};
  void settauc2(double newtauc2) {tauc2 = newtauc2;};
  double gettauc2(void) {return tauc2;};
  void setBCSQN0(double newBCSQN0) {BCSQN0 = newBCSQN0;};
  double getBCSQN0(void) {return BCSQN0;};

  void setvi(double newvi) {vi = newvi;}
  double getvi(void) {return vi;}
  void setvs(double newvs) {
    vs = newvs;
    for (int id = 0; id < n; id++) {
      cscp2[crupos[id]] = 1 / taup * vp[id] / vs;
    }
  }
  double getvs(void) {return vs;}

  double gettausi(void) {return tausi;}
  double settausi(double  newval) {tausi = newval; return tausi;}

  double gettautr(void) {return tautr;}
  double settautr(double  newval) {tautr = newval; return tautr;}






  //#ifdef ___SIGMOID
  void setKcp(double newKcp) {Kcp = newKcp;};
  double getKcp(void) {return Kcp;};
  void setpedk12(double newpedk12) {pedk12 = newpedk12;};
  double getpedk12(void) {return pedk12;};
  void setpedk43(double newpedk43) {pedk43 = newpedk43;};
  double getpedk43(void) {return pedk43;};
  //#endif

  double getMaxSR(void) {return MaxSR;};
  void setMaxSR(double newMaxSR) {MaxSR = newMaxSR;};
  double getMinSR(void) {return MinSR;};
  void setMinSR(double newMinSR) {MinSR = newMinSR;};
  double getec50SR(void) {return ec50SR;};
  void setec50SR(double newec50SR) {ec50SR = newec50SR;};
  double gethkosrca(void) {return hkosrca;};
  void sethkosrca(double newhkosrca) {hkosrca = newhkosrca;};


  //functions: diffusion
  double gettausT(void) {return tausT;}
  double settausT(double  newval) {tausT = newval; taumninv = 2 / tausL + 4 / tausT; return tausT;}
  double gettausL(void) {return tausL;}
  double settausL(double  newval) {tausL = newval; taumninv = 2 / tausL + 4 / tausT; return tausL;}
  double gettauiT(void) {return tauiT;}
  // double settauiT(double  newval){tauiT=newval;return tauiT;}
  double settauiT(double  newval) {tauiT /= newval; return tauiT;}
  double gettauiL(void) {return tauiL;}
  // double settauiL(double  newval){tauiL=newval;return tauiL;}
  double settauiL(double  newval) {tauiL /= newval; return tauiL;}
  double gettaunsrT(void) {return taunsrT;}
  double settaunsrT(double  newval) {taunsrT = newval; return taunsrT;}
  double gettaunsrL(void) {return taunsrL;}
  double settaunsrL(double  newval) {taunsrL = newval; return taunsrL;}

  double calcjcalave(void) {return icaave;};
  double calcjnacaave(void) {return incxave;};
  double calcicabkave(void) {return icabkave;};
  double calcislcapave(void) {return islcapave;};
  void srand(int sed = 0);
  void setboundary(int bcc = 3);
  void resetBuffer(void);




  void check_crupos() {

    for (int i = 0; i < n; ++i)
    {
      std::cout << i << " " << crupos[i] << std::endl;
    }

  }


  void set_lateral_Ttubule();
  void output_Ttubule_map();


  double ICaT_tot;
  iCaT icat_class;
  int NUM_Ttubule;

  void update_single_CRU(double V);


  double get_jSR_inst_buffering(double);
  double get_cytosol_cai_inst_buffering(double Cai);
  double get_cleft_caj_inst_buffering(double Caj);
  double get_submem_casl_inst_buffering(double casl);
  double calculate_dynamic_buffer_cytosol(int id, double Cai, double dt);   // unit of buffers: mM

template<typename U>
void output_map( U *map, const char outputfile[]);
void set_CRU_type();





  // static const double cao = 1.8; // [mM]
  // static const double ko = 5.4;
  // static const double nao = 136; //  [mM]
  static constexpr double F = 96.5; // [C/mmol]
  static constexpr double R = 8.314; // [J/mol*K]
  static constexpr double T = 308; // [K]
  static constexpr double rtf = R * T / F; //~26.5
  static constexpr double rtf2 = R * T / (2 * F);




  // BT, kon, koff, konEGTA, cati[id]

#ifdef ___EGTA
  double *caEGTAi, *caEGTAs;
  double setBEGTA(double newBEGTA) {BEGTA = newBEGTA; resetBuffer(); return BEGTA;}
  double getBEGTA(void) {return BEGTA;}
#endif

  

#ifdef ___DEBUG
  bool bSTOP;
#endif

};




// Definition of template function should be in .h or .hpp file


template<typename U>
void CSubcell::output_map(U *map, const char outputfile[]) {

  ofstream out_ci(outputfile);

  if (map == NULL) {
    std::cerr << "map==NULL in CSubcell::output_map()" << std::endl;
    std::exit(0);
  }
  out_ci <<  "# vtk DataFile Version 3.0" << std::endl;
  out_ci <<  "vtk output"  << std::endl;
  out_ci <<  "ASCII"  << std::endl;
  out_ci <<  "DATASET STRUCTURED_POINTS"  << std::endl;
  out_ci <<  "DIMENSIONS " << nx << " " << ny << " " << nz   << std::endl;
  out_ci <<  "SPACING 1 1 1"  << std::endl;
  out_ci <<  "ORIGIN 0 0 0"  << std::endl;
  out_ci <<  "POINT_DATA " <<  nx*ny*nz << std::endl;
  out_ci <<  "SCALARS HumanAtrium float 1"  << std::endl;
  out_ci <<  "LOOKUP_TABLE default"   << std::endl;

  for (int k = 0; k < nz; ++k)
    for (int i = 0; i < ny; ++i)
    {
      for (int j = 0; j < nx; ++j)
      {
        int id = j + i * nx + k * nx * ny;
        out_ci << map[id] << "\t";
      }
      out_ci << std::endl;

    }
  out_ci.close();
}


#endif