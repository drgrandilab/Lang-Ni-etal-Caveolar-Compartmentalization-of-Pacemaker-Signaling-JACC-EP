#include "subcell.hpp"
#include "ExplicitSolver.hpp"


#ifdef ___DEBUG
#include <iomanip>
#endif

//cell parameters
const double CSubcell::vjsr = 0.02 * 0.5/**0.5*0.5*/; // [um^3]  // 19:43:48, Thu, 13-December-2018, By Haibo times 0.5
const double CSubcell::vp_ave = 0.00126 * 2; // because the ratio between vp and vi is smaller than the shannon model (increased two fold here); // [um^3]
// const double CSubcell::Cmem = 1.1 * (1e-10); //  [F] *) xw: add Cm from GB model to convert the current and deliver tp "cell.cc"
const double CSubcell::Cmem = 0.25 * (1e-10); //  Haibo modified to 25 pF for SAN of mouse
// const double CSubcell::Cmem=3.1*(1e-10); //  [F]  *) xw: get this from Icabk equation


inline unsigned int xorshift(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
  unsigned int t = (*xx ^ (*xx << 11)); *xx = *yy; *yy = *zz; *zz = *ww;
  return ( *ww = (*ww ^ (*ww >> 19)) ^ (t ^ (t >> 8)) );
}


CSubcell::CSubcell(int sizex, int sizey, int sizez, int fmesh, double xii)
{
  if (fmesh == 1)
    dt = 0.1;
  else if (fmesh == 5)
    dt = 0.005;
  else if (fmesh == 3)
    dt = 0.01;
  else {
    cerr << "fine mesh incorrect !\n";
    exit(1);
  }


  ncx_scale = 1.0;

  nx = sizex;
  ny = sizey;
  nz = sizez;

  NUM_Ttubule = nx * ny * nz;
  finemesh = fmesh;
  xi = xii;
  n = nn = 0;
  layer = 2; // *) xw: the default # layer

  //parameters
  cao = 1.8; // [mM]
  nao = 136.0;  // [mM]
  vup = 0.3; //[uM/ms]  // original
  vup = 0.3 * 4; //[uM/ms] // 16:25:22, Tue, 13-August-2019, By Haibo
  vup = 0.3 * 4; //[uM/ms] // 17:35:32, Thu, 09-January-2020, By Haibo
  // vnaca=21.0;  //*) xw : [uM/ms];When using the Cmem as 310pF, vnaca shoule be 21.0[uM/ms] or 3.155 [A/F]; after using 110pF, it should be 7.452 [uM/ms]
  // vnaca = 7.452 * 4;
  // vnaca=21.0*0.7;  // *) xw: change the maxium NCX flux according to the Atrium -70% INCX;
  // vnaca=1.05*0.7; // *) xw: change the maxium NCX flux according to the Atrium -70% INCX; Then the paper indicate this value should be 1.05,  but it's a typo
#ifdef ___USE_ORG_PARAM
  Jmax = 0.0147;
  Ku = 3.8 * 0.0001;
  Kb = 5 * 0.00001;
  tauu = 125.00;
  taub = 5.0;
  tauc1 = 1.0;
  tauc2 = 1.0;
  hh = 23;
  KK = 850;
  gcabk = 0;
  qslcap = 0;
#else
  Jmax = 0.0147 * 18;
  Ku = 5.0;
  Kb = 0.005;
  tauu = 1250.0;
  taub = 0.5;
  tauc1 = 2.0;
  tauc2 = 0.3;
  hh = 10;
  KK = 1400;
  gcabk = 0.0002513 * 20;
  qslcap = 2.35 * 4;
#endif

  MaxSR = 15;
  MinSR = 1;
  ec50SR = 450;
  hkosrca = 2.5;



  kup = 0.123;
  // kup = 0.123*1.5; // 17:35:03, Thu, 09-January-2020, By Haibo
  KNSR = 1700;

  gca = 1.0;
  // gca=0;
  NoCaL = 4;

  gleak = 1.035 * 0.00001;


  taup = 0.022;
  BCSQN = 400;

  tautr = 5.0;
  // tautr = 40.0;  // 40 from Maltsev and Lakatta 2009

  Kc = 600.0;
  nM = 15;
  nD = 35;
  rhoinf = 5000;

  BCSQN0 = 400;

  //#ifdef ___SIGMOID
  Kcp = 100;
  pedk12 = 0.000001;
  pedk43 = 0.000001;
  //#endif

  seed = 0;
  initialized = false;
  bc = 0;
// #ifdef ___NCX
  // NCXalpha = 0.11;
  NCXalpha = 0.5;
  NCX_gamma = 0.0;

  if (NCXalpha + NCX_gamma > 1) {
    std::cerr << "NCXalpha + NCX_gamma > 1; NCX_gamma = " << NCX_gamma << std::endl;
    std::exit(0);
  }
  // NCXalpha = 0.8;

  BkCa_alpha = 1 * 0.11;
  Cap_alpha = 1 * 0.11;

  // CaT distribution
  // alpha -cleft; beta (SL) = 1 - alpha - gamma; gamma - cytosol
  alpha_CaT = 1 * 0.5;
  gamma_CaT = 0.0;

  if (alpha_CaT + gamma_CaT > 1) {
    std::cerr << "alpha_CaT + gamma_CaT > 1; gamma_CaT = " << gamma_CaT << std::endl;
    std::exit(0);
  }

  // NCXalpha = 1;
// #endif
#ifdef ___DEBUG
  bSTOP = false;
#endif

#ifdef __EGTA_
  BEGTA = 350.0;
#endif

  ica_stan = 0.0;
  incx_stan = 0.0;
  icabk_stan = 0.0;
  ipca_stan = 0.0;

  svr = 1.0;//  // 11:08:06, Fri, 07-December-2018, By Haibo 1.0 / 3.5;

#ifndef __UNIFORM_TUBULE
  hp_ryr = 1.0;
#endif


  num_open_ICaL = 0;
}



void CSubcell::init(double initci, double initcj, int LTCC_alpha_in, int LTCC_gamma_in, int rand_seed)
{
  nxny = nx * ny;
  n = nx * ny * nz;

  finemesh3 = finemesh * finemesh * finemesh;
  nnx = finemesh * nx;
  nny = finemesh * ny;
  nnz = finemesh * nz;
  nnxnny = nnx * nny;
  nn = nnx * nny * nnz;

  //cell parameters
  vi =  1.5 * 0.5 / finemesh3;
  vs = 0.025 / finemesh3;
  vnsr = 0.025 / finemesh3;


  tausT = /*0.2**/xi * 1.42 / (finemesh * finemesh);
  // tausL = xi * 3.4 / (finemesh * finemesh);
  tausL = /*0.2**/xi * 1.42 / (finemesh * finemesh);  // 17:14:34, Sun, 03-May-2020, By Haibo; this would accerate rate
  taumninv = 2 / tausL + 4 / tausT;


  tauiT = xi * 2.93 / (finemesh * finemesh);
  tauiL = xi * 2.32 / (finemesh * finemesh);
  taunsrT = xi * 7.2 / (finemesh * finemesh);
  taunsrL = xi * 24.0 / (finemesh * finemesh);


  // *)xw:test value
  // tauiT=xi*0.6/(finemesh*finemesh);
  // tauiL=xi*0.6/(finemesh*finemesh);
  // taunsrT=xi*15.0/(finemesh*finemesh);
  // taunsrL=xi*15.0/(finemesh*finemesh);
  // %%%%%%%%%%%%%%%%%%%%%

  tausi = 0.1 / (finemesh * finemesh);

  ci            = new double [nn];
  cs            = new double [nn];
  cp            = new double [n];
  cjsr          = new double [n];
  cnsr          = new double [nn];
  cati          = new double [nn];

  Tropc_vec     = new double [nn];
  CaM_vec       = new double [nn];
  Myosin_Ca_vec = new double [nn];
  Myosin_Mg_vec = new double [nn];
  SRB_vec       = new double [nn];


  for (int id = 0; id < nn; ++id)
  {
    Tropc_vec[id]     = 8.773191e-3;  //  // unit of buffers: mM
    CaM_vec[id]       = 2.911916e-4;  //  // unit of buffers: mM
    Myosin_Ca_vec[id] = 1.298754e-3;  //  // unit of buffers: mM
    Myosin_Mg_vec[id] = 1.381982e-1;  //  // unit of buffers: mM
    SRB_vec[id]       = 2.143165e-3;  //  // unit of buffers: mM
  }
  //  cats=new double [nn];

#ifdef ___DETERMINISTIC
  c1 = new double [n];
  c2 = new double [n];
  i1ca = new double [n];
  i1ba = new double [n];
  i2ca = new double [n];
  i2ba = new double [n];

  fryr1 = new double [n];
  fryr2 = new double [n];
  fryr3 = new double [n];
#else
  y = new int [n * NoCaL];
  ryr1 = new int [n];
  ryr2 = new int [n];
  ryr3 = new int [n];
  nryr = new int [n];
#endif

  tubule_flag = new double [n];

  CRU_type = new int [n];

  set_lateral_Ttubule();
  // 18:15:28, Tue, 11-December-2018, By Haibo


  LTCC_alpha = LTCC_alpha_in;
  LTCC_gamma = LTCC_gamma_in;

  LTCC_beta = 8 - (LTCC_alpha + LTCC_gamma);  // NUM of LTCC in SL region

  if (LTCC_alpha + LTCC_gamma > 8 or LTCC_alpha + LTCC_gamma < 0) {
    std::cerr << "LTCC_alpha + LTCC_gamma > 8; LTCC_gamma = " << LTCC_gamma << std::endl;
    std::cerr << "or LTCC_alpha + LTCC_gamma < 0; LTCC_gamma = " << LTCC_gamma << std::endl;
    std::exit(0);
  }


  // debug here.
  // for comparison, each space should have 8 LTCCs to ensure the random numbers are identical (strictly)
  // and need to scale the current of each compartment accordingly, for example, evenly (1/3) split among the 3 spaces
  /*LTCC_alpha = 8;
  LTCC_beta = 8;
  LTCC_gamma = 8;*/
  /*LTCC_alpha = 0;
  LTCC_beta = 4;
  LTCC_gamma = 0;*/

  for (int i = 0; i < n; ++i)
  {
    if (LTCC_alpha > 0)
      ical13_vec.push_back(ical13(LTCC_alpha, i, rand_seed)); // default, seed = 0;
    if (LTCC_beta > 0)
      ical13_vec_sl.push_back(ical13(LTCC_beta, i, rand_seed)); // default, seed = 0;
    if (LTCC_gamma > 0)
      ical13_vec_ci.push_back(ical13(LTCC_gamma, i, rand_seed)); // default, seed = 0;
    // random_num.push_back(xor_rand(0, i));

    if (tubule_flag[i] == 1.0)
      RyR_vec.push_back(RyR(100, i, rand_seed));
    else
      RyR_vec.push_back(RyR(1, i, rand_seed));
  }

  // *) xw: output currents matrix
  j_serca = new double [nn];
  j_ryr = new double [n];

  ica_array = new double [n];
  ica_array_sl = new double [n];
  ica_array_ci = new double [n];


  for (int i = 0; i < n; ++i)
  {
    ica_array[i] = 0;
    ica_array_sl[i] = 0;
    ica_array_ci[i] = 0;
  }
  ncx_array = new double [nn];

  for (int i = 0; i < nn; ++i)
  {
    ncx_array[i] = 0;
  }

  icabk_array = new double [nn];
  jpca_array = new double [nn];

  vp = new double [n];
  // tubule_flag=new double [n];

  Jmaxx = new double [n];
  cscp1 = new double [nn];
  cscp2 = new double [nn];

  Ici = new double [nn];
  Icnsr = new double [nn];
#ifndef __RYR_UNIFORM
  sigma = new double [n];
#endif

#ifdef ___NO_CS_BUFFER
  csmn = new double [nn];
#else
  Ics = new double [nn];
  Idps = new double [nn];
#endif

#ifdef ___EGTA
  caEGTAi = new double [nn];
  caEGTAs = new double [nn];
#endif
  // *) xw: add the output of sacrolemmal area NCX flux
#ifdef ___NCX
  ncx_array_sl = new double [n];

  for (int i = 0; i < n; ++i)
  {
    ncx_array_sl[i] = 0;
  }
#endif

  Itr = new double [nn];

  crupos = new int [n];

  //random number
  xsx = new unsigned int [n];
  xsy = new unsigned int [n];
  xsz = new unsigned int [n];
  xsw = new unsigned int [n];
#pragma ivdep
#pragma vector always
  for (int id = 0; id < n; id++)
  {
    xsx[id] = 123456789 + id + seed;
    xsy[id] = 362436069 + id * 100 + seed * 10;
    xsz[id] = 521288629 + id * 1000 + seed * 100;
    xsw[id] = 88675123 + id * 10000 + seed * 1000;
    for (int i = 0; i < 1000; i++)
      xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]);
  }

  //initial conditions

#pragma ivdep
#pragma vector always
  for (int id = 0; id < nn; id++)
  {
    ci[id] = initci;
    cs[id] = initci;
    cnsr[id] = initcj;
    Icnsr[id] = 0;

#ifdef ___NO_CS_BUFFER
    csmn[id] = 0;
#else
    Ics[id] = 0;
    Idps[id] = 0;
#endif
    // *) xw: j_serca
    j_serca[id] = 0;

    Ici[id] = 0;

    cscp1[id] = 0;
    cscp2[id] = 0;
    Itr[id] = 0;
  }
  resetBuffer();
#pragma ivdep
#pragma vector always
  for (int id = 0; id < n; id++)
  {
    Jmaxx[id] = Jmax;

    cp[id] = initci;
    cjsr[id] = initcj;
#ifdef ___DETERMINISTIC
    c1[id] = 1;
    c2[id] = 0;
    i1ca[id] = 0;
    i1ba[id] = 0;
    i2ca[id] = 0;
    i2ba[id] = 0;

    fryr1[id] = 0.03;
    fryr2[id] = 0;
    fryr3[id] = 0;

#else
    ryr1[id] = 0 + int(5.0 * xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]) / (double)(UINT_MAX));
    ryr2[id] = 0;
    ryr3[id] = 0;
    nryr[id] = 100;
    for (int j = 0; j < NoCaL; j++)
      y[id * NoCaL + j] = 2;
#endif

    // *) xw: j_ryr
    j_ryr[id] = 0;

#ifdef ___UNIFORM
    double r = 0; // *) xw: the vp is always vp_ave
#else
    double r = calcvp(0, 0.3, -0.8, 0.8, id); //Gaussian distribution (0,0.3) range(-0.8~0.8)
#endif
    vp[id] = vp_ave * (1 + r);


#ifndef __RYR_UNIFORM // *)xw: use the r parameter directlly 
    sigma[id] = 1 + calcvp(0, 0.2, -0.8, 0.8, id); // will be used on nryr
    // Jmaxx[id]=Jmax*(sigma[id]);
#endif


    // tubule_flag[id] = 1.0;

    int kk = id / nxny;
    kk = kk * finemesh + finemesh / 2;
    int modi = id % nxny;
    int jj = modi / nx;
    jj = jj * finemesh + finemesh / 2;
    int ii = modi % nx;
    ii = ii * finemesh + finemesh / 2;
    crupos[id] = ii + jj * nnx + kk * nnxnny;
    cscp2[crupos[id]] = 1 / taup * vp[id] / vs;
  }

  iupave = icaave = incxave = irave = ileakave = icabkave = islcapave = 0;
  initialized = true;
}



/* initialise random generator with a different sed */
void CSubcell::srand(int sed)
{
  seed = sed;
#pragma ivdep
#pragma vector always
  for (int id = 0; id < n; id++)
  {
    xsx[id] = 123456789 + id + seed;
    xsy[id] = 362436069 + id * 100 + seed * 10;
    xsz[id] = 521288629 + id * 1000 + seed * 100;
    xsw[id] = 88675123 + id * 10000 + seed * 1000;
  }
  for (int id = 0; id < n; id++)
    for (int i = 0; i < 1000; i++)
      xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]);
}

CSubcell::~CSubcell()
{
  if (initialized)
  {
    delarray();
  }
}


// Thu 25 Dec 2025 04:21:46 PM CST
// verify if all memory is properly handled at the end of simulations.
// Verify with Valgrind if possible.

void CSubcell::delarray(void)
{
  delete [] ci;
  delete [] cs;
  delete [] cp;
  delete [] cjsr;
  delete [] cnsr;

#ifdef ___DETERMINISTIC
  delete [] c1;
  delete [] c2;
  delete [] i1ca;
  delete [] i1ba;
  delete [] i2ca;
  delete [] i2ba;

  delete [] fryr1;
  delete [] fryr2;
  delete [] fryr3;
#else
  delete [] ryr1;
  delete [] ryr2;
  delete [] ryr3;
  delete [] nryr;
  delete [] y;
#endif

  // *) xw: clear the output currents matrix
  delete [] j_serca;
  delete [] j_ryr;
  delete [] ica_array;
  delete [] ica_array_ci;
  delete [] ica_array_sl;
  delete [] ncx_array;
  delete [] icabk_array;
  delete [] jpca_array;
  delete [] tubule_flag;

#ifdef ___NCX
  delete [] ncx_array_sl;
#endif

  delete [] cati;
  //      delete [] cats;

  delete [] Jmaxx;

#ifndef __RYR_UNIFORM
  delete[] sigma;
#endif

  delete [] vp;
  delete [] cscp1;
  delete [] cscp2;
  delete [] Itr;

  delete [] Ici;
  delete [] Icnsr;

#ifdef ___NO_CS_BUFFER
  delete [] csmn;
#else
  delete [] Ics;
  delete [] Idps;
#endif
  delete [] crupos;

  delete [] xsx;
  delete [] xsy;
  delete [] xsz;
  delete [] xsw;


  delete [] Tropc_vec     ;
  delete [] CaM_vec       ;
  delete [] Myosin_Ca_vec ;
  delete [] Myosin_Mg_vec ;
  delete [] SRB_vec       ;

}




// Thu 25 Dec 2025 04:18:22 PM CST; comment:
// THis function is no longer used and thus DEPERCATED since 2019
// please use with care!
// Avoid using operator = to copy instances whenever possible;
// use (smart) pointers instead
CSubcell& CSubcell::operator=(const CSubcell& sc)
{
  if (&sc == this)return (*this);
  if (initialized)
  {
    delarray();
  }
  //constructor
  dt = sc.dt;
  nx = sc.nx;
  ny = sc.ny;
  nz = sc.nz;
  finemesh = sc.finemesh;
  xi = sc.xi;
  // *) xw: set the layer from outside
  layer = sc.layer;

  cao = sc.cao;
  vup = sc.vup;
  kup = sc.kup;
  KNSR = sc.KNSR;
  // vnaca = sc.vnaca;
  Jmax = sc.Jmax;
  gca = sc.gca;

  gcabk = sc.gcabk;
  qslcap = sc.qslcap;
  gleak = sc.gleak;
  BCSQN = sc.BCSQN;

  Kc = sc.Kc;
  nM = sc.nM;
  nD = sc.nD;
  hh = sc.hh;
  KK = sc.KK;
  rhoinf = sc.rhoinf;

  Ku = sc.Ku;
  Kb = sc.Kb;
  tauu = sc.tauu;
  taub = sc.taub;
  tauc1 = sc.tauc1;
  tauc2 = sc.tauc2;
  BCSQN0 = sc.BCSQN0;

  //#ifdef ___SIGMOID
  Kcp = sc.Kcp;
  pedk12 = sc.pedk12;
  pedk43 = sc.pedk43;
  //#endif
  NoCaL = sc.NoCaL;
#ifdef ___NCX
  NCXalpha = sc.NCXalpha;
#endif

  MaxSR = sc.MaxSR;
  MinSR = sc.MinSR;
  ec50SR = sc.ec50SR;
  hkosrca = sc.hkosrca;


  init();

  //initial conditions
#pragma ivdep
#pragma vector always
  for (int id = 0; id < nn; id++)
  {
    ci[id] = sc.ci[id];
    cs[id] = sc.cs[id];
    cati[id] = sc.cati[id];
    //      cats[id]=sc.cats[id];
    cnsr[id] = sc.cnsr[id];
    Icnsr[id] = sc.Icnsr[id];
#ifdef ___NO_CS_BUFFER
    csmn[id] = sc.csmn[id];
#else
    Ics[id] = sc.Ics[id];
    Idps[id] = sc.Idps[id];
#endif
    // *) xw: j_serca and j_ryr and ica and ncx for ouput the data for each milisecond;
    // *) xw: tubule_flag is for cleft area and sacrolemmal area currents opening and closing
    j_serca[id] = sc.j_serca[id];
    ncx_array[id] = sc.ncx_array[id];
    icabk_array[id] = sc.icabk_array[id];
    jpca_array[id] = sc.jpca_array[id];

    cscp1[id] = sc.cscp1[id];
    cscp2[id] = sc.cscp2[id];
    Itr[id] = sc.Itr[id];
  }
#pragma ivdep
#pragma vector always
  for (int id = 0; id < n; id++)
  {
    Jmaxx[id] = sc.Jmaxx[id];
#ifndef __RYR_UNIFORM
    sigma[id] = sc.sigma[id]; // *)xw: define sc.sigma
#endif
    cp[id] = sc.cp[id];
    cjsr[id] = sc.cjsr[id];
#ifdef ___DETERMINISTIC
    c1[id] = sc.c1[id];
    c2[id] = sc.c2[id];
    i1ca[id] = sc.i1ca[id];
    i1ba[id] = sc.i1ba[id];
    i2ca[id] = sc.i2ca[id];
    i2ba[id] = sc.i2ba[id];

    fryr1[id] = sc.fryr1[id];
    fryr2[id] = sc.fryr2[id];
    fryr3[id] = sc.fryr3[id];
#else
    ryr1[id] = sc.ryr1[id];
    ryr2[id] = sc.ryr2[id];
    ryr3[id] = sc.ryr3[id];
    nryr[id] = sc.nryr[id];
    for (int j = 0; j < NoCaL; j++)
      y[id * NoCaL + j] = sc.y[id * NoCaL + j];
#endif

    // *) xw: j_serca and j_ryr and ica and ncx for ouput the data for each milisecond;
    // *) xw: tubule_flag is for cleft area and sacrolemmal area currents opening and closing
    j_ryr[id] = sc.j_ryr[id];
    ica_array[id] = sc.ica_array[id];

#ifndef __UNIFORM_TUBULE
    tubule_flag[id] = sc.tubule_flag[id];
    hp_ryr = sc.hp_ryr;
#endif

#ifdef ___NCX
    ncx_array_sl[id] = sc.ncx_array_sl[id];  // updated [n] -> [id]; based on reviewer's comment; Thu 25 Dec 2025 04:17:29 PM CST
#endif

    vp[id] = sc.vp[id];
    cscp2[crupos[id]] = sc.cscp2[crupos[id]];
    Ici[id] = sc.Ici[id];
  }

#pragma ivdep
#pragma vector always
  for (int id = 0; id < n; id++)
  {
    xsx[id] = sc.xsx[id];
    xsy[id] = sc.xsy[id];
    xsz[id] = sc.xsz[id];
    xsw[id] = sc.xsw[id];
  }
  return (*this);
}





void CSubcell::pace(double v, double nai)
{
  iupave = icaave = incxave = irave = ileakave = icabkave = islcapave = 0;

#ifndef ___NO_DIFFUSION
  //set diffusion terms
  computeIci();//diffusion ci
  computeIcnsr();//diffusion cnsr
#ifdef ___NO_CS_BUFFER
  computecsmn();//diffusion cs
#else
  computeIcs();//diffusion cs
#endif
#endif

  LTCC_unitary  ICaL_unitary(v);

  double sumica = 0;

  double sumir = 0;

  inaca NCX(nai, nao, v, cao); // this also calculates  NCX (ca independent part)
  int total_ICaL_open = 0;
  double sum_jnaca_j_flux = 0;
  double sum_j_jcabk = 0;
  double sum_j_jslcap = 0;

  ICaT_tot = ICaT_scale * icat_class.compute_ICaT(v, dt);  //// 9:34:09, Wed, 29-April-2020, By Haibo
  double ICaT_per_CRU = ICaT_tot / NUM_Ttubule;



  #pragma omp parallel for reduction(+: total_ICaL_open) schedule(auto)

  for (int id = 0; id < NUM_Ttubule; id++)
  {
    int NL = 0;

    int id_tt = tubule_ID_vec[id];
    // if (tubule_flag[id_tt] == 1.0) {
    if (LTCC_alpha > 0) {
      NL = ical13_vec[id_tt].update_ical13_stochastic_version_2(dt, v, cp[id_tt]);
      // // 9:35:13, Wed, 29-April-2020, By Haibo
      ica_array[id_tt] = ICaL_scale * gca * ICaL_unitary.compute_LTCC_unitary(cp[id_tt], cao) * NL;  // per Vp volume
      total_ICaL_open += NL;
    } else {
      ica_array[id_tt] = 0.0;
    }

    if (LTCC_beta > 0) {
      NL = ical13_vec_sl[id_tt].update_ical13_stochastic_version_2(dt, v, cs[id_tt]);
      // // 9:35:13, Wed, 29-April-2020, By Haibo
      ica_array_sl[id_tt] = ICaL_scale * gca * ICaL_unitary.compute_LTCC_unitary(cs[id_tt], cao) * NL;
      ica_array_sl[id_tt] = ica_array_sl[id_tt] * vp[id_tt] / vs;  // per (cleft) Vp volume converted to SL volume (vs)
      total_ICaL_open += NL;
    } else {
      ica_array_sl[id_tt] = 0.0;
    }


    if (LTCC_gamma > 0) {
      NL = ical13_vec_ci[id_tt].update_ical13_stochastic_version_2(dt, v, ci[id_tt]);
      // // 9:35:13, Wed, 29-April-2020, By Haibo
      ica_array_ci[id_tt] = ICaL_scale * gca * ICaL_unitary.compute_LTCC_unitary(ci[id_tt], cao) * NL;
      ica_array_ci[id_tt] = ica_array_ci[id_tt] * vp[id_tt] / vi;  // per (cleft) Vp volume converted to cytosolic volume (vi)
      total_ICaL_open += NL;
    } else {
      ica_array_ci[id_tt] = 0.0;
    }

    // ncx_array_sl [id_tt] = tubule_flag[id_tt] * NCXalpha * NCX.compute_NCX(nai, nao, cp[id_tt], cao)  * vs / vp[id_tt];

    // }
  }

  num_open_ICaL = total_ICaL_open;


  #pragma omp parallel for reduction(+: sumica, sumir,sum_jnaca_j_flux,sum_j_jcabk,sum_j_jslcap) schedule(auto)
// #pragma ivdep
// #pragma vector always
  for (int id = 0; id < n; id++)
  {
    //L-double Ca current

    double Ica = ica_array[id];

    // sumica+=Ica;
    // *) xw: Ica [umol/L/ms], here the vp varies around vp_ave, so use the Ica*vp[id] [10^15*umol/ms]
    sumica += Ica * vp[id];  // calculate flux, vp may be heterogeneous

    // *) xw: output the Ica [uM/ms] for each CRU and each milisecond
    // ica_array[id] = Ica; // *)xw: [uM/ms]

    //release current Ir
#ifdef ___DETERMINISTIC
    double Po = fryr2[id] + fryr3[id];
#else
    // double Po = (RyR_vec[id].RyR_2 + RyR_vec[id].RyR_3) / 100.0;
    double Po = RyR_vec[id].Update_RyR_stochastic(dt, cp[id], cjsr[id]);
    // double Po = (ryr2[id] + ryr3[id]) / 100.0;
#endif


    double Ir = RyR_release_scale * RyR_vec[id].Jmax * Po * (cjsr[id] - cp[id]) / vp[id];
    sumir += Ir;
    j_ryr[id] = Ir; // [uM/ms]
    //Diffusion from proximal space to submembrane space Idps
    double kr = RyR_release_scale * RyR_vec[id].Jmax * Po / vp[id];

#ifdef ___NCX




    double jnaca = ncx_scale * tubule_flag[id] * NCXalpha * NCX.compute_NCX_version_2(v, nai, nao, cp[id], cao)  * vs / vp[id];  // per vs volume change to vp (dyadic space)

    // double jnaca = ncx_array_sl[id]; // *)xw: add the ouput of ncx sl
    sum_jnaca_j_flux += jnaca * vp[id];


    // -------Icabk (bac3*kground SL Ca flux) following Shannon --------------
    double eca = rtf2 * log(cao * 1000 / cp[id]);


    double icabk = gcabk * (v - eca); // current [A/F]  negative current
    //Icabk[A/F]*(Cm[pF]/(65*27*11*fine3))/(z*F[C/mol]*(vs*10-9/fine3)[ul])=Icabk*3.3286  -> *) xw: Cm ~ 310pF, larger than 110pF
    // double jcabk=icabk*3.3286;//coming in
    double jcabk = tubule_flag[id] * (BkCa_alpha) * icabk * Cmem * (1e10) * vs / vp[id]; //  *) xw: change the Cmem to GB model

    // -------Islcap (SL Ca pump) following Shannon --------------
    const double vmax = 2.2 * 0.01;
    const double km = 0.5;
    const double h = 1.6;

    // change tubule_flag to something else if finemesh >1 // 18:57:54, Wed, 12-December-2018, By Haibo
    double islcap = tubule_flag[id] * qslcap * vmax / (1 + pow(km / cp[id], h)); //
    // double jslcap=islcap*3.3286;//going out
    double jslcap = tubule_flag[id] * Cap_alpha * islcap *  Cmem * (1e10) * vs / vp[id]; //  *) xw: change the Cmem to GB model

    sum_j_jcabk += jcabk * vp[id];
    sum_j_jslcap += jslcap * vp[id];

    //Leak current Ileak
    const double KJSR = 500;
    double cjsr2 = cjsr[id] * cjsr[id];
    double Ileak = gleak * (cjsr2 * (cjsr[id] - cp[id])) / (cjsr2 + KJSR * KJSR);


    double junc_CaT = tubule_flag[id] * (alpha_CaT) * ICaT_per_CRU * Cmem  / (2.0 * F) / (1e-15 * vp[id]);


    double newcp = (cs[crupos[id]] + taup * (kr * cjsr[id] - Ica - junc_CaT + jnaca - jcabk - jslcap)) / (1 + taup * Jmaxx[id] * Po / vp[id]);
#else
    double newcp = (cs[crupos[id]] + taup * (kr * cjsr[id] - Ica)) / (1 + taup * Jmaxx[id] * Po / vp[id]);
#endif
    if (newcp <= 0)newcp = 0.00001;

    //Diffusion from NSR to JSR
    Itr[crupos[id]] = (cnsr[crupos[id]] - cjsr[id]) / tautr;

    //Nearest-nighbor diffusive current Ici, Ics, IcNSR
    //RyR


    //double m=MM*(nM*BCSQN-cB);
    //      double rhocp=rhoinf*pow(cp[id],hh)/(pow(KK,hh)+pow(cp[id],hh));
    //      double MMcp=(sqrt(1+8*rhocp*BCSQN)-1)/(4*rhocp*BCSQN);

    // RyR_vec[id].Update_RyR_stochastic(dt, cp[id], cjsr[id]);
    //update
    double dcjsr = get_jSR_inst_buffering(cjsr[id]) * (Itr[crupos[id]] - Ir * (vp[id] / vjsr) - Ileak * (vi / vjsr));

#ifdef ___NO_CS_BUFFER
    cscp2[crupos[id]] = vp[id] / (vs * taup);
    cscp1[crupos[id]] = cp[id] * cscp2[crupos[id]];
    //      cscp1[crupos[id]]=cp[id]*vp[id]/vs/taup;
    //      cscp2[crupos[id]]=1/taup*vp[id]/vs;
#else
    Idps[crupos[id]] = (cp[id] - cs[crupos[id]]) / taup;
#endif





#ifdef ___CPDIFF




    //Instantaneous buffering functions
    /*    const double KCAM = 7.0;
        const double BCAM = 24.0;
        const double KSR = 0.6;
        const double BSR = 47.0;
        const double KMCa = 0.033;
        const double BMCa = 140.0;
        const double KMMg = 3.64;
        const double BMMg = 140.0;
        const double KSLH = 0.3;
        const double BSLH = 13.4;

        double CAM = BCAM * KCAM / ((cp[id] + KCAM) * (cp[id] + KCAM));
        double SR = BSR * KSR / ((cp[id] + KSR) * (cp[id] + KSR));
        double MCa = BMCa * KMCa / ((cp[id] + KMCa) * (cp[id] + KMCa));
        double MMg = BMMg * KMMg / ((cp[id] + KMMg) * (cp[id] + KMMg));
        double SLH = BSLH * KSLH / ((cp[id] + KSLH) * (cp[id] + KSLH)); // only for cs
        double Betap = 1 / (1 + CAM + SR + MCa + MMg + SLH);*/


#ifdef ___DEBUG
    if (isnan(Idps[crupos[id]])) //(Idsi != Idsi)
    {
      cout << setprecision(10) << id << "\t" << Idps[crupos[id]] << "\t cp=" << cp[id] << "\t cs=" << cs[crupos[id]] << endl;
      bSTOP = true;
    }
#endif

#ifdef ___NCX
    double dcp = get_cleft_caj_inst_buffering(cp[id])  * (Ir - Ica - junc_CaT + jnaca - jcabk - jslcap - Idps[crupos[id]] + Ileak * (vi / vp[id]));
#else
    double dcp = get_cleft_caj_inst_buffering(cp[id])  * (Ir - Ica - Idps[crupos[id]]);
#endif
    if (not CLAMP_Cai)
      cp[id] += dcp * dt;
#else
    if (not CLAMP_Cai)
      cp[id] = newcp;
#endif

    if (not CLAMP_Cai)
      cjsr[id] += dcjsr * dt;

  }


  double sumjup = 0;
  double sumjleak = 0;
  double sum_jnaca_sl = 0;
  double sum_jnaca_ci = 0;
  double sumjcabk = 0;
  double sumjslcap = 0;

  double sumica_sl = 0;
  double sumica_ci = 0;


  #pragma omp parallel for reduction(+: sumjup, sumjleak,sum_jnaca_sl,sumjcabk,sumjslcap,sum_jnaca_ci, sumica_ci, sumica_sl) schedule(auto)
// #pragma ivdep
// #pragma vector always
  for (int id = 0; id < nn; id++)
  {
    //SERCA Uptake current Iup
    const double H = 1.787;
    double Iup = SERCA_scale * vup * (pow(ci[id] / kup, H) - pow(cnsr[id] / KNSR, H)) / (1 + pow(ci[id] / kup, H) + pow(cnsr[id] / KNSR, H));

    j_serca[id] = Iup; //  [uM/ms]s

    //Leak current Ileak
    const double KJSR = 500;
    double cjsr2 = cnsr[id] * cnsr[id];
    double Ileak = 0;//  0*gleak * (cjsr2 * (cnsr[id] - ci[id])) / (cjsr2 + KJSR * KJSR);

    double j_CaT = tubule_flag[id] * (1 - alpha_CaT - gamma_CaT) * ICaT_per_CRU * Cmem  / (2.0 * F) / (1e-15 * vs);
    double j_CaT_ci = tubule_flag[id] * (gamma_CaT) * ICaT_per_CRU * Cmem  / (2.0 * F) / (1e-15 * vi);


    // double jnaca = tubule_flag[id] * (1 - NCXalpha) * NCX.compute_NCX(nai, nao, cs[id], cao);
    double jnaca = ncx_scale * tubule_flag[id] * (1 - NCXalpha - NCX_gamma) * NCX.compute_NCX_version_2(v, nai, nao, cs[id], cao);   // per cs volume
    double jnaca_ci = ncx_scale * tubule_flag[id] * NCX_gamma * NCX.compute_NCX_version_2(v, nai, nao, ci[id], cao) * vs / vi; // per cs volume convert to ci volume (cytosolic)


    // ICa calculated from previous steps,
    double ICa_sl = ica_array_sl[id];  // per SL volume
    double ICa_ci = ica_array_ci[id];  // per Cytosol volume



    ncx_array[id] = jnaca;  // need new array for jnanca_ci;
    // -------Icabk (bac3*kground SL Ca flux) following Shannon --------------
    double eca = rtf2 * log(cao * 1000 / cs[id]);


    double icabk = gcabk * (v - eca); // current [A/F]  negative current
    //Icabk[A/F]*(Cm[pF]/(65*27*11*fine3))/(z*F[C/mol]*(vs*10-9/fine3)[ul])=Icabk*3.3286  -> *) xw: Cm ~ 310pF, larger than 110pF
    // double jcabk=icabk*3.3286;//coming in
    double jcabk = tubule_flag[id] * (1 - BkCa_alpha) * icabk * Cmem * (1e10); //  *) xw: change the Cmem to GB model
    icabk_array[id] = icabk; // *) xw: output [A/F]
    // -------Islcap (SL Ca pump) following Shannon --------------
    const double vmax = 2.2 * 0.01;
    const double km = 0.5;
    const double h = 1.6;

    // change tubule_flag to something else if finemesh >1 // 18:57:54, Wed, 12-December-2018, By Haibo
    double islcap = tubule_flag[id] * qslcap * vmax / (1 + pow(km / cs[id], h)); //
    // double jslcap=islcap*3.3286;//going out
    double jslcap = tubule_flag[id] * (1 - Cap_alpha) * islcap *  Cmem * (1e10); //  *) xw: change the Cmem to GB model
    jpca_array[id] = jslcap; // *) xw: output [uM/ms]


    //Instantaneous buffering functions
    /* const double KCAM = 7.0;
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

     double CAM = BCAM * KCAM / ((ci[id] + KCAM) * (ci[id] + KCAM));
     double SR = BSR * KSR / ((ci[id] + KSR) * (ci[id] + KSR));
     double MCa = BMCa * KMCa / ((ci[id] + KMCa) * (ci[id] + KMCa));
     double MMg = BMMg * KMMg / ((ci[id] + KMMg) * (ci[id] + KMMg));
     double Betai = 1 / (1 + CAM + SR + MCa + MMg);*/

    const double BT = 70.0;
    const double kon = 0.0327;
    const double koff = 0.0196;

    double ITCi = kon * ci[id] * (BT - cati[id]) - koff * cati[id];

#ifdef ___EGTA
    const double konEGTA = 4E-3;
    const double koffEGTA = 2E-3;
    double IEGTAi = konEGTA * ci[id] * (BEGTA - caEGTAi[id]) - koffEGTA * caEGTAi[id];
    double IEGTAs = konEGTA * cs[id] * (BEGTA - caEGTAs[id]) - koffEGTA * caEGTAs[id];
#endif



    // double ITCs=kon*cs[id]*(BT-cats[id])-koff*cats[id];
    double ITCs = 0;

    //Diffusion from submembrane to myoplasm Idsi
    double Idsi = (cs[id] - ci[id]) / tausi;
#ifdef ___DEBUG
    if (isnan(Idsi)) //(Idsi != Idsi)
    {
      cout << setprecision(10) << id << "\t" << Idsi << "\t cs=" << cs[id] << "\t ci=" << ci[id] << endl;
      bSTOP = true;
    }
#endif


#ifdef ___NO_CS_BUFFER
#ifdef ___NO_DIFFUSION
    double newcs = (cscp1[id] + jnaca + ci[id] / tausi - ITCs + csmn[id] - jcabk - jslcap - j_CaT - ICa_sl) / (cscp2[id] + 1 / tausi);
#else
    double newcs = (cscp1[id] + jnaca + ci[id] / tausi - ITCs + csmn[id] - jcabk - jslcap - j_CaT - ICa_sl) / (cscp2[id] + 1 / tausi + taumninv);
#endif
    if (newcs <= 0)newcs = 0.00001; // to avoid div 0
#else


    /*    CAM = BCAM * KCAM / ((cs[id] + KCAM) * (cs[id] + KCAM));
        SR = BSR * KSR / ((cs[id] + KSR) * (cs[id] + KSR));
        MCa = BMCa * KMCa / ((cs[id] + KMCa) * (cs[id] + KMCa));
        MMg = BMMg * KMMg / ((cs[id] + KMMg) * (cs[id] + KMMg));
        double SLH = BSLH * KSLH / ((cs[id] + KSLH) * (cs[id] + KSLH)); // only for cs
        double Betas = 1 / (1 + CAM + SR + MCa + MMg + SLH);*/

#ifdef ___EGTA
    double dcs = get_submem_casl_inst_buffering(cs[id]) * (Idps[id] * vp[id] / vs + jnaca - Idsi - ITCs - IEGTAs + Ics[id] - jcabk - jslcap - j_CaT - ICa_sl);
#else
    double dcs = get_submem_casl_inst_buffering(cs[id]) * (Idps[id] * vp[id] / vs + jnaca - Idsi - ITCs + Ics[id] - jcabk - jslcap - j_CaT - ICa_sl);
#endif
#endif


#ifdef ___EGTA
    double dci = get_cytosol_cai_inst_buffering(ci[id]) * (Idsi * (vs / vi) - Iup + jnaca_ci + Ileak - ITCi - IEGTAi + Ici[id] - j_CaT_ci - ICa_ci);
#else
    double dci = get_cytosol_cai_inst_buffering(ci[id]) * (Idsi * (vs / vi) - Iup + jnaca_ci + Ileak - ITCi + Ici[id] - j_CaT_ci - ICa_ci);

    // double dci = Idsi * (vs / vi) - Iup + Ileak + Ici[id] - 1000 * calculate_dynamic_buffer_cytosol(id, ci[id] / 1000.0, dt);
#endif

    double dcnsr = ((Iup - Ileak) * (vi / vnsr) - Itr[id] * (vjsr / vnsr) + Icnsr[id]);


    if (not CLAMP_Cai) {

      ci[id] += dci * dt;
      cati[id] += ITCi * dt;
#ifdef ___NO_CS_BUFFER
      cs[id] = newcs;
#else
      cs[id] += dcs * dt;
#endif
      //      cats[id]+=ITCs*dt;
      cnsr[id] += dcnsr * dt;

    }
#ifdef ___EGTA
    caEGTAi[id] += IEGTAi * dt;
    caEGTAs[id] += IEGTAs * dt;
#endif

    // *)xw: the sum of flux on surface, because they are 0 inside of the myocyte
    sumjup += Iup;
    sumjleak += Ileak;
    sum_jnaca_sl += jnaca;
    sum_jnaca_ci += jnaca_ci;
    sumjcabk += jcabk;
    sumjslcap += jslcap;

    sumica_sl += ICa_sl;
    sumica_ci += ICa_ci;

  }

  irave = sumir / n;
  iupave = sumjup / nn;

  ica_stan = 0.001 * sumica * (1e-15) * 2 * F * 1000 / Cmem //  [pA/pF], cleft space  sumica has vp accounted for 
             + 0.001 * sumica_sl * vs * (1e-15) * 2 * F * 1000 / Cmem //  [pA/pF], SL space
             + 0.001 * sumica_ci * vi * (1e-15) * 2 * F * 1000 / Cmem; //  [pA/pF], cytosol space

  incx_stan = sum_jnaca_sl * 0.001 * vs * (1e-15) * F * 1000 / Cmem
              + sum_jnaca_ci * 0.001 * vi * (1e-15) * F * 1000 / Cmem
              + sum_jnaca_j_flux * 0.001 * (1e-15) * F * 1000 / Cmem; //  [pA/pF]  sum_jnaca_j_flux has vp considered (vp maybe variable)
  icabk_stan = sumjcabk * 0.001 * vs * (1e-15) * 2 * F * 1000 / Cmem
               + sum_j_jcabk * 0.001 * (1e-15) * 2 * F * 1000 / Cmem; //  [pA/pF]
  ipca_stan = sumjslcap * 0.001 * vs * (1e-15) * 2 * F * 1000 / Cmem
              + sum_j_jslcap * 0.001 * (1e-15) * 2 * F * 1000 / Cmem; //  [pA/pF]
}





// Thu 25 Dec 2025 04:18:22 PM CST; comment:
// THis function is no longer used and thus DEPERCATED since 2019
// please use with care!
// this function simulates the transition of RyR in the model/ stochastically .
int CSubcell::bino(double num, double p, int ii)
{
  int res;
  double lambda = num * p;
  if (lambda > 12)
  {
    //Gaussian
    double x1, x2, w;
    do
    {
      x1 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
      x2 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    w = sqrt((-2.0 * log(w)) / w);
    double y1 = x1 * w;
    //double y2=x2*w;
    res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
    res = int(res + 0.5); //round
  }
  else if (100 * p < 6.6 + 52 * pow(num, double(-0.5)))
  {
    //Poisson
    double L = exp(-lambda);
    double k = 0;
    double pp = 1;
    do
    {
      k++;
      double u = xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX);
      pp *= u;
    } while (pp >= L);
    res = k - 1;
  }
  else
  {
    //Gaussian
    double x1, x2, w;
    do
    {
      x1 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
      x2 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    w = sqrt((-2.0 * log(w)) / w);
    double y1 = x1 * w;
    //double y2=x2*w;
    res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
    res = int(res + 0.5); //round
  }
  if (res < 0)res = 0;

  return res;
}




// Thu 25 Dec 2025 04:18:22 PM CST; comment:
// THis function is no longer used and thus DEPERCATED since 2019
// please use with care!
double CSubcell::calcvp(double mean, double std, double lim1, double lim2, int ii)
{
  double res;
  do {
    //Gaussian
    double x1, x2, w;
    do
    {
      x1 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
      x2 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    w = sqrt((-2.0 * log(w)) / w );
    double y1 = x1 * w;
    //double y2 = x2 * w;
    res = y1 * std + mean;
  } while (res < lim1 || res > lim2);
  return res;
}


double CSubcell::compute_avg_ci(void)
{
  double sum = 0;
  #pragma omp parallel for reduction(+: sum) schedule (auto)
  for (int id = 0; id < nn; id++)
    sum += ci[id];
  return (sum / nn);
}
double CSubcell::compute_avg_cs(void)
{
  double sum = 0;
  #pragma omp parallel for reduction(+: sum) schedule (auto)
  for (int id = 0; id < nn; id++)
    sum +=    cs[id];
  return (sum / nn);
}
double CSubcell::compute_avg_cnsr(void)
{
  double sum = 0;
  #pragma omp parallel for reduction(+: sum) schedule (auto)
  for (int id = 0; id < nn; id++)
    sum += cnsr[id];
  return (sum / nn);
}


double CSubcell::compute_avg_cjsr(void)
{
  double sum = 0;
  #pragma omp parallel for reduction(+: sum) schedule (auto)
  for (int id = 0; id < nn; id++)
    sum += cjsr[id];
  return (sum / nn);
}

double CSubcell::compute_avg_cp(void) {
  double sum = 0;
  #pragma omp parallel for reduction(+: sum) schedule (auto)
  for (int id = 0; id < n; id++)
    sum +=  cp[id];  // assuming that vp[id] is identical  // // 11:21:50, Wed, 19-December-2018, By Haibo
  return (sum / n);
}
void CSubcell::setboundary(int bcc)
{
  if (bc > 0) //corner
  {
    Jmaxx[0 + 0 * nx + 0 * nxny] = 0;
    Jmaxx[0 + (ny - 1)*nx + (nz - 1)*nxny] = 0;
    Jmaxx[0 + (ny - 1)*nx + 0 * nxny] = 0;
    Jmaxx[0 + 0 * nx + (nz - 1)*nxny] = 0;
    Jmaxx[(nx - 1) + (ny - 1)*nx + 0 * nxny] = 0;
    Jmaxx[(nx - 1) + 0 * nx + (nz - 1)*nxny] = 0;
    Jmaxx[(nx - 1) + 0 * nx + 0 * nxny] = 0;
    Jmaxx[(nx - 1) + (ny - 1)*nx + (nz - 1)*nxny] = 0;
  }
  if (bc > 1) //edge
  {
    //x fixed
// #ifdef _OPENMP
    #pragma omp parallel for
// #endif
    for (int j = 1; j < ny - 1; j++)
    {
      Jmaxx[0 + j * nx + 0 * nxny] = 0;
      Jmaxx[(nx - 1) + j * nx + 0 * nxny] = 0;
      Jmaxx[0 + j * nx + (nz - 1)*nxny] = 0;
      Jmaxx[(nx - 1) + j * nx + (nz - 1)*nxny] = 0;
    }
    //y fixed
// #ifdef _OPENMP
    #pragma omp parallel for
// #endif
    for (int i = 1; i < (nx - 1); i++)
    {
      Jmaxx[i + 0 * nx + 0 * nxny] = 0;
      Jmaxx[i + (ny - 1)*nx + 0 * nxny] = 0;
      Jmaxx[i + 0 * nx + (nz - 1)*nxny] = 0;
      Jmaxx[i + (ny - 1)*nx + (nz - 1)*nxny] = 0;
    }
#pragma ivdep
#pragma vector always
    for (int k = 1; k < nz - 1; k++)
    {
      Jmaxx[0 + 0 * nx + k * nxny] = 0;
      Jmaxx[0 + (ny - 1)*nx + k * nxny] = 0;
      Jmaxx[(nx - 1) + 0 * nx + k * nxny] = 0;
      Jmaxx[(nx - 1) + (ny - 1)*nx + k * nxny] = 0;
    }
  }
  if (bc > 2) //surface
  {
    //x fixed
// #ifdef _OPENMP
    #pragma omp parallel for
// #endif
    for (int j = 1; j < ny - 1; j++)
    {
#pragma ivdep
#pragma vector always
      for (int k = 1; k < nz - 1; k++)
      {
        Jmaxx[0 + j * nx + k * nxny] = 0;
        Jmaxx[(nx - 1) + j * nx + k * nxny] = 0;
      }
    }
    //y fixed
// #ifdef _OPENMP
    #pragma omp parallel for
// #endif
    for (int i = 1; i < (nx - 1); i++)
    {
#pragma ivdep
#pragma vector always
      for (int k = 1; k < nz - 1; k++)
      {
        Jmaxx[i + 0 * nx + k * nxny] = 0;
        Jmaxx[i + (ny - 1)*nx + k * nxny] = 0;
      }
      //z fixed
#pragma ivdep
#pragma vector always
      for (int j = 1; j < ny - 1; j++)
      {
        Jmaxx[i + j * nx + 0 * nxny] = 0;
        Jmaxx[i + j * nx + (nz - 1)*nxny] = 0;
      }
    }
  }
  if (bc > 3) //corner more
  {
    Jmaxx[1 + 1 * nx + 1 * nxny] = 0;
    Jmaxx[1 + (ny - 2)*nx + (nz - 2)*nxny] = 0;
    Jmaxx[1 + (ny - 2)*nx + 1 * nxny] = 0;
    Jmaxx[1 + 1 * nx + (nz - 2)*nxny] = 0;
    Jmaxx[(nx - 2) + (ny - 2)*nx + 1 * nxny] = 0;
    Jmaxx[(nx - 2) + 1 * nx + (nz - 2)*nxny] = 0;
    Jmaxx[(nx - 2) + 1 * nx + 1 * nxny] = 0;
    Jmaxx[(nx - 2) + (ny - 2)*nx + (nz - 2)*nxny] = 0;
  }


}
void CSubcell::setJmax(double newJmax) {
  Jmax = newJmax;
  for (int id = 0; id < n; id++)Jmaxx[id] = Jmax;
  setboundary(bc);
}
void CSubcell::resetBuffer(void) {
  const double BT = 70.0;
  const double kon = 0.0327;
  const double koff = 0.0196;

  const double konEGTA = 4E-3;
  const double koffEGTA = 2E-3;

#pragma ivdep
#pragma vector always
  for (int id = 0; id < nn; id++)
  {
    cati[id] = kon * ci[id] * BT / (kon * ci[id] + koff);
    //      cats[id]=kon*cs[id]*BT/(kon*cs[id]+koff);
#ifdef ___EGTA
    caEGTAi[id] = konEGTA * ci[id] * BEGTA / (konEGTA * ci[id] + koffEGTA);
    caEGTAs[id] = konEGTA * cs[id] * BEGTA / (konEGTA * cs[id] + koffEGTA);
#endif
  }
}




void CSubcell::set_lateral_Ttubule() {
  int nxny = nx * ny;
  int nxyz = nx * ny * nz;

  int counter = 0;
  for (int id = 0; id < (nxyz); id++) {
    if ((id % nx < layer) || (id % nx > (nx - 1 - layer))) {
      set_Ttubule(id, 1.0); counter++;
    }
    else if ((((id % nxny - id % nx) / nx) < layer) || (((id % nxny - id % nx) / nx) > (ny - 1 - layer))) {
      set_Ttubule(id, 1.0); counter++;
    }
    else if ((((id - id % nxny) / nxny) < layer) || (((id - id % nxny) / nxny) > (nz - 1 - layer))) {
      set_Ttubule(id, 1.0); counter++;
    }
    else set_Ttubule(id, 0.0);
  }

  NUM_Ttubule = counter;


  for (int id = 0; id < (nxyz); id++) {
    if ( tubule_flag[id] == 1.0)
      tubule_ID_vec.push_back(id);
  }

}



void CSubcell::set_CRU_type() {


  for (int k = 0; k < nz; ++k)
    for (int i = 0; i < ny; ++i)
    {
      for (int j = 0; j < nx; ++j)
      {
        int id = j + i * nx + k * nx * ny;
        // out_ci << tubule_flag[id] << "\t";

        int type = 0;

        if (j == 0 or j == nx - 1 or  i == 0 or i == ny - 1 or k == 0 or k == nz - 1 )
          type = 2;
        else if (j == 1 or  j == nx - 2 or i == 1 or i == ny - 2 or k == 1 or k == nz - 2)
          type = 1;

        CRU_type[id] = type;
      }

    }

}




void CSubcell::output_Ttubule_map() {

  ofstream out_ci("Ttubule_map.vtk");

  if (tubule_flag == NULL) {
    std::cerr << "tubule_flag==NULL" << std::endl;
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
        out_ci << tubule_flag[id] << "\t";
      }
      out_ci << std::endl;

    }
  out_ci.close();

}



