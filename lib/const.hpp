  //// Model Parameters
  // Constants
  double Temp = 310; // 300 or 310 // [K]    delete because from protocol???
  double R = 8314;       // [J/kmol*K]
  double Frdy = 96485;   // [C/mol]
  double FoRT = Frdy/R/Temp;
  double Cmem = 1.1e-10; // [F] membrane capacitance 1.3810e-10;
  double Qpow = (Temp-310)/10;


  // Cell geometry
  double cellLength = 100;     // cell length [um] 113;%100
  double cellRadius = 10.25;   // cell radius [um] 12;%10.25
  double junctionLength = 160e-3;  // junc length [um]
  double junctionRadius = 15e-3;   // junc radius [um]
  double distSLcyto = 0.45;    // dist. SL to cytosol [um]
  double distJuncSL = 0.5;  // dist. junc to SL [um]
  double DcaJuncSL = 1.64e-6;  // Dca junc to SL [cm^2/sec]
  double DcaSLcyto = 1.22e-6; // Dca SL to cyto [cm^2/sec]
  double DnaJuncSL = 1.09e-5;  // Dna junc to SL [cm^2/sec]
  double DnaSLcyto = 1.79e-5;  // Dna SL to cyto [cm^2/sec]
  double Vcell = M_PI*pow((cellRadius),2)*cellLength*1e-15;    // [L]
  double Vmyo = 0.65*Vcell;
  double Vsr = 0.035*Vcell;
  double Vsl = 0.02*Vcell;
  double Vjunc = 1*0.0539*.01*Vcell;
  double SAjunc = 20150*M_PI*2*junctionLength*junctionRadius;  // [um^2]
  double SAsl = M_PI*2*cellRadius*cellLength;          // [um^2]
  //J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;// [L/msec] = 1.1074e-13
  //J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  // [L/msec] = 1.5714e-12
  //J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;// [L/msec] = 7.36e-13
  //J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  // [L/msec] = 2.3056e-11
  //J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;// [L/msec] = 9.9664e-014
  //J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  // [L/msec] = 1.7460e-012
  //J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;// [L/msec] = 6.6240e-013
  //J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  // [L/msec] = 2.5618e-011
  // tau's from c-code, not used here
  double J_ca_juncsl = 1/1.2134e12; // [L/msec] = 8.2413e-13
  double J_ca_slmyo = 1/2.68510e11; // [L/msec] = 3.2743e-12
  double J_na_juncsl = 1/(1.6382e12/3*100); // [L/msec] = 6.1043e-13
  double J_na_slmyo = 1/(1.8308e10/3*100);  // [L/msec] = 5.4621e-11

  // Fractional currents in compartments
  double Fjunc = 0.11;
  double Fsl = 1-Fjunc;
  double Fjunc_CaL = 0.9;
  double Fsl_CaL = 1-Fjunc_CaL;

  // Fixed ion concentrations
  double Cli = 15;   // Intracellular Cl  [mM]
  double Clo = 150;  // Extracellular Cl  [mM]
  double Ko = 5.4;   // Extracellular K   [mM]
  double Nao = 140; // 140 % Extracellular Na  [mM]
  double Cao = 1.8;  // Extracellular Ca  [mM]
  double Mgi = 1;    // Intracellular Mg  [mM]



  //// Na transport parameters
  double GNaB = 1*0.597e-3;    // [mS/uF]
  double IbarNaK = 1*1.26;     // [uA/uF]
  double KmKo = 1.5;         // [mM]1.5
  double Q10NaK = 1.63;
  double Q10KmNai = 1.39;

  //// K current parameters
  double pNaK = 0.01833;
  double gkp = 0.002;

  //// Cl current parameters
  double GClCa = 0.0548;   // [mS/uF]
  double GClB = 9e-3;        // [mS/uF]
  double KdClCa = 100e-3;    // [mM]
  double GClCFTR = 0; //4.9e-3*ISO; // [mS/uF]

  //// Ca transport parameters
  // I_ca parameteres
  double Q10CaL = 1.8;

  // I_cabk parameteres
  double GCaB = 6.0643e-4;   // [uA/uF] 3

  // NCX parameteres
  double KmCai = 3.59e-3;    // [mM]
  double KmCao = 1.3;        // [mM]
  double KmNai = 12.29;      // [mM]
  double KmNao = 87.5;       // [mM]
  double ksat = 0.27;        // [none]
  double nu = 0.35;          // [none]
  double Kdact = 0.384e-3;   // [mM] 0.256 rabbit384
  double Q10NCX = 1.57;      // [none]

  // I_pca parameteres
  double IbarSLCaP = 0.0471; // IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
  double KmPCa = 0.5e-3;     // [mM]
  double Q10SLCaP = 2.35;    // [none]

  // SR flux parameters
  double Q10SRCaP = 2.6;          // [none]
  double Vmax_SRCaP = 5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
  double Kmr = 1.7;               // [mM]L cytosol
  double hillSRCaP = 1.787;       // [mM]
  double ks = 25;                 // [1/ms]
  double kom = 0.06;              // [1/ms]
  double kiCa = 0.5;              // [1/mM/ms]
  double kim = 0.005;             // [1/ms]
  double ec50SR = 0.45;           // [mM]

  //// Buffering parameters
  // koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
  double Bmax_Naj = 7.561;       // [mM] // Na buffering
  double Bmax_Nasl = 1.65;       // [mM]
  double koff_na = 1e-3;         // [1/ms]
  double kon_na = 0.1e-3;        // [1/mM/ms]
  double Bmax_TnClow = 70e-3;    // [mM]                      // TnC low affinity
  double kon_tncl = 32.7;        // [1/mM/ms]
  double Bmax_TnChigh = 140e-3;  // [mM]                      // TnC high affinity
  double koff_tnchca = 0.032e-3; // [1/ms]
  double kon_tnchca = 2.37;      // [1/mM/ms]
  double koff_tnchmg = 3.33e-3;  // [1/ms]
  double kon_tnchmg = 3e-3;      // [1/mM/ms]
  double Bmax_CaM = 24e-3;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
  double koff_cam = 238e-3;      // [1/ms]
  double kon_cam = 34;           // [1/mM/ms]
  double Bmax_myosin = 140e-3;   // [mM]                      // Myosin buffering
  double koff_myoca = 0.46e-3;   // [1/ms]
  double kon_myoca = 13.8;       // [1/mM/ms]
  double koff_myomg = 0.057e-3;  // [1/ms]
  double kon_myomg = 0.0157;     // [1/mM/ms]
  double Bmax_SR = 19*.9e-3;     // [mM] (Bers text says 47e-3) 19e-3
  double koff_sr = 60e-3;        // [1/ms]
  double kon_sr = 100;           // [1/mM/ms]
  double Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl; // [mM]    // SL buffering
  double Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1; // [mM] //Fei *0.1!!! junction reduction factor
  double koff_sll = 1300e-3;     // [1/ms]
  double kon_sll = 100;          // [1/mM/ms]
  double Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl; // [mM]
  double Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1; // [mM] //Fei *0.1!!! junction reduction factor
  double koff_slh = 30e-3;       // [1/ms]
  double kon_slh = 100;          // [1/mM/ms]
  double Bmax_Csqn = 140e-3*Vmyo/Vsr; // [mM] // Bmax_Csqn = 2.6;      // Csqn buffering
  double koff_csqn = 65;         // [1/ms]
  double kon_csqn = 100;        // [1/mM/ms]





  //// I_Na: Voltage-Gated Na Current (NEW MARKOV MODEL)
  // INa Markov Model Parameters /////////////////////////////

  double P2a1=0.1027;
  double P3a1=2.5;
  double P4a1=17;
  double P5a1=0.20;
  double P6a1=150;
  double P4a2=15;
  double P5a2=0.23;
  double P4a3=12;
  double P5a3=0.25;
  double P1b1=0.1917;
  double P2b1=20.3;
  double P1b2=0.2;
  double P2b2=2.5;
  double P1b3=0.22;
  double P2b3=7.5;
  double P1a4=0.188495;
  double P2a4=16.6;
  double P3a4=0.393956;
  double P4a4=7;
  double P1a5=7e-7;
  double P2a5=7.2;
  double P1b5=0.0044;
  double P2b5=2e-5;
  double P1a6=100;
  double P1b6=8.9554e-7;
  double P2b6=11.3944;
  double P1a7=0.487e-4;
  double P2a7=23.2696;
  double P1b7=0.2868e-3;
  double P2b7=35.9898;
  double P1a8=0.1e-7;
  double P1b8=9.8e-3;

  double diffusion = 5500;       // Ranolazine        //drug = 1 * (1E-6); % (M)
  double pH = 7.4;
  double pKa = 7.2;
  double dd = -0.7;
  double kd0 = 100.5 * (1e-6);
  double kd0_b = 1.5012 * (1e-6); // bursting
  double k_off_0 = 400 * (1e-6);
  double ki_off_0 = 5.4 * (1e-6);
  double kc_off_0 = 800 * (1e-6);
  double Pa3_c = 3.6811;
  double Pa4_c = 6.8705e+04;
  double Pa5_c = 4.0832e-02;
  double Pb5_c = 1.7561e-01;
  double Pa6_c = 1*8;
  double Pb6_c = 1/4;
  double Pa7_c = 1;
  double Pb7_c = 1;
  double Pa3_n = 2.3570e+02;
  double Pa4_n = 2.1182e+02;
  double Pb5_n = 1.2197e-03;
  double Pa6_n = 1;
  double Pa7_n = 1;

  // INa Markov Model Transition Rates /////////////////////////////
  double Q10_INa = 2.1;
  double Tfactor_INa = pow(Q10_INa,(Temp-300)/10);


