#include "subcell.hpp"
// #include "ExplicitSolver.hpp"


#ifdef ___DEBUG
#include <iomanip>
#endif

void CSubcell::pace_new(double v, double nai)
{



  iupave = icaave = incxave = irave = ileakave = icabkave = islcapave = 0;

  LTCC_unitary  ICaL_unitary(v);

  double sumica = 0;
  double sumir = 0;

  inaca NCX(nai, nao, v, cao); // this also calculates  NCX (ca independent part)
  int total_ICaL_open = 0;
  double sum_jnaca_j_flux = 0;

  double sumjup = 0;
  double sumjleak = 0;
  double sum_jnaca_sl = 0;
  double sumjcabk = 0;
  double sumjslcap = 0;
  // #ifndef ___NO_DIFFUSION
  //set diffusion terms

  // #endif
  // #endif
  ICaT_tot = icat_class.compute_ICaT(v, dt);
  // }
  double ICaT_per_CRU = ICaT_tot / NUM_Ttubule;

  computeIci();//diffusion ci
  computeIcnsr();//diffusion cnsr
  // #ifdef ___NO_CS_BUFFER
  //   computecsmn();//diffusion cs
  // #else
  computeIcs();//diffusion cs
  #pragma omp parallel  default (shared) 
  {

    #pragma omp  for reduction(+: sumica, sumir,sum_jnaca_j_flux, total_ICaL_open)  /*schedule(dynamic)*/
    // #pragma ivdep
    // #pragma vector always
    for (int id = 0; id < n; id++)
    {
      //L-double Ca current
      int NL = 0;
      double Ica = 0.0;
      if (tubule_flag[id] == 1.0) {
        NL = ical13_vec[id].update_ical13_stochastic_version_2(dt, v, cp[id]);
        Ica = gca * ICaL_unitary.compute_LTCC_unitary(cp[id], cao)  * NL;
      }
      // sumica+=Ica;
      // *) xw: Ica [umol/L/ms], here the vp varies around vp_ave, so use the Ica*vp[id] [10^15*umol/ms]
      sumica += Ica * vp[id];  // calculate flux, vp may be heterogeneous
      total_ICaL_open += NL;
      // *) xw: output the Ica [uM/ms] for each CRU and each milisecond
      ica_array[id] = Ica; // *)xw: [uM/ms]

      //release current Ir
      // #ifdef ___DETERMINISTIC
      // double Po = fryr2[id] + fryr3[id];
      // #else
      // double Po = (RyR_vec[id].RyR_2 + RyR_vec[id].RyR_3) / 100.0;
      double Po = RyR_vec[id].Update_RyR_stochastic(dt, cp[id], cjsr[id]);
      // double Po = (ryr2[id] + ryr3[id]) / 100.0;
      // #endif


      double Ir = RyR_vec[id].Jmax * Po * (cjsr[id] - cp[id]) / vp[id];
      sumir += Ir;
      j_ryr[id] = Ir; // [uM/ms]
      //Diffusion from proximal space to submembrane space Idps
      double kr = RyR_vec[id].Jmax * Po / vp[id];

      // #ifdef ___NCX



      double jnaca = tubule_flag[id] * NCXalpha * NCX.compute_NCX(nai, nao, cp[id], cao)  * vs / vp[id];
      // double jnaca = tubule_flag[id] * NCXalpha * NCX.compute_NCX_version_2(v, nai, nao, cp[id], cao)  * vs / vp[id];

      sum_jnaca_j_flux += jnaca * vp[id];

      ncx_array_sl[id] = jnaca; // *)xw: add the ouput of ncx sl

      double newcp = (cs[crupos[id]] + taup * (kr * cjsr[id] - Ica + jnaca)) / (1 + taup * Jmaxx[id] * Po / vp[id]);
      // #else
      //     double newcp = (cs[crupos[id]] + taup * (kr * cjsr[id] - Ica)) / (1 + taup * Jmaxx[id] * Po / vp[id]);
      // #endif
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
      double dcjsr = get_jSR_inst_buffering(cjsr[id]) * (Itr[crupos[id]] - Ir * (vp[id] / vjsr));

      // #ifdef ___NO_CS_BUFFER
      //     cscp2[crupos[id]] = vp[id] / (vs * taup);
      //     cscp1[crupos[id]] = cp[id] * cscp2[crupos[id]];
      //     //      cscp1[crupos[id]]=cp[id]*vp[id]/vs/taup;
      //     //      cscp2[crupos[id]]=1/taup*vp[id]/vs;
      // #else
      Idps[crupos[id]] = (cp[id] - cs[crupos[id]]) / taup;
      // #endif





      // #ifdef ___CPDIFF




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


      // #ifdef ___DEBUG
      //     if (isnan(Idps[crupos[id]])) //(Idsi != Idsi)
      //     {
      //       cout << setprecision(10) << id << "\t" << Idps[crupos[id]] << "\t cp=" << cp[id] << "\t cs=" << cs[crupos[id]] << endl;
      //       bSTOP = true;
      //     }
      // #endif

      // #ifdef ___NCX
      double dcp = get_cleft_caj_inst_buffering(cp[id])  * (Ir - Ica + jnaca - Idps[crupos[id]]);
      // #else
      //     double dcp = get_cleft_caj_inst_buffering(cp[id])  * (Ir - Ica - Idps[crupos[id]]);
      // #endif
      cp[id] += dcp * dt;
      // #else
      //     cp[id] = newcp;
      // #endif
      cjsr[id] += dcjsr * dt;

    }


    #pragma omp master
    {
      num_open_ICaL = total_ICaL_open;
    }


    #pragma omp for reduction(+: sumjup, sumjleak,sum_jnaca_sl,sumjcabk,sumjslcap) /*schedule (dynamic)*/
    // #pragma ivdep
    // #pragma vector always
    for (int id = 0; id < nn; id++)
    {
      //SERCA Uptake current Iup
      const double H = 1.787;
      double Iup = vup * (pow(ci[id] / kup, H) - pow(cnsr[id] / KNSR, H)) / (1 + pow(ci[id] / kup, H) + pow(cnsr[id] / KNSR, H));

      j_serca[id] = Iup; //  [uM/ms]s

      //Leak current Ileak
      const double KJSR = 500;
      double cjsr2 = cnsr[id] * cnsr[id];
      double Ileak = gleak * (cjsr2 * (cnsr[id] - ci[id])) / (cjsr2 + KJSR * KJSR);

      double j_CaT = tubule_flag[id] * ICaT_per_CRU * Cmem  / (2.0 * F) / (1e-15 * vs);
      double jnaca = tubule_flag[id] * (1 - NCXalpha) * NCX.compute_NCX(nai, nao, cs[id], cao);
      // double jnaca = tubule_flag[id] * (1 - NCXalpha) * NCX.compute_NCX_version_2(v, nai, nao, cs[id], cao);


      ncx_array[id] = jnaca;
      // -------Icabk (bac3*kground SL Ca flux) following Shannon --------------
      double eca = rtf2 * log(cao * 1000 / cs[id]);


      double icabk = gcabk * (v - eca); // current [A/F]  negative current
      //Icabk[A/F]*(Cm[pF]/(65*27*11*fine3))/(z*F[C/mol]*(vs*10-9/fine3)[ul])=Icabk*3.3286  -> *) xw: Cm ~ 310pF, larger than 110pF
      // double jcabk=icabk*3.3286;//coming in
      double jcabk = tubule_flag[id] * icabk * Cmem * (1e10); //  *) xw: change the Cmem to GB model
      icabk_array[id] = icabk; // *) xw: output [A/F]
      // -------Islcap (SL Ca pump) following Shannon --------------
      const double vmax = 2.2 * 0.01;
      const double km = 0.5;
      const double h = 1.6;

      // change tubule_flag to something else if finemesh >1 // 18:57:54, Wed, 12-December-2018, By Haibo
      double islcap = tubule_flag[id] * qslcap * vmax / (1 + pow(km / cs[id], h)); //
      // double jslcap=islcap*3.3286;//going out
      double jslcap = tubule_flag[id] * islcap *  Cmem * (1e10); //  *) xw: change the Cmem to GB model
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

      // #ifdef ___EGTA
      //     const double konEGTA = 4E-3;
      //     const double koffEGTA = 2E-3;
      //     double IEGTAi = konEGTA * ci[id] * (BEGTA - caEGTAi[id]) - koffEGTA * caEGTAi[id];
      //     double IEGTAs = konEGTA * cs[id] * (BEGTA - caEGTAs[id]) - koffEGTA * caEGTAs[id];
      // #endif



      // double ITCs=kon*cs[id]*(BT-cats[id])-koff*cats[id];
      double ITCs = 0;

      //Diffusion from submembrane to myoplasm Idsi
      double Idsi = (cs[id] - ci[id]) / tausi;
      // #ifdef ___DEBUG
      //     if (isnan(Idsi)) //(Idsi != Idsi)
      //     {
      //       cout << setprecision(10) << id << "\t" << Idsi << "\t cs=" << cs[id] << "\t ci=" << ci[id] << endl;
      //       bSTOP = true;
      //     }
      // #endif


      // #ifdef ___NO_CS_BUFFER
      // #ifdef ___NO_DIFFUSION
      //     double newcs = (cscp1[id] + jnaca + ci[id] / tausi - ITCs + csmn[id] - jcabk - jslcap - j_CaT) / (cscp2[id] + 1 / tausi);
      // #else
      //     double newcs = (cscp1[id] + jnaca + ci[id] / tausi - ITCs + csmn[id] - jcabk - jslcap - j_CaT) / (cscp2[id] + 1 / tausi + taumninv);
      // #endif
      //     if (newcs <= 0 )newcs = 0.00001; // to avoid div 0
      // #else



      // #ifdef ___EGTA
      //     double dcs = get_submem_casl_inst_buffering(cs[id]) * (Idps[id] * vp[id] / vs + jnaca - Idsi - ITCs - IEGTAs + Ics[id] - jcabk - jslcap - j_CaT);
      // #else
      //     double dcs = get_submem_casl_inst_buffering(cs[id]) * (Idps[id] * vp[id] / vs + jnaca - Idsi - ITCs + Ics[id] - jcabk - jslcap - j_CaT);
      // #endif
      // #endif


      // #ifdef ___EGTA
      //     double dci = get_cytosol_cai_inst_buffering(ci[id]) * (Idsi * (vs / vi) - Iup + Ileak - ITCi - IEGTAi + Ici[id]);
      // #else
      //     double dci = get_cytosol_cai_inst_buffering(ci[id]) * (Idsi * (vs / vi) - Iup + Ileak - ITCi + Ici[id]);
      // #endif


      double dcs = get_submem_casl_inst_buffering(cs[id]) * (Idps[id] * vp[id] / vs + jnaca - Idsi - ITCs + Ics[id] - jcabk - jslcap - j_CaT);
      double dci = get_cytosol_cai_inst_buffering(ci[id]) * (Idsi * (vs / vi) - Iup + Ileak - ITCi + Ici[id]);

      double dcnsr = ((Iup - Ileak) * (vi / vnsr) - Itr[id] * (vjsr / vnsr) + Icnsr[id]);



      ci[id] += dci * dt;
      cati[id] += ITCi * dt;
      // #ifdef ___NO_CS_BUFFER
      //     cs[id] = newcs;
      // #else
      cs[id] += dcs * dt;
      // #endif
      //      cats[id]+=ITCs*dt;
      cnsr[id] += dcnsr * dt;

      // #ifdef ___EGTA
      //     caEGTAi[id] += IEGTAi * dt;
      //     caEGTAs[id] += IEGTAs * dt;

      // #endif

      // *)xw: the sum of flux on surface, because they are 0 inside of the myocyte
      sumjup += Iup;
      sumjleak += Ileak;
      sum_jnaca_sl += jnaca;
      sumjcabk += jcabk;
      sumjslcap += jslcap;

    }
  }

  irave = sumir / n;
  iupave = sumjup / nn;

  ica_stan = 0.001 * sumica * (1e-15) * 2 * F * 1000 / Cmem; //  [pA/pF]

  incx_stan = sum_jnaca_sl * 0.001 * vs * (1e-15) * F * 1000 / Cmem   + sum_jnaca_j_flux * 0.001 * (1e-15) * F * 1000 / Cmem; //  [pA/pF]
  icabk_stan = sumjcabk * 0.001 * vs * (1e-15) * 2 * F * 1000 / Cmem; //  [pA/pF]
  ipca_stan = sumjslcap * 0.001 * vs * (1e-15) * 2 * F * 1000 / Cmem; //  [pA/pF]

}


