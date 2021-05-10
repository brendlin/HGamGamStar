#include "HGamGamStar/MergedElectronID_v2F.h"
#include "PathResolver/PathResolver.h"

//______________________________________________________________________________
HG::MergedElectronID_v2F::MergedElectronID_v2F() :
  m_detaCorrectionTool()
{

}

//______________________________________________________________________________
HG::MergedElectronID_v2F::~MergedElectronID_v2F()
{

}






//______________________________________________________________________________

EL::StatusCode HG::MergedElectronID_v2F::initialize(Config &config)
{
  // These are enforced via the original MergedElectronID passPreselection cut!
  // m_PreselNPassBlayer = config.getInt("MergedElectrons.Preselection.NtracksPassingBlayer",1);
  // m_PreselRhad        = config.getNum("MergedElectrons.Preselection.RhadMin",0.10);
  // m_mergedElePtCut    = config.getNum("MergedElectrons.Selection.PtPreCutGeV",20.) * GeV;
  // m_mergedEleEtaCut   = config.getNum("MergedElectrons.Selection.MaxAbsEta"  ,2.37);
  m_deltaEtaHistName  = config.getStr("MergedElectrons.Selection.deltaEtaHistName","sdetaRun2") ;
  m_deltaEtaFileName  = config.getStr("MergedElectrons.Selection.deltaEtaFileName","HGamGamStar/DeltaEtaRun2.root") ;

  auto detaFile= TFile::Open( PathResolverFindCalibFile( m_deltaEtaFileName ).data() );
  m_sdetaCorr =  (TH2*) detaFile->Get(m_deltaEtaHistName.data());
  m_sdetaCorr->SetDirectory(0);
  detaFile->Close();
  delete detaFile;

  m_detaCorrectionTool.setTypeAndName("ElectronPhotonVariableCorrectionTool/kurtDetaTool");
  std::string configFile = "HGamGamStar/DetaMasterCorrectionTool.conf";
  m_detaCorrectionTool.setProperty("ConfigFile", configFile);
  if (m_detaCorrectionTool.initialize().isFailure())
    fatal("Could not initialize Deta correction tool.");

  return EL::StatusCode::SUCCESS;
}


int HG::MergedElectronID_v2F::getEtaBin( float eta )const{
  eta =  std::fabs(eta);
  int cutBinEta = -1;
  if( eta <= 0.8 ){
    cutBinEta = 0;
  } else if( eta <= 1.37 ){
    cutBinEta = 1;
  } else if( eta <= 1.52){
    cutBinEta = -1;
  } else if( eta <= 2.01){
    cutBinEta = 2;
  } else if( eta <= 2.37){
      cutBinEta = 3;
  } else {
    cutBinEta = -1;
  }
  return cutBinEta;
}

int HG::MergedElectronID_v2F::getPtBin( float pt )const {
  int cutBinPT;
  if( pt <= 30)
    cutBinPT = 0;
  else if( pt <= 40)
    cutBinPT = 1;
  else if( pt <= 50)
    cutBinPT = 2;
  else if( pt <= 60)
    cutBinPT = 3;
  else
    cutBinPT = 4;
  return cutBinPT;
}


//______________________________________________________________________________
bool HG::MergedElectronID_v2F::passPIDCut(const xAOD::Electron &ele, bool isMC) const{


  int trk1_index =  HG::EleAcc::vtxTrkIndex1(ele);
  int trk2_index =  HG::EleAcc::vtxTrkIndex2(ele);

  // Two tracks are found
  if(trk1_index < 0 || trk2_index < 0)
    return false;

  // Old MxAOD ( < p3877 ) that does not save the necessary information
  if ( !HG::EleAcc::vtxE.isAvailable(ele) )
    return false;

  // Vertex Fit failed
  if( HG::EleAcc::vtxE(ele) < 0)
    return false;

  // calculate shower shapes and other discriminating variables
  float trk_dEta1 =  0; // correctDeta1( ele, isMC);

//    double f1 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::f1);
//     double fSide = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::fracs1);

  auto ele_tp  =  ele.trackParticle(trk1_index);
  auto ele_tp2 =  ele.trackParticle(trk2_index);

  uint8_t shared;
  int sumOfSplitHits = 0;
  if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerSharedHits) )
    sumOfSplitHits += (int)shared;
  if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerSharedHits) )
    sumOfSplitHits += (int)shared;
  if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerSplitHits ) )
    sumOfSplitHits += (int)shared;
  if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerSplitHits) )
    sumOfSplitHits += (int)shared;

  if(sumOfSplitHits == 0){
    if( fabs( ele_tp2->z0() - ele_tp->z0() ) >  1.0 )
      return false;
  }


  float vtx_deta  = HG::EleAcc::vtxdEta(ele);
  float vtx_dphi  = HG::EleAcc::vtxdPhi(ele);
  
  float trk1_eta = ele_tp->eta();
  float trk2_eta = ele_tp2->eta();
  float trk_deta = fabs(trk1_eta-trk2_eta);

  float pte  = ele.pt()/1000.;
  float etae = ele.caloCluster()->etaBE(2);

  float Rhad = HG::EleAcc::RhadForPID(ele);
  float trk_TRT_PID1 = HG::EleAcc::vtxTrk1_TRT_PID_trans(ele);
  float trk_TRT_PID2 = HG::EleAcc::vtxTrk2_TRT_PID_trans(ele);
  float EoP =  ele.e() / HG::EleAcc::vtxE(ele); // The ele energy is calibrated earlier.


  float Eratio = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Eratio);
  //recovering efficiency at larger deta separation
  if(trk_deta>0.0031) Eratio =1.;
  
  float wTotS1 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::wtots1);

  float Reta = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Reta);
  float Rphi = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rphi);
  float wEta2 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::weta2);
  float f3 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::f3);


  std::vector<float> values{ vtx_dphi,  Rhad, Eratio, Rphi,  Reta, wEta2, wTotS1, trk_TRT_PID1, trk_TRT_PID2, EoP, trk_dEta1 };
  
  //enum ShowerShapeCuts { dphiVtx=0,  rHad,  eratio, rphi, reta, weta2, wstot, TRT1, TRT2, EoverPvtx, deltaEta };
  static const std::vector<float> lessOrGreater { -1 , -1,  1,   1,    1,    -1,    -1,    1,    1,    -1,        -1 };
                            //  0        1          2         3          4         5          6          7           8           9          10
  static const std::vector<std::vector<std::vector<float>>> cuts 
//Original values stored for ease of validation
/*{{{ 0.0300000,   0.0105940,   0.8329660,   0.5953560,   0.9363030,   0.0127560,   4.6969719,  -0.0756640,  -0.0549280,   2.6145539,   0.0019220,  },
{ 0.0300000,   0.0270270,   0.8584570,   0.3616200,   0.9180440,   0.0133440,   3.5941801,  -0.0818560,  -0.0466790,   1.7981640,   0.0023090,  },
{ 0.0300000,   0.0310100,   0.9212140,   0.6103670,   0.9252990,   0.0139590,   3.3785040,   0.0349330,  -0.3154060,   2.4994309,   0.0026530,  },
{ 0.0300000,   0.0267040,   0.9036450,   0.6786780,   0.9191900,   0.0139650,   1.4679290,  -100.0000000,  -100.0000000,   1.8899610,   0.0013590,  },
},
{{ 0.0200000,   0.0182030,   0.7504130,   0.3785450,   0.9394170,   0.0130860,   2.7477190,  -0.0711260,  -0.1188420,  17.9215679,   0.0011600,  },
{ 0.0200000,   0.0109030,   0.8207010,   0.3956080,   0.8988330,   0.0120400,   2.7920411,  -0.0991470,  -0.0494720,  14.4196806,   0.0014970,  },
{ 0.0200000,   0.0113850,   0.8909660,   0.5873610,   0.9289490,   0.0127910,   5.5125132,  -0.0885720,  -0.0746730,   2.2455759,   0.0024110,  },
{ 0.0200000,   0.0080610,   0.8713640,   0.6345500,   0.9168550,   0.0120360,   1.4877059,  -100.0000000,  -100.0000000,   2.0794089,   0.0018210,  },
},
{{ 0.0200000,   0.0080060,   0.8661140,   0.3620060,   0.9319200,   0.0115360,   2.8427091,  -0.0684320,  -0.1545530,   1.4550411,   0.0046400,  },
{ 0.0200000,   0.0117770,   0.8060710,   0.5083730,   0.9364480,   0.0117270,   3.9171829,  -0.0812850,  -0.1595100,   5.2027779,   0.0009850,  },
{ 0.0200000,   0.0114860,   0.9274730,   0.4355980,   0.9257420,   0.0123690,   3.6617560,  -0.0115630,  -0.4172020,   5.7887449,   0.0012500,  },
{ 0.0200000,   0.0078000,   0.9603420,   0.7435950,   0.9204130,   0.0123470,   1.4110889,  -100.0000000,  -100.0000000,   2.3693550,   0.0016520,  },
},
{{ 0.0200000,   0.0086540,   0.8951950,   0.6093500,   0.9418880,   0.0164810,   2.5004561,  -0.2912720,  -0.0871120,   1.7240840,   0.0018310,  },
{ 0.0200000,   0.0094340,   0.8274770,   0.5161480,   0.9294780,   0.0115400,   4.1216989,  -0.0026910,  -0.3502790,   3.5812171,   0.0014940,  },
{ 0.0200000,   0.0087940,   0.8065940,   0.7729230,   0.9274570,   0.0114050,   4.0711770,  -0.0279640,  -0.1921500,   2.5303569,   0.0025620,  },
{ 0.0200000,   0.0092970,   0.5140150,   0.8470310,   0.9276710,   0.0129360,   2.1915030,  -100.0000000,  -100.0000000,   2.1050539,   0.0035600,  },
},
{{ 0.0200000,   0.0079200,   0.8256560,   0.6857070,   0.9454740,   0.0105450,   4.8037701,  -0.0487090,  -0.2710820,   1.9960220,   0.0035800,  },
{ 0.0200000,   0.0102880,   0.9207660,   0.6600210,   0.9350340,   0.0126330,   3.1967671,  -0.1541700,  -0.1322650,   6.7192769,   0.0010380,  },
{ 0.0200000,   0.0066570,   0.9349810,   0.7924500,   0.9294110,   0.0120000,   3.2884941,  -0.2159610,  -0.1959750,   1.9640920,   0.0033550,  },
{ 0.0200000,   0.0050300,   0.9029250,   0.7422760,   0.9288400,   0.0155390,   1.7839900,  -100.0000000,  -100.0000000,   6.2313299,   0.0012500,  },
},
};*/
{{{ 0.0300000,   0.0092465,   0.8322660,   0.6051060,   0.9325530,   0.0128523,   4.7652221,  -0.0756640,  -0.0549280,   2.6145539,   0.0019220,  },
{ 0.0300000,   0.0258445,   0.8568570,   0.3651700,   0.9148189,   0.0133912,   3.6860552,  -0.0818560,  -0.0466790,   1.7981640,   0.0023090,  },
{ 0.0300000,   0.0305975,   0.9200140,   0.6137670,   0.9185490,   0.0142838,   3.5008166,   0.0349330,  -0.3154060,   2.4994309,   0.0026530,  },
{ 0.0300000,   0.0274190,   0.9034450,   0.6818280,   0.9113150,   0.0143045,   1.5946790,  -100.0000000,  -100.0000000,   1.8899610,   0.0013590,  },
},
{{ 0.0200000,   0.0174330,   0.7502130,   0.3884450,   0.9395670,   0.0131665,   2.8110940,  -0.0711260,  -0.1188420,  17.9215679,   0.0011600,  },
{ 0.0200000,   0.0105730,   0.8183010,   0.4004080,   0.8952330,   0.0120680,   2.8827410,  -0.0991470,  -0.0494720,  14.4196806,   0.0014970,  },
{ 0.0200000,   0.0117150,   0.8887660,   0.5917610,   0.9241490,   0.0131360,   5.6227632,  -0.0885720,  -0.0746730,   2.2455759,   0.0024110,  },
{ 0.0200000,   0.0083360,   0.8699640,   0.6381500,   0.9111550,   0.0123760,   1.6258309,  -100.0000000,  -100.0000000,   2.0794089,   0.0018210,  },
},
{{ 0.0200000,   0.0064110,   0.8669140,   0.3710060,   0.9290700,   0.0116760,   2.9077091,  -0.0684320,  -0.1545530,   1.4550411,   0.0046400,  },
{ 0.0200000,   0.0116670,   0.8086710,   0.5115730,   0.9355480,   0.0117095,   3.9984329,  -0.0812850,  -0.1595100,   5.2027779,   0.0009850,  },
{ 0.0200000,   0.0107260,   0.9256730,   0.4386980,   0.9216920,   0.0126405,   3.8042560,  -0.0115630,  -0.4172020,   5.7887449,   0.0012500,  },
{ 0.0200000,   0.0076900,   0.9581420,   0.7441950,   0.9147630,   0.0127070,   1.5492139,  -100.0000000,  -100.0000000,   2.3693550,   0.0016520,  },
},
{{ 0.0200000,   0.0077190,   0.8955950,   0.6147500,   0.9414380,   0.0165300,   2.5638311,  -0.2912720,  -0.0871120,   1.7240840,   0.0018310,  },
{ 0.0200000,   0.0088290,   0.8304770,   0.5200480,   0.9284280,   0.0115890,   4.2143240,  -0.0026910,  -0.3502790,   3.5812171,   0.0014940,  },
{ 0.0200000,   0.0082440,   0.8061940,   0.7783230,   0.9244570,   0.0116975,   4.1513019,  -0.0279640,  -0.1921500,   2.5303569,   0.0025620,  },
{ 0.0200000,   0.0098470,   0.5132150,   0.8497310,   0.9202210,   0.0132755,   2.3222530,  -100.0000000,  -100.0000000,   2.1050539,   0.0035600,  },
},
{{ 0.0200000,   0.0070217,   0.8261893,   0.6911070,   0.9446240,   0.0106500,   4.8712282,  -0.0487090,  -0.2710820,   1.9960220,   0.0035800,  },
{ 0.0200000,   0.0096830,   0.9226993,   0.6643210,   0.9335840,   0.0126773,   3.2731421,  -0.1541700,  -0.1322650,   6.7192769,   0.0010380,  },
{ 0.0200000,   0.0065470,   0.9341810,   0.7976500,   0.9258110,   0.0122888,   3.3715358,  -0.2159610,  -0.1959750,   1.9640920,   0.0033550,  },
{ 0.0200000,   0.0055067,   0.9025250,   0.7455760,   0.9213900,   0.0158988,   1.9091567,  -100.0000000,  -100.0000000,   6.2313299,   0.0012500,  },
},
};




  if(f3 > 0.03)
    return false;
  if( fabs(vtx_deta)> 0.010000)
    return false; 
  if(wTotS1<=-100)
    return false;

  int cutBinEta = getEtaBin( etae );
  if( cutBinEta < 0 )
    return false;
  int cutBinPT = getPtBin(pte);
  if( cutBinPT < 0 )
    return false;
  for( size_t i(0); i < values.size(); ++i){
    if( lessOrGreater[i] * values[i] < lessOrGreater[i] * cuts[cutBinPT][cutBinEta][i] ){
      return false;
    }
  }
  return true;
}

float HG::MergedElectronID_v2F::correctDeta1(const xAOD::Electron &ele, bool isMC) const{

  double corr =  0;
  if(!isMC){
    double eta = ele.caloCluster()->etaBE(2);
    double phi = ele.caloCluster()->phiBE(2);
    int    bin = m_sdetaCorr->FindBin( eta, phi );
    corr = m_sdetaCorr->GetBinContent(bin) * 1e-3;
  }
  return ele.trackCaloMatchValue(xAOD::EgammaParameters::TrackCaloMatchType::deltaEta1) - corr;

}

void HG::MergedElectronID_v2F::nilsDeta1(xAOD::Electron &ele, bool isMC) const{

  double corr =  0;
  if(!isMC){
    m_detaCorrectionTool->applyCorrection(ele);
  }
  return;
}

float HG::MergedElectronID_v2F::GetScaleFactor(const xAOD::Electron &ele, MergedSystematic sys) const{

  if (!HG::isMC()) return 1;

  float pt  = ele.pt()/1000.;
  float eta = fabs(ele.caloCluster()->etaBE(2));

  // This should really not happen:
  if (pt < 20) return 1;
  if (eta > 2.37) return 1;

  unsigned iPt = 0;
  if      (pt < 30.0) iPt = 0;
  else if (pt < 40.0) iPt = 1;
  else if (pt < 50.0) iPt = 2;
  else iPt = 3;

  unsigned iEta = 0;
  if      (eta < 0.8 ) iEta = 0;
  else if (eta < 1.37) iEta = 1;
  else if (eta < 1.52) iEta = 2; // crack will be excluded anyway, but we leave it in the matrix.
  else if (eta < 2.01) iEta = 3;
  else iEta = 4;

  static const std::vector<std::vector<float>> nominal_sf ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {      0.998,     0.976,   999,     0.940,     1.054}, // pT [20, 30] GeV
      {      0.999,     0.922,   999,     0.926,     0.925}, // pT [30, 40] GeV
      {      1.030,     1.003,   999,     0.944,     0.897}, // pT [40, 50] GeV
      {      1.077,     0.952,   999,     0.974,     0.832}  // pT [50, 100+] GeV
    });

  static const std::vector<std::vector<float>> stat_unc ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {   0.020,  0.020,   999,  0.023,  0.022}, // pT [20, 30] GeV
      {   0.023,  0.025,   999,  0.030,  0.035}, // pT [30, 40] GeV
      {   0.041,  0.045,   999,  0.060,  0.066}, // pT [40, 50] GeV
      {   0.034,  0.039,   999,  0.065,  0.065}  // pT [50, 100+] GeV
    });

  static const std::vector<std::vector<float>> syst_unc ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {   0.001,  0.009,   999,  0.034,  0.008}, // pT [20, 30] GeV
      {   0.027,  0.018,   999,  0.014,  0.010}, // pT [30, 40] GeV
      {   0.003,  0.039,   999,  0.065,  0.031}, // pT [40, 50] GeV
      {   0.083,  0.022,   999,  0.035,  0.033}  // pT [50, 100+] GeV
    });

  float sf = nominal_sf[iPt][iEta];
  if (sys == MERGEDUNC_STAT_UP  ) sf = sf + stat_unc[iPt][iEta];
  if (sys == MERGEDUNC_STAT_DOWN) sf = sf - stat_unc[iPt][iEta];
  // sys uncertainties are relative
  if (sys == MERGEDUNC_SYST_UP  ) sf = sf * (1. + syst_unc[iPt][iEta]);
  if (sys == MERGEDUNC_SYST_DOWN) sf = sf * (1. - syst_unc[iPt][iEta]);

  return sf;
}
