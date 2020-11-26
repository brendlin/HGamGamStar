#include "HGamGamStar/MergedElectronID_v3.h"
#include "PathResolver/PathResolver.h"

//______________________________________________________________________________
HG::MergedElectronID_v3::MergedElectronID_v3()
{

}

//______________________________________________________________________________
HG::MergedElectronID_v3::~MergedElectronID_v3()
{

}






//______________________________________________________________________________

EL::StatusCode HG::MergedElectronID_v3::initialize(Config &config)
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

  return EL::StatusCode::SUCCESS;
}


//______________________________________________________________________________
bool HG::MergedElectronID_v3::passPIDCut(const xAOD::Electron &ele, bool isMC) const{


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
    float trk_dEta1 =  correctDeta1( ele, isMC);

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
    
    float trk_dPhiIP = EleAcc::deltaPhiTrksIP(ele);

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


    // 60-800 GeV - at 60% WP

    if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.006560 && Eratio>0.887153 && Rphi>0.354076 && Reta>0.946109 && f3<0.03 && wEta2<0.011205 && wTotS1<3.158942 && trk_dPhiIP<0.150507 && trk_TRT_PID1>-0.211624 && trk_TRT_PID2>-0.262335 && EoP<1.414835 && trk_dEta1<0.039386)) return true;


    if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.007336 && Eratio>0.860833 && Rphi>0.622187 && Reta>0.926088 && f3<0.03 && wEta2<0.012049 && wTotS1<5.379566 && trk_dPhiIP<0.009770 && trk_TRT_PID1>-0.088449 && trk_TRT_PID2>-0.334621 && EoP<1.725958 && trk_dEta1<0.004401)) return true;


    if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.008747 && Eratio>0.901390 && Rphi>0.731188 && Reta>0.931949 && f3<0.03 && wEta2<0.011243 && wTotS1<5.905310 && trk_dPhiIP<0.013196 && trk_TRT_PID1>-0.628788 && trk_TRT_PID2>-0.082618 && EoP<1.981424 && trk_dEta1<0.004889)) return true;


    if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.007008 && Eratio>0.879733 && Rphi>0.789521 && Reta>0.927420 && f3<0.03 && wEta2<0.014103 && wTotS1<1.477728 && trk_dPhiIP<0.012838 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<1.904189 && trk_dEta1<0.004660)) return true;


    // 50-60 GeV - at 55% WP


    if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.010435 && Eratio>0.831915 && Rphi>0.576838 && Reta>0.947720 && f3<0.03 && wEta2<0.016926 && wTotS1<6.672210 && trk_dPhiIP<0.131297 && trk_TRT_PID1>-0.074202 && trk_TRT_PID2>-0.100878 && EoP<1.847073 && trk_dEta1<0.004698)) return true;


    if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.170988 && Eratio>0.315813 && Rphi>0.634531 && Reta>0.939027 && f3<0.03 && wEta2<0.012735 && wTotS1<3.461166 && trk_dPhiIP<0.014465 && trk_TRT_PID1>-0.351490 && trk_TRT_PID2>-0.120732 && EoP<1.621894 && trk_dEta1<0.001774)) return true;


    if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.014455 && Eratio>0.849625 && Rphi>0.714681 && Reta>0.928613 && f3<0.03 && wEta2<0.014186 && wTotS1<2.872270 && trk_dPhiIP<0.009929 && trk_TRT_PID1>-0.118859 && trk_TRT_PID2>-0.305445 && EoP<2.193178 && trk_dEta1<0.001412)) return true;


    if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.008910 && Eratio>0.839963 && Rphi>0.883823 && Reta>0.912895 && f3<0.03 && wEta2<0.014994 && wTotS1<2.186882 && trk_dPhiIP<0.010639 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<1.916158 && trk_dEta1<0.006709)) return true;


    // 40-50 GeV


    if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.007482 && Eratio>0.539577 && Rphi>0.344396 && Reta>0.937325 && f3<0.03 && wEta2<0.014920 && wTotS1<8.531023 && trk_dPhiIP<0.031725 && trk_TRT_PID1>-0.197248 && trk_TRT_PID2>-0.101549 && EoP<1.482727 && trk_dEta1<0.002771)) return true;


    if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.022826 && Eratio>0.902553 && Rphi>0.289044 && Reta>0.929602 && f3<0.03 && wEta2<0.012781 && wTotS1<3.812475 && trk_dPhiIP<0.047040 && trk_TRT_PID1>-0.059714 && trk_TRT_PID2>-0.162280 && EoP<1.690371 && trk_dEta1<0.001880)) return true;


    if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.010122 && Eratio>0.777696 && Rphi>0.689765 && Reta>0.931818 && f3<0.03 && wEta2<0.018233 && wTotS1<12.825915 && trk_dPhiIP<0.018945 && trk_TRT_PID1>-0.072239 && trk_TRT_PID2>-0.263620 && EoP<1.749146 && trk_dEta1<0.002461)) return true;


    if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.009852 && Eratio>0.923376 && Rphi>0.790327 && Reta>0.904579 && f3<0.03 && wEta2<0.015232 && wTotS1<1.949186 && trk_dPhiIP<0.004136 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<1.680171 && trk_dEta1<0.003559)) return true;


    // 30-40 GeV
    // changed to 45% WP

    if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.020486 && Eratio>0.832314 && Rphi>0.501272 && Reta>0.934665 && f3<0.03 && wEta2<0.015396 && wTotS1<3.728755 && trk_dPhiIP<0.010906 && trk_TRT_PID1>-0.066411 && trk_TRT_PID2>-0.141264 && EoP<6.930804 && trk_dEta1<0.001190)) return true;


    if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.954574 && Eratio>0.817603 && Rphi>0.325208 && Reta>0.933737 && f3<0.03 && wEta2<0.013057 && wTotS1<3.794917 && trk_dPhiIP<0.011270 && trk_TRT_PID1>-0.038766 && trk_TRT_PID2>-0.205847 && EoP<1.602289 && trk_dEta1<0.001598)) return true;


    if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.022812 && Eratio>0.713753 && Rphi>0.471539 && Reta>0.937620 && f3<0.03 && wEta2<0.013474 && wTotS1<3.824326 && trk_dPhiIP<0.019852 && trk_TRT_PID1>-0.123303 && trk_TRT_PID2>-0.096169 && EoP<2.765352 && trk_dEta1<0.001578)) return true;


    if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.005352 && Eratio>0.928801 && Rphi>0.567649 && Reta>0.912931 && f3<0.03 && wEta2<0.013158 && wTotS1<1.761505 && trk_dPhiIP<0.065222 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<1.580360 && trk_dEta1<0.002968)) return true;


    // 20-30 GeV
    // changed to 40% WP - smoothing the pt spectrum

    if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.032511 && Eratio>0.908455 && Rphi>0.342760 && Reta>0.909404 && f3<0.03 && wEta2<0.011093 && wTotS1<3.546919 && trk_dPhiIP<0.013798 && trk_TRT_PID1>-0.219264 && trk_TRT_PID2>-0.066326 && EoP<3.183485 && trk_dEta1<0.001190)) return true;


    if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.019496 && Eratio>0.599547 && Rphi>0.665738 && Reta>0.896621 && f3<0.03 && wEta2<0.014106 && wTotS1<4.649633 && trk_dPhiIP<0.008322 && trk_TRT_PID1>0.040271 && trk_TRT_PID2>-0.343206 && EoP<12.195064 && trk_dEta1<0.001898)) return true;


    if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.086947 && Eratio>0.821553 && Rphi>0.611141 && Reta>0.895740 && f3<0.03 && wEta2<0.013302 && wTotS1<35.360333 && trk_dPhiIP<0.004841 && trk_TRT_PID1>-0.074095 && trk_TRT_PID2>-0.099599 && EoP<2.413265 && trk_dEta1<0.001733)) return true;


    if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000 && trk_dPhiIP>-0.02)&&(Rhad<0.010985 && Eratio>0.462911 && Rphi>0.252542 && Reta>0.911919 && f3<0.03 && wEta2<0.013844 && wTotS1<2.959467 && trk_dPhiIP<0.008167 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<1.547559 && trk_dEta1<0.001967)) return true;

    
    return(false);
}

float HG::MergedElectronID_v3::correctDeta1(const xAOD::Electron &ele, bool isMC) const{

  double corr =  0;
  if(!isMC){
    double eta = ele.caloCluster()->etaBE(1);
    double phi = ele.caloCluster()->phiBE(2);
    int    bin = m_sdetaCorr->FindBin( eta, phi );
    corr = m_sdetaCorr->GetBinContent(bin) * 1e-3;
  }
  return ele.trackCaloMatchValue(xAOD::EgammaParameters::TrackCaloMatchType::deltaEta1) - corr;

}

float HG::MergedElectronID_v3::GetScaleFactor(const xAOD::Electron &ele, MergedSystematic sys) const{

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
  
  //SFs for v3
  static const std::vector<std::vector<float>> nominal_sf ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {      1.075,     0.882,   999,     1.063,     1.061}, // pT [20, 30] GeV
      {      1.056,     1.035,   999,     1.035,     1.006}, // pT [30, 40] GeV
      {      0.985,     1.134,   999,     0.846,     1.237}, // pT [40, 50] GeV
      {      1.115,     1.033,   999,     1.013,     0.848}  // pT [50, 100+] GeV
    });

  static const std::vector<std::vector<float>> stat_unc ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {   0.029,  0.032,   999,  0.036,  0.037}, // pT [20, 30] GeV
      {   0.031,  0.036,   999,  0.045,  0.046}, // pT [30, 40] GeV
      {   0.042,  0.055,   999,  0.073,  0.135}, // pT [40, 50] GeV
      {   0.036,  0.048,   999,  0.078,  0.089}  // pT [50, 100+] GeV
    });

  static const std::vector<std::vector<float>> syst_unc ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {   0.001,  0.008,   999,  0.034,  0.008}, // pT [20, 30] GeV
      {   0.027,  0.019,   999,  0.014,  0.010}, // pT [30, 40] GeV
      {   0.010,  0.036,   999,  0.064,  0.031}, // pT [40, 50] GeV
      {   0.083,  0.022,   999,  0.039,  0.035}  // pT [50, 100+] GeV
    });

  float sf = nominal_sf[iPt][iEta];
  if (sys == MERGEDUNC_STAT_UP  ) sf = sf + stat_unc[iPt][iEta];
  if (sys == MERGEDUNC_STAT_DOWN) sf = sf - stat_unc[iPt][iEta];
  // sys uncertainties are relative
  if (sys == MERGEDUNC_SYST_UP  ) sf = sf * (1. + syst_unc[iPt][iEta]);
  if (sys == MERGEDUNC_SYST_DOWN) sf = sf * (1. - syst_unc[iPt][iEta]);

  return sf;
}
