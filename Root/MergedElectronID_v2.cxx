#include "HGamGamStar/MergedElectronID_v2.h"
#include "PathResolver/PathResolver.h"

//______________________________________________________________________________
HG::MergedElectronID_v2::MergedElectronID_v2()
{

}

//______________________________________________________________________________
HG::MergedElectronID_v2::~MergedElectronID_v2()
{

}






//______________________________________________________________________________

EL::StatusCode HG::MergedElectronID_v2::initialize(Config &config)
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
bool HG::MergedElectronID_v2::passPIDCut(const xAOD::Electron &ele, bool isMC) const{


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

    unsigned iPt = 99;
    if      (pte < 30.) iPt = 0;
    else if (pte < 40.) iPt = 1;
    else if (pte < 50.) iPt = 2;
    else if (pte < 60.) iPt = 3;
    else iPt = 4;

    unsigned iEta = 99;
    if      (fabs(etae) < 0.80) iEta = 0;
    else if (fabs(etae) < 1.52) iEta = 1;
    else if (fabs(etae) < 2.01) iEta = 2;
    else iEta = 3;

    static const std::vector<std::vector<float>> cut_vtx_dphi ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.030000,  0.030000,  0.030000,  0.030000},
        { 0.020000,  0.020000,  0.020000,  0.020000},
        { 0.020000,  0.020000,  0.020000,  0.020000},
        { 0.020000,  0.020000,  0.020000,  0.020000},
        { 0.020000,  0.020000,  0.020000,  0.020000} });

    if (fabs(vtx_dphi) > cut_vtx_dphi[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_vtx_deta ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.010000,  0.010000,  0.010000,  0.010000},
        { 0.010000,  0.010000,  0.010000,  0.010000},
        { 0.010000,  0.010000,  0.010000,  0.010000},
        { 0.010000,  0.010000,  0.010000,  0.010000},
        { 0.010000,  0.010000,  0.010000,  0.010000} });

    if (fabs(vtx_deta) > cut_vtx_deta[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_Rhad ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.010594,  0.027027,  0.031010,  0.026704},
        { 0.018203,  0.010903,  0.011385,  0.008061},
        { 0.008006,  0.011777,  0.011486,  0.007800},
        { 0.008654,  0.009434,  0.008794,  0.009297},
        { 0.007920,  0.010288,  0.006657,  0.005030} });

    if (Rhad > cut_Rhad[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_Eratio ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.832966,  0.858457,  0.921214,  0.903645},
        { 0.750413,  0.820701,  0.890966,  0.871364},
        { 0.866114,  0.806071,  0.927473,  0.960342},
        { 0.895195,  0.827477,  0.806594,  0.514015},
        { 0.825656,  0.920766,  0.934981,  0.902925} });

    if (Eratio < cut_Eratio[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_Rphi ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.595356,  0.361620,  0.610367,  0.678678},
        { 0.378545,  0.395608,  0.587361,  0.634550},
        { 0.362006,  0.508373,  0.435598,  0.743595},
        { 0.609350,  0.516148,  0.772923,  0.847031},
        { 0.685707,  0.660021,  0.792450,  0.742276} });

    if (Rphi < cut_Rphi[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_Reta ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.936303,  0.918044,  0.925299,  0.919190},
        { 0.939417,  0.898833,  0.928949,  0.916855},
        { 0.931920,  0.936448,  0.925742,  0.920413},
        { 0.941888,  0.929478,  0.927457,  0.927671},
        { 0.945474,  0.935034,  0.929411,  0.928840} });

    if (Reta < cut_Reta[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_f3 ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.03,  0.03,  0.03,  0.03},
        { 0.03,  0.03,  0.03,  0.03},
        { 0.03,  0.03,  0.03,  0.03},
        { 0.03,  0.03,  0.03,  0.03},
        { 0.03,  0.03,  0.03,  0.03} });

    if (f3 > cut_f3[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_wEta2 ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.012756,  0.013344,  0.013959,  0.013965},
        { 0.013086,  0.012040,  0.012791,  0.012036},
        { 0.011536,  0.011727,  0.012369,  0.012347},
        { 0.016481,  0.011540,  0.011405,  0.012936},
        { 0.010545,  0.012633,  0.012000,  0.015539} });

    if (wEta2 > cut_wEta2[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_wTotS1 ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 4.696972,  3.594180,  3.378504,  1.467929},
        { 2.747719,  2.792041,  5.512513,  1.487706},
        { 2.842709,  3.917183,  3.661756,  1.411089},
        { 2.500456,  4.121699,  4.071177,  2.191503},
        { 4.803770,  3.196767,  3.288494,  1.783990} });

    if (wTotS1 > cut_wTotS1[iPt][iEta]) return false;
    if (wTotS1 < -100.) return false;

    static const std::vector<std::vector<float>> cut_trk_TRT_PID1 ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { -0.075664,  -0.081856,   0.034933,  -100.000000},
        { -0.071126,  -0.099147,  -0.088572,  -100.000000},
        { -0.068432,  -0.081285,  -0.011563,  -100.000000},
        { -0.291272,  -0.002691,  -0.027964,  -100.000000},
        { -0.048709,  -0.154170,  -0.215961,  -100.000000} });

    if (trk_TRT_PID1 < cut_trk_TRT_PID1[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_trk_TRT_PID2 ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { -0.054928,  -0.046679,  -0.315406,  -100.000000},
        { -0.118842,  -0.049472,  -0.074673,  -100.000000},
        { -0.154553,  -0.159510,  -0.417202,  -100.000000},
        { -0.087112,  -0.350279,  -0.192150,  -100.000000},
        { -0.271082,  -0.132265,  -0.195975,  -100.000000} });

    if (trk_TRT_PID2 < cut_trk_TRT_PID2[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_EoP ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        {  2.614554,  1.798164,  2.499431,  1.889961},
        { 17.921568, 14.419681,  2.245576,  2.079409},
        {  1.455041,  5.202778,  5.788745,  2.369355},
        {  1.724084,  3.581217,  2.530357,  2.105054},
        {  1.996022,  6.719277,  1.964092,  6.231330} });

    if (EoP > cut_EoP[iPt][iEta]) return false;

    static const std::vector<std::vector<float>> cut_trk_dEta1 ({
        // 0-0.8,    to1.37,    to2.01,    to2.37
        { 0.001922,  0.002309,  0.002653,  0.001359},
        { 0.001160,  0.001497,  0.002411,  0.001821},
        { 0.004640,  0.000985,  0.001250,  0.001652},
        { 0.001831,  0.001494,  0.002562,  0.003560},
        { 0.003580,  0.001038,  0.003355,  0.001250} });

    if (trk_dEta1 > cut_trk_dEta1[iPt][iEta]) return false;

    return true;
}

float HG::MergedElectronID_v2::correctDeta1(const xAOD::Electron &ele, bool isMC) const{

  double corr =  0;
  if(!isMC){
    double eta = ele.caloCluster()->etaBE(1);
    double phi = ele.caloCluster()->phiBE(2);
    int    bin = m_sdetaCorr->FindBin( eta, phi );
    corr = m_sdetaCorr->GetBinContent(bin) * 1e-3;
  }
  return ele.trackCaloMatchValue(xAOD::EgammaParameters::TrackCaloMatchType::deltaEta1) - corr;

}

float HG::MergedElectronID_v2::GetScaleFactor(const xAOD::Electron &ele, MergedSystematic sys) const{

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
      {      1.011,     0.983,   999,     0.945,     1.121}, // pT [20, 30] GeV
      {      1.002,     0.941,   999,     0.933,     0.876}, // pT [30, 40] GeV
      {      1.018,     0.966,   999,     0.927,     0.963}, // pT [40, 50] GeV
      {      1.159,     0.959,   999,     0.954,     0.993}  // pT [50, 100+] GeV
    });

  static const std::vector<std::vector<float>> stat_unc ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {   0.021,  0.020,   999,  0.025,  0.030}, // pT [20, 30] GeV
      {   0.023,  0.027,   999,  0.033,  0.045}, // pT [30, 40] GeV
      {   0.042,  0.047,   999,  0.066,  0.110}, // pT [40, 50] GeV
      {   0.036,  0.040,   999,  0.068,  0.093}  // pT [50, 100+] GeV
    });

  static const std::vector<std::vector<float>> syst_unc ({
      // 0-0.8,to1.37,to1.52,to2.01,to2.37
      {   0.002,  0.009,   999,  0.033,  0.009}, // pT [20, 30] GeV
      {   0.023,  0.017,   999,  0.015,  0.006}, // pT [30, 40] GeV
      {   0.030,  0.037,   999,  0.060,  0.029}, // pT [40, 50] GeV
      {   0.080,  0.029,   999,  0.039,  0.034}  // pT [50, 100+] GeV
    });

  float sf = nominal_sf[iPt][iEta];
  if (sys == MERGEDUNC_STAT_UP  ) sf = sf + stat_unc[iPt][iEta];
  if (sys == MERGEDUNC_STAT_DOWN) sf = sf - stat_unc[iPt][iEta];
  // sys uncertainties are relative
  if (sys == MERGEDUNC_SYST_UP  ) sf = sf * (1. + syst_unc[iPt][iEta]);
  if (sys == MERGEDUNC_SYST_DOWN) sf = sf * (1. - syst_unc[iPt][iEta]);

  return sf;
}
