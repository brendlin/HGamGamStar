#include "HGamGamStar/ZyUnfoldingInput.h"
#include "HGamAnalysisFramework/HGamCommon.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>

// this is needed to distribute the algorithm to the workers
ClassImp(ZyUnfoldingInput)

ZyUnfoldingInput::ZyUnfoldingInput(const char *name)
  : HgammaAnalysis(name) {}

// Here you delete any memory you allocated during your analysis.
ZyUnfoldingInput::~ZyUnfoldingInput() { }

EL::StatusCode ZyUnfoldingInput::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  return EL::StatusCode::SUCCESS;
}

TH1F *ZyUnfoldingInput::createAndRegisterTH1F(TString name, int Nbins, double min, double max, TString title)
{
  histoStore()->createTH1F(name, Nbins, min, max, title);
  wk()->addOutput(histoStore()->getTH1F(name));
  return histoStore()->getTH1F(name);
}

TH1F *ZyUnfoldingInput::createAndRegisterTH1F(TString name, const std::vector<double>& xbins, TString title)
{
  histoStore()->createTH1F(name, xbins, title);
  wk()->addOutput(histoStore()->getTH1F(name));
  return histoStore()->getTH1F(name);
}

TH2F *ZyUnfoldingInput::createAndRegisterTH2F(TString name, const std::vector<double>& xbins, const std::vector<double>& ybins, TString title)
{
  histoStore()->createTH2F(name, xbins, ybins, title);
  wk()->addOutput(histoStore()->getTH2F(name));
  return histoStore()->getTH2F(name);
}

EL::StatusCode ZyUnfoldingInput::execute()
{
  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  if (!HG::isMAOD()) { HG::fatal("This code only works on MxAOD input."); }

  int id = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
  //hardcoding PDF weights since MetaData missing in MxAOD
  int meta_dsid =366146; //take MetaData info from DSID366146 as default (Sherpa 2.2.4 sample)
  if (id>=345775 && id<=345782) meta_dsid=345775; //or Madgraph
  #include "data/meta_366146.h"
  #include "data/meta_345775.h"
  if (HG::isMC()){
    if (meta_dsid==366146){
      for (int pdf=261001;pdf<=261100;++pdf) {CreateHists(Form("PDF%i",pdf));}
      CreateHists("PDF269000"); //NNPDF30_nnlo_as_0117
      CreateHists("PDF270000"); //NNPDF30_nnlo_as_0119
      CreateHists("PDF13000"); //CT14nnlo
      CreateHists("PDF25300"); //MMHT2014nnlo68cl
    }
    else if (meta_dsid==345775){
      for (int pdf=260001;pdf<=260100;++pdf) {CreateHists(Form("PDF%i",pdf));}
    }
    for (int is=0; is<6; is++) {CreateHists(Form("SCALE%i",is));}
  }

  // Start loop over systematics
  for (auto sys : getSystematics()) {
    //std::cout<< "sys " << sys.name() << " " << isSystematicAvailable(sys) << std::endl;
    CreateHists(sys.name());

    // Some MxAODs don't have systematics, skip this!
    if (!isSystematicAvailable(sys)) {
      continue;
    }

    CP_CHECK("HGamExample", applySystematicVariation(sys));

    // Check whether basic event selections are passed
    //if (!eventHandler()->pass()) { continue; } // GRL, trig, etc.

    // Check whether event passes HGam kinematic selections
    if (HG::isMAOD()) { // MxAODs store flags
      //if (!pass()) { continue; } // kinematic selection

      // This will update pileup weights, even if no objects are passed
      setSelectedObjects();
      doSelAndSaveHists(sys.name());
    }

    //if (sys.name()!="") break; //don't run systematics
  }

  // Return to nominal for truth event weight()
  CP_CHECK("HGamExample", applySystematicVariation(CP::SystematicSet()));

  // Fill histograms using nominal information
  ATH_MSG_DEBUG("Using nominal to fill histograms");

  m_nEvents[id]++;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZyUnfoldingInput::doSelAndSaveHists(TString sysname)
{

  TString prefix = HG::isData() ? "data" : HG::mcType()+getMCSampleName();
  if (sysname!="") prefix += "_" + sysname;
 
  const int nVariables=7;
  const char *variables[nVariables] = {"pty", "etay", "mlly", "dphilly", "ptlly", "ptmlly", "etcone20"};

  // final weight
  double w = HG::isData() ? 1.0 : lumiXsecWeight() * weight(); //lumiXsecWeight contains sum of weights from cutflow histogram
 
  bool isPassedZy = eventHandler()->getVar<char>("isPassedZy");
  bool isISR = eventHandler()->getVar<float>("m_ll")+eventHandler()->getVar<float>("m_lly")>182e3;

  double pt_y=-99, eta_y=-99, m_lly=-99, pt_lly=-99, dphi_lly=-99, ptm_lly=-99;
  pt_y     = eventHandler()->getVar<float>("pT_y1")/ HG::GeV;
  eta_y    = std::abs(eventHandler()->getVar<float>("eta_y1"));
  m_lly    = eventHandler()->getVar<float>("m_lly")/ HG::GeV;
  dphi_lly = std::abs(eventHandler()->getVar<float>("deltaPhi_ll_y"));
  pt_lly   = eventHandler()->getVar<float>("pt_lly")/ HG::GeV;
  ptm_lly   = pt_lly/m_lly;

  vals[0]=pt_y;
  vals[1]=eta_y;
  vals[2]=m_lly;
  vals[3]=dphi_lly;
  vals[4]=pt_lly;
  vals[5]=ptm_lly;
  vals[6]=0;

  int meta_dsid =366146;
  int dsid = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
  if (dsid>=345775 && dsid<=345782) meta_dsid=345775;

  if (HG::isMC() && sysname=="") { //truth
     
    truew =  lumiXsecWeight() * eventHandler()->mcWeight();

    const xAOD::TruthEventContainer* truthEvent = 0;
    if( ! event()->retrieve( truthEvent, "TruthEvents").isSuccess() ){
      Error("execute()", "Failed to retrieve TruthEvent. Exiting." );
      return EL::StatusCode::FAILURE;
    }
    xAOD::TruthEventContainer::const_iterator itr = truthEvent->begin();
    m_pdfw.clear();
    m_scalew.clear();
    double pdfw, scalew;
    int index;
    int dsid = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
    if (dsid>=345775 && dsid<=345782) meta_dsid=345775;
    if (meta_dsid==366146){
      for (int id=0;id<=104;++id) {
        int pdf = 261000+id;
        if (id==101) pdf = 269000;
        else if (id==102) pdf = 270000;
        else if (id==103) pdf = 13000;
        else if (id==104) pdf = 25300;
        index = m_variations[meta_dsid][Form("MUR1_MUF1_PDF%i", pdf)];
        pdfw = ((*itr)->weights())[index];
        m_pdfw.push_back(pdfw);
      }
      index = m_variations[meta_dsid]["MUR1_MUF0.5_PDF261000"]  ; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["MUR1_MUF2_PDF261000"]    ; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["MUR0.5_MUF1_PDF261000"]  ; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["MUR0.5_MUF0.5_PDF261000"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["MUR2_MUF1_PDF261000"]    ; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["MUR2_MUF2_PDF261000"]    ; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
    }
    else if (meta_dsid==345775) {
      for (int id=1;id<=100;++id) {
        index = m_variations[meta_dsid]["muR=0.10000E+01 muF=0.10000E+01"];
        pdfw = ((*itr)->weights())[index];
        m_pdfw.push_back(pdfw);
      }
      index = m_variations[meta_dsid]["muR=0.10000E+01 muF=0.10000E+01"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["muR=0.10000E+01 muF=0.20000E+01"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["muR=0.10000E+01 muF=0.50000E+00"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["muR=0.20000E+01 muF=0.10000E+01"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["muR=0.20000E+01 muF=0.20000E+01"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["muR=0.50000E+00 muF=0.10000E+01"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
      index = m_variations[meta_dsid]["muR=0.50000E+00 muF=0.50000E+00"]; scalew = ((*itr)->weights())[index]; m_scalew.push_back(scalew);
    }

    xAOD::TruthParticleContainer truth_photons   = truthHandler()->getPhotons();
    xAOD::TruthParticleContainer truth_electrons = truthHandler()->getElectrons();
    xAOD::TruthParticleContainer truth_muons     = truthHandler()->getMuons();
    xAOD::TruthParticleContainer fid_photons   = truthHandler()->applyPhotonSelection   (truth_photons);
    xAOD::TruthParticleContainer fid_electrons = xAOD::TruthParticleContainer(SG::VIEW_ELEMENTS);
    xAOD::TruthParticleContainer fid_muons = xAOD::TruthParticleContainer(SG::VIEW_ELEMENTS);
    double truem_ll=-99, truem_lly=-99, truept_y=-99, trueeta_y=-99, truedphi_lly=-99, truept_lly=-99, trueetcone20=-99, trueptm_lly=-99;
    isFiducial=false;

    for (auto muon : truth_muons){
      if (muon->auxdata<float>("pt_dressed")<25.*HG::GeV) continue;
      if (std::abs(muon->auxdata<float>("eta_dressed"))>2.47) continue;
      fid_muons.push_back(muon);
    }
    for (auto electron : truth_electrons){
      if (electron->auxdata<float>("pt_dressed")<25.*HG::GeV) continue;
      if (std::abs(electron->auxdata<float>("eta_dressed"))>2.47) continue;
      fid_electrons.push_back(electron);
    }
    if (fid_muons.size()>2 || fid_electrons.size()>2){std::cout << "more than 2 leptons" << std::endl;} 

    TLorentzVector lep1(0., 0., 0., 0.), lep2(0., 0., 0., 0.);
    if (fid_muons.size()>=2 && (fid_muons[0]->auxdata<float>("pt_dressed")>30.*HG::GeV || fid_muons[1]->auxdata<float>("pt_dressed")>30.*HG::GeV)){
      lep1.SetPtEtaPhiE(fid_muons[0]->auxdata<float>("pt_dressed"),fid_muons[0]->auxdata<float>("eta_dressed"),fid_muons[0]->auxdata<float>("phi_dressed"),fid_muons[0]->auxdata<float>("e_dressed")); 
      lep2.SetPtEtaPhiE(fid_muons[1]->auxdata<float>("pt_dressed"),fid_muons[1]->auxdata<float>("eta_dressed"),fid_muons[1]->auxdata<float>("phi_dressed"),fid_muons[1]->auxdata<float>("e_dressed"));
    }else if (fid_electrons.size()>=2 && (fid_electrons[0]->auxdata<float>("pt_dressed")>30.*HG::GeV || fid_electrons[1]->auxdata<float>("pt_dressed")>30.*HG::GeV)){
      lep1.SetPtEtaPhiE(fid_electrons[0]->auxdata<float>("pt_dressed"),fid_electrons[0]->auxdata<float>("eta_dressed"),fid_electrons[0]->auxdata<float>("phi_dressed"),fid_electrons[0]->auxdata<float>("e_dressed"));
      lep2.SetPtEtaPhiE(fid_electrons[1]->auxdata<float>("pt_dressed"),fid_electrons[1]->auxdata<float>("eta_dressed"),fid_electrons[1]->auxdata<float>("phi_dressed"),fid_electrons[1]->auxdata<float>("e_dressed"));
    }
    TLorentzVector gamma(0., 0., 0., 0.);
    //double UE;
    for (auto pho : fid_photons){
      if (std::abs(pho->eta())>2.37) continue;
      if (pho->p4().DeltaR(lep1)<0.4 || pho->p4().DeltaR(lep2)<0.4) continue; 
      if (pho->pt()<15*HG::GeV) continue;
      trueetcone20 =pho->auxdata<float>("etcone20")/pho->pt();
      if (pho->auxdata<float>("etcone20")/pho->pt()>0.07) continue;
      gamma = pho->p4();
      //if (std::abs(pho->eta())<1.5) UE = truthHandler()->centralEventShapeDensity()*(TMath::Pi()*0.2*0.2 - 5*7*0.025*TMath::Pi()/128);
      //else UE = truthHandler()->forwardEventShapeDensity()*(TMath::Pi()*0.2*0.2 - 5*7*0.025*TMath::Pi()/128);
      //std::cout << "etcone20 " << pho->auxdata<float>("etcone20") << " UE-subtracted etcone20 " << pho->auxdata<float>("etcone20")-UE << std::endl;
      //trueetcone20 = (pho->auxdata<float>("etcone20")-UE)/pho->pt(); 
      break;
    }
    if (lep1.Pt()>0 && gamma.Pt()>0){
      if (lep1.Pt()==0 || lep1.Eta()==0 || lep1.Phi()==0 || lep1.E()==0 ) std::cout << "zero in lepton four momentum" << std::endl;
      truept_y = gamma.Pt() / HG::GeV;
      trueeta_y = std::abs(gamma.Eta());
      truem_ll = (lep1+lep2).M() /HG::GeV;
      truem_lly = (lep1+lep2+gamma).M() /HG::GeV;
      truedphi_lly = std::abs(gamma.DeltaPhi(lep1+lep2));
      truept_lly = (lep1+lep2+gamma).Pt() / HG::GeV;
      trueptm_lly = truept_lly/truem_lly;
    }
    if (truem_ll>40 && truem_ll+truem_lly>182) isFiducial=true;

    truevals[0]=truept_y;
    truevals[1]=trueeta_y;
    truevals[2]=truem_lly;
    truevals[3]=truedphi_lly;
    truevals[4]=truept_lly;
    truevals[5]=trueptm_lly;
    truevals[6]=trueetcone20;
  }

  for (int ivar=0; ivar<nVariables; ivar++){
    bool isPassed = isPassedZy && isISR;
    bool isFid = isFiducial;
    if (ivar!=0){
      isPassed = isPassedZy && isISR && eventHandler()->getVar<float>("pT_y1")>30e3;
      isFid = isFiducial && truevals[0]>30;
    }
    if (isFid==true && isPassed==true){
       histoStore()->fillTH1F(prefix + Form("_fidcorrnum_%s", variables[ivar]), vals[ivar], w);
       histoStore()->fillTH1F(prefix + Form("_effcorrnum_%s", variables[ivar]), truevals[ivar], w);
       histoStore()->fillTH2F(prefix + Form("_migration_matrix_%s", variables[ivar]), vals[ivar], truevals[ivar], w);
    }
    if (isPassed==true){
       histoStore()->fillTH1F(prefix + Form("_A_%s", variables[ivar]), vals[ivar], w);
       histoStore()->fillTH1F(prefix + Form("_fidcorrden_%s", variables[ivar]), vals[ivar], w);
    }
    if (isFid==true){
      histoStore()->fillTH1F(prefix + Form("_effcorrden_%s", variables[ivar]), truevals[ivar], truew);
      histoStore()->fillTH1F(prefix + Form("_truth_%s", variables[ivar]), truevals[ivar], truew);
      if (sysname=="") {
        for (int id=1;id<=104;++id) {
          int pdf = 261000+id;
          if (id==101) pdf = 269000;
          else if (id==102) pdf = 270000;
          else if (id==103) pdf = 13000;
          else if (id==104) pdf = 25300;
          if (meta_dsid==345775){
            pdf = 260000+id;
            if (id>100) break;
          }
          double pdftruew = truew*m_pdfw[id]/m_pdfw[0];
          histoStore()->fillTH1F(prefix + Form("_PDF%i", pdf) + Form("_truth_%s", variables[ivar]), truevals[ivar], pdftruew);
        }
        for (int is=0; is<6; is++){
          double scaletruew = truew*m_scalew[is]/m_pdfw[0];
          histoStore()->fillTH1F(prefix + Form("_SCALE%i", is) + Form("_truth_%s", variables[ivar]), truevals[ivar], scaletruew);
        }
      }
    }
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZyUnfoldingInput::CreateHists(TString sysname) 
{

  TString prefix = HG::isData() ? "data" : HG::mcType()+getMCSampleName();
  if (sysname!="") prefix += "_" + sysname;

  const std::vector<double> ptybins = { 15., 20., 30., 40., 50., 60., 70., 80., 100., 120., 150., 200., 300., 500., 1200.};
  const std::vector<double> etaybins = {0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};
  const std::vector<double> mllybins = {95., 120., 130., 140., 150., 160., 170., 190., 220., 250., 300., 400., 500., 700., 2500.};
  const std::vector<double> dphillybins = {0., 0.1*TMath::Pi(), 0.2*TMath::Pi(), 0.3*TMath::Pi(), 0.4*TMath::Pi(), 0.5*TMath::Pi(), 0.6*TMath::Pi(), 0.7*TMath::Pi(), 0.8*TMath::Pi(), 0.9*TMath::Pi(), 0.95*TMath::Pi(), TMath::Pi()};
  const std::vector<double> ptllybins = {0., 5., 10., 15., 20., 30., 40., 50., 70., 100., 200., 1500.};
  const std::vector<double> ptmllybins = {0., 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 5.0};
  const std::vector<double> etconebins = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 1};
  const int nVariables=7;
  const std::vector<double> bins[nVariables] = {ptybins, etaybins, mllybins, dphillybins, ptllybins, ptmllybins, etconebins};

  const char *variables[nVariables] = {"pty", "etay", "mlly", "dphilly", "ptlly", "ptmlly", "etcone20"};
  const char *histtitles[nVariables] = {"#it{p}_{T#gamma} [GeV]", "|#it{#eta}_{#gamma}|", "#it{m}_{ll#gamma} [GeV]", "#Delta #phi(ll, #gamma)", "#it{p}_{Tll#gamma} [GeV]", "#it{p}_{Tll#gamma}/#it{m}_{ll#gamma}", "E_{T}^{cone20}/p_{T#gamma}"};
  int nbinsx[nVariables]={240, 25, 500, 20, 300, 100, 100};
  double xup[nVariables]={1200, 2.5, 2500, TMath::Pi(), 1500, 5, 1};

  int id = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
  if (m_nEvents[id] == 0) {
    // New sample!
    for (int ivar=0; ivar<nVariables; ivar++){
      if (ivar==0) std::cout << "creating " << prefix << std::endl;
      createAndRegisterTH1F(prefix + Form("_A_%s", variables[ivar]), nbinsx[ivar], 0, xup[ivar], prefix +";"+ histtitles[ivar]);
      createAndRegisterTH1F(prefix + Form("_fidcorrden_%s", variables[ivar]), nbinsx[ivar], 0, xup[ivar], prefix + ";"+ histtitles[ivar]);
      if (HG::isMC()) {
        createAndRegisterTH1F(prefix + Form("_truth_%s", variables[ivar]), nbinsx[ivar], 0, xup[ivar], prefix +";"+ histtitles[ivar]);
        createAndRegisterTH1F(prefix + Form("_fidcorrnum_%s", variables[ivar]), nbinsx[ivar], 0, xup[ivar], prefix + ";"+ histtitles[ivar]);
        createAndRegisterTH1F(prefix + Form("_effcorrnum_%s", variables[ivar]), nbinsx[ivar], 0, xup[ivar], prefix +";"+ histtitles[ivar]); 
        createAndRegisterTH1F(prefix + Form("_effcorrden_%s", variables[ivar]), nbinsx[ivar], 0, xup[ivar], prefix + ";"+ histtitles[ivar]);
        createAndRegisterTH2F(prefix + Form("_migration_matrix_%s", variables[ivar]), bins[ivar], bins[ivar], prefix +";"+ histtitles[ivar]);
      }
    }
  }
  return EL::StatusCode::SUCCESS;
}
