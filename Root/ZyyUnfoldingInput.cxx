#include "HGamGamStar/ZyyUnfoldingInput.h"
#include "HGamAnalysisFramework/HGamCommon.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>

// this is needed to distribute the algorithm to the workers
ClassImp(ZyyUnfoldingInput)

ZyyUnfoldingInput::ZyyUnfoldingInput(const char *name)
  : HgammaAnalysis(name) {}

// Here you delete any memory you allocated during your analysis.
ZyyUnfoldingInput::~ZyyUnfoldingInput() { }

EL::StatusCode ZyyUnfoldingInput::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  return EL::StatusCode::SUCCESS;
}

TH1F *ZyyUnfoldingInput::createAndRegisterTH1F(TString name, int Nbins, double min, double max, TString title)
{
  histoStore()->createTH1F(name, Nbins, min, max, title);
  wk()->addOutput(histoStore()->getTH1F(name));
  return histoStore()->getTH1F(name);
}

TH1F *ZyyUnfoldingInput::createAndRegisterTH1F(TString name, const std::vector<double>& xbins, TString title)
{
  histoStore()->createTH1F(name, xbins, title);
  wk()->addOutput(histoStore()->getTH1F(name));
  return histoStore()->getTH1F(name);
}

TH2F *ZyyUnfoldingInput::createAndRegisterTH2F(TString name, const std::vector<double>& xbins, const std::vector<double>& ybins, TString title)
{
  histoStore()->createTH2F(name, xbins, ybins, title);
  wk()->addOutput(histoStore()->getTH2F(name));
  return histoStore()->getTH2F(name);
}

EL::StatusCode ZyyUnfoldingInput::execute()
{
  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  if (!HG::isMAOD()) { HG::fatal("This code only works on MxAOD input."); }

  int id = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
  
  //hardcoding PDF weights since MetaData missing in MxAOD
  //int meta_dsid =364865; //take MetaData info from DSID366146 as default (Sherpa 2.2.4 sample)
  //if (id>=345775 && id<=345782) meta_dsid=345775; //or Madgraph
  #include "data/meta_364865.h"
  //#include "data/meta_345775.h"
  if (HG::isMC()){
    
    for (int pdf=261001;pdf<=261100;++pdf) {CreateHists(Form("PDF%i",pdf));}
    CreateHists("PDF269000"); //NNPDF30_nnlo_as_0117
    CreateHists("PDF270000"); //NNPDF30_nnlo_as_0119
    CreateHists("PDF13000"); //CT14nnlo
    CreateHists("PDF25300"); //MMHT2014nnlo68cl
    
    for (int is=0; is<6; is++) {CreateHists(Form("SCALE%i",is));}
  }
  

  // Start loop over systematics
  for (auto sys : getSystematics()) {
    //std::cout<< "sys " << sys.name() << " " << isSystematicAvailable(sys) << std::endl;
    //if (sys.name()!="") break; //don't run systematics

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

    
  }
  

  
  // Return to nominal for truth event weight()
  CP_CHECK("HGamExample", applySystematicVariation(CP::SystematicSet()));

  // Fill histograms using nominal information
  ATH_MSG_DEBUG("Using nominal to fill histograms");

  m_nEvents[id]++;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZyyUnfoldingInput::doSelAndSaveHists(TString sysname)
{

  TString prefix = HG::isData() ? "data" : HG::mcType()+getMCSampleName();
  if (sysname!="") prefix += "_" + sysname;
 
  const int nVariables=6;
  const char *variables[nVariables] = {"pty1", "pty2", "ptll", "ptllyy", "myy", "mllyy"};

  // final weight
  double w = HG::isData() ? 1.0 : lumiXsecWeight() * weight(); //lumiXsecWeight contains sum of weights from cutflow histogram

  //std::cout << "Weight " << w << std::endl;
 
  bool isPassedZyy = eventHandler()->getVar<char>("isPassedZyy");
  bool isPassedpT = eventHandler()->getVar<float>("pT_y1")>20e3 && eventHandler()->getVar<float>("pT_y2")>20e3;

  double pt_y1=-99, pt_y2=-99, pt_ll=-99, pt_llyy=-99, m_yy=-99, m_llyy=-99;
  pt_y1     = eventHandler()->getVar<float>("pT_y1")/ HG::GeV;
  pt_y2     = eventHandler()->getVar<float>("pT_y2")/ HG::GeV;
  pt_ll     = eventHandler()->getVar<float>("pt_ll")/ HG::GeV;
  pt_llyy   = eventHandler()->getVar<float>("pt_llyy")/ HG::GeV;
  m_yy      = eventHandler()->getVar<float>("m_yy")/ HG::GeV;
  m_llyy    = eventHandler()->getVar<float>("m_llyy")/ HG::GeV;

  vals[0]=pt_y1;
  vals[1]=pt_y2;
  vals[2]=pt_ll;
  vals[3]=pt_llyy;
  vals[4]=m_yy;
  vals[5]=m_llyy;

  /*int meta_dsid =366146;
  int dsid = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
  if (dsid>=345775 && dsid<=345782) meta_dsid=345775;*/

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
    //int dsid = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
    //if (dsid>=345775 && dsid<=345782) meta_dsid=345775;
    int meta_dsid = 364865;

    if (meta_dsid==364865){
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

    
    /*else if (meta_dsid==345775) {
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
      }*/
    

    xAOD::TruthParticleContainer truth_photons   = truthHandler()->getPhotons();
    xAOD::TruthParticleContainer truth_electrons = truthHandler()->getElectrons();
    xAOD::TruthParticleContainer truth_muons     = truthHandler()->getMuons();
    xAOD::TruthParticleContainer fid_photons   = truthHandler()->applyPhotonSelection(truth_photons);
    xAOD::TruthParticleContainer fid_electrons = xAOD::TruthParticleContainer(SG::VIEW_ELEMENTS);
    xAOD::TruthParticleContainer fid_muons = xAOD::TruthParticleContainer(SG::VIEW_ELEMENTS);
    xAOD::TruthParticleContainer fid_sel_photons = xAOD::TruthParticleContainer(SG::VIEW_ELEMENTS);
    double truept_y1=-99, truept_y2=-99, truept_ll=-99, truept_llyy=-99, truem_yy=-99, truem_llyy=-99, truem_ll=-99, truem_lly1=-99, truem_lly2=-99;
    isFiducial=false;

    for (auto muon : truth_muons){
      if (muon->auxdata<float>("pt_dressed")<20.*HG::GeV) continue;
      if (std::abs(muon->auxdata<float>("eta_dressed"))>2.47) continue;
      fid_muons.push_back(muon);
    }
    for (auto electron : truth_electrons){
      if (electron->auxdata<float>("pt_dressed")<20.*HG::GeV) continue;
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
    TLorentzVector gamma1(0., 0., 0., 0.), gamma2(0., 0., 0., 0.);
    //double UE;
    for (auto pho : fid_photons){
      if (std::abs(pho->eta())>2.37) continue;
      if (pho->p4().DeltaR(lep1)<0.4 || pho->p4().DeltaR(lep2)<0.4) continue; 
      if (pho->pt()<20.*HG::GeV) continue;
      //trueetcone20 =pho->auxdata<float>("etcone20")/pho->pt();
      if (pho->auxdata<float>("etcone20")/pho->pt()>0.07) continue;
      fid_sel_photons.push_back(pho);
      //if (std::abs(pho->eta())<1.5) UE = truthHandler()->centralEventShapeDensity()*(TMath::Pi()*0.2*0.2 - 5*7*0.025*TMath::Pi()/128);
      //else UE = truthHandler()->forwardEventShapeDensity()*(TMath::Pi()*0.2*0.2 - 5*7*0.025*TMath::Pi()/128);
      //std::cout << "etcone20 " << pho->auxdata<float>("etcone20") << " UE-subtracted etcone20 " << pho->auxdata<float>("etcone20")-UE << std::endl;
      //trueetcone20 = (pho->auxdata<float>("etcone20")-UE)/pho->pt(); 
    }

    if(fid_sel_photons.size()>=2){
      gamma1 = fid_sel_photons[0]->p4();
      gamma2 = fid_sel_photons[1]->p4();
    
    
      if(lep1.Pt()>0 && lep2.Pt()>0 && gamma1.Pt()>0 && gamma2.Pt()>0){
	truept_y1 = std::max(gamma1.Pt(),gamma2.Pt())/HG::GeV;
	truept_y2 = std::min(gamma1.Pt(),gamma2.Pt())/HG::GeV;
	truept_ll = (lep1+lep2).Pt()/HG::GeV;
	truept_llyy = (lep1+lep2+gamma1+gamma2).Pt()/HG::GeV;
	truem_yy = (gamma1+gamma2).M()/HG::GeV;
	truem_llyy = (lep1+lep2+gamma1+gamma2).M()/HG::GeV;
	
	truem_ll = (lep1+lep2).M() /HG::GeV;
	truem_lly1 = (lep1+lep2+gamma1).M() /HG::GeV;
	truem_lly2 = (lep1+lep2+gamma2).M() /HG::GeV;
	
	
	

	if(gamma1.DeltaR(gamma2)>0.4 && truem_ll>40 && truem_ll+std::min(truem_lly1,truem_lly2)>182)isFiducial=true;
	if(isFiducial)std::cout << "is fid" << std::endl;
      }
    }
    /*if (lep1.Pt()>0 && gamma.Pt()>0){
      if (lep1.Pt()==0 || lep1.Eta()==0 || lep1.Phi()==0 || lep1.E()==0 ) std::cout << "zero in lepton four momentum" << std::endl;
      truept_y = gamma.Pt() / HG::GeV;
      trueeta_y = std::abs(gamma.Eta());
      truem_ll = (lep1+lep2).M() /HG::GeV;
      truem_lly = (lep1+lep2+gamma).M() /HG::GeV;
      truedphi_lly = std::abs(gamma.DeltaPhi(lep1+lep2));
      truept_lly = (lep1+lep2+gamma).Pt() / HG::GeV;
      trueptm_lly = truept_lly/truem_lly;
      }*/
    //if (truem_ll>40 && truem_ll+truem_lly>182) isFiducial=true;

    truevals[0]=truept_y1;
    truevals[1]=truept_y2;
    truevals[2]=truept_ll;
    truevals[3]=truept_llyy;
    truevals[4]=truem_yy;
    truevals[5]=truem_llyy;
    
  }

  for (int ivar=0; ivar<nVariables; ivar++){
    bool isPassed = isPassedZyy && isPassedpT;
    bool isFid = isFiducial;
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
          /*if (meta_dsid==345775){
            pdf = 260000+id;
            if (id>100) break;
	    }*/
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

EL::StatusCode ZyyUnfoldingInput::CreateHists(TString sysname) 
{

  TString prefix = HG::isData() ? "data" : HG::mcType()+getMCSampleName();
  if (sysname!="") prefix += "_" + sysname;

  const std::vector<double> pty1bins = {20., 40., 60., 100., 2000.};
  const std::vector<double> pty2bins = {20., 25., 30., 35., 50., 1200.};
  const std::vector<double> ptllbins = {0., 30., 60., 120., 1200.};
  const std::vector<double> ptllyybins = {0., 10., 30., 50., 100., 1200.};
  const std::vector<double> myybins = {0., 50., 75., 100., 150., 3000.};
  const std::vector<double> mllyybins = {150., 200., 250., 350., 500., 3000.};
  
  const int nVariables=6;
  const std::vector<double> bins[nVariables] = {pty1bins, pty2bins, ptllbins, ptllyybins, myybins, mllyybins};

  const char *variables[nVariables] = {"pty1", "pty2", "ptll", "ptllyy", "myy", "mllyy"};
  const char *histtitles[nVariables] = {"#it{p}_{T}^{#gamma1} [GeV]", "#it{p}_{T}^{#gamma2} [GeV]", "#it{p}_{T}^{ll} [GeV]", "#it{p}_{T}^{ll#gamma#gamma} [GeV]", "#it{m}_{#gamma#gamma} [GeV]", "#it{m}_{ll#gamma#gamma} [GeV]"};
  //int nbinsx[nVariables]={240, 25, 500, 20, 300, 100, 100};
  //double xup[nVariables]={1200, 2.5, 2500, TMath::Pi(), 1500, 5, 1};

  int id = HG::isData() ? 0 : eventInfo()->mcChannelNumber();
  if (m_nEvents[id] == 0) {
    // New sample!
    for (int ivar=0; ivar<nVariables; ivar++){
      if (ivar==0) std::cout << "creating " << prefix << std::endl;
      createAndRegisterTH1F(prefix + Form("_A_%s", variables[ivar]), bins[ivar], prefix +";"+ histtitles[ivar]);
      createAndRegisterTH1F(prefix + Form("_fidcorrden_%s", variables[ivar]), bins[ivar], prefix + ";"+ histtitles[ivar]);
      if (HG::isMC()) {
        createAndRegisterTH1F(prefix + Form("_truth_%s", variables[ivar]), bins[ivar], prefix +";"+ histtitles[ivar]);
        createAndRegisterTH1F(prefix + Form("_fidcorrnum_%s", variables[ivar]), bins[ivar], prefix + ";"+ histtitles[ivar]);
        createAndRegisterTH1F(prefix + Form("_effcorrnum_%s", variables[ivar]), bins[ivar], prefix +";"+ histtitles[ivar]); 
        createAndRegisterTH1F(prefix + Form("_effcorrden_%s", variables[ivar]), bins[ivar], prefix + ";"+ histtitles[ivar]);
        createAndRegisterTH2F(prefix + Form("_migration_matrix_%s", variables[ivar]), bins[ivar], bins[ivar], prefix +";"+ histtitles[ivar]);
      }
    }
  }
  return EL::StatusCode::SUCCESS;
}
