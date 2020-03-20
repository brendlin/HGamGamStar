#include "HGamGamStar/ZyCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>

#include "PhotonVertexSelection/PhotonPointingTool.h"
#include "ZMassConstraint/ConstraintFit.h"

#include "HGamGamStar/HggStarVariables.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZyCutflowAndMxAOD)

ZyCutflowAndMxAOD::ZyCutflowAndMxAOD(const char *name)
: MxAODTool(name), m_goodFakeComb(false) { }

ZyCutflowAndMxAOD::~ZyCutflowAndMxAOD() {}


EL::StatusCode ZyCutflowAndMxAOD::createOutput()
{
  // Read the output branch names - add option to make this configurable in future ?
  m_photonContainerName = "HGam"+config()->getStr("PhotonHandler.ContainerName");
  m_jetContainerName    = "HGam"+config()->getStr("JetHandler.ContainerName");
  m_elecContainerName   = "HGam"+config()->getStr("ElectronHandler.ContainerName");
  m_muonContainerName   = "HGam"+config()->getStr("MuonHandler.ContainerName");
  m_evtInfoName         = "EventInfo";
  m_truthEvtsName       = "TruthEvents";

  // What selection to require. For defintion, see CutEnum in the header file
  // Negative values means no event selection
  m_skimCut = config()->getInt("SkimmingCut", -1);

  // Whether or not to apply systematic variations
  // If true, all variables will be stored
  m_applySystematics = config()->getBool("ApplySystematicVariations",false);

  if (m_applySystematics && m_skimCut >= 1)
    HG::fatal("Running over systematics and applying a skimming cut is currently not possible.");

  // Whether to allow more than two good photons
  m_allowMoreThanTwoPhotons = config()->getBool("AllowMoreThanTwoPhotons",false);

  // Whether to save objects (photons, jets ...)
  m_saveObjects = config()->getBool("SaveObjects",false);

  // Whether to save the list of differential variables
  m_saveDetailed = config()->getBool("SaveDetailedVariables",false);

  // Whether to save the truth objects and differential variables
  m_saveTruthObjects = HG::isMC() && config()->getBool("SaveTruthObjects",false);
  m_saveTruthVars    = HG::isMC() && config()->getBool("SaveTruthVariables",false);

  //Save fake photon combinations
  m_enableFakePhotons = HG::isMC() && config()->getBool("SaveFakePhotonCombinations", false);

  //select opposite flavour instead of same flavour
  m_checkemu = config()->getBool("SelectOppFlavour", false);

  // Ignore photon selection and save inclusive Z events
  m_saveAllZ = config()->getBool("SaveAllZ",false);

  // do VBS selection
  m_isVBSsel = config()->getBool("isVBSselection",false);

  // Whether we are running with yybb-tool in detailed mode.
  m_detailedHHyybb = config()->getBool("HHyybb.DetailedInfo",false);

  // Temporary hack for large PhotonAllSys samples
  m_photonAllSys = config()->getStr("PhotonHandler.Calibration.decorrelationModel") == "FULL_v1";

  //Triboson selection
  m_isWZysel = config()->getBool("WZyConfig",false);
  m_isZZysel = config()->getBool("ZZyConfig",false);

  // a. Event variables
  StrV ignore = {};
  if (HG::isData()) ignore = {".mcChannelNumber", ".mcEventWeights", ".RandomRunNumber", ".truthCentralEventShapeDensity", ".truthForwardEventShapeDensity"};

  StrV trigs = config()->getStrV("EventHandler.RequiredTriggers");
  StrV extra = {};
  for (auto trig: trigs) extra.push_back(".passTrig_"+trig);

  declareOutputVariables(m_evtInfoName,"MxAOD.Variables.EventInfo", extra, ignore);

  // a.2 TruthEvents variables
  ignore = {};
  extra = {};
  declareOutputVariables(m_truthEvtsName,"MxAOD.Variables.TruthEvents", extra, ignore);

  // a.3 BTagging variables - If 20.7 switch decoration names.
  TString mv2_tagger = "MV2c10";
  #ifdef __Rel20p1__
    mv2_tagger = "MV2c20";
    TString vars = config()->getStr("MxAOD.Variables.Jet");
    vars.ReplaceAll("MV2c10","MV2c20");
    config()->setValue("MxAOD.Variables.Jet",vars);
  #endif

  // b. Selected objects

  // If we have a detailed run then append a list of extra variables to add to jets
  if(m_detailedHHyybb)
  {
    TString yybb_detailedJetVars = config()->getStr("MxAOD.Variables.Jet")+config()->getStr("MxAOD.yybb-Detailed.Jet");
    yybb_detailedJetVars.ReplaceAll(" ",""); //Remove spaces in lists...
    config()->setValue("MxAOD.Variables.Jet",yybb_detailedJetVars.Data());
  }

  if (HG::isData()) ignore = {".isEMTight_nofudge", ".isTight_nofudge", ".topoetcone20_DDcorrected", ".topoetcone40_DDcorrected", ".truthOrigin", ".truthType", ".truthConvRadius", ".scaleFactor", ".truthLink", ".parentPdgId", ".pdgId"};
  declareOutputVariables(m_photonContainerName, "MxAOD.Variables.Photon"  , {}, ignore);
  declareOutputVariables("HGamPhotonsWithFakes","MxAOD.Variables.Photon"  , {}, ignore);
  if (HG::isData()) ignore = {".SF_"+mv2_tagger+"_FixedCutBEff_60", ".SF_"+mv2_tagger+"_FixedCutBEff_70", ".SF_"+mv2_tagger+"_FixedCutBEff_77", ".SF_"+mv2_tagger+"_FixedCutBEff_85", ".Eff_"+mv2_tagger+"_FixedCutBEff_60", ".Eff_"+mv2_tagger+"_FixedCutBEff_70", ".Eff_"+mv2_tagger+"_FixedCutBEff_77", ".Eff_"+mv2_tagger+"_FixedCutBEff_85", ".InEff_"+mv2_tagger+"_FixedCutBEff_60", ".InEff_"+mv2_tagger+"_FixedCutBEff_70", ".InEff_"+mv2_tagger+"_FixedCutBEff_77", ".InEff_"+mv2_tagger+"_FixedCutBEff_85", ".HadronConeExclTruthLabelID"};
  declareOutputVariables(m_jetContainerName   , "MxAOD.Variables.Jet"     , {}, ignore);
  if (HG::isData()) ignore = {".scaleFactor", ".truthLink"};
  declareOutputVariables(m_elecContainerName  , "MxAOD.Variables.Electron", {}, ignore);
  if (HG::isData()) ignore = {".scaleFactor"};
  declareOutputVariables(m_muonContainerName  , "MxAOD.Variables.Muon"    , {}, ignore);
  declareOutputVariables("HGamMuonsInJets"    , "MxAOD.Variables.Muon"    , {}, ignore);

  // c. Truth objects
  if (HG::isMC()) {
    m_photonTruthContainerName = "HGam"+config()->getStr("TruthHandler.PhotonContainerName");
    m_elecTruthContainerName   = "HGam"+config()->getStr("TruthHandler.ElectronContainerName");
    m_muonTruthContainerName   = "HGam"+config()->getStr("TruthHandler.MuonContainerName");
    m_jetTruthContainerName    = "HGam"+config()->getStr("TruthHandler.JetContainerName");

    declareOutputVariables(m_photonTruthContainerName  , "MxAOD.Variables.TruthPhotons"    );
    declareOutputVariables(m_elecTruthContainerName    , "MxAOD.Variables.TruthElectrons"  );
    declareOutputVariables(m_muonTruthContainerName    , "MxAOD.Variables.TruthMuons"      );
    declareOutputVariables(m_jetTruthContainerName     , "MxAOD.Variables.TruthJets"       );
    declareOutputVariables("HGam"+config()->getStr("TruthHandler.HiggsBosonContainerName"), "MxAOD.Variables.TruthHiggsBosons");
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZyCutflowAndMxAOD::execute()
{
  // Needed for all underlaying tools to be working corectly!
  HgammaAnalysis::execute();

  // Handle File Metadata (need to put this here because we need sample ID to define cutflow histo
  if (m_newFileMetaData) {

    // Initialize cutflow histograms by calling them.
    getCutFlowHisto();
    if (HG::isMC()) {
      getCutFlowWeightedHisto();
    }

    // Fill the AOD and DAOD entries of the cutflow histograms.
    if (MxAODTool::fillCutFlowWithBookkeeperInfo() == EL::StatusCode::FAILURE)
    {
      return EL::StatusCode::FAILURE;
    }

    m_newFileMetaData = false;
  }

  // apply cuts. Returned value will be the last passed cut
  m_cutFlow = cutflow();

  // fill the cut-flow histograms up to tight selection
  double wi = weightInitial();
  for (int cut=ALLEVTS;cut<m_cutFlow;++cut) {
    if (cut <= GAM_TIGHTID)
      fillCutFlow(CutEnum(cut),wi);
  }

  // if desired, apply skimming.
  // That is, only write events that pass a given cut to the output
  if (m_skimCut >= 1 && m_cutFlow <= m_skimCut)
    return EL::StatusCode::SUCCESS;

  // Selects the objects, does overlap removal, and calculate all
  // variables that will be saved
  doReco();

  // fill the cut-flow histograms after tight selection including SF weights
  double w = weight();
  for (int cut=ALLEVTS;cut<m_cutFlow;++cut) {
    if (cut > GAM_TIGHTID)
      fillCutFlow(CutEnum(cut),w);
  }

  // check if we should apply systematics or not
  if (m_applySystematics) {
    for (auto sys: getSystematics()) {
      // ignore nominal case, already done!
      if (sys.name() == "") continue;

      // apply the systmeatic variation and calculate the outupt
      CP_CHECK("ZyCutflowAndMxAOD::execute()", applySystematicVariation(sys));
      m_cutFlow = cutflow();
      doReco(true);
    }
  }

  if ( HG::isMC() && (m_saveTruthObjects || m_saveTruthVars))
    doTruth();

  // Write nominal EventInfo to output
  eventHandler()->writeEventInfo();

  // Save all information written to output
  event()->fill();

  return EL::StatusCode::SUCCESS;
}

// Returns value of the last cut passed in the cut sequence
ZyCutflowAndMxAOD::CutEnum ZyCutflowAndMxAOD::cutflow()
{
  m_preSelPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selElectrons = xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
  m_selMuons = xAOD::MuonContainer(SG::VIEW_ELEMENTS);
  m_allJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);
  m_selJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);

  m_isZ = false;
  bool c_oneLooseGam=true, c_ambiguity=true, c_onePhotonPostor=true;

  //Check if there are two good fakes. Needed so we dont slim the event at trigger.
  m_goodFakeComb = false;
  if(HG::isMC() && m_enableFakePhotons){
      double weightFakePhotons = 1;
      xAOD::PhotonContainer photonsWithFakes = getFakePhotons(weightFakePhotons);
      m_goodFakeComb = photonsWithFakes.size()>1 ? true : false;
      if(m_goodFakeComb) return PASSALL;
  }

  //==== CUT 0 : Remove duplicate events (only for data) ====
  static bool checkDuplicates = config()->getBool("EventHandler.CheckDuplicates");
  if ( checkDuplicates && eventHandler()->isDuplicate() ) return DUPLICATE;

  //==== CUT 1 : GRL ====
  static bool requireGRL = config()->getBool("EventHandler.CheckGRL");
  if ( requireGRL && HG::isData() && !eventHandler()->passGRL(eventInfo()) ) return GRL;

  //==== CUT 2 : Require trigger ====
  static bool requireTrigger = config()->getBool("EventHandler.CheckTriggers");
  if ( requireTrigger && !eventHandler()->passTriggers() ) return TRIGGER;

  //==== CUT 3 : Detector quality ====
  if ( !(eventHandler()->passLAr (eventInfo()) &&
         eventHandler()->passTile(eventInfo()) &&
         eventHandler()->passSCT (eventInfo()) ) )
    return DQ;

  //==== CUT 4 : Require a vertex ====
  if ( !eventHandler()->passVertex(eventInfo()) ) return VERTEX;

  // retrieve electrons, muons
  m_allElectrons = electronHandler()->getCorrectedContainer();
  m_preSelElectrons = electronHandler()->applySelection(m_allElectrons);

  m_allMuons = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer dirtyMuons = muonHandler()->applySelection(m_allMuons);

  //==== CUT 5 : 2 SF leptons (before OR) ====
  if (m_checkemu && (m_preSelElectrons.size()==0 || dirtyMuons.size()==0)) return TWO_SF_LEPTONS;
  else if (!m_checkemu && m_preSelElectrons.size() < 2 && dirtyMuons.size() < 2) return TWO_SF_LEPTONS;

  // Get object containers
  m_allPhotons = photonHandler()->getCorrectedContainer();
  m_preSelPhotons = photonHandler()->applyPreSelection(m_allPhotons);
  if (m_preSelPhotons.size()  ) m_selPhotons.push_back(m_preSelPhotons[0]);
  if (m_preSelPhotons.size() > 1) { m_selPhotons.push_back(m_preSelPhotons[1]); }

  int nloose=0, namb=0, nHV=0;
  for (auto gam:m_allPhotons) {
    if (photonHandler()->passOQCut(gam)       &&
        photonHandler()->passCleaningCut(gam) &&
        photonHandler()->passPtEtaCuts(gam)   &&
        photonHandler()->passPIDCut(gam,egammaPID::PhotonIDLoose))
      ++nloose;

    if (photonHandler()->passOQCut(gam)       &&
        photonHandler()->passCleaningCut(gam) &&
        photonHandler()->passPtEtaCuts(gam)   &&
        photonHandler()->passPIDCut(gam,egammaPID::PhotonIDLoose) &&
        photonHandler()->passAmbCut(gam))
      ++namb;

    if (photonHandler()->passOQCut(gam)       &&
        photonHandler()->passCleaningCut(gam) &&
        photonHandler()->passPtEtaCuts(gam)   &&
        photonHandler()->passPIDCut(gam,egammaPID::PhotonIDLoose) &&
        photonHandler()->passAmbCut(gam)/*       &&
        photonHandler()->passHVCut(gam)*/)
      ++nHV;
  }

  //==== CUT 6 : Require one loose photons, pT>15 GeV ====
  if (nloose<1 && !m_saveAllZ) return ONE_LOOSE_GAM;
  if (nloose<1 && m_saveAllZ) c_oneLooseGam = false;

  //==== CUT 7 : Ambiguity
  // - Require two loose photons that also pass e-gamma ambiguity ====
  static bool requireAmbiguity = config()->getBool("PhotonHandler.Selection.ApplyAmbiguityCut", false);
  if (requireAmbiguity && namb<1 && !m_saveAllZ) return AMBIGUITY;
  if (requireAmbiguity && namb<1 && m_saveAllZ) c_ambiguity = false;

  // static bool requireHV = config()->getBool("PhotonHandler.Selection.ApplyHVCut", false);
  // if (requireHV && nHV<1) return AMBIGUITY;


  //==== CUT 8 : Two OSSF leptons (after OR) ====
  m_allJets = jetHandler()->getCorrectedContainer();
  m_selJets = jetHandler()->applySelection(m_allJets);

  // Removes overlap with candidate diphoton system, and any additional tight photons (if option set)
  overlapHandler()->removeOverlap(m_selPhotons, m_selJets, m_preSelElectrons, dirtyMuons);

  //above doesn't have option to remove photon overlapping with lepton
  overlapHandler()->removeOverlap(m_selPhotons, m_preSelElectrons, 0.4);
  overlapHandler()->removeOverlap(m_selPhotons, dirtyMuons, 0.4);
  overlapHandler()->removeOverlap(m_preSelElectrons, dirtyMuons, 0.2);

  // Muon cleaning should be done after overlap removal
  m_preSelMuons = muonHandler()->applyCleaningSelection(dirtyMuons);

  // Select Z candidate after overlap removal.
  // choose leading OSSF pair
  int nOSSFpair=0;

  // Container for all possible Z candidates with mll > 10GeV
  std::vector<std::pair<int,int>> lepPairs;
  int nElecPairs=0, nMuonPairs=0;

  if(!m_checkemu && m_preSelElectrons.size()>=2){
    for(int ilepton1=0; ilepton1<((int)m_preSelElectrons.size()-1); ilepton1++){
      for(int ilepton2=ilepton1+1; ilepton2< (int)m_preSelElectrons.size(); ilepton2++){
        if(m_preSelElectrons[ilepton1]->charge() + m_preSelElectrons[ilepton2]->charge() == 0){
          nOSSFpair++;
          m_selElectrons.push_back(m_preSelElectrons[ilepton1]);
          m_selElectrons.push_back(m_preSelElectrons[ilepton2]);

	  if((m_preSelElectrons[ilepton1]->p4()+m_preSelElectrons[ilepton2]->p4()).M()>10e3){
	    lepPairs.push_back(std::make_pair(ilepton1,ilepton2));
	    nElecPairs++;
	  }

        }
      }
    }
  }


  if(!m_checkemu && m_preSelMuons.size()>=2){
    for(int ilepton1=0; ilepton1<((int)m_preSelMuons.size()-1); ilepton1++){
      for(int ilepton2=ilepton1+1; ilepton2< (int)m_preSelMuons.size(); ilepton2++){
        if(m_preSelMuons[ilepton1]->charge() + m_preSelMuons[ilepton2]->charge() == 0){
          nOSSFpair++;
          m_selMuons.push_back(m_preSelMuons[ilepton1]);
          m_selMuons.push_back(m_preSelMuons[ilepton2]);

	  if((m_preSelMuons[ilepton1]->p4()+m_preSelMuons[ilepton2]->p4()).M()>10e3){
	    lepPairs.push_back(std::make_pair(ilepton1,ilepton2));
	    nMuonPairs++;
	  }
        }
      }
    }
  }


  // For triboson config, need multiple lepton pairs so save all preselected leptons
  SG::AuxElement::Accessor<int> isZLepton("isZLepton");

  if(m_isWZysel || m_isZZysel){
    m_selElectrons.clear();
    m_selElectrons = m_preSelElectrons;
    for(int ielec=0;ielec<m_selElectrons.size();ielec++) isZLepton(*m_selElectrons[ielec]) = -1;
    m_selMuons.clear();
    m_selMuons = m_preSelMuons;
    for(int imuon=0;imuon<m_selMuons.size();imuon++) isZLepton(*m_selMuons[imuon]) = -1;
  }
  
  //emu pair
  int nOSOFpair=0;
  if(m_checkemu && m_preSelMuons.size()>=1 && m_preSelElectrons.size()>=1){
    for(int ilepton1=0; ilepton1<((int)m_preSelMuons.size()); ilepton1++){
      for(int ilepton2=0; ilepton2< (int)m_preSelElectrons.size(); ilepton2++){
        if(m_preSelMuons[ilepton1]->charge() + m_preSelElectrons[ilepton2]->charge() == 0){
           nOSOFpair++;
           m_selMuons.push_back(m_preSelMuons[ilepton1]);
           m_selElectrons.push_back(m_preSelElectrons[ilepton2]);
        }
      }
    }
  }

  double  mZ=91187;
  double m_ll1=-99, m_ll2=-99;
  int Z_pair1=-1, Z_pair2=-1;
  int nQuad=0;

  // WZy Z selection
  if(m_isWZysel && (nElecPairs+nMuonPairs)>=1){
    // Loop over all possible lepton pairs and find one with mll closest to mZ
    for(int ipair=0; ipair<(nElecPairs+nMuonPairs); ipair++){
      double m_pair = 0.;
      //For muon pairs
      if(ipair>=nElecPairs){
	m_pair = (m_selMuons[lepPairs[ipair].first]->p4()+m_selMuons[lepPairs[ipair].second]->p4()).M(); 
      }
      //For electron pairs
      else{
	m_pair = (m_selElectrons[lepPairs[ipair].first]->p4()+m_selElectrons[lepPairs[ipair].second]->p4()).M(); 
      }
      //Choose pair with closest inv mass to mZ
      if(abs(m_pair-mZ)<abs(m_ll1-mZ)){
	m_ll1 = m_pair;
	Z_pair1 = ipair;
      }
    }
  }
  
  // ZZy Z selection
  if(m_isZZysel && (nElecPairs+nMuonPairs)>=2){
    // Loop over all combinations of lepton pairs (quads)
    for(int ipair=0; ipair<(nElecPairs+nMuonPairs-1); ipair++){
      for(int jpair=ipair+1; jpair<(nElecPairs+nMuonPairs); jpair++){
	TLorentzVector vi,vj;
	// If either of the two pairs contain a common lepton, continue (only applicable for comparing like flavour pairs)
	if((ipair<nElecPairs && jpair<nElecPairs) || (ipair>=nElecPairs && jpair>=nElecPairs)){
	  if(lepPairs[ipair].first==lepPairs[jpair].first || lepPairs[ipair].first==lepPairs[jpair].second || lepPairs[ipair].second==lepPairs[jpair].first || lepPairs[ipair].second==lepPairs[jpair].second)continue;
	}
	nQuad++;
	//Muon pair
	if(ipair>=nElecPairs){
	  vi = m_selMuons[lepPairs[ipair].first]->p4()+m_selMuons[lepPairs[ipair].second]->p4(); 
	}
	//Electron pair
	else{
	  vi = m_selElectrons[lepPairs[ipair].first]->p4()+m_selElectrons[lepPairs[ipair].second]->p4(); 
	}
	//Muon pair
	if(jpair>=nElecPairs){
	  vj = m_selMuons[lepPairs[jpair].first]->p4()+m_selMuons[lepPairs[jpair].second]->p4(); 
	}
	//Electrons pair
	else{
	  vj = m_selElectrons[lepPairs[jpair].first]->p4()+m_selElectrons[lepPairs[jpair].second]->p4(); 
	}
	// Choose quad with closest combined inv mass to 2mZ
	if( (abs(vi.M()-mZ)+abs(vj.M()-mZ)) < (abs(m_ll1-mZ)+abs(m_ll2-mZ)) ){
	  // Store selected pairs in order of pT, i.e. Z_pair1 is the harder Z
	  if(vi.Pt()>vj.Pt()){
	    m_ll1 = vi.M(); Z_pair1 = ipair;
	    m_ll2 = vj.M(); Z_pair2 = jpair;
	  }
	  else{
	    m_ll1 = vj.M(); Z_pair1 = jpair;
	    m_ll2 = vi.M(); Z_pair2 = ipair;
	  }
	}
      }
    }
  }

 
  // Decorate lepton containers with selected "Z" leptons 
  // isZlepton = -1 not a Z lepton
  // isZlepton = 1 leading Z lepton
  // isZlepton = 2 subleading Z lepton
  if(m_isWZysel || m_isZZysel){
    for(int ipair=0;ipair<nElecPairs+nMuonPairs;ipair++){
      if(ipair==Z_pair1){
	if(ipair>=nElecPairs){
	  isZLepton(*m_selMuons[lepPairs[Z_pair1].first])=1;
	  isZLepton(*m_selMuons[lepPairs[Z_pair1].second])=1;
	}
	else{
	  isZLepton(*m_selElectrons[lepPairs[Z_pair1].first])=1;
	  isZLepton(*m_selElectrons[lepPairs[Z_pair1].second])=1;
	}
      }
      if(ipair==Z_pair2){
	if(ipair>=nElecPairs){
	  isZLepton(*m_selMuons[lepPairs[Z_pair2].first])=2;
	  isZLepton(*m_selMuons[lepPairs[Z_pair2].second])=2;
	}
	else{
	  isZLepton(*m_selElectrons[lepPairs[Z_pair2].first])=2;
	  isZLepton(*m_selElectrons[lepPairs[Z_pair2].second])=2;
	}
      }
    }
  }

  
  if (m_checkemu && nOSOFpair==0) return TWO_SF_LEPTONS_POSTOR; 
  else if (m_isZZysel && nQuad==0) return TWO_SF_LEPTONS_POSTOR;
  else if (m_isWZysel && (m_selElectrons.size()+m_selMuons.size()<3 || nOSSFpair==0)) return TWO_SF_LEPTONS_POSTOR;
  else if (!m_isZZysel && !m_isWZysel && !m_checkemu && nOSSFpair==0) return TWO_SF_LEPTONS_POSTOR;


  //==== CUT 9 : 1 photon after OR ====
  if (m_selPhotons.size()==0 && !m_saveAllZ) return ONE_PHOTON_POSTOR;
  if (m_selPhotons.size()==0 && m_saveAllZ) c_onePhotonPostor=false;

  //==== CUT 10 : Trigger matching ====
  //trigger matching
  static bool requireTriggerMatch = config()->getBool("EventHandler.CheckTriggerMatching", true);
  //if ( requireTriggerMatch && !passTriggerMatch(NULL, &m_selElectrons, &m_selMuons, NULL) ) return TRIG_MATCH; //doesn't work
  if ( requireTriggerMatch){
    StrV m_requiredTriggers = config()->getStrV("EventHandler.RequiredTriggers");
    int itrigmatch=0;
    for (auto trig: m_requiredTriggers) {
      if (passTriggerMatch(trig, NULL, &m_selElectrons, &m_selMuons, NULL) ) itrigmatch++;
    }
    if (itrigmatch==0) return TRIG_MATCH;
  }

  //==== CUT 11 : 30 GeV cut on leading lepton ====
  if (!m_isVBSsel && !((m_selElectrons.size()>0 && m_selElectrons[0]->pt()>30*HG::GeV) || (m_selMuons.size()>0 && m_selMuons[0]->pt()>30*HG::GeV))) return LEADLEPTON_PT;

  //==== CUT 12 : MLL>40 GeV ====
  double m_ll=-99, m_emu=-99;
  if (!m_isWZysel && !m_isZZysel && m_selElectrons.size()>=2 ) m_ll = (m_selElectrons[0]->p4() + m_selElectrons[1]->p4()).M();
  else if (!m_isWZysel && !m_isZZysel && m_selMuons.size()>=2 ) m_ll = (m_selMuons[0]->p4() + m_selMuons[1]->p4()).M();


  if (m_selElectrons.size()>=1 && m_selMuons.size()>=1) m_emu = (m_selElectrons[0]->p4() + m_selMuons[0]->p4()).M();
  

  if (m_checkemu && m_emu<40*HG::GeV ) return MASSCUT;
  else if (m_isWZysel && m_ll1<40*HG::GeV) return MASSCUT;
  else if (m_isZZysel && (m_ll1<40*HG::GeV || m_ll2<40*HG::GeV)) return MASSCUT;
  else if (!m_isWZysel && !m_isZZysel && !m_checkemu && m_ll<40*HG::GeV ) return MASSCUT;

  
  // isZ flag set properly whilst not influincing the original cutflow decision
  m_isZ = true;
  if(m_saveAllZ){
    if(!c_oneLooseGam) return ONE_LOOSE_GAM;
    if(!c_ambiguity) return AMBIGUITY;
    if(!c_onePhotonPostor) return ONE_PHOTON_POSTOR;
  }
  //return GAM_TIGHTID;
  //==== CUT 13 : tight ID for photon ====

  if (!photonHandler()->passPIDCut(m_selPhotons[0])) return GAM_TIGHTID;

  //==== CUT 14 : pass FixedCutLoose isolation for photon ====

  if (!photonHandler()->passIsoCut(m_selPhotons[0], HG::Iso::FixedCutLoose)) return GAM_ISOLATION;

  //==== CUT 7 : Require two loose photons to pass trigger matching
  // static bool requireTriggerMatch = config()->getBool("EventHandler.CheckTriggerMatching", true);
  // if ( requireTriggerMatch && !passTriggerMatch(&loosePhotons) ) return TRIG_MATCH;

  // Our *Higgs candidate photons* are the two, leading pre-selected photons
  // These are also the photons used to define the primary vertex
  // xAOD::Photon* gam1 = loosePhotons[0];

  //==== CUT 8 : Require both photons to pass photon ID (isEM) ====
  // Do we really want to require the highest-pt photon to pass tight ID? Can we ask for lower-pt gam?
  // static bool requireTight = config()->getBool("PhotonHandler.Selection.ApplyPIDCut", true);
  // if (requireTight && (!photonHandler()->passPIDCut(gam1)) ) return GAM_TIGHTID;

  //==== CUT 9 : Require both photons to fulfill the isolation criteria ===
  // static bool requireIso = config()->getBool("PhotonHandler.Selection.ApplyIsoCut", true);
  // if (requireIso && (!photonHandler()->passIsoCut(gam1))) return GAM_ISOLATION;

  //==== CUT 11 : Mass window cut ====
  // if ( !passMyyWindowCut(loosePhotons) ) return MASSCUT;

  return PASSALL;
}

EL::StatusCode  ZyCutflowAndMxAOD::doReco(bool isSys){
  // Do anything you missed in cutflow, and save the objects.

  // Rebuild MET using selected objects
  m_allMET = etmissHandler()->getCorrectedContainer(&m_selPhotons    ,
                                                    &m_allJets     ,
                                                    &m_selElectrons,
                                                    &m_selMuons    );
  m_selMET = etmissHandler()->applySelection(m_allMET);

  xAOD::JetContainer allPFlowJets = jetHandlerPFlow()->getCorrectedContainer();
  xAOD::JetContainer selPFlowJets = jetHandlerPFlow()->applySelection(allPFlowJets);

  // Save JVT weight (needs special overlap removal)
  m_jvtJets = jetHandler()->applySelectionNoJvt(m_allJets);
  xAOD::ElectronContainer jvtElecs = m_selElectrons;
  xAOD::MuonContainer jvtMuons = m_selMuons;
  overlapHandler()->removeOverlap(m_selPhotons, m_jvtJets, jvtElecs, jvtMuons);

  // Special overlap removal for PFlow jets
  xAOD::ElectronContainer pflowElecs = m_selElectrons;
  xAOD::MuonContainer pflowMuons = m_selMuons;
  overlapHandler()->removeOverlap(m_selPhotons, selPFlowJets, pflowElecs, pflowMuons);
  
  
  // Adds event weights and catgory to TStore
  // Also sets pointer to photon container, etc., which is used by var's
  setSelectedObjects(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets, &m_selMET, &m_jvtJets);

  // add lepton SF and trigger SF weight to total weight
  if (HG::isMC()) {
    if ( m_selElectrons.size()>=2||m_selMuons.size()>=2) var::weightTrigSF.setValue(eventHandler()->triggerScaleFactor(&m_selElectrons,&m_selMuons));
    double myweight = var::weightSF();
    static SG::AuxElement::Accessor<float> scaleFactor("scaleFactor");
    if (m_selMuons.size()>=2) myweight *= scaleFactor(*m_selMuons[0])*scaleFactor(*m_selMuons[1]);
    else if (m_selElectrons.size()>=2)  myweight *= scaleFactor(*m_selElectrons[0])*scaleFactor(*m_selElectrons[1]);
    var::weightSF.setValue(myweight*var::weightTrigSF());
    var::weight.setValue(weightInitial()*var::weightSF());
  }
  if (photonHandler()->getPointingTool()->updatePointingAuxdata(m_selPhotons).isFailure())
  { HG::fatal("Couldn't retrieve PrimaryVertices, exiting!"); }
  SG::AuxElement::Accessor<float> cluster_time("cluster_time");
  for (auto photon: m_selPhotons) {
    cluster_time(*photon) = photon->caloCluster()->time();
  }

  if (not m_photonAllSys) {
    // Must come before writing out containers (Detailed mode decorates objects for studying)
    // Decorate MET information to HGamEventInfo
    metCatTool()->saveCategoryAndWeight(m_selPhotons, m_selElectrons, m_selMuons, m_selJets, m_selMET);
  }

  // Adds event-level variables to TStore
  if (m_photonAllSys)
    writePhotonAllSys(isSys);
  else
    writeNominalAndSystematic(isSys);

  if (not isSys && not m_photonAllSys) {
    writeNominalOnly();

    if (m_saveDetailed)
      writeDetailed();

    if (m_saveObjects) {
      CP_CHECK("execute()", photonHandler  ()->writeContainer(m_selPhotons  ));
      CP_CHECK("execute()", electronHandler()->writeContainer(m_selElectrons));
      CP_CHECK("execute()", jetHandler     ()->writeContainer(m_selJets     ));
      CP_CHECK("execute()", jetHandlerPFlow()->writeContainer(selPFlowJets));
      CP_CHECK("execute()", muonHandler    ()->writeContainer(m_selMuons    ));
      CP_CHECK("execute()", etmissHandler  ()->writeContainer(m_selMET      ));
    }
  }

  for (auto photon: m_selPhotons) {
    if (photon->pt() != photon->pt()) {
      HG::fatal("Photons have NaN pT?");
    }
  }

  // Adds all event variables (weight, category, isPassed, and pT_yy etc.)
  // to the TEvent output stream
  HG::VarHandler::getInstance()->write();

  return EL::StatusCode::SUCCESS;
}

void ZyCutflowAndMxAOD::writePhotonAllSys(bool isSys)
{
  // Basic event selection flags
  var::isPassedBasic.setValue(m_goodFakeComb ? true : eventHandler()->pass());
  var::isPassed.setValue(m_goodFakeComb ? true : eventHandler()->pass() && pass(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets));
  var::cutFlow.setValue(m_cutFlow);

  if (!isSys) {
    int Nloose = m_preSelPhotons.size();
    eventHandler()->storeVar<char>("isPassedPreselection",Nloose>=2);
  }

  // Add MC only variables
  if (HG::isMC()) {
    if (config()->isDefined(TString::Format("CrossSection.%d", eventInfo()->mcChannelNumber()))) {
      double xs = getCrossSection(), kf = 1.0, ge = 1.0;
      if (config()->isDefined(TString::Format("kFactor.%d", eventInfo()->mcChannelNumber())))
        kf = getKFactor();
      if (config()->isDefined(TString::Format("GeneratorEfficiency.%d", eventInfo()->mcChannelNumber())))
        ge = getGeneratorEfficiency();
      eventHandler()->storeVar<float>("crossSectionBRfilterEff", xs*kf*ge);
    } else {
      eventHandler()->storeVar<float>("crossSectionBRfilterEff", -1);
    }
  }

  writePhotonAllSysVars();
}

void ZyCutflowAndMxAOD::writePhotonAllSysVars(bool /*truth*/)
{

}

void ZyCutflowAndMxAOD::writeNominalAndSystematic(bool isSys)
{
  // Basic event selection flags
   var::isPassedBasic.setValue(m_goodFakeComb ? true : eventHandler()->pass());
  // var::isPassed.setValue(m_goodFakeComb ? true : var::isPassedBasic() && pass(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets));
  var::cutFlow.setValue(m_cutFlow);

  passJetEventCleaning();

  // Basic event weights
  eventHandler()->pileupWeight();
  if (HG::isMC()) eventHandler()->vertexWeight();

  if (!isSys) {
    // Make sure every trigger is checked, and decorated to EventInfo
    eventHandler()->getPassedTriggers();
  }

  // Additional variables useful for non-framework analysis
  int Nloose = m_preSelPhotons.size();
  eventHandler()->storeVar<int>("NLoosePhotons",Nloose);
  eventHandler()->storeVar<float>("met_hardVertexTST", m_selMET["hardVertexTST"] ? m_selMET["hardVertexTST"]->met() : m_selMET["TST"]->met());

  //store passAll flag
  bool passPID = false, passIso = false, passPre = false, passAll = false, passVBSPre = false;
  if (m_selPhotons.size()>0){
    xAOD::Photon *y1 = m_selPhotons[0];
    passPID = photonHandler()->passPIDCut(y1);
    passIso = photonHandler()->passIsoCut(y1, HG::Iso::FixedCutLoose);
  }
  passPre = var::cutFlow()>15; //all except ID and iso
  passAll = passPre && passPID && passIso;
  passVBSPre = passPre && var::N_j()>=2 && var::m_jj_50()>150 * HG::GeV && var::Dy_j_j()>1.0;

  eventHandler()->storeVar<char>("isPassedZyPreSel", passPre);
  eventHandler()->storeVar<char>("isPassedZy", passAll);
  eventHandler()->storeVar<char>("isPassedZyVBSPreSel", passVBSPre);

  if(m_saveAllZ){
    eventHandler()->storeVar<char>("isPassedZ",m_isZ);
  }

  var::isPassed.setValue(passAll);

  
  writeNominalAndSystematicVars();
}

void ZyCutflowAndMxAOD::writeNominalAndSystematicVars(bool truth)
{
  var::m_yy.addToStore(truth);
  var::pT_y1.addToStore(truth);
  var::m_lly.addToStore(truth);
  var::m_ll.addToStore(truth);
  var::m_l1y.addToStore(truth);
  var::m_l2y.addToStore(truth);
  var::pt_lly.addToStore(truth);
  var::pt_ll.addToStore(truth);
  var::pt_llyy.addToStore(truth);
  var::m_llyy.addToStore(truth);
  var::deltaPhi_ll_y.addToStore(truth);
  var::eta_y1.addToStore(truth);
  var::m_emu.addToStore(truth);
  var::m_emuy.addToStore(truth);
  var::N_mu.addToStore(truth);
  var::N_e.addToStore(truth);
  var::N_j.addToStore(truth);
  var::N_j_central.addToStore(truth);
  var::N_j_btag30.addToStore(truth);
  var::N_j_btag.addToStore(truth);
  var::m_jj_50.addToStore(truth);
  var::Dy_j_j_50.addToStore(truth);
  var::Zy_centrality.addToStore(truth);
  var::m_jj.addToStore(truth);
  var::DRmin_y_j.addToStore(truth);
  var::Dy_j_j.addToStore(truth);
  var::Dphi_lly_jj.addToStore(truth);
  var::DR_Zy_jj.addToStore(truth);
  var::deltaR_ll.addToStore(truth);
  var::m_lly2.addToStore(truth);
  if (!truth) {
    var::met_TST  .addToStore(truth);
    var::sumet_TST.addToStore(truth);
  }
}


void ZyCutflowAndMxAOD::writeNominalOnly()
{
  eventHandler()->mu();
  eventHandler()->runNumber();
  eventHandler()->storeVar<unsigned long long>("eventNumber", eventInfo()->eventNumber());
  eventHandler()->centralEventShapeDensity();
  eventHandler()->forwardEventShapeDensity();
  if (HG::isMC()) truthHandler()->centralEventShapeDensity();
  if (HG::isMC()) truthHandler()->forwardEventShapeDensity();
 
  // Additional cut flow granularity
  int Nloose = m_preSelPhotons.size();
  bool passTrigMatch = passTriggerMatch(&m_preSelPhotons);
  bool passIso = false, passPID = false;
  if (Nloose>=2) {
    xAOD::Photon* y1 = m_preSelPhotons[0], *y2 = m_preSelPhotons[1];
    passIso = m_goodFakeComb ? true : photonHandler()->passIsoCut(y1) && photonHandler()->passIsoCut(y2);
    passPID = m_goodFakeComb ? true : photonHandler()->passPIDCut(y1) && photonHandler()->passPIDCut(y2);
  }
  eventHandler()->storeVar<char>("isPassedPreselection",Nloose>=2);
  eventHandler()->storeVar<char>("isPassedTriggerMatch",passTrigMatch);
  eventHandler()->storeVar<char>("isPassedPID",passPID);
  eventHandler()->storeVar<char>("isPassedIsolation",passIso);
  eventHandler()->storeVar<char>("isPassedMassCut",passMyyWindowCut(m_preSelPhotons));

  // Vertex information
  eventHandler()->numberOfPrimaryVertices();
  eventHandler()->selectedVertexZ();
  eventHandler()->hardestVertexZ();
  eventHandler()->pileupVertexSumPt2(); // also sets pileupVertexZ internally

  if (HG::isMC()) truthHandler()->vertexZ();

  const xAOD::VertexContainer* vertices = nullptr;
  if (event()->contains<xAOD::VertexContainer>("PrimaryVertices")) {
    if (event()->retrieve(vertices,"PrimaryVertices").isFailure())
      HG::fatal("Error retrieving PrimaryVertices, exiting");

    //std::vector<float> verticesZ;
    //std::vector<float> verticesScore;
    //static SG::AuxElement::ConstAccessor<float> vertexScore("vertexScore");

    //if (vertices->size() == 1) {
    //  // DxAODs in p2669 have issues with only dummy vertex and vertex score decoration
    //  verticesZ.push_back(-999);
    //  verticesScore.push_back(-99);
    //} else {
    //  for (auto vertex : *vertices){
    //    verticesZ.push_back(vertex->z());
    //    verticesScore.push_back(vertexScore(*vertex));
    //  }
    //}

    //eventHandler()->storeVar<std::vector<float> >("PrimaryVerticesZ"    , verticesZ    ); 
    //eventHandler()->storeVar<std::vector<float> >("PrimaryVerticesScore", verticesScore); 
  }

  // Bunch train information
  eventHandler()->bunchDistanceFromFront();
  eventHandler()->bunchGapBeforeTrain();

  // Add MC only variables
  if (HG::isMC()) {
    truthHandler()->catCoup();

    if (config()->isDefined(TString::Format("CrossSection.%d", eventInfo()->mcChannelNumber()))) {
      double xs = getCrossSection(), kf = 1.0, ge = 1.0;
      if (config()->isDefined(TString::Format("kFactor.%d", eventInfo()->mcChannelNumber())))
        kf = getKFactor();
      if (config()->isDefined(TString::Format("GeneratorEfficiency.%d", eventInfo()->mcChannelNumber())))
        ge = getGeneratorEfficiency();
      eventHandler()->storeVar<float>("crossSectionBRfilterEff", xs*kf*ge);
    } else {
      eventHandler()->storeVar<float>("crossSectionBRfilterEff", -1);
    }
  }

  writeNominalOnlyVars();

}

void ZyCutflowAndMxAOD::writeNominalOnlyVars(bool /*truth*/)
{

}

void ZyCutflowAndMxAOD::writeDetailed()
{
  writeDetailedVars();
}

void ZyCutflowAndMxAOD::writeDetailedVars(bool /*truth*/)
{

}

EL::StatusCode  ZyCutflowAndMxAOD::doTruth()
{
  // Truth particles
  xAOD::TruthParticleContainer all_photons   = truthHandler()->getPhotons();
  xAOD::TruthParticleContainer all_electrons = truthHandler()->getElectrons();
  xAOD::TruthParticleContainer all_muons     = truthHandler()->getMuons();
  xAOD::JetContainer           all_jets      = truthHandler()->getJets();
  xAOD::MissingETContainer     all_met       = truthHandler()->getMissingET();
  xAOD::TruthParticleContainer all_higgs     = truthHandler()->getHiggsBosons();

  // Apply fiducial selections to all containers
  xAOD::TruthParticleContainer photons   = truthHandler()->applyPhotonSelection   (all_photons);
  xAOD::TruthParticleContainer electrons = truthHandler()->applyElectronSelection (all_electrons);
  xAOD::TruthParticleContainer muons     = truthHandler()->applyMuonSelection     (all_muons);
  xAOD::JetContainer           jets      = truthHandler()->applyJetSelection      (all_jets);
  xAOD::JetContainer           bjets     = truthHandler()->applyBJetSelection     (jets);
  xAOD::MissingETContainer     met       = truthHandler()->applyMissingETSelection(all_met);

  // remove truth jets that are from electrons or photons
  truthHandler()->removeOverlap(photons, jets, electrons, muons);

  const xAOD::TruthParticleContainer *truthMuons = nullptr;
  if (event()->retrieve(truthMuons, "TruthMuons").isFailure()) {
    HG::fatal("Can't access TruthMuons Container");
  }
  const xAOD::TruthParticleContainer *truthElectrons = nullptr;
  if (event()->retrieve(truthElectrons, "TruthElectrons").isFailure()) {
    HG::fatal("Can't access TruthElectrons Container");
  }
  HG::DecorateLeptonDressing(all_muons, *truthMuons);
  HG::DecorateLeptonDressing(all_electrons, *truthElectrons);

  //lepton particle isolation
  const xAOD::TruthParticleContainer *truthParts = nullptr;
  if (event()->retrieve(truthParts, "TruthParticles").isFailure()) {
    HG::fatal("Can't access TruthParticles Container");
  }
  SG::AuxElement::Accessor<float> etcone20("etcone20");
  SG::AuxElement::Accessor<float> ptcone20("ptcone20");
  static std::vector<int> ignorePdgIds = { 13, 12, 14, 16, 18 }; // mu, nus
  for (auto part : all_electrons) {
    etcone20(*part) = HG::getTruthIsolation(part, truthParts, 0.2, false, ignorePdgIds);
    ptcone20(*part) = HG::getTruthIsolation(part, truthParts, 0.2, true, {}, 1.0 * HG::GeV);
  }
  for (auto part : all_muons) {
    etcone20(*part) = HG::getTruthIsolation(part, truthParts, 0.2, false, ignorePdgIds);
    ptcone20(*part) = HG::getTruthIsolation(part, truthParts, 0.2, true, {}, 1.0 * HG::GeV);
  }

  // Save truth containers, if configured
  if (m_saveTruthObjects) {
    truthHandler()->writePhotons    (all_photons  );
    truthHandler()->writeElectrons  (all_electrons);
    truthHandler()->writeMuons      (all_muons    );
    truthHandler()->writeJets       (all_jets     );
    truthHandler()->writeMissingET  (met          );
    truthHandler()->writeHiggsBosons(all_higgs    );
    truthHandler()->writeTruthEvents(             );

    addTruthLinks(m_photonContainerName.Data(), m_photonTruthContainerName.Data());
    addTruthLinks(m_elecContainerName.Data()  , m_elecTruthContainerName.Data());
  }

  HG::VarHandler::getInstance()->setTruthContainers(&all_photons, &electrons, &muons, &jets);
  HG::VarHandler::getInstance()->setHiggsBosons(&all_higgs);

  // Adds event-level variables to TStore (this time using truth containers)
  bool truth = true;
  if (m_saveTruthVars) {

    if (m_photonAllSys) {
      writePhotonAllSysVars(truth);
      HG::VarHandler::getInstance()->writeTruth();
      return EL::StatusCode::SUCCESS;
    }

    writeNominalAndSystematicVars(truth);
    writeNominalOnlyVars(truth);
    if (m_saveDetailed)
      writeDetailedVars(truth);

    //var::pT_h1.addToStore(truth);
    //var::y_h1.addToStore(truth);
    //var::m_h1.addToStore(truth);

    // High mass fiducial variables
    static SG::AuxElement::Accessor<float> etcone40("etcone40");

  }

  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}



void ZyCutflowAndMxAOD::fillCutFlow(CutEnum cut, double w) {
  getCutFlowHisto()->Fill(cut);
  if (HG::isData()) return;
  getCutFlowWeightedHisto()->Fill(cut,w);
  return;
}


EL::StatusCode ZyCutflowAndMxAOD::finalize() {
  printf("\nEvent selection cut flow:\n");
  printCutFlowHistos();

  // Write the output to file
  HgammaAnalysis::finalize();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZyCutflowAndMxAOD::fileExecute() {
  // Things that you need to process for each individual input file, even those containing no events.

  HgammaAnalysis::fileExecute();

  // Signals to the code (in execute) that we have to book cutflow information.
  m_newFileMetaData = true;

  return EL::StatusCode::SUCCESS;
}
