#include "HGamGamStar/HiggsGamGamStarCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>

#include "PhotonVertexSelection/PhotonPointingTool.h"
#include "ZMassConstraint/ConstraintFit.h"

#include "HGamGamStar/HggStarVariables.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HiggsGamGamStarCutflowAndMxAOD)

HiggsGamGamStarCutflowAndMxAOD::HiggsGamGamStarCutflowAndMxAOD(const char *name)
: MxAODTool(name), m_goodFakeComb(false) { }

HiggsGamGamStarCutflowAndMxAOD::~HiggsGamGamStarCutflowAndMxAOD() {}

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::createOutput()
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

  // Whether we are running with yybb-tool in detailed mode.
  m_detailedHHyybb = config()->getBool("HHyybb.DetailedInfo",false);

  // Temporary hack for large PhotonAllSys samples
  m_photonAllSys = config()->getStr("PhotonHandler.Calibration.decorrelationModel") == "FULL_v1";

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

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::execute()
{
  // Needed for all underlaying tools to be working corectly!
  HgammaAnalysis::execute();

  // Handle File Metadata (need to put this here because we need sample ID to define cutflow histo
  if (m_newFileMetaData) {

    // Initialize cutflow histograms by calling them.
    getCutFlowHisto();
    if (HG::isMC()) {
      getCutFlowHisto(true);
      getCutFlowWeightedHisto();
      getCutFlowWeightedHisto(true);
    }

    // Fill the AOD and DAOD entries of the cutflow histograms.
    if (MxAODTool::fillCutFlowWithBookkeeperInfo() == EL::StatusCode::FAILURE)
    {
      return EL::StatusCode::FAILURE;
    }

    m_crossSectionBRfilterEff = -1;

    if (HG::isMC()) {
      // We want the code to fail if the cross section is not defined.
      m_crossSectionBRfilterEff = getCrossSection();

      if (config()->isDefined(Form("kFactor.%d", eventInfo()->mcChannelNumber())))
      { m_crossSectionBRfilterEff *= getKFactor(); }

      if (config()->isDefined(Form("GeneratorEfficiency.%d", eventInfo()->mcChannelNumber())))
      { m_crossSectionBRfilterEff *= getGeneratorEfficiency(); }
    }

    m_newFileMetaData = false;
  }

  m_isNonHyyStarHiggs = false;

  // flag current event as a MC Dalitz event
  // (needed for cut-flow histograms)
  if (HG::isMC()) {
    const xAOD::TruthParticleContainer *truthParticles = nullptr;
    TString truthPartContName = config()->getStr("TruthParticles.ContainerName", "TruthParticle");
    if (event()->retrieve(truthParticles, truthPartContName.Data()).isFailure())
    { HG::fatal("Can't access TruthParticleContainer"); }
    m_isNonHyyStarHiggs = HG::isMC() && HG::eventIsNonHyyStarHiggs(truthParticles);
    var::isNonHyyStarHiggs.setTruthValue(m_isNonHyyStarHiggs);
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
      CP_CHECK("HiggsGamGamStarCutflowAndMxAOD::execute()", applySystematicVariation(sys));
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
HiggsGamGamStarCutflowAndMxAOD::CutEnum HiggsGamGamStarCutflowAndMxAOD::cutflow()
{

  m_preSelPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selElectrons = xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
  m_selMuons = xAOD::MuonContainer(SG::VIEW_ELEMENTS);
  m_allJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);
  m_selJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);

  //Check if there are two good fakes. Needed so we dont slim the event at trigger.
  m_goodFakeComb = false;
  if(HG::isMC() && m_enableFakePhotons){
      double weightFakePhotons = 1;
      xAOD::PhotonContainer photonsWithFakes = getFakePhotons(weightFakePhotons);
      m_goodFakeComb = photonsWithFakes.size()>1 ? true : false;
      if(m_goodFakeComb) return PASSALL;
  }

  //==== CUT "-1" : Remove truth Higgs events that are not leptonic Dalitz ====
  if ( m_isNonHyyStarHiggs ) return HIGGS_LEP_DALITZ;

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

  //==== CUT 4 : 2 SF leptons (before OR) ====
  if (m_preSelElectrons.size() < 2 && dirtyMuons.size() < 2) return TWO_SF_LEPTONS;

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

  //==== CUT 5 : Require one loose photon, pT > "PtPreCutGeV" GeV ====
  if (nloose<1) return ONE_LOOSE_GAM;

  //==== CUT 6 : Preselection
  // - Require two loose photons that also pass e-gamma ambiguity ====
  static bool requireAmbiguity = config()->getBool("PhotonHandler.Selection.ApplyAmbiguityCut", false);
  if (requireAmbiguity && namb<1) return AMBIGUITY;

  // static bool requireHV = config()->getBool("PhotonHandler.Selection.ApplyHVCut", false);
  // if (requireHV && nHV<1) return AMBIGUITY;

  m_allJets = jetHandler()->getCorrectedContainer();
  m_selJets = jetHandler()->applySelection(m_allJets);

  // Removes overlap with candidate photon, and any additional tight photons (if option set)
  overlapHandler()->removeOverlap(m_selPhotons, m_selJets, m_preSelElectrons, dirtyMuons);

  //above doesn't have option to remove photon overlapping with lepton
  overlapHandler()->removeOverlap(m_selPhotons, m_preSelElectrons, 0.4);
  overlapHandler()->removeOverlap(m_selPhotons, dirtyMuons, 0.4);

  // Muon cleaning should be done after overlap removal
  m_preSelMuons = muonHandler()->applyCleaningSelection(dirtyMuons);

  // Select Z candidate after overlap removal.
  // Find the highest-mll pair closest to (13 TeV)

  int sel_muon1 = -1, sel_muon2 = -1;
  double return_mmumu = -1;
  HG::AssignZbosonIndices(m_preSelMuons,sel_muon1,sel_muon2,return_mmumu,13000.*HG::GeV);

  int sel_ele1 = -1, sel_ele2 = -1;
  double return_mee = -1;
  HG::AssignZbosonIndices(m_preSelElectrons,sel_ele1,sel_ele2,return_mee,13000.*HG::GeV);

  if (return_mmumu < 0 && return_mee < 0) return TWO_SF_LEPTONS_POSTOR;

  if (return_mmumu > return_mee) {
    m_selMuons.push_back(m_preSelMuons[sel_muon1]);
    m_selMuons.push_back(m_preSelMuons[sel_muon2]);
  } else {
    m_selElectrons.push_back(m_preSelElectrons[sel_ele1]);
    m_selElectrons.push_back(m_preSelElectrons[sel_ele2]);
  }

  if (m_selPhotons.size()==0) return ONE_PHOTON_POSTOR;

  //trigger matching
  static bool requireTriggerMatch = config()->getBool("EventHandler.CheckTriggerMatching", true);

  if ( requireTriggerMatch ){
    StrV m_requiredTriggers = config()->getStrV("EventHandler.RequiredTriggers");
    int itrigmatch=0;
    for (auto trig: m_requiredTriggers) {
      if (passTriggerMatch(trig, NULL, &m_selElectrons, &m_selMuons, NULL) ) itrigmatch++;
    }
    if (itrigmatch==0) return TRIG_MATCH;
  }


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

EL::StatusCode  HiggsGamGamStarCutflowAndMxAOD::doReco(bool isSys){
  // Do anything you missed in cutflow, and save the objects.

  // Rebuild MET using selected objects
  m_allMET = etmissHandler()->getCorrectedContainer(&m_selPhotons    ,
                                                    &m_allJets     ,
                                                    &m_selElectrons,
                                                    &m_selMuons    );
  m_selMET = etmissHandler()->applySelection(m_allMET);

  // Save JVT weight (needs special overlap removal)
  m_jvtJets = jetHandler()->applySelectionNoJvt(m_allJets);
  xAOD::ElectronContainer jvtElecs = m_selElectrons;
  xAOD::MuonContainer jvtMuons = m_selMuons;
  overlapHandler()->removeOverlap(m_selPhotons, m_jvtJets, jvtElecs, jvtMuons);

  // Adds event weights and catgory to TStore
  // Also sets pointer to photon container, etc., which is used by var's
  setSelectedObjects(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets, &m_selMET, &m_jvtJets);

  if (not m_photonAllSys) {
    // Must come before writing out containers (Detailed mode decorates objects for studying)
    // Decorate MET information to HGamEventInfo
    metCatTool()->saveCategoryAndWeight(m_selPhotons, m_selElectrons, m_selMuons, m_selJets, m_selMET);
  }

  // Adds event-level variables to TStore
  if (m_photonAllSys)
  {
    // Write when using FULL_v1 photon systematics (separated from the other systs)
    writePhotonAllSys(isSys);
  }
  else
  {
    // Write in the nominal and systematics loops
    writeNominalAndSystematic();

    // Write only in the nominal loop
    if (not isSys)
    {
      writeNominalOnly();

      if (m_saveDetailed) { writeDetailed(); }

      if (m_saveObjects) {
        CP_CHECK("execute()", photonHandler  ()->writeContainer(m_selPhotons  ));
        CP_CHECK("execute()", electronHandler()->writeContainer(m_selElectrons));
        CP_CHECK("execute()", jetHandler     ()->writeContainer(m_selJets     ));
        CP_CHECK("execute()", muonHandler    ()->writeContainer(m_selMuons    ));
        CP_CHECK("execute()", etmissHandler  ()->writeContainer(m_selMET      ));
      }
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

void HiggsGamGamStarCutflowAndMxAOD::writePhotonAllSys(bool isSys)
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
    eventHandler()->storeVar<float>("crossSectionBRfilterEff", m_crossSectionBRfilterEff);
  }

  writePhotonAllSysVars();
}

void HiggsGamGamStarCutflowAndMxAOD::writePhotonAllSysVars(bool /*truth*/)
{

}

void HiggsGamGamStarCutflowAndMxAOD::writeNominalAndSystematic()
{
  // Basic event selection flags
  // var::isPassedBasic.setValue(m_goodFakeComb ? true : eventHandler()->pass());
  // var::isPassed.setValue(m_goodFakeComb ? true : var::isPassedBasic() && pass(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets));
  var::cutFlow.setValue(m_cutFlow);
  passJetEventCleaning();

  // Basic event weights
  eventHandler()->pileupWeight();
  eventHandler()->vertexWeight();

  // Additional variables useful for non-framework analysis
  int Nloose = m_preSelPhotons.size();
  eventHandler()->storeVar<int>("NLoosePhotons",Nloose);

  writeNominalAndSystematicVars();
}

void HiggsGamGamStarCutflowAndMxAOD::writeNominalAndSystematicVars(bool truth)
{
  // var::m_yy.addToStore(truth);
  var::m_lly.addToStore(truth);
  var::m_ll.addToStore(truth);
  var::pt_lly.addToStore(truth);
  var::pt_ll.addToStore(truth);

  var::N_mu   .addToStore(truth);
  var::N_e    .addToStore(truth);
}


void HiggsGamGamStarCutflowAndMxAOD::writeNominalOnly()
{
  eventHandler()->mu();
  eventHandler()->runNumber();

  // Make sure every trigger is checked, and decorated to EventInfo
  eventHandler()->getPassedTriggers();

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
    eventHandler()->storeVar<float>("crossSectionBRfilterEff", m_crossSectionBRfilterEff);
  }

  writeNominalOnlyVars();

}

void HiggsGamGamStarCutflowAndMxAOD::writeNominalOnlyVars(bool /*truth*/)
{

}

void HiggsGamGamStarCutflowAndMxAOD::writeDetailed()
{
  writeDetailedVars();
}

void HiggsGamGamStarCutflowAndMxAOD::writeDetailedVars(bool /*truth*/)
{

}

EL::StatusCode  HiggsGamGamStarCutflowAndMxAOD::doTruth()
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

    // Add truth-only variables here.
    var::pT_h1.addToStore(truth);
    var::y_h1.addToStore(truth);
    var::m_h1.addToStore(truth);

    var::pT_l1_h1.addToStore(truth);
    var::pT_l2_h1.addToStore(truth);
    var::deltaR_l1l2_h1.addToStore(truth);
    var::ystar_pdg_flavor.addToStore(truth);
    var::pT_yDirect_h1.addToStore(truth);
    var::m_yStar_undressed_h1.addToStore(truth);

    // High mass fiducial variables
    static SG::AuxElement::Accessor<float> etcone40("etcone40");

  }

  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}




void HiggsGamGamStarCutflowAndMxAOD::fillCutFlow(CutEnum cut, double w) {
  getCutFlowHisto()->Fill(cut);
  if (HG::isData()) return;
  getCutFlowWeightedHisto()->Fill(cut,w);
  if (m_isNonHyyStarHiggs) return;
  getCutFlowHisto(true)->Fill(cut);
  getCutFlowWeightedHisto(true)->Fill(cut,w);
}


EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::finalize() {
  printf("\nEvent selection cut flow:\n");
  printCutFlowHistos();

  // Write the output to file
  HgammaAnalysis::finalize();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::fileExecute() {
  // Things that you need to process for each individual input file, even those containing no events.

  HgammaAnalysis::fileExecute();

  // Signals to the code (in execute) that we have to book cutflow information.
  m_newFileMetaData = true;

  return EL::StatusCode::SUCCESS;
}
