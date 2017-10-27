#include "HGamGamStar/HggStarCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

#include "PhotonVertexSelection/PhotonPointingTool.h"
#include "ZMassConstraint/ConstraintFit.h"

#include "HGamGamStar/HggStarVariables.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HggStarCutflowAndMxAOD)

HggStarCutflowAndMxAOD::HggStarCutflowAndMxAOD(const char *name)
: HgammaAnalysis(name), m_goodFakeComb(false), m_N_xAOD(0), m_N_DxAOD(0),
  m_sumw_xAOD(0.0), m_sumw2_xAOD(0.0), m_sumw_DxAOD(0.0), m_sumw2_DxAOD(0.0) { }

HggStarCutflowAndMxAOD::~HggStarCutflowAndMxAOD() {}

void HggStarCutflowAndMxAOD::declareOutputVariables(TString outName, TString configKey, StrV extra, StrV ignore) {
  if (config()->isDefined(configKey)) {
    TString vars = config()->getStr(configKey).Data();

    if (m_saveDetailed) {
      TString detailKey = configKey.ReplaceAll("Variables", "DetailedVariables");
      if (config()->isDefined(detailKey)) {
        TString detailed = config()->getStr(detailKey);
        vars += "." + detailed;
      }
    }

    for (TString val: extra)
      vars += val;

    for (TString val: ignore)
      vars = vars.ReplaceAll(val, "");

    event()->setAuxItemList((outName+"Aux.").Data(), vars.Data());
  }
  else HG::fatal("Cannot find "+configKey);
}

EL::StatusCode HggStarCutflowAndMxAOD::createOutput()
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

EL::StatusCode HggStarCutflowAndMxAOD::execute()
{
  // Needed for all underlaying tools to be working corectly!
  HgammaAnalysis::execute();

  // if it's a new file, update the book-keeper
  if (m_newFile && HG::isDAOD() ) {
    addBookKeeping(getCutFlowHisto(),m_N_xAOD,m_N_DxAOD);
    if (HG::isMC()) {
      addBookKeeping(getCutFlowHisto(true),m_N_xAOD,m_N_DxAOD);
      addBookKeeping(getCutFlowWeightedHisto(),m_sumw_xAOD,m_sumw_DxAOD,m_sumw2_xAOD,m_sumw2_DxAOD);
      addBookKeeping(getCutFlowWeightedHisto(true),m_sumw_xAOD,m_sumw_DxAOD,m_sumw2_xAOD,m_sumw2_DxAOD);
    }
    m_newFile=false;
  }

  // apply cuts. Returned value will be the last passed cut
  m_cutFlow = cutflow();

  // flag current event as a MC Dalitz event
  // (needed for cut-flow histograms)
  m_isDalitz = HG::isMC() && eventHandler()->isDalitz();

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
      CP_CHECK("HggStarCutflowAndMxAOD::execute()", applySystematicVariation(sys));
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
HggStarCutflowAndMxAOD::CutEnum HggStarCutflowAndMxAOD::cutflow()
{
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
  m_selElectrons = electronHandler()->applySelection(m_allElectrons);

  m_allMuons = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer dirtyMuons = muonHandler()->applySelection(m_allMuons);

  //==== CUT 4 : 2 SF leptons (before OR) ====
  if (m_selElectrons.size() < 2 && dirtyMuons.size() < 2) return TWO_SF_LEPTONS;

  // Get object containers
  m_allPhotons = photonHandler()->getCorrectedContainer();
  m_preSelPhotons = photonHandler()->applyPreSelection(m_allPhotons);
  m_selPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  if (m_preSelPhotons.size()  ) m_selPhotons.push_back(m_preSelPhotons[0]);

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

  //==== CUT 5 : Require two loose photons, pT>25 GeV ====
  if (nloose<1) return ONE_LOOSE_GAM;

  //==== CUT 6 : Preselection
  // - Require two loose photons that also pass e-gamma ambiguity ====
  static bool requireAmbiguity = config()->getBool("PhotonHandler.Selection.ApplyAmbiguityCut", false);
  if (requireAmbiguity && namb<1) return AMBIGUITY;

  // static bool requireHV = config()->getBool("PhotonHandler.Selection.ApplyHVCut", false);
  // if (requireHV && nHV<1) return AMBIGUITY;

  m_allJets = jetHandler()->getCorrectedContainer();
  m_selJets = jetHandler()->applySelection(m_allJets);

  // Removes overlap with candidate diphoton system, and any additional tight photons (if option set)
  overlapHandler()->removeOverlap(m_selPhotons, m_selJets, m_selElectrons, dirtyMuons);

  // Muon cleaning should be done after overlap removal
  m_selMuons = muonHandler()->applyCleaningSelection(dirtyMuons);

  // Select Z candidate after overlap removal.

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

EL::StatusCode  HggStarCutflowAndMxAOD::doReco(bool isSys){
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

void HggStarCutflowAndMxAOD::writePhotonAllSys(bool isSys)
{
  // Basic event selection flags
  var::isPassedBasic.setValue(m_goodFakeComb ? true : eventHandler()->pass());
  var::isPassed.setValue(m_goodFakeComb ? true : eventHandler()->pass() && pass(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets));
  var::cutFlow.setValue(m_cutFlow);
  if (HG::isMC()) var::isDalitzEvent.setValue(m_isDalitz);

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

void HggStarCutflowAndMxAOD::writePhotonAllSysVars(bool /*truth*/)
{

}

void HggStarCutflowAndMxAOD::writeNominalAndSystematic(bool isSys)
{
  // Basic event selection flags
  // var::isPassedBasic.setValue(m_goodFakeComb ? true : eventHandler()->pass());
  // var::isPassed.setValue(m_goodFakeComb ? true : var::isPassedBasic() && pass(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets));
  var::cutFlow.setValue(m_cutFlow);
  if (HG::isMC()) var::isDalitzEvent.setValue(m_isDalitz);
  passJetEventCleaning();

  // Basic event weights
  eventHandler()->pileupWeight();
  eventHandler()->vertexWeight();

  if (!isSys) {
    // Make sure every trigger is checked, and decorated to EventInfo
    eventHandler()->getPassedTriggers();
  }

  // Additional variables useful for non-framework analysis
  int Nloose = m_preSelPhotons.size();
  eventHandler()->storeVar<int>("NLoosePhotons",Nloose);

  writeNominalAndSystematicVars();
}

void HggStarCutflowAndMxAOD::writeNominalAndSystematicVars(bool truth)
{
  // var::m_yy.addToStore(truth);
  var::m_lly.addToStore(truth);
  var::m_ll.addToStore(truth);

  var::N_mu   .addToStore(truth);
  var::N_e    .addToStore(truth);
}


void HggStarCutflowAndMxAOD::writeNominalOnly()
{
  eventHandler()->mu();
  eventHandler()->runNumber();

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

    std::vector<float> verticesZ;
    std::vector<float> verticesScore;
    static SG::AuxElement::ConstAccessor<float> vertexScore("vertexScore");

    if (vertices->size() == 1) {
      // DxAODs in p2669 have issues with only dummy vertex and vertex score decoration
      verticesZ.push_back(-999);
      verticesScore.push_back(-99);
    } else {
      for (auto vertex : *vertices){
        verticesZ.push_back(vertex->z());
        verticesScore.push_back(vertexScore(*vertex));
      }
    }

    eventHandler()->storeVar<std::vector<float> >("PrimaryVerticesZ"    , verticesZ    ); 
    eventHandler()->storeVar<std::vector<float> >("PrimaryVerticesScore", verticesScore); 
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

void HggStarCutflowAndMxAOD::writeNominalOnlyVars(bool /*truth*/)
{

}

void HggStarCutflowAndMxAOD::writeDetailed()
{
  writeDetailedVars();
}

void HggStarCutflowAndMxAOD::writeDetailedVars(bool /*truth*/)
{

}

EL::StatusCode  HggStarCutflowAndMxAOD::doTruth()
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

    var::pT_h1.addToStore(truth);
    var::y_h1.addToStore(truth);
    var::m_h1.addToStore(truth);

    // High mass fiducial variables
    static SG::AuxElement::Accessor<float> etcone40("etcone40");

  }

  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}



TH1F* HggStarCutflowAndMxAOD::makeCutFlowHisto(int id, TString suffix) {
  int Ncuts = s_cutDescs.size();

  // create meaningful name of the cutflow histo
  TString name(Form("CutFlow_%s%d",HG::isData()?"Run":"MC",std::abs(id)));

  bool hasMCname = HG::isMC() && config()->isDefined(Form("SampleName.%d",std::abs(id)));

  if(hasMCname){
    name = Form("CutFlow_%s",getMCSampleName(std::abs(id)).Data());
  }

  if (HG::isMC()&&!hasMCname&&config()->getStr("SampleName","sample")!="sample")
    name="CutFlow_"+config()->getStr("SampleName");
  name+=suffix;

  // maybe should make sure we don't switch directory?
  TDirectory *dir = gDirectory;
  TFile *file = wk()->getOutputFile("MxAOD");
  TH1F *h = new TH1F(name,name,Ncuts,0,Ncuts);
  h->SetDirectory(file); // probably not needed
  for (int bin=1;bin<=Ncuts;++bin)
    h->GetXaxis()->SetBinLabel(bin,s_cutDescs[bin-1]);
  dir->cd();
  return h;
}

void HggStarCutflowAndMxAOD::fillCutFlow(CutEnum cut, double w) {
  getCutFlowHisto()->Fill(cut);
  if (HG::isData()) return;
  getCutFlowWeightedHisto()->Fill(cut,w);
  if (!m_isDalitz) return;
  getCutFlowHisto(true)->Fill(cut);
  getCutFlowWeightedHisto(true)->Fill(cut,w);
}


EL::StatusCode HggStarCutflowAndMxAOD::finalize() {
  printf("\nEvent selection cut flow:\n");
  printCutFlowHistos();

  // Write the output to file
  HgammaAnalysis::finalize();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HggStarCutflowAndMxAOD::fileExecute() {
  HgammaAnalysis::fileExecute();

  // Tell the code a new file has just been opened
  m_newFile=true;

  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (MetaData == nullptr)
    HG::fatal("Couldn't find MetaData TTree in event, is this a proper xAOD file? Exiting.");

  // The isDAOD flag set in HgammaAnalysis::changeInput is not available here, so we get it maually.
  bool tmp_isDAOD = false;
  for (auto branch : * MetaData->GetListOfBranches())
  {
    if (TString(branch->GetName()).Contains("StreamDAOD")) {
      tmp_isDAOD = true;
      break;
    }
  }

  if ( tmp_isDAOD ) {
    // If we get here, this is a DxAOD

    // 1. Check if there if the incomplete book-keeper object has entreies (that would be bad!)
    const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
    EL_CHECK( "fileExecute()", wk()->xaodEvent()->retrieveMetaInput(incompleteCBC,"IncompleteCutBookkeepers") );
    if ( incompleteCBC->size() != 0 ) {
      // fatal or error? let's start with a hard fatal!
      // HG::fatal("Issue with DxAOD book keeper. It's incomplete. File corrupted?");
      Warning("fileExecute()",
	      "Issue with DxAOD book keeper. It's incomplete. File corrupted? %s %s",
	      "If this is data, this is known to happen (but not understood).",
              "If this is MC, this is expected to happen, and can probably be ignored.");
    }

    // 2. now get the actual bookkeeper
    const xAOD::CutBookkeeperContainer* completeCBC = nullptr;
    EL_CHECK( "fileExecute()", wk()->xaodEvent()->retrieveMetaInput(completeCBC,"CutBookkeepers") );

    int maxAcycle = -1, maxDcycle = -1;
    for ( auto cbk : *completeCBC ) {
      Info("fileExecute()","  Book keeper name=%s, inputStream=%s, cycle=%d, nAcceptedEvents=%d", cbk->name().c_str(), cbk->inputStream().c_str(), cbk->cycle(), (int)cbk->nAcceptedEvents());

      if ( cbk->name().empty() ) continue;

      // Use the DxAOD numbers from the largest cycle
      if (TString(cbk->name()).Contains("HIGG1D1") && cbk->inputStream() == "StreamAOD" && cbk->cycle() > maxDcycle) {
        maxDcycle     = cbk->cycle();
	m_N_DxAOD     = cbk->nAcceptedEvents();
	m_sumw_DxAOD  = cbk->sumOfEventWeights();
	m_sumw2_DxAOD = cbk->sumOfEventWeightsSquared();
      }


      // Use the xAOD numbers from the largest cycle
      if (cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD" && cbk->cycle() > maxAcycle) {
        maxAcycle    = cbk->cycle();
	m_N_xAOD     = cbk->nAcceptedEvents();
	m_sumw_xAOD  = cbk->sumOfEventWeights();
	m_sumw2_xAOD = cbk->sumOfEventWeightsSquared();
      }
    }
    Info("fileExecute()","Book keeper anticipates %i events in current input file (%i in parent xAOD)",m_N_DxAOD,m_N_xAOD);
  }
  return EL::StatusCode::SUCCESS;
}

void HggStarCutflowAndMxAOD::printCutFlowHistos() {
  for ( auto entry : m_cFlowHistos ) {
    printf("\n%s %d cut-flow%s\n",HG::isMC()?"MC sample":"Data run",
           std::abs(entry.first),entry.first>0?", all events":", only Dalitz events");
    printCutFlowHisto(entry.second,0);
  }
  for ( auto entry : m_cFlowHistosWeighted ) {
    printf("\n%s %d cut-flow, weighted events%s\n",HG::isMC()?"MC sample":"Data run",
           std::abs(entry.first),entry.first>0?", all events":", only Dalitz events");
    printCutFlowHisto(entry.second,2);
  }
  printf("\n");
}

void HggStarCutflowAndMxAOD::printCutFlowHisto(TH1F *h, int Ndecimals) {
  TString format("  %-24s%10."); format+=Ndecimals; format+="f%11.2f%%%11.2f%%\n";
  int all_bin = h->FindBin(ALLEVTS);
  printf("  %-24s%10s%12s%12s\n","Event selection","Nevents","Cut rej.","Tot. eff.");
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    double ntot=h->GetBinContent(all_bin), n=h->GetBinContent(bin), nprev=h->GetBinContent(bin-1);
    TString cutName(h->GetXaxis()->GetBinLabel(bin));
    cutName.ReplaceAll("#it{m}_{#gamma#gamma}","m_yy");
    if (bin==1||nprev==0||n==nprev)
      printf(format.Data(),cutName.Data(),n,-1e-10,n/ntot*100);
    else // if the cut does something, print more information
      printf(format.Data(),cutName.Data(),n,(n-nprev)/nprev*100,n/ntot*100);
  }
}
