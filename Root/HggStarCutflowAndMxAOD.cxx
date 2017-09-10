#include "HGamGamStar/HggStarCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

#include "PhotonVertexSelection/PhotonPointingTool.h"

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
  m_saveTruthObjects = isMC() && config()->getBool("SaveTruthObjects",false);
  m_saveTruthVars    = isMC() && config()->getBool("SaveTruthVariables",false);

  //Save fake photon combinations
  m_enableFakePhotons = isMC() && config()->getBool("SaveFakePhotonCombinations", false);

  // Whether we are running with yybb-tool in detailed mode.
  m_detailedHHyybb = config()->getBool("HHyybb.DetailedInfo",false);

  // Temporary hack for large PhotonAllSys samples
  m_photonAllSys = config()->getStr("PhotonHandler.Calibration.decorrelationModel") == "FULL_v1";

  // a. Event variables
  StrV ignore = {};
  if (isData()) ignore = {".mcChannelNumber", ".mcEventWeights", ".RandomRunNumber", ".truthCentralEventShapeDensity", ".truthForwardEventShapeDensity"};

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

  if (isData()) ignore = {".isEMTight_nofudge", ".isTight_nofudge", ".topoetcone20_DDcorrected", ".topoetcone40_DDcorrected", ".truthOrigin", ".truthType", ".truthConvRadius", ".scaleFactor", ".truthLink", ".parentPdgId", ".pdgId"};
  declareOutputVariables(m_photonContainerName, "MxAOD.Variables.Photon"  , {}, ignore);
  declareOutputVariables("HGamPhotonsWithFakes","MxAOD.Variables.Photon"  , {}, ignore);
  if (isData()) ignore = {".SF_"+mv2_tagger+"_FixedCutBEff_60", ".SF_"+mv2_tagger+"_FixedCutBEff_70", ".SF_"+mv2_tagger+"_FixedCutBEff_77", ".SF_"+mv2_tagger+"_FixedCutBEff_85", ".Eff_"+mv2_tagger+"_FixedCutBEff_60", ".Eff_"+mv2_tagger+"_FixedCutBEff_70", ".Eff_"+mv2_tagger+"_FixedCutBEff_77", ".Eff_"+mv2_tagger+"_FixedCutBEff_85", ".InEff_"+mv2_tagger+"_FixedCutBEff_60", ".InEff_"+mv2_tagger+"_FixedCutBEff_70", ".InEff_"+mv2_tagger+"_FixedCutBEff_77", ".InEff_"+mv2_tagger+"_FixedCutBEff_85", ".HadronConeExclTruthLabelID"};
  declareOutputVariables(m_jetContainerName   , "MxAOD.Variables.Jet"     , {}, ignore);
  if (isData()) ignore = {".scaleFactor", ".truthLink"};
  declareOutputVariables(m_elecContainerName  , "MxAOD.Variables.Electron", {}, ignore);
  if (isData()) ignore = {".scaleFactor"};
  declareOutputVariables(m_muonContainerName  , "MxAOD.Variables.Muon"    , {}, ignore);
  declareOutputVariables("HGamMuonsInJets"    , "MxAOD.Variables.Muon"    , {}, ignore);

  // c. Truth objects
  if (isMC()) {
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
  if (m_newFile && isDAOD() ) {
    addBookKeeping(getCutFlowHisto(),m_N_xAOD,m_N_DxAOD);
    if (isMC()) {
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
  m_isDalitz = isMC() && eventHandler()->isDalitz();

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

  if ( isMC() && (m_saveTruthObjects || m_saveTruthVars))
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
  if(isMC() && m_enableFakePhotons){
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
  if ( requireGRL && isData() && !eventHandler()->passGRL(eventInfo()) ) return GRL;

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
  if (isMC()) var::isDalitzEvent.setValue(m_isDalitz);

  if (!isSys) {
    int Nloose = m_preSelPhotons.size();
    eventHandler()->storeVar<char>("isPassedPreselection",Nloose>=2);
  }

  // Cross-section fiducial regions
  xAOD::JetContainer jets30(SG::VIEW_ELEMENTS);
  double weightBJet30 = 1.0;
  static SG::AuxElement::ConstAccessor<float> SF_bjet("SF_MV2c10_FixedCutBEff_70");
  for (auto jet: m_selJets) {
    if (jet->pt() < 30.0*HG::GeV) continue;
    jets30.push_back(jet);
    weightBJet30 *= SF_bjet(*jet);
  }

  xAOD::JetContainer bjets30 = jetHandler()->applyBJetSelection(jets30);
  int nbjet30 = bjets30.size();

  eventHandler()->storeVar<int>("N_j_btag30", nbjet30);

  char xsec_ttHsemi = var::N_lep_15() >= 1 && var::N_j_30() >= 3 && nbjet30 >= 1;
  char xsec_ttHhad = var::N_lep_15() == 0 && var::N_j_30() >= 4 && nbjet30 >= 1;
  eventHandler()->storeVar<char>("catXS_ttH", xsec_ttHsemi || xsec_ttHhad);
  eventHandler()->storeVar<float>("weightCatXS_ttH", weightBJet30*var::weightN_lep_15());

  // Mass measurement variables
  xAOD::EventInfo *ei = HG::VarHandler::getInstance()->getEventInfoFromStore();
  static SG::AuxElement::Accessor<float> eta_y1("eta_y1"), eta_y2("eta_y2");
  static SG::AuxElement::Accessor<float> etas2_y1("etas2_y1"), etas2_y2("etas2_y2");
  static SG::AuxElement::Accessor<int> conversionType_y1("conversionType_y1"), conversionType_y2("conversionType_y2");

  if (m_selPhotons.size() > 0) {
    eta_y1(*ei) = m_selPhotons[0]->eta();
    etas2_y1(*ei) = m_selPhotons[0]->caloCluster()->etaBE(2);
    conversionType_y1(*ei) = m_selPhotons[0]->conversionType();
  } else {
    eta_y1(*ei) = -99;
    etas2_y1(*ei) = -99;
    conversionType_y1(*ei) = -99;
  }

  if (m_selPhotons.size() > 1) {
    eta_y2(*ei) = m_selPhotons[1]->eta();
    etas2_y2(*ei) = m_selPhotons[1]->caloCluster()->etaBE(2);
    conversionType_y2(*ei) = m_selPhotons[1]->conversionType();
  } else {
    eta_y2(*ei) = -99;
    etas2_y2(*ei) = -99;
    conversionType_y2(*ei) = -99;
  }

  // Add MC only variables
  if (isMC()) {
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

void HggStarCutflowAndMxAOD::writePhotonAllSysVars(bool truth)
{
  var::m_yy.addToStore(truth);
  var::N_j.addToStore(truth);
  var::pT_y1.addToStore(truth);

  // Differential variables
  var::N_j_30.addToStore(truth);
  var::N_j_50.addToStore(truth);
  var::pT_j1_30.addToStore(truth);
  var::pT_j2_30.addToStore(truth);
  var::pT_j3_30.addToStore(truth);
  var::yAbs_j1_30.addToStore(truth);
  var::yAbs_j2_30.addToStore(truth);
  var::HT_30.addToStore(truth);
  var::HTall_30.addToStore(truth);
  var::m_jj_30.addToStore(truth);
  var::Dy_j_j_30.addToStore(truth);
  var::Dphi_j_j_30.addToStore(truth);
  var::Dphi_j_j_30_signed.addToStore(truth);

  if (!truth) {
    var::weightN_lep_15.addToStore(truth);
    var::met_TST.addToStore(truth);
  }
}

void HggStarCutflowAndMxAOD::writeNominalAndSystematic(bool isSys)
{
  // Basic event selection flags
  var::isPassedBasic.setValue(m_goodFakeComb ? true : eventHandler()->pass());
  var::isPassed.setValue(m_goodFakeComb ? true : var::isPassedBasic() && pass(&m_selPhotons, &m_selElectrons, &m_selMuons, &m_selJets));
  var::cutFlow.setValue(m_cutFlow);
  if (isMC()) var::isDalitzEvent.setValue(m_isDalitz);
  passJetEventCleaning();

  // Basic event weights
  eventHandler()->pileupWeight();
  eventHandler()->vertexWeight();

  // Define the signal jets from the no-jvt cut collection
  // Jets with pT > 30/50 GeV are used for pT_j1_30, N_j_50, etc.
  xAOD::JetContainer jvtJets30(SG::VIEW_ELEMENTS), jvtJets50(SG::VIEW_ELEMENTS);
  for (auto jet: m_jvtJets) {
    if (jet->pt() < 30.0*HG::GeV)
      continue;
    jvtJets30.push_back(jet);

    if (jet->pt() < 50.0*HG::GeV)
      continue;
    jvtJets50.push_back(jet);
  }

  eventHandler()->storeVar<float>("weightJvt", HG::JetHandler::multiplyJvtWeights(&m_jvtJets));
  eventHandler()->storeVar<float>("weightJvt_30", HG::JetHandler::multiplyJvtWeights(&jvtJets30));
  eventHandler()->storeVar<float>("weightJvt_50", HG::JetHandler::multiplyJvtWeights(&jvtJets50));

  // Default b-jet information for people outside the framework
  xAOD::JetContainer bjets = jetHandler()->applyBJetSelection(m_selJets);
  eventHandler()->storeVar<int>("N_j_btag", bjets.size());

  xAOD::JetContainer jets30(SG::VIEW_ELEMENTS);
  double weightBJet30 = 1.0;
  static SG::AuxElement::ConstAccessor<float> SF_bjet("SF_MV2c10_FixedCutBEff_70");
  for (auto jet: m_selJets) {
    if (jet->pt() < 30.0*HG::GeV) continue;
    jets30.push_back(jet);
    weightBJet30 *= SF_bjet(*jet);
  }

  xAOD::JetContainer bjets30 = jetHandler()->applyBJetSelection(jets30);
  int nbjet30 = bjets30.size();

  eventHandler()->storeVar<int>("N_j_btag30", nbjet30);

  char xsec_ttHsemi = var::N_lep_15() >= 1 && var::N_j_30() >= 3 && nbjet30 >= 1;
  char xsec_ttHhad = var::N_lep_15() == 0 && var::N_j_30() >= 4 && nbjet30 >= 1;
  eventHandler()->storeVar<char>("catXS_ttH", xsec_ttHsemi || xsec_ttHhad);
  eventHandler()->storeVar<float>("weightCatXS_ttH", weightBJet30*var::weightN_lep_15());

  //Fake photon weight to be included into final weight. If we have it enabled.
  if (m_goodFakeComb)
    var::weight.setValue(weight()*eventHandler()->getVar<float>("weightFakePhotons"));

  if (!isSys) {
    // Make sure every trigger is checked, and decorated to EventInfo
    eventHandler()->getPassedTriggers();
  }

  // Additional variables useful for non-framework analysis
  int Nloose = m_preSelPhotons.size();
  eventHandler()->storeVar<float>("m_yy_resolution",diphotonMassResolution(m_selPhotons));
  eventHandler()->storeVar<int>("NLoosePhotons",Nloose);
  eventHandler()->storeVar<float>("met_hardVertexTST"  , m_selMET["hardVertexTST"] ? m_selMET["hardVertexTST"]->met  () : m_selMET["TST"]->met  ());
  eventHandler()->storeVar<float>("sumet_hardVertexTST", m_selMET["hardVertexTST"] ? m_selMET["hardVertexTST"]->sumet() : m_selMET["TST"]->sumet());
  eventHandler()->storeVar<float>("phi_hardVertexTST"  , m_selMET["hardVertexTST"] ? m_selMET["hardVertexTST"]->phi  () : m_selMET["TST"]->phi  ());

  // High mass variables
  bool passPID = false, passIsoMyy = false, passRelMyy = false, passIsoExot = false, passPtExot = false;
  if (Nloose>=2) {
    xAOD::Photon* y1 = m_preSelPhotons[0], *y2 = m_preSelPhotons[1];

    passPID = m_goodFakeComb ? true : photonHandler()->passPIDCut(y1) && photonHandler()->passPIDCut(y2);
    passIsoMyy = photonHandler()->passIsoCut(y1, HG::Iso::FixedCutTight) &&
                 photonHandler()->passIsoCut(y2, HG::Iso::FixedCutTight);
    passRelMyy = y1->pt()/var::m_yy() >= 0.4 && y2->pt()/var::m_yy() >= 0.3;

    passIsoExot = photonHandler()->passIsoCut(y1, HG::Iso::FixedCutLooseCaloOnly) &&
                  photonHandler()->passIsoCut(y2, HG::Iso::FixedCutLooseCaloOnly);
    passPtExot  = y1->pt() >= 55.0*HG::GeV && y2->pt() >= 55.0*HG::GeV;
  }

  bool passTrigMatch = passTriggerMatch(&m_preSelPhotons);

  // Spin-0 selection
  eventHandler()->storeVar<char>("isPassedIsolationLowHighMyy", passIsoMyy);
  eventHandler()->storeVar<char>("isPassedRelPtCutsLowHighMyy", passRelMyy);
  eventHandler()->storeVar<char>("isPassedLowHighMyy"         , var::isPassedBasic() && Nloose >= 2 && passTrigMatch && passPID && passIsoMyy && passRelMyy);

  // Spin-2 selection
  eventHandler()->storeVar<char>("isPassedIsolationExotic", passIsoExot);
  eventHandler()->storeVar<char>("isPassedlPtCutsExotic"  , passPtExot);
  eventHandler()->storeVar<char>("isPassedExotic"         , var::isPassedBasic() && Nloose >= 2 && passTrigMatch && passPID && passIsoExot && passPtExot);
  eventHandler()->storeVar<char>("isPassedExoticTight"    , var::isPassedBasic() && Nloose >= 2 && passTrigMatch && passPID && passIsoMyy && passPtExot);

  writeNominalAndSystematicVars();
}

void HggStarCutflowAndMxAOD::writeNominalAndSystematicVars(bool truth)
{
  var::m_yy.addToStore(truth);

  // Differential variables
  var::N_j_30.addToStore(truth);
  var::N_j_50.addToStore(truth);
  var::pT_j1_30.addToStore(truth);
  var::pT_j2_30.addToStore(truth);
  var::pT_j3_30.addToStore(truth);
  var::yAbs_j1_30.addToStore(truth);
  var::yAbs_j2_30.addToStore(truth);
  var::HT_30.addToStore(truth);
  var::HTall_30.addToStore(truth);
  var::m_jj_30.addToStore(truth);
  var::Dy_j_j_30.addToStore(truth);
  var::Dphi_j_j_30.addToStore(truth);
  var::Dphi_j_j_30_signed.addToStore(truth);

  var::pT_hard.addToStore(truth);
  var::N_mu   .addToStore(truth);
  var::N_e    .addToStore(truth);
  if (!truth) {
    var::weightN_lep_15.addToStore(truth);
    var::met_TST  .addToStore(truth);
    var::sumet_TST.addToStore(truth);
    var::phi_TST  .addToStore(truth);
  }
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
  eventHandler()->storeVar<char>("isPassedRelPtCuts",passRelativePtCuts(m_preSelPhotons));
  eventHandler()->storeVar<char>("isPassedMassCut",passMyyWindowCut(m_preSelPhotons));

  // Masses with different PV definitions
  const CP::PhotonPointingTool *pointingTool = photonHandler()->getPointingTool();
  if (pointingTool) {
    if (m_cutFlow > AMBIGUITY) {
      xAOD::PhotonContainer leadPhotons = m_preSelPhotons;
      leadPhotons.resize(2);

      // Store m_yy using hardest vertex
      eventHandler()->storeVar<float>("m_yy_hardestVertex",  pointingTool->getCorrectedMass(leadPhotons, eventHandler()->hardestVertexZ()));
      if (isMC())
        eventHandler()->storeVar<float>("m_yy_truthVertex", pointingTool->getCorrectedMass(leadPhotons, truthHandler()->truthVertexZ()) );

      // Store m_yy using zCommon
      eventHandler()->storeVar<float>("m_yy_zCommon",  pointingTool->getCorrectedMass(leadPhotons, xAOD::PVHelpers::getZCommonAndError(eventInfo(), &leadPhotons).first ) );
      var::zCommon.setValue(xAOD::PVHelpers::getZCommonAndError(eventInfo(), &leadPhotons).first);
    } else {
      eventHandler()->storeVar<float>("m_yy_hardestVertex", -99);
      if (isMC())
        eventHandler()->storeVar<float>("m_yy_truthVertex", -99);
      eventHandler()->storeVar<float>("m_yy_zCommon", -99);
    }
  }

  // Vertex information
  eventHandler()->numberOfPrimaryVertices();
  eventHandler()->selectedVertexZ();
  eventHandler()->hardestVertexZ();
  eventHandler()->pileupVertexSumPt2(); // also sets pileupVertexZ internally

  if (isMC()) truthHandler()->truthVertexZ();

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
  if (isMC()) {
    truthHandler()->truthCategory();
    truthHandler()->isVyyOverlap();

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

void HggStarCutflowAndMxAOD::writeNominalOnlyVars(bool truth)
{
  // Truth and reco vars
  var::passMeyCut   .addToStore(truth);
  var::pT_y1        .addToStore(truth);
  var::E_y1         .addToStore(truth);
  var::Zepp         .addToStore(truth);

  var::massTrans    .addToStore(truth);
  var::pTlepMET     .addToStore(truth);

  var::N_j          .addToStore(truth);
  var::N_j_central  .addToStore(truth);
  var::N_j_central30.addToStore(truth);
  var::pT_j1        .addToStore(truth);
  var::pT_j2        .addToStore(truth);
  var::pT_jj        .addToStore(truth);
  var::m_jj         .addToStore(truth);
  var::Dy_j_j       .addToStore(truth);
  var::Dphi_j_j     .addToStore(truth);

  var::DRmin_y_j    .addToStore(truth);

  var::N_e          .addToStore(truth);
  var::N_mu         .addToStore(truth);
  var::N_lep        .addToStore(truth);
  var::m_ee         .addToStore(truth);
  var::m_mumu       .addToStore(truth);

  if (not truth) {
    var::Deta_j_j     .addToStore(truth);
    var::DRmin_y_j_2  .addToStore(truth);
    var::m_alljet     .addToStore(truth);
    var::m_alljet_30  .addToStore(truth);
  }

}

void HggStarCutflowAndMxAOD::writeDetailed()
{
  // Just calling these adds the variables to the TStore
  eventHandler()->selectedVertexSumPt2();
  eventHandler()->hardestVertexSumPt2();
#if __RELEASE__ < 2100
  eventHandler()->eventShapeDensity();
#endif
  eventHandler()->centralEventShapeDensity();
  eventHandler()->forwardEventShapeDensity();

  writeDetailedVars();
}

void HggStarCutflowAndMxAOD::writeDetailedVars(bool truth)
{
  var::Dphi_y_y     .addToStore(truth);
  var::yAbs_j1      .addToStore(truth);
  var::yAbs_j2      .addToStore(truth);
  var::pT_yyj       .addToStore(truth);
  var::Dy_yy_jj     .addToStore(truth);
  var::m_yyjj       .addToStore(truth);

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

    truthHandler()->passFiducial(&all_photons); // calls passFiducialKinOnly internally
    truthHandler()->centralEventShapeDensity();
    truthHandler()->forwardEventShapeDensity();

    var::pT_h1.addToStore(truth);
    var::pT_h2.addToStore(truth);
    var::y_h1.addToStore(truth);
    var::y_h2.addToStore(truth);
    var::m_h1.addToStore(truth);
    var::m_h2.addToStore(truth);

    eventHandler()->storeTruthVar<int>("N_j_btag30", bjets.size());

    char xsec_ttHsemi = var::N_lep_15.truth() >= 1 && var::N_j_30.truth() >= 3 && bjets.size() >= 1;
    char xsec_ttHhad = var::N_lep_15.truth() == 0 && var::N_j_30.truth() >= 4 && bjets.size() >= 1;
    eventHandler()->storeTruthVar<char>("catXS_ttH", xsec_ttHsemi || xsec_ttHhad);

    eventHandler()->storeTruthVar<float>("met_NonInt"  , met["NonInt"]->met()); // MET from neutrinos
    eventHandler()->storeTruthVar<float>("sumet_Int"   , met["Int"   ]->sumet()); // SumET from hadrons, etc.
    eventHandler()->storeTruthVar<float>("met_NonHad"  , truthHandler()->getMissingET_NonHad());

    // High mass fiducial variables
    static SG::AuxElement::Accessor<float> etcone40("etcone40");

    bool isFiducialLowHighMyy = false, isFiducialExotic = false;
    if (all_photons.size() > 1) {
      const xAOD::TruthParticle *gam1 = all_photons[0], *gam2 = all_photons[1];

      isFiducialLowHighMyy = isFiducialExotic = true;

      // Eta cut
      if (fabs(gam1->eta()) >= 2.37 || fabs(gam2->eta()) >= 2.37)
        isFiducialLowHighMyy = isFiducialExotic = false;
      // Isolation cut
      if (etcone40(*gam1)/(gam1->pt() + 120e3) >= 0.05 ||
          etcone40(*gam2)/(gam2->pt() + 120e3) >= 0.05 )
        isFiducialLowHighMyy = isFiducialExotic = false;

      // Scalar relative pT cut
      if (gam1->pt()/var::m_yy.truth() < 0.4 || gam2->pt()/var::m_yy.truth() < 0.3)
        isFiducialLowHighMyy = false;

      // Exotic pT cut
      if (gam1->pt() < 55.0*HG::GeV || gam2->pt() < 55.0*HG::GeV)
        isFiducialExotic = false;

    }

    eventHandler()->storeTruthVar<char>("isFiducialLowHighMyy", isFiducialLowHighMyy);
    eventHandler()->storeTruthVar<char>("isFiducialExotic", isFiducialExotic);

  }

  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}



TH1F* HggStarCutflowAndMxAOD::makeCutFlowHisto(int id, TString suffix) {
  int Ncuts = s_cutDescs.size();

  // create meaningful name of the cutflow histo
  TString name(Form("CutFlow_%s%d",isData()?"Run":"MC",std::abs(id)));

  bool hasMCname = isMC() && config()->isDefined(Form("SampleName.%d",std::abs(id)));

  if(hasMCname){
    name = Form("CutFlow_%s",getMCSampleName(std::abs(id)).Data());
  }

  if (isMC()&&!hasMCname&&config()->getStr("SampleName","sample")!="sample")
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
  if (isData()) return;
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
    printf("\n%s %d cut-flow%s\n",isMC()?"MC sample":"Data run",
           std::abs(entry.first),entry.first>0?", all events":", only Dalitz events");
    printCutFlowHisto(entry.second,0);
  }
  for ( auto entry : m_cFlowHistosWeighted ) {
    printf("\n%s %d cut-flow, weighted events%s\n",isMC()?"MC sample":"Data run",
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
