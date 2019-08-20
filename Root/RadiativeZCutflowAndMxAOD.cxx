#include "HGamGamStar/RadiativeZCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>
#include "HGamAnalysisFramework/TruthUtils.h"

#include "HGamGamStar/ExtraHggStarObjects.h"
#include "HGamGamStar/TrackElectronMap.h"
#include "xAODTracking/VertexAuxContainer.h"



// #include "PhotonVertexSelection/PhotonPointingTool.h"
// #include "ZMassConstraint/ConstraintFit.h"

// this is needed to distribute the algorithm to the workers
ClassImp(RadiativeZCutflowAndMxAOD)

RadiativeZCutflowAndMxAOD::RadiativeZCutflowAndMxAOD(const char *name)
: MxAODTool(name)
  , m_trackHandler(nullptr)
{ }

RadiativeZCutflowAndMxAOD::~RadiativeZCutflowAndMxAOD() {}

EL::StatusCode RadiativeZCutflowAndMxAOD::initialize()
{
  // Make sure that all of our tools here are initialized before calling HgammaAnalysis::initialize()
  // in order to correctly fillSystematicsList()
  HgammaAnalysis::initialize();

  HG::ExtraHggStarObjects::getInstance()->setEventAndStore(event(), store());

  for (auto trig: eventHandler()->getRequiredTriggers()) {
    if (!config()->isDefined("EventHandler.TriggerMatchType." + trig))
      HG::fatal("You must define a EventHandler.TriggerMatchType for trigger "+trig);
    if (!config()->isDefined("EventHandler.RunNumbers." + trig))
      HG::fatal("You must define EventHandler.RunNumbers (years of validity) for trigger "+trig);
  }

  m_trackHandler = new HG::TrackHandler("TrackHandler", event(), store());
  ANA_CHECK(m_trackHandler->initialize(*config()));

  m_mergedElectronID = new HG::MergedElectronID();
  ANA_CHECK(m_mergedElectronID->initialize(*config()));

  m_mergedElectronID_v2 = new HG::MergedElectronID_v2();
  ANA_CHECK(m_mergedElectronID_v2->initialize(*config()));

  // Resolved electron ID preselection. Applied in FindZboson_ElectronChannelAware.
  m_eleIDPreselection = config()->getStr("ResolvedElectrons.Preselection.PID","VeryLoose");

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RadiativeZCutflowAndMxAOD::createOutput()
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

  // Whether to save objects (photons, jets ...)
  m_saveObjects = config()->getBool("SaveObjects",false);
  m_skipElectronObjects = false; // to be filled in execute, when the mcChannelNumber is known
  m_skipMuonObjects = false; // to be filled in execute, when the mcChannelNumber is known

  // Whether to save the list of differential variables
  m_saveDetailed = config()->getBool("SaveDetailedVariables",false);

  // Whether to save the truth objects and differential variables
  m_saveTruthObjects = HG::isMC() && config()->getBool("SaveTruthObjects",false);
  m_saveTruthVars    = HG::isMC() && config()->getBool("SaveTruthVariables",false);

  // a. Event variables
  StrV ignore = {};
  if (HG::isData()) ignore = {".mcChannelNumber", ".mcEventWeights", ".RandomRunNumber", ".truthCentralEventShapeDensity", ".truthForwardEventShapeDensity"};

  StrV trigs = config()->getStrV("EventHandler.RequiredTriggers");
  StrV extra = {};

  if (!m_applySystematics) {
    for (auto trig: trigs) extra.push_back(".passTrig_"+trig);
  }

  declareOutputVariables(m_evtInfoName,"MxAOD.Variables.EventInfo", extra, ignore);

  if (config()->isDefined("MxAOD.Variables.HGamEventInfo")) {
    TString HGamEvtInfo_str = "HGam"+m_evtInfoName;
    declareOutputVariables(HGamEvtInfo_str,"MxAOD.Variables.HGamEventInfo", extra, ignore);
  }

  // a.2 TruthEvents variables
  ignore = {};
  extra = {};
  if (m_saveTruthVars) {
    declareOutputVariables(m_truthEvtsName,"MxAOD.Variables.TruthEvents", extra, ignore);
  }

  // b. Selected objects

  if (HG::isData()) ignore = {".isEMTight_nofudge", ".isTight_nofudge", ".topoetcone20_DDcorrected", ".topoetcone40_DDcorrected", ".truthOrigin", ".truthType", ".truthConvRadius", ".scaleFactor", ".truthLink", ".parentPdgId", ".pdgId"};
  declareOutputVariables(m_photonContainerName, "MxAOD.Variables.Photon"  , {}, ignore);
  declareOutputVariables("HGamPhotonsWithFakes","MxAOD.Variables.Photon"  , {}, ignore);

  if (HG::isData()) ignore = {".SF_MV2c10_FixedCutBEff_60", ".SF_MV2c10_FixedCutBEff_70",
                              ".SF_MV2c10_FixedCutBEff_77", ".SF_MV2c10_FixedCutBEff_85",
                              ".Eff_MV2c10_FixedCutBEff_60", ".Eff_MV2c10_FixedCutBEff_70",
                              ".Eff_MV2c10_FixedCutBEff_77", ".Eff_MV2c10_FixedCutBEff_85",
                              ".InEff_MV2c10_FixedCutBEff_60", ".InEff_MV2c10_FixedCutBEff_70",
                              ".InEff_MV2c10_FixedCutBEff_77", ".InEff_MV2c10_FixedCutBEff_85",
                              ".HadronConeExclTruthLabelID"};
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

    if (m_saveTruthObjects) {
      declareOutputVariables(m_photonTruthContainerName  , "MxAOD.Variables.TruthPhotons"    );
      declareOutputVariables(m_elecTruthContainerName    , "MxAOD.Variables.TruthElectrons"  );
      declareOutputVariables(m_muonTruthContainerName    , "MxAOD.Variables.TruthMuons"      );
      declareOutputVariables(m_jetTruthContainerName     , "MxAOD.Variables.TruthJets"       );
      declareOutputVariables("HGam"+config()->getStr("TruthHandler.HiggsBosonContainerName"), "MxAOD.Variables.TruthHiggsBosons");
    }
  }

  if (m_applySystematics) {
    for (auto sys: getSystematics()) {
      // ignore nominal case, already done!
      if (sys.name() == "") continue;

      // Copied from VarHandler.cxx
      TString sysName = sys.name() == "" ? "" : ("_" + sys.name()).c_str();
      sysName.ReplaceAll(" ", "_");

      TString evtInfoNameSys = TString("HGam") + m_evtInfoName + sysName;

      ignore = {};
      extra = {};
      declareOutputVariables(evtInfoNameSys,"MxAOD.Variables.EventInfoSyst", extra, ignore);

    }
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RadiativeZCutflowAndMxAOD::execute()
{
  // Needed for all underlaying tools to be working corectly!
  HgammaAnalysis::execute();

  // Clear containers which point to objects from previous event
  HG::ExtraHggStarObjects::getInstance()->clearContainers();

  // Handle File Metadata (need to put this here because we need sample ID to define cutflow histo
  if (m_newFileMetaData) {

    // Initialize cutflow histograms by calling them.
    getCutFlowHisto(false /*weighted*/);

    // per-channel cutflows (mainly for debugging)
    getCutFlowHisto(false, HG::DIMUON);
    getCutFlowHisto(false, HG::RESOLVED_DIELECTRON);

    if (HG::isMC()) {
      getCutFlowHisto(true);

      getCutFlowHisto(true, HG::DIMUON);
      getCutFlowHisto(true, HG::RESOLVED_DIELECTRON);

    }

    // Fill the AOD and DAOD entries of the cutflow histograms.
    if (MxAODTool::fillCutFlowWithBookkeeperInfo() == EL::StatusCode::FAILURE)
    {
      return EL::StatusCode::FAILURE;
    }

    m_newFileMetaData = false;
  }

  // Things to do once per file
  if (m_newFileLoaded) {
    m_crossSectionBRfilterEff = -1;

    if (HG::isMC()) {
      // We want the code to fail if the cross section is not defined.
      m_crossSectionBRfilterEff = getCrossSection();

      if (config()->isDefined(Form("kFactor.%d", eventInfo()->mcChannelNumber())))
      { m_crossSectionBRfilterEff *= getKFactor(); }

      if (config()->isDefined(Form("GeneratorEfficiency.%d", eventInfo()->mcChannelNumber())))
      { m_crossSectionBRfilterEff *= getGeneratorEfficiency(); }

      StrV skipMuonsV = config()->getStrV("SkipSavingMuonObjects",{});
      for (auto mcChanNum : skipMuonsV) {
        if ( std::atoi(mcChanNum.Data()) == (int)eventInfo()->mcChannelNumber() )
        { m_skipMuonObjects = true; }
      }

      StrV skipElectronsV = config()->getStrV("SkipSavingElectronObjects",{});
      for (auto mcChanNum : skipElectronsV) {
        if ( std::atoi(mcChanNum.Data()) == (int)eventInfo()->mcChannelNumber() )
        { m_skipElectronObjects = true; }
      }

    }

    m_newFileLoaded = false;
  }

  // Sets the truth channel,  the truth channel in order to be able to use it in the cutflow
  SetTruthZBosonInformation();

  // Set this for every event, just in case.
  var::ZyChannel.setValue(HG::CHANNELUNKNOWN);

  // apply cuts. Returned value will be the last passed cut
  m_cutFlow = cutflow();

  // fill the cut-flow histograms. For each event, fills the bins up to (excluding) the cut it fails.
  double wi = weightInitial();
  for (int cut=ALLEVTS; cut < m_cutFlow; ++cut) {

    // Starting at the choice of the two leptons, multiply by ID/IP SFs
    // (affects all subsequent entries)
    if (cut == ZBOSON_ASSIGNMENT) wi *= m_lepWeight;

    // Starting at photon tight ID
    if (cut == GAM_TIGHTID) wi *= m_phIDWeight;

    // Starting at photon tight iso
    if (cut == GAM_ISOLATION) wi *= m_phIsoWeight;

    // Trigger SF?

    fillCutFlow(CutEnum(cut),wi);
  }

  // if desired, apply skimming.
  // That is, only write events that pass a given cut to the output
  if (m_skimCut >= 1 && m_cutFlow <= m_skimCut)
    return EL::StatusCode::SUCCESS;

  // Selects the objects, does overlap removal, and calculate all
  // variables that will be saved
  doReco();

  // check if we should apply systematics or not
  if (m_applySystematics) {
    for (auto sys: getSystematics()) {
      // ignore nominal case, already done!
      if (sys.name() == "") continue;

      // apply the systmeatic variation and calculate the outupt
      CP_CHECK("RadiativeZCutflowAndMxAOD::execute()", applySystematicVariation(sys));

      // Set this for every event, just in case.
      var::ZyChannel.setValue(HG::CHANNELUNKNOWN);

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
RadiativeZCutflowAndMxAOD::CutEnum RadiativeZCutflowAndMxAOD::cutflow()
{
  m_lepWeight   = 1.0;
  m_phIDWeight  = 1.0;
  m_phIsoWeight = 1.0;

  m_preSelPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selElectrons = xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
  m_selMuons = xAOD::MuonContainer(SG::VIEW_ELEMENTS);
  m_allJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);
  m_selJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);

  //==== CUT 4 : Remove duplicate events (only for data) ====
  static bool checkDuplicates = config()->getBool("EventHandler.CheckDuplicates");
  if ( checkDuplicates && eventHandler()->isDuplicate() ) return DUPLICATE;

  //==== CUT 5 : GRL ====
  static bool requireGRL = config()->getBool("EventHandler.CheckGRL");
  if ( requireGRL && HG::isData() && !eventHandler()->passGRL(eventInfo()) ) return GRL;

  //==== CUT 6 : Require trigger ====
  static bool requireTrigger = config()->getBool("EventHandler.CheckTriggers");
  // passTrigger() will impose the RunNumbers restriction, if specified via EventHandler.RunNumbers.TRIG
  if ( requireTrigger && !eventHandler()->passTriggers() ) return TRIGGER;

  //==== CUT 7 : Detector quality ====
  if ( !(eventHandler()->passLAr (eventInfo()) &&
         eventHandler()->passTile(eventInfo()) &&
         eventHandler()->passSCT (eventInfo()) ) )
    return DQ;

  //==== CUT 8 : Require a vertex ====
  if ( !eventHandler()->passVertex(eventInfo()) ) return VERTEX;

  // Apply electron preselection.
  // HGamCore does not have an electron preselection step, so we make our own here:
  m_allElectrons = electronHandler()->getCorrectedContainer();

  AddElectronDecorations(m_allElectrons);

  static bool ele_requireID = config()->getBool("ElectronHandler.Selection.ApplyPIDCut", true);
  static bool ele_requireIP = config()->getBool("ElectronHandler.Selection.ApplyIPCuts", true);
  static bool ele_requireIso = config()->getBool("ElectronHandler.Selection.ApplyIsoCut", true);

  xAOD::ElectronContainer m_preSelElectrons(SG::VIEW_ELEMENTS);
  for (auto electron : m_allElectrons) {
    if (!electronHandler()->passOQCut(electron)) { continue; }
    if (!electronHandler()->passPtEtaCuts(electron)) { continue; }
    if (!electronHandler()->passHVCut(electron)) { continue; }
    if (ele_requireID  && !electronHandler()->passPIDCut(electron)) { continue; }
    if (ele_requireIP  && !electronHandler()->passIPCuts(electron)) { continue; }
    if (ele_requireIso && !electronHandler()->passIsoCut(electron)) { continue; }
    m_preSelElectrons.push_back(electron);
  }

  // Apply muon preselection.
  // HGamCore does not have a muon preselection step, so we make our own here:
  m_allMuons = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer m_preSelMuons(SG::VIEW_ELEMENTS);

  static bool mu_requireID = config()->getBool("MuonHandler.Selection.ApplyPIDCut", true);
  static bool mu_requireIP = config()->getBool("MuonHandler.Selection.ApplyIPCuts", true);
  static bool mu_requireIso = config()->getBool("MuonHandler.Selection.ApplyIsoCut", true);

  for (auto muon : m_allMuons) {
    if (!muonHandler()->passPtCuts(muon)) { continue; }
    if (mu_requireID  && !muonHandler()->passPIDCut(muon)) { continue; } // This includes MaxEta cut
    if (mu_requireIP  && !muonHandler()->passIPCuts(muon)) { continue; }
    if (mu_requireIso && !muonHandler()->passIsoCut(muon)) { continue; }
    m_preSelMuons.push_back(muon);
  }

  //==== CUT 9 : 2 SF leptons (before OR) ====
  if (m_preSelElectrons.size() < 2 && m_preSelMuons.size() < 2) return TWO_SF_LEPTONS;

  // Get object containers
  m_allPhotons = photonHandler()->getCorrectedContainer();
  AddPhotonDecorations(m_allPhotons);

  // This section is just for the cutflow purposes.
  int nquality=0, namb=0, nHV=0;
  for (auto gam:m_allPhotons) {
    bool passQuality = (photonHandler()->passOQCut(gam)       &&
                        photonHandler()->passCleaningCut(gam) &&
                        photonHandler()->passPtEtaCuts(gam));

    if (!passQuality) continue;
    ++nquality;

    if (!photonHandler()->passAmbCut(gam)) continue;
    ++namb;

    if (!photonHandler()->passHVCut(gam)) continue;
    ++nHV;

    m_preSelPhotons.push_back(gam);
  }

  //==== CUT 10 : Require one loose photon, pT > "PtPreCutGeV" GeV ====
  if (nquality<1) return ONE_RECO_GAM;

  //==== CUT 11 : ambiguity / HV
  // - Require two loose photons that also pass e-gamma ambiguity ====
  static bool requireAmbiguity = config()->getBool("PhotonHandler.Selection.ApplyAmbiguityCut", false);
  if (requireAmbiguity && namb<1) return AMBIGUITY;

  // sneak HV requirement into AMBIGUITY bit
  static bool requireHV = config()->getBool("PhotonHandler.Selection.ApplyHVCut", false);
  if (requireHV && nHV<1) return HVCUT;

  // Our *Higgs candidate photon* is the leading pre-selected (reco-level, quality) photon
  if ( m_preSelPhotons.size() ) m_selPhotons.push_back(m_preSelPhotons[0]);

  // Select Z candidate after overlap removal.
  // Find the highest-pt pair closest to (13 TeV)
  int sel_muon1 = -1, sel_muon2 = -1;
  double return_mmumu = -1;

  float lead_pt_cut = config()->getNum("MuonHandler.Selection.PtLeadCutGeV",11);
  HG::AssignZbosonIndices(m_preSelMuons,sel_muon1,sel_muon2,return_mmumu,
                          /*sort by Z pt? (false is sort-by-mass)*/ false,
                          /*Find value (pt or mass) closest to X (MeV) */ 13000.*HG::GeV,
                          /*leading lepton pt cut (GeV)*/ lead_pt_cut);

  xAOD::ElectronContainer candElectrons(SG::VIEW_ELEMENTS);

  // HG::AssignZbosonIndices( electrons )
  int sel_ele1 = -1, sel_ele2 = -1;
  double return_mee = -1;
  HG::AssignZbosonIndices(m_preSelElectrons,sel_ele1,sel_ele2,return_mee,
                          /*sort by Z pt? (false is sort-by-mass)*/ false,
                          /*Find value (pt or mass) closest to X (MeV)*/ 13000.*HG::GeV);

  //==== CUT 12 : Whether SF leptons survive OR
  if (return_mmumu < 0 && return_mee < 0) return ZBOSON_ASSIGNMENT;

  double m_lly = -999, m_ll = -999;

  if (return_mmumu > 0) {
    xAOD::Muon* mu0 = m_preSelMuons[sel_muon1];
    xAOD::Muon* mu1 = m_preSelMuons[sel_muon2];

    m_selMuons.push_back(mu0);
    m_selMuons.push_back(mu1);
    m_ll = return_mmumu;
    m_lly = (mu0->p4() + mu1->p4() + m_selPhotons[0]->p4()).M();
    var::ZyChannel.setValue(HG::DIMUON);

    if (mu_requireID) {
      // ID scale factors for cutflow
      m_lepWeight *= (HG::MuonHandler::effSF(*mu0) * HG::MuonHandler::effSFTTVA(*mu0));
      m_lepWeight *= (HG::MuonHandler::effSF(*mu1) * HG::MuonHandler::effSFTTVA(*mu1));
    }

    if (mu_requireIso) {
      // Iso scale factors for cutflow
      m_lepWeight *= HG::MuonHandler::effSFIso(*mu0);
      m_lepWeight *= HG::MuonHandler::effSFIso(*mu1);
    }

  }

  else {

    xAOD::Electron* ele0 = m_preSelElectrons[sel_ele1];
    xAOD::Electron* ele1 = m_preSelElectrons[sel_ele2];

    m_selElectrons.push_back(ele0);
    m_selElectrons.push_back(ele1);
    m_ll = return_mee;
    m_lly = (ele0->p4() + ele1->p4() + m_selPhotons[0]->p4()).M();
    var::ZyChannel.setValue(HG::RESOLVED_DIELECTRON);

    if (ele_requireID) {
      // ID / reco scale factors for cutflow
      m_lepWeight *= (HG::ElectronHandler::effIDSF(*ele0) * HG::ElectronHandler::effRecoSF(*ele0));
      m_lepWeight *= (HG::ElectronHandler::effIDSF(*ele1) * HG::ElectronHandler::effRecoSF(*ele1));
    }

    if (ele_requireIso) {
      // Iso scale factors for cutflow
      m_lepWeight *= HG::ElectronHandler::effIsoSF(*ele0);
      m_lepWeight *= HG::ElectronHandler::effIsoSF(*ele1);
    }

  }

  // Need to convert this into photons
  // m_mergedElectronID->decorateMergedVariables(*ele,*trk0,*trk1);
  // HG::EleAcc::passPID(*ele) = m_mergedElectronID->passPIDCut(*ele,*trk0,*trk1);
  // HG::EleAcc::passTMVAPID(*ele) = m_mergedElectronID_v2->passPIDCut(*ele);

  m_allJets = jetHandler()->getCorrectedContainer();
  m_selJets = jetHandler()->applySelection(m_allJets);

  // Removes overlap with candidate photon, and any additional tight photons (if option set)
  overlapHandler()->removeOverlap(m_selPhotons, m_selJets, m_selElectrons, m_selMuons);

  //==== CUT 13 : Whether SF leptons survive OR
  if (m_selElectrons.size() == 0 && m_selMuons.size() < 2) return TWO_SF_LEPTONS_POSTOR;
  if (m_selMuons.size() == 0 && m_selElectrons.size() < 2) return TWO_SF_LEPTONS_POSTOR;

  //==== CUT 14 : Muon cleaning should be done after overlap removal. Cleaning removes isBad muons.
  for (auto muon : m_selMuons) {
    SG::AuxElement::Accessor<char> isBad("isBad");
    if (isBad.isAvailable(*muon) && isBad(*muon)) return BAD_MUON;
  }

  //==== CUT 15 : Photon gets lost in overlap removal
  if (m_selPhotons.size()==0) return ONE_PHOTON_POSTOR;

  //==== CUT 16: Trigger matching (Set in config file using EventHandler.CheckTriggerMatching)
  if ( eventHandler()->doTrigMatch() ){
    bool isTrigMatched = false;
    for (auto trig: eventHandler()->getRequiredTriggers()) {
      // You need passTrigger() here in order to make sure the RunNumbers restriction is imposed on trigs.
      if (eventHandler()->passTrigger(trig) &&
          eventHandler()->passTriggerMatch(trig, &m_selPhotons, &m_selElectrons, &m_selMuons, NULL))
      {
        isTrigMatched = true;
        break;
      }
    }
    if (!isTrigMatched) return TRIG_MATCH;
  }

  xAOD::Photon* selPhoton  = m_selPhotons[0];

  //==== CUT 20 : Require both photons to pass photon ID (isEM) ====
  // Do we really want to require the highest-pt photon to pass tight ID? Can we ask for lower-pt gam?
  static bool ph_requireTight = config()->getBool("PhotonHandler.Selection.ApplyPIDCut", false);
  if (ph_requireTight)
  {
    if (!photonHandler()->passPIDCut(selPhoton)) return GAM_TIGHTID;
    m_phIDWeight *= HG::PhotonHandler::effSF(*selPhoton);
  }

  //==== CUT 21 : Require both photons to fulfill the isolation criteria ===
  static bool ph_requireIso = config()->getBool("PhotonHandler.Selection.ApplyIsoCut", false);
  if (ph_requireIso) {
    if (!photonHandler()->passIsoCut(selPhoton)) return GAM_ISOLATION;
    m_phIsoWeight *= HG::PhotonHandler::isoSF(*selPhoton);
  }

  //==== CUT 22 : Z Mass window cut ====
  if ( m_ll < 10.*HG::GeV ) return ZMASSCUT;

  //==== CUT 23 : lly window cut ====
  if ( 80.*HG::GeV > m_lly || m_lly < 100.*HG::GeV ) return LLGMASSCUT;
  
  return PASSALL;
}

EL::StatusCode  RadiativeZCutflowAndMxAOD::doReco(bool isSys){
  // Do anything you missed in cutflow, and save the objects.

  // Adds event weights and catgory to TStore
  // Also sets pointer to photon container, etc., which is used by var's
  HG::VarHandler::getInstance()->setContainers(&m_selPhotons,&m_selElectrons,&m_selMuons,&m_selJets);

  // Weights
  // total weight (mc, prw, vtx, sf)
  // Reset every time, since weightInitial() might be updated
  // for new pileup weighting
  double myweight = weightInitial();

  // Only apply all scale factors in the case where all objects pass. For all other instances,
  // you have to figure out exactly what weights you need!
  if (m_cutFlow > GAM_ISOLATION) {
    myweight *= m_lepWeight;
    myweight *= m_phIDWeight;
    myweight *= m_phIsoWeight;
    // trigger?
  }

  var::weight.setValue(myweight);

  // Adds event-level variables to TStore
  // Write in the nominal and systematics loops
  writeNominalAndSystematic();
  writeNominalAndSystematicVars();

  // Write only in the nominal loop
  if (not isSys)
  {
    writeNominalOnly();
    writeNominalOnlyVars();

    if (m_saveDetailed) {
      writeDetailed();
      writeDetailedVars();
    }

    if (m_saveObjects) {
      CP_CHECK("execute()", photonHandler  ()->writeContainer(m_selPhotons  ));
      CP_CHECK("execute()", jetHandler     ()->writeContainer(m_selJets     ));
      if (!m_skipElectronObjects) {
        CP_CHECK("execute()", electronHandler()->writeContainer(m_selElectrons));
      }
      if (!m_skipMuonObjects) {
        CP_CHECK("execute()", muonHandler    ()->writeContainer(m_selMuons    ));
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

void RadiativeZCutflowAndMxAOD::writePhotonAllSys(bool isSys)
{
  (void)isSys;

  // Basic event selection flags
  var::cutFlow.setValue(m_cutFlow);

  // Add MC only variables
  if (HG::isMC()) {
    eventHandler()->storeVar<float>("crossSectionBRfilterEff", m_crossSectionBRfilterEff);
  }

}

void RadiativeZCutflowAndMxAOD::writePhotonAllSysVars(bool /*truth*/)
{

}

void RadiativeZCutflowAndMxAOD::writeNominalAndSystematic()
{
  // Basic event selection flags
  var::cutFlow.setValue(m_cutFlow);

  // Additional variables useful for non-framework analysis
  eventHandler()->storeVar<char>("isPassedEventSelection",m_cutFlow >= PASSALL);

}

void RadiativeZCutflowAndMxAOD::writeNominalAndSystematicVars(bool truth)
{
  // Put here all of the HGamVariables that you want to save in the nominal loop and
  // the systematics loops.
  // In the truth case, there is no "systematic" case, so they are saved only once.

  var::m_lly.addToStore(truth);
  var::m_lly_gev.addToStore(truth);

}

void RadiativeZCutflowAndMxAOD::writeTruthOnlyVars()
{
  // Put here the truth-only HGamVariables that you want to save to the MxAOD.

  // bool truth = true;
  // var::deltaR_l1l2_h1.addToStore(truth);
}


void RadiativeZCutflowAndMxAOD::writeNominalOnly()
{
  // Put here the things that you want to save only in the nominal loop (not the systematics loops).

  eventHandler()->mu();
  eventHandler()->runNumber();

  // Make sure every trigger is checked, and decorated to EventInfo
  eventHandler()->getPassedTriggers();

  // Add some convenience variables
  eventHandler()->storeVar<char>("isPassedObjPreselection",m_cutFlow > TRIG_MATCH);
  eventHandler()->storeVar<char>("isPassedObjSelection",m_cutFlow > GAM_ISOLATION);

  // Vertex information
  eventHandler()->numberOfPrimaryVertices();
  eventHandler()->selectedVertexZ();
  eventHandler()->hardestVertexZ();
  eventHandler()->pileupVertexSumPt2(); // also sets pileupVertexZ internally

  if (HG::isMC()) truthHandler()->vertexZ();

  // Bunch train information
  eventHandler()->bunchDistanceFromFront();
  eventHandler()->bunchGapBeforeTrain();

  // Add MC only variables
  if (HG::isMC()) {

    // Basic event weights
    eventHandler()->pileupWeight();
    eventHandler()->vertexWeight();

    truthHandler()->catCoup();
    eventHandler()->storeVar<float>("crossSectionBRfilterEff", m_crossSectionBRfilterEff);
  }

}

void RadiativeZCutflowAndMxAOD::writeNominalOnlyVars(bool truth)
{
  // Put here all of the HGamVariables that you want to save in the nominal loop only
  // (not the systematics loops).

  var::m_ll.addToStore(truth);
  var::deltaR_ll.addToStore(truth);
  var::pt_lly.addToStore(truth);
  var::pt_ll.addToStore(truth);

  var::pTt_lly.addToStore(truth);

  if (!truth)
  {
    // var::m_lly_track4mom.addToStore(false);
  }

  return;
}

void RadiativeZCutflowAndMxAOD::writeDetailed()
{
  // Put here all of the things that you want to save in the case where
  // "SaveDetailedVariables" is set to TRUE in the config file.

}

void RadiativeZCutflowAndMxAOD::writeDetailedVars(bool /*truth*/)
{
  // Put here all of the HGamVariables that you want to save in the case where
  // "SaveDetailedVariables" is set to TRUE in the config file.

}

void RadiativeZCutflowAndMxAOD::SetTruthZBosonInformation(void)
{

  var::ZyChannel.setTruthValue( (int)HG::CHANNELUNKNOWN );

  if (!HG::isMC())
    return;

  // Get all truth particles
  const xAOD::TruthParticleContainer* truthParticles = truthHandler()->getTruthParticles();

  HG::TruthPtcls promptleps(SG::VIEW_ELEMENTS);

  // Works with Sherpa
  for (auto part : *truthParticles) {

    // isGoodTruthElectron and Muon exclude e's and mu's from taus
    bool isGoodE   = (part->absPdgId() == 11 && HG::isGoodTruthElectron(part));
    bool isGoodMu  = (part->absPdgId() == 13 && HG::isGoodTruthMuon(part));
    // We ignore taus for now.

    if (isGoodE || isGoodMu)
    {
      promptleps.push_back(part);
    }
  }

  bool isDilep = (promptleps.size() == 2 && (promptleps[0]->absPdgId() == promptleps[1]->absPdgId()) );
  if (!isDilep)
  {
    return;
  }

  if (promptleps[0]->absPdgId() == 11) {
    var::ZyChannel.setTruthValue( (int)HG::RESOLVED_DIELECTRON );
  }

  if (promptleps[0]->absPdgId() == 13) {
    var::ZyChannel.setTruthValue( (int)HG::DIMUON );
  }

  return;
}

EL::StatusCode  RadiativeZCutflowAndMxAOD::doTruth()
{
  // Truth particles
  xAOD::TruthParticleContainer all_photons   = truthHandler()->getPhotons();
  xAOD::TruthParticleContainer all_electrons = truthHandler()->getElectrons();
  xAOD::TruthParticleContainer all_muons     = truthHandler()->getMuons();
  xAOD::JetContainer           all_jets      = truthHandler()->getJets();
  xAOD::TruthParticleContainer all_higgs     = truthHandler()->getHiggsBosons();

  // Apply fiducial selections to all containers
  xAOD::TruthParticleContainer photons   = truthHandler()->applyPhotonSelection   (all_photons);
  xAOD::TruthParticleContainer electrons = truthHandler()->applyElectronSelection (all_electrons);
  xAOD::TruthParticleContainer muons     = truthHandler()->applyMuonSelection     (all_muons);
  xAOD::JetContainer           jets      = truthHandler()->applyJetSelection      (all_jets);
  xAOD::JetContainer           bjets     = truthHandler()->applyBJetSelection     (jets);

  // remove truth jets that are from electrons or photons
  truthHandler()->removeOverlap(photons, jets, electrons, muons);

  // Save truth containers, if configured
  if (m_saveTruthObjects) {
    truthHandler()->writePhotons    (all_photons  );
    if (!m_skipElectronObjects) truthHandler()->writeElectrons  (all_electrons);
    if (!m_skipMuonObjects)     truthHandler()->writeMuons      (all_muons    );
    truthHandler()->writeJets       (all_jets     );
    truthHandler()->writeHiggsBosons(all_higgs    );
    truthHandler()->writeTruthEvents(             );

    addTruthLinks(m_photonContainerName.Data(), m_photonTruthContainerName.Data());
    if (!m_skipElectronObjects)
    { addTruthLinks(m_elecContainerName.Data()  , m_elecTruthContainerName.Data()); }
  }

  HG::VarHandler::getInstance()->setTruthContainers(&all_photons, &electrons, &muons, &jets);

  // Adds event-level variables to TStore (this time using truth containers)
  bool truth = true;
  if (m_saveTruthVars) {

    writeNominalAndSystematicVars(truth);
    writeNominalOnlyVars(truth);
    writeTruthOnlyVars();
    if (m_saveDetailed)
    {
      writeDetailedVars(truth);
    }

  }

  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}


void RadiativeZCutflowAndMxAOD::fillCutFlow(CutEnum cut, double w) {

  getCutFlowHisto(false /*weighted*/)->Fill(cut);

  if (HG::isData()) return;

  getCutFlowHisto(true)->Fill(cut,w);

  // per-channel cutflows (mainly for debugging)
  HG::ChannelEnum truthChan = (HG::ChannelEnum)var::ZyChannel.truth();

  // unweighted MC, per channel:
  if (truthChan == HG::DIMUON) getCutFlowHisto(false, HG::DIMUON)->Fill(cut);
  if (truthChan == HG::RESOLVED_DIELECTRON) getCutFlowHisto(false, HG::RESOLVED_DIELECTRON)->Fill(cut);

  // weighted MC
  getCutFlowHisto(true)->Fill(cut,w);

  if (truthChan == HG::DIMUON) getCutFlowHisto(true, HG::DIMUON)->Fill(cut,w);
  if (truthChan == HG::RESOLVED_DIELECTRON) getCutFlowHisto(true, HG::RESOLVED_DIELECTRON)->Fill(cut,w);

  return;
}


EL::StatusCode RadiativeZCutflowAndMxAOD::finalize() {
  printf("\nEvent selection cut flow:\n");
  printCutFlowHistos();

  // Write the output to file
  HgammaAnalysis::finalize();

  SafeDelete(m_trackHandler);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RadiativeZCutflowAndMxAOD::fileExecute() {
  // Things that you need to process for each individual input file, even those containing no events.

  HgammaAnalysis::fileExecute();

  // Signals to the code (in execute) that we have to book cutflow information.
  m_newFileMetaData = true;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RadiativeZCutflowAndMxAOD::changeInput(bool firstFile) {
  // Here you do everything you need to do when we change input files.

  HgammaAnalysis::changeInput(firstFile);

  // Signals to the code (in execute) that we have to reset crossSectionBRfilterEff
  m_newFileLoaded = true;

  return EL::StatusCode::SUCCESS;
}

void RadiativeZCutflowAndMxAOD::AddElectronDecorations(xAOD::ElectronContainer& electrons) {

  // Make a dummy vertex container (does not need to be fancy because we do not save it
  // but it does need to be safe for running over systematics (i.e. do not make many copies of it)
  // (to avoid "Trying to overwrite object with key" errors)

  xAOD::VertexContainer* outVertices = nullptr;
  if (!store()->contains<xAOD::VertexContainer>("TempVertices")) {
    outVertices = new xAOD::VertexContainer();
    xAOD::VertexAuxContainer* outVerticesAux = new xAOD::VertexAuxContainer();
    outVertices->setStore(outVerticesAux);
    store()->record( outVertices, "TempVertices" );
    store()->record( outVerticesAux, "TempVerticesAux." );
  }
  else if (store()->retrieve( outVertices , "TempVertices" ).isFailure())
  {
    HG::fatal("Could not retrieve TempVertices. Exiting.");
  }

  for (auto electron : electrons) {

    // NEED TO initialize merged electron ID variables here!
    HG::EleAcc::EOverP0P1(*electron) = -999;
    HG::EleAcc::dRExtrapTrk12(*electron) = -999;
    HG::EleAcc::dRExtrapTrk12_LM(*electron) = -999;
    HG::EleAcc::delta_z0sinTheta_tracks(*electron) = -999;
    HG::EleAcc::delta_z0_tracks(*electron) = -999;
    HG::EleAcc::passPID(*electron) = false;
    HG::EleAcc::passTMVAPID(*electron) = false;

    // HG::EleAcc::dRbetweenTracks_LM_L1(*electron) = -999;
    // HG::EleAcc::dRbetweenTracks_LM_L2(*electron) = -999;
    // HG::EleAcc::dRbetweenTracks_P_L1(*electron) = -999;
    // HG::EleAcc::dRbetweenTracks_P_L2(*electron) = -999;

    // Decorate Rhad
    double feta = fabs(electron->eta());
    HG::EleAcc::RhadForPID(*electron) = (0.8 < feta && feta < 1.37) ?
      electron->showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rhad) :
      electron->showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rhad1);



    int index1 = HG::EleAcc::vtxTrkIndex1(*electron);
    int index2 = HG::EleAcc::vtxTrkIndex2(*electron);
    HG::EleAcc::vtxTrk1_TRT_PID_trans(*electron) = index1<0 ? -999 : trackHandler()->calculateTRT_PID(*electron->trackParticle(index1)) ;
    HG::EleAcc::vtxTrk2_TRT_PID_trans(*electron) = index2<0 ? -999 : trackHandler()->calculateTRT_PID(*electron->trackParticle(index2)) ;

    xAOD::Photon* photon = HG::createPhotonFromElectron(electron);


    if(photon)
    {
      HG::setPhotonConversionVertex( electron, photon, 20, outVertices);
      HG::EleAcc::calibratedPhotonEnergy(*electron) = photon->e();

      delete photon;
    } else {
      HG::EleAcc::calibratedPhotonEnergy(*electron) = -999;
    }

  }

    return;
}

void RadiativeZCutflowAndMxAOD::AddPhotonDecorations(xAOD::PhotonContainer& photons) {


  for (auto photon : photons) {

    double feta = fabs(photon->eta());

    HG::PhAcc::RhadForPID(*photon) = (0.8 < feta && feta < 1.37) ?
      photon->showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rhad) :
      photon->showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rhad1);

    HG::PhAcc::ambiguousE_deltaEta1(*photon) = -999;

    if (photon->ambiguousObject()) {
      const xAOD::Electron* el = dynamic_cast<const xAOD::Electron*>( photon->ambiguousObject() );
      HG::PhAcc::ambiguousE_deltaEta1(*photon) = el->trackCaloMatchValue(xAOD::EgammaParameters::deltaEta1);
    }

    Amg::Vector3D vtxMom = xAOD::EgammaHelpers::momentumAtVertex(photon);
    HG::PhAcc::vtxE(*photon) = sqrt(vtxMom.x()*vtxMom.x() + vtxMom.y()*vtxMom.y() + vtxMom.z()*vtxMom.z());

    HG::PhAcc::vtxTrk1_TRT_PID_trans(*photon) = -999;
    HG::PhAcc::vtxTrk2_TRT_PID_trans(*photon) = -999;

    const xAOD::Vertex *vxPh = photon->vertex();

    if (vxPh){
      const xAOD::TrackParticle *trk1 = ( vxPh->nTrackParticles() ? vxPh->trackParticle(0) : 0 );
      const xAOD::TrackParticle *trk2 = ( vxPh->nTrackParticles() > 1 ? vxPh->trackParticle(1) : 0 );

      if (trk1) {
        HG::PhAcc::vtxTrk1_TRT_PID_trans(*photon) = trackHandler()->calculateTRT_PID(*trk1);
      }
      if (trk2) {
        HG::PhAcc::vtxTrk2_TRT_PID_trans(*photon) = trackHandler()->calculateTRT_PID(*trk2);
      }
    }

  }

  return;
}

void RadiativeZCutflowAndMxAOD::printCutFlowHistos() {
  for ( auto entry : m_cFlowHistos ) {

    TString isMC = HG::isMC()?"MC sample":"Data run";

    int cutFlowID = entry.first;

    TString typeOfEvents = "all events";

    // Weighted / unweighted
    TString weighted = "";
    int nDecimals = 0; // for unweighted
    if (cutFlowID > (int)1e8) {
      weighted = "(weighted)";
      nDecimals = 2;
      cutFlowID -= (int)1e8;
    }

    // figure out channel
    TString truthChannel = "All channels";
    for (auto chan : {HG::RESOLVED_DIELECTRON,HG::DIMUON}) {
      if (cutFlowID > (int)(1e6 * (int)chan))
      {
        truthChannel = GetChannelName(chan) + TString(" truth channel");
        cutFlowID -= (int)(1e6 * (int)chan);
        break;
      }
    }

    printf("\n%s %d, %s, %s %s\n",isMC.Data(),cutFlowID,typeOfEvents.Data(),
           truthChannel.Data(),weighted.Data());
    printCutFlowHisto(entry.second,nDecimals);
  }
}
