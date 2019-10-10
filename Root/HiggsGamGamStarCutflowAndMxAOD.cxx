#include "HGamGamStar/HiggsGamGamStarCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>
#include "HGamAnalysisFramework/TruthUtils.h"

#include "HGamGamStar/ExtraHggStarObjects.h"
#include "HGamGamStar/TrackElectronMap.h"
#include "xAODTracking/VertexAuxContainer.h"



// #include "PhotonVertexSelection/PhotonPointingTool.h"
// #include "ZMassConstraint/ConstraintFit.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HiggsGamGamStarCutflowAndMxAOD)

HiggsGamGamStarCutflowAndMxAOD::HiggsGamGamStarCutflowAndMxAOD(const char *name)
: MxAODTool(name)
  , m_trackHandler(nullptr)
{ }

HiggsGamGamStarCutflowAndMxAOD::~HiggsGamGamStarCutflowAndMxAOD() {}

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::initialize()
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

  // We will evaluate the close-by correction for every working point, provided that it falls
  // into the Resovled category.
  StrV eleIsoWPs = config()->getStrV("ElectronHandler.Selection.IsoCriteria");

  // Working points used for cutflow
  m_eleMergedIsoWP = electronHandler()->getIsoType(config()->getStr("MergedElectrons.IsoCriteria"));
  m_eleResolvedIsoWP = electronHandler()->getIsoType(config()->getStr("ResolvedElectrons.IsoCriteria"));

  for (auto isoStr : eleIsoWPs) {
    HG::Iso::IsolationType iso = electronHandler()->getIsoType(isoStr);
    ToolHandle<CP::IIsolationSelectionTool> isoTool = electronHandler()->getIsoTool(iso);
    m_isoCloseByTools_Ele[iso] = new CP::IsolationCloseByCorrectionTool(("CBT_ele_"+isoStr).Data());
    m_isoCloseByTools_Ele[iso]->setProperty("IsolationSelectionTool", isoTool);
    m_isoCloseByTools_Ele[iso]->setProperty("BackupPrefix", "original");
    m_isoCloseByTools_Ele[iso]->initialize();

    m_eleIsoAccCorr[iso] = new SG::AuxElement::Accessor<char>(("isIsoWithCorr" + isoStr).Data());
  }

  StrV muIsoWPs = config()->getStrV("MuonHandler.Selection.IsoCriteria");
  m_muonIsoWP = muonHandler()->getIsoType(muIsoWPs[0]); // default WP is used for cut

  for (auto isoStr : muIsoWPs) {
    HG::Iso::IsolationType iso = muonHandler()->getIsoType(isoStr);
    ToolHandle<CP::IIsolationSelectionTool> isoTool = muonHandler()->getIsoTool(iso);
    m_isoCloseByTools_Muon[iso] = new CP::IsolationCloseByCorrectionTool(("CBT_mu_"+isoStr).Data());
    m_isoCloseByTools_Muon[iso]->setProperty("IsolationSelectionTool", isoTool);
    m_isoCloseByTools_Muon[iso]->setProperty("BackupPrefix", "original");
    m_isoCloseByTools_Muon[iso]->initialize();

    m_muIsoAccCorr[iso] = new SG::AuxElement::Accessor<char>(("isIsoWithCorr" + isoStr).Data());
  }

  // Resolved electron ID preselection. Applied in FindZboson_ElectronChannelAware.
  m_eleIDPreselection = config()->getStr("ResolvedElectrons.Preselection.PID","VeryLoose");
  
  //configurable overlap removal cone sizes
  m_OR_e_DR_y = config()->getNum("OverlapRemoval.Electron_DR_Photon", 0.4);
  m_OR_jet_DR_y = config()->getNum("OverlapRemoval.Jet_DR_Photon", 0.4);
  m_OR_jet_DR_e = config()->getNum("OverlapRemoval.Jet_DR_Electron", 0.2);
  m_OR_e_DR_jet = config()->getNum("OverlapRemoval.Electron_DR_Jet", 0.4);
  m_OR_mu_DR_y = config()->getNum("OverlapRemoval.Muon_DR_Photon", 0.4);
  m_OR_mu_DR_jet = config()->getNum("OverlapRemoval.Muon_DR_Jet", 0.4);
  
  m_passTriggers = false;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::createOutput()
{
  // Read the output branch names - add option to make this configurable in future ?
  m_photonContainerName = "HGam"+config()->getStr("PhotonHandler.ContainerName");
  m_jetContainerName    = "HGam"+config()->getStr("JetHandler.ContainerName");
  m_elecContainerName   = "HGam"+config()->getStr("ElectronHandler.ContainerName");
  m_muonContainerName   = "HGam"+config()->getStr("MuonHandler.ContainerName");
  m_trackContainerName  = "HGam"+config()->getStr("TrackHandler.ContainerName");
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

  ignore = {};
  declareOutputVariables(m_trackContainerName , "MxAOD.Variables.Track"   , {}, ignore);

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

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::execute()
{
  // Needed for all underlaying tools to be working corectly!
  HgammaAnalysis::execute();

  // Clear containers which point to objects from previous event
  HG::ExtraHggStarObjects::getInstance()->clearContainers();

  // Handle File Metadata (need to put this here because we need sample ID to define cutflow histo
  if (m_newFileMetaData) {

    // Initialize cutflow histograms by calling them.
    getCutFlowHisto(false /*weighted*/, false /*onlyDalitz*/);

    // per-channel cutflows (mainly for debugging)
    getCutFlowHisto(false, false, HG::DIMUON);
    getCutFlowHisto(false, false, HG::RESOLVED_DIELECTRON);
    getCutFlowHisto(false, false, HG::MERGED_DIELECTRON);

    if (HG::isMC()) {
      getCutFlowHisto(true, false);
      getCutFlowHisto(true, true );

      getCutFlowHisto(true, true , HG::DIMUON);
      getCutFlowHisto(true, true , HG::RESOLVED_DIELECTRON);
      getCutFlowHisto(true, true , HG::MERGED_DIELECTRON);

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
  // m_isNonHyyStarHiggs is set in SetTruthHiggsInformation().
  SetTruthHiggsInformation();

  // Set this for every event, just in case.
  var::yyStarChannel.setValue(HG::CHANNELUNKNOWN);

  // apply cuts. Returned value will be the last passed cut
  m_cutFlow = cutflow();

  // fill the cut-flow histograms. For each event, fills the bins up to (excluding) the cut it fails.
  double wi = weightInitial();
  for (int cut=ALLEVTS; cut < m_cutFlow; ++cut) {

    // Starting at the pass-leptonID, multiply by ID/IP SFs
    // (affects all subsequent entries)
    if (cut == LEP_MEDID) wi *= m_lepIDWeight;

    // Starting at the pass-leptonIso level, multiply by the Iso SF.
    if (cut == LEP_ISO) wi *= m_lepIsoWeight;

    // Starting at photon tight ID
    if (cut == GAM_TIGHTID) wi *= HG::PhotonHandler::effSF(*(m_selPhotons[0]));

    // Starting at photon tight iso
    if (cut == GAM_ISOLATION) wi *= HG::PhotonHandler::isoSF(*(m_selPhotons[0]));

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
      CP_CHECK("HiggsGamGamStarCutflowAndMxAOD::execute()", applySystematicVariation(sys));

      // Set this for every event, just in case.
      var::yyStarChannel.setValue(HG::CHANNELUNKNOWN);

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
  m_passTriggers = false;
  m_lepIDWeight = 1.0;
  m_lepIsoWeight = 1.0;

  m_preSelPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  m_selElectrons = xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
  m_selTracks = xAOD::TrackParticleContainer(SG::VIEW_ELEMENTS);
  m_selMuons = xAOD::MuonContainer(SG::VIEW_ELEMENTS);
  m_allJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);
  m_selJets = xAOD::JetContainer(SG::VIEW_ELEMENTS);

  //==== CUT 3 : Remove truth Higgs events that are not leptonic Dalitz ====
  if ( m_isNonHyyStarHiggs ) return HIGGS_LEP_DALITZ;

  //==== CUT 4 : Remove duplicate events (only for data) ====
  static bool checkDuplicates = config()->getBool("EventHandler.CheckDuplicates");
  if ( checkDuplicates && eventHandler()->isDuplicate() ) return DUPLICATE;

  //==== CUT 5 : GRL ====
  static bool requireGRL = config()->getBool("EventHandler.CheckGRL");
  if ( requireGRL && HG::isData() && !eventHandler()->passGRL(eventInfo()) ) return GRL;

  //==== CUT 6 : Require trigger ====
  m_passTriggers = eventHandler()->passTriggers();
  static bool requireTrigger = config()->getBool("EventHandler.CheckTriggers");
  // passTrigger() will impose the RunNumbers restriction, if specified via EventHandler.RunNumbers.TRIG
  if ( requireTrigger && !m_passTriggers ) return TRIGGER;

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

  xAOD::ElectronContainer m_preSelElectrons(SG::VIEW_ELEMENTS);
  for (auto electron : m_allElectrons) {
    if (!electronHandler()->passOQCut(electron)) { continue; }
    if (!electronHandler()->passPtEtaCuts(electron)) { continue; }
    if (!electronHandler()->passHVCut(electron)) { continue; }
    m_preSelElectrons.push_back(electron);
  }

  m_allTracks = trackHandler()->getCorrectedContainer();
  HG::TrackElectronMap trkElectronMap;
  m_preSelTracks = trackHandler()->findTracksFromElectrons(m_allTracks,m_preSelElectrons,trkElectronMap);
  m_preSelTracks.sort(HG::TrackHandler::comparePt);

  // Apply muon preselection.
  // HGamCore does not have a muon preselection step, so we make our own here:
  m_allMuons = muonHandler()->getCorrectedContainer();

  AddMuonDecorations(m_allMuons);

  xAOD::MuonContainer m_preSelMuons(SG::VIEW_ELEMENTS);
  for (auto muon : m_allMuons) {
    if (!muonHandler()->passPtCuts(muon)) { continue; }
    if (!muonHandler()->passPIDCut(muon)) { continue; } // This includes MaxEta cut
    m_preSelMuons.push_back(muon);
  }

  //==== CUT 9 : 2 SF leptons (before OR) ====
  if (m_preSelTracks.size() < 2 && m_preSelMuons.size() < 2) return TWO_SF_LEPTONS;

  // Get object containers
  m_allPhotons = photonHandler()->getCorrectedContainer();

  // photon applyPreSelection applies OQ, cleaning, PtEta, Loose PID, ambiguity and HV cuts
  m_preSelPhotons = photonHandler()->applyPreSelection(m_allPhotons);

  // Our *Higgs candidate photon* is the leading pre-selected (Loose) photon
  if (m_preSelPhotons.size()  ) m_selPhotons.push_back(m_preSelPhotons[0]);


  // This section is just for the cutflow purposes.
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
        photonHandler()->passAmbCut(gam)       &&
        photonHandler()->passHVCut(gam))
      ++nHV;
  }

  //==== CUT 10 : Require one loose photon, pT > "PtPreCutGeV" GeV ====
  if (nloose<1) return ONE_LOOSE_GAM;

  //==== CUT 11 : ambiguity / HV
  // - Require two loose photons that also pass e-gamma ambiguity ====
  static bool requireAmbiguity = config()->getBool("PhotonHandler.Selection.ApplyAmbiguityCut", false);
  if (requireAmbiguity && namb<1) return AMBIGUITY;

  // sneak HV requirement into AMBIGUITY bit
  static bool requireHV = config()->getBool("PhotonHandler.Selection.ApplyHVCut", false);
  if (requireHV && nHV<1) return AMBIGUITY;

  // Muon cleaning should be done after overlap removal. Cleaning removes isBad muons.
  // m_preSelMuons = muonHandler()->applyCleaningSelection(dirtyMuons);

  // Select Z candidate after overlap removal.
  // Find the highest-pt pair closest to (13 TeV)
  int sel_muon1 = -1, sel_muon2 = -1;
  double return_mmumu = -1;
  float lead_pt_cut = config()->getNum("MuonHandler.Selection.PtLeadCutGeV",11);
  HG::AssignZbosonIndices(m_preSelMuons,sel_muon1,sel_muon2,return_mmumu,/*sortby_pt*/ true,13000.*HG::GeV,lead_pt_cut);

  xAOD::TrackParticle* sel_trk1 = nullptr;
  xAOD::TrackParticle* sel_trk2 = nullptr;
  double return_mtrktrk = -1;

  xAOD::ElectronContainer candElectrons(SG::VIEW_ELEMENTS);
  HG::ChannelEnum echan = FindZboson_ElectronChannelAware(&m_preSelTracks,sel_trk1,sel_trk2,return_mtrktrk,
                                                          trkElectronMap,&m_preSelElectrons,&candElectrons);

  //==== CUT 12 : Whether SF leptons survive OR
  if (return_mmumu < 0 && return_mtrktrk < 0) return ZBOSON_ASSIGNMENT;

  double m_lly = -999, m_ll = -999, pt_ll = -999;

  if (return_mmumu > 0) {
    m_selMuons.push_back(m_preSelMuons[sel_muon1]);
    m_selMuons.push_back(m_preSelMuons[sel_muon2]);
    m_ll = return_mmumu;
    pt_ll = (m_selMuons[0]->p4() + m_selMuons[1]->p4()).Pt();
    m_lly = (m_selMuons[0]->p4() + m_selMuons[1]->p4() + m_selPhotons[0]->p4()).M();
    var::yyStarChannel.setValue(HG::DIMUON);
  } else {
    if (!sel_trk1 || !sel_trk2) HG::fatal("Pointer error. Check code.");
    m_selTracks.push_back(sel_trk1);
    m_selTracks.push_back(sel_trk2);

    for (auto ele : candElectrons) {
      m_selElectrons.push_back(ele);
    }
    m_selElectrons.sort(HG::ElectronHandler::comparePt);
    var::yyStarChannel.setValue(echan);

    if (m_selElectrons.size() == 1) {
      TLorentzVector merged = HG::MergedEleTLV(*m_selTracks[0],*m_selTracks[1],*m_selElectrons[0]);
      m_ll = merged.M();
      pt_ll = merged.Pt();
      m_lly = (merged + m_selPhotons[0]->p4()).M();
      // Decorate merged variables as soon as you find out the channel is merged
      xAOD::Electron* ele = m_selElectrons[0];
      xAOD::TrackParticle* trk0 = m_selTracks[0];
      xAOD::TrackParticle* trk1 = m_selTracks[1];
      m_mergedElectronID->decorateMergedVariables(*ele,*trk0,*trk1);
      // Need to cut on lead track pt in merged case when processing MxAODs, for that need a variable which is defined for each event (i.e. also need it in muon/resolved events where the value is -99)
      HG::EleAcc::passPID(*ele) = m_mergedElectronID->passPIDCut(*ele,*trk0,*trk1);
      HG::EleAcc::passTMVAPID(*ele) = m_mergedElectronID_v2->passPIDCut(*ele, HG::isMC() );
    }
    else if (m_selElectrons.size() == 2) {
      m_ll = (m_selElectrons[0]->p4() + m_selElectrons[1]->p4()).M();
      pt_ll = (m_selElectrons[0]->p4() + m_selElectrons[1]->p4()).Pt();
      m_lly = (m_selElectrons[0]->p4() + m_selElectrons[1]->p4() + m_selPhotons[0]->p4()).M();
    }

  }

  decorateCorrectedIsoCut(m_selElectrons, m_selMuons);

  m_allJets = jetHandler()->getCorrectedContainer();
  m_selJets = jetHandler()->applySelection(m_allJets);

  unsigned int electrons_preOR = m_selElectrons.size();

  // Removes overlap with candidate photon, and any additional tight photons (if option set)
//   overlapHandler()->removeOverlap(m_selPhotons, m_selJets, m_selElectrons, m_selMuons);
  //Not using HGam code directly (line above) to be able to prioritize muon in muon-jet OR (otherwise same)

  // 1. remove electrons overlapping with photons
  overlapHandler()->removeOverlap(m_selElectrons, m_selPhotons,  m_OR_e_DR_y);
  // 2. jets
  // 2.a remove jets overlapping with photons
  overlapHandler()->removeOverlap(m_selJets, m_selPhotons, m_OR_jet_DR_y);
  // 2.b remove jets overlapping with electrons
  overlapHandler()->removeOverlap(m_selJets, m_selElectrons, m_OR_jet_DR_e);
  // 3. remove electrons too close to jets (usually 0.4)
  overlapHandler()->removeOverlap(m_selElectrons, m_selJets, m_OR_e_DR_jet);
  // 4. remove muons overlapping photons and **remove jets** overlapping with muons
  overlapHandler()->removeOverlap(m_selMuons,  m_selPhotons, m_OR_mu_DR_y);
  overlapHandler()->removeOverlap(m_selJets, m_selMuons, m_OR_mu_DR_jet);

  //==== CUT 13 : Whether SF leptons survive OR
  if (m_selElectrons.size() == 0 && m_selMuons.size() < 2) return TWO_SF_LEPTONS_POSTOR;
  if (m_selMuons.size() == 0 && m_selElectrons.size() != electrons_preOR) return TWO_SF_LEPTONS_POSTOR;

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


  if(var::yyStarChannel()==HG::DIMUON){

    xAOD::Muon* mu0 = m_selMuons[0];
    xAOD::Muon* mu1 = m_selMuons[1];

    //==== CUT 17: Require muons to pass medium PID
    static bool requireMedium = config()->getBool("MuonHandler.Selection.ApplyPIDCut", true);
    if (requireMedium) {
      if (!muonHandler()->passPIDCut(mu0) || !muonHandler()->passPIDCut(mu1)) return LEP_MEDID;

      // ID scale factors for cutflow
      m_lepIDWeight *= (HG::MuonHandler::effSF(*mu0) * HG::MuonHandler::effSFTTVA(*mu0));
      m_lepIDWeight *= (HG::MuonHandler::effSF(*mu1) * HG::MuonHandler::effSFTTVA(*mu1));
    }

    //==== CUT 18: Require muons to pass IP
    static bool requireIP = config()->getBool("MuonHandler.Selection.ApplyIPCuts", true);
    if (requireIP) {
      if (!muonHandler()->passIPCuts(mu0) || !muonHandler()->passIPCuts(mu1)) return LEP_IP;
    }

    //==== CUT 19: Require muons to pass isolation
    static bool requireIso = config()->getBool("MuonHandler.Selection.ApplyIsoCut", true);
    if(requireIso){
      static bool correctIsolation = config()->getBool("MuonHandler.Selection.UseCorrectedIso", false);
      if(correctIsolation){ //isolation cut taking into account close-by objects
        if ( !(*m_muIsoAccCorr[m_muonIsoWP])(*mu0) ||
             !(*m_muIsoAccCorr[m_muonIsoWP])(*mu1)) return LEP_ISO;
      }
      else {
        if ( !muonHandler()->passIsoCut(mu0) ||
             !muonHandler()->passIsoCut(mu1)) return LEP_ISO;
      }

      // Iso scale factors for cutflow
      m_lepIsoWeight *= HG::MuonHandler::effSFIso(*mu0);
      m_lepIsoWeight *= HG::MuonHandler::effSFIso(*mu1);

    }

  }

  else if(var::yyStarChannel()==HG::RESOLVED_DIELECTRON){

    xAOD::Electron* ele0 = m_selElectrons[0];
    xAOD::Electron* ele1 = m_selElectrons[1];

    //==== CUT 17: Require electrons to pass medium PID
    static bool requireMedium = config()->getBool("ElectronHandler.Selection.ApplyPIDCut", true);
    if (requireMedium) {
      if (!electronHandler()->passPIDCut(ele0) || !electronHandler()->passPIDCut(ele1)) return LEP_MEDID;

      // ID / reco scale factors for cutflow
      m_lepIDWeight *= (HG::ElectronHandler::effIDSF(*ele0) * HG::ElectronHandler::effRecoSF(*ele0));
      m_lepIDWeight *= (HG::ElectronHandler::effIDSF(*ele1) * HG::ElectronHandler::effRecoSF(*ele1));
    }

    //==== CUT 18: Require electrons to pass IP
    static bool requireIP = config()->getBool("ElectronHandler.Selection.ApplyIPCuts", true);
    if (requireIP) {
      if (!electronHandler()->passIPCuts(ele0) || !electronHandler()->passIPCuts(ele1)) return LEP_IP;
    }

    //==== CUT 19: Require electrons to pass isolation
    static bool requireIso = config()->getBool("ElectronHandler.Selection.ApplyIsoCut", true);
    if(requireIso){
      static bool correctIsolation = config()->getBool("ElectronHandler.Selection.UseCorrectedIso", false);
      if(correctIsolation){ //isolation cut taking into account close-by objects
        if ( !(*m_eleIsoAccCorr[m_eleResolvedIsoWP])(*ele0) ||
             !(*m_eleIsoAccCorr[m_eleResolvedIsoWP])(*ele1)) return LEP_ISO;
      }
      else {
        if ( !electronHandler()->passIsoCut(ele0,m_eleResolvedIsoWP) ||
             !electronHandler()->passIsoCut(ele1,m_eleResolvedIsoWP)) return LEP_ISO;
      }

      // Iso scale factors for cutflow
      m_lepIsoWeight *= HG::ElectronHandler::effIsoSF(*ele0);
      m_lepIsoWeight *= HG::ElectronHandler::effIsoSF(*ele1);

    }
  }

  else if(var::yyStarChannel()==HG::MERGED_DIELECTRON){

    //==== CUT 17: Require electrons to pass merged PID
    static bool requireMerged = config()->getBool("ElectronHandler.Selection.ApplyPIDCut", true);
    if (requireMerged && !HG::EleAcc::passTMVAPID(*m_selElectrons[0]) ) return LEP_MEDID;

    //==== CUT 18: Require electrons to pass IP
    static bool requireIP = config()->getBool("ElectronHandler.Selection.ApplyIPCuts", true);
    if (requireIP && (!trackHandler()->passIPCuts(*m_selTracks[0]) || !trackHandler()->passIPCuts(*m_selTracks[1])) ) return LEP_IP;

    //==== CUT 19: Require melectrons to pass isolation
    static bool requireIso = config()->getBool("ElectronHandler.Selection.ApplyIsoCut", true);
    if (requireIso && (!electronHandler()->passIsoCut(m_selElectrons[0],m_eleMergedIsoWP)) ) return LEP_ISO;
  }
  else {
    HG::fatal("Unknown channel categorization - please check!");
  }

  //==== CUT 20 : Require both photons to pass photon ID (isEM) ====
  // Do we really want to require the highest-pt photon to pass tight ID? Can we ask for lower-pt gam?
  static bool requireTight = config()->getBool("PhotonHandler.Selection.ApplyPIDCut", true);
  if (requireTight && (!photonHandler()->passPIDCut(m_selPhotons[0])) ) return GAM_TIGHTID;

  //==== CUT 21 : Require both photons to fulfill the isolation criteria ===
  static bool requireIso = config()->getBool("PhotonHandler.Selection.ApplyIsoCut", true);
  if (requireIso && (!photonHandler()->passIsoCut(m_selPhotons[0]))) return GAM_ISOLATION;

  //==== CUT 22 : Z Mass window cut ====
  if ( m_ll > 45.*HG::GeV ) return ZMASSCUT;

  //==== CUT 23 : lly window cut ====
  if ( 105.*HG::GeV > m_lly || m_lly > 160.*HG::GeV ) return LLGMASSCUT;
  
  //==== CUT 24 : mll window (veto J/Psi and Y) cut ====
  if(var::yyStarChannel()==HG::DIMUON){
    if ( (m_ll > 2.9*HG::GeV && m_ll < 3.3*HG::GeV) || (m_ll > 9.1*HG::GeV && m_ll < 10.6*HG::GeV) ) return LLMASSCUT;
  }
  else{
    if ( (m_ll > 2.5*HG::GeV && m_ll < 3.5*HG::GeV) || (m_ll > 8.0*HG::GeV && m_ll < 11.0*HG::GeV) ) return LLMASSCUT;
  }
  
  //==== CUT 25 : pt_ll fraction ggF cut ==== 
  if (pt_ll/m_lly < 0.3) return DILEP_PT_FRAC;

  //==== CUT 26 : pt_y fraction ggF cut ==== 
  if (m_selPhotons[0]->pt()/m_lly < 0.3) return GAM_PT_FRAC;

  return PASSALL;
}

EL::StatusCode  HiggsGamGamStarCutflowAndMxAOD::doReco(bool isSys){
  // Do anything you missed in cutflow, and save the objects.

  // Adds event weights and catgory to TStore
  // Also sets pointer to photon container, etc., which is used by var's
  HG::VarHandler::getInstance()->setContainers(&m_selPhotons,&m_selElectrons,&m_selMuons,&m_selJets);
  HG::ExtraHggStarObjects::getInstance()->setElectronTrackContainer(&m_selTracks);

  // Weights
  // total weight (mc, prw, vtx, sf)
  // Reset every time, since weightInitial() might be updated
  // for new pileup weighting
  double myweight = weightInitial();

  // Only apply all scale factors in the case where all objects pass. For all other instances,
  // you have to figure out exactly what weights you need!
  if (m_cutFlow > GAM_ISOLATION) {
    myweight *= m_lepIDWeight;
    myweight *= m_lepIsoWeight;
    myweight *= HG::PhotonHandler::effSF(*(m_selPhotons[0]));
    myweight *= HG::PhotonHandler::isoSF(*(m_selPhotons[0]));
    // trigger?
  }

  var::weight.setValue(myweight);

  // Set the Merged electron TLV
  if (var::yyStarChannel() == HG::MERGED_DIELECTRON) {
    HG::ExtraHggStarObjects::getInstance()->setMergedElectronTLV(*m_selTracks[0],*m_selTracks[1],*m_selElectrons[0]);
  }

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
      CP_CHECK("execute()", trackHandler   ()->writeContainer(m_selTracks   ));
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

void HiggsGamGamStarCutflowAndMxAOD::writePhotonAllSys(bool isSys)
{
  (void)isSys;

  // Basic event selection flags
  var::cutFlow.setValue(m_cutFlow);

  // Add MC only variables
  if (HG::isMC()) {
    eventHandler()->storeVar<float>("crossSectionBRfilterEff", m_crossSectionBRfilterEff);
  }

}

void HiggsGamGamStarCutflowAndMxAOD::writePhotonAllSysVars(bool /*truth*/)
{

}

void HiggsGamGamStarCutflowAndMxAOD::writeNominalAndSystematic()
{
  // Basic event selection flags
  var::cutFlow.setValue(m_cutFlow);

  // Additional variables useful for non-framework analysis
  eventHandler()->storeVar<char>("isPassedEventSelection",m_cutFlow >= PASSALL);

}

void HiggsGamGamStarCutflowAndMxAOD::writeNominalAndSystematicVars(bool truth)
{
  // Put here all of the HGamVariables that you want to save in the nominal loop and
  // the systematics loops.
  // In the truth case, there is no "systematic" case, so they are saved only once.

  var::m_lly.addToStore(truth);
  var::m_lly_gev.addToStore(truth);
  var::yyStarCategory.addToStore(truth);

}

void HiggsGamGamStarCutflowAndMxAOD::writeTruthOnlyVars()
{
  // Put here the truth-only HGamVariables that you want to save to the MxAOD.

  bool truth = true;

  var::pT_h1.addToStore(truth);
  var::y_h1.addToStore(truth);
  var::m_h1.addToStore(truth);

  var::pT_l1_h1.addToStore(truth);
  var::pT_l2_h1.addToStore(truth);
  var::deltaR_l1l2_h1.addToStore(truth);
  var::ystar_pdg_flavor.addToStore(truth);
  var::pT_yDirect_h1.addToStore(truth);
  var::m_yStar_undressed_h1.addToStore(truth);

}


void HiggsGamGamStarCutflowAndMxAOD::writeNominalOnly()
{
  // Put here the things that you want to save only in the nominal loop (not the systematics loops).

  eventHandler()->mu();
  eventHandler()->runNumber();

  // Make sure every trigger is checked, and decorated to EventInfo
  eventHandler()->getPassedTriggers();
  eventHandler()->storeVar<char>("isPassedTriggers",m_passTriggers);

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

void HiggsGamGamStarCutflowAndMxAOD::writeNominalOnlyVars(bool truth)
{
  // Put here all of the HGamVariables that you want to save in the nominal loop only
  // (not the systematics loops).

  var::m_ll.addToStore(truth);
  var::deltaR_ll.addToStore(truth);
  var::pt_lly.addToStore(truth);
  var::pt_ll.addToStore(truth);

  var::m_jj.addToStore(truth);
  var::Deta_j_j.addToStore(truth);
  var::Dphi_lly_jj.addToStore(truth);
  var::Zepp_lly.addToStore(truth);
  var::pTt_lly.addToStore(truth);
  var::pT_llyjj.addToStore(truth);
  var::DRmin_y_ystar_2jets.addToStore(truth);
  var::DRmin_y_leps_2jets.addToStore(truth);

  if (!truth)
  {
    var::m_lly_track4mom.addToStore(false);
    var::m_ll_track4mom.addToStore(false);
    var::Resolved_dRExtrapTrk12.addToStore(false);
    var::Resolved_deltaPhiRescaled2.addToStore(false);
    var::Resolved_deltaEta2.addToStore(false);
    var::trk_lead_pt.addToStore(false);
  }

  return;
}

void HiggsGamGamStarCutflowAndMxAOD::writeDetailed()
{
  // Put here all of the things that you want to save in the case where
  // "SaveDetailedVariables" is set to TRUE in the config file.

}

void HiggsGamGamStarCutflowAndMxAOD::writeDetailedVars(bool /*truth*/)
{
  // Put here all of the HGamVariables that you want to save in the case where
  // "SaveDetailedVariables" is set to TRUE in the config file.

}

void HiggsGamGamStarCutflowAndMxAOD::SetTruthHiggsInformation(void)
{
  // For data this is automatically false:
  m_isNonHyyStarHiggs = false;

  if (!HG::isMC())
    return;

  xAOD::TruthParticleContainer all_higgs = truthHandler()->getHiggsBosons();
  HG::VarHandler::getInstance()->setHiggsBosons(&all_higgs);

  // Set the truth decay product containers in ExtraHggStarObjects
  const xAOD::TruthParticleContainer* all_particles = truthHandler()->getTruthParticles();
  HG::ExtraHggStarObjects::getInstance()->setTruthHiggsDecayProducts(all_particles);

  // Code for the Truth Channel determination
  const xAOD::ElectronContainer all_reco_elecs = electronHandler()->getCorrectedContainer();
  auto childleps = *(HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons());

  HG::ChannelEnum truthChannel = HG::truthChannel(childleps,all_reco_elecs);
  var::yyStarChannel.setTruthValue( (int)truthChannel );

  // flag current event as a MC Dalitz event
  // (needed for cut-flow histograms)
  m_isNonHyyStarHiggs = HG::eventIsNonHyyStarHiggs(all_particles);
  var::isNonHyyStarHiggs.setTruthValue(m_isNonHyyStarHiggs);

  return;
}

EL::StatusCode  HiggsGamGamStarCutflowAndMxAOD::doTruth()
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

  // Now done in HiggsGamGamStarCutflowAndMxAOD::SetTruthHiggsParticles()
  // HG::VarHandler::getInstance()->setHiggsBosons(&all_higgs);

  // Set the truth decay product containers in ExtraHggStarObjects
  // Now done in HiggsGamGamStarCutflowAndMxAOD::SetTruthHiggsParticles()
  //
  //const xAOD::TruthParticleContainer* all_particles = truthHandler()->getTruthParticles();
  //HG::ExtraHggStarObjects::getInstance()->setTruthHiggsDecayProducts(all_particles);

  // Adds event-level variables to TStore (this time using truth containers)
  bool truth = true;
  if (m_saveTruthVars) {

    writeNominalAndSystematicVars(truth);
    writeNominalOnlyVars(truth);
    writeTruthOnlyVars();
    if (m_saveDetailed)
    { writeDetailedVars(truth); }
    //
    const xAOD::TruthParticleContainer *higgsLeptons = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
    const xAOD::TruthParticleContainer *higgsPhotons = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsPhotons();

    bool isFiducial = true;
    int nFidLeptons=0;
    int nFidPhotons=0;

    if( higgsLeptons->size() !=2 || higgsPhotons->size()!=1 ){
      isFiducial = false;
    } else {
      for( const auto& part: *higgsLeptons){
        if( abs(part->pdgId()) == 11 ){
          if( part->pt() < 4.5 * HG::GeV )
            continue;
          if( fabs(part->eta() ) < 2.5 )
            continue;
        } else if ( abs(part->pdgId()) == 13 ){
          if( part->pt() < 3. * HG::GeV )
            continue;
          if( fabs(part->eta() ) < 2.7 )
            continue;
        } else {
          continue;
        }
        ++nFidLeptons;
      }

      for( const auto& part: *higgsPhotons){
        if( part->pt() < 30. * HG::GeV  )
          continue;
        if( fabs(part->eta() ) < 2.37 )
          continue;
        ++nFidPhotons;
      }
    }

    if( nFidPhotons !=1 || nFidLeptons !=2 ){
      isFiducial = false;
    }
    eventHandler()->storeTruthVar<char>("isFiducial", isFiducial);
  }


  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}


void HiggsGamGamStarCutflowAndMxAOD::fillCutFlow(CutEnum cut, double w) {

  getCutFlowHisto(false /*weighted*/, false /*onlyDalitz*/)->Fill(cut);

  if (HG::isData()) return;

  getCutFlowHisto(true, false /*onlyDalitz*/)->Fill(cut,w);

  // per-channel cutflows (mainly for debugging)
  HG::ChannelEnum truthChan = (HG::ChannelEnum)var::yyStarChannel.truth();

  // unweighted MC, per channel:
  if (truthChan == HG::DIMUON) getCutFlowHisto(false, false, HG::DIMUON)->Fill(cut);
  if (truthChan == HG::RESOLVED_DIELECTRON) getCutFlowHisto(false, false, HG::RESOLVED_DIELECTRON)->Fill(cut);
  if (truthChan == HG::MERGED_DIELECTRON) getCutFlowHisto(false, false, HG::MERGED_DIELECTRON)->Fill(cut);

  if (m_isNonHyyStarHiggs) return;

  getCutFlowHisto(true, true /*onlyDalitz*/)->Fill(cut,w);

  if (truthChan == HG::DIMUON) getCutFlowHisto(true, true , HG::DIMUON)->Fill(cut,w);
  if (truthChan == HG::RESOLVED_DIELECTRON) getCutFlowHisto(true, true , HG::RESOLVED_DIELECTRON)->Fill(cut,w);
  if (truthChan == HG::MERGED_DIELECTRON) getCutFlowHisto(true, true , HG::MERGED_DIELECTRON)->Fill(cut,w);

  return;
}


EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::finalize() {
  printf("\nEvent selection cut flow:\n");
  printCutFlowHistos();

  // Write the output to file
  HgammaAnalysis::finalize();

  SafeDelete(m_trackHandler);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::fileExecute() {
  // Things that you need to process for each individual input file, even those containing no events.

  HgammaAnalysis::fileExecute();

  // Signals to the code (in execute) that we have to book cutflow information.
  m_newFileMetaData = true;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsGamGamStarCutflowAndMxAOD::changeInput(bool firstFile) {
  // Here you do everything you need to do when we change input files.

  HgammaAnalysis::changeInput(firstFile);

  // Signals to the code (in execute) that we have to reset crossSectionBRfilterEff
  m_newFileLoaded = true;

  return EL::StatusCode::SUCCESS;
}

void HiggsGamGamStarCutflowAndMxAOD::decorateCorrectedIsoCut(xAOD::ElectronContainer & electrons, xAOD::MuonContainer & muons){

  if(var::yyStarChannel()==HG::DIMUON){
    std::vector<const xAOD::IParticle*> muonsVec;
    for(auto muon: muons) muonsVec.push_back((const xAOD::IParticle*) muon);
    for(auto muon: muons)
    {
      // Make the correction / decisions for every working point
      for (auto dec : m_muIsoAccCorr)
      {
        (*dec.second)(*muon) = m_isoCloseByTools_Muon[dec.first]->acceptCorrected(*muon, muonsVec);
        // this actually modifies isolation of individual objects (muons)
        m_isoCloseByTools_Muon[dec.first]->getCloseByIsoCorrection(nullptr, &muons);
      }
    }
  }
  else if(var::yyStarChannel()==HG::RESOLVED_DIELECTRON){
    std::vector<const xAOD::IParticle*> electronsVec;
    for(auto electron: electrons) electronsVec.push_back((const xAOD::IParticle*) electron);
    for(auto ele: electrons)
    {
      for (auto dec : m_eleIsoAccCorr)
      {
        // Make the correction / decisions for every working point
        (*dec.second)(*ele) = m_isoCloseByTools_Ele[dec.first]->acceptCorrected(*ele, electronsVec);
        //this actually modifies isolation of individual objects (electrons)
        m_isoCloseByTools_Ele[dec.first]->getCloseByIsoCorrection(&electrons);
      }
    }
  }
  //don't care about merged ele channel, since correction would not do anything there
}

HG::ChannelEnum HiggsGamGamStarCutflowAndMxAOD::FindZboson_ElectronChannelAware(xAOD::TrackParticleContainer* inTracks,
                                                                                xAOD::TrackParticle*& sel_trk1,
                                                                                xAOD::TrackParticle*& sel_trk2,
                                                                                double& return_mll,
                                                                                const HG::TrackElectronMap& trkEleMap,
                                                                                xAOD::ElectronContainer* inEleCont,
                                                                                xAOD::ElectronContainer* outEleCont
                                                                                )
{

  HG::ChannelEnum return_chan = HG::CHANNELUNKNOWN;
  double max_pt = -1;

  for (auto tracki : *inTracks) {
    for (auto trackj : *inTracks) {
      if (tracki == trackj) continue;
      if (tracki->pt() < trackj->pt()) continue;
      if (tracki->charge() == trackj->charge()) continue;

      xAOD::ElectronContainer tmp_eles(SG::VIEW_ELEMENTS);
      HG::ChannelEnum tmp_chan = HG::ClassifyElectronChannelsByBestMatch(tracki,trackj,trkEleMap,inEleCont,&tmp_eles);

      if (tmp_chan != HG::RESOLVED_DIELECTRON &&
          tmp_chan != HG::MERGED_DIELECTRON) continue;

      // apply preselection cuts
      double tmp_pt = (tracki->p4() + trackj->p4()).Pt();
      double tmp_m = -1;

      if (tmp_chan == HG::RESOLVED_DIELECTRON) {
        if (tmp_eles.size() != 2) HG::fatal("Z boson assignment error (resolved).");

        // Resolved preselection
        if (!m_eleIDPreselection.IsNull()) {
          if (!electronHandler()->passPIDCut(tmp_eles[0],m_eleIDPreselection)) continue;
          if (!electronHandler()->passPIDCut(tmp_eles[1],m_eleIDPreselection)) continue;
        }
        
        float lead_pt_cut = config()->getNum("ElectronHandler.Selection.PtLeadCutGeV", 13);
        if (tmp_eles[0]->pt() < lead_pt_cut*HG::GeV && tmp_eles[1]->pt() < lead_pt_cut*HG::GeV) continue;

        // In order to use standard electron ID, these tracks must be primary tracks!
        bool both_prim = (HG::MapHelpers::getMatchingTrackIndex(tmp_eles[0],tracki) == 0 &&
                          HG::MapHelpers::getMatchingTrackIndex(tmp_eles[1],trackj) == 0);
        both_prim = both_prim | (HG::MapHelpers::getMatchingTrackIndex(tmp_eles[1],tracki) == 0 &&
                                 HG::MapHelpers::getMatchingTrackIndex(tmp_eles[0],trackj) == 0);
        if (!both_prim) continue;

        TLorentzVector tmp_tlv = tmp_eles[0]->p4() + tmp_eles[1]->p4();
        tmp_m  = tmp_tlv.M();
      }
      else if (tmp_chan == HG::MERGED_DIELECTRON) {
        if (tmp_eles.size() != 1) HG::fatal("Z boson assignment error (merged).");

        // apply resolved preselection cuts
        if (!m_mergedElectronID->passPreselection(*tmp_eles[0],*tracki,*trackj)) continue;

        TLorentzVector merged = HG::MergedEleTLV(*tracki,*trackj,*tmp_eles[0]);
        tmp_m  = merged.M();
      }

      // Sort by max (vector sum) pt of the di-track system.
      // Other things tried:
      // - [prefer resolved channel always -> bad]
      // - [use calorimeter pt, with trk-pt tiebreaker -> ok, but not better.]
      if (tmp_pt > max_pt) {
        max_pt = tmp_pt;
        sel_trk1 = tracki;
        sel_trk2 = trackj;
        return_mll = tmp_m;
      }
    }
  }

  if (return_mll > 0)
  {
    // Re-do, but fill electron containers
    return_chan = HG::ClassifyElectronChannelsByBestMatch(sel_trk1,sel_trk2,trkEleMap,inEleCont,outEleCont);
  }

  return return_chan;
}

void HiggsGamGamStarCutflowAndMxAOD::AddMuonDecorations(xAOD::MuonContainer& muons) {

  //set corrected iso decision same as non-corrected by default
  for(auto muon: muons){
    for (auto dec : m_muIsoAccCorr){
      (*dec.second)(*muon) = muonHandler()->passIsoCut(muon,dec.first);
    }
  }

  return;
}

void HiggsGamGamStarCutflowAndMxAOD::AddElectronDecorations(xAOD::ElectronContainer& electrons) {

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

    //set corrected iso decision same as non-corrected by default
    for (auto dec : m_eleIsoAccCorr){
      (*dec.second)(*electron) = electronHandler()->passIsoCut(electron,dec.first);
    }

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
      photonHandler()->getCalibrationAndSmearingTool()->applyCorrection(*photon, *eventInfo());
      HG::EleAcc::calibratedPhotonEnergy(*electron) = photon->e();

      delete photon;
    } else {
      HG::EleAcc::calibratedPhotonEnergy(*electron) = -999;
    }

  }

    return;
}

void HiggsGamGamStarCutflowAndMxAOD::printCutFlowHistos() {
  for ( auto entry : m_cFlowHistos ) {

    TString isMC = HG::isMC()?"MC sample":"Data run";

    int cutFlowID = entry.first;

    TString typeOfEvents = "all events";
    if (cutFlowID > (int)1e9) {
      typeOfEvents = "only Dalitz events";
      cutFlowID -= (int)1e9;
    }

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
    for (auto chan : {HG::MERGED_DIELECTRON,HG::RESOLVED_DIELECTRON,HG::DIMUON}) {
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
