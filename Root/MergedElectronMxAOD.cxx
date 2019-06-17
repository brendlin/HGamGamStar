#include "HGamGamStar/MergedElectronMxAOD.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include <EventLoop/Worker.h>
#include "HGamAnalysisFramework/TruthUtils.h"
#include "xAODEgamma/EgammaxAODHelpers.h"

#include "HGamGamStar/ExtraHggStarObjects.h"
#include "HGamGamStar/TrackElectronMap.h"
#include "HGamGamStar/SimpleVertexFit.h"

#include <InDetTrackSelectionTool/IInDetTrackSelectionTool.h>
#include <TrackVertexAssociationTool/ITrackVertexAssociationTool.h>
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "xAODTracking/VertexAuxContainer.h"
#include "PhotonVertexSelection/PhotonVertexHelpers.h"


// #include "PhotonVertexSelection/PhotonPointingTool.h"
// #include "ZMassConstraint/ConstraintFit.h"

//template <class T >
std::pair<bool, int> findInVector(const DataVector<xAOD::IParticle>& vecOfElements,const xAOD::IParticle* element)
{
  std::pair<bool, int> result;

  // Find given element in vector
  result.first = false;
  result.second = -1;
  for(unsigned int i(0); i < vecOfElements.size();++i){
    if(HG::MapHelpers::getTheOriginalPointer(*vecOfElements[i]) == element){
      result.second = i;
      result.first = true;
    }
  }

  return result;
}


// this is needed to distribute the algorithm to the workers
ClassImp(MergedElectronMxAOD)

MergedElectronMxAOD::MergedElectronMxAOD(const char *name)
: MxAODTool(name)
  , m_trackHandler(nullptr)
  , m_trkselTool()
  , m_ttvaTool()
{

 //        m_trkselTool.declarePropertyFor(this, "TrackSelectionTool", "TrackSelectionTool to select tracks which made it actually into the isolation"); // Makes the track selection tool a settable property of this tool
 //         m_ttvaTool.declarePropertyFor(this, "TTVASelectionTool", "TTVASelectionTool to correct for the pile-up robust WPs");

}

MergedElectronMxAOD::~MergedElectronMxAOD() {}

EL::StatusCode MergedElectronMxAOD::initialize()
{
  HgammaAnalysis::initialize();

  if (!m_trkselTool.isUserConfigured()) {
    m_trkselTool.setTypeAndName("InDet::InDetTrackSelectionTool/TrackParticleSelectionTool");
    ATH_MSG_INFO("No TrackSelectionTool provided, so I will create and configure my own, called: " << m_trkselTool.name());
    // The z0 cut is checked in any case either by the
    // track to vertex association tool or by the tracking tool
    ANA_CHECK(m_trkselTool.setProperty("maxZ0SinTheta", 3.));
    // The minimum Pt requirement is lowered to 500 MeV because
    // the Loose ttva cone variables accept very low-pt tracks
    // https://gitlab.cern.ch/atlas/athena/blob/21.2/Reconstruction/RecoAlgs/IsolationAlgs/python/IsoUpdatedTrackCones.py#L21
    ANA_CHECK(m_trkselTool.setProperty("minPt", 500.));
    ANA_CHECK(m_trkselTool.setProperty("CutLevel", "Loose"));
  }
  if (!m_ttvaTool.isUserConfigured()){
    m_ttvaTool.setTypeAndName("CP::TrackVertexAssociationTool/ttva_selection_tool");
    ANA_CHECK(m_ttvaTool.setProperty("WorkingPoint", "Loose"));
  }
  ANA_CHECK(m_trkselTool.retrieve());
  ANA_CHECK(m_ttvaTool.retrieve())

  HG::ExtraHggStarObjects::getInstance()->setEventAndStore(event(), store());

  m_trackHandler = new HG::TrackHandler("TrackHandler", event(), store());
  ANA_CHECK(m_trackHandler->initialize(*config()));

  m_mergedElectronID = new HG::MergedElectronID();
  ANA_CHECK(m_mergedElectronID->initialize(*config()));

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

 // Resolved electron ID preselection. Applied in FindZboson_ElectronChannelAware.
  m_eleIDPreselection = config()->getStr("ResolvedElectrons.Preselection.PID","VeryLoose");

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode MergedElectronMxAOD::createOutput()
{
  // Read the output branch names - add option to make this configurable in future ?
  m_photonContainerName = "HGam"+config()->getStr("PhotonHandler.ContainerName");
  m_jetContainerName    = "HGam"+config()->getStr("JetHandler.ContainerName");
  m_elecContainerName   = "HGam"+config()->getStr("ElectronHandler.ContainerName");
  m_muonContainerName   = "HGam"+config()->getStr("MuonHandler.ContainerName");
  m_trackContainerName  = "HGam"+config()->getStr("TrackHandler.ContainerName");
  m_evtInfoName         = "EventInfo";
  m_truthEvtsName       = "TruthEvents";


  // Whether or not to apply systematic variations
  // If true, all variables will be stored
  m_applySystematics = config()->getBool("ApplySystematicVariations",false);

  // Whether to save objects (photons, jets ...)
  m_saveObjects = config()->getBool("SaveObjects",false);

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
  for (auto trig: trigs) extra.push_back(".passTrig_"+trig);

  declareOutputVariables(m_evtInfoName,"MxAOD.Variables.EventInfo", extra, ignore);

  // a.2 TruthEvents variables
  ignore = {};
  extra = {};
  declareOutputVariables(m_truthEvtsName,"MxAOD.Variables.TruthEvents", extra, ignore);

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

    declareOutputVariables(m_photonTruthContainerName  , "MxAOD.Variables.TruthPhotons"    );
    declareOutputVariables(m_elecTruthContainerName    , "MxAOD.Variables.TruthElectrons"  );
    declareOutputVariables(m_muonTruthContainerName    , "MxAOD.Variables.TruthMuons"      );
    declareOutputVariables(m_jetTruthContainerName     , "MxAOD.Variables.TruthJets"       );
    declareOutputVariables("HGam"+config()->getStr("TruthHandler.HiggsBosonContainerName"), "MxAOD.Variables.TruthHiggsBosons");
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MergedElectronMxAOD::execute()
{
  // Needed for all underlaying tools to be working corectly!
  HgammaAnalysis::execute();

  // Clear containers which point to objects from previous event
  HG::ExtraHggStarObjects::getInstance()->clearContainers();

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
    }

    m_newFileLoaded = false;
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

    HG::ExtraHggStarObjects::getInstance()->setTruthHiggsDecayProducts(truthParticles);
  }

  // Set this for every event, just in case.
  var::yyStarChannel.setValue(HG::CHANNELUNKNOWN);

  // apply cuts. Returned value will be the last passed cut
  m_cutFlow = cutflow();

  // fill the cut-flow histograms up to tight selection
  double wi = weightInitial();
  for (int cut=ALLEVTS;cut<m_cutFlow;++cut) {
    fillCutFlow(CutEnum(cut),wi);
  }

  // if desired, apply skimming.
  // That is, only write events that pass a given cut to the output
  if ( m_cutFlow < PASSALL)
    return EL::StatusCode::SUCCESS;

  // Selects the objects, does overlap removal, and calculate all
  // variables that will be saved
  doReco();

  // fill the cut-flow histograms after tight selection including SF weights


  // check if we should apply systematics or not
  if (m_applySystematics) {
    for (auto sys: getSystematics()) {
      // ignore nominal case, already done!
      if (sys.name() == "") continue;

      // apply the systmeatic variation and calculate the outupt
      CP_CHECK("MergedElectronMxAOD::execute()", applySystematicVariation(sys));
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
MergedElectronMxAOD::CutEnum MergedElectronMxAOD::cutflow()
{

  m_selElectrons = xAOD::ElectronContainer(SG::VIEW_ELEMENTS);


  //==== CUT 4 : Remove duplicate events (only for data) ====
  static bool checkDuplicates = config()->getBool("EventHandler.CheckDuplicates");
  if ( checkDuplicates && eventHandler()->isDuplicate() ) return DUPLICATE;

  //==== CUT 5 : GRL ====
  static bool requireGRL = config()->getBool("EventHandler.CheckGRL");
  if ( requireGRL && HG::isData() && !eventHandler()->passGRL(eventInfo()) ) return GRL;

  //==== CUT 7 : Detector quality ====
  if ( !(eventHandler()->passLAr (eventInfo()) &&
         eventHandler()->passTile(eventInfo()) &&
         eventHandler()->passSCT (eventInfo()) ) )
    return DQ;

  //==== CUT 8 : Require a vertex ====
  if ( !eventHandler()->passVertex(eventInfo()) ) return VERTEX;

  // Apply electron preselection.
  // HGamCore does not have an electron preselection step, so we make our own here:
  xAOD::ElectronContainer allElectrons = electronHandler()->getCorrectedContainer();
  AddElectronDecorations(allElectrons);




  for (auto electron : allElectrons) {
    if (!electronHandler()->passOQCut(electron)) { continue; }
    if (!electronHandler()->passPtEtaCuts(electron)) { continue; }
    if (!electronHandler()->passHVCut(electron)) { continue; }
    //Count the number of Si tracks matching the electron
    int nSiTrack(0);
    for( unsigned int trk_i(0); trk_i < electron->nTrackParticles(); ++trk_i){
      auto ele_tp =  electron->trackParticle(trk_i);
      if(!ele_tp){
        continue;
      }

      int nSiHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfSiliconHitsAndDeadSensors(ele_tp);
      if(nSiHitsPlusDeadSensors >= 7)
        ++nSiTrack;
    }
    //If 2 or more the electron is selected
    if(nSiTrack>1)
      m_selElectrons.push_back(electron);
  }


  decorateCorrectedIsoCut(m_selElectrons);

  //==== CUT 13 : Whether SF leptons survive OR
  if (m_selElectrons.size() == 0 )
    return ELECTRON;


  return PASSALL;


}

EL::StatusCode  MergedElectronMxAOD::doReco(bool /*isSys*/){
  // Do anything you missed in cutflow, and save the objects.


  // Adds event weights and catgory to TStore
  // Also sets pointer to photon container, etc., which is used by var's
  xAOD::PhotonContainer noPhotons = xAOD::PhotonContainer(SG::VIEW_ELEMENTS);
  xAOD::MuonContainer noMuons = xAOD::MuonContainer(SG::VIEW_ELEMENTS);

  setSelectedObjects(&noPhotons, &m_selElectrons, &noMuons, nullptr, nullptr, nullptr);

  xAOD::TrackParticleContainer noTracks = xAOD::TrackParticleContainer(SG::VIEW_ELEMENTS);
  HG::ExtraHggStarObjects::getInstance()->setElectronTrackContainer(&noTracks);

  // Adds event-level variables to TStore
  // Write in the nominal and systematics loops
  writeNominalAndSystematic();
  writeNominalAndSystematicVars();

  writeNominalOnly();
  writeNominalOnlyVars();


  if (m_saveDetailed) {
    writeDetailed();
    writeDetailedVars();
  }

  CP_CHECK("execute()", electronHandler()->writeContainer(m_selElectrons));


  // Adds all event variables (weight, category, isPassed, and pT_yy etc.)
  // to the TEvent output stream
  HG::VarHandler::getInstance()->write();


  return EL::StatusCode::SUCCESS;
}

void MergedElectronMxAOD::writeNominalAndSystematic()
{
  // Basic event selection flags
  var::cutFlow.setValue(m_cutFlow);

  if (HG::isMC()) {

    // Basic event weights
    eventHandler()->pileupWeight();
    eventHandler()->vertexWeight();
  }

  // Additional variables useful for non-framework analysis
  eventHandler()->storeVar<char>("isPassedEventSelection",m_cutFlow >= PASSALL);

}

void MergedElectronMxAOD::writeNominalAndSystematicVars(bool truth)
{
  // Put here all of the HGamVariables that you want to save in the nominal loop and
  // the systematics loops.
  // In the truth case, there is no "systematic" case, so they are saved only once.

  var::m_ll.addToStore(truth);
  var::deltaR_ll.addToStore(truth);
  var::pt_ll.addToStore(truth);

  if (!truth)
  {
    var::m_ll_track4mom.addToStore(false);
    var::Resolved_dRExtrapTrk12.addToStore(false);
    var::Resolved_deltaPhiRescaled2.addToStore(false);
    var::Resolved_deltaEta2.addToStore(false);
  }

  var::N_e    .addToStore(truth);
}


void MergedElectronMxAOD::writeTruthOnlyVars()
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


void MergedElectronMxAOD::writeNominalOnly()
{
  // Put here the things that you want to save only in the nominal loop (not the systematics loops).

  eventHandler()->mu();
  eventHandler()->runNumber();

  // Make sure every trigger is checked, and decorated to EventInfo
  eventHandler()->getPassedTriggers();

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
    truthHandler()->catCoup();
    eventHandler()->storeVar<float>("crossSectionBRfilterEff", m_crossSectionBRfilterEff);
  }

}

void MergedElectronMxAOD::writeNominalOnlyVars(bool /*truth*/)
{
  // Put here all of the HGamVariables that you want to save in the nominal loop only
  // (not the systematics loops).

}

void MergedElectronMxAOD::writeDetailed()
{
  // Put here all of the things that you want to save in the case where
  // "SaveDetailedVariables" is set to TRUE in the config file.

}

void MergedElectronMxAOD::writeDetailedVars(bool /*truth*/)
{
  // Put here all of the HGamVariables that you want to save in the case where
  // "SaveDetailedVariables" is set to TRUE in the config file.

}

EL::StatusCode  MergedElectronMxAOD::doTruth()
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
    truthHandler()->writeElectrons  (all_electrons);
    truthHandler()->writeHiggsBosons(all_higgs    );
    //truthHandler()->writePhotons    (all_photons  );
    //truthHandler()->writeMuons      (all_muons    );
    //truthHandler()->writeJets       (all_jets     );
    //truthHandler()->writeTruthEvents(             );
    addTruthLinks(m_elecContainerName.Data()  , m_elecTruthContainerName.Data());
  }

  HG::VarHandler::getInstance()->setTruthContainers(&all_photons, &electrons, &muons, &jets);
  HG::VarHandler::getInstance()->setHiggsBosons(&all_higgs);

  // Set the truth decay product containers in ExtraHggStarObjects
  //const xAOD::TruthParticleContainer* all_particles = truthHandler()->getTruthParticles();
  //HG::ExtraHggStarObjects::getInstance()->setTruthHiggsDecayProducts(all_particles);

  // Set the truth decay product containers in ExtraHggStarObjects
  var::yyStarChannel.setTruthValue( (int) HG::CHANNELUNKNOWN );
  var::yyStarChannel.setTruthValue( (int) truthClass() );


  // Adds event-level variables to TStore (this time using truth containers)
  bool truth = true;
  if (m_saveTruthVars) {

    writeNominalAndSystematicVars(truth);
    writeNominalOnlyVars(truth);
    writeTruthOnlyVars();
    if (m_saveDetailed)
    { writeDetailedVars(truth); }
    //
  }


  // Adds all event variables to the TEvent output stream
  HG::VarHandler::getInstance()->writeTruth();

  return EL::StatusCode::SUCCESS;
}


void MergedElectronMxAOD::fillCutFlow(CutEnum cut, double w) {
  getCutFlowHisto()->Fill(cut);
  if (HG::isData()) return;
  getCutFlowWeightedHisto()->Fill(cut,w);
  if (m_isNonHyyStarHiggs) return;
  getCutFlowHisto(true)->Fill(cut);
  getCutFlowWeightedHisto(true)->Fill(cut,w);
}


EL::StatusCode MergedElectronMxAOD::finalize() {
  printf("\nEvent selection cut flow:\n");
  printCutFlowHistos();

  // Write the output to file
  HgammaAnalysis::finalize();

  SafeDelete(m_trackHandler);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MergedElectronMxAOD::fileExecute() {
  // Things that you need to process for each individual input file, even those containing no events.

  HgammaAnalysis::fileExecute();

  // Signals to the code (in execute) that we have to book cutflow information.
  m_newFileMetaData = true;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MergedElectronMxAOD::changeInput(bool firstFile) {
  // Here you do everything you need to do when we change input files.

  HgammaAnalysis::changeInput(firstFile);

  // Signals to the code (in execute) that we have to reset crossSectionBRfilterEff
  m_newFileLoaded = true;

  return EL::StatusCode::SUCCESS;
}

void MergedElectronMxAOD::decorateCorrectedIsoCut(xAOD::ElectronContainer & electrons){

  //set corrected iso decision same as non-corrected by default
  for(auto electron: electrons){
    for (auto dec : m_eleIsoAccCorr){
      (*dec.second)(*electron) = electronHandler()->passIsoCut(electron,dec.first);
    }
  }
}




void MergedElectronMxAOD::AddElectronDecorations(xAOD::ElectronContainer& electrons) {
  const xAOD::VertexContainer* vertices = nullptr;
  if (event()->contains<xAOD::VertexContainer>("PrimaryVertices"))
    if (event()->retrieve(vertices,"PrimaryVertices").isFailure())
      HG::fatal("Error retrieving PrimaryVertices, exiting");

  const xAOD::Vertex *primaryVertex = xAOD::PVHelpers::getHardestVertex(vertices);

  xAOD::VertexContainer* outVerticies = new xAOD::VertexContainer();
  xAOD::VertexAuxContainer* outVerticiesAux = new xAOD::VertexAuxContainer();
  outVerticies->setStore(outVerticiesAux);
  eventHandler()->evtStore()->record( outVerticies, "TempVerticies" );
  eventHandler()->evtStore()->record( outVerticiesAux, "TempVerticiesAux." );

  xAOD::TrackParticleContainer all_tracks = trackHandler()->getCorrectedContainer();
  // Truth-track map
  HG::TruthTrackMap trkTruthMap = trackHandler()->MakeTruthTrackMapFromElectronContainer(electrons);
  // Track-electron map
  HG::TrackElectronMap trkEleMap;
  trackHandler()->findTracksFromElectrons(all_tracks,electrons,trkEleMap,true);
  // Get the track assoicated to the two leptons
  const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();

  std::vector<const xAOD::IParticle*> leptonTracks;
  for(const auto& lepton: *childleps){
    auto leptonTrack  =  HG::MapHelpers::getTrackMatchingTruth( lepton, trkTruthMap );
    if(leptonTrack)
      leptonTracks.push_back(leptonTrack);
  }

  for(const auto el: electrons){
    HG::EleAcc::truthTrackIndexA(*el) = -999;
    HG::EleAcc::truthTrackIndexB(*el) = -999;
    int nfound =0 ;
    for( unsigned int trk_i(0); trk_i < el->nTrackParticles(); ++trk_i)
    {

      auto found = std::find(leptonTracks.begin(),leptonTracks.end(),el->trackParticle(trk_i));
      if( found !=leptonTracks.end() )
      {
        if(nfound==0)
          HG::EleAcc::truthTrackIndexA(*el) = trk_i;
        else if (nfound==1)
          HG::EleAcc::truthTrackIndexB(*el) = trk_i;
        else
          std::cout << "Found too many tracks " <<  std::endl;
        ++nfound;
      }

    }

    float ambiR  = -999;
    int   ambiCT = -999;
    if( el->ambiguousObject() ){
      auto ambiPhoton = dynamic_cast<const xAOD::Photon*>( el->ambiguousObject() );
      if(ambiPhoton){
        ambiR = ambiPhoton->conversionRadius();
        ambiCT= ambiPhoton->conversionType();
      }
    }
    HG::EleAcc::ambiguousPhotonR(*el) = ambiR;
    HG::EleAcc::ambiguousPhotonCT(*el) = ambiCT;



    xAOD::Photon* photon = createPhotonFromElectron(el, outVerticies);


    if(photon)
    {
      //std::cout << "Photon created " <<  std::endl;
      photonHandler()->getCalibrationAndSmearingTool()->applyCorrection(*photon, *eventInfo());
      //std::cout << "Photon calibrated " <<  std::endl;

      HG::EleAcc::calibratedPhotonEnergy(*el) = photon->e();
    } else {
      HG::EleAcc::calibratedPhotonEnergy(*el) = -999;
    }
    delete photon;
    HG::EleAcc::calibratedElectronEnergy(*el) = el->e();

    //std::cout << "E " <<  HG::EleAcc::calibratedElectronEnergy(*el)  << " ph " << HG::EleAcc::calibratedPhotonEnergy(*el) << std::endl;

  }
  //delete outVerticies;

  for (auto electron : electrons) {
    std::vector<int>   passTTVA;
    std::vector<float> trackPT;
    std::vector<float> trackD0;
    std::vector<float> trackZ0;
    std::vector<float> trackP;
    std::vector<float> trackD0Sig;
    std::vector<float> trackZ0Sig;
    std::vector<float> trackTRT_PID_trans;
    std::vector<int>   trackNPix;
    std::vector<int>   trackNSCT;
    std::vector<int>   trackPassBL;
    std::vector<int>   trackNIBL;
    std::vector<int>   trackNBL;
    std::vector<int>   trackSharedIBL;
    std::vector<int>   trackSharedBL;
    std::vector<int>   trackSplitIBL;
    std::vector<int>   trackSplitBL;
    std::vector<int>   trackPdgID;
    std::vector<int>   trackBarcode;
    std::vector<float> trackTruthE;
    std::vector<int>   trackFromHiggs;
    int   isTrueMergedE = 0;
    float trueEnergy = -999;
    float trueMass = -999;
    float trueEta = -999;
    float truePhi = -999;


    int nFromHiggs(0);
    TLorentzVector sumOfTruthProducts;

    std::vector <const xAOD::TrackParticle*> trksToFit;
    std::vector <int>  indexToFit;

    for( unsigned int trk_i(0); trk_i < electron->nTrackParticles(); ++trk_i){
      auto ele_tp =  electron->trackParticle(trk_i);

      passTTVA.push_back(0);
      trackPT.push_back(-999);
      trackD0.push_back(-999);
      trackZ0.push_back(-999);
      trackP.push_back(-999);
      trackD0Sig.push_back(-999);
      trackZ0Sig.push_back(-999);
      trackTRT_PID_trans.push_back(-999);
      trackNPix.push_back(-999);
      trackNSCT.push_back(-999);
      trackPassBL.push_back(-999);
      trackNBL.push_back(-999);
      trackNIBL.push_back(-999);
      trackSharedIBL.push_back(-999);
      trackSharedBL.push_back(-999);
      trackSplitIBL.push_back(-999);
      trackSplitBL.push_back(-999);
      trackPdgID.push_back(-999);
      trackBarcode.push_back(-999);
      trackTruthE.push_back(-999);
      trackFromHiggs.push_back(0);
      if(!ele_tp)
        continue;

      //auto cov = ele_tp->definingParametersCovMatrix();
      //std::cout << cov << std::endl;

      trackPT.back() = ele_tp->pt();
      trackD0.back() = ele_tp->d0();
      trackZ0.back() = ele_tp->z0();
      trackP.back() = 1./(ele_tp->qOverP());
      trackD0Sig.back() = xAOD::TrackingHelpers::d0significance(ele_tp);
      trackZ0Sig.back() = xAOD::TrackingHelpers::z0significance(ele_tp);
      trackTRT_PID_trans.back() = trackHandler()->calculateTRT_PID(*ele_tp);

      uint8_t shared;
      if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerSharedHits) )
        trackSharedIBL.back() = (int)shared;
      if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerSharedHits) )
        trackSharedBL.back() = (int)shared;
      if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerSplitHits ) )
        trackSplitIBL.back() = (int)shared;
      if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerSplitHits) )
        trackSplitBL.back() = (int)shared;

      if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerHits) )
        trackNIBL.back() = (int)shared;
      if( trackNIBL.back() == 0 && ele_tp->summaryValue(shared, xAOD::expectInnermostPixelLayerHit) )
        trackNIBL.back() = -(int)shared;

      if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerHits) )
        trackNBL.back() = (int)shared;
      if( trackNBL.back() == 0 && ele_tp->summaryValue(shared, xAOD::expectNextToInnermostPixelLayerHit) )
        trackNBL.back() = -(int)shared;


      trackNPix.back() = ElectronSelectorHelpers::numberOfPixelHitsAndDeadSensors(ele_tp);
      trackNSCT.back() = ElectronSelectorHelpers::numberOfSCTHitsAndDeadSensors(ele_tp);
      trackPassBL.back() = ElectronSelectorHelpers::passBLayerRequirement(ele_tp);
      const xAOD::TruthParticle* truthPart = xAOD::TruthHelpers::getTruthParticle(*ele_tp);
      if( trackNPix.back() + trackNSCT.back()  > 7 ){
        if(trksToFit.size() == 0 )
        {
          trksToFit.push_back( ele_tp );
          indexToFit.push_back( trk_i );
        }else {
          if(trksToFit.size() == 1 && trksToFit[0]->charge() !=ele_tp->charge()){
            trksToFit.push_back( ele_tp );
            indexToFit.push_back( trk_i );
          }
        }
      }

      if(truthPart){
        trackPdgID.back()   =  truthPart->pdgId();
        trackBarcode.back() =  truthPart->barcode();
        trackTruthE.back()  =  truthPart->p4().E();
        if(HG::isFromHiggs(truthPart))
        {
          trackFromHiggs.back() = 1;
          sumOfTruthProducts += truthPart->p4();
          ++nFromHiggs;
        }
      }
      if(nFromHiggs==2){
        isTrueMergedE=1;
        trueEnergy = sumOfTruthProducts.E();
        trueMass = sumOfTruthProducts.M();
        trueEta = sumOfTruthProducts.Eta();
        truePhi = sumOfTruthProducts.Phi();
      }


      if( !m_trkselTool->accept(*ele_tp, primaryVertex) )
        continue;
      if( !m_ttvaTool->isCompatible(*ele_tp,*primaryVertex) )
        continue;
      passTTVA.back()=1;
    }

    if( trksToFit.size() == 2 ){
      SimpleVertexFit svf;
      AmgVector(3) bs( -0.5,
                     -0.5,
                      0  );
      Vertex vtx = svf.fitVertex( trksToFit, bs );

      HG::EleAcc::standAloneVertexR(*electron) = vtx.position().perp();
      HG::EleAcc::standAloneIndexA(*electron) = indexToFit[0];
      HG::EleAcc::standAloneIndexB(*electron) = indexToFit[1];

    } else {
      HG::EleAcc::standAloneVertexR(*electron) = -999;
      HG::EleAcc::standAloneIndexA(*electron) = -999;
      HG::EleAcc::standAloneIndexB(*electron) = -999;
    }

    HG::EleAcc::passTTVA(*electron)      = passTTVA;
    HG::EleAcc::trackPT(*electron)       = trackPT;
    HG::EleAcc::trackD0(*electron)       = trackD0;
    HG::EleAcc::trackZ0(*electron)       = trackZ0;
    HG::EleAcc::trackP(*electron)        = trackP;
    HG::EleAcc::trackD0Sig(*electron)    = trackD0Sig;
    HG::EleAcc::trackZ0Sig(*electron)    = trackZ0Sig;
    HG::EleAcc::trackTRT_PID_trans(*electron) = trackTRT_PID_trans;
    HG::EleAcc::trackNPix(*electron)     = trackNPix;
    HG::EleAcc::trackNSCT(*electron)     = trackNSCT;
    HG::EleAcc::trackPassBL(*electron)   = trackPassBL;
    HG::EleAcc::trackNBL(*electron)      = trackNBL;
    HG::EleAcc::trackNIBL(*electron)     = trackNIBL;
    HG::EleAcc::trackSharedBL(*electron)  = trackSharedBL;
    HG::EleAcc::trackSharedIBL(*electron) = trackSharedIBL;
    HG::EleAcc::trackSplitBL(*electron)   = trackSplitBL;
    HG::EleAcc::trackSplitIBL(*electron)  = trackSplitIBL;

    HG::EleAcc::isTrueMergedE(*electron) = isTrueMergedE;
    HG::EleAcc::trueEnergy(*electron)    = trueEnergy;
    HG::EleAcc::trueMass(*electron)      = trueMass;
    HG::EleAcc::trueEta(*electron)       = trueEta;
    HG::EleAcc::truePhi(*electron)       = truePhi;
    HG::EleAcc::trackPdgID(*electron)    = trackPdgID;
    HG::EleAcc::trackBarcode(*electron)  = trackBarcode;
    HG::EleAcc::trackFromHiggs(*electron)= trackFromHiggs;


    if( HG::EleAcc::vtxTrkIndex1.isAvailable(*electron) && HG::EleAcc::vtxTrkIndex2.isAvailable(*electron)  ){
      int index = HG::EleAcc::vtxTrkIndex1(*electron);
      bool hasIndex = index >= 0;
      HG::EleAcc::vtxTrk1_PT (*electron)= hasIndex ? trackPT[index] : -999;
      HG::EleAcc::vtxTrk1_P (*electron) = hasIndex ? trackP[index] : -999;
      HG::EleAcc::vtxTrk1_D0 (*electron) = hasIndex ? trackD0[index] : -999;
      HG::EleAcc::vtxTrk1_D0Sig (*electron) = hasIndex ? trackD0Sig[index] : -999;
      HG::EleAcc::vtxTrk1_Z0 (*electron)  = hasIndex ? trackZ0[index] : -999;
      HG::EleAcc::vtxTrk1_Z0Sig (*electron) = hasIndex ? trackD0Sig[index] : -999;
      HG::EleAcc::vtxTrk1_TRT_PID_trans(*electron) = hasIndex ? trackTRT_PID_trans[index] : -999;
      HG::EleAcc::vtxTrk1_NPix(*electron) = hasIndex ? trackNPix[index] : -999;
      HG::EleAcc::vtxTrk1_NSCT(*electron) = hasIndex ? trackNSCT[index] : -999;
      HG::EleAcc::vtxTrk1_PassBL(*electron) = hasIndex ? trackPassBL[index] : -999;
      HG::EleAcc::vtxTrk1_NBL(*electron) = hasIndex ? trackNBL[index] : -999;
      HG::EleAcc::vtxTrk1_NIBL(*electron) = hasIndex ? trackNIBL[index] : -999;
      HG::EleAcc::vtxTrk1_SplitBL(*electron) = hasIndex ? trackSplitBL[index] : -999;
      HG::EleAcc::vtxTrk1_SharedBL(*electron) = hasIndex ? trackSharedBL[index] : -999;
      HG::EleAcc::vtxTrk1_SplitIBL(*electron) = hasIndex ? trackSplitIBL[index] : -999;
      HG::EleAcc::vtxTrk1_SharedIBL(*electron) = hasIndex ? trackSharedIBL[index] : -999;
      HG::EleAcc::vtxTrk1_PdgID(*electron) = hasIndex ? trackPdgID[index] : -999;
      HG::EleAcc::vtxTrk1_Barcode(*electron) = hasIndex ? trackBarcode[index] : -999;
      HG::EleAcc::vtxTrk1_TruthE(*electron) = hasIndex ? trackTruthE[index] : -999;
      HG::EleAcc::vtxTrk1_FromHiggs(*electron) = hasIndex ? trackFromHiggs[index] : -999;
      HG::EleAcc::vtxTrk1_dEta2_P(*electron) = hasIndex ? HG::EleAcc::TrackMatchingP_dEta2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dEta1_P(*electron) = hasIndex ? HG::EleAcc::TrackMatchingP_dEta1(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dPhi2_P(*electron) = hasIndex ? HG::EleAcc::TrackMatchingP_dPhi2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dEta2_LM(*electron) = hasIndex ? HG::EleAcc::TrackMatchingLM_dEta2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dEta1_LM(*electron) = hasIndex ? HG::EleAcc::TrackMatchingLM_dEta1(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dPhi2_LM(*electron) = hasIndex ? HG::EleAcc::TrackMatchingLM_dPhi2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dEta2_T(*electron) = hasIndex ? HG::EleAcc::TrackMatchingTrue_dEta2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk1_dPhi2_T(*electron) = hasIndex ? HG::EleAcc::TrackMatchingTrue_dPhi2(*electron)[index] : -999;

      index = HG::EleAcc::vtxTrkIndex2(*electron);
      hasIndex = index >= 0;
      HG::EleAcc::vtxTrk2_PT (*electron)= hasIndex ? trackPT[index] : -999;
      HG::EleAcc::vtxTrk2_P (*electron) = hasIndex ? trackP[index] : -999;
      HG::EleAcc::vtxTrk2_D0 (*electron) = hasIndex ? trackD0[index] : -999;
      HG::EleAcc::vtxTrk2_D0Sig (*electron) = hasIndex ? trackD0Sig[index] : -999;
      HG::EleAcc::vtxTrk2_Z0 (*electron)  = hasIndex ? trackZ0[index] : -999;
      HG::EleAcc::vtxTrk2_Z0Sig (*electron) = hasIndex ? trackD0Sig[index] : -999;
      HG::EleAcc::vtxTrk2_TRT_PID_trans(*electron) = hasIndex ? trackTRT_PID_trans[index] : -999;
      HG::EleAcc::vtxTrk2_NPix(*electron) = hasIndex ? trackNPix[index] : -999;
      HG::EleAcc::vtxTrk2_NSCT(*electron) = hasIndex ? trackNSCT[index] : -999;
      HG::EleAcc::vtxTrk2_PassBL(*electron) = hasIndex ? trackPassBL[index] : -999;
      HG::EleAcc::vtxTrk2_NBL(*electron) = hasIndex ? trackNBL[index] : -999;
      HG::EleAcc::vtxTrk2_NIBL(*electron) = hasIndex ? trackNIBL[index] : -999;
      HG::EleAcc::vtxTrk2_SplitBL(*electron) = hasIndex ? trackSplitBL[index] : -999;
      HG::EleAcc::vtxTrk2_SharedBL(*electron) = hasIndex ? trackSharedBL[index] : -999;
      HG::EleAcc::vtxTrk2_SplitIBL(*electron) = hasIndex ? trackSplitIBL[index] : -999;
      HG::EleAcc::vtxTrk2_SharedIBL(*electron) = hasIndex ? trackSharedIBL[index] : -999;
      HG::EleAcc::vtxTrk2_PdgID(*electron) = hasIndex ? trackPdgID[index] : -999;
      HG::EleAcc::vtxTrk2_Barcode(*electron) = hasIndex ? trackBarcode[index] : -999;
      HG::EleAcc::vtxTrk2_TruthE(*electron) = hasIndex ? trackTruthE[index] : -999;
      HG::EleAcc::vtxTrk2_FromHiggs(*electron) = hasIndex ? trackFromHiggs[index] : -999;
      HG::EleAcc::vtxTrk2_dEta2_P(*electron) = hasIndex ? HG::EleAcc::TrackMatchingP_dEta2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dEta1_P(*electron) = hasIndex ? HG::EleAcc::TrackMatchingP_dEta1(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dPhi2_P(*electron) = hasIndex ? HG::EleAcc::TrackMatchingP_dPhi2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dEta2_LM(*electron) = hasIndex ? HG::EleAcc::TrackMatchingLM_dEta2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dEta1_LM(*electron) = hasIndex ? HG::EleAcc::TrackMatchingLM_dEta1(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dPhi2_LM(*electron) = hasIndex ? HG::EleAcc::TrackMatchingLM_dPhi2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dEta2_T(*electron) = hasIndex ? HG::EleAcc::TrackMatchingTrue_dEta2(*electron)[index] : -999;
      HG::EleAcc::vtxTrk2_dPhi2_T(*electron) = hasIndex ? HG::EleAcc::TrackMatchingTrue_dPhi2(*electron)[index] : -999;
    }


    HG::EleAcc::trueType(*electron) = truthType( electron );



    // NEED TO initialize merged electron ID variables here!
    HG::EleAcc::EOverP0P1(*electron)               = -999;
    HG::EleAcc::dRExtrapTrk12(*electron)           = -999;
    HG::EleAcc::dRExtrapTrk12_LM(*electron)        = -999;
    HG::EleAcc::delta_z0sinTheta_tracks(*electron) = -999;
    HG::EleAcc::delta_z0_tracks(*electron)         = -999;
    // HG::EleAcc::dRbetweenTracks_LM_L1(*electron) = -999;
    // HG::EleAcc::dRbetweenTracks_LM_L2(*electron) = -999;
    // HG::EleAcc::dRbetweenTracks_P_L1(*electron) = -999;
    // HG::EleAcc::dRbetweenTracks_P_L2(*electron) = -999;

    // Decorate Rhad
    double feta = fabs(electron->eta());
    HG::EleAcc::RhadForPID(*electron) = (0.8 < feta && feta < 1.37) ?
      electron->showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rhad) :
      electron->showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rhad1);

  }

  return;
}


MergedElectronMxAOD::ElectronTruthType MergedElectronMxAOD::truthType( const xAOD::Electron* el ) const
{
  size_t nTracks = HG::EleAcc::trackPassBL(*el).size();
  int nEl(0);
  int nTot(0);
  int nHiggs(0);
  for( size_t trk(0); trk <  nTracks; ++trk){
    if( HG::EleAcc::trackFromHiggs(*el)[trk] > 0 )
    {
      if(nTot<2)
        ++nHiggs;
    }
    if(HG::EleAcc::trackPassBL(*el)[trk]>0)
    {
      if(abs(HG::EleAcc::trackPdgID(*el)[trk]) == 11)
      {
        ++nEl;
      }
      ++nTot;
    }
  }

  if(HG::EleAcc::isTrueMergedE(*el) )
  {

    int truthIndexA = HG::EleAcc::truthTrackIndexA(*el);
    int truthIndexB = HG::EleAcc::truthTrackIndexB(*el);

    if( truthIndexA != -999 &&  truthIndexB != -999 &&
       abs(HG::EleAcc::trackPdgID(*el)[truthIndexA]) == 11 &&
       abs(HG::EleAcc::trackPdgID(*el)[truthIndexB]) == 11 &&
       abs(HG::EleAcc::trackBarcode(*el)[truthIndexA]) < 200000 &&
       abs(HG::EleAcc::trackBarcode(*el)[truthIndexB]) < 200000 )
    {
      return SignalGood;
    }
    return SignalCompromised;
  }

  if(nTot == 1 && nEl==0)
    return BackgroundHad;

  if(nEl == 1)
    return BackgroundElHad;

  if(nEl >= 2)
    return BackgroundElEl;

  return BackgroundHadHad;
}



HG::ChannelEnum MergedElectronMxAOD::truthClass()
{
  // Get Higgs Final State Decay Products

  if (!HG::isMC())
    return HG::CHANNELUNKNOWN;

  // Get Leptons from Higgs decay products
  const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
  if (childleps->size() != 2)
    return HG::OTHER;

  // Check for out-of-acceptance
  for(const auto& lepton: *childleps){
    if (fabs(lepton->pdgId()) == 11) {
      if (lepton->pt()/1000. < 0.3) return HG::OUT_OF_ACCEPTANCE;
      if (fabs(lepton->eta()) > 2.5) return HG::OUT_OF_ACCEPTANCE;
    }
    else if (fabs(lepton->pdgId()) == 13) {
      if (lepton->pt()/1000. < 3.0) return HG::OUT_OF_ACCEPTANCE;
      if (fabs(lepton->eta()) > 2.7) return HG::OUT_OF_ACCEPTANCE;
    }
  }

  // Check if there are electrons in the decay
  bool isElectron = true;
  bool isMuon =  false;
  for(const auto& lepton: *childleps){
    if( fabs(lepton->pdgId()) != 11 )
      isElectron =  false;
    if( fabs(lepton->pdgId()) == 13 )
      isMuon  = true;
  }

  if(isMuon)
    return HG::DIMUON;

  if(!isElectron)
    return HG::OTHER;

  // Fill maps linking truth particle and track and track and electron
  // To be used to work out how many times a truth particle matches an electron candidate
  xAOD::ElectronContainer all_elecs = electronHandler()->getCorrectedContainer();
  xAOD::TrackParticleContainer all_tracks = trackHandler()->getCorrectedContainer();

  // Truth-track map
  HG::TruthTrackMap trkTruthMap = trackHandler()->MakeTruthTrackMapFromElectronContainer(all_elecs);

  // Track-electron map
  HG::TrackElectronMap trkEleMap;
  trackHandler()->findTracksFromElectrons(all_tracks,all_elecs,trkEleMap,true);

  // Get the track assoicated to the two leptons
  std::vector<const xAOD::TrackParticle*> leptonTracks;
  for(const auto& lepton: *childleps){
    auto leptonTrack  =  HG::MapHelpers::getTrackMatchingTruth( lepton, trkTruthMap );
    if(leptonTrack)
      leptonTracks.push_back(leptonTrack);
  }

  // If the the two tracks are not reconstructed then the truth matching has failed, exit
  if(leptonTracks.size()!=2)
    return HG::FAILEDTRKELECTRON;


  SimpleVertexFit svf;
  AmgVector(3) bs( -0.5,
                   -0.5,
                   0  );
  Vertex vtx = svf.fitVertex(leptonTracks, bs );
  var::vertexTruthFitRadius.setValue( vtx.position().perp() );
  var::vertexTruthFitRadius.setTruthValue( vtx.position().perp() );


  // Electron disambiguation is continued below in (reco-only) function
  return ClassifyElectronChannelsByBestMatch(leptonTracks[0],leptonTracks[1],trkEleMap);
}


HG::ChannelEnum MergedElectronMxAOD::ClassifyElectronChannelsByBestMatch(const xAOD::TrackParticle* trk0,
                                                                        const xAOD::TrackParticle* trk1,
                                                                        const HG::TrackElectronMap& trkEleMap)
{

  // Find the electrons associated to the tracks
  auto Trk0_Electrons = HG::MapHelpers::getElectronsMatchingTrack( trk0, trkEleMap );
  auto Trk1_Electrons = HG::MapHelpers::getElectronsMatchingTrack( trk1, trkEleMap );

  // Count the number of electrons a track matches to
  int Trk0_nElectron = Trk0_Electrons.size();
  int Trk1_nElectron = Trk1_Electrons.size();

  // If each track only matches to 1 electron each then it is quite simple
  if( Trk0_nElectron == 1 &&  Trk1_nElectron == 1){
    // Match the same electron --  Merged
    if( Trk0_Electrons.front() == Trk1_Electrons.front() )
    {
      return HG::MERGED_DIELECTRON;
    }
    else
    {
      return HG::RESOLVED_DIELECTRON;
    }
  }

  // Tracks match more than one electron each - lets see if they are the best match for any electron
  // Determine the match ranking for the track to each electron
  std::vector<int> Trk0_TrackNo = HG::MapHelpers::getMatchingTrackIndices(Trk0_Electrons,trk0);

  // Check if the track is ever the primary track
  int Trk0_PrimaryE(-1);
  for( unsigned int i(0); i < Trk0_TrackNo.size(); ++i){
    if( Trk0_TrackNo[i] == 0 ){
      // If the track is the primary track for multiple electrons
      // choose the one with higher pT
      if( Trk0_PrimaryE > -1 ){
        if( Trk0_Electrons[i]->pt() > Trk0_Electrons[Trk0_PrimaryE]->pt() )
          Trk0_PrimaryE = i;
      }else{
        Trk0_PrimaryE = i;
      }
    }
  }

  // Same again for the other electron
  std::vector<int> Trk1_TrackNo = HG::MapHelpers::getMatchingTrackIndices(Trk1_Electrons,trk1);

  int Trk1_PrimaryE(-1);
  for( unsigned int i(0); i < Trk1_TrackNo.size(); ++i){
    if( Trk1_TrackNo[i] == 0 ){
      // If the track is the primary track for multiple electrons
      // choose the one with higher pT
      if( Trk1_PrimaryE > -1 ){
        if( Trk1_Electrons[i]->pt() > Trk1_Electrons[Trk1_PrimaryE]->pt() )
          Trk1_PrimaryE = i;
      }else{
        Trk1_PrimaryE = i;
      }
    }
  }

  // If both tracks are the primary track for an electron the it is resolved
  if( Trk0_PrimaryE > -1 && Trk1_PrimaryE > -1 ){
    // Get to the correct pair iterator
    // This is a pair of TrackParticle , Electron
    auto el0 = Trk0_Electrons[Trk0_PrimaryE];
    auto el1 = Trk1_Electrons[Trk1_PrimaryE];

    // Compare if the electrons are the same
    if( el0 == el1 )
    {
      HG::fatal("Electron classification truth error!");
      return HG::AMBIGUOUS_DIELECTRON;
    }
    return HG::RESOLVED_DIELECTRON;
  }

  // If either are primary
  if( Trk0_PrimaryE > -1 || Trk1_PrimaryE > -1 ){
    const xAOD::Electron*  el = nullptr;
    //const xAOD::Electron*  elOther = nullptr;
    const xAOD::TrackParticle* otherTrack = nullptr;
    int nEleOther = 0;
    if( Trk0_PrimaryE >= 0 ){
      el = Trk0_Electrons[Trk0_PrimaryE];
      otherTrack = trk1;
      nEleOther = Trk1_nElectron;
    }else{
      el = Trk1_Electrons[Trk1_PrimaryE];
      otherTrack = trk0;
      nEleOther = Trk0_nElectron;
    }

    // Search for the other track in the electron
    for( unsigned int trk_i(0); trk_i < el->nTrackParticles(); ++trk_i){
      if( el->trackParticle(trk_i) == otherTrack )
      {
        return HG::MERGED_DIELECTRON;
      }
    }
    //If the other track is only matched to one electron and its not the primary track
    //We assume the reco has made a mistake so we will call it resolved
    if(nEleOther == 1)
    {
      return HG::RESOLVED_DIELECTRON;
    }
  }

  //Tracks are not primary for any electron candidate, tracks match to mutiple candiates
  // --  generally the reco has made a mess.
  // Might want to keep looking for candidates


  return HG::AMBIGUOUS_DIELECTRON;
}



xAOD::Photon*  MergedElectronMxAOD::createPhotonFromElectron (const xAOD::Electron* el, xAOD::VertexContainer* vertexContainer) const
{


  int index1 = -999;
  int index2 = -999;

  if( HG::EleAcc::vtxTrkIndex1.isAvailable(*el) && HG::EleAcc::vtxTrkIndex2.isAvailable(*el)  ){
    index1 = HG::EleAcc::vtxTrkIndex1(*el);
    index2 = HG::EleAcc::vtxTrkIndex2(*el);
  }
  if(index1 < 0 || index2 < 0 ){
    return 0;
  }
  //std::cout << "Index 1/2  " << index1  << " " << index2 << std::endl;

  xAOD::Photon* photon = new xAOD::Photon();
  photon->makePrivateStore();

  if( el->ambiguousObject() ){
    //std::cout << "Copying photon" <<  std::endl;
    auto ambiPhoton = dynamic_cast<const xAOD::Photon*>( el->ambiguousObject() );
    photon->Photon_v1::operator=(*ambiPhoton);
    photon->setCaloClusterLinks(el->caloClusterLinks());
  } else {
    //std::cout << "Creating photon" <<  std::endl;
    photon = new xAOD::Photon();
    photon->Egamma_v1::operator=(*el);
    photon->setCaloClusterLinks(el->caloClusterLinks());
  }

  //std::cout << "Setting Topo" <<  std::endl;

  static SG::AuxElement::Decorator<int> nClu("numTopoClusters") ;
  auto clusterVec = xAOD::EgammaHelpers::getAssociatedTopoClusters( photon->caloCluster() );
  nClu(*el)= clusterVec.size();


  static SG::AuxElement::Accessor<float> vtxPhi("vtxPhi") ;
  static SG::AuxElement::Accessor<float> vtxZ("vtxZ") ;

  float vtxR = 20; //
  float vtxX = vtxR * cos(vtxPhi(*el));
  float vtxY = vtxR * sin(vtxPhi(*el));

  //std::cout << "Creating Vtx " <<  std::endl;

  xAOD::Vertex*  vertex = new xAOD::Vertex();
  vertexContainer->push_back(vertex);

  vertex->setZ(vtxZ(*el));
  vertex->setX(vtxX);
  vertex->setY(vtxY);

  //std::cout << "Decorating Vtx " <<  std::endl;

  // decorate with pt1, pt2
  vertex->auxdata<float>("pt1") = el->trackParticle(index1)->pt();
  vertex->auxdata<float>("pt2") = el->trackParticle(index2)->pt();


  //std::cout << "Setting Trk element links Vtx " <<  std::endl;
  //links to tracks
  std::vector<ElementLink<xAOD::TrackParticleContainer>> links_tracks;
  links_tracks.push_back( el->trackParticleLinks()[index1] );
  links_tracks.push_back( el->trackParticleLinks()[index2] );

  //set vertex - track links
  //vertex->setTrackParticleLinks(links_tracks);

  //std::cout << "Setting Vtx element links Photon " <<  std::endl;

  std::vector<ElementLink<xAOD::VertexContainer>> links_verticies;
  //std::cout << "-- Setting Vtx element links Photon " <<  std::endl;
  links_verticies.push_back(ElementLink<xAOD::VertexContainer>(vertex, *vertexContainer));
  //std::cout << "---- Setting Vtx element links Photon " <<  std::endl;
  photon->setVertexLinks(links_verticies);

  //std::cout << "Setting deltas " <<  std::endl;

  static SG::AuxElement::Accessor<float> vtxdEta("vtxdEta") ;
  static SG::AuxElement::Accessor<float> vtxdPhi("vtxdPhi") ;


  float dEta = vtxdEta(*el);
  float dPhi = vtxdPhi(*el);
  photon->setVertexCaloMatchValue( dEta, xAOD::EgammaParameters::convMatchDeltaEta1 );
  photon->setVertexCaloMatchValue( dEta, xAOD::EgammaParameters::convMatchDeltaEta2 );
  photon->setVertexCaloMatchValue( dPhi, xAOD::EgammaParameters::convMatchDeltaPhi1 );
  photon->setVertexCaloMatchValue( dPhi, xAOD::EgammaParameters::convMatchDeltaPhi2 );


  return photon;
}
