#ifndef HGamGamStar_HiggsGamGamStarCutflowAndMxAOD_H
#define HGamGamStar_HiggsGamGamStarCutflowAndMxAOD_H

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamGamStar/MxAODTool.h"
/* #include "FsrUtils/FsrPhotonTool.h" */
#include "AsgTools/ToolHandle.h"
/* #include "ZMassConstraint/IConstraintFit.h" */
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/TrackHandler.h"
#include "HGamGamStar/MergedElectronID.h"
#include "HGamGamStar/MergedElectronID_v2.h"

#include "IsolationSelection/IsolationCloseByCorrectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"

class HiggsGamGamStarCutflowAndMxAOD : public MxAODTool
{

private:

  // Cut-flow - need to keep the same order!
  enum CutEnum {
    // Covered in MxAODTool: NxAOD=0, NDxAOD=1, ALLEVTS=2
    HIGGS_LEP_DALITZ=3, DUPLICATE=4, GRL=5, TRIGGER=6, DQ=7, VERTEX=8,
    TWO_SF_LEPTONS=9,ONE_LOOSE_GAM=10, AMBIGUITY=11,
    ZBOSON_ASSIGNMENT=12,TWO_SF_LEPTONS_POSTOR=13,BAD_MUON=14,ONE_PHOTON_POSTOR=15,
    TRIG_MATCH=16, LEP_MEDID=17, LEP_IP=18, LEP_ISO=19,
    GAM_TIGHTID=20, GAM_ISOLATION=21,
    ZMASSCUT=22, LLGMASSCUT=23, LLMASSCUT=24, DILEP_PT_FRAC=25, GAM_PT_FRAC=26, PASSALL=27
  };

  // names of all cuts (do not includ "pass all")
  const std::vector<TString> s_cutDescs =
    {"Lepton Dalitz truth","No duplicates","GRL","Pass trigger","Detector DQ","Has PV",
     "2 same-flavor leptons","1 loose photon","e-#gamma ambiguity / HV",
     "Z-boson assignment","2 same-flavor leptons (post-OR)","Bad muon","1 loose photon (post-OR)",
     "Trigger match","lepton ID","lepton impact parameter","lepton isolation",
     "photon tight ID","photon isolation",
     "#it{m}_{ll} < 45 GeV",
     "#it{m}_{ll#gamma} #in [105,160] GeV",
     "#it{m}_{ll} J/#Psi/#Upsilon window",
     "p^{ll}_{T}/#it{m}_{ll#gamma} > 0.3",
     "p^{#gamma}_{T}/#it{m}_{ll#gamma} > 0.3" };

  /// value of cut that fail selection: PASSALL if all cuts passed
  CutEnum m_cutFlow;

  // Lepton ID and isolation scale factors, used to apply the right weights at the right time in the
  // cutflow.
  float m_lepIDWeight;
  float m_lepIsoWeight;

  // names of the output containers
  TString m_photonContainerName, m_jetContainerName, m_elecContainerName, m_muonContainerName, m_trackContainerName;
  TString m_photonTruthContainerName, m_jetTruthContainerName, m_elecTruthContainerName, m_muonTruthContainerName;
  TString m_evtInfoName, m_truthEvtsName;

  HG::Iso::IsolationType m_eleMergedIsoWP;
  HG::Iso::IsolationType m_eleResolvedIsoWP;
  HG::Iso::IsolationType m_muonIsoWP;

  TString m_eleIDPreselection;

  // what skimming to apply
  int m_skimCut;

  // whether it's a Dalitz event
  bool m_isNonHyyStarHiggs;

  // store the crossSectionBRfilterEff value once per file
  float m_crossSectionBRfilterEff;

  // whether to apply systematics, save the differential variables and the truth
  bool m_applySystematics, m_saveObjects, m_saveTruthObjects, m_saveTruthVars;
  bool m_skipElectronObjects, m_skipMuonObjects;

  // Containers
  xAOD::PhotonContainer m_allPhotons; //!
  xAOD::PhotonContainer m_preSelPhotons; //!
  xAOD::PhotonContainer m_selPhotons; //!

  xAOD::JetContainer m_allJets; //!
  xAOD::JetContainer m_selJets; //!
  xAOD::JetContainer m_jvtJets; //!

  xAOD::ElectronContainer m_allElectrons; //!
  xAOD::ElectronContainer m_selElectrons; //!
  // xAOD::ElectronContainer m_preSelElectrons; //!

  xAOD::TrackParticleContainer m_allTracks; //!
  xAOD::TrackParticleContainer m_preSelTracks; //!
  xAOD::TrackParticleContainer m_selTracks; //!

  xAOD::MuonContainer m_allMuons; //!
  xAOD::MuonContainer m_selMuons; //!
  xAOD::MuonContainer m_preSelMuons; //!

  std::map<HG::Iso::IsolationType, CP::IsolationCloseByCorrectionTool*> m_isoCloseByTools_Ele; //!
  std::map<HG::Iso::IsolationType, CP::IsolationCloseByCorrectionTool*> m_isoCloseByTools_Muon; //!

  std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_eleIsoAccCorr; //!
  std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_muIsoAccCorr; //!

private:
  /// helper methods to create, fill and print the cut flow

  TH1F* getCutFlowHisto(bool weighted,
                        bool onlyDalitz,
                        HG::ChannelEnum truthChannel=HG::CHANNELUNKNOWN) {

    // 2147483647
    // 1000000000 onlyDalitz
    //  100000000 weighted
    //   XX000000 channel (up to 99 channels apparently)
    //     345678 sample ID

    // We are just looking to get a unique integer here, to keep track of individual histos.
    int ID = getSampleID();

    int map_i = ID;
    if (onlyDalitz)   map_i += (int)1e9;
    if (weighted)     map_i += (int)1e8;
    if (truthChannel) map_i += (int)(1e6 * (int)truthChannel);

    if (TH1F *h = m_cFlowHistos[map_i]) return h;

    TString suffix = "";
    if (onlyDalitz) suffix += "_onlyDalitz";
    if (truthChannel != HG::CHANNELUNKNOWN) suffix += TString("_") + GetChannelName(truthChannel);
    if (weighted) suffix += "_weighted";

    m_cFlowHistos[map_i] = makeCutFlowHisto(ID, s_cutDescs, suffix);
    return m_cFlowHistos[map_i];
  }

  void printCutFlowHistos(void);

  /// \brief fill the cut flow histograms
  void fillCutFlow(CutEnum cut, double w);

  // apply cut flow
  // returns enum corresponding to failed cut - or PASSALL if all cuts passed
  CutEnum cutflow();
  EL::StatusCode doReco(bool isSys = false);
  EL::StatusCode doTruth();

  // Need to be able to set the truth channel earlier, for cutflow purposes.
  void SetTruthHiggsInformation(void);

  void decorateCorrectedIsoCut(xAOD::ElectronContainer & electrons, xAOD::MuonContainer & muons);
  void AddElectronDecorations(xAOD::ElectronContainer& electrons);
  void AddMuonDecorations(xAOD::MuonContainer& muons);

  HG::ChannelEnum FindZboson_ElectronChannelAware(xAOD::TrackParticleContainer* inTracks,
                                                  xAOD::TrackParticle*& sel_trk1,
                                                  xAOD::TrackParticle*& sel_trk2,
                                                  double& return_mll,
                                                  const HG::TrackElectronMap& trkEleMap,
                                                  xAOD::ElectronContainer* inEleCont,
                                                  xAOD::ElectronContainer* outEleCont);
  //configurable overlap removal cone sizes
  float m_OR_e_DR_y;
  float m_OR_jet_DR_y;
  float m_OR_jet_DR_e;
  float m_OR_e_DR_jet;
  float m_OR_mu_DR_y;
  float m_OR_mu_DR_jet;

private:
#ifndef __CINT__
  HG::TrackHandler *m_trackHandler; //!
  HG::MergedElectronID * m_mergedElectronID; //!
  HG::MergedElectronID_v2 * m_mergedElectronID_v2; //!
#endif // __CINT__

protected:
  inline virtual HG::TrackHandler *trackHandler() { return m_trackHandler; }

public:

  // this is a standard constructor
  HiggsGamGamStarCutflowAndMxAOD() { }
  HiggsGamGamStarCutflowAndMxAOD(const char *name);
  virtual ~HiggsGamGamStarCutflowAndMxAOD();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode initialize();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode finalize();
  virtual EL::StatusCode fileExecute();
  virtual EL::StatusCode changeInput(bool firstFile);

  // Variables to write when using FULL_v1 photon systematics
  // (separated from the other systs due to file size)
  void writePhotonAllSys(bool isSys);

  // Variables to write to nominal file and (non-photon) systematics file
  void writeNominalAndSystematic();

  // Variables to write only to the nominal file
  void writeNominalOnly();

  // Detailed variables to write if "SaveDetailedVariables" is set (only called for nominal case)
  void writeDetailed();

  // Functions for writing variables (E.g. "writeBlahVars" is called by "writeBlah")
  // Use these for variables for which a truth and a reco variable exists.
  void writePhotonAllSysVars(bool truth = false);
  void writeNominalAndSystematicVars(bool truth = false);
  void writeNominalOnlyVars(bool truth = false);
  void writeDetailedVars(bool truth = false);
  void writeTruthOnlyVars();

  // this is needed to distribute the algorithm to the workers
  ClassDef(HiggsGamGamStarCutflowAndMxAOD, 1);
};

#endif // HGamGamStar_HiggsGamGamStarCutflowAndMxAOD_H
