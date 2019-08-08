#ifndef HGamGamStar_RadiativeZCutflowAndMxAOD_H
#define HGamGamStar_RadiativeZCutflowAndMxAOD_H

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamGamStar/MxAODTool.h"
#include "AsgTools/ToolHandle.h"
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/TrackHandler.h"
#include "HGamGamStar/MergedElectronID.h"
#include "HGamGamStar/MergedElectronID_v2.h"

class RadiativeZCutflowAndMxAOD : public MxAODTool
{

private:

  // Cut-flow - need to keep the same order!
  enum CutEnum {
    // Covered in MxAODTool: NxAOD=0, NDxAOD=1, ALLEVTS=2
    DUPLICATE=3, GRL=4, TRIGGER=5, DQ=6, VERTEX=7,
    TWO_SF_LEPTONS=8,ONE_RECO_GAM=9, AMBIGUITY=10, HVCUT=11,
    ZBOSON_ASSIGNMENT=12,TWO_SF_LEPTONS_POSTOR=13,BAD_MUON=14,ONE_PHOTON_POSTOR=15,
    TRIG_MATCH=16,
    GAM_TIGHTID=17, GAM_ISOLATION=18,
    ZMASSCUT=19, LLGMASSCUT=20, PASSALL=21
  };

  // names of all cuts (do not includ "pass all")
  const std::vector<TString> s_cutDescs =
    {"No duplicates","GRL","Pass trigger","Detector DQ","Has PV",
     "2 same-flavor leptons","1 reco photon","e-#gamma ambiguity","HV cut",
     "Z-boson assignment","2 same-flavor leptons (post-OR)","Bad muon","1 reco photon (post-OR)",
     "Trigger match",
     "photon tight ID","photon isolation",
     "#it{m}_{ll} > 10 GeV",
     "#it{m}_{ll#gamma} #in [80,100] GeV",};

  /// value of cut that fail selection: PASSALL if all cuts passed
  CutEnum m_cutFlow;

  // Lepton ID and isolation scale factors, used to apply the right weights at the right time in the
  // cutflow.
  float m_lepWeight;
  float m_phIDWeight;
  float m_phIsoWeight;

  // names of the output containers
  TString m_photonContainerName, m_jetContainerName, m_elecContainerName, m_muonContainerName, m_trackContainerName;
  TString m_photonTruthContainerName, m_jetTruthContainerName, m_elecTruthContainerName, m_muonTruthContainerName;
  TString m_evtInfoName, m_truthEvtsName;

  TString m_eleIDPreselection;

  // what skimming to apply
  int m_skimCut;

  // store the crossSectionBRfilterEff value once per file
  float m_crossSectionBRfilterEff;

  // whether to apply systematics, save the differential variables and the truth
  bool m_applySystematics, m_saveObjects, m_saveTruthObjects, m_saveTruthVars;

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

  xAOD::MuonContainer m_allMuons; //!
  xAOD::MuonContainer m_selMuons; //!
  xAOD::MuonContainer m_preSelMuons; //!

private:
  /// helper methods to create, fill and print the cut flow

  TH1F* getCutFlowHisto(bool weighted,
                        HG::ChannelEnum truthChannel=HG::CHANNELUNKNOWN) {

    // 2147483647
    //  100000000 weighted
    //   XX000000 channel (up to 99 channels apparently)
    //     345678 sample ID

    // We are just looking to get a unique integer here, to keep track of individual histos.
    int ID = getSampleID();

    int map_i = ID;
    if (weighted)     map_i += (int)1e8;
    if (truthChannel) map_i += (int)(1e6 * (int)truthChannel);

    if (TH1F *h = m_cFlowHistos[map_i]) return h;

    TString suffix = "";
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

  // Space for truth cuts, eventually
  void SetTruthZBosonInformation(void);

  void AddElectronDecorations(xAOD::ElectronContainer& electrons);

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
  RadiativeZCutflowAndMxAOD() { }
  RadiativeZCutflowAndMxAOD(const char *name);
  virtual ~RadiativeZCutflowAndMxAOD();

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
  ClassDef(RadiativeZCutflowAndMxAOD, 1);
};

#endif // HGamGamStar_RadiativeZCutflowAndMxAOD_H
