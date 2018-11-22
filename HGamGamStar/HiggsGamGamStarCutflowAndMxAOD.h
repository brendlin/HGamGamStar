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
    TRIG_MATCH=16, LEP_MEDID=17, LEP_IP=18, LEP_ISO=19, GAM_TIGHTID=20, GAM_ISOLATION=21, ZMASSCUT=22, LLGMASSCUT=23, PASSALL=24
  };

  // names of all cuts (do not includ "pass all")
  const std::vector<TString> s_cutDescs =
    {"Lepton Dalitz truth","No duplicates","GRL","Pass trigger","Detector DQ","Has PV",
     "2 same-flavor leptons","1 loose photon","e-#gamma ambiguity","Z-boson assignment",
     "2 same-flavor leptons (post-OR)","Bad muon","1 loose photon (post-OR)",
     "Trigger match","tight ID","isolation","#it{m}_{ll} < 45 GeV",
     "#it{m}_{ll#gamma} #in [105,160] GeV"};

  /// value of cut that fail selection: PASSALL if all cuts passed
  CutEnum m_cutFlow;

  // names of the output containers
  TString m_photonContainerName, m_jetContainerName, m_elecContainerName, m_muonContainerName, m_trackContainerName;
  TString m_photonTruthContainerName, m_jetTruthContainerName, m_elecTruthContainerName, m_muonTruthContainerName;
  TString m_evtInfoName, m_truthEvtsName;
  
  TString m_electronIsoWP;
  TString m_muonIsoWP;

  TString m_eleIDPreselection;

  // what skimming to apply
  int m_skimCut;

  // whether it's a Dalitz event
  bool m_isNonHyyStarHiggs;

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

  xAOD::TrackParticleContainer m_allTracks; //!
  xAOD::TrackParticleContainer m_preSelTracks; //!
  xAOD::TrackParticleContainer m_selTracks; //!

  xAOD::MuonContainer m_allMuons; //!
  xAOD::MuonContainer m_selMuons; //!
  xAOD::MuonContainer m_preSelMuons; //!
  
  CP::IsolationCloseByCorrectionTool* m_isoCloseByTool_Electron; //!
  CP::IsolationSelectionTool* m_isoSelTool_Electron; //!
  
  CP::IsolationCloseByCorrectionTool* m_isoCloseByTool_Muon; //!
  CP::IsolationSelectionTool* m_isoSelTool_Muon; //!

private:
  /// helper methods to create, fill and print the cut flow

  TH1F* getCutFlowHisto(bool onlyDalitz=false) {
    int ID = getSampleID()*(onlyDalitz?-1:1);
    if (TH1F *h = m_cFlowHistos[ID]) return h;
    m_cFlowHistos[ID] = makeCutFlowHisto(ID, s_cutDescs, onlyDalitz?"_onlyDalitz":"");
    return m_cFlowHistos[ID];
  }

  TH1F* getCutFlowWeightedHisto(bool onlyDalitz=false) {
    int ID = getSampleID()*(onlyDalitz?-1:1);
    if (TH1F *h = m_cFlowHistosWeighted[ID]) return h;
    m_cFlowHistosWeighted[ID] = makeCutFlowHisto(ID, s_cutDescs,
                                                 onlyDalitz?"_onlyDalitz_weighted":"_weighted");
    return m_cFlowHistosWeighted[ID];
  }

  /// \brief fill the cut flow histograms
  void fillCutFlow(CutEnum cut, double w);

  // apply cut flow
  // returns enum corresponding to failed cut - or PASSALL if all cuts passed
  CutEnum cutflow();
  EL::StatusCode doReco(bool isSys = false);
  EL::StatusCode doTruth();
  
  void decorateCorrectedIsoCut(xAOD::ElectronContainer & electrons, xAOD::MuonContainer & muons);

  HG::ChannelEnum truthClass();
  HG::ChannelEnum ClassifyElectronChannelsByBestMatch(const xAOD::TrackParticle* trk0,
                                                      const xAOD::TrackParticle* trk1,
                                                      const HG::TrackElectronMap& trkEleMap,
                                                      xAOD::ElectronContainer* inEleCont=nullptr,
                                                      xAOD::ElectronContainer* outEleCont=nullptr);

  HG::ChannelEnum ClassifyElectronsOld(xAOD::TrackParticle* trk0,
                                       xAOD::TrackParticle* trk1,
                                       const HG::TrackElectronMap& trkEleMap,
                                       xAOD::ElectronContainer* inEleCont=nullptr,
                                       xAOD::ElectronContainer* outEleCont=nullptr);

private:
#ifndef __CINT__
  HG::TrackHandler *m_trackHandler; //!
  HG::MergedElectronID * m_mergedElectronID; //!
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
  void writePhotonAllSysVars(bool truth = false);
  void writeNominalAndSystematicVars(bool truth = false);
  void writeNominalOnlyVars(bool truth = false);
  void writeDetailedVars(bool truth = false);
  void writeTruthOnlyVars();

  // this is needed to distribute the algorithm to the workers
  ClassDef(HiggsGamGamStarCutflowAndMxAOD, 1);
};

#endif // HGamGamStar_HiggsGamGamStarCutflowAndMxAOD_H
