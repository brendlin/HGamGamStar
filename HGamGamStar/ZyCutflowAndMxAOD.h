#ifndef HGamGamStar_ZyCutflowAndMxAOD_H
#define HGamGamStar_ZyCutflowAndMxAOD_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamGamStar/MxAODTool.h"
#include "FsrUtils/FsrPhotonTool.h"
#include "AsgTools/ToolHandle.h"
#include "ZMassConstraint/IConstraintFit.h"

class ZyCutflowAndMxAOD : public MxAODTool
{

private:

  // Cut-flow - need to keep the same order!
  enum CutEnum {
    // Covered in MxAODTool: NxAOD=0, NDxAOD=1, ALLEVTS=2
    DUPLICATE=3, GRL=4, TRIGGER=5, DQ=6, VERTEX=7,
    TWO_SF_LEPTONS=8,ONE_LOOSE_GAM=9, AMBIGUITY=10,
    TWO_SF_LEPTONS_POSTOR=11,ONE_PHOTON_POSTOR=12,
    TRIG_MATCH=13, GAM_TIGHTID=14, GAM_ISOLATION=15, RELPTCUTS=16, MASSCUT=17, PASSALL=18 };

  // names of all cuts (do not includ "pass all")
  const std::vector<TString> s_cutDescs =
    {"No duplicates","GRL","Pass trigger","Detector DQ","Has PV",
     "2 same-flavor leptons","1 loose photon","e-#gamma ambiguity",
     "2 same-flavor leptons (post-OR)","1 loose photon (post-OR)",
     "Trigger match","tight ID","isolation","rel. #it{p}_{T} cuts",
     "#it{m}_{#gamma#gamma} #in [105,160] GeV"};

  /// value of cut that fail selection: PASSALL if all cuts passed
  CutEnum m_cutFlow;

  // names of the output containers
  TString m_photonContainerName, m_jetContainerName, m_elecContainerName, m_muonContainerName;
  TString m_photonTruthContainerName, m_jetTruthContainerName, m_elecTruthContainerName, m_muonTruthContainerName;
  TString m_evtInfoName, m_truthEvtsName;

  // what skimming to apply
  int m_skimCut;

  // whether to apply systematics, save the differential variables and the truth
  bool m_applySystematics, m_saveObjects, m_saveDetailed, m_saveTruthObjects, m_saveTruthVars;
  bool m_allowMoreThanTwoPhotons;

  // whether to save fake photon combinations
  bool m_enableFakePhotons;
  //If we have two good fakes then we need to pass the slimming
  bool m_goodFakeComb;

  //Whether we are running yybb-tool in detailed mode
  bool m_detailedHHyybb;

  // Temporary flag for photon all sys
  bool m_photonAllSys;

  // Tools
  ToolHandle<ZMassConstraint::IConstraintFit> m_massConstraint;

  // Containers
  xAOD::PhotonContainer m_allPhotons; //!
  xAOD::PhotonContainer m_preSelPhotons; //!
  xAOD::PhotonContainer m_selPhotons; //!

  xAOD::JetContainer m_allJets; //!
  xAOD::JetContainer m_selJets; //!
  xAOD::JetContainer m_jvtJets; //!

  xAOD::ElectronContainer m_allElectrons; //!
  xAOD::ElectronContainer m_selElectrons; //!
  xAOD::ElectronContainer m_preSelElectrons; //!

  xAOD::MuonContainer m_allMuons; //!
  xAOD::MuonContainer m_selMuons; //!
  xAOD::MuonContainer m_preSelMuons; //!

  xAOD::MissingETContainer m_allMET; //!
  xAOD::MissingETContainer m_selMET; //!

private:
  /// helper methods to create, fill and print the cut flow

  TH1F* getCutFlowHisto() {
    int ID = getSampleID();
    if (TH1F *h = m_cFlowHistos[ID]) return h;
    m_cFlowHistos[ID] = makeCutFlowHisto(ID, s_cutDescs);
    return m_cFlowHistos[ID];
  }

  TH1F* getCutFlowWeightedHisto() {
    int ID = getSampleID();
    if (TH1F *h = m_cFlowHistosWeighted[ID]) return h;
    m_cFlowHistosWeighted[ID] = makeCutFlowHisto(ID, s_cutDescs, "_weighted");
    return m_cFlowHistosWeighted[ID];
  }

  /// \brief fill the cut flow histograms
  void fillCutFlow(CutEnum cut, double w);

  // apply cut flow
  // returns enum corresponding to failed cut - or PASSALL if all cuts passed
  CutEnum cutflow();
  EL::StatusCode doReco(bool isSys = false);
  EL::StatusCode doTruth();

public:

  // this is a standard constructor
  ZyCutflowAndMxAOD() { }
  ZyCutflowAndMxAOD(const char *name);
  virtual ~ZyCutflowAndMxAOD();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode finalize();
  virtual EL::StatusCode fileExecute();

  // Functions for saving information
  void writePhotonAllSys(bool isSys);
  void writeNominalAndSystematic(bool isSys);
  void writeNominalOnly();
  void writeDetailed();

  // Functions for writting variables
  void writePhotonAllSysVars(bool truth = false);
  void writeNominalAndSystematicVars(bool truth = false);
  void writeNominalOnlyVars(bool truth = false);
  void writeDetailedVars(bool truth = false);

  // this is needed to distribute the algorithm to the workers
  ClassDef(ZyCutflowAndMxAOD, 1);
};

#endif // HGamGamStar_ZyCutflowAndMxAOD_H
