#ifndef HGamGamStar_MergedElectronMxAOD_H
#define HGamGamStar_MergedElectronMxAOD_H

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamGamStar/MxAODTool.h"
/* #include "FsrUtils/FsrPhotonTool.h" */
#include "AsgTools/ToolHandle.h"
#include <AsgTools/AnaToolHandle.h>


/* #include "ZMassConstraint/IConstraintFit.h" */
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/TrackHandler.h"
#include "HGamGamStar/MergedElectronID.h"

#include "IsolationSelection/IsolationCloseByCorrectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"

namespace InDet {
    class IInDetTrackSelectionTool;
}

namespace CP {
    class ITrackVertexAssociationTool;
}

class MergedElectronMxAOD : public MxAODTool
{

private:

  // Cut-flow - need to keep the same order!
  enum CutEnum {
    // Covered in MxAODTool: NxAOD=0, NDxAOD=1, ALLEVTS=2
    DUPLICATE=3, GRL=4,  DQ=5, VERTEX=6, ELECTRON=7, PASSALL=8
  };

  // names of all cuts (do not includ "pass all")
  const std::vector<TString> s_cutDescs =  {"No duplicates","GRL","Detector DQ",
                                            "Has PV","Electron","Pass"};

  /// value of cut that fail selection: PASSALL if all cuts passed
  CutEnum m_cutFlow;

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

  // Containers

  xAOD::ElectronContainer m_selElectrons; //!


  std::map<HG::Iso::IsolationType, CP::IsolationCloseByCorrectionTool*> m_isoCloseByTools_Ele; //!
  std::map<HG::Iso::IsolationType, CP::IsolationCloseByCorrectionTool*> m_isoCloseByTools_Muon; //!

  std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_eleIsoAccCorr; //!
  std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_muIsoAccCorr; //!

private:
  /// helper methods to create, fill and print the cut flow

  TH1F* getCutFlowHisto(bool onlyDalitz=false) {
    int ID = getSampleID()*(onlyDalitz?-1:1);
    if (TH1F *h = m_cFlowHistos[ID]) return h;
    m_cFlowHistos[ID] = makeCutFlowHisto(ID, s_cutDescs, "");
    return m_cFlowHistos[ID];
  }

  TH1F* getCutFlowWeightedHisto(bool onlyDalitz=false) {
    int ID = getSampleID()*(onlyDalitz?-1:1);
    if (TH1F *h = m_cFlowHistosWeighted[ID]) return h;
    m_cFlowHistosWeighted[ID] = makeCutFlowHisto(ID, s_cutDescs,"_weighted");
    return m_cFlowHistosWeighted[ID];
  }

  /// \brief fill the cut flow histograms
  void fillCutFlow(CutEnum cut, double w);

  // apply cut flow
  // returns enum corresponding to failed cut - or PASSALL if all cuts passed
  CutEnum cutflow();
  EL::StatusCode doReco(bool isSys = false);
  EL::StatusCode doTruth();

  void decorateCorrectedIsoCut(xAOD::ElectronContainer & electrons);
  void AddElectronDecorations(xAOD::ElectronContainer& electrons);


private:
#ifndef __CINT__
  HG::TrackHandler *m_trackHandler; //!
  HG::MergedElectronID * m_mergedElectronID; //!

  asg::AnaToolHandle<InDet::IInDetTrackSelectionTool> m_trkselTool; //!
  asg::AnaToolHandle<CP::ITrackVertexAssociationTool> m_ttvaTool; //!
#endif // __CINT__

protected:
  inline virtual HG::TrackHandler *trackHandler() { return m_trackHandler; }

public:

  // this is a standard constructor
  MergedElectronMxAOD() { }
  MergedElectronMxAOD(const char *name);
  virtual ~MergedElectronMxAOD();

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
  void writeNominalAndSystematicVars(bool truth = false);
  void writeNominalOnlyVars(bool truth = false);
  void writeDetailedVars(bool truth = false);
  void writeTruthOnlyVars();

  // this is needed to distribute the algorithm to the workers
  ClassDef(MergedElectronMxAOD, 1);
};

#endif // HGamGamStar_MergedElectronMxAOD_H
