#ifndef HGamGamStar_HggStarCutflowAndMxAOD_H
#define HGamGamStar_HggStarCutflowAndMxAOD_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "FsrUtils/FsrPhotonTool.h"
#include "AsgTools/ToolHandle.h"
#include "ZMassConstraint/IConstraintFit.h"

class HggStarCutflowAndMxAOD : public HgammaAnalysis
{

private:

  // Cut-flow - need to keep the same order!
  enum CutEnum {
    NxAOD=0, NDxAOD=1, ALLEVTS=2, DUPLICATE=3, GRL=4, TRIGGER=5, DQ=6, VERTEX=7,
    TWO_SF_LEPTONS=8,ONE_LOOSE_GAM=9, AMBIGUITY=10,
    TWO_SF_LEPTONS_POSTOR=11,ONE_PHOTON_POSTOR=12,
    TRIG_MATCH=13, GAM_TIGHTID=14, GAM_ISOLATION=15, RELPTCUTS=16, MASSCUT=17, PASSALL=18 };
  // names of all cuts (do not includ "pass all")
  const std::vector<TString> s_cutDescs =
    {"#it{N}_{xAOD}","#it{N}_{DxAOD}","All events","No duplicates","GRL","Pass trigger","Detector DQ","Has PV",
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

  // whether it's a Dalitz event
  bool m_isDalitz;

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

  // cut-flow histograms and book keeper yields
  std::map<int,TH1F*> m_cFlowHistos, m_cFlowHistosWeighted;
  int m_N_xAOD, m_N_DxAOD;
  double m_sumw_xAOD, m_sumw2_xAOD, m_sumw_DxAOD, m_sumw2_DxAOD;
  bool m_newFile;

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

  /// \brief create a new cut-flow histogram
  TH1F *makeCutFlowHisto(int id, TString suffix="");

  void addBookKeeping(TH1F *cutflow, double N_xAOD, double N_DxAOD,
		      double sumw2_xAOD=-1, double sumw2_DxAOD=-1) {
    int bin_xAOD=cutflow->FindBin(NxAOD), bin_DxAOD=cutflow->FindBin(NDxAOD);
    cutflow->AddBinContent(bin_xAOD,N_xAOD);
    cutflow->AddBinContent(bin_DxAOD,N_DxAOD);
    if (sumw2_xAOD>0.0) {
      cutflow->SetBinError(bin_xAOD, sqrt(pow(cutflow->GetBinError(bin_xAOD),2)  + sumw2_xAOD ));
      cutflow->SetBinError(bin_DxAOD,sqrt(pow(cutflow->GetBinError(bin_DxAOD),2) + sumw2_DxAOD));
    }
  }

  /// \brief get the "sample ID", meaning run number for data and MC channel number for MC
  inline int getSampleID() { return HG::isMC() ? eventInfo()->mcChannelNumber() : eventInfo()->runNumber(); }

  TH1F* getCutFlowHisto(bool onlyDalitz=false) {
    int ID = getSampleID()*(onlyDalitz?-1:1);
    if (TH1F *h = m_cFlowHistos[ID]) return h;
    m_cFlowHistos[ID] = makeCutFlowHisto(ID,onlyDalitz?"":"_onlyDalitz");
    return m_cFlowHistos[ID];
  }

  TH1F* getCutFlowWeightedHisto(bool onlyDalitz=false) {
    int ID = getSampleID()*(onlyDalitz?-1:1);
    if (TH1F *h = m_cFlowHistosWeighted[ID]) return h;
    m_cFlowHistosWeighted[ID] = makeCutFlowHisto(ID,onlyDalitz?"_weighted":"_onlyDalitz_weighted");
    return m_cFlowHistosWeighted[ID];
  }

  /// \brief fill the cut flow histograms
  void fillCutFlow(CutEnum cut, double w);

  /// \brief print a given cut flow histogram
  void printCutFlowHisto(TH1F *h, int Ndecimals=0);

  /// \brief print all cut flow histograms stored in the maps
  void printCutFlowHistos();

  // apply cut flow
  // returns enum corresponding to failed cut - or PASSALL if all cuts passed
  CutEnum cutflow();
  EL::StatusCode doReco(bool isSys = false);
  EL::StatusCode doTruth();

  /// Declares list of output variables to be written.
  /// configKey defines the key that specifies the list of variable names
  /// the name in the output file will be of the form outName+"Aux."+VARNAME
  void declareOutputVariables(TString outName, TString configKey, StrV extra = {}, StrV ignore = {});

  double diphotonMassResolution(xAOD::PhotonContainer &gams) {
    static SG::AuxElement::Accessor<float> dE("relEreso");
    if (gams.size()<2) return -99.0;
    // m_yy = sqrt(E1*E2*A), where the angular term A=2*(1-cos(alpha)) = 2*(cosh(Deta)-cos(Dphi))
    // error propagation gives
    //   Dm_yy / m_yy = 1/2*DE1/E1 (+) 1/2*DE2/E2 (+) dA/A
    // ignoring angular resolution, and with dE = DE/E =>
    return sqrt( pow(0.5*dE(*gams[0]),2) + pow(0.5*dE(*gams[1]),2) );
  }

public:

  // this is a standard constructor
  HggStarCutflowAndMxAOD() { }
  HggStarCutflowAndMxAOD(const char *name);
  virtual ~HggStarCutflowAndMxAOD();

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
  ClassDef(HggStarCutflowAndMxAOD, 1);
};

#endif // HGamGamStar_HggStarCutflowAndMxAOD_H
