#ifndef HGamGamStar_MxAODTool_H
#define HGamGamStar_MxAODTool_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class MxAODTool : public HgammaAnalysis
{

public:
  // this is a standard constructor
  MxAODTool() { }
  MxAODTool(const char *name);
  virtual ~MxAODTool();

 protected:

  enum CutEnum { NxAOD=0, NDxAOD=1, ALLEVTS=2 };

  const std::vector<TString> s_cutDescs_base = {"#it{N}_{xAOD}","#it{N}_{DxAOD}","All events"};

  bool m_saveDetailed;
  bool m_newFileMetaData;

  /// Declares list of output variables to be written.
  /// configKey defines the key that specifies the list of variable names
  /// the name in the output file will be of the form outName+"Aux."+VARNAME
  virtual void declareOutputVariables(TString outName, TString configKey, StrV extra = {}, StrV ignore = {});

  // cut-flow histograms and book keeper yields
  std::map<int,TH1F*> m_cFlowHistos, m_cFlowHistosWeighted;

  // Add details to bookkeeping
  void addBookKeeping(TH1F *cutflow, double N_xAOD, double N_DxAOD,
                      double sumw2_xAOD=-1, double sumw2_DxAOD=-1);

  EL::StatusCode fillCutFlowWithBookkeeperInfo();

  /// \brief create a new cut-flow histogram
  TH1F *makeCutFlowHisto(int id, const std::vector<TString> otherCuts, TString suffix="");

  void printCutFlowHistos();
  void printCutFlowHisto(TH1F *h, int Ndecimals);

  /// \brief get the "sample ID", meaning run number for data and MC channel number for MC
  inline int getSampleID()
  {
    return HG::isMC() ? eventInfo()->mcChannelNumber() : eventInfo()->runNumber();
  }

  // this is needed to distribute the algorithm to the workers
  ClassDef(MxAODTool, 1);
};

#endif // HGamGamStar_MxAODTool_H
