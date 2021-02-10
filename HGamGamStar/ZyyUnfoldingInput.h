#ifndef HGamGamStar_ZyyUnfoldingInput_H
#define HGamGamStar_ZyyUnfoldingInput_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class ZyyUnfoldingInput : public HgammaAnalysis {
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  std::map<int, std::map<std::string, int> > m_variations;

  std::map<int, int> m_nEvents;
  TH1F *createAndRegisterTH1F(TString name, int Nbins, double min, double max, TString title = "");
  TH1F *createAndRegisterTH1F(TString name, const std::vector<double>& xbins, TString title = "");
  TH2F *createAndRegisterTH2F(TString name, const std::vector<double>& xbins, const std::vector<double>& ybins, TString title = "");

public:
  // this is a standard constructor
  ZyyUnfoldingInput() { }
  ZyyUnfoldingInput(const char *name);
  virtual ~ZyyUnfoldingInput();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode doSelAndSaveHists(TString sysname);
  virtual EL::StatusCode CreateHists(TString sysname);

  bool isFiducial=false;
  double truew;
  double vals[7];
  double truevals[7];
  std::vector<double> m_pdfw, m_scalew;

  // this is needed to distribute the algorithm to the workers
  ClassDef(ZyyUnfoldingInput, 1);
};

#endif // HGamGamStar_ZyyUnfoldingInput_H
