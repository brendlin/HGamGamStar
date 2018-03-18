#ifndef HGamGamStar_MxAODTool_H
#define HGamGamStar_MxAODTool_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class MxAODTool : public HgammaAnalysis
{

public:

private:

public:
  // this is a standard constructor
  MxAODTool() { }
  MxAODTool(const char *name);
  virtual ~MxAODTool();

 protected:

  bool m_saveDetailed;

  /// Declares list of output variables to be written.
  /// configKey defines the key that specifies the list of variable names
  /// the name in the output file will be of the form outName+"Aux."+VARNAME
  virtual void declareOutputVariables(TString outName, TString configKey, StrV extra = {}, StrV ignore = {});

  // this is needed to distribute the algorithm to the workers
  ClassDef(MxAODTool, 1);
};

#endif // HGamGamStar_MxAODTool_H
