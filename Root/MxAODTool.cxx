#include "HGamGamStar/MxAODTool.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

// this is needed to distribute the algorithm to the workers
ClassImp(MxAODTool)



MxAODTool::MxAODTool(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



MxAODTool::~MxAODTool()
{
  // Here you delete any memory you allocated during your analysis.
}



void MxAODTool::declareOutputVariables(TString outName, TString configKey, StrV extra, StrV ignore) {
  if (config()->isDefined(configKey)) {
    TString vars = config()->getStr(configKey).Data();

    if (m_saveDetailed) {
      TString detailKey = configKey.ReplaceAll("Variables", "DetailedVariables");
      if (config()->isDefined(detailKey)) {
        TString detailed = config()->getStr(detailKey);
        vars += "." + detailed;
      }
    }

    for (TString val: extra)
      vars += val;

    for (TString val: ignore)
      vars = vars.ReplaceAll(val, "");

    event()->setAuxItemList((outName+"Aux.").Data(), vars.Data());
  }
  else HG::fatal("Cannot find "+configKey);
}

