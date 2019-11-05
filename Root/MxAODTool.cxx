#include "HGamGamStar/MxAODTool.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

// this is needed to distribute the algorithm to the workers
ClassImp(MxAODTool)



MxAODTool::MxAODTool(const char *name)
: HgammaAnalysis(name),
  m_newFileMetaData(false)
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



EL::StatusCode MxAODTool::fillCutFlowWithBookkeeperInfo()
{
  // Call this from execute(), with m_newFileMetaData = true;

  // cutflow bookkeeper yields
  int N_xAOD=0, N_DxAOD=0;
  double sumw_xAOD=0.0, sumw2_xAOD=0.0, sumw_DxAOD=0.0, sumw2_DxAOD=0.0;

  asg::AsgMetadataTool amdt("MetaDataTool");
  const xAOD::FileMetaData *fmd = nullptr;

  if (amdt.inputMetaStore()->retrieve(fmd, "FileMetaData").isFailure())
  { HG::fatal("Cannot find FileMetaData in the input file."); }

  // Get the dataType (e.g. "StreamDAOD_HIGG1D1")
  std::string dataType;
  fmd->value(xAOD::FileMetaData::dataType, dataType);

  // Get the DAOD name (e.g. "HIGG1D1")
  TString str_daod = dataType;
  str_daod.ReplaceAll("StreamDAOD_", "");

  if ( !HG::isDAOD() )
  {
    return EL::StatusCode::FAILURE;
  }

  // If we get here, this is a DxAOD

  // 1. Check if there if the incomplete book-keeper object has entreies (that would be bad!)
  const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
  EL_CHECK( "fileExecute()", wk()->xaodEvent()->retrieveMetaInput(incompleteCBC,"IncompleteCutBookkeepers") );
  if ( incompleteCBC->size() != 0 ) {
    // fatal or error? let's start with a hard fatal!
    // HG::fatal("Issue with DxAOD book keeper. It's incomplete. File corrupted?");
    Warning("fileExecute()",
            "Issue with DxAOD book keeper. It's incomplete. File corrupted? %s %s",
            "If this is data, this is known to happen (but not understood).",
            "If this is MC, this is expected to happen, and can probably be ignored.");
  }

  // 2. now get the actual bookkeeper
  const xAOD::CutBookkeeperContainer* completeCBC = nullptr;
  EL_CHECK( "fileExecute()", wk()->xaodEvent()->retrieveMetaInput(completeCBC,"CutBookkeepers") );

  int maxAcycle = -1, maxDcycle = -1;
  for ( auto cbk : *completeCBC ) {
    if (TString(cbk->name()).Contains(str_daod))
      Info("fileExecute()","  Book keeper name=%s, inputStream=%s, cycle=%d, nAcceptedEvents=%d", cbk->name().c_str(), cbk->inputStream().c_str(), cbk->cycle(), (int)cbk->nAcceptedEvents());

    if ( cbk->name().empty() ) continue;

    // Use the DxAOD numbers from the largest cycle
    // str_daod is filled above. E.g. Use "HIGG1D1" or "HIGG1D2" or ...
    if (TString(cbk->name()).Contains(str_daod) && cbk->inputStream() == "StreamAOD" && cbk->cycle() > maxDcycle) {
      maxDcycle     = cbk->cycle();
      N_DxAOD     = cbk->nAcceptedEvents();
      sumw_DxAOD  = cbk->sumOfEventWeights();
      sumw2_DxAOD = cbk->sumOfEventWeightsSquared();
    }


    // Use the xAOD numbers from the largest cycle
    if (cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD" && cbk->cycle() > maxAcycle) {
      maxAcycle    = cbk->cycle();
      N_xAOD     = cbk->nAcceptedEvents();
      sumw_xAOD  = cbk->sumOfEventWeights();
      sumw2_xAOD = cbk->sumOfEventWeightsSquared();
    }
  }
  Info("fileExecute()",
       "Book keeper anticipates %i events in current input file (%i in parent xAOD)",N_DxAOD,N_xAOD);

  // Update the cutflow histos
  for (auto cf : m_cFlowHistos )
  {
    addBookKeeping(cf.second,N_xAOD,N_DxAOD);
  }

  // Update the weighted cutflow histos
  for (auto cf : m_cFlowHistosWeighted )
  {
    addBookKeeping(cf.second,sumw_xAOD,sumw_DxAOD,sumw2_xAOD,sumw2_DxAOD);
  }

  return EL::StatusCode::SUCCESS;
}

void MxAODTool::addBookKeeping(TH1F *cutflow, double N_xAOD, double N_DxAOD,
                               double sumw2_xAOD, double sumw2_DxAOD) {
  int bin_xAOD=cutflow->FindBin(NxAOD), bin_DxAOD=cutflow->FindBin(NDxAOD);
  cutflow->AddBinContent(bin_xAOD,N_xAOD);
  cutflow->AddBinContent(bin_DxAOD,N_DxAOD);
  if (sumw2_xAOD>0.0) {
    cutflow->SetBinError(bin_xAOD, sqrt(pow(cutflow->GetBinError(bin_xAOD),2)  + sumw2_xAOD ));
    cutflow->SetBinError(bin_DxAOD,sqrt(pow(cutflow->GetBinError(bin_DxAOD),2) + sumw2_DxAOD));
  }
  return;
}

void MxAODTool::printCutFlowHistos() {
  for ( auto entry : m_cFlowHistos ) {
    printf("\n%s %d cut-flow%s\n",HG::isMC()?"MC sample":"Data run",
           std::abs(entry.first),entry.first>0?", all events":", only Dalitz events");
    printCutFlowHisto(entry.second,0);
  }
  for ( auto entry : m_cFlowHistosWeighted ) {
    printf("\n%s %d cut-flow, weighted events%s\n",HG::isMC()?"MC sample":"Data run",
           std::abs(entry.first),entry.first>0?", all events":", only Dalitz events");
    printCutFlowHisto(entry.second,2);
  }
  printf("\n");
}

void MxAODTool::printCutFlowHisto(TH1F *h, int Ndecimals) {
  TString format("  %-24s%10."); format+=Ndecimals; format+="f%11.2f%%%11.2f%%\n";
  int all_bin = h->FindBin(ALLEVTS);
  printf("  %-24s%10s%12s%12s\n","Event selection","Nevents","Cut rej.","Tot. eff.");
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    double ntot=h->GetBinContent(all_bin), n=h->GetBinContent(bin), nprev=h->GetBinContent(bin-1);
    TString cutName(h->GetXaxis()->GetBinLabel(bin));
    cutName.ReplaceAll("#it{m}_{#gamma#gamma}","m_yy");
    if (bin==1||nprev==0||n==nprev)
      printf(format.Data(),cutName.Data(),n,-1e-10,n/ntot*100);
    else // if the cut does something, print more information
      printf(format.Data(),cutName.Data(),n,(n-nprev)/nprev*100,n/ntot*100);
  }
}

TH1F* MxAODTool::makeCutFlowHisto(int id, const std::vector<TString> otherCuts, TString suffix) {
  std::vector<TString> tmp_cutDescs;

  for (auto cut : s_cutDescs_base)
  { tmp_cutDescs.push_back(cut); }

  for (auto cut : otherCuts)
  { tmp_cutDescs.push_back(cut); }

  int Ncuts = tmp_cutDescs.size();

  // create meaningful name of the cutflow histo
  TString name(Form("CutFlow_%s%d",HG::isData()?"Run":"MC",std::abs(id)));

  bool hasMCname = HG::isMC() && config()->isDefined(Form("SampleName.%d",std::abs(id)));

  if(hasMCname){
    name = Form("CutFlow_%s",getMCSampleName(std::abs(id)).Data());
  }

  if (HG::isMC()&&!hasMCname&&config()->getStr("SampleName","sample")!="sample")
    name="CutFlow_"+config()->getStr("SampleName");
  name+=suffix;

  // maybe should make sure we don't switch directory?
  TDirectory *dir = gDirectory;
  TFile *file = wk()->getOutputFile("MxAOD");
  TH1F *h = new TH1F(name,name,Ncuts,0,Ncuts);
  h->SetDirectory(file); // probably not needed
  for (int bin=1;bin<=Ncuts;++bin)
    h->GetXaxis()->SetBinLabel(bin,tmp_cutDescs[bin-1]);
  dir->cd();
  return h;
}
