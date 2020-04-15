#include "HGamGamStar/HiggsGamGamStarCutflowAndMxAOD.h"
#include "HGamGamStar/MergedElectronMxAOD.h"

#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  // HiggsGamGamStarCutflowAndMxAOD *alg = new  HiggsGamGamStarCutflowAndMxAOD("HiggsGamGamStarCutflowAndMxAOD");
  auto *alg = new MergedElectronMxAOD("HiggsGamGamStarCutflowAndMxAOD");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
