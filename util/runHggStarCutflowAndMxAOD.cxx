#include "HGamGamStar/HggStarCutflowAndMxAOD.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  HggStarCutflowAndMxAOD *alg = new HggStarCutflowAndMxAOD("HggStarCutflowAndMxAOD");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
