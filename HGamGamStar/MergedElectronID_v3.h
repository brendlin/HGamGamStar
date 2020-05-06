#ifndef HGamGamStar_MergedElectronID_v3
#define HGamGamStar_MergedElectronID_v3

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/AngularPosition.h"
#include "HGamGamStar/TrackModel.h"
#include "HGamGamStar/TrackHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"
#include <string>
#include <TH2.h>

#include "HGamGamStar/MergedElectronID_v2.h"
//enum MergedSystematic {
//  MERGEDUNC_NOMINAL=0,
//  MERGEDUNC_STAT_UP=1,
//  MERGEDUNC_STAT_DOWN=2,
//  MERGEDUNC_SYST_UP=3,
//  MERGEDUNC_SYST_DOWN=4
//};

namespace HG {

  class MergedElectronID_v3 {

  private:

    // Preselection cuts
    // These are applied via the original MergedElectronID passPreselection function!
    /* int m_PreselNPassBlayer; */
    /* float m_PreselRhad; */
    /* flaot m_mergedElePtCut; */

    bool m_isMC;
    std::string m_deltaEtaFileName;
    std::string m_deltaEtaHistName;
    TH2* m_sdetaCorr;

    float correctDeta1(const xAOD::Electron &ele, bool isMC) const;
  

  public:

    /// constructor
    MergedElectronID_v3();

    /// destructor
    virtual ~MergedElectronID_v3();

    virtual EL::StatusCode initialize(Config &config);

    bool passPIDCut(const xAOD::Electron &ele, bool isMC) const ;
    float GetScaleFactor(const xAOD::Electron &ele, MergedSystematic sys=MERGEDUNC_NOMINAL) const;

  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID_v3
