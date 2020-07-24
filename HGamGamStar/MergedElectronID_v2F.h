#ifndef HGamGamStar_MergedElectronID_v2FF
#define HGamGamStar_MergedElectronID_v2FF

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/AngularPosition.h"
#include "HGamGamStar/TrackModel.h"
#include "HGamGamStar/TrackHandler.h"
#include "HGamGamStar/MergedElectronID_v2.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"
#include <string>
#include <TH2.h>


namespace HG {
  //* The F is for Fudge :D*/
  
  class MergedElectronID_v2F {

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
    MergedElectronID_v2F();

    /// destructor
    virtual ~MergedElectronID_v2F();

    virtual EL::StatusCode initialize(Config &config);

    int getEtaBin( float eta ) const;
    int getPtBin( float pt ) const;

    bool passPIDCut(const xAOD::Electron &ele, bool isMC) const ;
    float GetScaleFactor(const xAOD::Electron &ele, MergedSystematic sys=MERGEDUNC_NOMINAL) const;

  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID_v2FF
