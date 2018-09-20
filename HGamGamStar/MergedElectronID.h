#ifndef HGamGamStar_MergedElectronID
#define HGamGamStar_MergedElectronID

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamGamStar/HggStarVariables.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"


namespace HG {

  class MergedElectronID {

  private:

    bool m_doTrqCuts;
    int m_nSiMin;
    int m_nPixMin;

    double  m_etaCut;
    double  m_ptCut;

  public:

    /// constructor
    MergedElectronID();

    /// destructor
    virtual ~MergedElectronID();

    virtual EL::StatusCode initialize(Config &config);
    
    bool passPIDCut(xAOD::Electron *ele,xAOD::TrackParticle *trk1,xAOD::TrackParticle *trk2);

    xAOD::TrackParticle m_trk1;
    xAOD::TrackParticle m_trk2;
    xAOD::Electron m_ele;

  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID
