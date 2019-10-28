#ifndef HGamGamStar_MergedElectronID
#define HGamGamStar_MergedElectronID

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


namespace HG {

  class MergedElectronID {

  private:
    
    unsigned getPtBin(const xAOD::Electron * const el) const;
    unsigned getEtaBin(const xAOD::Electron * const el) const;
    
    enum class extrapolationStartPositionEnum {
      Perigee,
      FirstMeasurement,
      LastMeasurement
    };

    // Preselection cuts
    int m_PreselNPassBlayer;
    float m_PreselRhad;
    float m_mergedElePtCut;
    float m_mergedEleEtaCut;
    
    AngularPosition getExtrapolatedTrackPosition(
	const xAOD::TrackParticle * track,
	const extrapolationStartPositionEnum extrapolationStartPosition,
	const bool kalmanUpdate,
	const bool verbose,
	const bool printTrajectory);
    



  public:

    /// constructor
    MergedElectronID();

    /// destructor
    virtual ~MergedElectronID();

    virtual EL::StatusCode initialize(Config &config);
    
    bool passPIDCut(const xAOD::Electron &ele,const xAOD::TrackParticle &trk1,const xAOD::TrackParticle &trk2);

    void decorateMergedVariables(xAOD::Electron &ele,
                                 xAOD::TrackParticle &trk1,
                                 xAOD::TrackParticle &trk2);

    bool passPreselection(const xAOD::Electron &ele,
                          const xAOD::TrackParticle &trk1,
                          const xAOD::TrackParticle &trk2);
    
  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID
