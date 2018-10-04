#ifndef HGamGamStar_MergedElectronID
#define HGamGamStar_MergedElectronID

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/AngularPosition.h"
#include "HGamGamStar/TrackModel.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"


namespace HG {

  class MergedElectronID {

  private:
    
    bool passCut(const float obsValue, const std::string cutString);
    unsigned getPtBin(const xAOD::Electron * const el) const;
    unsigned getEtaBin(const xAOD::Electron * const el) const;
    
    enum class extrapolationStartPositionEnum {
      Perigee,
      FirstMeasurement,
      LastMeasurement
    };
    
    AngularPosition getExtrapolatedTrackPosition(
	const xAOD::TrackParticle * track,
	const extrapolationStartPositionEnum extrapolationStartPosition,
	const bool kalmanUpdate,
	const bool verbose,
	const bool printTrajectory);
    
    extrapolationStartPositionEnum m_electron_trk_ex_origin;



  public:

    /// constructor
    MergedElectronID();

    /// destructor
    virtual ~MergedElectronID();

    virtual EL::StatusCode initialize(Config &config);
    
    bool passPIDCut(xAOD::Electron &ele,xAOD::TrackParticle &trk1,xAOD::TrackParticle &trk2);
    
  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID
