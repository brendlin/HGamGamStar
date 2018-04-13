#ifndef HGamGamStar_TrackHandler
#define HGamGamStar_TrackHandler

#include "HGamAnalysisFramework/HgammaHandler.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

namespace HG {

  class TrackHandler : public HgammaHandler<xAOD::TrackParticle, xAOD::TrackParticleContainer, xAOD::TrackParticleAuxContainer> {

  private:

  public:

    /// constructor
    TrackHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);

    /// destructor
    virtual ~TrackHandler();

    virtual EL::StatusCode initialize(Config &config);

  };

} // namespace HG

#endif // HGamGamStar_TrackHandler
