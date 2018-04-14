#ifndef HGamGamStar_TrackHandler
#define HGamGamStar_TrackHandler

#include "HGamAnalysisFramework/HgammaHandler.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

namespace HG {

  class TrackHandler : public HgammaHandler<xAOD::TrackParticle, xAOD::TrackParticleContainer, xAOD::TrackParticleAuxContainer> {

  private:

    bool m_doTrqCuts;
    int m_nSiMin;
    int m_nPixMin;

    double  m_etaCut;
    double  m_ptCut;

    // pt, eta cuts on underlying electron
    double  m_eleEtaCut;
    double  m_elePtCut;

  public:

    /// constructor
    TrackHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);

    /// destructor
    virtual ~TrackHandler();

    virtual EL::StatusCode initialize(Config &config);

    virtual xAOD::TrackParticleContainer getCorrectedContainer();
    virtual xAOD::TrackParticleContainer applySelection(xAOD::TrackParticleContainer &container);
    virtual CP::SystematicCode    applySystematicVariation(const CP::SystematicSet &sys);

    xAOD::TrackParticleContainer findTracksFromElectrons(xAOD::TrackParticleContainer& container,
                                                         const xAOD::ElectronContainer& elecs);

    xAOD::ElectronContainer GetElecsAssociatedToTracks(const xAOD::TrackParticle& trk1,
                                                       const xAOD::TrackParticle& trk2,
                                                       xAOD::ElectronContainer& preSelElecs);

  };

} // namespace HG

#endif // HGamGamStar_TrackHandler
