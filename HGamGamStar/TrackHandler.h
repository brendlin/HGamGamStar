#ifndef HGamGamStar_TrackHandler
#define HGamGamStar_TrackHandler

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaHandler.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "HGamGamStar/TrackElectronMap.h"

namespace HG {

  class TrackHandler : public HgammaHandler<xAOD::TrackParticle, xAOD::TrackParticleContainer, xAOD::TrackParticleAuxContainer> {

  private:

    int m_nSiMin;
    int m_nPixMin;

    int m_truth_nSiMin;

    double  m_etaCut;
    double  m_ptCut;

    double m_d0BySigd0Cut;
    double m_z0Cut;

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
                                                         const xAOD::ElectronContainer& elecs,
                                                         TrackElectronMap& trkEleMap,
                                                         bool doTruthClassify=false);

    xAOD::ElectronContainer GetElecsAssociatedToTracks(xAOD::TrackParticle& trk1,
                                                       xAOD::TrackParticle& trk2,
                                                       xAOD::ElectronContainer& preSelElecs);

    TruthTrackMap MakeTruthTrackMapFromGSFContainer(xAOD::TrackParticleContainer& tracks);
    TruthTrackMap MakeTruthTrackMapFromElectronContainer(const xAOD::ElectronContainer& elecs);

    bool passIPCuts(xAOD::TrackParticle& trk);

    float calculateIPSig(const xAOD::TrackParticle& trk) const;
    void decorateIPCut(xAOD::TrackParticle& trk);
    float calculateTRT_PID(const xAOD::TrackParticle& trk)const;

    void decorateAdditionalCuts(xAOD::TrackParticle& trk);

    size_t nMatchedElectrons(const xAOD::TrackParticle& trk) const;

  };

} // namespace HG

#endif // HGamGamStar_TrackHandler
