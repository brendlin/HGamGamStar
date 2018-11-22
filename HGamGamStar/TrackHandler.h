#ifndef HGamGamStar_TrackHandler
#define HGamGamStar_TrackHandler

#include "HGamAnalysisFramework/HgammaHandler.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "HGamGamStar/TrackElectronMap.h"

namespace HG {

  class TrackHandler : public HgammaHandler<xAOD::TrackParticle, xAOD::TrackParticleContainer, xAOD::TrackParticleAuxContainer> {

  private:

    bool m_doTrqCuts;
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
    void decorateIPCut(xAOD::TrackParticle& trk);
    void decorateTRT_PID(xAOD::TrackParticle& trk);

    static SG::AuxElement::Accessor< std::vector<int> > MatchedElectrons;
    static SG::AuxElement::Accessor<char>  passIPCut;
    static SG::AuxElement::Accessor<float>  d0significance;
    static SG::AuxElement::Accessor<float>  z0sinTheta;
    static SG::AuxElement::Accessor<char>  isTrueHiggsElectron;
    static SG::AuxElement::Accessor<float> TRT_PID_trans;

    size_t nMatchedElectrons(const xAOD::TrackParticle& trk) const;

  };

} // namespace HG

#endif // HGamGamStar_TrackHandler
