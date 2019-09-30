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

    // Perform the index-based track selection, as in the DAOD
    bool m_doIndexBasedTrackSelection;

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
                                                         xAOD::ElectronContainer& elecs,
                                                         TrackElectronMap& trkEleMap);

    TruthTrackMap MakeTruthTrackMapFromGSFContainer(xAOD::TrackParticleContainer& tracks);
    TruthTrackMap MakeTruthTrackMapFromElectronContainer(const xAOD::ElectronContainer& elecs);

    bool passTrackPreselection(const xAOD::TrackParticle* trk,bool doTruthClassify) const;
    void fillMergedIndices(xAOD::Electron* ele,int& tmp_vtxTrkIndex1,int& tmp_vtxTrkIndex2) const;

    bool passIPCuts(xAOD::TrackParticle& trk);

    float calculateIPSig(const xAOD::TrackParticle& trk) const;
    void decorateIPCut(xAOD::TrackParticle& trk);
    float calculateTRT_PID(const xAOD::TrackParticle& trk)const;

    void decorateAdditionalCuts(xAOD::TrackParticle& trk);

    size_t nMatchedElectrons(const xAOD::TrackParticle& trk) const;

  };

} // namespace HG

#endif // HGamGamStar_TrackHandler
