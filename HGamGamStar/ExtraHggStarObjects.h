#ifndef HGamGamStar_ExtraHggStarObjects_H
#define HGamGamStar_ExtraHggStarObjects_H

//
// This class is simply meant to save some extra, easily accessible variables
// so they do not have to be calculated again (such as the selected tracks, or the
// true uh stuff.
//

#include "xAODBase/IParticleContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTracking/TrackParticleAuxContainer.h"
#include "xAODEgamma/ElectronContainer.h"

namespace xAOD {
  class TEvent;
  class TStore;
}

namespace HG {

  class ExtraHggStarObjects {
  private:
    static ExtraHggStarObjects *m_ptr;

    xAOD::TEvent *m_event;
    xAOD::TStore *m_store;

    xAOD::IParticleContainer m_tracks;
    bool m_tracksAvail;

    xAOD::TruthParticleContainer m_higgsLeps;
    bool m_higgsLepsAvail;

    xAOD::TruthParticleContainer m_higgsPhotons;
    bool m_higgsPhotonsAvail;

    xAOD::TruthParticleContainer m_fsrPhotons;
    bool m_fsrPhotonsAvail;

  public:
    /// Get instance of singleton class
    static ExtraHggStarObjects *getInstance();

    void setElectronTrackContainer(const xAOD::IParticleContainer* tracks);

    /// Run the code to find the photon and leptons from the Higgs
    void setTruthHiggsDecayProducts(const xAOD::TruthParticleContainer* all_particles);

    /// Reset containers to null pointers to avoid carry-over from previous event
    void clearContainers();

    /// Get pointer to collection
    const xAOD::IParticleContainer *getElectronTracks(bool truth = false) const;

    const xAOD::TruthParticleContainer *getTruthHiggsLeptons() const;
    const xAOD::TruthParticleContainer *getTruthHiggsPhotons() const;
    const xAOD::TruthParticleContainer *getTruthHiggsFSRPhotons() const;

    /// Set TEvent and TStore
    void setEventAndStore(xAOD::TEvent *event, xAOD::TStore *store);

  private:
    /// Default constructor - note that it is private
    ExtraHggStarObjects();

    /// Default detructor - note that is is private
    ~ExtraHggStarObjects();

  }; // class ExtraHggStarObjects

}

#endif // HGamGamStar_ExtraHggStarObjects_H
