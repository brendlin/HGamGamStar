#ifndef HGamGamStar_ExtraHggStarObjects_H
#define HGamGamStar_ExtraHggStarObjects_H

//
// This class is simply meant to save some extra, easily accessible variables
// so they do not have to be calculated again (such as the selected tracks, or the
// true uh stuff.
//

#include "xAODBase/IParticleContainer.h"

namespace HG {

  class ExtraHggStarObjects {
  private:
    static ExtraHggStarObjects *m_ptr;

    xAOD::IParticleContainer m_tracks;
    bool m_tracksAvail;

  public:
    /// Get instance of singleton class
    static ExtraHggStarObjects *getInstance();

    void setElectronTrackContainer(const xAOD::IParticleContainer* tracks);

    /// Reset containers to null pointers to avoid carry-over from previous event
    void clearContainers();

    /// Get pointer to collection
    const xAOD::IParticleContainer *getElectronTracks(bool truth = false) const;

  private:
    /// Default constructor - note that it is private
    ExtraHggStarObjects();

    /// Default detructor - note that is is private
    ~ExtraHggStarObjects();

  }; // class ExtraHggStarObjects

}

#endif // HGamGamStar_ExtraHggStarObjects_H
