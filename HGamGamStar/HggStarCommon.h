#ifndef HGamGamStar_HggStarCommon_H
#define HGamGamStar_HggStarCommon_H

#include "AthContainers/AuxElement.h"

namespace HG {

  enum ChannelEnum {
    CHANNELUNKNOWN=0,
    DIMUON=1,
    RESOLVED_DIELECTRON=2,
    MERGED_DIELECTRON=3,
    AMBIGUOUS_DIELECTRON=4,
    FAILEDTRKELECTRON=5, //// Tracking Failed Electron Decay
    OTHER=6, //// Other Decay
    OUT_OF_ACCEPTANCE=7
  };

  namespace TrkAcc {
    static SG::AuxElement::Accessor< std::vector<int> > MatchedElectrons("MatchedElectrons");
    static SG::AuxElement::Accessor<char>  passIPCut("passIPCut");
    static SG::AuxElement::Accessor<float> d0significance("d0significance");
    static SG::AuxElement::Accessor<float> z0pv("z0pv");
    static SG::AuxElement::Accessor<float> z0sinTheta("z0sinTheta");
    static SG::AuxElement::Accessor<char>  isTrueHiggsElectron("isTrueHiggsElectron");
    static SG::AuxElement::Accessor<float> TRT_PID_trans("TRT_PID_trans");
    static SG::AuxElement::Accessor<char>  passBLayerRequirement("passBLayerRequirement");
    static SG::AuxElement::Accessor<float> pt("pt");
    static SG::AuxElement::Accessor<int> mergedTrackParticleIndex("mergedTrackParticleIndex");

  } // namespace TrkAcc

  namespace EleAcc {

    //
    // New Accessors for output MxAOD:
    //

    // New ID variables
    static SG::AuxElement::Accessor<float> RhadForPID("RhadForPID");
    static SG::AuxElement::Accessor<float> EOverP0P1("EOverP0P1");
    static SG::AuxElement::Accessor<float> delta_z0_tracks("delta_z0_tracks");
    static SG::AuxElement::Accessor<float> delta_z0sinTheta_tracks("delta_z0sinTheta_tracks");

    // Locally calculated, using the perigee parameters, extrapolating to the calorimeter.
    static SG::AuxElement::Accessor<float> dRExtrapTrk12("dRExtrapTrk12");
    // Locally calculated, using the last measurement, extrapolating to the calorimeter.
    static SG::AuxElement::Accessor<float> dRExtrapTrk12_LM("dRExtrapTrk12_LM");

    // The "official" track-matching deltaR variables, e.g. propagating to the calorimeter from
    // the perigee (P) or last measurement (LM), to layer 1 (2)
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_LM_L1("dRbetweenTracks_LM_L1"); */
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_LM_L2("dRbetweenTracks_LM_L2"); */
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_P_L1("dRbetweenTracks_P_L1"); */
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_P_L2("dRbetweenTracks_P_L2"); */

    //
    // DAOD Accessors:
    //
    static SG::AuxElement::Accessor<uint8_t> ambiguityType("ambiguityType");
    static SG::AuxElement::Accessor<float> conversionRadius("conversionRadius");

    // Accessors for deltaEta / deltaPhi (perigee) - the index is the Nth-matched xAOD::Electron track
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingP_dEta1("TrackMatchingP_dEta1");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingP_dEta2("TrackMatchingP_dEta2");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingP_dPhi1("TrackMatchingP_dPhi1");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingP_dPhi2("TrackMatchingP_dPhi2");

    // Accessors for deltaEta / deltaPhi (last measurement) - the index is the Nth-matched xAOD::Electron track
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingLM_dEta1("TrackMatchingLM_dEta1");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingLM_dEta2("TrackMatchingLM_dEta2");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingLM_dPhi1("TrackMatchingLM_dPhi1");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingLM_dPhi2("TrackMatchingLM_dPhi2");

    // Sub-cluster accessors
    static SG::AuxElement::Accessor<std::vector<float>> subCluster_E("SubCluster_E");
    static SG::AuxElement::Accessor<std::vector<float>> subCluster_dEta("SubCluster_dEta");
    static SG::AuxElement::Accessor<std::vector<float>> subCluster_dPhi("SubCluster_dPhi");

  } // namespace EleAcc

} // namespace HG

#endif // HGamGamStar_HggStarCommon_H
