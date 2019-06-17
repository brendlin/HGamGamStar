#ifndef HGamGamStar_HggStarCommon_H
#define HGamGamStar_HggStarCommon_H

#include "AthContainers/AuxElement.h"
#include <vector>

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

    //Tracks matched to the electron variables
    static SG::AuxElement::Accessor< std::vector<int> >   passTTVA ("trackPassTTVA");
    static SG::AuxElement::Accessor< std::vector<float> > trackPT ("trackPT");
    static SG::AuxElement::Accessor< std::vector<float> > trackP ("trackP");
    static SG::AuxElement::Accessor< std::vector<float> > trackD0 ("trackD0");
    static SG::AuxElement::Accessor< std::vector<float> > trackD0Sig ("trackD0Sig");
    static SG::AuxElement::Accessor< std::vector<float> > trackZ0 ("trackZ0");
    static SG::AuxElement::Accessor< std::vector<float> > trackZ0Sig ("trackZ0Sig");
    static SG::AuxElement::Accessor< std::vector<float> > trackTRT_PID_trans("trackTRT_PID_trans");

    static SG::AuxElement::Accessor< std::vector<int> >   trackNPix("trackNPix");
    static SG::AuxElement::Accessor< std::vector<int> >   trackNSCT("trackNSCT");
    static SG::AuxElement::Accessor< std::vector<int> >   trackPassBL("trackPassBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackNBL("trackNBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackNIBL("trackNIBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackSplitBL("trackSplitBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackSharedBL("trackSharedBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackSplitIBL("trackSplitIBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackSharedIBL("trackSharedIBL");
    static SG::AuxElement::Accessor< std::vector<int> >   trackPdgID("trackPdgID");
    static SG::AuxElement::Accessor< std::vector<int> >   trackBarcode("trackBarcode");
    static SG::AuxElement::Accessor< std::vector<float> > trackTruthE("trackTruthE");
    static SG::AuxElement::Accessor< std::vector<int> >   trackFromHiggs("trackFromHiggs");
    static SG::AuxElement::Accessor< int >   isTrueMergedE("isTrueMergedE");
    static SG::AuxElement::Accessor< int >   trueType("trueType");
    static SG::AuxElement::Accessor< float > trueEnergy("trueEnergy");
    static SG::AuxElement::Accessor< float > trueMass("trueMass");
    static SG::AuxElement::Accessor< float > trueEta("trueEta");
    static SG::AuxElement::Accessor< float > truePhi("truePhi");
    static SG::AuxElement::Accessor< float > standAloneVertexR("standAloneVertexR");
    static SG::AuxElement::Accessor< int >   standAloneIndexA("standAloneIndexA");
    static SG::AuxElement::Accessor< int >   standAloneIndexB("standAloneIndexB");

    static SG::AuxElement::Accessor< float > ambiguousPhotonR("ambiguousPhotonR");
    static SG::AuxElement::Accessor< int >   ambiguousPhotonCT("ambiguousPhotonCT");
    static SG::AuxElement::Accessor< float > calibratedPhotonEnergy("calibratedPhotonEnergy");
    static SG::AuxElement::Accessor< float > calibratedElectronEnergy("calibratedElectronEnergy");


    static SG::AuxElement::Accessor< int >   truthTrackIndexA("truthTrackIndexA");
    static SG::AuxElement::Accessor< int >   truthTrackIndexB("truthTrackIndexB");


    static SG::AuxElement::Accessor< int >  vtxTrkIndex1("vtxTrkParticleIndex1");
    static SG::AuxElement::Accessor< int >  vtxTrkIndex2("vtxTrkParticleIndex2");


    static SG::AuxElement::Accessor< int >   vtxTrk1_d0("vtxTrk1_d0");
    static SG::AuxElement::Accessor< float > vtxTrk1_PT ("vtxTrk1_PT");
    static SG::AuxElement::Accessor< float > vtxTrk1_P ("vtxTrk1_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_D0 ("vtxTrk1_D0");
    static SG::AuxElement::Accessor< float > vtxTrk1_D0Sig ("vtxTrk1_D0Sig");
    static SG::AuxElement::Accessor< float > vtxTrk1_Z0 ("vtxTrk1_Z0");
    static SG::AuxElement::Accessor< float > vtxTrk1_Z0Sig ("vtxTrk1_Z0Sig");
    static SG::AuxElement::Accessor< float > vtxTrk1_TRT_PID_trans("vtxTrk1_TRT_PID_trans");
    static SG::AuxElement::Accessor< int >   vtxTrk1_NPix("vtxTrk1_NPix");
    static SG::AuxElement::Accessor< int >   vtxTrk1_NSCT("vtxTrk1_NSCT");
    static SG::AuxElement::Accessor< int >   vtxTrk1_PassBL("vtxTrk1_PassBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_NBL("vtxTrk1_NBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_NIBL("vtxTrk1_NIBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_SplitBL("vtxTrk1_SplitBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_SharedBL("vtxTrk1_SharedBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_SplitIBL("vtxTrk1_SplitIBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_SharedIBL("vtxTrk1_SharedIBL");
    static SG::AuxElement::Accessor< int >   vtxTrk1_PdgID("vtxTrk1_PdgID");
    static SG::AuxElement::Accessor< int >   vtxTrk1_Barcode("vtxTrk1_Barcode");
    static SG::AuxElement::Accessor< float > vtxTrk1_TruthE("vtxTrk1_TruthE");
    static SG::AuxElement::Accessor< int >   vtxTrk1_FromHiggs("vtxTrk1_FromHiggs");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta2_P("vtxTrk1_dEta2_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta1_P("vtxTrk1_dEta1_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_dPhi2_P("vtxTrk1_dPhi2_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta2_LM("vtxTrk1_dEta2_LM");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta1_LM("vtxTrk1_dEta1_LM");
    static SG::AuxElement::Accessor< float > vtxTrk1_dPhi2_LM("vtxTrk1_dPhi2_LM");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta2_T("vtxTrk1_dEta2_T");
    static SG::AuxElement::Accessor< float > vtxTrk1_dPhi2_T("vtxTrk1_dPhi2_T");

    static SG::AuxElement::Accessor< int >   vtxTrk2_d0("vtxTrk2_d0");
    static SG::AuxElement::Accessor< float > vtxTrk2_PT ("vtxTrk2_PT");
    static SG::AuxElement::Accessor< float > vtxTrk2_P ("vtxTrk2_P");
    static SG::AuxElement::Accessor< float > vtxTrk2_D0 ("vtxTrk2_D0");
    static SG::AuxElement::Accessor< float > vtxTrk2_D0Sig ("vtxTrk2_D0Sig");
    static SG::AuxElement::Accessor< float > vtxTrk2_Z0 ("vtxTrk2_Z0");
    static SG::AuxElement::Accessor< float > vtxTrk2_Z0Sig ("vtxTrk2_Z0Sig");
    static SG::AuxElement::Accessor< float > vtxTrk2_TRT_PID_trans("vtxTrk2_TRT_PID_trans");
    static SG::AuxElement::Accessor< int >   vtxTrk2_NPix("vtxTrk2_NPix");
    static SG::AuxElement::Accessor< int >   vtxTrk2_NSCT("vtxTrk2_NSCT");
    static SG::AuxElement::Accessor< int >   vtxTrk2_PassBL("vtxTrk2_PassBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_NBL("vtxTrk2_NBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_NIBL("vtxTrk2_NIBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_SplitBL("vtxTrk2_SplitBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_SharedBL("vtxTrk2_SharedBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_SplitIBL("vtxTrk2_SplitIBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_SharedIBL("vtxTrk2_SharedIBL");
    static SG::AuxElement::Accessor< int >   vtxTrk2_PdgID("vtxTrk2_PdgID");
    static SG::AuxElement::Accessor< int >   vtxTrk2_Barcode("vtxTrk2_Barcode");
    static SG::AuxElement::Accessor< float > vtxTrk2_TruthE("vtxTrk2_TruthE");
    static SG::AuxElement::Accessor< int >   vtxTrk2_FromHiggs("vtxTrk2_FromHiggs");
    static SG::AuxElement::Accessor< float > vtxTrk2_dEta2_P("vtxTrk2_dEta2_P");
    static SG::AuxElement::Accessor< float > vtxTrk2_dEta1_P("vtxTrk2_dEta1_P");
    static SG::AuxElement::Accessor< float > vtxTrk2_dPhi2_P("vtxTrk2_dPhi2_P");
    static SG::AuxElement::Accessor< float > vtxTrk2_dEta2_LM("vtxTrk2_dEta2_LM");
    static SG::AuxElement::Accessor< float > vtxTrk2_dEta1_LM("vtxTrk2_dEta1_LM");
    static SG::AuxElement::Accessor< float > vtxTrk2_dPhi2_LM("vtxTrk2_dPhi2_LM");
    static SG::AuxElement::Accessor< float > vtxTrk2_dEta2_T("vtxTrk2_dEta2_T");
    static SG::AuxElement::Accessor< float > vtxTrk2_dPhi2_T("vtxTrk2_dPhi2_T");

    // The "official" track-matching deltaR variables, e.g. propagating to the calorimeter from
    // the perigee (P) or last measurement (LM), to layer 1 (2)
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_LM_L1("dRbetweenTracks_LM_L1"); */
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_LM_L2("dRbetweenTracks_LM_L2"); */
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_P_L1("dRbetweenTracks_P_L1"); */
    /* static SG::AuxElement::Accessor<float> dRbetweenTracks_P_L2("dRbetweenTracks_P_L2"); */

    //
    // DAOD Accessors:
    //

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

    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingTrue_dEta2("TrackMatchingTrue_dEta2");
    static SG::AuxElement::Accessor<std::vector<float>> TrackMatchingTrue_dPhi2("TrackMatchingTrue_dPhi2");

    // Sub-cluster accessors
    static SG::AuxElement::Accessor<std::vector<float>> subCluster_E("SubCluster_E");
    static SG::AuxElement::Accessor<std::vector<float>> subCluster_dEta("SubCluster_dEta");
    static SG::AuxElement::Accessor<std::vector<float>> subCluster_dPhi("SubCluster_dPhi");

  } // namespace EleAcc

} // namespace HG

#endif // HGamGamStar_HggStarCommon_H
