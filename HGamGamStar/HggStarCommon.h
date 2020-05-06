#ifndef HGamGamStar_HggStarCommon_H
#define HGamGamStar_HggStarCommon_H


#include "AthContainers/AuxElement.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTracking/VertexContainer.h"

#include "HGamGamStar/TrackElectronMap.h"

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

  TString GetChannelName(ChannelEnum channel);

  enum CategoryEnum {
    CATEGORYUNKNOWN=0,
    GGF_DIMUON=1,
    GGF_RESOLVED_DIELECTRON=2,
    GGF_MERGED_DIELECTRON=3,
    VBF_DIMUON=4,
    VBF_RESOLVED_DIELECTRON=5,
    VBF_MERGED_DIELECTRON=6,
    HIPTT_DIMUON=7,
    HIPTT_RESOLVED_DIELECTRON=8,
    HIPTT_MERGED_DIELECTRON=9,
  };

  namespace PhAcc {
    static SG::AuxElement::Accessor<float> RhadForPID("RhadForPID");
    static SG::AuxElement::Accessor<float> ambiguousE_deltaEta1("ambiguousE_deltaEta1");
    static SG::AuxElement::Accessor<float> vtxE("vtxE");
    static SG::AuxElement::Accessor< float > vtxTrk1_TRT_PID_trans("vtxTrk1_TRT_PID_trans");
    static SG::AuxElement::Accessor< float > vtxTrk2_TRT_PID_trans("vtxTrk2_TRT_PID_trans");
  }

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

    static SG::AuxElement::Accessor<float> p ("p");
    static SG::AuxElement::Accessor<int> nPixHitsAndDeadSens("nPixHitsAndDeadSens");
    static SG::AuxElement::Accessor<int> nSCTHitsAndDeadSens("nSCTHitsAndDeadSens");
    static SG::AuxElement::Accessor<int> PassTTVA ("PassTTVA");

    static SG::AuxElement::Accessor<int> NBL("NBL");
    static SG::AuxElement::Accessor<int> NIBL("NIBL");
    static SG::AuxElement::Accessor<int> SplitBL("SplitBL");
    static SG::AuxElement::Accessor<int> SharedBL("SharedBL");
    static SG::AuxElement::Accessor<int> SplitIBL("SplitIBL");
    static SG::AuxElement::Accessor<int> SharedIBL("SharedIBL");

    static SG::AuxElement::Accessor<int> PdgID("PdgID");
    static SG::AuxElement::Accessor<int> Barcode("Barcode");
    static SG::AuxElement::Accessor<float> TruthE("TruthE");

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
    static SG::AuxElement::Accessor< int >   isTrueMergedE("isTrueMergedE");
    static SG::AuxElement::Accessor< int >   trueType("trueType");
    static SG::AuxElement::Accessor< float > trueEnergy("trueEnergy");
    static SG::AuxElement::Accessor< float > trueMass("trueMass");
    static SG::AuxElement::Accessor< float > trueEta("trueEta");
    static SG::AuxElement::Accessor< float > truePhi("truePhi");
    static SG::AuxElement::Accessor< float > standAloneVertexR("standAloneVertexR");
    static SG::AuxElement::Accessor< int >   standAloneIndexA("standAloneIndexA");
    static SG::AuxElement::Accessor< int >   standAloneIndexB("standAloneIndexB");

    // Variables to connect the output MxAOD TrackParticle container with the electron container
    static SG::AuxElement::Accessor< std::vector<int> >  trkParticleIndex_MxAOD("trkParticleIndex_MxAOD");
    static SG::AuxElement::Accessor< int >  vtxTrkParticleIndex1_MxAOD("vtxTrkParticleIndex1_MxAOD");
    static SG::AuxElement::Accessor< int >  vtxTrkParticleIndex2_MxAOD("vtxTrkParticleIndex2_MxAOD");

    static SG::AuxElement::Accessor< float > ambiguousPhotonR("ambiguousPhotonR");
    static SG::AuxElement::Accessor< int >   ambiguousPhotonCT("ambiguousPhotonCT");
    static SG::AuxElement::Accessor< float > calibratedPhotonEnergy("calibratedPhotonEnergy");
    static SG::AuxElement::Accessor< float > calibratedPhotonEnergy50("calibratedPhotonEnergy50");
    static SG::AuxElement::Accessor< float > calibratedPhotonEnergy100("calibratedPhotonEnergy100");
    static SG::AuxElement::Accessor< float > calibratedPhotonEnergy200("calibratedPhotonEnergy200");
    static SG::AuxElement::Accessor< float > calibratedPhotonEnergy400("calibratedPhotonEnergy400");
    static SG::AuxElement::Accessor< float > calibratedElectronEnergy("calibratedElectronEnergy");


    static SG::AuxElement::Accessor< int >   truthTrackIndexA("truthTrackIndexA");
    static SG::AuxElement::Accessor< int >   truthTrackIndexB("truthTrackIndexB");


    static SG::AuxElement::Accessor<int>     passTMVAPIDv2("passTMVAPIDv2") ;
    static SG::AuxElement::Accessor<int>     passTMVAPID("passTMVAPID") ;
    static SG::AuxElement::Accessor<int>     passPID("passPID") ;

    static SG::AuxElement::Accessor< float > vtxTrk1_TRT_PID_trans("vtxTrk1_TRT_PID_trans");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta2_P("vtxTrk1_dEta2_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta1_P("vtxTrk1_dEta1_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_dPhi2_P("vtxTrk1_dPhi2_P");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta2_LM("vtxTrk1_dEta2_LM");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta1_LM("vtxTrk1_dEta1_LM");
    static SG::AuxElement::Accessor< float > vtxTrk1_dPhi2_LM("vtxTrk1_dPhi2_LM");
    static SG::AuxElement::Accessor< float > vtxTrk1_dEta2_T("vtxTrk1_dEta2_T");
    static SG::AuxElement::Accessor< float > vtxTrk1_dPhi2_T("vtxTrk1_dPhi2_T");

    static SG::AuxElement::Accessor< float > vtxTrk2_TRT_PID_trans("vtxTrk2_TRT_PID_trans");
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

    // Stored in the DAOD
    static SG::AuxElement::Accessor< int >  vtxTrkIndex1("vtxTrkParticleIndex1");
    static SG::AuxElement::Accessor< int >  vtxTrkIndex2("vtxTrkParticleIndex2");

    static SG::AuxElement::Accessor<float>   vtxdEta("vtxdEta") ;
    static SG::AuxElement::Accessor<float>   vtxdPhi("vtxdPhi") ;
    static SG::AuxElement::Accessor<float>   vtxPhi("vtxPhi") ;
    static SG::AuxElement::Accessor<float>   vtxEta("vtxEta") ;
    static SG::AuxElement::Accessor<float>   vtxZ("vtxZ") ;
    static SG::AuxElement::Accessor<float>   vtxR("vtxR") ;
    static SG::AuxElement::Accessor<float>   vtxE("vtxE") ;
    static SG::AuxElement::Accessor<float>   vtxM("vtxM") ;

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
  
  namespace MuAcc {

    //
    // New Accessors for output MxAOD:
    //
    static SG::AuxElement::Accessor<float> d0significance("d0significance");
    
  } // namespace MuAcc

  xAOD::Photon*  createPhotonFromElectron (const xAOD::Electron* el) ;
  void setPhotonConversionVertex( const xAOD::Electron* el,
                             xAOD::Photon* ph, float vtxR,
                             xAOD::VertexContainer* vertexContainer ) ;

  ChannelEnum truthChannel(const xAOD::TruthParticleContainer& childleps,
                           const xAOD::ElectronContainer& all_elecs);

  ChannelEnum ClassifyElectronChannelsByBestMatch(const xAOD::TrackParticle* trk0,
                                                  const xAOD::TrackParticle* trk1,
                                                  const HG::TrackElectronMap& trkEleMap,
                                                  xAOD::ElectronContainer* inEleCont=nullptr,
                                                  xAOD::ElectronContainer* outEleCont=nullptr);

} // namespace HG

#endif // HGamGamStar_HggStarCommon_H
