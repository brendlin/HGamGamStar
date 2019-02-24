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
    static SG::AuxElement::Accessor<float> z0sinTheta("z0sinTheta");
    static SG::AuxElement::Accessor<char>  isTrueHiggsElectron("isTrueHiggsElectron");
    static SG::AuxElement::Accessor<float> TRT_PID_trans("TRT_PID_trans");
    static SG::AuxElement::Accessor<char>  passBLayerRequirement("passBLayerRequirement");
    static SG::AuxElement::Accessor<float> pt("pt");

  } // namespace TrkAcc

  namespace EleAcc {

    static SG::AuxElement::Accessor<float> RhadForPID("RhadForPID");
    static SG::AuxElement::Accessor<float> EOverP0P1("EOverP0P1");
    static SG::AuxElement::Accessor<float> dRExtrapTrk12("dRExtrapTrk12");

  } // namespace EleAcc

} // namespace HG

#endif // HGamGamStar_HggStarCommon_H
