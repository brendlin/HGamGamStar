#include "HGamGamStar/MergedElectronID_v2.h"

//______________________________________________________________________________
HG::MergedElectronID_v2::MergedElectronID_v2()
{

}

//______________________________________________________________________________
HG::MergedElectronID_v2::~MergedElectronID_v2()
{

}

//______________________________________________________________________________

EL::StatusCode HG::MergedElectronID_v2::initialize(Config &config)
{
  m_PreselNPassBlayer = config.getInt("MergedElectrons.Preselection.NtracksPassingBlayer",1);
  m_PreselRhad        = config.getNum("MergedElectrons.Preselection.RhadMin",0.10);
  m_mergedElePtCut    = config.getNum("MergedElectrons.Selection.PtPreCutGeV",20.) * GeV;

  return EL::StatusCode::SUCCESS;
}


//______________________________________________________________________________
bool HG::MergedElectronID_v2::passPIDCut(const xAOD::Electron &ele) const{


    int trk1_index =  HG::EleAcc::vtxTrkIndex1(ele);
    int trk2_index =  HG::EleAcc::vtxTrkIndex2(ele);

    // Two tracks are found
    if(trk1_index < 0 || trk2_index < 0)
      return false;

    // Vertex Fit failed
    if( HG::EleAcc::vtxE(ele) < 0)
      return false;

    // calculate shower shapes and other discriminating variables
    float trk_dEta1 = ele.trackCaloMatchValue(xAOD::EgammaParameters::TrackCaloMatchType::deltaEta1);

//    double f1 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::f1);
//     double fSide = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::fracs1);

    auto ele_tp  =  ele.trackParticle(trk1_index);
    auto ele_tp2 =  ele.trackParticle(trk2_index);

    uint8_t shared;
    int sumOfSplitHits = 0;
    if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerSharedHits) )
      sumOfSplitHits += (int)shared;
    if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerSharedHits) )
      sumOfSplitHits += (int)shared;
    if( ele_tp->summaryValue(shared, xAOD::numberOfInnermostPixelLayerSplitHits ) )
      sumOfSplitHits += (int)shared;
    if( ele_tp->summaryValue(shared, xAOD::numberOfNextToInnermostPixelLayerSplitHits) )
      sumOfSplitHits += (int)shared;

    if(sumOfSplitHits == 0){
      if( fabs( ele_tp2->z0() - ele_tp->z0() ) >  1.0 )
        return false;
    }


    float vtx_deta  = HG::EleAcc::vtxdEta(ele);
    float vtx_dphi  = HG::EleAcc::vtxdPhi(ele);

    float pte  = ele.pt();
    float etae = ele.eta();

    float Rhad = HG::EleAcc::RhadForPID(ele);
    float trk_TRT_PID1 = HG::EleAcc::vtxTrk1_TRT_PID_trans(ele);
    float trk_TRT_PID2 = HG::EleAcc::vtxTrk2_TRT_PID_trans(ele);
    float EoP =  HG::EleAcc::calibratedPhotonEnergy(ele) / HG::EleAcc::vtxE(ele);


    float Eratio = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Eratio);
    float wTotS1 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::wtots1);

    float Reta = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Reta);
    float Rphi = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rphi);
    float wEta2 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::weta2);
    float f3 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::f3);





    // 60-800 GeV - at 60% WP

     //signal efficiency: 20180/32870=0.613934
     //bkg efficiency: 84/101448=0.00082801
     if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.007920 && Eratio>0.825656 && Rphi>0.685707 && Reta>0.945474 && f3<0.098034 && wEta2<0.010545 && wTotS1<4.803770 && trk_TRT_PID1>-0.048709 && trk_TRT_PID2>-0.271082 && EoP<1.996022 && trk_dEta1<0.003580)) return true;

     //signal efficiency: 12242/19675=0.622211
     //bkg efficiency: 80/70384=0.00113662
     if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.010288 && Eratio>0.920766 && Rphi>0.660021 && Reta>0.935034 && f3<0.023042 && wEta2<0.012633 && wTotS1<3.196767 && trk_TRT_PID1>-0.154170 && trk_TRT_PID2>-0.132265 && EoP<6.719277 && trk_dEta1<0.001038)) return true;

     //signal efficiency: 6552/10641=0.615732
     //bkg efficiency: 95/54701=0.00173671
    // if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.014723 && Eratio>0.918744 && Rphi>0.662869 && Reta>0.934536 && f3<0.033993 && wEta2<0.014340 && wTotS1<3.331633 && trk_TRT_PID1>-0.126811 && trk_TRT_PID2>-0.490945 && EoP<2.980618 && trk_dEta1<0.002115)) return true;
    //signal efficiency: 5280/10641=0.496194
    //bkg efficiency: 63/54701=0.00115172
    if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.006657 && Eratio>0.934981 && Rphi>0.792450 && Reta>0.929411 && f3<0.017305 && wEta2<0.012000 && wTotS1<3.288494 && trk_TRT_PID1>-0.215961 && trk_TRT_PID2>-0.195975 && EoP<1.964092 && trk_dEta1<0.003355)) return true;

     //signal efficiency: 4095/6649=0.615882
     //bkg efficiency: 113/40809=0.002769
    // if((wTotS1>-100. && pte>60.000000 && pte<800.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.007715 && Eratio>0.911347 && Rphi>0.591312 && Reta>0.920665 && f3<0.037214 && wEta2<0.012723 && wTotS1<1.402248 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<2.475802 && trk_dEta1<0.004850)) return true;
    //signal efficiency: 3437/6649=0.51692
    //bkg efficiency: 87/40809=0.00213188
    if((wTotS1>-100. && pte>60.000000 && pte<8000.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.005030 && Eratio>0.902925 && Rphi>0.742276 && Reta>0.928840 && f3<0.035636 && wEta2<0.015539 && wTotS1<1.783990 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<6.231330 && trk_dEta1<0.001250)) return true;


    // 50-60 GeV - at 55% WP

     //signal efficiency: 8009/14562=0.549993
     //bkg efficiency: 53/77572=0.000683236
     if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.008654 && Eratio>0.895195 && Rphi>0.609350 && Reta>0.941888 && f3<0.042226 && wEta2<0.016481 && wTotS1<2.500456 && trk_TRT_PID1>-0.291272 && trk_TRT_PID2>-0.087112 && EoP<1.724084 && trk_dEta1<0.001831)) return true;

     //signal efficiency: 5765/9928=0.580681
     //bkg efficiency: 78/53404=0.00146056
     if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.009434 && Eratio>0.827477 && Rphi>0.516148 && Reta>0.929478 && f3<0.187417 && wEta2<0.011540 && wTotS1<4.121699 && trk_TRT_PID1>-0.002691 && trk_TRT_PID2>-0.350279 && EoP<3.581217 && trk_dEta1<0.001494)) return true;

     //signal efficiency: 3266/5885=0.55497
     //bkg efficiency: 67/42756=0.00156703
    // if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.007913 && Eratio>0.465984 && Rphi>0.592916 && Reta>0.928148 && f3<0.049200 && wEta2<0.013342 && wTotS1<3.791015 && trk_TRT_PID1>-0.072917 && trk_TRT_PID2>-0.485202 && EoP<2.708062 && trk_dEta1<0.002058)) return true;
    //signal efficiency: 2786/5885=0.473407
    //bkg efficiency: 55/42756=0.00128637
    if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.008794 && Eratio>0.806594 && Rphi>0.772923 && Reta>0.927457 && f3<0.059607 && wEta2<0.011405 && wTotS1<4.071177 && trk_TRT_PID1>-0.027964 && trk_TRT_PID2>-0.192150 && EoP<2.530357 && trk_dEta1<0.002562)) return true;

     //signal efficiency: 2215/4011=0.552231
     //bkg efficiency: 88/32942=0.00267136
    // if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.006816 && Eratio>0.943300 && Rphi>0.716474 && Reta>0.921586 && f3<0.036596 && wEta2<0.018066 && wTotS1<1.819416 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<2.179533 && trk_dEta1<0.007114)) return true;
    //signal efficiency: 1992/4011=0.496634
    //bkg efficiency: 70/32942=0.00212495
    if((wTotS1>-100. && pte>50.000000 && pte<60.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.009297 && Eratio>0.514015 && Rphi>0.847031 && Reta>0.927671 && f3<0.077539 && wEta2<0.012936 && wTotS1<2.191503 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<2.105054 && trk_dEta1<0.003560)) return true;


    // 40-50 GeV

     //signal efficiency: 4665/8951=0.521171
     //bkg efficiency: 177/486142=0.000364091
     if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.008006 && Eratio>0.866114 && Rphi>0.362006 && Reta>0.931920 && f3<0.010005 && wEta2<0.011536 && wTotS1<2.842709 && trk_TRT_PID1>-0.068432 && trk_TRT_PID2>-0.154553 && EoP<1.455041 && trk_dEta1<0.004640)) return true;

     //signal efficiency: 3288/6610=0.497428
     //bkg efficiency: 153/337512=0.000453317
     if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.011777 && Eratio>0.806071 && Rphi>0.508373 && Reta>0.936448 && f3<0.011271 && wEta2<0.011727 && wTotS1<3.917183 && trk_TRT_PID1>-0.081285 && trk_TRT_PID2>-0.159510 && EoP<5.202778 && trk_dEta1<0.000985)) return true;

     //signal efficiency: 2299/4254=0.540433
     //bkg efficiency: 205/264079=0.000776283
    // if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.011242 && Eratio>0.689382 && Rphi>0.671466 && Reta>0.921320 && f3<0.031406 && wEta2<0.013469 && wTotS1<3.021417 && trk_TRT_PID1>-0.119441 && trk_TRT_PID2>-0.447686 && EoP<1.975685 && trk_dEta1<0.002366)) return true;
    //signal efficiency: 1761/4254=0.413963
    //bkg efficiency: 125/264079=0.000473343
    if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.011486 && Eratio>0.927473 && Rphi>0.435598 && Reta>0.925742 && f3<0.010574 && wEta2<0.012369 && wTotS1<3.661756 && trk_TRT_PID1>-0.011563 && trk_TRT_PID2>-0.417202 && EoP<5.788745 && trk_dEta1<0.001250)) return true;

     //signal efficiency: 1614/2940=0.54898
     //bkg efficiency: 165/200137=0.000824435
    // if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.009445 && Eratio>0.855301 && Rphi>0.832559 && Reta>0.915871 && f3<0.028457 && wEta2<0.012498 && wTotS1<1.668140 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<3.852275 && trk_dEta1<0.001508)) return true;
    //signal efficiency: 1231/2940=0.418707
    //bkg efficiency: 80/200137=0.000399726
    if((wTotS1>-100. && pte>40.000000 && pte<50.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.007800 && Eratio>0.960342 && Rphi>0.743595 && Reta>0.920413 && f3<0.038835 && wEta2<0.012347 && wTotS1<1.411089 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<2.369355 && trk_dEta1<0.001652)) return true;


    // 30-40 GeV
    // changed to 45% WP

     //signal efficiency: 1992/4389=0.453862
     //bkg efficiency: 93/80643=0.00115323
     if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.018203 && Eratio>0.750413 && Rphi>0.378545 && Reta>0.939417 && f3<0.007743 && wEta2<0.013086 && wTotS1<2.747719 && trk_TRT_PID1>-0.071126 && trk_TRT_PID2>-0.118842 && EoP<17.921568 && trk_dEta1<0.001160)) return true;

     //signal efficiency: 1537/3508=0.438141
     //bkg efficiency: 85/56006=0.00151769
     if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.010903 && Eratio>0.820701 && Rphi>0.395608 && Reta>0.898833 && f3<0.009788 && wEta2<0.012040 && wTotS1<2.792041 && trk_TRT_PID1>-0.099147 && trk_TRT_PID2>-0.049472 && EoP<14.419681 && trk_dEta1<0.001497)) return true;

     //signal efficiency: 1067/2490=0.428514
     //bkg efficiency: 114/46102=0.00247278
    // if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.009782 && Eratio>0.494125 && Rphi>0.578891 && Reta>0.929619 && f3<0.021705 && wEta2<0.012530 && wTotS1<5.977871 && trk_TRT_PID1>-0.072906 && trk_TRT_PID2>-0.283956 && EoP<3.785885 && trk_dEta1<0.001317)) return true;
    //signal efficiency: 948/2490=0.380723
    //bkg efficiency: 72/46102=0.00156175
    if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.011385 && Eratio>0.890966 && Rphi>0.587361 && Reta>0.928949 && f3<0.076013 && wEta2<0.012791 && wTotS1<5.512513 && trk_TRT_PID1>-0.088572 && trk_TRT_PID2>-0.074673 && EoP<2.245576 && trk_dEta1<0.002411)) return true;

     //signal efficiency: 839/1872=0.448184
     //bkg efficiency: 113/36215=0.00312025
    // if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.009827 && Eratio>0.286698 && Rphi>0.690578 && Reta>0.924369 && f3<0.032053 && wEta2<0.013278 && wTotS1<1.627911 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<4.401609 && trk_dEta1<0.001417)) return true;
    //signal efficiency: 715/1872=0.381944
    //bkg efficiency: 74/36215=0.00204335
    if((wTotS1>-100. && pte>30.000000 && pte<40.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.020000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.008061 && Eratio>0.871364 && Rphi>0.634550 && Reta>0.916855 && f3<0.034073 && wEta2<0.012036 && wTotS1<1.487706 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<2.079409 && trk_dEta1<0.001821)) return true;


    // 20-30 GeV
    // changed to 40% WP - smoothing the pt spectrum

     //signal efficiency: 782/2053=0.380906
     //bkg efficiency: 280/264932=0.00105687
     if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>0.000000 && fabs(etae)<0.800000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.010594 && Eratio>0.832966 && Rphi>0.595356 && Reta>0.936303 && f3<0.114127 && wEta2<0.012756 && wTotS1<4.696972 && trk_TRT_PID1>-0.075664 && trk_TRT_PID2>-0.054928 && EoP<2.614554 && trk_dEta1<0.001922)) return true;

     //signal efficiency: 636/1599=0.397749
     //bkg efficiency: 303/185142=0.00163658
     if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>0.800000 && fabs(etae)<1.370000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.027027 && Eratio>0.858457 && Rphi>0.361620 && Reta>0.918044 && f3<0.227751 && wEta2<0.013344 && wTotS1<3.594180 && trk_TRT_PID1>-0.081856 && trk_TRT_PID2>-0.046679 && EoP<1.798164 && trk_dEta1<0.002309)) return true;

     //signal efficiency: 464/1272=0.36478
     //bkg efficiency: 317/155927=0.002033
    // if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.008386 && Eratio>0.581808 && Rphi>0.584478 && Reta>0.921533 && f3<0.027292 && wEta2<0.015991 && wTotS1<4.975387 && trk_TRT_PID1>-0.009200 && trk_TRT_PID2>-0.233873 && EoP<2.976654 && trk_dEta1<0.002188)) return true;
    //signal efficiency: 401/1272=0.315252
    //bkg efficiency: 246/155927=0.00157766
    if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>1.520000 && fabs(etae)<2.010000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.031010 && Eratio>0.921214 && Rphi>0.610367 && Reta>0.925299 && f3<0.018186 && wEta2<0.013959 && wTotS1<3.378504 && trk_TRT_PID1>0.034933 && trk_TRT_PID2>-0.315406 && EoP<2.499431 && trk_dEta1<0.002653)) return true;

     //signal efficiency: 394/986=0.399594
     //bkg efficiency: 395/123543=0.00319727
    // if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.013270 && Eratio>0.933329 && Rphi>0.672153 && Reta>0.914604 && f3<0.103046 && wEta2<0.013576 && wTotS1<1.546700 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<2.026565 && trk_dEta1<0.001696)) return true;
    //signal efficiency: 320/986=0.324544
    //bkg efficiency: 253/123543=0.00204787
    if((wTotS1>-100. && pte>20.000000 && pte<30.000000 && fabs(etae)>2.010000 && fabs(etae)<2.470000 && fabs(vtx_dphi)<0.030000 && fabs(vtx_deta)<0.010000)&&(Rhad<0.026704 && Eratio>0.903645 && Rphi>0.678678 && Reta>0.919190 && f3<0.007630 && wEta2<0.013965 && wTotS1<1.467929 && trk_TRT_PID1>-100.000000 && trk_TRT_PID2>-100.000000 && EoP<1.889961 && trk_dEta1<0.001359)) return true;

    return(false);
}
