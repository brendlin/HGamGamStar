#include "HGamGamStar/HggStarCommon.h"
#include "xAODEgamma/EgammaxAODHelpers.h"
#include "xAODTruth/xAODTruthHelpers.h"


#include "HGamAnalysisFramework/HGamCommon.h"

TString HG::GetChannelName(ChannelEnum channel){
  if (channel == DIMUON              ) return "Dimuon";
  if (channel == RESOLVED_DIELECTRON ) return "ResolvedDielectron";
  if (channel == MERGED_DIELECTRON   ) return "MergedDielectron";
  if (channel == AMBIGUOUS_DIELECTRON) return "AmbiguousDielectron";
  if (channel == FAILEDTRKELECTRON   ) return "FailedTrackElectron";
  if (channel == OTHER               ) return "Other";
  if (channel == OUT_OF_ACCEPTANCE   ) return "OutOfAcceptance";
  return "";
}

xAOD::Photon*  HG::createPhotonFromElectron (const xAOD::Electron* el)
{


  int index1 = -999;
  int index2 = -999;

  if( HG::EleAcc::vtxTrkIndex1.isAvailable(*el) && HG::EleAcc::vtxTrkIndex2.isAvailable(*el)  ){
    index1 = HG::EleAcc::vtxTrkIndex1(*el);
    index2 = HG::EleAcc::vtxTrkIndex2(*el);
  }
  if(index1 < 0 || index2 < 0 ){
    return 0;
  }
  //std::cout << "Index 1/2  " << index1  << " " << index2 << std::endl;

  xAOD::Photon* photon = new xAOD::Photon();
  photon->makePrivateStore();

  if( el->ambiguousObject() ){
    //std::cout << "Copying photon" <<  std::endl;
    auto ambiPhoton = dynamic_cast<const xAOD::Photon*>( el->ambiguousObject() );
    photon->Photon_v1::operator=(*ambiPhoton);
    photon->setCaloClusterLinks(el->caloClusterLinks());
  } else {
    //std::cout << "Creating photon" <<  std::endl;
    photon = new xAOD::Photon();
    photon->Egamma_v1::operator=(*el);
    photon->setCaloClusterLinks(el->caloClusterLinks());
  }

  //std::cout << "Setting Topo" <<  std::endl;

  static SG::AuxElement::Decorator<int> nClu("numTopoClusters") ;
  auto clusterVec = xAOD::EgammaHelpers::getAssociatedTopoClusters( photon->caloCluster() );
  nClu(*el)= clusterVec.size();


  return photon;
}

void   HG::setPhotonConversionVertex( const xAOD::Electron* el,
                                       xAOD::Photon* photon, float vtxR,
                                     xAOD::VertexContainer* vertexContainer )
{

  int index1 = -999;
  int index2 = -999;

  if( HG::EleAcc::vtxTrkIndex1.isAvailable(*el) && HG::EleAcc::vtxTrkIndex2.isAvailable(*el)  ){
    index1 = HG::EleAcc::vtxTrkIndex1(*el);
    index2 = HG::EleAcc::vtxTrkIndex2(*el);
  }
  if(index1 < 0 || index2 < 0 ){
    return;
  }

  // If the merged ID acceptors do not exist, then stop (pre-p3877)
  if (!HG::EleAcc::vtxPhi.isAvailable(*el)) {
    return;
  }

  float vtxX = vtxR * cos(HG::EleAcc::vtxPhi(*el));
  float vtxY = vtxR * sin(HG::EleAcc::vtxPhi(*el));

  //std::cout << "Creating Vtx " <<  std::endl;

  xAOD::Vertex*  vertex = new xAOD::Vertex();
  vertexContainer->push_back(vertex);

  vertex->setZ(HG::EleAcc::vtxZ(*el));
  vertex->setX(vtxX);
  vertex->setY(vtxY);

  //std::cout << "Decorating Vtx " <<  std::endl;

  // decorate with pt1, pt2
  vertex->auxdata<float>("pt1") = el->trackParticle(index1)->pt();
  vertex->auxdata<float>("pt2") = el->trackParticle(index2)->pt();


  //std::cout << "Setting Trk element links Vtx " <<  std::endl;
  //links to tracks
  std::vector<ElementLink<xAOD::TrackParticleContainer>> links_tracks;
  links_tracks.push_back( el->trackParticleLinks()[index1] );
  links_tracks.push_back( el->trackParticleLinks()[index2] );

  //set vertex - track links
  vertex->setTrackParticleLinks(links_tracks);

  //std::cout << "Setting Vtx element links Photon " <<  std::endl;

  std::vector<ElementLink<xAOD::VertexContainer>> links_verticies;
  //std::cout << "-- Setting Vtx element links Photon " <<  std::endl;
  links_verticies.push_back(ElementLink<xAOD::VertexContainer>(vertex, *vertexContainer));
  //std::cout << "---- Setting Vtx element links Photon " <<  std::endl;
  photon->setVertexLinks(links_verticies);

  //std::cout << "Setting deltas " <<  std::endl;


  float dEta = HG::EleAcc::vtxdEta(*el);
  float dPhi = HG::EleAcc::vtxdPhi(*el);
  photon->setVertexCaloMatchValue( dEta, xAOD::EgammaParameters::convMatchDeltaEta1 );
  photon->setVertexCaloMatchValue( dEta, xAOD::EgammaParameters::convMatchDeltaEta2 );
  photon->setVertexCaloMatchValue( dPhi, xAOD::EgammaParameters::convMatchDeltaPhi1 );
  photon->setVertexCaloMatchValue( dPhi, xAOD::EgammaParameters::convMatchDeltaPhi2 );

  return;
}

HG::ChannelEnum HG::truthChannel(const xAOD::TruthParticleContainer& childleps,
                                 const xAOD::ElectronContainer& all_elecs)
{
  // Get Higgs Final State Decay Products

  if (!HG::isMC())
    return HG::CHANNELUNKNOWN;

  if (childleps.size() != 2)
    return HG::OTHER;

  // Check for out-of-acceptance
  TLorentzVector allLeps;

  for(const auto& lepton: childleps){
    /*
    if (fabs(lepton->pdgId()) == 11) {
      if (lepton->pt()/1000. < 0.3) return HG::OUT_OF_ACCEPTANCE;
      if (fabs(lepton->eta()) > 2.5) return HG::OUT_OF_ACCEPTANCE;
    }
    else if (fabs(lepton->pdgId()) == 13) {
      if (lepton->pt()/1000. < 3.0) return HG::OUT_OF_ACCEPTANCE;
      if (fabs(lepton->eta()) > 2.7) return HG::OUT_OF_ACCEPTANCE;
    }
    */
    allLeps += lepton->p4();
  }

  if( allLeps.M() > 50e3 )
    return HG::OUT_OF_ACCEPTANCE;

  // Check if there are electrons in the decay
  bool isElectron = true;
  bool isMuon =  false;
  for(const auto& lepton: childleps){
    if( fabs(lepton->pdgId()) != 11 )
      isElectron =  false;
    if( fabs(lepton->pdgId()) == 13 )
      isMuon  = true;
  }

  if(isMuon)
    return HG::DIMUON;

  if(!isElectron)
    return HG::OTHER;

  // Make a track-truth map (one-to-one map)
  HG::TruthTrackMap trkTruthMap;

  // trkEleMap (multimap) for channel classification
  HG::TrackElectronMap trkEleMap;

  for (auto electron : all_elecs) {
    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {
      const xAOD::TrackParticle* trkParticle = electron->trackParticle(i);

      int nSi = xAOD::EgammaHelpers::numberOfSiHits( trkParticle );
      if( nSi < 3 ) continue;

      const xAOD::TruthParticle* truthPart = xAOD::TruthHelpers::getTruthParticle(*trkParticle);
      if (!truthPart) continue;

      // Find tracks that match to our two leptons
      for(const auto& lepton: childleps){

        // Need to match the original pointers
        auto lepton_p = (xAOD::TruthParticle*)MapHelpers::getTheOriginalPointer(*lepton);

        if (truthPart == lepton_p) {

          // Add the lepton (provided its truth match probability beats out the others)
          MapHelpers::AddTruthTrackMapEntry(trkParticle,trkTruthMap);

          // Add track and electron to map (it is okay if we have extras here)
          MapHelpers::AddTrackElectronMapEntry(trkParticle,electron,trkEleMap);
        }
      }
    }
  }

  if (trkTruthMap.size() != 2)
    return HG::FAILEDTRKELECTRON;

  const xAOD::TrackParticle* trk0 = trkTruthMap.begin()->second; // get the first track
  const xAOD::TrackParticle* trk1 = trkTruthMap.rbegin()->second; // get the second track

  // Electron disambiguation is continued below in (reco-only) function
  return ClassifyElectronChannelsByBestMatch(trk0,trk1,trkEleMap);
}

HG::ChannelEnum HG::ClassifyElectronChannelsByBestMatch(const xAOD::TrackParticle* trk0,
                                                        const xAOD::TrackParticle* trk1,
                                                        const HG::TrackElectronMap& trkEleMap,
                                                        xAOD::ElectronContainer* inEleCont,
                                                        xAOD::ElectronContainer* outEleCont)
{

  // Find the electrons associated to the tracks
  auto Trk0_Electrons = HG::MapHelpers::getElectronsMatchingTrack( trk0, trkEleMap );
  auto Trk1_Electrons = HG::MapHelpers::getElectronsMatchingTrack( trk1, trkEleMap );

  // Count the number of electrons a track matches to
  int Trk0_nElectron = Trk0_Electrons.size();
  int Trk1_nElectron = Trk1_Electrons.size();

  // If each track only matches to 1 electron each then it is quite simple
  if( Trk0_nElectron == 1 &&  Trk1_nElectron == 1){
    // Match the same electron --  Merged
    if( Trk0_Electrons.front() == Trk1_Electrons.front() )
    {
      if (inEleCont)
        outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,Trk0_Electrons.front()));
      return HG::MERGED_DIELECTRON;
    }
    else
    {
      // Match different electrons  -- Resolved
      if (inEleCont) {
        outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,Trk0_Electrons.front()));
        outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,Trk1_Electrons.front()));
      }
      return HG::RESOLVED_DIELECTRON;
    }
  }

  // Tracks match more than one electron each - lets see if they are the best match for any electron
  // Determine the match ranking for the track to each electron
  std::vector<int> Trk0_TrackNo = HG::MapHelpers::getMatchingTrackIndices(Trk0_Electrons,trk0);

  // Check if the track is ever the primary track
  int Trk0_PrimaryE(-1);
  for( unsigned int i(0); i < Trk0_TrackNo.size(); ++i){
    if( Trk0_TrackNo[i] == 0 ){
      // If the track is the primary track for multiple electrons
      // choose the one with higher pT
      if( Trk0_PrimaryE > -1 ){
        if( Trk0_Electrons[i]->pt() > Trk0_Electrons[Trk0_PrimaryE]->pt() )
          Trk0_PrimaryE = i;
      }else{
        Trk0_PrimaryE = i;
      }
    }
  }

  // Same again for the other electron
  std::vector<int> Trk1_TrackNo = HG::MapHelpers::getMatchingTrackIndices(Trk1_Electrons,trk1);

  int Trk1_PrimaryE(-1);
  for( unsigned int i(0); i < Trk1_TrackNo.size(); ++i){
    if( Trk1_TrackNo[i] == 0 ){
      // If the track is the primary track for multiple electrons
      // choose the one with higher pT
      if( Trk1_PrimaryE > -1 ){
        if( Trk1_Electrons[i]->pt() > Trk1_Electrons[Trk1_PrimaryE]->pt() )
          Trk1_PrimaryE = i;
      }else{
        Trk1_PrimaryE = i;
      }
    }
  }

  // If both tracks are the primary track for an electron the it is resolved
  if( Trk0_PrimaryE > -1 && Trk1_PrimaryE > -1 ){
    // Get to the correct pair iterator
    // This is a pair of TrackParticle , Electron
    auto el0 = Trk0_Electrons[Trk0_PrimaryE];
    auto el1 = Trk1_Electrons[Trk1_PrimaryE];

    // Compare if the electrons are the same
    if( el0 == el1 )
    {
      HG::fatal("Electron classification truth error!");
      return HG::AMBIGUOUS_DIELECTRON;
    }
    if (inEleCont) {
      outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,el0));
      outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,el1));
    }
    return HG::RESOLVED_DIELECTRON;
  }

  // If either are primary
  if( Trk0_PrimaryE > -1 || Trk1_PrimaryE > -1 ){
    const xAOD::Electron*  el = nullptr;
    const xAOD::Electron*  elOther = nullptr;
    const xAOD::TrackParticle* otherTrack = nullptr;
    int nEleOther = 0;
    if( Trk0_PrimaryE >= 0 ){
      el = Trk0_Electrons[Trk0_PrimaryE];
      otherTrack = trk1;
      nEleOther = Trk1_nElectron;
      if (nEleOther == 1) elOther = Trk1_Electrons[0];
    }else{
      el = Trk1_Electrons[Trk1_PrimaryE];
      otherTrack = trk0;
      nEleOther = Trk0_nElectron;
      if (nEleOther == 1) elOther = Trk0_Electrons[0];
    }

    // Search for the other track in the electron
    for( unsigned int trk_i(0); trk_i < el->nTrackParticles(); ++trk_i){
      if( el->trackParticle(trk_i) == otherTrack )
      {
        if (inEleCont)
          outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,el));
        return HG::MERGED_DIELECTRON;
      }
    }
    //If the other track is only matched to one electron and its not the primary track
    //We assume the reco has made a mistake so we will call it resolved
    if(nEleOther == 1)
    {
      if (inEleCont) {
        outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,el));
        outEleCont->push_back(HG::MapHelpers::FindElectron(inEleCont,elOther));
      }
      return HG::RESOLVED_DIELECTRON;
    }
  }

  //Tracks are not primary for any electron candidate, tracks match to mutiple candiates
  // --  generally the reco has made a mess.
  // Might want to keep looking for candidates

  return HG::AMBIGUOUS_DIELECTRON;
}
