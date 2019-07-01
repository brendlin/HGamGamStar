#include "HGamGamStar/HggStarCommon.h"
#include "xAODEgamma/EgammaxAODHelpers.h"



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
