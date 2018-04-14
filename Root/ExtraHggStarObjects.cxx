#include "HGamGamStar/ExtraHggStarObjects.h"

//____________________________________________________________________________
namespace HG {
  ExtraHggStarObjects *ExtraHggStarObjects::m_ptr = nullptr;

  ExtraHggStarObjects *ExtraHggStarObjects::getInstance()
  {
    if (m_ptr == nullptr)
    { m_ptr = new ExtraHggStarObjects(); }

    return m_ptr;
  }

}

//____________________________________________________________________________
HG::ExtraHggStarObjects::ExtraHggStarObjects()
{

}

//____________________________________________________________________________
HG::ExtraHggStarObjects::~ExtraHggStarObjects()
{

}
