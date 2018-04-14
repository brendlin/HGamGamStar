#ifndef HGamGamStar_ExtraHggStarObjects_H
#define HGamGamStar_ExtraHggStarObjects_H

//
// This class is simply meant to save some extra, easily accessible variables
// so they do not have to be calculated again (such as the selected tracks, or the
// true uh stuff.
//

namespace HG {

  class ExtraHggStarObjects {
  private:
    static ExtraHggStarObjects *m_ptr;

  public:
    /// Get instance of singleton class
    static ExtraHggStarObjects *getInstance();

  private:
    /// Default constructor - note that it is private
    ExtraHggStarObjects();

    /// Default detructor - note that is is private
    ~ExtraHggStarObjects();

  }; // class ExtraHggStarObjects

}

#endif // HGamGamStar_ExtraHggStarObjects_H
