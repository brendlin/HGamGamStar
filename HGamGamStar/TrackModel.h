/*  ============================================================
Filename: TrackModel.h
Purpose:  Class for analysis of track model
============================================================ */

#ifndef HGamGamStar_TrackModel
#define HGamGamStar_TrackModel

#include <string>
#include <vector>

class TH3D;
class TGraph;


double sign(const double x);


class TrackModel
{

public:
  
  //!< Public structures
  struct HelixPosition
  { 
  
    HelixPosition()
    : x(0),
      y(0),
      z(0)
    {
      init();
    }
    
    HelixPosition(double x1,double y1, double z1)
    : x(x1),
      y(y1),
      z(z1)
    {
      init();
    }
  
    HelixPosition operator+(const HelixPosition& m) const
    {
      HelixPosition b(x + m.x , y + m.y, z + m.z);
      b.init(); 
      return b;
    }

    HelixPosition& operator+=(const HelixPosition& m)
    {
      x += m.x; 
      y += m.y;
      z += m.z;
      init(); 
      return *this;
    }
  
    HelixPosition& operator*=(const double& d) 
    {
      x *= d; 
      y *= d;
      z *= d;
      init();
      return *this;
    }
   
    HelixPosition operator*(const double& d) const
    {
      HelixPosition b (x * d, y * d, z * d);
      b.init(); 
      return b;
    }
  
  
    double x;
    double y;
    double z;

    double r;
    double phi;
    double theta;
    double eta;
     
    void init()
    {
      calc_r();
      calc_phi();
      calc_theta();
      calc_eta();
    } 
     
  private:
    

    
    void calc_r(){
      r =  sqrt( x*x + y*y );
    }
    
    void calc_phi(){
      phi =  atan2(y,x);
    }
    
    void calc_theta(){
      theta = atan2(r,z);
    }   

    void calc_eta(){
      eta = -log(tan(theta*0.5));
    }   

  };
  
  //!< Public structures
  struct HelixMomentum
  { 
  
    HelixMomentum()
    : x(0),
      y(0),
      z(0)
    {
      init();
    }
    
    HelixMomentum(double x1,double y1, double z1)
    : x(x1),
      y(y1),
      z(z1)
    {
      init();
    }

    HelixMomentum(double x1,double y1, double z1, double invP)
    : x(x1*invP),
      y(y1*invP),
      z(z1*invP)
    {
      init();
    }
  
  
  
    HelixMomentum operator+(const HelixMomentum& m) const
    {
      HelixMomentum b(x + m.x , y + m.y, z + m.z);
      b.init(); 
      return b;
    }

    HelixMomentum& operator+=(const HelixMomentum& m)
    {
      x += m.x; 
      y += m.y;
      z += m.z;
      init(); 
      return *this;
    }
  
    HelixMomentum& operator*=(const double& d)
    {
      x *= d; 
      y *= d;
      z *= d;
      init();
      return *this;
    }
   
    HelixMomentum operator*(const double& d) const
    {
      HelixMomentum b (x * d, y * d, z * d);
      b.init(); 
      return b;
    }
  
  
  
    double x;
    double y;
    double z;

    double pT;
    double phi;
    double theta;
     
    double p(){
      return sqrt( x*x + y*y + z*z);
    }
     
  private:
    
    void init()
    {
      calc_pT();
      calc_phi();
      calc_theta();
    }
    
    void calc_pT(){
      pT =  sqrt( x*x + y*y );
    }
    
    void calc_phi(){
      phi =  atan2(y,x);
    }
    
    void calc_theta(){
      theta = atan2(pT,z);
    }   
  };
  
  
  struct HelixState
  {

    HelixState(){ 
      position = HelixPosition();
      momentum = HelixMomentum(); 
      qOverPT  = 0.001;
    }   

    HelixState( HelixPosition x, HelixMomentum p, double qonpT){ 
      position = x;
      momentum = p; 
      qOverPT  = qonpT;
    }   


    HelixState operator+(const HelixState& m) const
    {
      HelixMomentum hm = momentum + m.momentum;
      HelixState b( position + m.position , hm, sign(qOverPT)/hm.pT );
      return b;
    }

    HelixState& operator+=(const HelixState& m)
    {
      position += m.position;
      momentum += m.momentum;
      qOverPT = sign(qOverPT)/momentum.pT;
      return *this;
    }
  
    HelixState& operator*=(const double& d)
    {
      position *= d;
      momentum *= d;
      qOverPT /= d;
      return *this;
    }
   
    HelixState operator*(const double& d) const
    {
      HelixState b( position *d , momentum * d, qOverPT / d );
      return b;
    }


    HelixPosition position;
    HelixMomentum momentum;
    double qOverPT;  
  };
  
  
  struct HelixParameters
  {
  
    //!< Default constructor
    HelixParameters()
    :
    x0(0.),
    y0(0.),
    z0(0.),
    phi0(0.),
    cotTheta(0.),
    qOverPT(0.)
    {}
  
   //!< Constructor with parameters
   HelixParameters( double theX0, double theY0, double theZ0, double thePhi0, double theCotTheta, double theQoverPT )
    :
    x0(theX0),
    y0(theY0),
    z0(theZ0),
    phi0(thePhi0),
    cotTheta(theCotTheta),
    qOverPT(theQoverPT)
    {}
  
    double x0;
    double y0;
    double z0;
    double phi0;
    double cotTheta;
    double qOverPT;
    
  };

  struct ParameterisedHelix
  {
    
    TH3D* helixHistogram;
    
  };


  TrackModel( const HelixParameters&);
  
  //!< Default constructor
  TrackModel();
  
  //!< Destructor
  ~TrackModel();
  
  //!< Calculate the position at some R or Z boundary 
  TrackModel::HelixState  stateAtBoundaryRZ( const double maxR, const double  maxZ, const bool printTrajectory) const;

  
  //!< run method
  void run( double x0, double y0, double z0, double phi0, double cotTheta, double qOverPT ) const;
  
  //!< method to determine helix position from given pathlength and parameters
  const HelixPosition determineHelix( double pathLength, const HelixParameters& parameters ) const;
  
  //!< method to determine helix from a given phi and parameters
  const HelixPosition determineHelixFromPhi( double phi, const HelixParameters& parameters ) const;
  
  //<! method to get the trajectory histogram
  TH3D* trajectoryHistogram( const HelixParameters&, std::string instanceName ) const;
  
  //!< method to generate an r / phi graph
  TGraph* phiVsR( const TrackModel::HelixParameters& parameters ) const;
  
  //!< Method to generate a cotTheta / phi graph
  TGraph* phiVsCotTheta( const TrackModel::HelixParameters& parameters ) const;
  
  //!< Method to generate a cot Theta / rho graph
  TGraph* cotThetaVsRho( const TrackModel::HelixParameters& parameters ) const;
  
  //!< method to generate an r / eta graph
  TGraph* etaVsR( const TrackModel::HelixParameters& parameters ) const;

  //!< method to generate an r / eta graph
  TGraph* etaVsZ( const TrackModel::HelixParameters& parameters ) const;

  //!< method to generate an r / z graph
  TGraph* RVsZ( const TrackModel::HelixParameters& parameters ) const;



  //<! method to get the trajectory histogram
  TH3D* trajectoryHistogram( std::string instanceName ) const;
  
  //!< method to generate an r / phi graph
  TGraph* phiVsR() const;
  
  //!< Method to generate a cotTheta / phi graph
  TGraph* phiVsCotTheta() const;
  
  //!< Method to generate a cot Theta / rho graph
  TGraph* cotThetaVsRho() const;
  
  //!< method to generate an r / eta graph
  TGraph* etaVsR() const;

  //!< method to generate an r / eta graph
  TGraph* etaVsZ() const;

  //!< method to generate an r / z graph
  TGraph* RVsZ() const;
  
  
  std::vector<HelixPosition> positionAtR(std::vector<double>&  radii );

  
  int RungeKuttaNystroem(const double step, const TrackModel::HelixState& vIn, TrackModel::HelixState& vOut) const;

  TrackModel::HelixState determineInitialHelixState(const TrackModel::HelixParameters& parameters ) const;


private:
  
  void BuildTrajectory() const;
  void BuildTrajectoryOLD() const;
  

  //!< Private members
  mutable std::vector<double> m_cotThetaVector;
  mutable std::vector<double> m_phiVector;
  mutable std::vector<double> m_rVector;
  mutable std::vector<double> m_xVector;
  mutable std::vector<double> m_yVector;
  mutable std::vector<double> m_zVector;
  mutable std::vector<double> m_sVector;
  mutable std::vector<double> m_etaVector;
    
  mutable const TrackModel::HelixParameters* m_cachedparameters;  
  
  mutable double m_theta;
  mutable double m_pT;
  mutable double m_radiusOfCurvature;
  mutable double m_phi;
    
};


#endif
