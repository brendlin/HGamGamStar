#include <fstream>
#include <iostream>
#include <math.h>

#include "HGamGamStar/TrackModel.h"

#include "TGraph.h"
#include "TH3D.h"
#include "TF1.h"

#include <assert.h>

double sign(const double x) 
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}


static const int    s_numberOfSteps = 30000;
static const double s_stepSize      = 0.001;
static const double s_pi            = 3.141592654;
static const double s_ec   = 2.9979251E-1;  // This constant is for units M,GEV/C and T
static const double s_field = 2.; //Field strength in Tesla

TrackModel::TrackModel( const HelixParameters& hp )
{
  m_theta=0;
  m_pT=0;
  m_radiusOfCurvature=0;
  m_cachedparameters= &hp;
  // BuildTrajectory(); // this doesn't do anything anyways
}


TrackModel::TrackModel()
{
  m_theta=0;
  m_pT=0;
  m_radiusOfCurvature=0;
  m_cachedparameters=0;
}

TrackModel::~TrackModel()
{}

void TrackModel::run( double x0, double y0, double z0, double phi0, double cotTheta0, double qOverPT ) const
{

  
  TrackModel::HelixParameters parameters( x0, y0, z0, phi0, cotTheta0, qOverPT );

  TH3D* trajectory = trajectoryHistogram( parameters, "test" );
  
  trajectory->Draw();
  
}

const TrackModel::HelixPosition TrackModel::determineHelix( double pathLength, const TrackModel::HelixParameters& parameters ) const
{
  if (!m_pT || m_cachedparameters != &parameters ) {
    m_theta = atan( 1. / parameters.cotTheta );
    m_pT = 1. / parameters.qOverPT;
    m_radiusOfCurvature =  m_pT / ( s_field * s_ec ) ;
    m_cachedparameters = &parameters;
  }
  
  double phi = parameters.phi0 + pathLength * sin( m_theta ) / m_radiusOfCurvature;
  
  return determineHelixFromPhi( phi, parameters );  
  
}

const TrackModel::HelixPosition TrackModel::determineHelixFromPhi( double phi, const TrackModel::HelixParameters& parameters ) const
{
  
  TrackModel::HelixPosition position;
    
  position.x = parameters.x0 + m_radiusOfCurvature * ( sin( phi ) - sin( parameters.phi0 ) );
  position.y = parameters.y0 - m_radiusOfCurvature * ( cos( phi ) - cos( parameters.phi0 ) );
  position.z = parameters.z0 + m_radiusOfCurvature * parameters.cotTheta * ( phi - parameters.phi0 );
  
  //position.x = parameters.z0 + radiusOfCurvature * parameters.cotTheta * ( phi - parameters.phi0 );
  //position.y = -1. * ( parameters.x0 + radiusOfCurvature * ( sin( phi ) - sin( parameters.phi0 ) ) );
  //position.z = parameters.y0 - radiusOfCurvature * ( cos( phi ) - cos( parameters.phi0 ) );
  
  return position;
  
}

TrackModel::HelixState TrackModel::determineInitialHelixState(const TrackModel::HelixParameters& parameters ) const
{
  

  if (!m_pT || m_cachedparameters != &parameters ) {
    m_theta = atan( 1. / parameters.cotTheta );
    m_pT = fabs(1. / parameters.qOverPT);
    m_radiusOfCurvature =  m_pT / ( 2. * s_ec ) ;
    m_cachedparameters = &parameters;
  }
  
  
  double sinPhi0 = sin( parameters.phi0 );
  double cosPhi0 = cos( parameters.phi0 );
  
  double x = parameters.x0;// + m_radiusOfCurvature * ( 0  - sinPhi0 );
  double y = parameters.y0;// - m_radiusOfCurvature * ( 1. - cosPhi0 );
  double z = parameters.z0;// + m_radiusOfCurvature * parameters.cotTheta * ( - parameters.phi0 );

  TrackModel::HelixPosition position(x,y,z);
  
  double px = m_pT * cosPhi0;
  double py = m_pT * sinPhi0;
  double pz = m_pT * parameters.cotTheta;

  TrackModel::HelixMomentum momentum(px,py,pz);
  
  return TrackModel::HelixState( position, momentum, parameters.qOverPT);
  
}



void TrackModel::BuildTrajectoryOLD() const{
  
  if (!m_cachedparameters) {
    std::cout << "ERROR no parameters" << std::endl;
    return;
  }
  
  m_pT = 0;
  m_cotThetaVector.clear();
  m_phiVector.clear();
  m_rVector.clear();
  m_sVector.clear();
  m_xVector.clear();
  m_yVector.clear();
  m_zVector.clear();
  m_etaVector.clear();

  double s(0.);
  
  int index(0);
  
  double currentRadius(0.);
  
  while ( currentRadius <= 1.2 && index < s_numberOfSteps ){
    
    TrackModel::HelixPosition position = determineHelix( s, *m_cachedparameters );
    //if(index == 1)
    //  std::cout << position.z << " " <<  position.x << " " <<  position.y << std::endl;
      
    if ( position.x || position.y ){
      
      double phi      = atan2( position.y , position.x );
      double r        = sqrt( position.x * position.x + position.y * position.y );
      double cotTheta = position.z / r;
      
      m_cotThetaVector.push_back( cotTheta );
      m_phiVector.push_back( phi );
      m_rVector.push_back( r );
      m_xVector.push_back( position.x );
      m_yVector.push_back( position.y );
      m_zVector.push_back( position.z );
      m_sVector.push_back( s );
      
      double theta  = atan2(r, position.z);
      m_etaVector.push_back( -log(tan(theta*0.5))); //atanh( position.z /  sqrt( position.x * position.x + position.y * position.y  + position.z * position.z  ) ) );
            
    }
    
    currentRadius = sqrt( position.x * position.x + position.y * position.y );
        
    s += s_stepSize;
    
    ++index;
    
  }


}


void TrackModel::BuildTrajectory() const{
  
  if (!m_cachedparameters) {
    std::cout << "ERROR no parameters" << std::endl;
    return;
  }
  
  m_pT = 0;
  m_cotThetaVector.clear();
  m_phiVector.clear();
  m_rVector.clear();
  m_sVector.clear();
  m_xVector.clear();
  m_yVector.clear();
  m_zVector.clear();
  m_etaVector.clear();


  double s(0.);
  
  int index(0);
  
  double currentRadius(0.);
  double currentZ(0.);
  
  TrackModel::HelixState initState = determineInitialHelixState( *m_cachedparameters );
  
  
  
  while ( currentRadius <= 1.2 && fabs(currentZ) <= 3 && index < s_numberOfSteps ){
    
    TrackModel::HelixState nextState;
      
    int RkStatus = RungeKuttaNystroem( s_stepSize, initState, nextState  );
    if(RkStatus !=0 ) break;
    
    double cotTheta = nextState.position.z / nextState.position.r;
    
    m_cotThetaVector.push_back( cotTheta );
    m_phiVector.push_back( nextState.position.phi );
    m_rVector.push_back( nextState.position.r );
    m_xVector.push_back( nextState.position.x );
    m_yVector.push_back( nextState.position.y );
    m_zVector.push_back( nextState.position.z );
    m_sVector.push_back( s );
    m_etaVector.push_back( -log(tan(nextState.position.theta*0.5))); //atanh( position.z /  sqrt( position.x * position.x + position.y * position.y  + position.z * position.z  ) ) );
            
    
    currentRadius = nextState.position.r;
    currentZ = nextState.position.z;    
    s += s_stepSize;
    
    ++index;
    
    initState = nextState;
    
  }
  
  return;

}




int TrackModel::RungeKuttaNystroem(const double step, const HelixState& vIn, HelixState& vOut ) const
{
  //
  // Runge-Kutta method for tracking a particle through a magnetic field.
  // Uses NYSTROEM algorithm (See HANDBOOK NAT. BUR. OF STANDARDS, PROCEDURE 25.5.20)
  // Input parameters:  charge  - particle charge
  //      step  - step size
  //      vIn - initial co-ords, direction , momentum
  //
  // Output parameters: vOut  - output co-ords, direction , momentum
  //
  // Return values: 0 - OK
  //      8 - underflow of output direction 
  //
  // Authors (original Fortran code): R.Brun, M.Hansroul
  //

  // const double ec   = 2.9979251E-4;  // This constant is for units CM,GEV/C and KGAUSS
  //const double ec   = 2.9979251E-1; // This constant is for units M,GEV/C and T
  const double pinv = s_ec * vIn.qOverPT;

  const double h    = step;
  const double h2  = 0.5 * h;
  const double h4  = 0.25 * h;
  const double ph   = pinv * h;
  const double ph2 = 0.5 * ph;

  double xyz[3], xyzt[3];
  double a, d, c, at, dt, ct;
  double secxs[4], secys[4], seczs[4];
  double field[3] = {0,0,2.};

  /* First intermediate point */
  xyz[0] = vIn.position.x;
  xyz[1] = vIn.position.y;
  xyz[2] = vIn.position.z;
  
    
  a = vIn.momentum.x;
  d = vIn.momentum.y;
  c = vIn.momentum.z;
  
  double adcTot =sqrt(a*a + d*d + c*c);
  double invadcTot = 1./adcTot;
  a*=invadcTot; 
  d*=invadcTot; 
  c*=invadcTot; 
  
  secxs[0] = (d * field[2] - c * field[1]) * ph2;
  secys[0] = (c * field[0] - a * field[2]) * ph2;
  seczs[0] = (a * field[1] - d * field[0]) * ph2;
  xyzt[0] = xyz[0] + h2 * a + h4 * secxs[0];
  xyzt[1] = xyz[1] + h2 * d + h4 * secys[0];
  xyzt[2] = xyz[2] + h2 * c + h4 * seczs[0];

  /* Second intermediate point */
  at = a + secxs[0];
  dt = d + secys[0];
  ct = c + seczs[0];
  secxs[1] = (dt * field[2] - ct * field[1]) * ph2;
  secys[1] = (ct * field[0] - at * field[2]) * ph2;
  seczs[1] = (at * field[1] - dt * field[0]) * ph2;
  at = a + secxs[1];
  dt = d + secys[1];
  ct = c + seczs[1];
  secxs[2] = (dt * field[2] - ct * field[1]) * ph2;
  secys[2] = (ct * field[0] - at * field[2]) * ph2;
  seczs[2] = (at * field[1] - dt * field[0]) * ph2;
  xyzt[0] = xyz[0] + h * (a + secxs[2]);
  xyzt[1] = xyz[1] + h * (d + secys[2]);
  xyzt[2] = xyz[2] + h * (c + seczs[2]);
  at = a + 2. * secxs[2];
  dt = d + 2. * secys[2];
  ct = c + 2. * seczs[2];

  /* Final point */
  double oneThird = 1./3.;
  xyz[0] += (a + (secxs[0] + secxs[1] + secxs[2]) * oneThird) * h;
  xyz[1] += (d + (secys[0] + secys[1] + secys[2]) * oneThird) * h;
  xyz[2] += (c + (seczs[0] + seczs[1] + seczs[2]) * oneThird) * h;
  secxs[3] = (dt * field[2] - ct * field[1]) * ph2;
  secys[3] = (ct * field[0] - at * field[2]) * ph2;
  seczs[3] = (at * field[1] - dt * field[0]) * ph2;
  a += (secxs[0] + secxs[3] + 2. * (secxs[1] + secxs[2])) * oneThird;
  d += (secys[0] + secys[3] + 2. * (secys[1] + secys[2])) * oneThird;
  c += (seczs[0] + seczs[3] + 2. * (seczs[1] + seczs[2])) * oneThird;

  /* Final processing and storing of output */
  double adcTot2 =sqrt(a*a + d*d + c*c);
  if (adcTot2 < 1e-37)
      return 8;
  
  a*=adcTot;
  d*=adcTot;
  c*=adcTot;
  
  //std::cout << "a,d,c " << a << "  "<<  d << "  "<< c << " " <<  xyz[0]<< " "<<  xyz[1] << " "<< xyz[2]<< std::endl; 
  
  vOut.position = HelixPosition( xyz[0], xyz[1], xyz[2] );
  vOut.momentum = HelixMomentum( a, d, c);
  vOut.qOverPT  = vIn.qOverPT;

  return 0;

}


std::vector<TrackModel::HelixPosition> TrackModel::positionAtR(std::vector<double>&  radii )
{
  //Radii must be in assending order;
  
  std::vector<TrackModel::HelixPosition> positions;
  
  unsigned int nr = radii.size();
  int r(0);
  
  if(nr < 2)
    return positions;
  
  // Ensure target radius is larger than the current radius 
  double dr = m_rVector[1] - m_rVector[0];
  if (dr > 0) {
    while(radii[r] < m_rVector[0]){
      ++r;
    } 
  } 

  double previousDr = 0;
  for( unsigned int i(1); i < m_rVector.size()-1; ++i){

    previousDr = dr;
    dr = m_rVector[i+1] - m_rVector[i];
    //std::cout << radii[r] <<"  "<<   m_rVector[i]  << " " << dr <<  std::endl;
    if (dr * previousDr <= 0){
      // turning point look for layer
      //std::cout << "Sign Flip" <<  std::endl;
      if (dr > 0){
        ++r;
      } else if ( dr < 0) {
        --r;
      }
      
      if(r < 0 ||  r >= (int)nr){
        std::cout << "ERROR  r out of range" << std::endl;
        r = 0;
      }
      
    }
    // CHeck if we have crossed a boundry;
    if(dr > 0 &&  m_rVector[i]  > radii[r]){
      
      double diff    = 1./(m_rVector[i] - m_rVector[i-1]);
      double weight1 = (m_rVector[i]-radii[r]) * diff;
      double weight2 = (radii[r]-m_rVector[i-1]) * diff;
      
      //std::cout << "POS " << i << " " << radii[r] <<"  "<<   m_rVector[i]  << " " << weight1*m_xVector[i]+weight2*m_xVector[i-1] <<  std::endl;
      
      positions.push_back(TrackModel::HelixPosition( weight2*m_xVector[i]+weight1*m_xVector[i-1], 
                                         weight2*m_yVector[i]+weight1*m_yVector[i-1],
                                         weight2*m_zVector[i]+weight1*m_zVector[i-1] ) );
      // Increase the next search radius                                   
      if(r < (int)nr-1 )++r;
      
                                         
    } else if (dr < 0 &&  m_rVector[i]  < radii[r]) {

      double diff    = 1./(m_rVector[i] - m_rVector[i-1]);
      double weight1 = (m_rVector[i]-radii[r]) * diff;
      double weight2 = (radii[r]-m_rVector[i-1]) * diff;

      //std::cout << "NEG " << i << " " << radii[r] <<"  "<<   m_rVector[i]  <<" " << m_rVector[i-1] << " " << weight2*m_rVector[i]+weight1*m_rVector[i-1] <<  std::endl;
      
      positions.push_back(TrackModel::HelixPosition( weight2*m_xVector[i]+weight1*m_xVector[i-1], 
                                         weight2*m_yVector[i]+weight1*m_yVector[i-1],
                                         weight2*m_zVector[i]+weight1*m_zVector[i-1] ) );
      // Decrease the next search radius as dr is neg                                                                          
      if (r > 0) --r;                                   
    }
  }
  
  return positions;
}




TrackModel::HelixState  TrackModel::stateAtBoundaryRZ( const double maxR, const double  maxZ, const bool printTrajectory) const
{  
  
  TrackModel::HelixState nextState;

  if (!m_cachedparameters) {
    std::cout << "ERROR no parameters" << std::endl;
    return nextState;
  }
  

  double s(0.);
  
  int index(0);
  
  double currentRadius(0.);
  double currentZ(0.);
  
  TrackModel::HelixState startState = determineInitialHelixState( *m_cachedparameters );
  TrackModel::HelixState initState = startState;
  
  // safety check: both coordinates outside the boundaries should not happen!
  // just one coordinate outside the boundaries may happen because
  // the trackParameters(LastMeasurement) is an extrapolation and not
  // really a measurement
  if (initState.position.r > maxR && fabs(initState.position.z) > maxZ){
    std::cout << "initial position outside of maxR and maxZ" << std::endl;
    std::cout << "  initial position is" << std::endl;
    std::cout << Form("  - (x, y, z)       = (%2.5f, %2.5f, %2.5f)", initState.position.x, initState.position.y, initState.position.z) << std::endl;
    std::cout << Form("  - (r, theta, phi) = (%2.5f, %2.5f, %2.5f)", initState.position.r, initState.position.theta, initState.position.phi) << std::endl;
    std::cout << Form("  boundaries  are (r, z) = (%2.5f, %2.5f)", maxR, maxZ) << std::endl << std::flush;
    assert(false);
  }

  std::ofstream ofs;
  if (printTrajectory) ofs.open("/nfs/dust/atlas/user/rauch/Hyy/software/ATLAS-Hystary-ll-trunk/trajectory.dat", std::ofstream::out);

  while ( index < s_numberOfSteps ){

    int RkStatus = RungeKuttaNystroem( s_stepSize, initState, nextState  );
    if (printTrajectory) ofs << nextState.position.x << "   " << nextState.position.y << "   " << nextState.position.z << std::endl;
    if(RkStatus !=0 ) break;
    
    double cotTheta = nextState.position.z / nextState.position.r;

    currentRadius = nextState.position.r;
    currentZ = nextState.position.z;

    ++index;

    if( currentRadius > maxR || fabs(currentZ) > maxZ ) break;

    initState = nextState;

  }

  if (printTrajectory) ofs.close();

  // safety check
  if (index == s_numberOfSteps){
    std::cout << "failed to capture boundary - RKN terminated by step number (index = " << index << ")" << std::endl;
    std::cout << "  - start" << std::endl;
    std::cout << "    - position was" << std::endl;
    std::cout << Form("      - (x, y, z)        = (%2.5f, %2.5f, %2.5f)", startState.position.x, startState.position.y, startState.position.z) << std::endl;
    std::cout << Form("      - (r, theta, phi)  = (%2.5f, %2.5f, %2.5f)", startState.position.r, startState.position.theta, startState.position.phi) << std::endl;
    std::cout << "    - momentum was" << std::endl;
    std::cout << Form("      - (x, y, z)        = (%2.5f, %2.5f, %2.5f)", startState.momentum.x, startState.momentum.y, startState.momentum.z) << std::endl;
    std::cout << Form("      - (pT, theta, phi) = (%2.5f, %2.5f, %2.5f)", startState.momentum.pT, startState.momentum.theta, startState.momentum.phi) << std::endl;
    std::cout << "  - current" << std::endl;
    std::cout << "    - position is" << std::endl;
    std::cout << Form("      - (x, y, z)        = (%2.5f, %2.5f, %2.5f)", nextState.position.x, nextState.position.y, nextState.position.z) << std::endl;
    std::cout << Form("      - (r, theta, phi)  = (%2.5f, %2.5f, %2.5f)", nextState.position.r, nextState.position.theta, nextState.position.phi) << std::endl;
    std::cout << "    - momentum is" << std::endl;
    std::cout << Form("      - (x, y, z)        = (%2.5f, %2.5f, %2.5f)", nextState.momentum.x, nextState.momentum.y, nextState.momentum.z) << std::endl;
    std::cout << Form("      - (pT, theta, phi) = (%2.5f, %2.5f, %2.5f)", nextState.momentum.pT, nextState.momentum.theta, nextState.momentum.phi) << std::endl;
    std::cout << Form("  boundaries  are (r, z) = (%2.5f, %2.5f)", maxR, maxZ) << std::endl << std::flush;
    assert(false);
  }

  double weightNext = 0.;
  if(fabs(currentZ) > maxZ){
    weightNext = (fabs(nextState.position.z) - maxZ)/(fabs(nextState.position.z - initState.position.z));
  } else if( currentRadius > maxR ) {
    weightNext = (nextState.position.r - maxR)/(fabs(nextState.position.r - initState.position.r));
  } else {
    std::cout << "FAILED to extrpolate to (" << maxR << ", " << maxZ << ")" << std::endl; 
    return nextState; 
  }
  double weightInit = 1.-weightNext;

  nextState = nextState * weightNext + initState * weightInit;
  
  return nextState; 

}





TH3D* TrackModel::trajectoryHistogram( const TrackModel::HelixParameters& parameters, std::string instanceName ) const
{

  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }
      
  std::string histogramName = "trajectory_" + instanceName;
  
  TH3D* trajectory = new TH3D( histogramName.c_str(), "", 200, -0.0, 2, 120, -0.6, 0.6, 160, -0.4, 1.2  );
  

  int vectorSize = (int) m_phiVector.size();
  int index(0);
  
  for ( ; index < vectorSize ; ++index ){


    trajectory->Fill( m_zVector[index], m_xVector[index], m_yVector[index]);
  } 
  trajectory->SetMarkerStyle(20);
  trajectory->SetMarkerSize(0.5);
  
  trajectory->GetXaxis()->SetTitle("z [m]");
  trajectory->GetXaxis()->SetTitleOffset(1.9);
  trajectory->GetXaxis()->SetNdivisions( 305 );
    
  trajectory->GetYaxis()->SetTitle("x [m]");
  trajectory->GetYaxis()->SetTitleOffset(1.9);
  trajectory->GetYaxis()->SetNdivisions( 305 );
  
  trajectory->GetZaxis()->SetTitle("y [m]");
  trajectory->GetZaxis()->SetNdivisions( 305 );
  //trajectory->GetZaxis()->CenterTitle();
    
  return trajectory;
  
}
  

TGraph* TrackModel::phiVsR( const TrackModel::HelixParameters& parameters ) const
{
  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }
  
  
  int vectorSize = (int) m_phiVector.size();
  
  double phiArray[vectorSize];
  double rArray[vectorSize];
  
  std::vector<double>::const_iterator phi      = m_phiVector.begin();
  std::vector<double>::const_iterator r        = m_rVector.begin();
  
  int index(0);
  
  for ( ; index < vectorSize ; ++index, ++phi, ++r ){
  
    phiArray[index]      = *phi;
    rArray[index]        = *r;
    
  }
  
  TGraph* graph = new TGraph( vectorSize, rArray, phiArray );
  
  graph->SetTitle("");
  
  graph->GetXaxis()->SetTitle("R [m]");
  //graph->GetXaxis()->SetTitleOffset(1.9);
    
  graph->GetYaxis()->SetTitle("#phi [rad]");
  //graph->GetYaxis()->SetTitleOffset(1.9);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(4);
  graph->SetLineWidth(2);
  
  //TF1* fit = new TF1( "fit", "[0] + [1]*x", 0., 1.2 );
  //fit->SetLineWidth(1);
  //fit->SetLineColor(8);
  //graph->Fit("fit", "Q");
  
  //double slope = fit->GetParameter(1);
  
  //double slope(0.);
  
  //std::cout << "Slope normalised to q/p: " << slope / parameters.qOverPT << std::endl;
      
  return graph;
    
}

TGraph* TrackModel::phiVsCotTheta( const TrackModel::HelixParameters& parameters ) const
{
  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }

  
  int vectorSize = (int) m_cotThetaVector.size();

  double cotThetaArray[vectorSize];
  double phiArray[vectorSize];

  std::vector<double>::const_iterator cotTheta = m_cotThetaVector.begin();
  std::vector<double>::const_iterator phi      = m_phiVector.begin();

  int index(0);

  for ( ; index < vectorSize ; ++index, ++cotTheta, ++phi ){
  
    cotThetaArray[index] = *cotTheta;
    phiArray[index]      = *phi;
    
  }


  TGraph* graph = new TGraph( vectorSize, cotThetaArray, phiArray );
  graph->GetXaxis()->SetLimits( parameters.cotTheta - 0.2 * parameters.cotTheta, parameters.cotTheta + 0.2 * parameters.cotTheta );

  graph->SetTitle("");
  
  graph->GetXaxis()->SetTitle( "cot( #theta_{t} )" );
  graph->GetXaxis()->SetTitleOffset(1.9);
  
  graph->GetYaxis()->SetTitle("#varphi_{t}");
  graph->GetYaxis()->SetTitleOffset(1.9);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(6);
  graph->SetLineWidth(2);
    
  return graph;
  
}

TGraph* TrackModel::cotThetaVsRho( const TrackModel::HelixParameters& parameters ) const
{

  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }

  
  int vectorSize = (int) m_cotThetaVector.size();
  
  double cotThetaArray[vectorSize];
  double rArray[vectorSize];
  
  std::vector<double>::const_iterator cotTheta = m_cotThetaVector.begin();
  std::vector<double>::const_iterator r        = m_rVector.begin();
  
  int index(0);
  
  for ( ; index < vectorSize ; ++index, ++cotTheta, ++r ){
    
    cotThetaArray[index] = *cotTheta;
    rArray[index]        = *r;
  
    //std::cout << cotThetaArray[index] << '\t' << rArray[index] << std::endl;
    
  }
  
  TGraph* graph = new TGraph( vectorSize, cotThetaArray, rArray );
  graph->GetXaxis()->SetLimits( parameters.cotTheta - 0.2 * parameters.cotTheta, parameters.cotTheta + 0.2 * parameters.cotTheta );
  
  graph->SetTitle("");
  
  graph->GetXaxis()->SetTitle( "cot( #theta_{t} )" );
  graph->GetXaxis()->SetTitleOffset(1.9);
  
  graph->GetYaxis()->SetTitle("#rho");
  graph->GetYaxis()->SetTitleOffset(1.9);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(6);
  graph->SetLineWidth(1);
    
  return graph;
    
}


TGraph* TrackModel::etaVsR( const TrackModel::HelixParameters& parameters ) const
{

  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }

  
  int vectorSize = (int) m_cotThetaVector.size();
  
  double etaArray[vectorSize];
  double rArray[vectorSize];
  
  std::vector<double>::const_iterator eta = m_etaVector.begin();
  std::vector<double>::const_iterator r   = m_rVector.begin();
  
  int index(0);
  
  for ( ; index < vectorSize ; ++index, ++eta, ++r ){
    
    etaArray[index]      = *eta;
    rArray[index]        = *r;
  
    //std::cout << cotThetaArray[index] << '\t' << rArray[index] << std::endl;
    
  }
  
  TGraph* graph = new TGraph( vectorSize, etaArray, rArray );
  
  graph->SetTitle("");
  
  graph->GetXaxis()->SetTitle( "#eta_{t}" );
  //graph->GetXaxis()->SetTitleOffset(1.9);
  
  graph->GetYaxis()->SetTitle("R [m]");
  //graph->GetYaxis()->SetTitleOffset(1.9);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(6);
  graph->SetLineWidth(2);
    
  return graph;
    
}


TGraph* TrackModel::etaVsZ( const TrackModel::HelixParameters& parameters ) const
{

  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }

  
  int vectorSize = (int) m_cotThetaVector.size();
  
  double etaArray[vectorSize];
  double rArray[vectorSize];
  
  std::vector<double>::const_iterator eta = m_etaVector.begin();
  std::vector<double>::const_iterator r   = m_zVector.begin();
  
  int index(0);
  
  for ( ; index < vectorSize ; ++index, ++eta, ++r ){
    
    etaArray[index]      = *eta;
    rArray[index]        = *r;
  
    //std::cout << cotThetaArray[index] << '\t' << rArray[index] << std::endl;
    
  }
  
  TGraph* graph = new TGraph( vectorSize, etaArray, rArray );
  
  graph->SetTitle("");
  
  graph->GetXaxis()->SetTitle( "#eta_{t}" );
  //graph->GetXaxis()->SetTitleOffset(1.9);
  
  graph->GetYaxis()->SetTitle("z [m]");
  //graph->GetYaxis()->SetTitleOffset(1.9);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(6);
  graph->SetLineWidth(2);
    
  return graph;
    
}

TGraph* TrackModel::RVsZ( const TrackModel::HelixParameters& parameters ) const
{

  if (m_cachedparameters != &parameters ) {
    std::cout << "Recalculating" << std::endl;
    m_cachedparameters = &parameters;
    BuildTrajectory();
  }

  
  int vectorSize = (int) m_cotThetaVector.size();
  
  double zArray[vectorSize];
  double rArray[vectorSize];
  
  std::vector<double>::const_iterator z = m_zVector.begin();
  std::vector<double>::const_iterator r = m_rVector.begin();
  
  int index(0);
  
  for ( ; index < vectorSize ; ++index, ++z, ++r ){
    
    zArray[index]  = *z;
    rArray[index]  = *r;
  
    //std::cout << cotThetaArray[index] << '\t' << rArray[index] << std::endl;
    
  }
  
  TGraph* graph = new TGraph( vectorSize, zArray, rArray );
  
  graph->SetTitle("");
  
  graph->GetXaxis()->SetTitle( "Z [m]" );
  //graph->GetXaxis()->SetTitleOffset(1.9);
  
  graph->GetYaxis()->SetTitle("R [m]");
  //graph->GetYaxis()->SetTitleOffset(1.9);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(6);
  graph->SetLineWidth(2);
    
  return graph;
    
}



  //<! method to get the trajectory histogram
  TH3D* TrackModel::trajectoryHistogram( std::string instanceName )const{
    if(!m_cachedparameters)
      return 0;
    return trajectoryHistogram(*m_cachedparameters, instanceName);
  }
  
  //!< method to generate an r / phi graph
  TGraph* TrackModel::phiVsR() const{
    if(!m_cachedparameters)
      return 0;
    return phiVsR(*m_cachedparameters);
  }
  
  
  
  //!< Method to generate a cotTheta / phi graph
  TGraph* TrackModel::phiVsCotTheta() const{
    if(!m_cachedparameters)
      return 0;
    return phiVsCotTheta(*m_cachedparameters);
  }
  
  //!< Method to generate a cot Theta / rho graph
  TGraph* TrackModel::cotThetaVsRho() const{
    if(!m_cachedparameters)
      return 0;
    return cotThetaVsRho(*m_cachedparameters);
  }
  
  //!< method to generate an r / eta graph
  TGraph* TrackModel::etaVsR() const{
    if(!m_cachedparameters)
      return 0;
    return etaVsR(*m_cachedparameters);
  }

  //!< method to generate an r / eta graph
  TGraph* TrackModel::etaVsZ() const{
    if(!m_cachedparameters)
      return 0;
    return etaVsZ(*m_cachedparameters);
  }

  //!< method to generate an r / z graph
  TGraph* TrackModel::RVsZ() const{
    if(!m_cachedparameters)
      return 0;
    return RVsZ(*m_cachedparameters);
  }
