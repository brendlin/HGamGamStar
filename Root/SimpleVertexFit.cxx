#include "HGamGamStar/SimpleVertexFit.h"



Vertex SimpleVertexFit::fitVertex( std::vector<const xAOD::TrackParticle*>& tracks, const AmgVector(3) & initialPosition){

  std::vector<LinearTrackParameters> linTracks;
  for(const auto& track:  tracks)
   linTracks.push_back(lineariseTrack( *track, initialPosition));

  AmgSymMatrix(3) cov;
  cov <<100,0,0,
        0,100,0,
        0,0,100;
  Vertex vertex( initialPosition, cov);
  vertex.setFitQuality(0,0);

  for(const auto trkLin:  linTracks){
    vertex = updateVertex( vertex, trkLin, 1, 1);
  }

  return vertex;
}

Vertex SimpleVertexFit::updateVertex( const Vertex& vtx, const LinearTrackParameters& trk, double trackWeight,int sign ){

   //  std::cout<<"Position update called " <<std::endl;
   //  std::cout<<"Position update: position Jacobian: "<<trk->positionJacobian()<<std::endl;
   //  std::cout<<"Position update: momentum Jacobian: "<<trk->momentumJacobian()<<std::endl;
   //  std::cout<<"Expected parameters: "<<trk->expectedParametersAtPCA()<<std::endl;

   // linearized track information
   const AmgMatrix(5,3) A = trk.positionJacobian;
   const AmgMatrix(5,3) B = trk.momentumJacobian;
   const AmgVector(5)  trackParameters = trk.parametersAtPCA;
   const AmgVector(5)  constantTerm = trk.constantTerm;
   const AmgSymMatrix(5)  trackParametersWeight  = trk.covarianceAtPCA.inverse();
   // const double trackWeight = trk->weight();

   //vertex to be updated
   const AmgVector(3)  old_pos = vtx.position();
   const AmgSymMatrix(3) old_vrt_weight = vtx.covariance().inverse();
   //ATH_MSG_VERBOSE ("Old vertex weight " << old_vrt_weight);

   //making the intermediate quantities:
   //W_k = (B^T*G*B)^(-1)
   //  std::cout << "B matrix calculated: " << B << std::endl;

   AmgSymMatrix(3) S = B.transpose()*(trackParametersWeight*B);
   //  std::cout << "S matrix calculated " << S << std::endl;
   S = S.inverse().eval();

   //  std::cout << "Main inversion passed successfully" << std::endl;
   //  std::cout << "S matrix inverted" << S << std::endl;
   //G_b = G_k - G_k*B_k*W_k*B_k^(T)*G_k
   AmgSymMatrix(5) gB = trackParametersWeight - trackParametersWeight*(B*(S*B.transpose()))*trackParametersWeight.transpose();
   //  std::cout << "Gain factor obtained: " << trackParametersWeight*(B*(S*B.transpose()))*trackParametersWeight.transpose() << std::endl;
   //  std::cout << "Resulting Gain Matrix: " << gB << std::endl;

   //new vertex weight matrix
   AmgSymMatrix(3) new_vrt_weight_later_cov = old_vrt_weight + trackWeight * sign * A.transpose()*(gB*A);
   //  std::cout << "Track weight: " << trackWeight << std::endl;
   //  std::cout << "Gain similarity: " << A.transpose()*(gB*A) << std::endl;

   //  std::cout << "New vertex weight obtained: " << new_vrt_weight << std::endl;

   new_vrt_weight_later_cov = new_vrt_weight_later_cov.inverse().eval();

   //  std::cout << "New vertex covariance obtained: " << new_vrt_cov << std::endl;
   //new vertex position
   AmgVector(3) new_vrt_position =  new_vrt_weight_later_cov*(old_vrt_weight * old_pos + trackWeight * sign * A.transpose() * gB *(trackParameters - constantTerm) );
   //  std::cout << "New vertex position obtained: " << new_vrt_position << std::endl;

   // return a vertex which has an updated position and covariance but with no tracks
   Vertex r_vtx(new_vrt_position,new_vrt_weight_later_cov);
   double chi2 = vtx.chi2() + sign * ( vertexPositionChi2(vtx, r_vtx) + trackWeight* trackParametersChi2( r_vtx, trk )  );
   double dof = vtx.dof() +  sign * trackWeight * (2.0);
   r_vtx.setFitQuality( chi2, dof );
   // r_vtx was a RecVertex before vertex EDM migration and instantiated using RecVertex(pos,cov)

   //ATH_MSG_VERBOSE ( "Maths done, returning a valid xAOD::Vertex.");

   return r_vtx;
}


//xAOD chi2 increment method
double SimpleVertexFit::trackParametersChi2(const Vertex& new_vtx, const LinearTrackParameters& trk)
{

  const AmgVector(3)&  new_vrt_pos  = new_vtx.position();

  // track information
  const AmgMatrix(5,3)& A = trk.positionJacobian;
  const AmgMatrix(5,3) & B = trk.momentumJacobian;
  const AmgVector(5) & trkParameters = trk.parametersAtPCA;
  const AmgSymMatrix(5) & trkParametersWeight = trk.covarianceAtPCA.inverse();

  //calculation of S matrix
  AmgSymMatrix(3) Sm = B.transpose()*(trkParametersWeight*B);

  Sm = Sm.inverse().eval();
  AmgVector(5) theResidual = trk.constantTerm;

  //refitted track momentum
  AmgVector(3) newTrackMomentum = Sm*B.transpose()*trkParametersWeight*(trkParameters - theResidual - A*new_vrt_pos);

  //refitted track parameters
  AmgVector(5) refTrackParameters = theResidual + A * new_vrt_pos + B * newTrackMomentum;

  //parameters difference
  AmgVector(5) paramDifference = trkParameters - refTrackParameters;

  double chi2 = paramDifference.transpose() * ( trkParametersWeight * paramDifference );

  //here the momentum part of the chi2 is calculated;
  //making the vertex position chi2

  return chi2;

}//end of chi2 increment method

double SimpleVertexFit::vertexPositionChi2(const Vertex& old_vtx, const Vertex& new_vtx)
{
  AmgSymMatrix(3)  old_wrt_weight = old_vtx.covariance().inverse();
  AmgVector(3) posDifference = new_vtx.position() - old_vtx.position();
  double chi2 = posDifference.transpose()*(old_wrt_weight*posDifference);
  return chi2;
}//end of vertex position chi2


LinearTrackParameters SimpleVertexFit::lineariseTrack(const xAOD::TrackParticle& track, const AmgVector(3)&  linPoint)
{
  AmgVector(5)    paramsAtPCA = track.definingParameters();
  AmgSymMatrix(5) parCovarianceAtPCA =  track.definingParametersCovMatrix();
  AmgVector(3)   positionAtPCA =  linPoint;

  // phiV and functions
  double phiV    = paramsAtPCA[2];
  double sinPhiV = std::sin(phiV);
  double cosPhiV = std::cos(phiV);

  positionAtPCA[0] += paramsAtPCA[0] * sinPhiV;
  positionAtPCA[1] += paramsAtPCA[0] * cosPhiV;
  positionAtPCA[2] += paramsAtPCA[1];

  // theta and functions
  double th    = paramsAtPCA[3];
  double sinTh = std::sin(th);
  double tanTh = std::tan(th);

  // q over p
  double qOvP = paramsAtPCA[4];
  double sgnH = (qOvP < 0.) ? -1 : 1;

  AmgVector(3) momentumAtPCA(phiV, th, qOvP);

  // get B-field z-component at current position
  double Bz = 2;

  double rho;
  // Curvature is infinite w/o b field
  if (Bz == 0. || std::abs(qOvP) < 1.e-15) {
    rho = std::numeric_limits<double>::max();
  } else {
    // signed(!) rho
    rho = sinTh * 0.299792 / (qOvP * Bz);
  }

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X  = positionAtPCA(0) - linPoint.x() + rho * sinPhiV;
  double Y  = positionAtPCA(1) - linPoint.y() - rho * cosPhiV;
  double S2 = (X * X + Y * Y);
  double S  = std::sqrt(S2);

  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  AmgVector(5) predParamsAtPCA;

  int sgnX = (X < 0.) ? -1 : 1;
  int sgnY = (Y < 0.) ? -1 : 1;

  double phiAtPCA;
  if (std::abs(X) > std::abs(Y)) {
    phiAtPCA = sgnH * sgnX * std::acos(-sgnH * Y / S);
  } else {
    phiAtPCA = std::asin(sgnH * X / S);
    if ((sgnH * sgnY) > 0) {
      phiAtPCA = sgnH * sgnX * M_PI - phiAtPCA;
    }
  }

  // Eq. 5.33 in Ref(1) (see .hpp)
  predParamsAtPCA[0] = rho - sgnH * S;
  predParamsAtPCA[1]
      = positionAtPCA[2] - linPoint.z() + rho * (phiV - phiAtPCA) / tanTh;
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = th;
  predParamsAtPCA[4] = qOvP;

  // Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  AmgMatrix(5, 3) positionJacobian;
  positionJacobian.setZero();
  // First row
  positionJacobian(0, 0) = -sgnH * X / S;
  positionJacobian(0, 1) = -sgnH * Y / S;

  // Second row
  positionJacobian(1, 0) = rho * Y / (tanTh * S2);
  positionJacobian(1, 1) = -rho * X / (tanTh * S2);
  positionJacobian(1, 2) = 1.;

  // Third row
  positionJacobian(2, 0) = -Y / S2;
  positionJacobian(2, 1) = X / S2;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  AmgMatrix(5, 3) momentumJacobian;
  momentumJacobian.setZero();

  double R    = X * cosPhiV + Y * sinPhiV;
  double Q    = X * sinPhiV - Y * cosPhiV;
  double dPhi = phiAtPCA - phiV;

  // First row
  momentumJacobian(0, 0) = -sgnH * rho * R / S;

  double qOvSred = 1 - sgnH * Q / S;

  momentumJacobian(0, 1) = qOvSred * rho / tanTh;
  momentumJacobian(0, 2) = -qOvSred * rho / qOvP;

  // Second row
  momentumJacobian(1, 0) = (1 - rho * Q / S2) * rho / tanTh;
  momentumJacobian(1, 1) = (dPhi + rho * R / (S2 * tanTh * tanTh)) * rho;
  momentumJacobian(1, 2) = (dPhi - rho * R / S2) * rho / (qOvP * tanTh);

  // Third row
  momentumJacobian(2, 0) = rho * Q / S2;
  momentumJacobian(2, 1) = -rho * R / (S2 * tanTh);
  momentumJacobian(2, 2) = rho * R / (qOvP * S2);

  // Last two rows:
  momentumJacobian(3, 1) = 1.;
  momentumJacobian(4, 2) = 1.;

  // const term F(V_0, p_0) in Talyor expansion
  AmgVector(5) constTerm = predParamsAtPCA - positionJacobian * positionAtPCA
      - momentumJacobian * momentumAtPCA;

  LinearTrackParameters ltp;
  ltp.parametersAtPCA = paramsAtPCA;
  ltp.covarianceAtPCA = parCovarianceAtPCA;
  ltp.linearizationPoint = linPoint;
  ltp.positionJacobian =  positionJacobian;
  ltp.momentumJacobian = momentumJacobian;
  ltp.positionAtPCA = positionAtPCA;
  ltp.momentumAtPCA = momentumAtPCA;
  ltp.constantTerm = constTerm;

  return ltp;

}
