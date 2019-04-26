#include <vector>
#include <xAODTracking/TrackParticleContainer.h>

struct LinearTrackParameters
{
  AmgVector(5)    parametersAtPCA;
  AmgSymMatrix(5) covarianceAtPCA;
  AmgVector(3)     linearizationPoint;
  AmgMatrix(5, 3) positionJacobian;
  AmgMatrix(5, 3)  momentumJacobian;
  AmgVector(3)     positionAtPCA;
  AmgVector(3)     momentumAtPCA;
  AmgVector(5)   constantTerm;
};

class Vertex
{
 public:
  Vertex(const AmgVector(3)& position, const AmgSymMatrix(3)& covariance):
  m_position(position),
  m_covariance(covariance)
  {};

  AmgVector(3) position() const {return m_position;};
  AmgSymMatrix(3) covariance() const {return m_covariance;};
  void setFitQuality( double chi2, double dof ){ m_chi2 = chi2; m_dof = dof;};
  double chi2()const{return m_chi2;};
  double dof()const{return m_dof;};

 private:
  AmgVector(3) m_position;
  AmgSymMatrix(3) m_covariance;

  double m_chi2;
  double m_dof;
};


class SimpleVertexFit
{
public:
  Vertex fitVertex( std::vector<const xAOD::TrackParticle*>& , const AmgVector(3)& initialPosition);
private:
  Vertex updateVertex(const Vertex& vertexInit,const  LinearTrackParameters&, double trackWeight,int sign );
  double trackParametersChi2(const Vertex& new_vtx, const LinearTrackParameters& trk) ;
  double vertexPositionChi2(const Vertex& old_vtx, const Vertex& new_vtx);

  LinearTrackParameters lineariseTrack(const xAOD::TrackParticle&, const AmgVector(3)& );
};
