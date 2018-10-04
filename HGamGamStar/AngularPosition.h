#ifndef HGamGamStar_AngularPosition
#define HGamGamStar_AngularPosition

#include <iostream>
#include <math.h>
#include <string.h>

#include <TString.h>
#include <TMath.h>


struct AngularPosition {

public:
  double eta, phi;

  AngularPosition()
    : eta(0.0), phi(0.0) {}

  AngularPosition(double thisEta, double thisPhi)
    : eta(thisEta), phi(thisPhi) {}


  double deltaEta(AngularPosition other){
    return(fabs(other.eta - eta));
  }

  double deltaEtaSigned(AngularPosition other){
    return(other.eta - eta);
  }

  double deltaPhi(AngularPosition other){
    double dPhi = fabs(other.phi - phi);
    while (dPhi >= TMath::Pi()) dPhi = fabs(dPhi - 2.0*TMath::Pi());
    return(dPhi);
  }

  double deltaPhiSigned(AngularPosition other){
    double dPhi = other.phi - phi;
    while (dPhi > +TMath::Pi()) dPhi -= 2.0*TMath::Pi();
    while (dPhi < -TMath::Pi()) dPhi += 2.0*TMath::Pi();
    return(dPhi);
  }

  double deltaR(AngularPosition other){
    double dE = deltaEta(other);
    double dP = deltaPhi(other);
    return(sqrt(dE*dE + dP*dP));
  }

  std::string print(){
    return(Form("(%f, %f)", eta, phi));
  }

};


#endif
 
