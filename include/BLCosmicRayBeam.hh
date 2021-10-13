#ifndef COSMICRAYBEAM_HH
#define COSMICRAYBEAM_HH

#include <string>
#include <iostream>


class CosmicRayBeam {
  public:
	/// Constructor.
	CosmicRayBeam();
    ~CosmicRayBeam(){};
	/// init() will initialize internal variables.
	void init(double,double,double,double,double);

	/// generateReferenceParticle() generates the reference particle.

	/// nextBeamEvent() generates the next beam event.
	void nextBeamEvent(double& x,double& y,double& z,double& dx, double& dy,double& dz,double& ekin);

	/// cosmicRayMuonMomentum() returns a random momentum value distributed
	/// like the muons from cosmic rays. Average momentum is 3.0*GeV/c, 
	/// cutoff is 120*GeV/c. Return value is in geant4 units (MeV).
	/// Data from Kremer et al, Phys. Rev. Lett., Vol 83 no 21, p4241 (1999)
	double cosmicRayMuonMomentum();

	/// cosmicRayMuonAngle() returns a random angle distributed like the 
	/// polar angle of cosmic ray muons.
	/// The return value is the polar angle from vertical (radians). It is
	/// cut off at 70 degrees.
	/// Note this is the distribution for muons of ~3 GeV/c, which is the
	/// average momentum. In fact, lower energy muons have a steeper
	/// distribution and higher ones have a flatter distribution.
	/// But this is a reasonable approximation.
	/// Particle Data Group, Review of Particle Properties, 2002.
	/// Section 23.3.1.
	double cosmicRayMuonAngle();
    
  private:

    double meanMomentum;
    double beamY,boxX,boxY;
    double beamWidth;
    double beamHeight;
    double beamLength;
	double boxWidth;
	double boxLength;
    double boxHeight;
    double radius;
    double sterradians;
    double hitsPerM2PerSecPerSterrad;
    std::string particle;
};

#endif
