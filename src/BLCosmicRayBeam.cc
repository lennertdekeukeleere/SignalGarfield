#include "BLCosmicRayBeam.hh"
#include "Garfield/Random.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <assert.h>

CosmicRayBeam::CosmicRayBeam()
{
	meanMomentum = 3.0;
	beamY = 0.0;
	boxX =0.0;
	boxY=0.0;
	beamWidth = 0.0;
	beamHeight = 0.0;
	beamLength = 0.0;
	boxWidth=0.0;
	boxLength=0.0;
    boxHeight=0.0;
	radius = 0.0;
	sterradians = 0.0;
	hitsPerM2PerSecPerSterrad = 0.0;
	particle = "mu+";
}


void CosmicRayBeam::init(double gasboxX,double gasboxY,double gasboxZ,double gasboxOffset, double muonPlane)
{
	beamWidth = gasboxX;
	beamLength = gasboxY;
    beamHeight = gasboxZ;
    beamY = gasboxOffset;
	
	
	boxWidth=beamWidth;

    boxHeight=beamHeight;
    boxLength=0;
    boxY=muonPlane;
	
	radius = sqrt(beamWidth*beamWidth/4+beamLength*beamLength/4+beamHeight*beamHeight/4);
}


void CosmicRayBeam::nextBeamEvent(double& xIni,double& yIni,double& zIni,double& dxIni, double& dyIni,double& dzIni,double& ekinIni)
{
    std::cout << "CosmicRayBeam::nextBeamEvent" << std::endl;
	double mass = 0.105658367*std::pow(10,9);
	double position[3];
	double direction[3];
	double time = 0.0;
	double momentum = 0.0;

	// work in local coordinates (relative to the beam box)
	double x = boxWidth*Garfield::RndmUniform() - boxWidth/2.0+boxX;
	double z = boxHeight*Garfield::RndmUniform() - boxHeight/2.0;
	double y = boxLength*Garfield::RndmUniform() - boxLength/2.0+boxY;
    std::cout << "cosmicRayMuonMomentum" << std::endl;
	momentum = cosmicRayMuonMomentum();
    std::cout << "cosmicRayMuonAngle" << std::endl;
	double theta = cosmicRayMuonAngle();
	double phi = 2.0*CLHEP::pi*Garfield::RndmUniform();
	direction[0] = momentum*sin(theta)*cos(phi);
	direction[2] = momentum*sin(theta)*sin(phi);
	direction[1] = -momentum*cos(theta);


//	if(fabs(x) < beamWidth/2.0 && fabs(z) < beamHeight/2.0)
//		++hits;

	double X[5];
	double Y[5];
	double Z[5];
	double R[5];
 
	X[1] = beamWidth/2.;
	X[2] = -beamWidth/2.;
	Y[0] = beamLength/2.+beamY;
	Z[3] = beamHeight/2.;
	Z[4] = -beamHeight/2.;
	
	Y[1] = y+(X[1]-x)*direction[1]/direction[0];
	Z[1] = z+(Y[1]-y)*direction[2]/direction[1];
	
	Y[2] = y+(X[2]-x)*direction[1]/direction[0];
	Z[2] = z+(Y[2]-y)*direction[2]/direction[1];
	
	X[0] = x+(Y[0]-y)*direction[0]/direction[1];
	Z[0] = z+(Y[0]-y)*direction[2]/direction[1];

	Y[3] = y+(Z[3]-z)*direction[1]/direction[2];
	X[3] = x+(Y[3]-y)*direction[0]/direction[1];

	Y[4] = y+(Z[4]-z)*direction[1]/direction[2];
	X[4] = x+(Y[4]-y)*direction[0]/direction[1];		
	
	
	if(fabs(X[0])<=beamWidth/2. && fabs(Z[0])<=beamHeight/2.){
		position[0] = X[0];
		position[2] = Z[0];
		position[1] = Y[0];
	}
	else{
		y=Y[0]-beamLength;
		int ind=1;
		for(int i=1;i<5;i++){	
			if(Y[i]>=y && Y[i]<Y[0] && fabs(X[i])<=beamWidth/2. && fabs(Z[i])<=beamHeight/2.){
				y=Y[i];
				ind=i;
			}
		}
		position[0] = X[ind];
		position[2] = Z[ind];
		position[1] = Y[ind];
	}
		
	
    
    double ke = sqrt(momentum*momentum + mass*mass) - mass;
    
    std::cout << "Muon generator" << std::endl;
    std::cout << "position: " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;
    std::cout << "time: " << time << std::endl;
    std::cout << "energy: " << ke << std::endl;
    std::cout << "direction: " << direction[0] << ", " << direction[1] << ", " << direction[2] << std::endl;

    xIni = position[0];
    yIni = position[1];
    zIni = position[2];

    double dirMag = std::sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
    dxIni = direction[0]/dirMag;
    dyIni = direction[1]/dirMag;
    dzIni = direction[2]/dirMag;

    ekinIni = ke;
    
}


double CosmicRayBeam::cosmicRayMuonMomentum()
{
	// Data from Kremer et al, Phys. Rev. Lett., Vol 83 no 21, p4241 (1999).
	// values are lower bin edge, bin average, mu+ rate, mu- rate
	// (laid out this silly way so verification with the paper is easy)
	// NOTE: units are GeV/c, and counts/(GeV/c m^2 sr s)
	static double vals[] = {
		0.0,	0.0,	0.0,	0.0,
		0.2,	0.25,	14.0,	11.0,
		0.3,	0.35,	16.8,	13.6,
		0.4,	0.47,	17.2,	14.4,
		0.55,	0.62,	16.6,	13.5,
		0.70,	0.78,	15.6,	13.3,
		0.85,	0.92,	14.8,	12.1,
		1.0,	1.1,	13.0,	11.0,
		1.2,	1.3,	12.0,	10.1,
		1.4,	1.5,	10.2,	8.7,
		1.6,	1.84,	9.1,	7.3,
		2.1,	2.49,	6.6,	5.2,
		2.94,	3.49,	4.12,	3.38,
		4.12,	4.78,	2.53,	1.98,
		5.5,	6.21,	1.61,	1.25,
		7.0,	8.37,	0.90,	0.69,
		10.0,	12.42,	0.389,	0.309,
		15.5,	18.85,	0.138,	0.108,
		23.0,	26.68,	0.063,	0.046,
		31.1,	36.69,	0.028,	0.019,
		43.6,	51.47,	0.0099,	0.0071,
		61.1,	72.08,	0.0036,	0.0030,
		85.6,	100.96,	0.0014,	0.0012,
		120.0,	120.0,	0.0,	0.0}; // cutoff at 120 GeV/c
	const int nvals = sizeof(vals)/sizeof(vals[0]);
	const int nbins = nvals/4 - 1;
	const int npdf=256;
	static double pdf[npdf];
	static double pmax = vals[4*nbins];
	static bool init=true;

	if(init) {
		// RandGeneral needs equal-sized bins for pdf[]
		// it returns a value in the range [0,1)
		hitsPerM2PerSecPerSterrad = 0.0;
		for(int i=0,ibin=0; i<npdf; ++i) {
			double p = (i+0.5)*pmax/npdf;
			while(p >= vals[4*ibin+5]) ++ibin;
			assert(ibin <= nbins);
			double f = (p - vals[4*ibin+1]) /
						(vals[4*ibin+5]-vals[4*ibin+1]);
			assert(0.0 <= f && f <= 1.0);
			pdf[i] = (1.0-f)*(vals[4*ibin+2]+vals[4*ibin+3]) +
					f*(vals[4*ibin+6]+vals[4*ibin+7]);
			hitsPerM2PerSecPerSterrad += pdf[i] * pmax/npdf;
		}
		init = false;
	}

	CLHEP::RandGeneral generator(pdf,npdf); // BUG in RandGeneral - cannot use new
	return generator.shoot() * pmax * std::pow(10,9);
}

double CosmicRayBeam::cosmicRayMuonAngle()
{
	const int npdf=128;
	static double pdf[npdf];
	const double thetamax = 70.0*(CLHEP::pi/180.0);
	static bool init = true;

	if(init) {
		// RandGeneral needs equal-sized bins for pdf[]
		// it returns a value in the range [0,1)
		sterradians = 0.0;
		double dtheta = thetamax / npdf;
		for(int i=0; i<npdf; ++i) {
			// Particle Data Group, Review of Particle Properties,
			// 2002. Section 23.3.1.
			double c = cos(dtheta*i);
			pdf[i] = c*c;
			sterradians += 2.0*CLHEP::pi*c*c*sin(dtheta*i)*dtheta;
		}
		init = false;
	}

	CLHEP::RandGeneral generator(pdf,npdf); // BUG in RandGeneral - cannot use new
	return generator.shoot() * thetamax;
}
