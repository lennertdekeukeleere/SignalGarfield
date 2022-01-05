#ifndef SIGNALGENERATOR_H
#define SIGNALGENERATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>

#include "Garfield/Sensor.hh"
#include "Garfield/MediumMagboltz.hh"

#include "Garfield/GeometrySimple.hh"		//Geometry
#include "Garfield/SolidBox.hh"			//Geometry
#include "Garfield/ComponentAnalyticField.hh"	//Garfield field
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewCell.hh"			//Visualization
#include "TCanvas.h"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/TrackHeed.hh"


class AnalysisManager;
class CosmicRayBeam;

class SignalGenerator {
  public:
    SignalGenerator(AnalysisManager*, std::string);
    ~SignalGenerator();

    void Initialise();
    void Run();

    inline void SetGasFile(std::string s) { gasFile = s;};
    inline void SetIonMobilityFile(std::string s) { ionMobFile = s; };

    inline void SetLoadComsolField(bool comsol){loadComsolfield = comsol;};
    inline void SetBinWidth(double binwidth){fBinWidth = binwidth;};
    inline void SetNBins(int nbins){fNbins = nbins;};
    inline double GetRiseTime(){return fRiseTime;};
    inline double GetAvalancheSize(){return fAvalancheSize/fNumOfElectrons;};
    inline void SetTrackMicroscopic(bool b){trackMicro=b;};
    inline void SetCreateAvalancheMC(bool b){createAval=b;};
    inline void SetVisualizeChamber(bool b){visualizeChamber = b;};
    inline void SetVisualizeSignals(bool b){visualizeSignal = b;};
    inline void SetVisualizeField(bool b){visualizeField = b;};
    inline void SetDriftRKF(bool b){driftRKF=b;};

    void SetTransferFile(std::string name);
    void SetComsolFile(std::string name){comsolFile=name;};

    inline double GetGasBoxSizeX(){return gasboxhalfX*2;};
    inline double GetGasBoxSizeY(){return gasboxhalfY*2;};
    inline double GetGasBoxSizeZ(){return gasboxhalfZ*2;};
    inline double GetGasBoxOffset(){return offset;};

  private:
    void ProcessEvent(float,float,int,int);
    void GetSignalWirePosition(int,double&,double&);

    CosmicRayBeam* muonGen;
    AnalysisManager* anManager;
    Garfield::MediumMagboltz* fMediumMagboltz;
    Garfield::Sensor* fSensor;
    //  Garfield::TrackHeed* fTrackHeed;
    Garfield::TrackHeed* fTrackHeed;
    Garfield::GeometrySimple* geo;
    Garfield::SolidBox* box; //Full chamber with 5 planes
    Garfield::ComponentAnalyticField* comp;
    Garfield::AvalancheMC* fDrift;
    Garfield::DriftLineRKF* fDriftRKF;
    Garfield::AvalancheMicroscopic* fAvalanche;
    TCanvas* fChamber;
    TCanvas* fSignal;
    TCanvas* fField;
    Garfield::ViewCell* cellView;
    Garfield::ViewDrift* viewDrift;
    Garfield::ViewSignal* viewSignal;
    Garfield::ViewField* viewField;

    std::vector<double> transferTimeVector;
    std::vector<double> transferValuesVector;

    void makeGas();
    void buildBox();
    void BuildCompField();
    void BuildSensor();
    void SetTracking();
    void CreateChamberView();
    void CreateSignalView();
    void CreateFieldView();
    void ReadConfigFile();
    bool WireIsHit(int,double,double,double);

    bool loadComsolfield;
    bool trackMicro;
    bool createAval;
    bool driftIons;
    bool fullSignal;
    bool visualizeChamber;
    bool visualizeSignal;
    bool visualizeField;
    bool driftRKF;
    bool tracks;
    bool signalCFD;
    bool fullChamber;
    bool muonDist;

    double fBinWidth;
    int fNbins;
    double sigWireVoltage;
    double sigWireVoltageUpper;
    double sigWireVoltageLower;
    double fieldWireVoltage;
    double primE;
    int nSim;
    double fRiseTime;
    double fAvalancheSize;
    int fNumOfElectrons;
    int fNumOfIons;
    int avSize;
    double signalThreshold;
    std::string particleType;
    std::string ionMobFile;
    std::string gasFile;
    std::string comsolFile;
    std::string inFileName;
    std::string transferFilename;
    int nrofPlanes;
    int nrofWiresPlane;
    double Swr = 12.5e-4;
    double Fwr = 50.e-4;
    double cellR = sqrt(3)*0.5;
    double planespacing = 0.5;
    double planeOffset = 2.95;
    double spacerThickness = 0.3;
    double innerEndCapThickness = 4;
    double firstWireXpos;
    double offset;
    double gasboxhalfX;
    double gasboxhalfY;
    double gasboxhalfZ;
    double muonPlane;

    double signalWireVoltages[10];

    int cellNr;

};

double transfer(double t);
#endif /* GARFIELDMODELCONFIG_HH_ */
