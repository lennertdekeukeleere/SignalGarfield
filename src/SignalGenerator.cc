#include "SignalGenerator.hh"
#include "AnalysisManager.hh"
#include "OptionContainer.hh"
#include "BLCosmicRayBeam.hh"
#include "G4SystemOfUnits.hh"
#include "Garfield/Medium.hh"

#include <iostream>
#include <sstream>
#include <math.h>
#include <fstream>

SignalGenerator::SignalGenerator(AnalysisManager* an, std::string in) : anManager(an), inFileName(in) {
  fMediumMagboltz=0;
  fSensor=0;
  fTrackHeed=0;
  trackMicro=false;
  createAval=false;
  driftRKF=false;
  visualizeChamber = false;
  visualizeSignal = false;
  debugging = false;
  fBinWidth = 0.;
  fNbins = 0;
  comp=0;
  fAvalanche=0;
  fDrift=0;
  fDriftRKF=0;
  gasFile="";
  nrofPlanes = 10;
  nrofWiresPlane=8;
  transferFilename="";
 }


SignalGenerator::~SignalGenerator() {
  delete fTrackHeed;
  delete fMediumMagboltz;
  delete fSensor;
  delete comp;
  delete fAvalanche;
  delete fDrift;
  delete fDriftRKF;

}


void SignalGenerator::Initialise() {

	ReadConfigFile();

	makeGas();

	buildBox();


	BuildCompField();

	BuildSensor();

	SetTracking();

	if(visualizeChamber) CreateChamberView();
	if(visualizeSignal) CreateSignalView();
	if(visualizeField) CreateFieldView();

	anManager->Initialise(fSensor,fullSignal);
	if(muonDist){
		muonGen = new CosmicRayBeam();
		muonGen->init(gasboxhalfX*2,gasboxhalfY*2,gasboxhalfZ*2,offset,muonPlane);
	}

}

void SignalGenerator::ReadConfigFile(){
	fBinWidth = OptionContainer::GetInstance().GetOption<double>("Signal.Resolution");
	fNbins = OptionContainer::GetInstance().GetOption<double>("Signal.Length")/fBinWidth;
	signalThreshold = OptionContainer::GetInstance().GetOption<double>("Signal.Threshold");
	fullSignal = OptionContainer::GetInstance().GetOption<bool>("Signal.Full");
	visualizeSignal = OptionContainer::GetInstance().GetOption<bool>("Signal.View");
        transferFilename = OptionContainer::GetInstance().GetOption<std::string>("Signal.Transfer");
	fullChamber = OptionContainer::GetInstance().GetOption<bool>("Tracking.FullChamber");
	tracks = OptionContainer::GetInstance().GetOption<bool>("Tracking.Track");
	driftRKF = OptionContainer::GetInstance().GetOption<bool>("Tracking.DriftInt");
	createAval = OptionContainer::GetInstance().GetOption<bool>("Tracking.Avalanche");
	trackMicro = OptionContainer::GetInstance().GetOption<bool>("Tracking.Microscopic");
	driftIons = OptionContainer::GetInstance().GetOption<bool>("Tracking.DriftIons");
	visualizeChamber = OptionContainer::GetInstance().GetOption<bool>("Tracking.View");
	visualizeField = OptionContainer::GetInstance().GetOption<bool>("Field.View");
	if(visualizeField)
		cellNr = OptionContainer::GetInstance().GetOption<int>("Field.PlotProfileCell");
	sigWireVoltage = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV");
	sigWireVoltageUpper = OptionContainer::GetInstance().GetOption<double>("Field.SigWireVUpper");
	sigWireVoltageLower = OptionContainer::GetInstance().GetOption<double>("Field.SigWireVLower");
	nrofPlanes =  OptionContainer::GetInstance().GetOption<int>("Chamber.NumberOfPlanes");
	if(fullChamber){
		signalWireVoltages[0] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV1");
		signalWireVoltages[1] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV2");
		signalWireVoltages[2] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV3");
		signalWireVoltages[3] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV4");
		signalWireVoltages[4] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV5");
		signalWireVoltages[5] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV6");
		signalWireVoltages[6] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV7");
		signalWireVoltages[7] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV8");
		signalWireVoltages[8] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV9");
		signalWireVoltages[9] = OptionContainer::GetInstance().GetOption<double>("Field.SigWireV10");
	}
	fieldWireVoltage = OptionContainer::GetInstance().GetOption<double>("Field.FieldWireV");
	primE = OptionContainer::GetInstance().GetOption<double>("Particle.Energy");
	muonDist = OptionContainer::GetInstance().GetOption<bool>("Particle.CustomEnergy");
	particleType = OptionContainer::GetInstance().GetOption<std::string>("Particle.Type");
	nSim = OptionContainer::GetInstance().GetOption<int>("Particle.Number");
	avSize = OptionContainer::GetInstance().GetOption<int>("Tracking.AvalancheSize");
  debugging = OptionContainer::GetInstance().GetOption<bool>("Tracking.Debugging");
}

void SignalGenerator::makeGas(){
	fMediumMagboltz = new Garfield::MediumMagboltz();
	double pressure = OptionContainer::GetInstance().GetOption<int>("Mixture.Pres");

	double isobutanePerc = OptionContainer::GetInstance().GetOption<int>("Mixture.Iso");
	fMediumMagboltz->SetComposition("isobutane",isobutanePerc,"helium",
					100.-isobutanePerc);
//	fMediumMagboltz->SetComposition("NEO-PENTANE");
//	fMediumMagboltz->SetComposition("4He",100.-isobutanePerc,"CO2",isobutanePerc);
	fMediumMagboltz->SetTemperature(288.15);
	fMediumMagboltz->SetPressure(pressure*0.76);
	fMediumMagboltz->EnableDebugging();
	fMediumMagboltz->Initialise(true);
  fMediumMagboltz->DisableDebugging();
  // Set the Penning transfer efficiency.
	const double rPenning = 1;
	const double lambdaPenning = 0.;
	fMediumMagboltz->EnablePenningTransfer(rPenning, lambdaPenning, "helium");
//	fMediumMagboltz->EnablePenningTransfer(rPenning, lambdaPenning, "NEO-PENTANE");
// Load the ion mobilities.


	std::string gasdata = OptionContainer::GetInstance().GetOption<std::string>("gasdata");
	std::ostringstream s;
	double hePerc = 100. - isobutanePerc;
	s << gasdata << "iso_" << isobutanePerc << ".0%_he_" << hePerc << ".0%_" << pressure << ".0mbar_penning.gas";
//	s << gasdata << "4He_" << hePerc << ".0%_CO2_" << isobutanePerc << ".0%_" << pressure << ".0mbar_penning.gas";
//	s << gasdata << "NEO-PENTANE_100.0%_" << pressure << ".0mbar_penning.gas";
	std::string gasFile = s.str();
	std::cout << gasFile << std::endl;
	if(gasFile!="")
    	fMediumMagboltz->LoadGasFile(gasFile.c_str());
    if(driftIons){
  		const std::string path = getenv("GARFIELD_HOME");
  		fMediumMagboltz->LoadIonMobility(path + "/Data/IonMobility_He+_He.txt");
  	}
}

void SignalGenerator::buildBox(){
	geo = new Garfield::GeometrySimple();
	if(fullChamber){
		double positionFirstSignalPlane = planeOffset + 2. * planespacing;

		double baseUpperEndCapPos =
			(planeOffset + planespacing * 2. +
		       	(nrofPlanes - 1) * planespacing*3 + planespacing * 2. +
         		spacerThickness);
		double baseLowerEndCapPos = -(planeOffset + spacerThickness);

		  // face of the scintillator is 32 mm away from the aluminium endcaps
		  // gasBoxY = 2.0*(baseUpperEndCapPos + 33.0*mm);
		double gasBoxY = baseUpperEndCapPos - baseLowerEndCapPos + 2*innerEndCapThickness;
		offset = (baseUpperEndCapPos + baseLowerEndCapPos)/2.;

		gasboxhalfX = 12.;
		gasboxhalfY = gasBoxY*0.5;
		gasboxhalfZ = 12.;
		muonPlane = gasBoxY*0.5 + offset;
                firstWireXpos = 16*cellR- 6.47;
	}
	else{
		gasboxhalfX = 3.;
		gasboxhalfY = 3.;
		gasboxhalfZ = 12.;
		offset = 0.;
		muonPlane = gasboxhalfY;
	}
	box = new Garfield::SolidBox(0., offset, 0.,gasboxhalfX,gasboxhalfY,gasboxhalfZ);
	geo->AddSolid(box, fMediumMagboltz);

}

void SignalGenerator::BuildCompField(){
	comp = new Garfield::ComponentAnalyticField();
	comp->SetGeometry(geo);
	if(fullChamber){
		int cn = 0;
		double x,y = 0;
		for (int j = 0; j < nrofPlanes; j++) {
			for (int i = 0; i < nrofWiresPlane; i++) {
				std::stringstream cn_stream;
				cn_stream << cn;
				GetSignalWirePosition(i + j * nrofWiresPlane,x,y);
				comp->AddWire(x,y,2 * Swr, signalWireVoltages[j], cn_stream.str());
				comp->AddReadout(cn_stream.str());
				cn++;
				if(i==0){
  				x = x + 2*cellR;
					comp->AddWire(x,y,2 * Swr, signalWireVoltages[j], "s_extra");
				}
				if(i==nrofWiresPlane-1){
          x = x - 2*cellR;
          comp->AddWire(x,y,2 * Swr, signalWireVoltages[j], "s_extra");
        }
			}
		}
		for (int j = 0; j < nrofPlanes + 1; j++) {
			for (int i = 0; i < nrofWiresPlane + 3; i++) {
				double wirePositionX = -(-firstWireXpos + (i-1) * cellR*2 - (j % 2) * cellR);
				double wirePositionY = planeOffset + j * planespacing*3;

				comp->AddWire(wirePositionX, wirePositionY,2 * Fwr , 0., "f");

				wirePositionX=-(-firstWireXpos + (i-1) * cellR*2 - !(j % 2) * cellR);
				wirePositionY=planeOffset + j * planespacing*3 + planespacing;

				comp->AddWire(wirePositionX, wirePositionY,2 * Fwr , 0., "f");

			}
		}
	}
	else{
		comp->AddWire(0., 0., 2 * Swr, sigWireVoltage, "s");
		comp->AddReadout("s"); //This calculates a weighting field for the signal wires
		comp->AddWire(0., 1., 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(0., -1., 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(cellR, planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(cellR, -planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-cellR, planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-cellR, -planespacing, 2 * Fwr, fieldWireVoltage, "f");
		//cell on the right
		comp->AddWire(2*cellR, 0., 2 * Swr, sigWireVoltage, "f");
		comp->AddWire(2*cellR, 1., 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(2*cellR, -1., 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(3*cellR, planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(3*cellR, -planespacing, 2 * Fwr, fieldWireVoltage, "f");
		//cell on the left
		comp->AddWire(-2*cellR, 0., 2 * Swr, sigWireVoltage, "f");
		comp->AddWire(-2*cellR, 1., 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-2*cellR, -1., 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-3*cellR, planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-3*cellR, -planespacing, 2 * Fwr, fieldWireVoltage, "f");
		//cell upper right
		comp->AddWire(cellR, 3*planespacing, 2 * Swr, sigWireVoltageUpper, "f");
		comp->AddWire(cellR, 5*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(0., 4*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(2*cellR, 4*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		//cell upper left
		comp->AddWire(-cellR, 3*planespacing, 2 * Swr, sigWireVoltageUpper, "f");
		comp->AddWire(-cellR, 5*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-2*cellR, 4*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		//cell lower right
		comp->AddWire(cellR, -3*planespacing, 2 * Swr, sigWireVoltageLower, "f");
		comp->AddWire(cellR, -5*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(0., -4*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(2*cellR, -4*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		//cell lower left
		comp->AddWire(-cellR, -3*planespacing, 2 * Swr, sigWireVoltageLower, "f");
		comp->AddWire(-cellR, -5*planespacing, 2 * Fwr, fieldWireVoltage, "f");
		comp->AddWire(-2*cellR, -4*planespacing, 2 * Fwr, fieldWireVoltage, "f");
	}


}

void SignalGenerator::GetSignalWirePosition(int cellNo,double &x,double &y) {
  int column = cellNo % nrofWiresPlane;
  int planeNo = cellNo / nrofWiresPlane;
  x=-column * cellR * 2 + firstWireXpos - (planeNo % 2) * cellR;
  y=planeOffset + 2 * planespacing + planeNo * planespacing*3;
}

void SignalGenerator::BuildSensor(){
	fSensor = new Garfield::Sensor();

	if(transferFilename == "")
		fSensor->SetTransferFunction(transfer);
	else
    SetTransferFile(transferFilename);

	fSensor->AddComponent(comp);

	fSensor->AddElectrode(comp,"s");

	fSensor->SetTimeWindow(0.,fBinWidth,fNbins); //Lowest time [ns], time bins [ns], number of bins
}

void SignalGenerator::SetTracking(){
	if(driftRKF){
		fDriftRKF = new Garfield::DriftLineRKF();
		fDriftRKF->SetSensor(fSensor);
		fDriftRKF->EnableDebugging();
	}
	else{
		fAvalanche = new Garfield::AvalancheMicroscopic();
		fAvalanche->SetSensor(fSensor);
		fAvalanche->EnableSignalCalculation();
    if (avSize>0)
		  fAvalanche->EnableAvalancheSizeLimit(avSize);
		std::cout << fAvalanche->GetAvalancheSizeLimit() << std::endl;
		fDrift = new Garfield::AvalancheMC();
		fDrift->SetSensor(fSensor);
		fDrift->EnableSignalCalculation();
		fDrift->SetDistanceSteps(1.e-3);
	    // fDrift->EnableDebugging();
		if(createAval) fDrift->EnableAttachment();
		else fDrift->DisableAttachment();
	}
	fTrackHeed = new Garfield::TrackHeed();
  fTrackHeed->SetSensor(fSensor);
  fTrackHeed->SetParticle(particleType);
  if (debugging) fTrackHeed->EnableDebugging();
	fTrackHeed->EnableDeltaElectronTransport();

}

void SignalGenerator::CreateChamberView(){
	fChamber = new TCanvas("c", "Chamber View", 700, 700);
	cellView = new Garfield::ViewCell();
	cellView->SetComponent(comp);
	cellView->SetCanvas(fChamber);
	std::cout << "Canvas Set" << std::endl;
	cellView->Plot2d();
	fChamber->Update();
	fChamber->Print("chamber_configuration.pdf");
//	gSystem->ProcessEvents();
	std::cout << "CreateCellView()" << std::endl;

	viewDrift = new Garfield::ViewDrift();
	viewDrift->SetCanvas(fChamber);
	if(!driftRKF){
		fAvalanche->EnablePlotting(viewDrift);
		fDrift->EnablePlotting(viewDrift);
	}
	if(tracks) fTrackHeed->EnablePlotting(viewDrift);

}

void SignalGenerator::CreateSignalView(){
	fSignal = new TCanvas("c2", "Signal on the wire", 700, 700);
	viewSignal = new Garfield::ViewSignal();
	viewSignal->SetSensor(fSensor);
	viewSignal->SetCanvas(fSignal);

}

void SignalGenerator::CreateFieldView(){
	fField = new TCanvas("c3", "Weightingfield", 700, 700);
	viewField = new Garfield::ViewField();
	viewField->SetCanvas(fField);
	viewField->SetComponent(comp);
//	viewField->SetArea(-1.1,-1.1,1.1,1.1);
	viewField->SetNumberOfContours(40);
	// viewField->PlotContour("e");
  viewField->Plot("e","Colz");
//	viewField->SetWeightingFieldRange(0.,500);
//	viewField->SetNumberOfContours(40);
//	viewField->PlotContourWeightingField("20","e");
	fField->Update();
	fField->Print("Efield.pdf");

	if(cellNr>=0){
		double x,y;
    cellNr = 0;
		GetSignalWirePosition(cellNr,x,y);
		fField->Clear();
		viewField->PlotProfile(x,y,0.,x+1.,y,0.,"e");
		fField->Update();
		// fField->Print("Efield_R_0deg.pdf");
		// fField->Clear();
		viewField->PlotProfile(x,y,0.,x+1.*cos(30*CLHEP::pi/180),y+1.*sin(30*CLHEP::pi/180),0.,"e");
		fField->Update();
    cellNr = 27;
		GetSignalWirePosition(cellNr,x,y);
    viewField->PlotProfile(x,y,0.,x+1.,y,0.,"e");
		fField->Update();
		viewField->PlotProfile(x,y,0.,x+1.*cos(30*CLHEP::pi/180),y+1.*sin(30*CLHEP::pi/180),0.,"e");
		fField->Update();
		fField->Print("Efield_profiles.pdf");
		int rSize = 500;
		int tSize = 2;
		double rMax=1.;
		double tMax = 30.;
		double rStart = Swr+0.000001;
		double b = log(rMax/rStart)/(rMax-rStart);
		double a = rMax / exp(rMax*b) ;
		std::cout << a << " " << b << std::endl;
		double r,t;
		double ex,ey,ez;
		double vx=0.,vy=0.,vz=0.,dl=0.,dt=0.,alpha=0.;
		int status;
		Garfield::Medium* medium = nullptr;
		std::ofstream File;
		char buffer[20];
		sprintf(buffer,"field_cell_%d.txt",cellNr);
		File.open(buffer);
		std::ofstream FileField;
		FileField.open("Param_vs_Field_long_range.txt");
		for(int i=0;i<100;i++){
			double field = (double)(1200*i);
                	bool vel = fMediumMagboltz->ElectronVelocity(field,0.,0.,0,0,0,vx,vy,vz);
                        bool diff = fMediumMagboltz->ElectronDiffusion(field,0.,0.,0.,0.,0.,dl,dt);
                        bool townsend = fMediumMagboltz->ElectronTownsend(field,0.,0.,0.,0.,0.,alpha);
                        FileField << field << " " << vx << " " << dl << " " << dt << " " << alpha << std::endl;
                }
		FileField.close();
		FileField.open("Param_vs_Field.txt");
		for(int i=0;i<100;i++){
			double field = (double)(100*i);
                	bool vel = fMediumMagboltz->ElectronVelocity(field,0.,0.,0,0,0,vx,vy,vz);
                        bool diff = fMediumMagboltz->ElectronDiffusion(field,0.,0.,0.,0.,0.,dl,dt);
                        bool townsend = fMediumMagboltz->ElectronTownsend(field,0.,0.,0.,0.,0.,alpha);
                        FileField << field << " " << vx << " " << dl << " " << dt << " " << alpha << std::endl;
                }
		FileField.close();

		for(int i=0;i<rSize;i++){
			for(int j=0;j<tSize;j++){
				//r = a*exp(b*(rStart+(rMax-rStart)*i/(rSize-1)));
				r = 0.002+i*0.002;
				t = tMax*j/tSize;
				comp->ElectricField(x+r*cos(t*CLHEP::pi/180),y+r*sin(t*CLHEP::pi/180),0,ex,ey,ez,medium,status);
				bool vel = fMediumMagboltz->ElectronVelocity(ex,ey,ez,0,0,0,vx,vy,vz);
				bool diff = fMediumMagboltz->ElectronDiffusion(ex,ey,ez,0.,0.,0.,dl,dt);
				bool townsend = fMediumMagboltz->ElectronTownsend(ex,ey,ez,0.,0.,0.,alpha);
				File << r << " " << t << " " << ex << " " << ey << " " << vx << " "  << vy << " " << dl << " " << dt << " " << alpha << std::endl;
			}
		}
		File.close();
    int cell_nrs[2] = {0,27};
    for(int k=0;k<2;k++){
      GetSignalWirePosition(cell_nrs[k],x,y);
      for(int j=0;j<tSize;j++){
        	t = tMax*j/(tSize-1);
        char buffer2[20];
        sprintf(buffer2,"field_cell_%d_theta_%d.txt",cell_nrs[k],(int)t);
        std::ofstream field_profiles;
        field_profiles.open(buffer2);
        for(int i=0;i<rSize;i++){
          r = 0.002+i*0.002;
  				comp->ElectricField(x+r*cos(t*CLHEP::pi/180),y+r*sin(t*CLHEP::pi/180),0,ex,ey,ez,medium,status);
          field_profiles << r << " " << ex << " " << ey << std::endl;
        }
        field_profiles.close();
      }
    }
	}
}



void SignalGenerator::Run() {
	const char* myfilename;
	myfilename = inFileName.c_str();

	std::ifstream input(myfilename);
  	std::string line;
  	double x,y,z;
  	if(tracks){
  		std::cout << "Full tracking" << std::endl;
  		double dx, dy,dz, ekin=0.;
  		if(muonDist){
  			for(int j=0;j<nSim;j++){
  				muonGen->nextBeamEvent(x,y,z,dx,dy,dz,ekin);
  				anManager->SetInitParam(x,y,z,dx,dy,dz,ekin);
  				fTrackHeed->SetKineticEnergy(1e9);
  				float ne=0, ni=0;
  				if (j % 10 == 0) std::cout << "Event " << j << std::endl;
  				fSensor->ClearSignal();
	  			fTrackHeed->NewTrack(x,y,z,0,dx,dy,dz);
	  			double x_pos,y_pos,z_pos,t,energy,extra;
	  			int nc;
	  			int nInitEl=0;
	  			int nDetEl=0;
	  			while(fTrackHeed->GetCluster(x_pos, y_pos, z_pos, t, nc, energy, extra)) { //Loop over all created clusters
					while (nc) {
						nc--;
				  		double xe, ye, ze, te;
				  		double ee, dxe, dye, dze;
				  		fTrackHeed->GetElectron(nc, xe, ye, ze, te, ee, dxe, dye, dze);
				  		if(std::sqrt(xe*xe+ye*ye)>2) continue;
				  		if(driftRKF) {
	                		fDriftRKF->DriftElectron(xe,ye,ze,te);
	                		if(driftIons) fDriftRKF->DriftIon(xe,ye,ze,te);
			            }
						else if(trackMicro){
							fAvalanche->AvalancheElectron(xe, ye, ze, te, ee, dxe, dye, dze);
							int neTemp=0,niTemp=0;
							fAvalanche->GetAvalancheSize(neTemp, niTemp);
							ne+=(float)neTemp;
							ni+=(float)niTemp;
							if(driftIons){
								const int np = fAvalanche->GetNumberOfElectronEndpoints();
	    						double xe1, ye1, ze1, te1, e1;
	    						double xe2, ye2, ze2, te2, e2;
	    						int status;
							    for (int k = np; k--;) {
							      fAvalanche->GetElectronEndpoint(k, xe1, ye1, ze1, te1, e1,
							                                   xe2, ye2, ze2, te2, e2, status);
							      fDrift->DriftIon(xe1, ye1, ze1, te1);
							    }
							}
						}
						else {
							double xe1, ye1, ze1, te1;
    						double xe2, ye2, ze2, te2;
    						int status;
							fDrift->DriftElectron(xe,ye,ze,te);
							fDrift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);
							if(std::sqrt(xe*xe+ye*ye)<1){
								nInitEl++;
								if(WireIsHit(status,xe2,ye2,ze2)) nDetEl++;
							}
							unsigned int neTemp=0,niTemp=0;
							fDrift->GetAvalancheSize(neTemp, niTemp);
							ne+=(float)neTemp;
							ni+=(float)niTemp;
							if(driftIons) fDrift->DriftIon(xe,ye,ze,te);
						}
					}
		  		}
		  		ProcessEvent(ne,ni,nInitEl,nDetEl);
  			}
  		}
  		else{
	  		while(getline(input,line)){
	  			std::istringstream stream(line);
	  			stream >> x >> y >> z >> dx >> dy >> dz >> ekin;
	  			anManager->SetInitParam(x,y,z,dx,dy,dz,ekin);
	  			fTrackHeed->SetKineticEnergy(ekin*1000000.);
	  			for(int j=0;j<nSim;j++){
	  				float ne=0, ni=0;
	  				if (j % 10 == 0) std::cout << "Event " << j << std::endl;
	  				fSensor->ClearSignal();
		  			fTrackHeed->NewTrack(x,y,z,0,dx,dy,dz);
		  			double x_pos,y_pos,z_pos,t,energy,extra;
		  			int nc;
		  			int nInitEl=0;
		  			int nDetEl=0;
		  			while(fTrackHeed->GetCluster(x_pos, y_pos, z_pos, t, nc, energy, extra)) { //Loop over all created clusters
              while (nc) {
                nc--;
                double xe, ye, ze, te;
                double ee, dxe, dye, dze;
                bool wire_hit=false;
                fTrackHeed->GetElectron(nc, xe, ye, ze, te, ee, dxe, dye, dze);
                if(std::sqrt(xe*xe+ye*ye)>1.5) continue;
                if(driftRKF) {
                  fDriftRKF->DriftElectron(xe,ye,ze,te);
                  if(driftIons) fDriftRKF->DriftIon(xe,ye,ze,te);
                }
                else if(trackMicro){
                  fAvalanche->AvalancheElectron(xe, ye, ze, te, ee, dxe, dye, dze);
                  int neTemp=0,niTemp=0;
                  fAvalanche->GetAvalancheSize(neTemp, niTemp);
                  ne+=(float)neTemp;
                  ni+=(float)niTemp;
                  if(driftIons){
                    const int np = fAvalanche->GetNumberOfElectronEndpoints();
                    double xe1, ye1, ze1, te1, e1;
                    double xe2, ye2, ze2, te2, e2;
                    int status;
                    for (int j = np; j--;) {
                      fAvalanche->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
                                                xe2, ye2, ze2, te2, e2, status);
                      if(WireIsHit(status,xe2,ye2,ze2)) {
                        nDetEl++;
                      }
                      fDrift->DriftIon(xe1, ye1, ze1, te1);
                    }
                    nInitEl++;
                  }
                }
                else {
                  double xe1, ye1, ze1, te1;
                  double xe2, ye2, ze2, te2;
                  int status;
                  fDrift->DriftElectron(xe,ye,ze,te);
                  fDrift->GetElectronEndpoint(fDrift->GetNumberOfElectronEndpoints()-1, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);
                  if(WireIsHit(status,xe2,ye2,ze2)) {
                    nDetEl++;
                  }
                  nInitEl++;
                  unsigned int neTemp=0,niTemp=0;
                  fDrift->GetAvalancheSize(neTemp, niTemp);
                  ne+=(float)neTemp;
                  ni+=(float)niTemp;
                  if(driftIons) fDrift->DriftIon(xe,ye,ze,te);
                }
              }
            }
			  		ProcessEvent(ne,ni,nInitEl,nDetEl);
			  	}
	  		}
	  	}
  	}
  	else {
  		while(getline(input,line)){
  			std::cout << "Loop over input file lines" << std::endl;
  			std::istringstream stream(line);
  			stream >> x >> y >> z;
  			anManager->SetInitParam(x,y,z,0,0,0,0);
  			for(int j=0;j<nSim;j++){
  				float ne=0, ni=0;
  				if (j % 10 == 0) std::cout << "Event " << j << std::endl;
  				fSensor->ClearSignal();
		  		if(driftRKF) {
	        		fDriftRKF->DriftElectron(x,y,0,0);
                                double neTemp=0,niTemp=0;
                                fDriftRKF->GetAvalancheSize(neTemp, niTemp);
                                ne+=(float)neTemp;
                                ni+=(float)niTemp;
	        		if(driftIons) fDriftRKF->DriftIon(x,y,0,0);
	            }
				else if(trackMicro){
					fAvalanche->AvalancheElectron(x, y, 0, 0, 0, 0, 0, 0);
        				int neTemp=0,niTemp=0;
					fAvalanche->GetAvalancheSize(neTemp, niTemp);
                                        ne+=(float)neTemp;
 					ni+=(float)niTemp;
					if(driftIons){
						const int np = fAvalanche->GetNumberOfElectronEndpoints();
						double xe1, ye1, ze1, te1, e1;
						double xe2, ye2, ze2, te2, e2;
						double xi1, yi1, zi1, ti1;
							double xi2, yi2, zi2, ti2;
						int status;
					    for (int j = np; j--;) {
					      fAvalanche->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
					                                   xe2, ye2, ze2, te2, e2, status);
					      fDrift->DriftIon(xe1, ye1, ze1, te1);
					    }
					}
				}
				else {
					std::cout << "Drift electrons" << std::endl;
					fDrift->DriftElectron(x,y,0,0);
					unsigned int neTemp=0,niTemp=0;
					fDrift->GetAvalancheSize(neTemp, niTemp);
					ne+=(float)neTemp;
					ni+=(float)niTemp;
					if(driftIons) fDrift->DriftIon(x,y,0,0);
				}
				ProcessEvent(ne,ni,1,0);
			}
  		}
  	}
}


void SignalGenerator::ProcessEvent(float ne, float ni, int nInitEl, int nDetEl){
	anManager->ProcessEvent(fSensor,ne,ni,nInitEl,nDetEl,signalThreshold);
	if(visualizeChamber){
		viewDrift->Plot(true,false);
		fChamber->Update();
		fChamber->Print("chamberview.pdf");
	}
	if(visualizeSignal){
		viewSignal->PlotSignal("s");
		fSignal->Update();
		fSignal->Print("signal.pdf");
	}
}




void SignalGenerator::SetTransferFile(std::string name){
	std::cout << "SignalGenerator::SetTransferFile" << std::endl;
	std::string transferData = OptionContainer::GetInstance().GetOption<std::string>("transferdata") + name;
	std::cout << transferData << std::endl;
	const char * filename;
	filename = transferData.c_str();
	std::cout << filename << std::endl;
	std::ifstream input(filename);
  std::string line;
  std::string item;
  double value;
  bool first=true;
  bool firstline=true;
  double firstvalue;
  while(getline(input,line)){
      std::istringstream streamline(line);
      first=true;
      while(getline(streamline,item,',')){
          std::istringstream streamitem(item);
          streamitem >> value;
          if(firstline){
              firstvalue=value;
              firstline=false;
          }
          if(first) {
              transferTimeVector.push_back((value-firstvalue)*1.e9-200);
//              std::cout << transferTimeVector[transferTimeVector.size()-1] << " " ;
              first=false;
          }
          else{
              transferValuesVector.push_back(value);
//               std::cout << value << std::endl;
          }
      }
  }
  fSensor->SetTransferFunction(transferTimeVector,transferValuesVector);
}


double transfer(double t){
	const double k = -100.;
	const double tau1 = 330.0;
	const double tau2 = 100.0;
	double Vout = k*tau1/(tau1-tau2) * (exp(-t/tau1) - exp(-t/tau2));

/*	double Vout = (-1.794*std::pow(10,4)*std::pow(t,3) + 3.736*std::pow(10,9)*std::pow(t,2) + 1.261*std::pow(10,15)*t - 3.736*std::pow(10,12)) /
					(5.155*std::pow(t,3) + 3.15*std::pow(10,5)*std::pow(t,2) + 1.4*std::pow(10,9)*t - 3.108*std::pow(10,7));
*/
	//double Vout = (/*1.486*std::pow(10,2)*std::pow(t,3) + */4.757*std::pow(10,7)*std::pow(t,2) + 1.936*std::pow(10,11)*t + 1.985*std::pow(10,9)) /
		//			(/*2.756*std::pow(t,3)*/ + 9.782*std::pow(10,5)*std::pow(t,2) + 6.267*std::pow(10,10)*t + 2.414*std::pow(10,14));
	return Vout;
}

bool SignalGenerator::WireIsHit(int status,double x,double y, double z){
  // std::cout << "Status: " << status << ", (x,y,z): (" << x << "," << y << "," << z << ")" << std::endl;
	if(status != Garfield::StatusLeftDriftMedium) return false;
	if(std::sqrt(x*x+y*y)<2*Swr) return true;
	return false;
}
