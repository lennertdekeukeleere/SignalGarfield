#include "AnalysisManager.hh"
#include "Garfield/Sensor.hh"

AnalysisManager::AnalysisManager(std::string out) : outFileName(out){

}

void AnalysisManager::Initialise(Garfield::Sensor* sens, bool full){
	double tStart, tStep;
	unsigned int nBins;
	sens->GetTimeWindow(tStart,tStep,nBins);
	std::string rootName = outFileName + ".root";
	f = new TFile(rootName.c_str(),"recreate", "Results", 9);
	totalSignalHist = new TH1D("tSig","Total Signal",nBins,tStart,nBins*tStep+tStart);
	electronSignalHist = new TH1D("eSig","Electron Signal",nBins,tStart,nBins*tStep+tStart);
	ionSignalHist = new TH1D("iSig","Ion Signal",nBins,tStart,nBins*tStep+tStart);

	tree = new TTree("Data","Signal data");
	// if(full){
	// 	// tree->Branch("totalHist","TH1D",&totalSignalHist,32000,0);
	// 	// tree->Branch("electronHist","TH1D",&electronSignalHist,32000,0);
	// 	// tree->Branch("ionHist","TH1D",&ionSignalHist,32000,0);
	//
	// }
	tree->Branch("tdcTotal",&tdcTotal);
	tree->Branch("adcTotal",&adcTotal);
	tree->Branch("signalHeight",&signalH);
	tree->Branch("signalHeight@100ns",&signalHat100ns);
	tree->Branch("signalHeight@200ns",&signalHat200ns);
	tree->Branch("tdcElectron",&tdcElectron);
	tree->Branch("adcElectron",&adcElectron);
	tree->Branch("tdcIon",&tdcIon);
	tree->Branch("adcIon",&adcIon);
	tree->Branch("NumOfEl",&nEl);
	tree->Branch("NumOfIon",&nIon);
	tree->Branch("NumOfInitEl",&nInitE);
	tree->Branch("NumOfDetEl",&nDetE);
	tree->Branch("Xinit",&x);
	tree->Branch("Yinit",&y);
	tree->Branch("Zinit",&z);
	tree->Branch("PXinit",&dx);
	tree->Branch("PYinit",&dy);
	tree->Branch("PZinit",&dz);
	tree->Branch("EKinit",&ekin);

	//char buffer[20];
	//sprintf(buffer,"%s.dat",outFileName);
	std::string datName = outFileName + ".dat";
	outfile.open(datName.c_str());



}

void AnalysisManager::SetInitParam(double x0,double y0, double z0, double dx0, double dy0, double dz0, double e0){
	x=x0;
	y=y0;
	z=z0;
	dx=dx0;
	dy=dy0;
	dz=dz0;
	ekin=e0;
	outfile << x << " " << y;
}

void AnalysisManager::ProcessEvent(Garfield::Sensor* sens,float ne, float ni, int nInitEl,int nDetEl, double threshold){
	outfile << " " << ne << std::endl;
	sens->ConvoluteSignal();
	double tStart, tStep;
	unsigned int nBins;
	sens->GetTimeWindow(tStart,tStep,nBins);

	totalSignalHist->Reset("ICESM");
	electronSignalHist->Reset("ICESM");
	ionSignalHist->Reset("ICESM");

	for(int i=0;i<nBins;i++){
		totalSignalHist->Fill(i*tStep,sens->GetSignal("s",i));
		electronSignalHist->Fill(i*tStep,sens->GetElectronSignal("s",i));
		ionSignalHist->Fill(i*tStep,sens->GetIonSignal("s",i));
	}
	tdcTotal = GetTDC(totalSignalHist,threshold,true);
	adcTotal = GetADC(totalSignalHist);
	tdcElectron = GetTDC(electronSignalHist,-1,false);
	adcElectron = GetADC(electronSignalHist);
	tdcIon = GetTDC(ionSignalHist,-1,false);
	adcIon = GetADC(ionSignalHist);
	nEl=ne;
	nIon=ni;
	nInitE=nInitEl;
	nDetE=nDetEl;
	tree->Fill();

}

void AnalysisManager::Write(){
	f->Write();
	outfile.close();
}

double AnalysisManager::GetTDC(TH1D* hist, double threshold, bool calcHeight){
	int thresh_cross = 0;
	double thresh_time = 0., thresh_level = 0.;
	bool thresh_rise = true;

	// Establish the range.
	double vMin = hist->GetBinContent(1);
	double vMax = hist->GetBinContent(1);
	double sig;
	int nbins = hist->GetSize()-2;

	for (int i = nbins; i>0; i--) {
		sig = hist->GetBinContent(i);
		if (sig > vMax) vMax = sig;
	}
	if(calcHeight){
		signalH=vMax;
	if(threshold<0)
		thresh_level = vMin +(vMax-vMin)*0.2;
	else
		thresh_level = threshold;

	double tdc = -1.0;
	for(int i = 1; i<= nbins; i++){
		if(hist->GetBinContent(i)>thresh_level) {
			tdc = hist->GetBinLowEdge(i);
			if (calcHeight){
				signalHat100ns = hist->GetBinContent(i);
				for (int j=i;j<i+100;j++){
					sig = hist->GetBinContent(j);
					if(sig > signalHat100ns) signalHat100ns = sig;
				}
				signalHat200ns = signalHat100ns;
				for (int j=i+100; j<i+200;j++){
					sig = hist->GetBinContent(j);
					if (sig > signalHat200ns) signalHat200ns = sig;
				}
			}
			break;
	}

	return tdc;
}

double AnalysisManager::GetADC(TH1D* hist){
	double adc = 0.;
	int nbins = hist->GetSize()-2;
	for(int i=1; i<nbins; i++){
//		G4cout << "Wire: " << s << "Signalheight at " << i << " : " <<  fSensor->GetSignal(s,i)*fBinWidth << G4endl;
		adc+= hist->GetBinContent(i)*hist->GetBinWidth(i);
	}
	return adc;
}
