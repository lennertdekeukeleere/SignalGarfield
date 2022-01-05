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
	tree->Branch("signalHeightElectron",&signalHElectron);
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
	double sig_h, sig_h1, sig_h2;
	double thresh_level = GetADC(totalSignalHist,true,1.,adcTotal,signalH,signalHat100ns,signalHat200ns);
	if(threshold>0) thresh_level=threshold;
	tdcTotal = GetTDC(totalSignalHist,1.,thresh_level);
	thresh_level = GetADC(electronSignalHist,true,-1.,adcElectron,signalHElectron,sig_h1,sig_h2);
	tdcElectron = GetTDC(electronSignalHist,-1.,thresh_level);
	thresh_level = GetADC(ionSignalHist,true,-1.,adcIon,sig_h,sig_h1,sig_h2);
	tdcIon = GetTDC(ionSignalHist,-1.,thresh_level);
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

double AnalysisManager::GetTDC(TH1D* hist, double factor, double thresh_level){
	// std::cout << "Threshold: " << thresh_level << std::endl;
	for(int i = 1; i<=hist->GetSize()-2; i++){
		std::cout << hist->GetBinContent(i)*factor << std::endl;
		if(hist->GetBinContent(i)*factor>thresh_level)
			return hist->GetBinLowEdge(i);
	}
	return -1.;
}

double AnalysisManager::GetADC(TH1D* hist, bool calc_height, double factor, double& adc, double& sig_height,double& sig_height_1, double& sig_height_2){
	double vMin = hist->GetBinContent(1)*factor;
	double vMax = hist->GetBinContent(1)*factor;
	double sig;
	int nbins = hist->GetSize()-2;

	for (int i = nbins; i>0; i--) {
		sig = hist->GetBinContent(i)*factor;
		if (sig > vMax) vMax = sig;
	}
	double thresh_level = vMin +(vMax-vMin)*0.2;

	adc = 0.;
	for(int i=1; i<nbins; i++){
		adc+= hist->GetBinContent(i)*factor*hist->GetBinWidth(i);
	}

	sig_height = -1;
	sig_height_1 = -1;
	sig_height_2 = -1;

	if (calc_height){
		for (int j=0;j<400;j++){
			sig = hist->GetBinContent(j)*factor;
			if(sig > sig_height_1) sig_height_1 = sig;
		}
		sig_height_2 = sig_height_1;
		for (int j=400; j<800;j++){
			sig = hist->GetBinContent(j)*factor;
			if (sig > sig_height_2) sig_height_2 = sig;
		}
		sig_height = vMax;
	}

	return thresh_level;
}
