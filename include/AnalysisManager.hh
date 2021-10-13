#ifndef ANALYSISMANAGER_H
#define ANALYSISMANAGER_H

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF2.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include <vector>

//Basic stuff
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>

#include "Garfield/Sensor.hh"


class AnalysisManager{
	public:
	AnalysisManager(std::string);
	~AnalysisManager(){};
	void ProcessEvent(Garfield::Sensor*, float ne, float ni, int nInitEl,int nDetEl, double threshold);
	void Initialise(Garfield::Sensor*, bool full);
        void SetInitParam(double,double,double,double,double);
        void Write();

	private:
	double GetTDC(TH1D*,double threshold,bool);
    	double GetADC(TH1D*);

        std::string outFileName;

    	TH1D* totalSignalHist;
    	TH1D* electronSignalHist;
    	TH1D* ionSignalHist;
        TFile* f;
        TTree* tree;
    	double tdcTotal,adcTotal;
        double tdcElectron,adcElectron;
        double tdcIon,adcIon;
        float nEl, nIon;
	int nInitE,nDetE;
        double x,y,dx,dy,ekin;
	double signalH;
	std::ofstream outfile;	

};


#endif
