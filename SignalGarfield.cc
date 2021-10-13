#include "OptionContainer.hh"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>

#include "SignalGenerator.hh"
#include "AnalysisManager.hh"

namespace po = boost::program_options;
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::ofstream;
using std::ifstream;

void ShowIntro() {
  double version = 0.1;


  cout << "------------------------------" << endl;
  cout << "-  SignalGarfield version " << version << "      -" << endl;
  cout << "-           June 27, 2018    -" << endl;
  cout << "------------------------------\n" << endl;
}

int main(int argc, char* argv[]) {
    ShowIntro(); 
    po::options_description cmdOptions("Commandline options");

  const char * filename;
  std::string outfilename, infilename, configName;
  std::string COMSOL_file; //This is the only COMSOL file that exists for now...

  cmdOptions.add_options()
    ("help,h", "Produce help message")
    ("output,o", po::value<std::string>(&outfilename)->default_value("output"), "Output file")
    ("input,I", po::value<std::string>(&infilename)->default_value("input.txt"), "Input file")
    ("config,c", po::value<std::string>(&configName)->default_value("config.txt"), "Config file")
    ("field,f", po::value<std::string>(&COMSOL_file)->default_value("1500V_1cell_wMed.garf"), "Comsol File");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOptions), vm);
  po::notify(vm);

  OptionContainer::GetInstance(configName);

  if(vm.count("help")) {
    cout << cmdOptions;
    cout << OptionContainer::GetConfigOptions();
    cout << OptionContainer::GetEnvOptions();
  }
  ///////////////
  //
  // This is where the main action happens...
  //
  else {

  	time_t start=time(0);
  	AnalysisManager* an = new AnalysisManager(outfilename);
  	SignalGenerator* gen = new SignalGenerator(an, infilename);
  	gen->Initialise();
  	gen->Run();
  	an->Write();
  	double duration = difftime(time(0),start);
	cout << "Simulation time: " << duration << endl;
    return 0;
  }
    
}
