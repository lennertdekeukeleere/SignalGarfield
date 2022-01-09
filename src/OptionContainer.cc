#include "OptionContainer.hh"

po::options_description OptionContainer::configOptions("\nConfiguration options");
po::options_description OptionContainer::envOptions("\nEnvironment variables");
po::variables_map OptionContainer::vm;

OptionContainer::OptionContainer(std::string configName) {
  configOptions.add_options()
    ("Signal.Resolution", po::value<double>()->default_value(1.), "Resolution of signal window (ns/bin)")
    ("Signal.Length", po::value<double>()->default_value(2000.), "Length of signal window (ns)")
    ("Signal.Threshold", po::value<double>()->default_value(-1.), "Use cfd(<0) or constant fraction discrimination with trigger threshold")
    ("Signal.Full", po::value<bool>()->default_value(false), "Output the full raw signal (true) or just ADC and TDC (false)")
    ("Signal.View", po::value<bool>()->default_value(false), "view signal (true) or not (false)")
    ("Signal.Transfer", po::value<std::string>()->default_value(""), "Path of the transferfunction file. If not given, analytical calculation is done.")
    ("Mixture.Iso", po::value<int>()->default_value(50), "Percentage isobutane (out of 100)")
    ("Mixture.Pres", po::value<int>()->default_value(600), "Pressure (in mbar)")
    ("Tracking.FullChamber", po::value<bool>()->default_value(false), "Full chamber configuration (true) or only adjacent cells (p9w5)(false)")
    ("Tracking.Track", po::value<bool>()->default_value(false), "Full tracking (true) or drift only (false)")
    ("Tracking.Avalanche", po::value<bool>()->default_value(false), "enable avalanche creation(true) or not(false)")
    ("Tracking.DriftInt", po::value<bool>()->default_value(false), "Integration of drift lines (true) or via Monte Carlo (false) / microscopic")
    ("Tracking.Microscopic", po::value<bool>()->default_value(false), "Microscopic tracking or via Monte Carlo (false)")
    ("Tracking.DriftIons", po::value<bool>()->default_value(false), "Drift only electrons (false) of drift ions too (true)")
    ("Tracking.View", po::value<bool>()->default_value(false), "View drift lines (true) or not (false)")
    ("Tracking.AvalancheSize", po::value<int>()->default_value(1000), "Set limit on the avalanche size")
    ("Tracking.Debugging", po::value<bool>()->default_value(false), "Enable debugging (true) or not (false)")
    ("Field.COMSOL", po::value<bool>()->default_value(false), "Field from COMSOL (true) or Garfield (false)")
    ("Field.View", po::value<bool>()->default_value(false), "Draw a fieldview (true) or not (false)")
    ("Field.PlotProfileCell", po::value<int>()->default_value(0), "Draw the profile of the electric field for cell with number x")
    ("Field.SigWireV", po::value<double>()->default_value(1500.0), "Magnitude of signal wire  potential (V) (plane 9)")
    ("Field.SigWireV1", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 1 potential (V) (plane 1)")
    ("Field.SigWireV2", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 2 potential (V) (plane 2)")
    ("Field.SigWireV3", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 3 potential (V) (plane 3)")
    ("Field.SigWireV4", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 4 potential (V) (plane 4)")
    ("Field.SigWireV5", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 5 potential (V) (plane 5)")
    ("Field.SigWireV6", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 6 potential (V) (plane 6)")
    ("Field.SigWireV7", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 7 potential (V) (plane 7)")
    ("Field.SigWireV8", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 8 potential (V) (plane 8)")
    ("Field.SigWireV9", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 9 potential (V) (plane 9)")
    ("Field.SigWireV10", po::value<double>()->default_value(1500.0), "Magnitude of signal wire 10 potential (V) (plane 10)")
    ("Field.SigWireVUpper", po::value<double>()->default_value(1500.0), "Magnitude of signal wire potential (V) of the upper plane (plane 10)")
    ("Field.SigWireVLower", po::value<double>()->default_value(1500.0), "Magnitude of signal wire potential (V) of the lower plane (plane 8)")
    ("Field.FieldWireV", po::value<double>()->default_value(0.0), "Magnitude of field wire potential (V)")
    ("Particle.Energy", po::value<double>()->default_value(1.0), "Beta particle energy (MeV)")
    ("Particle.CustomEnergy", po::value<bool>()->default_value(false), "read initial values from file or use random generator(only for muons at the moment)")
    ("Particle.Type", po::value<std::string>()->default_value("e-"), "Beta particle type (e+ or e-) or muons (mu)")
    ("Particle.Number", po::value<int>()->default_value(1), "Number of simulations")
    ("Chamber.NumberOfPlanes", po::value<int>()->default_value(8), "Number of planes")
 ;

  std::ifstream configStream(configName.c_str());
  if(!configStream.is_open()) {
    std::cerr << "ERROR: Config.txt cannot be found. Aborting.\n\n" << std::endl;
  }

  po::store(po::parse_config_file(configStream, configOptions), vm);
  po::notify(vm);

  envOptions.add_options()
    ("fielddata", po::value<std::string>()->default_value("/home/svn/svs_project042/MiniGarf/trunk/Fields/"), "Folder containing the COMSOL field data")
    ("gasdata", po::value<std::string>()->default_value("/home/svn/svs_project042/MiniGarf/trunk/GasesNew/"), "Folder of the mixture (gas, pressure) data")
    ("transferdata", po::value<std::string>()->default_value("/home/svn/svs_project042/MiniGarf/trunk/"), "Folder of the transfer function data")
    ("rootdata", po::value<std::string>()->default_value("/home/svn/svs_project042/MiniGarf/trunk/"), "Folder of the transfer function data")
  ;

  po::store(po::parse_environment(envOptions, "MiniGarf_"), vm);
  po::notify(vm);
}
