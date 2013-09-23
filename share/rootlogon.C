{
  G__loadfile("vector");
  G__loadfile("sstream");
  G__loadfile("bitset");
  G__loadfile("iomanip");
  
  std::cout << " Loading library \"Avec.so\"...";
  gSystem -> Load ( "Avec.so" );
  std::cout << " done" << endl;

  std::cout << " Loading library \"LAS_analysis.so\"...";
  gSystem -> Load ( "LAS_analysis.so" );
  std::cout << " done" << endl;

//   std::cout << "Loading LAS_vectorfloat_tools.C ...";
//   gROOT->ProcessLine(".L LAS_vectorfloat_tools.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_basic_tools.C ...";
//   gROOT->ProcessLine(".L LAS_basic_tools.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_globaldata_tools.C ...";
//   gROOT->ProcessLine(".L LAS_globaldata_tools.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_alpar.C ...";
//   gROOT->ProcessLine(".L LAS_alpar.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_Tec_Reconstruction.C ...";
//   gROOT->ProcessLine(".L LAS_Tec_Reconstruction.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_control_plots.C ...";
//   gROOT->ProcessLine(".L LAS_control_plots.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_RDC_tools.C ...";
//   gROOT->ProcessLine(".L LAS_RDC_tools.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_histo.C ...";
//   gROOT->ProcessLine(".L LAS_histo.C+");
//   std::cout << " done" << endl;

//   std::cout << "Loading LAS_stabil.C ...";
//   gROOT->ProcessLine(".L LAS_stabil.C+");
//   std::cout << " done" << endl;

  cout << "rootlogon loaded" << endl;
}
