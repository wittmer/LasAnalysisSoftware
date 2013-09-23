{
  cout << "******************************" << endl;
  cout << "*  Welcome to ROOT v" << gROOT->GetVersion() << "  *" << endl;
  cout << "******************************" << endl;
  cout << endl;

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

  std::cout << "Loading LAS_data_processing.C ...";
  gROOT->ProcessLine(".L LAS_data_processing.C+");
  std::cout << " done" << endl;

  gStyle -> SetCanvasBorderMode ( 0 );
  gStyle -> SetPadBorderMode ( 0 );
  gStyle -> SetCanvasColor ( 0 );
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetLabelSize(0.05);
  gStyle -> SetOptFit ( 1011 );
  
  cout << "rootlogon loaded" << endl;
}
