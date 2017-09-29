{
  //    gSystem->Load("libFWCoreFWLite.so");
  //    gSystem->Load("libDataFormatsFWLite.so"); 
    //   AutoLibraryLoader::enable();

 //  std::cout << std::endl << "Welcome to my rootlogon.C" << std::endl;
 // std::cout << "reading ~/setTDRStyle.C" << std::endl;
 // std::cout << "and some personal modifications." << std::endl << std::endl;

  //gSystem->Load("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/CfgManager/lib/CfgManagerDict.so");
  //gSystem->Load("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/lib/H4Dict.so");
  gROOT->ProcessLine(".L ./TDRStyle.C");
  setTDRStyle();
  tdrStyle->SetOptTitle(1);
  tdrStyle->SetPadTopMargin(0.06);
  tdrStyle->SetPadBottomMargin(0.14);
  //tdrStyle->SetPadLeftMargin(0.13);
  //tdrStyle->SetPadRightMargin(0.10); // per la paletta !
  tdrStyle->SetPadLeftMargin(0.10);
  tdrStyle->SetPadRightMargin(0.13); // per la paletta !
  tdrStyle->SetTitleXOffset(1.0);
  //tdrStyle->SetTitleYOffset(2.2);
  tdrStyle->SetTitleYOffset(1.0);
  tdrStyle->SetNdivisions(505, "X");
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPalette(1,0);

///////// pretty palette ///////////

   const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gROOT->ProcessLine("tdrStyle->SetNumberContours(NCont)");
  
/////////////////////////////////////

  gROOT->ForceStyle();
  //gROOT->SetStyle("tdrStyle");
}
