int ReadData(char* FileName, float x[], float xscale, float y[], float yscale )
{
 
 ifstream in;
 in.open(FileName);
 int nlines = 0; 

while (1) {
      in >> x[nlines] >> y[nlines];
      if (!in.good()) break;
      x[nlines]=xscale*x[nlines];
      y[nlines]=yscale*y[nlines];
      nlines++;
   }
 cout<<"*** In total "<<nlines<<" data points taken from "<<FileName<<" and rescaled"<<endl;
 cout<<"    First point: x="<<x[0]<<"  y="<<y[0]<<endl;
 cout<<"    Last point:  x="<<x[nlines-1]<<"  y="<<y[nlines-1]<<endl;
 return nlines;
}

void pseudorapidity(){
gROOT->Reset();
gROOT->ProcessLine("#include <vector>");
gROOT->LoadMacro("ReadData.cc");
  bool graypalette = 0;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(62);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
// Make a canva and a pad
TCanvas* depo = new TCanvas("depo","",200,10,1000,800);
TPad* pad = new TPad("pad","The pad with the historgram",0.02,0.02,0.98,0.98,21);
pad->SetBorderMode(0); pad->SetFillColor(kWhite); pad->SetTickx(); pad->SetTicky(); pad->Draw();
pad->cd();
  
TH2D* etaHistoBack = new TH2D("","",1000,5.5,20,1000,0,20);
  etaHistoBack->SetXTitle("#eta");
  etaHistoBack->SetYTitle("#frac{1}{N_{ev.}} #frac{dn_{Z}}{d#eta}");
  etaHistoBack->GetXaxis()->SetTitleOffset(0.8);
  etaHistoBack->GetYaxis()->SetTitleOffset(1.1);
  etaHistoBack->GetXaxis()->SetNdivisions(210);
  etaHistoBack->GetXaxis()->SetLabelOffset(0.007);
  etaHistoBack->SetStats(kFALSE);
etaHistoBack->Draw();

  TH1D* etaHisto = new TH1D("momentum_histo","", 40+7*5,5,20);
  etaHisto->SetLineColor(kBlack);
  etaHisto->SetLineWidth(3.5);
  etaHisto->SetLineStyle(kSolid);
  TH1D* etaHistoHe = new TH1D("momentum_histo","", 40+7*5 ,5,20);
  etaHistoHe->SetLineColor(kBlue);
  etaHistoHe->SetLineWidth(3.5);
  etaHistoHe->SetLineStyle(kDashed);
  TH1D* etaHistoFr = new TH1D("momentum_histo","", 40+7*5 ,5,20);
  etaHistoFr->SetLineColor(kRed+1);
  etaHistoFr->SetLineWidth(3.5);
  etaHistoFr->SetLineStyle(9);

std::vector<double>* px = 0;
std::vector<double>* py = 0;
std::vector<double>* pz = 0;
std::vector<double>* Z  = 0;
std::vector<double>* A  = 0;

TFile* ReadFile = new TFile("../dcmqgsm_in_aamcc_format_sps.root");
TTree *tree = (TTree*) ReadFile->Get("Glauber");
tree->SetBranchAddress("Z_on_A", &Z);
tree->SetBranchAddress("A_on_A", &A);
tree->SetBranchAddress("pX_on_A", &px);
tree->SetBranchAddress("pY_on_A", &py);
tree->SetBranchAddress("pZ_on_A", &pz);

double norm = tree->GetEntries();


int nlines;
float xVal[22];
float yVal[22];
float xErr[22];
float yErr[22];



nlines = ReadData("cherry1998_single_charge.dat", xVal, 1., yVal, 1.);
TGraphErrors* cherry1998 = new TGraphErrors(22, xVal, yVal, xErr, yErr);
nlines = ReadData("cherry1998_double_charge.dat", xVal, 1., yVal, 1.);
TGraphErrors* cherry1998He = new TGraphErrors(22, xVal, yVal, xErr, yErr);
nlines = ReadData("cherry1998_frags.dat", xVal, 1., yVal, 1.);
TGraphErrors* cherry1998Fr = new TGraphErrors(22, xVal, yVal, xErr, yErr);

  cherry1998->SetMarkerStyle(8);
  cherry1998->SetMarkerColor(kBlack);
  cherry1998->SetLineColor(kBlue);
  cherry1998->SetLineStyle(kSolid);
  cherry1998->SetMarkerSize(1.5);
  cherry1998->SetLineWidth(0);
    
  cherry1998He->SetMarkerStyle(8);
  cherry1998He->SetMarkerColor(kBlue);
  cherry1998He->SetLineColor(kBlue);
  cherry1998He->SetLineStyle(kSolid);
  cherry1998He->SetMarkerSize(1.5);
  cherry1998He->SetLineWidth(0);
  
  cherry1998Fr->SetMarkerStyle(8);
  cherry1998Fr->SetMarkerColor(kRed+1);
  cherry1998Fr->SetLineColor(kBlue);
  cherry1998Fr->SetLineStyle(kSolid);
  cherry1998Fr->SetMarkerSize(1.5);
  cherry1998Fr->SetLineWidth(0);
  

	for(int k = 1; k < tree->GetEntries(); k++){
		tree->GetEntry(k); 
		
		double pX = 0;
		double pY = 0;
		double pZ = 0;
		double p = 0;
		double tan_theta = 0;
		double cos_theta = 0;
		double theta = 0;
		double eta = 0;			           

		for(int k = 0; k< px->size(); k++ ){pX=px->at(k);pY=py->at(k);pZ=pz->at(k);
		                                   tan_theta = std::sqrt(pX*pX + pY*pY)/pZ;
						   theta = atan(tan_theta);
						   p = std::sqrt(pX*pX + pY*pY + pZ*pZ);
						   cos_theta = pZ/p;
						   //eta = 0.5*log((p+pZ)/(p-pZ));
						   eta = -log(std::sqrt((1 - cos_theta)/(1 + cos_theta)));
              // eta = atan(pZ/p);
						   if(Z->at(k) == 0) etaHisto->Fill(eta);
						   if(Z->at(k) == 2) etaHistoHe->Fill(eta);
						   if(Z->at(k) > 2) etaHistoFr->Fill(eta);
						   }
						   

	}


TLegend* leg = new TLegend(0.336,0.561,0.812,0.867,"Pb-Pb 158A GeV");
leg->SetTextSize(0.05);
leg->SetFillColor(0);
leg->AddEntry(etaHisto,"DCM-QGSM, Z = 1", "lp");
leg->AddEntry(etaHistoHe,"DCM-QGSM, Z = 2", "lp");
leg->AddEntry(etaHistoFr,"DCM-QGSM, Z #geq 3", "lp");
leg->AddEntry(cherry1998,"Cherry et al. (1998), Z = 1","lp");
leg->AddEntry(cherry1998He,"Cherry et al. (1998), Z = 2","lp");
leg->AddEntry(cherry1998Fr,"Cherry et al. (1998), Z #geq 3","lp");
leg->Draw();

cherry1998->Draw("P");
cherry1998He->Draw("P SAME");
cherry1998Fr->Draw("P SAME");

etaHisto->Scale(1./(0.1*norm));
etaHisto->Draw("HIST SAME");

etaHistoHe->Scale(1./(0.1*norm));
etaHistoHe->Draw("HIST SAME");

etaHistoFr->Scale(1./(0.1*norm));
etaHistoFr->Draw("HIST SAME");

pad->Update();

delete ReadFile;
}
