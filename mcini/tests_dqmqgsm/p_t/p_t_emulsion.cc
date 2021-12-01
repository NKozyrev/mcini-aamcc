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

void p_t_emulsion(){

gROOT->Reset();
gROOT->ProcessLine("#include <vector>");
//gROOT->LoadMacro("ReadData.cc");
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
  
TH2D* ptHistoBack = new TH2D("","",1000,0,350,1000,0,0.42);
  ptHistoBack->SetXTitle("p_{T}, MeV/#it{c}");
  ptHistoBack->SetYTitle("1/N_{Z = 2} #times dn_{Z = 2}/dp_{T}");
  ptHistoBack->GetXaxis()->SetTitleOffset(1.);
  ptHistoBack->GetYaxis()->SetTitleOffset(1.1);
  ptHistoBack->GetXaxis()->SetNdivisions(210);
  ptHistoBack->GetXaxis()->SetLabelOffset(0.007);
  ptHistoBack->SetStats(kFALSE);
ptHistoBack->Draw();

  TH1D* ptHistoAg = new TH1D("momentum_histo","", 400 ,0,10000);
  TH1D* ptHistoBr = new TH1D("momentum_histo","", 400 ,0,10000);

double sigma_Ag;
std::vector<double>* px_Ag = 0;
std::vector<double>* py_Ag = 0;
std::vector<double>* pz_Ag = 0;
std::vector<double>* Z_Ag = 0;
std::vector<double>* A_Ag = 0;

double sigma_Br;
std::vector<double>* px_Br = 0;
std::vector<double>* py_Br = 0;
std::vector<double>* pz_Br = 0;
std::vector<double>* Z_Br = 0;
std::vector<double>* A_Br = 0;


TFile* ReadFile1 = new TFile("../au2ag10_7AGeV_20k_latest.root");
TTree *tree1 = (TTree*) ReadFile1->Get("Glauber");
tree1->SetBranchAddress("Z_on_A", &Z_Ag);
tree1->SetBranchAddress("A_on_A", &A_Ag);
tree1->SetBranchAddress("pX_on_A", &px_Ag);
tree1->SetBranchAddress("pY_on_A", &py_Ag);
tree1->SetBranchAddress("pZ_on_A", &pz_Ag);

TFile* ReadFile2 = new TFile("../au2br10_7AGeV_20k_latest.root");
TTree *tree2 = (TTree*) ReadFile2->Get("Glauber");
tree2->SetBranchAddress("Z_on_A", &Z_Br);
tree2->SetBranchAddress("A_on_A", &A_Br);
tree2->SetBranchAddress("pX_on_A", &px_Br);
tree2->SetBranchAddress("pY_on_A", &py_Br);
tree2->SetBranchAddress("pZ_on_A", &pz_Br);


TTree *init_par_Ag = (TTree*) ReadFile1->Get("Conditions");
init_par_Ag->SetBranchAddress("Xsect_total", &sigma_Ag);
init_par_Ag->GetEntry(0);

TTree *init_par_Br = (TTree*) ReadFile2->Get("Conditions");
init_par_Br->SetBranchAddress("Xsect_total", &sigma_Br);
init_par_Br->GetEntry(0);

int nlines;
float xVal[100];
float yVal[100];
float xErr[100];
float yErr[100];

double norm1 = tree1->GetEntries();
double norm2 = tree2->GetEntries();

double nHe = 0;

	for(int k = 0; k < tree1->GetEntries(); k++){
		tree1->GetEntry(k);
		
		double pX = 0;
		double pY = 0;
		double pZ = 0;
		double pT = 0;
		double tan_theta = 0;
		double theta = 0;
		double eta = 0;			           

		for(int k = 0; k< px_Ag->size(); k++ ){pX=px_Ag->at(k);pY=py_Ag->at(k);pZ=pz_Ag->at(k);
		                                   tan_theta = std::sqrt(pX*pX + pY*pY)/pZ;
						   theta = atan(tan_theta);
						   //pT = std::sqrt(pX*pX + pY*pY)/A_Ag->at(k);
						   pT = pZ*sin(theta)/A_Ag->at(k);
						   if(Z_Ag->at(k) ==2 && A_Ag->at(k) > 2){ptHistoAg->Fill(pT);if(pT < 200) nHe +=1;}
						   
						   }

	}
	
	
	for(int k = 0; k < tree2->GetEntries(); k++){
		tree2->GetEntry(k);
		
		double pX = 0;
		double pY = 0;
		double pZ = 0;
		double pT = 0;
		double tan_theta = 0;
		double theta = 0;
		double eta = 0;			           

		for(int k = 0; k< px_Br->size(); k++ ){pX=px_Br->at(k);pY=py_Br->at(k);pZ=pz_Br->at(k);
		                                   tan_theta = std::sqrt(pX*pX + pY*pY)/pZ;
						   theta = atan(tan_theta);
						   //pT = std::sqrt(pX*pX + pY*pY)/A_Br->at(k);
						   pT = pZ*sin(theta)/A_Br->at(k);
						   if(Z_Br->at(k) ==2 && A_Br->at(k) > 2){ptHistoBr->Fill(pT); if(pT < 200) nHe +=sigma_Br/sigma_Ag;}
						   
						   }

	}
//nHe = nHe/(norm1+norm2*sigma_Br/sigma_Ag);

nlines = ReadData("cherry1994_p_t.dat", xVal, 1., yVal, 1.);


for(int i=0; i<100; i++) {
  xErr[i] = 0.025*xVal[i];
  yErr[i] = sqrt(461.*yVal[i])/461.;
}

double Int; for(int k = 1; k < nlines; k++){ Int += 0.5*(yVal[k]+yVal[k-1]);}
for(int k = 0; k < nlines; k++){yVal[k]/=Int;}
TGraphErrors* cherry1994 = new TGraphErrors(nlines-1, xVal, yVal, xErr, yErr);

  cherry1994->SetMarkerStyle(8);
  cherry1994->SetMarkerColor(kBlue);
  cherry1994->SetLineColor(kBlue);
  cherry1994->SetLineStyle(kSolid);
  cherry1994->SetMarkerSize(2.5);
  cherry1994->SetLineWidth(1);
  

//ptHistoAg->Scale(0.5/norm1);
//ptHistoBr->Scale(0.5/norm2);
ptHistoBr->SetLineColor(kBlack);
ptHistoAg->SetLineColor(kBlue);
delete ReadFile1;
delete ReadFile2;
TH1D* ptHistoEm = (TH1D*)ptHistoAg->Clone();
ptHistoEm->Add(ptHistoBr,sigma_Br/sigma_Ag);
ptHistoEm->Scale(1./nHe);

ptHistoEm->SetLineColor(kRed+1);
ptHistoEm->SetLineWidth(3.);
cherry1994->Draw("P");

//ptHistoAg->Draw("HIST SAME");
//ptHistoBr->Draw("HIST SAME");
ptHistoEm->Draw("HIST SAME");

TLegend* leg = new TLegend(0.24,0.7,0.79,0.84,"");
leg->SetTextSize(0.05);
leg->SetFillColor(0);
//leg->AddEntry(ptHistoAg,"Au+Ag at 10.7 GeV/nucleon", "lp");
//leg->AddEntry(ptHistoBr,"Au+Br at 10.7 GeV/nucleon", "lp");
leg->AddEntry(ptHistoEm,"Au-Em at 10.7 GeV/nucleon", "lp");
leg->AddEntry(cherry1994,"Cherry et al. (1994)","lp");
leg->Draw();


TLatex* latex = new TLatex( 40, 0.263, "Au+Em, T = 10.7 GeV/nucleon");
//latex->Draw();
pad->Update();

}
