{
std::vector<double> MassOnSideA;
std::vector<double> ChargeOnSideA;
std::vector<double> pXonSideA;
std::vector<double> pYonSideA;
std::vector<double> pZonSideA;
float btree;
TFile* ReadFile = new TFile("dcmqgsm_in_aamcc_format_sps.root","RECREATE","Demo ROOT file with histograms & trees");
TTree* tree = new TTree("Glauber","Events from glauber modeling");

tree->Branch("A_on_A", "std::vector" ,&MassOnSideA);
tree->Branch("Z_on_A", "std::vector" ,&ChargeOnSideA);
tree->Branch("pX_on_A", "std::vector" ,&pXonSideA);
tree->Branch("pY_on_A", "std::vector" ,&pYonSideA);
tree->Branch("pZ_on_A", "std::vector" ,&pZonSideA);
tree->Branch("impact_parameter", &btree, "impact_parameter/f");

TChain* fChain=new TChain("events");
for (int i = 0; i <= 49; i++) {fChain->Add(Form("%s/dcmqgsm_%d.root", "/media/sf_Ubuntu_share/sps_dcm", i));}
    //TFile* ReadFile = new TFile("../input/QA_dcmqgsm.root");

cout<<"Numer of entries equal "<<fChain->GetEntries()<<endl;
UEvent* fEvent = new UEvent;
EventInitialState* fIniState = new EventInitialState;
fChain->SetBranchAddress("event", &fEvent); 
fChain->SetBranchAddress("iniState", &fIniState); 

Double_t px_beam=0., py_beam=0., pz_beam=159.0, E_beam=159.003;
TLorentzVector Beam(px_beam, py_beam, pz_beam, E_beam);
TLorentzVector Target(0., 0., 0, 0.938);
TLorentzRotation rTransform;
rTransform.Boost((Beam + Target).BoostVector());
    
Long64_t lNEvents = fChain->GetEntries();
Long64_t fNpa;
UParticle* fParticle;
                    
  for (long i = 0; i < lNEvents; i++)
  {
        fChain->GetEntry(i);
        fNpa = fEvent->GetNpa();
        btree = fEvent->GetB();
        for (int j=0;j<fNpa;j++)
        {
            fParticle = fEvent->GetParticle(j);
            TLorentzVector fMomentum = (fParticle->GetMomentum()).Transform(rTransform);
            //TLorentzVector fMomentum = fParticle->GetMomentum();
            if (fParticle -> T() == 1)
            {
              if (fParticle->GetPdg()>1e9) {
                MassOnSideA.push_back(fParticle->GetPdg()/10%1000);
                ChargeOnSideA.push_back(fParticle->GetPdg()/10000%1000);
                pXonSideA.push_back(fMomentum.Px()*1000);
                pYonSideA.push_back(fMomentum.Py()*1000);
                pZonSideA.push_back(fMomentum.Pz()*1000);
              }
              else if (fParticle->GetPdg()==2212) {
                MassOnSideA.push_back(1.);
                ChargeOnSideA.push_back(1.);
                pXonSideA.push_back(fMomentum.Px()*1000);
                pYonSideA.push_back(fMomentum.Py()*1000);
                pZonSideA.push_back(fMomentum.Pz()*1000);
              }
              else if (fParticle->GetPdg()==2112) {
                MassOnSideA.push_back(1.);
                ChargeOnSideA.push_back(0.);
                pXonSideA.push_back(fMomentum.Px()*1000);
                pYonSideA.push_back(fMomentum.Py()*1000);
                pZonSideA.push_back(fMomentum.Pz()*1000);
              }
            }
        }
        tree->Fill();
        MassOnSideA.clear();
        ChargeOnSideA.clear();
        pXonSideA.clear();
        pYonSideA.clear();
        pZonSideA.clear();
    }
  ReadFile->Write();
  ReadFile->Close();
}