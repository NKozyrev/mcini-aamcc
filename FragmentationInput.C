void FragmentationInput(char* input_dirname, char* output_dirname, Int_t files, Int_t A, Int_t P_beam)
{
	TChain* fChain=new TChain("events");
    fChain->Add(Form("%s/mcini_Au2_11GeV_50k.root", input_dirname));
    cout<<"Numer of entries equal "<<fChain->GetEntries()<<endl;
    UEvent* fEvent = new UEvent;
    EventInitialState* fIniState = new EventInitialState;
	fChain->SetBranchAddress("event", &fEvent); 
    fChain->SetBranchAddress("iniState", &fIniState); 
    
    TH1F* hImpactParameter =          new TH1F("B",";B, fm;counts",200,0,20);
    TH1F* hNpart =                    new TH1F("Npart",";Npart;counts",200,0,200);
    TH1F* hNspec =                    new TH1F("Nspec",";Nspec;counts",200,0,200);
    TH1F* hNcoll =                    new TH1F("Ncoll",";Ncoll;counts",1000,0,1000);
    TH2F* hNpart_vs_ImpactParameter = new TH2F("Npart_VS_B",";Npart;B, fm",200,0,200,200,0,20);
    TH2F* hNspec_vs_ImpactParameter = new TH2F("Nspec_VS_B",";Nspec;B, fm",200,0,200,200,0,20);
    TH2F* hNcoll_vs_ImpactParameter = new TH2F("Ncoll_VS_B",";Ncoll;B, fm",1000,0,1000,200,0,20);
    TH2F* hNspec_vs_Npart =           new TH2F("Nspec_VS_Npart",";Nspec;Npart",200,0,200,200,0,200);
    TH2F* hNspec_vs_Ncoll =           new TH2F("Nspec_VS_Ncoll",";Nspec;Ncoll",200,0,200,1000,0,1000);
    TH2F* hNpart_vs_Ncoll =           new TH2F("Npart_VS_Ncoll",";Npart;Ncoll",200,0,200,1000,0,1000);

    TH2F* hEcc1_vs_B =                new TH2F ("Ecc1_VS_B", ";Ecc1;B, fm",200,0,1,200,0,20);
    TH2F* hEcc2_vs_B =                new TH2F ("Ecc2_VS_B", ";Ecc2;B, fm",200,0,1,200,0,20);
    TH2F* hEcc3_vs_B =                new TH2F ("Ecc3_VS_B", ";Ecc3;B, fm",200,0,1,200,0,20);
    TH2F* hEcc4_vs_B =                new TH2F ("Ecc4_VS_B", ";Ecc4;B, fm",200,0,1,200,0,20);
    
    TH1F* hNfrag =                              new TH1F("A_Fragments",";Mass number of fragments;counts",200,0,200);
    TH1F* hEnergy =                             new TH1F("Energy",";E, GeV;counts",250,0,3500);
    TH1F* hP =                                  new TH1F("P_nucleon",";P, GeV;counts",200,0,20);
    TH1F* hRapidity =                           new TH1F("Y_nucleon",";Rapidity;counts",500,0,5);
    TH1F* hPseudoRapidity =                     new TH1F("Eta_nucleon",";Pseudorapidity;counts",500,0,10);
    TH2F* hNprotons_vs_Nneutrons =              new TH2F("Nprotons_VS_Nneutrons",";Nprotons;Nneutrons",200,0,200,200,0,200);
    TH2F* hNfrag_vs_ImpactParameter =           new TH2F("A_VS_B",";Mass number of fragments;B, fm",200,0,200,200,0,20);
    TH2F* hNfrag_vs_Energy =                    new TH2F("A_VS_E",";Mass number of fragments;E, GeV",200,0,200,250,0,3500);
    TH2F* hNfrag_vs_Energy_scaled =             new TH2F("A_VS_E_scaled",";Mass number of fragments;E-E_{beam}(A-1), GeV",200,0,200,400,-100,200);
    TH2F* hNfrag_vs_Rapidity =                  new TH2F("A_VS_Y",";Mass number of fragments;Rapidity",200,0,200,500,0,5);
    TH2F* hImpactParameter_vs_Energy =          new TH2F("B_VS_Energy",";E, GeV;B, fm",250,0,3500,200,0,20);
    TH2F* hImpactParameter_vs_Energy_Fragment = new TH2F("B_VS_E",";E, GeV;B, fm",250,0,3500,200,0,20);
    TH2F* hImpactParameter_vs_Rapidity =        new TH2F("B_VS_Y",";B, fm;Rapidity",200,0,20,500,0,5);
    TH2F* hRapidity_vs_Energy =                 new TH2F("Y_VS_E",";E, GeV;Rapidity",250,0,3500,500,0,5);
    TH2F* hPseudoRapidity_vs_A =                new TH2F("Eta_VS_A",";Pseudorapidity;A",500,0,10,200,0,200);
    TH2F* hPx_vs_Py =                           new TH2F("Py_VS_Px",";Px, GeV;Py, GeV",200,0,5,200,0,5);
    
    TH3F* hImpactParameter_Nfrag_vs_Energy_scaled = new TH3F("B_VS_A_VS_E_scaled",";B, fm;Mass number of fragments;E-E_{beam}(A-1), GeV",200,0,20,200,0,200,200,0,20);
        
	Long64_t lNEvents = fChain->GetEntries();
    Long64_t fNpa;
    UParticle* fParticle;
    
    Double_t fEnergy=0;
    Double_t px_beam, py_beam, pz_beam, E_beam;
    if (P_beam == 3) px_beam=0., py_beam=0., pz_beam=3.300000, E_beam=3.430721;
    else if (P_beam == 12) px_beam=0., py_beam=0., pz_beam=12.0, E_beam=12.036604;
    else px_beam=0., py_beam=0., pz_beam=11.0, E_beam=11.04;

            
	for (long i = 0; i < lNEvents; i++)
	{
		if (i % 100 == 0) cout<<i<<endl;
        fChain->GetEntry(i);
        fEnergy=0;
        fNpa = fEvent->GetNpa();

        for (int j=0;j<fNpa;j++)
        {
            fParticle = fEvent->GetParticle(j);
            TLorentzVector fMomentum = fParticle->GetMomentum();
            
            if (fParticle->GetPdg()>1e9) 
            {
                hNprotons_vs_Nneutrons                  ->Fill(fParticle->GetPdg()/10000%1000, fParticle->GetPdg()/10%1000-fParticle->GetPdg()/10000%1000);
                hNfrag                                  ->Fill(fParticle->GetPdg()/10%1000);
                hNfrag_vs_ImpactParameter               ->Fill(fParticle->GetPdg()/10%1000, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(fParticle->GetPdg()/10%1000, fMomentum.E());
                hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), fParticle->GetPdg()/10%1000, fMomentum.E()-E_beam*(fParticle->GetPdg()/10%1000-1));
                hNfrag_vs_Energy_scaled                 ->Fill(fParticle->GetPdg()/10%1000, fMomentum.E()-E_beam*(fParticle->GetPdg()/10%1000-1));
                hNfrag_vs_Rapidity                      ->Fill(fParticle->GetPdg()/10%1000, fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fMomentum.E(), fEvent->GetB());
                hImpactParameter_vs_Rapidity            ->Fill(fEvent->GetB(), fMomentum.Rapidity());
                hPx_vs_Py                               ->Fill(fMomentum.Px(), fMomentum.Py());
                hRapidity_vs_Energy                     ->Fill(fMomentum.E(), fMomentum.Rapidity());
                hP                                      ->Fill(fMomentum.P()-P_beam*(fParticle->GetPdg()/10%1000-1));
                hRapidity                               ->Fill(fMomentum.Rapidity());
                hPseudoRapidity                         ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_vs_A                    ->Fill(fMomentum.PseudoRapidity(),fParticle->GetPdg()/10%1000);
                fEnergy=fEnergy+fMomentum.E();
            }       
            else if (fParticle->GetPdg()==2212) 
            {
                hNprotons_vs_Nneutrons                  ->Fill(1, 0.1);
                hNfrag                                  ->Fill(1);
                hNfrag_vs_ImpactParameter               ->Fill(1, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(1,  fMomentum.E());
                hNfrag_vs_Rapidity                      ->Fill(1, fMomentum.Rapidity());
                hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), 1, fMomentum.E());
                hNfrag_vs_Energy_scaled                 ->Fill(1, fMomentum.E());
                hImpactParameter_vs_Rapidity            ->Fill(fEvent->GetB(), fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fEvent->GetB(), fMomentum.E());
                hPx_vs_Py                               ->Fill(fMomentum.Px(), fMomentum.Py());
                hRapidity_vs_Energy                     ->Fill(fMomentum.E(), fMomentum.Rapidity());                
                hP                                      ->Fill(fMomentum.P());
                hRapidity                               ->Fill(fMomentum.Rapidity());
                hPseudoRapidity                         ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_vs_A                    ->Fill(fMomentum.PseudoRapidity(),1);

                fEnergy=fEnergy+fMomentum.E();
            }
            else if (fParticle->GetPdg()==2112) 
            {
                hNprotons_vs_Nneutrons                  ->Fill(0.1, 1);
                hNfrag                                  ->Fill(1);
                hNfrag_vs_ImpactParameter               ->Fill(1, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(1, fMomentum.E());
                hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), 1, fMomentum.E());
                hNfrag_vs_Energy_scaled                 ->Fill(1, fMomentum.E());
                hNfrag_vs_Rapidity                      ->Fill(1, fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fEvent->GetB(), fMomentum.E());
                hImpactParameter_vs_Rapidity            ->Fill(fEvent->GetB(), fMomentum.Rapidity());
                hPx_vs_Py                               ->Fill(fMomentum.Px(), fMomentum.Py());
                hRapidity_vs_Energy                     ->Fill(fMomentum.E(), fMomentum.Rapidity());
                hP                                      ->Fill(fMomentum.P());
                hRapidity                               ->Fill(fMomentum.Rapidity());
                hPseudoRapidity                         ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_vs_A                    ->Fill(fMomentum.PseudoRapidity(),1);
                fEnergy=fEnergy+fMomentum.E();
            }
        }
        
        
        hEnergy                       ->Fill(fEnergy);
        hImpactParameter_vs_Energy    ->Fill(fEnergy, fEvent->GetB());
        hImpactParameter              ->Fill(fEvent->GetB());
        hNpart                        ->Fill(fIniState->getNPart());
        hNspec                        ->Fill(A-fIniState->getNPart());
        hNcoll                        ->Fill(fIniState->getNColl());
        hNpart_vs_ImpactParameter     ->Fill(fIniState->getNPart(), fEvent->GetB());
        hNspec_vs_ImpactParameter     ->Fill(A-fIniState->getNPart(), fEvent->GetB());
        hNcoll_vs_ImpactParameter     ->Fill(fIniState->getNColl(), fEvent->GetB());
        hNspec_vs_Npart               ->Fill(A-fIniState->getNPart(), fIniState->getNPart());
        hNspec_vs_Ncoll               ->Fill(A-fIniState->getNPart(), fIniState->getNColl());
        hNpart_vs_Ncoll               ->Fill(fIniState->getNPart(), fIniState->getNColl());

        const Float_t * Ecc = fEvent -> GetEcc();
        hEcc1_vs_B                    ->Fill(Ecc[1], fEvent->GetB());
        hEcc2_vs_B                    ->Fill(Ecc[2], fEvent->GetB());
        hEcc3_vs_B                    ->Fill(Ecc[3], fEvent->GetB());
        hEcc4_vs_B                    ->Fill(Ecc[4], fEvent->GetB());
		
	}
    TCanvas* canvas = new TCanvas("canvas");
    hImpactParameter           ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf(");
    hNpart                     ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNspec                     ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNcoll                     ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNpart_vs_ImpactParameter  ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNspec_vs_ImpactParameter  ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNcoll_vs_ImpactParameter  ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNspec_vs_Npart            ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNspec_vs_Ncoll            ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNpart_vs_Ncoll            ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hEcc1_vs_B                 ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hEcc2_vs_B                 ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hEcc3_vs_B                 ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hEcc4_vs_B                 ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    
    gPad->SetLogy();
    hNfrag                               ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    gPad->SetLogy(0);
    hEnergy                              ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hP                                   ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hRapidity                            ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNprotons_vs_Nneutrons               ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNfrag_vs_ImpactParameter            ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNfrag_vs_Energy                     ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNfrag_vs_Energy_scaled              ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hNfrag_vs_Rapidity                   ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hImpactParameter_vs_Energy           ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hImpactParameter_vs_Energy_Fragment  ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hImpactParameter_vs_Rapidity         ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hRapidity_vs_Energy                  ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hPx_vs_Py                            ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hPseudoRapidity                      ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    hPseudoRapidity_vs_A                 ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf");
    
    hImpactParameter_Nfrag_vs_Energy_scaled ->Draw();
    canvas->Print("../output/histos/QA_Au2_11GeV_50k.pdf)");    

    TFile *fout = new TFile(output_dirname, "recreate");
    
    hImpactParameter           ->Write();
    hNpart                     ->Write();
    hNspec                     ->Write();
    hNcoll                     ->Write();
    hNpart_vs_ImpactParameter  ->Write();
    hNspec_vs_ImpactParameter  ->Write();
    hNcoll_vs_ImpactParameter  ->Write();
    hNspec_vs_Npart            ->Write();
    hNspec_vs_Ncoll            ->Write();
    hNpart_vs_Ncoll            ->Write();
    hEcc1_vs_B                 ->Write();
    hEcc2_vs_B                 ->Write();
    hEcc3_vs_B                 ->Write();
    hEcc4_vs_B                 ->Write();
    
    gPad->SetLogy();
    hNfrag                               ->Write();
    gPad->SetLogy(0);
    hEnergy                              ->Write();
    hP                                   ->Write();
    hRapidity                            ->Write();
    hNprotons_vs_Nneutrons               ->Write();
    hNfrag_vs_ImpactParameter            ->Write();
    hNfrag_vs_Energy                     ->Write();
    hNfrag_vs_Energy_scaled              ->Write();
    hNfrag_vs_Rapidity                   ->Write();
    hImpactParameter_vs_Energy           ->Write();
    hImpactParameter_vs_Energy_Fragment  ->Write();
    hImpactParameter_vs_Rapidity         ->Write();
    hRapidity_vs_Energy                  ->Write();
    hPx_vs_Py                            ->Write();
    hPseudoRapidity                      ->Write();
    hPseudoRapidity_vs_A                 ->Write();
    
    hImpactParameter_Nfrag_vs_Energy_scaled ->Write();
    
	fout->Close();
}
