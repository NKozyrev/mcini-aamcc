void FragmentationInput_com(char* input_dirname, char* output_dirname, Int_t files, Int_t A, Int_t P_beam)
{
	TChain* fChain=new TChain("events");
    fChain->Add(Form("%s/mcini_Au2_4.936GeV_20k_withMST.root", input_dirname));

    TFile* ReadFile = new TFile("../input/Au2_4.936GeV_20k_withMST.root");

    cout<<"Numer of entries equal "<<fChain->GetEntries()<<endl;
    UEvent* fEvent = new UEvent;
    EventInitialState* fIniState = new EventInitialState;
	fChain->SetBranchAddress("event", &fEvent); 
    fChain->SetBranchAddress("iniState", &fIniState); 
    gStyle -> SetOptStat(1000000001);
    gStyle -> SetStatX (0.3);

    
    TH1F* hImpactParameter =          new TH1F("b",";b, fm;counts",200,0,20);
    TH1F* hNpart =                    new TH1F("Npart",";Npart;counts",200,0,200);
    TH1F* hNspec =                    new TH1F("Nspec",";Nspec;counts",200,0,200);
    TH1F* hNcoll =                    new TH1F("Ncoll",";Ncoll;counts",1000,0,1000);
    TH2F* hNpart_vs_ImpactParameter = new TH2F("Npart_VS_B",";Npart;b, fm",200,0,200,200,0,20);
    TH2F* hNspec_vs_ImpactParameter = new TH2F("Nspec_VS_B",";Nspec;b, fm",200,0,200,200,0,20);
    TH2F* hNcoll_vs_ImpactParameter = new TH2F("Ncoll_VS_B",";Ncoll;b, fm",1000,0,1000,200,0,20);
    TH2F* hNspec_vs_Npart =           new TH2F("Nspec_VS_Npart",";Nspec;Npart",200,0,200,200,0,200);
    TH2F* hNspec_vs_Ncoll =           new TH2F("Nspec_VS_Ncoll",";Nspec;Ncoll",200,0,200,1000,0,1000);
    TH2F* hNpart_vs_Ncoll =           new TH2F("Npart_VS_Ncoll",";Npart;Ncoll",200,0,200,1000,0,1000);

    TH2F* hEcc1_vs_B =                new TH2F ("Ecc1_VS_B", ";#epsilon_{1};B, fm",200,0,1,200,0,20);
    TH2F* hEcc2_vs_B =                new TH2F ("Ecc2_VS_B", ";#epsilon_{2};B, fm",200,0,1,200,0,20);
    TH2F* hEcc3_vs_B =                new TH2F ("Ecc3_VS_B", ";#epsilon_{3};B, fm",200,0,1,200,0,20);
    TH2F* hEcc4_vs_B =                new TH2F ("Ecc4_VS_B", ";#epsilon_{4};B, fm",200,0,1,200,0,20);
    
    TH1F* hNfrag =                              new TH1F("A_Fragments",";Mass number of fragments;counts",200,0,200);
    TH1F* hEnergy =                             new TH1F("Energy",";E_{event}, GeV;counts",250,0,800);
    TH1F* hP =                                  new TH1F("P_nucleon",";P, GeV;counts",400,0,8);
    TH1F* hRapidity =                           new TH1F("Y_nucleon",";Rapidity;counts",500,0,3);
    TH1F* hPseudoRapidity =                     new TH1F("Eta_nucleon",";Pseudorapidity;counts",180,0,15);
    TH2F* hNprotons_vs_Nneutrons =              new TH2F("Nprotons_VS_Nneutrons",";N in fragment;Z in fragment",120,0,120,80,0,80);
    TH2F* hNfrag_vs_ImpactParameter =           new TH2F("A_VS_B",";Mass number of fragments;b, fm",200,0,200,200,0,20);
    TH2F* hNfrag_vs_Energy =                    new TH2F("A_VS_E",";Mass number of fragments;E_{frag}, GeV",200,0,200,150,0,800);
    TH2F* hNfrag_vs_Energy_scaled =             new TH2F("A_VS_E_scaled",";Mass number of fragments;E_{frag}-E_{beam}(A-1), GeV",200,0,200,150,-100,200);
    TH2F* hNfrag_vs_Rapidity =                  new TH2F("A_VS_Y",";Mass number of fragments;Rapidity",200,0,200,500,0,5);
    TH2F* hImpactParameter_vs_Energy_Fragment = new TH2F("B_VS_E",";b, fm;E_{event}, GeV",200,0,20,250,0,800);
    TH2F* hImpactParameter_vs_Rapidity =        new TH2F("B_VS_Y",";b, fm;Rapidity",200,0,20,250,0,3);
    TH2F* hRapidity_vs_Energy =                 new TH2F("Y_VS_E",";E, GeV;Rapidity",250,0,800,250,0,3);
    TH2F* hPseudoRapidity_vs_A =                new TH2F("Eta_VS_A",";Pseudorapidity;A",180,0,15,200,0,200);
    TH2F* hPx_vs_Py =                           new TH2F("Py_VS_Px",";Px, GeV;Py, GeV",200,0,5,200,0,5);
    
    //TH3F* hImpactParameter_Nfrag_vs_Energy_scaled = new TH3F("B_VS_A_VS_E_scaled",";b, fm;Mass number of fragments;E-E_{beam}(A-1), GeV",200,0,20,200,0,200,200,0,20);
        
	Long64_t lNEvents = fChain->GetEntries();
    Long64_t fNpa;
    UParticle* fParticle;
    
    Double_t fEnergy=0;
    Double_t px_beam, py_beam, pz_beam, E_beam;
    if (P_beam == 3) px_beam=0., py_beam=0., pz_beam=3.300000, E_beam=3.430721;
    else if (P_beam == 12) px_beam=0., py_beam=0., pz_beam=12.0, E_beam=12.036604;
    else px_beam=0., py_beam=0., pz_beam=1.0, E_beam=1.372;

            
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
                hNprotons_vs_Nneutrons                  ->Fill(fParticle->GetPdg()/10%1000-fParticle->GetPdg()/10000%1000,fParticle->GetPdg()/10000%1000);
                hNfrag                                  ->Fill(fParticle->GetPdg()/10%1000);
                hNfrag_vs_ImpactParameter               ->Fill(fParticle->GetPdg()/10%1000, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(fParticle->GetPdg()/10%1000, fMomentum.E());
                //hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), fParticle->GetPdg()/10%1000, fMomentum.E()-E_beam*(fParticle->GetPdg()/10%1000-1));
                hNfrag_vs_Energy_scaled                 ->Fill(fParticle->GetPdg()/10%1000, fMomentum.E()-E_beam*(fParticle->GetPdg()/10%1000-1));
                hNfrag_vs_Rapidity                      ->Fill(fParticle->GetPdg()/10%1000, fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fEvent->GetB(), fMomentum.E());
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
                hNfrag_vs_Energy                        ->Fill(1, fMomentum.E());
                hNfrag_vs_Rapidity                      ->Fill(1, fMomentum.Rapidity());
                //hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), 1, fMomentum.E());
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
                //hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), 1, fMomentum.E());
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
        //hImpactParameter_vs_Energy    ->Fill(fEnergy, fEvent->GetB());
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

    std::vector<Double_t>* MassOnSideA = 0;
    std::vector<Double_t>* MassOnSideB = 0;
    std::vector<Double_t>* ChargeOnSideA = 0;
    std::vector<Double_t>* ChargeOnSideB = 0;
    std::vector<Double_t>* pXonSideA = 0;
    std::vector<Double_t>* pYonSideA = 0;
    std::vector<Double_t>* pZonSideA = 0;
    std::vector<Double_t>* pXonSideB = 0;
    std::vector<Double_t>* pYonSideB = 0;
    std::vector<Double_t>* pZonSideB = 0;
    Float_t b;
    Double_t EnergyProj, EnergyTarg, Energy_frag;
    Int_t sumZ_proj, sumZ_targ;

    TH2F* hEnergyA_vs_ImpactParameter_proj =          new TH2F("EperA_VS_B_proj",";b, fm;E_{frag}/A, GeV",200,0,20,200,0,8);
    TH2F* hEnergyA_vs_ImpactParameter_targ =          new TH2F("EperA_VS_B_targ",";b, fm;E_{frag}/A, GeV",200,0,20,200,0,8);
    TH2F* hNspect_vs_Espect_proj =          new TH2F("Espect_VS_Nspect_proj",";Nspect;Espect, GeV",100,0,100,150,0,800);
    TH2F* hNspect_vs_Espect_targ =          new TH2F("Espect_VS_Nspect_targ",";Nspect;Espect, GeV",100,0,100,150,0,800);
    TH2F* hNspect_vs_sumZ_proj =          new TH2F("sumZ_VS_Nspect_proj",";sumZ;Nspect",80,0,80,110,0,110);
    TH2F* hNspect_vs_sumZ_targ =          new TH2F("sumZ_VS_Nspect_targ",";sumZ;Nspect",80,0,80,110,0,110);
    
    TH2F* hNspectT_vs_NspectP =          new TH2F("Nspect_proj_VS_Nspect_targ",";Nspect(target);Nspect(projectile)",110,0,110,110,0,110);
    TH2F* hNsumZT_vs_NsumZP =          new TH2F("sumZ_proj_VS_sumZ_targ",";sumZ(target);sumZ(projectile)",80,0,80,80,0,80);
    TH2F* hNnucl_vs_Nfrag_proj =          new TH2F("Nnucl_vs_Nfrag_proj",";Nfrag;Nnucl",30,0,30,90,0,90);
    TH2F* hNnucl_vs_Nfrag_targ =          new TH2F("Nnucl_vs_Nfrag_targ",";Nfrag;Nnucl",30,0,30,90,0,90);

    TH2F* hEnergyE_vs_sumZ_proj =          new TH2F("EnergyE_VS_sumZ_proj",";sumZ;E_{event}, GeV",80,0,80,175,0,800);
    TH2F* hEnergyE_vs_Nnucl_proj =          new TH2F("EnergyE_VS_Nnucl_proj",";Nnucl;E_{event}, GeV",90,0,90,175,0,800);
    TH2F* hEnergyE_vs_Nimf_proj =          new TH2F("EnergyE_VS_Nimf_proj",";N_{IMF};E_{event}, GeV",15,0,15,30,0,800);

    TH2F* hEnergyProj_vs_EnergyTarg        =          new TH2F("E_proj_VS_E_targ",";E_{event} (target), GeV;E_{event} (projectile), GeV",150,0,800,150,0,800);
    TH2F* hImpactParameter_vs_Energy_targ  =       new TH2F("B_VS_Energy_targ",";b, fm;E_{event} (target), GeV",200,0,20,250,0,800);
    TH2F* hImpactParameter_vs_Energy_proj  =       new TH2F("B_VS_Energy_proj",";b, fm;E_{event} (projectile), GeV",200,0,20,250,0,800);



    TTree *tree = (TTree*) ReadFile->Get("Glauber");  
    tree->SetBranchAddress("A_on_A", &MassOnSideA);
    tree->SetBranchAddress("A_on_B", &MassOnSideB);
    tree->SetBranchAddress("Z_on_A", &ChargeOnSideA); 
    tree->SetBranchAddress("Z_on_B", &ChargeOnSideB);
    tree->SetBranchAddress("pX_on_A", &pXonSideA);
    tree->SetBranchAddress("pY_on_A", &pYonSideA);
    tree->SetBranchAddress("pZ_on_A", &pZonSideA);
    tree->SetBranchAddress("pX_on_B", &pXonSideB);
    tree->SetBranchAddress("pY_on_B", &pYonSideB);
    tree->SetBranchAddress("pZ_on_B", &pZonSideB);
    tree->SetBranchAddress("impact_parameter", &b);

    for(int k = 0; k < tree->GetEntries(); k++){
        tree->GetEntry(k);
        EnergyProj=0; EnergyTarg=0; sumZ_proj =0; sumZ_targ=0;
        int Nucl_T=0, Frag_T=0,Nucl_P=0, Frag_P=0, IMF_P=0, IMF_T=0;
        for (Int_t ipart = 0; ipart < (MassOnSideA->size()); ipart++){
            Energy_frag = pow(pow(pXonSideA->at(ipart) / 1000, 2) + pow(pYonSideA->at(ipart) / 1000, 2) + pow(pZonSideA->at(ipart) / 1000, 2) + pow(0.9395654 * (MassOnSideA->at(ipart) - ChargeOnSideA->at(ipart)) + 0.9382721 * ChargeOnSideA->at(ipart), 2), 0.5);
            hEnergyA_vs_ImpactParameter_proj -> Fill(b, Energy_frag/MassOnSideA->at(ipart));
            if (MassOnSideA->at(ipart) == 1) Nucl_P+=1;
            else Frag_P+=1;
            if (ChargeOnSideA->at(ipart) >= 3 && ChargeOnSideA->at(ipart) <= 30) IMF_P+=1;
            if (Energy_frag/MassOnSideA->at(ipart) > 17) cout << "A = " << MassOnSideA->at(ipart) << "; Z = " << ChargeOnSideA->at(ipart) << "; Pz[MeV] = " << pZonSideA->at(ipart) <<"; Pt[MeV]= " << pow (pow (pYonSideA->at(ipart),2) + pow (pXonSideA->at(ipart),2), 0.5) <<endl;
            EnergyProj+=Energy_frag;
            sumZ_proj+=ChargeOnSideA->at(ipart);

        }
        for (Int_t ipart = 0; ipart < (MassOnSideB->size()); ipart++){
            hEnergyA_vs_ImpactParameter_targ -> Fill(b, pow(pow(pXonSideB->at(ipart) / 1000, 2) + pow(pYonSideB->at(ipart) / 1000, 2) + pow(pZonSideB->at(ipart) / 1000, 2) + pow(0.9395654 * (MassOnSideB->at(ipart) - ChargeOnSideB->at(ipart)) + 0.9382721 * ChargeOnSideB->at(ipart), 2), 0.5)/MassOnSideB->at(ipart));
            if (MassOnSideB->at(ipart) == 1) Nucl_T+=1;
            else Frag_T+=1;
            if (ChargeOnSideB->at(ipart) >= 3 && ChargeOnSideB->at(ipart) <= 30) IMF_T+=1;
            EnergyTarg+=pow(pow(pXonSideB->at(ipart) / 1000, 2) + pow(pYonSideB->at(ipart) / 1000, 2) + pow(pZonSideB->at(ipart) / 1000, 2) + pow(0.9395654 * (MassOnSideB->at(ipart) - ChargeOnSideB->at(ipart)) + 0.9382721 * ChargeOnSideB->at(ipart), 2), 0.5);
            sumZ_targ+=ChargeOnSideB->at(ipart);
        }
        hEnergyProj_vs_EnergyTarg -> Fill (EnergyTarg, EnergyProj);
        hNspect_vs_Espect_proj -> Fill (MassOnSideA->size(), EnergyProj);
        hNspect_vs_Espect_targ -> Fill (MassOnSideB->size(), EnergyTarg);
        hNspect_vs_sumZ_proj -> Fill (sumZ_proj, MassOnSideA->size());
        hNspect_vs_sumZ_targ -> Fill (sumZ_targ, MassOnSideB->size());
        hImpactParameter_vs_Energy_proj -> Fill (b, EnergyProj);
        hImpactParameter_vs_Energy_targ -> Fill (b, EnergyTarg);

        hEnergyE_vs_sumZ_proj -> Fill (sumZ_proj, EnergyProj);
        hEnergyE_vs_Nnucl_proj -> Fill (Nucl_P, EnergyProj);
        hEnergyE_vs_Nimf_proj -> Fill (IMF_P, EnergyProj);

        hNspectT_vs_NspectP -> Fill (MassOnSideB->size(), MassOnSideA->size());
        hNsumZT_vs_NsumZP -> Fill (sumZ_targ, sumZ_proj);
        hNnucl_vs_Nfrag_proj -> Fill (Frag_P, Nucl_P);
        hNnucl_vs_Nfrag_targ -> Fill (Frag_T, Nucl_T);
    }



    TCanvas* canvas = new TCanvas("canvas");
    gPad->SetTickx();
    gPad->SetTicky();
    hImpactParameter           ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf(");

    gPad->SetLogz();

   
    hEnergyE_vs_sumZ_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEnergyE_vs_Nnucl_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEnergyE_vs_Nimf_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

    gStyle -> SetOptStat(1000100001);
     hEnergyA_vs_ImpactParameter_targ ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEnergyA_vs_ImpactParameter_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEnergyProj_vs_EnergyTarg ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    
    hNspect_vs_sumZ_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNspect_vs_sumZ_targ ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
 
    
    hNspect_vs_Espect_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNspect_vs_Espect_targ ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

    gStyle -> SetOptStat(1000000001);
    
    hImpactParameter_vs_Energy_targ ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hImpactParameter_vs_Energy_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

        hNspectT_vs_NspectP ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
        hNsumZT_vs_NsumZP ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
        hNnucl_vs_Nfrag_proj ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
        hNnucl_vs_Nfrag_targ ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
   


    /*hNpart                     ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNspec                     ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNcoll                     ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNpart_vs_ImpactParameter  ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNspec_vs_ImpactParameter  ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNcoll_vs_ImpactParameter  ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNspec_vs_Npart            ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNspec_vs_Ncoll            ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNpart_vs_Ncoll            ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEcc1_vs_B                 ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEcc2_vs_B                 ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEcc3_vs_B                 ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hEcc4_vs_B                 ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");*/
    
    gPad->SetLogy();
    hNfrag                               ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    gPad->SetLogy(0);
    hEnergy                              ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    gPad->SetLogy();
    hP                                   ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hRapidity                            ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    gPad->SetLogy(0);
    hNprotons_vs_Nneutrons               ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hNfrag_vs_ImpactParameter            ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

    hNfrag_vs_Energy                     ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    //hNfrag_vs_Energy_scaled              ->Draw("COLZ");
    //canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

    hNfrag_vs_Rapidity                   ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

    //hImpactParameter_vs_Energy           ->Draw("COLZ");
    //canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");

    hImpactParameter_vs_Energy_Fragment  ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hImpactParameter_vs_Rapidity         ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hRapidity_vs_Energy                  ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hPx_vs_Py                            ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hPseudoRapidity                      ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf");
    hPseudoRapidity_vs_A                 ->Draw("COLZ");
    canvas->Print("../output/histos/Au2_4.936GeV_20k_withMST.pdf)");
       

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
    //hNfrag_vs_Energy_scaled              ->Write();
    hNfrag_vs_Rapidity                   ->Write();
    //hImpactParameter_vs_Energy           ->Write();
    hImpactParameter_vs_Energy_Fragment  ->Write();
    hImpactParameter_vs_Rapidity         ->Write();
    hRapidity_vs_Energy                  ->Write();
    hPx_vs_Py                            ->Write();
    hPseudoRapidity                      ->Write();
    hPseudoRapidity_vs_A                 ->Write();

    hEnergyA_vs_ImpactParameter_proj -> Write();
    hEnergyA_vs_ImpactParameter_targ -> Write();
     hNspect_vs_Espect_proj -> Write();
   hNspect_vs_Espect_targ -> Write();
    hNspect_vs_sumZ_proj -> Write();
   hNspect_vs_sumZ_targ -> Write();
    hEnergyProj_vs_EnergyTarg        -> Write();
    hImpactParameter_vs_Energy_targ   -> Write();
   hImpactParameter_vs_Energy_proj   -> Write();

    hNspectT_vs_NspectP -> Write();
        hNsumZT_vs_NsumZP -> Write();;
        hNnucl_vs_Nfrag_proj -> Write();
        hNnucl_vs_Nfrag_targ -> Write();

        hEnergyE_vs_sumZ_proj -> Write();
        hEnergyE_vs_Nnucl_proj -> Write();
        hEnergyE_vs_Nimf_proj -> Write();
    
    //hImpactParameter_Nfrag_vs_Energy_scaled ->Write();
    
	fout->Close();
}
