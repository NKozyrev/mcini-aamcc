void FragmentationInput_mcini_pb_lhc(char* input_dirname, char* output_dirname, Int_t files, Int_t A, Int_t P_beam)
{
	TChain* fChain=new TChain("events");
    //fChain->Add(Form("%s/mcini_auau12AGeV_dcmqgsm.root", input_dirname));
    for (int i = 0; i <= files; i++) {fChain->Add(Form("%s/dcmqgsm_%d.root", input_dirname, i));}
    //TFile* ReadFile = new TFile("../input/QA_dcmqgsm_lhc.root");

    cout<<"Numer of entries equal "<<fChain->GetEntries()<<endl;
    UEvent* fEvent = new UEvent;
    EventInitialState* fIniState = new EventInitialState;
	fChain->SetBranchAddress("event", &fEvent); 
    fChain->SetBranchAddress("iniState", &fIniState); 
    gStyle -> SetOptStat(1000000001);
    gStyle -> SetStatX (0.5);
    TGaxis::SetMaxDigits(3);

    
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
    
    TH1F* hNfrag =                              new TH1F("A_Fragments",";Mass number of fragments;counts",210,0,210);
    TH1F* hEnergy =                             new TH1F("Energy",";E_{event}, GeV;counts",250,0,530000);
    TH1F* hP =                                  new TH1F("P",";P/A, GeV;counts",400,2000,3000);
    TH1F* hPTnucl =                             new TH1F("Pt_nucl",";P_{t}, GeV;counts",100,0,1);
    TH2F* hPT_vs_A =                            new TH2F("Pt_vs_A",";P_{t}, GeV;A",100,0,3,208,0,208);
    TH1F* hPnucl =                              new TH1F("P_z_nucleons",";P/A, GeV;counts",400,2000,3000);
    TH1F* hPimf =                               new TH1F("P_z_imf",";P/A, GeV;counts",400,2000,3000);
    TH1F* hPheavy =                             new TH1F("P_z_heavy",";P/A, GeV;counts",400,2000,3000);
    TH1F* hRapidity =                           new TH1F("Y_nucleon",";Rapidity;counts",200,8,15);
    TH1F* hPseudoRapidity =                     new TH1F("Eta",";Pseudorapidity;counts",180,0,18);
    TH1F* hPseudoRapidity_nucl =                new TH1F("Eta_nucl",";Pseudorapidity (Z=1);counts",180,0,18);
    TH2F* hNprotons_vs_Nneutrons =              new TH2F("Nprotons_VS_Nneutrons",";N in fragment;Z in fragment",126,0,126,83,0,83);
    TH2F* hNfrag_vs_ImpactParameter =           new TH2F("A_VS_B",";Mass number of fragments;b, fm",210,0,210,200,0,20);
    TH2F* hNfrag_vs_Energy =                    new TH2F("A_VS_E",";Mass number of fragments;E_{frag}, GeV",200,0,200,250,0,530000);
    TH2F* hNfrag_vs_Energy_scaled =             new TH2F("A_VS_E_scaled",";Mass number of fragments;E_{frag}-E_{beam}(A-1), GeV",200,0,200,150,-100,200);
    TH2F* hNfrag_vs_Rapidity =                  new TH2F("A_VS_Y",";Mass number of fragments;Rapidity",200,0,200,500,8,15);
    TH2F* hImpactParameter_vs_Energy_Fragment = new TH2F("B_VS_E",";b, fm;E_{event}, GeV",200,0,20,250,0,530000);
    TH2F* hImpactParameter_vs_Rapidity =        new TH2F("B_VS_Y",";b, fm;Rapidity",200,0,20,250,8,15);
    TH2F* hRapidity_vs_Energy =                 new TH2F("Y_VS_E",";E, GeV;Rapidity",360,0,530000,250,8,15);
    TH2F* hPseudoRapidity_vs_A =                new TH2F("Eta_VS_A",";Pseudorapidity;A",180,0,18,210,0,210);
    TH2F* hPx_vs_Py =                           new TH2F("Py_VS_Px",";Px, GeV;Py, GeV",200,0,5,200,0,5);

    TH2F* hEnergyA_vs_ImpactParameter_proj =          new TH2F("EperA_VS_B_proj",";b, fm;E_{frag}/A, GeV",200,0,20,250,0,3000);
    //TH2F* hEnergyA_vs_ImpactParameter_targ =          new TH2F("EperA_VS_B_targ",";b, fm;E_{frag}/A, GeV",200,0,20,200,0,4);
    TH2F* hNspect_vs_Espect_proj =          new TH2F("Espect_VS_Nspect_proj",";Nspect;Espect, GeV",100,0,100,250,0,530000);
    //TH2F* hNspect_vs_Espect_targ =          new TH2F("Espect_VS_Nspect_targ",";Nspect;Espect, GeV",100,0,100,150,0,300);
    TH2F* hNspect_vs_sumZ_proj =          new TH2F("sumZ_VS_Nspect_proj",";sumZ;Nspect",83,0,83,110,0,110);
    TH2F* hNspect_vs_Zb2 =          new TH2F("Zb3_VS_Nspect_proj",";Zb3;Nspect",83,0,83,100,0,100);
    //TH2F* hNspect_vs_sumZ_targ =          new TH2F("sumZ_VS_Nspect_targ",";sumZ;Nspect",80,0,80,110,0,110);
    
    //TH2F* hNspectT_vs_NspectP =          new TH2F("Nspect_proj_VS_Nspect_targ",";Nspect(target);Nspect(projectile)",110,0,110,110,0,110);
    //TH2F* hNsumZT_vs_NsumZP =          new TH2F("sumZ_proj_VS_sumZ_targ",";sumZ(target);sumZ(projectile)",80,0,80,80,0,80);
    TH2F* hNnucl_vs_Nfrag_proj =          new TH2F("Nnucl_vs_Nfrag_proj",";Nfrag;Nnucl",40,0,40,90,0,90);
    //TH2F* hNnucl_vs_Nfrag_targ =          new TH2F("Nnucl_vs_Nfrag_targ",";Nfrag;Nnucl",30,0,30,90,0,90);

    TH2F* hEnergyE_vs_sumZ_proj =          new TH2F("EnergyE_VS_sumZ_proj",";sumZ;E_{event}, GeV",83,0,83,250,0,530000);
    TH2F* hEnergyE_vs_Nnucl_proj =          new TH2F("EnergyE_VS_Nnucl_proj",";Nnucl;E_{event}, GeV",90,0,90,250,0,530000);
    TH2F* hEnergyE_vs_Nimf_proj =          new TH2F("EnergyE_VS_Nimf_proj",";N_{IMF};E_{event}, GeV",15,0,15,53,0,530000);

    //TH2F* hEnergyProj_vs_EnergyTarg        =          new TH2F("E_proj_VS_E_targ",";E_{event} (target), GeV;E_{event} (projectile), GeV",150,0,300,150,0,530000);
    //TH2F* hImpactParameter_vs_Energy_targ  =       new TH2F("B_VS_Energy_targ",";b, fm;E_{event} (target), GeV",200,0,20,250,0,300);
    TH2F* hImpactParameter_vs_Energy  =       new TH2F("B_VS_Energy_proj",";b, fm;E_{event} (projectile), GeV",200,0,20,250,0,530000);

    
    //TH3F* hImpactParameter_Nfrag_vs_Energy_scaled = new TH3F("B_VS_A_VS_E_scaled",";b, fm;Mass number of fragments;E-E_{beam}(A-1), GeV",200,0,20,200,0,200,200,0,20);
        
	Long64_t lNEvents = fChain->GetEntries();
    Long64_t fNpa;
    UParticle* fParticle;
    
    Double_t fEnergy=0;
    Double_t AverNeut=0, AverProt=0, AverDeut=0;
    Double_t px_beam, py_beam, pz_beam, E_beam;
    if (P_beam == 3) px_beam=0., py_beam=0., pz_beam=3.53000000, E_beam=3.430721;
    else if (P_beam == 12) px_beam=0., py_beam=0., pz_beam=12.0, E_beam=12.036604;
    else if (P_beam == 159) px_beam=0., py_beam=0., pz_beam=159.0, E_beam=159.003;
    else px_beam=0., py_beam=0., pz_beam=1.0, E_beam=1.372;
    TLorentzVector Beam(px_beam, py_beam, pz_beam, E_beam);
    TLorentzVector Target(0., 0., 0, 0.938);
    TLorentzRotation rTransform;
    rTransform.Boost((Beam + Target).BoostVector()); 
            
	for (long i = 0; i < lNEvents; i++)
	{
		if (i % 100 == 0) cout<<i<<endl;
        fChain->GetEntry(i);
        fEnergy=0;
        Int_t sumZ = 0, Zb2 = 0, fNnucl = 0, fNfrag = 0, fNimf = 0, fNspect = 0;
        fNpa = fEvent->GetNpa();

        for (int j=0;j<fNpa;j++)
        {
            fParticle = fEvent->GetParticle(j);
            TLorentzVector fMomentum = fParticle->GetMomentum();
            //TLorentzVector fMomentum = (fParticle->GetMomentum()).Transform(rTransform);
            if (fParticle -> T() == 1)
            {
            if (fParticle->GetPdg()>1e9) 
            {
                hNprotons_vs_Nneutrons                  ->Fill(fParticle->GetPdg()/10%1000-fParticle->GetPdg()/10000%1000,fParticle->GetPdg()/10000%1000);
                hNfrag                                  ->Fill(fParticle->GetPdg()/10%1000);
                hNfrag_vs_ImpactParameter               ->Fill(fParticle->GetPdg()/10%1000, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(fParticle->GetPdg()/10%1000, fMomentum.E());
                hEnergyA_vs_ImpactParameter_proj        ->Fill(fEvent->GetB(), fMomentum.E()/(fParticle->GetPdg()/10%1000));
                //hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), fParticle->GetPdg()/10%1000, fMomentum.E()-E_beam*(fParticle->GetPdg()/10%1000-1));
                hNfrag_vs_Energy_scaled                 ->Fill(fParticle->GetPdg()/10%1000, fMomentum.E()-E_beam*(fParticle->GetPdg()/10%1000-1));
                hNfrag_vs_Rapidity                      ->Fill(fParticle->GetPdg()/10%1000, fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fEvent->GetB(), fMomentum.E());
                hImpactParameter_vs_Rapidity            ->Fill(fEvent->GetB(), fMomentum.Rapidity());
                hPx_vs_Py                               ->Fill(fMomentum.Px(), fMomentum.Py());
                hRapidity_vs_Energy                     ->Fill(fMomentum.E(), fMomentum.Rapidity());
                hP                                      ->Fill(fMomentum.P()/(fParticle->GetPdg()/10%1000));
                hPT_vs_A                                -> Fill (pow(fMomentum.Px()*fMomentum.Px() + fMomentum.Py()*fMomentum.Py(),0.5), fParticle->GetPdg()/10%1000);
                hRapidity                               ->Fill(fMomentum.Rapidity());
                hPseudoRapidity                         ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_vs_A                    ->Fill(fMomentum.PseudoRapidity(),fParticle->GetPdg()/10%1000);
                fEnergy=fEnergy+fMomentum.E();
                sumZ+=fParticle->GetPdg()/10000%1000;
                //if (fEvent->GetB())
                if (fParticle->GetPdg()/10000%1000 > 2) Zb2+=fParticle->GetPdg()/10000%1000;
                if (fParticle->GetPdg()/10000%1000 < 83 && fParticle->GetPdg()/10000%1000 < 209) fNfrag+=1;
                if (fParticle->GetPdg()/10000%1000 < 83 && fParticle->GetPdg()/10000%1000 < 209) fNspect+=1;
                if ( fParticle->GetPdg()/10000%1000 <= 30 && fParticle->GetPdg()/10000%1000 >=3){
                    fNimf+=1;
                    hPimf->Fill(fMomentum.P()/(fParticle->GetPdg()/10%1000));
                }
                if (fParticle->GetPdg()/10000%1000 > 30) hPheavy->Fill(fMomentum.P()/(fParticle->GetPdg()/10%1000)); 
            }       
            else if (fParticle->GetPdg()==2212) 
            {
                hNprotons_vs_Nneutrons                  ->Fill(1, 0.1);
                hNfrag                                  ->Fill(1);
                hNfrag_vs_ImpactParameter               ->Fill(1, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(1, fMomentum.E());
                hNfrag_vs_Rapidity                      ->Fill(1, fMomentum.Rapidity());
                hEnergyA_vs_ImpactParameter_proj        ->Fill(fEvent->GetB(), fMomentum.E());
                //hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), 1, fMomentum.E());
                hNfrag_vs_Energy_scaled                 ->Fill(1, fMomentum.E());
                hImpactParameter_vs_Rapidity            ->Fill(fEvent->GetB(), fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fEvent->GetB(), fMomentum.E());
                hPx_vs_Py                               ->Fill(fMomentum.Px(), fMomentum.Py());
                hRapidity_vs_Energy                     ->Fill(fMomentum.E(), fMomentum.Rapidity());                
                hP                                      ->Fill(fMomentum.P());
                hPTnucl     -> Fill (pow( fMomentum.Px()*fMomentum.Px() + fMomentum.Py()*fMomentum.Py(),0.5));
                hPT_vs_A                                -> Fill (pow(fMomentum.Px()*fMomentum.Px() + fMomentum.Py()*fMomentum.Py(),0.5), fParticle->GetPdg()/10%1000);

                hPnucl                                  ->Fill(fMomentum.P());
                hRapidity                               ->Fill(fMomentum.Rapidity());
                hPseudoRapidity                         ->Fill(fMomentum.PseudoRapidity());
                //hPseudoRapidity_nucl                    ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_vs_A                    ->Fill(fMomentum.PseudoRapidity(),1);

                fEnergy=fEnergy+fMomentum.E();
                sumZ+=1;
                fNnucl +=1;
                fNspect+=1;
            }
            else if (fParticle->GetPdg()==2112) 
            {
                hNprotons_vs_Nneutrons                  ->Fill(0.1, 1);
                hNfrag                                  ->Fill(1);
                hNfrag_vs_ImpactParameter               ->Fill(1, fEvent->GetB());
                hNfrag_vs_Energy                        ->Fill(1, fMomentum.E());
                hEnergyA_vs_ImpactParameter_proj        ->Fill(fEvent->GetB(), fMomentum.E());
                //hImpactParameter_Nfrag_vs_Energy_scaled ->Fill(fEvent->GetB(), 1, fMomentum.E());
                hNfrag_vs_Energy_scaled                 ->Fill(1, fMomentum.E());
                hNfrag_vs_Rapidity                      ->Fill(1, fMomentum.Rapidity());
                hImpactParameter_vs_Energy_Fragment     ->Fill(fEvent->GetB(), fMomentum.E());
                hImpactParameter_vs_Rapidity            ->Fill(fEvent->GetB(), fMomentum.Rapidity());
                hPx_vs_Py                               ->Fill(fMomentum.Px(), fMomentum.Py());
                hRapidity_vs_Energy                     ->Fill(fMomentum.E(), fMomentum.Rapidity());
                hP                                      ->Fill(fMomentum.P());
                hPTnucl     -> Fill (pow( fMomentum.Px()*fMomentum.Px() + fMomentum.Py()*fMomentum.Py(),0.5));
                hPT_vs_A                                -> Fill (pow(fMomentum.Px()*fMomentum.Px() + fMomentum.Py()*fMomentum.Py(),0.5), fParticle->GetPdg()/10%1000);

                hPnucl                                  ->Fill(fMomentum.P());
                hRapidity                               ->Fill(fMomentum.Rapidity());
                hPseudoRapidity                         ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_nucl                    ->Fill(fMomentum.PseudoRapidity());
                hPseudoRapidity_vs_A                    ->Fill(fMomentum.PseudoRapidity(),1);
                fEnergy=fEnergy+fMomentum.E();
                fNnucl+=1;
                fNspect+=1;
            }
            }
        }
        
        
        hNspect_vs_Espect_proj        ->Fill(fNspect, fEnergy);
        hNspect_vs_sumZ_proj          ->Fill(sumZ, fNspect);
        hNspect_vs_Zb2                ->Fill(Zb2, fNspect);
        hNnucl_vs_Nfrag_proj          ->Fill (fNfrag, fNnucl);
        hEnergyE_vs_sumZ_proj         ->Fill (sumZ, fEnergy);
        hEnergyE_vs_Nnucl_proj        ->Fill (fNnucl, fEnergy);
        hEnergyE_vs_Nimf_proj         ->Fill (fNimf, fEnergy);

        hEnergy                       ->Fill(fEnergy);
        hImpactParameter_vs_Energy    ->Fill(fEvent->GetB(), fEnergy);
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
		
	}


    TCanvas* canvas = new TCanvas("canvas");
    gPad->SetTickx();
    gPad->SetTicky();
    hImpactParameter           ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf(");

    gPad->SetLogz();

    
    hEnergyE_vs_sumZ_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hEnergyE_vs_Nnucl_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hEnergyE_vs_Nimf_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    gStyle -> SetOptStat(1000100001);
    hEnergyA_vs_ImpactParameter_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    
    hNspect_vs_sumZ_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hNspect_vs_Zb2 ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    
    
    //hNspect_vs_Espect_targ ->Draw("COLZ");
    //canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    gStyle -> SetOptStat(1000000001);
    hNspect_vs_Espect_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

        hNnucl_vs_Nfrag_proj ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

   
    
    gPad->SetLogy();
    hNfrag                               ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    gPad->SetLogy(0);
    hEnergy                              ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    gPad->SetLogy();
    hP                                   ->Draw("");
    hPnucl -> SetLineColor(2);
    hPnucl                                ->Draw("SAME");
    hPimf -> SetLineColor(4);
    hPimf                                ->Draw("SAME");
    hPheavy -> SetLineColor(6);
    hPheavy                               ->Draw("SAME");

    auto legend = new TLegend(0.75,0.61,0.96,0.96);
    legend -> AddEntry(hP, "all");
    legend -> AddEntry(hPnucl, "nucleons");
    legend -> AddEntry(hPimf, "IMF");
    legend -> AddEntry(hPheavy, "Z>30");
    //legend -> SetTextSize(0.045);
    legend -> Draw();


    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hPTnucl                            ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hRapidity                            ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hPseudoRapidity                      ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hPseudoRapidity_nucl                 ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    gPad->SetLogy(0);
        hPT_vs_A                           ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hNprotons_vs_Nneutrons               ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hNfrag_vs_ImpactParameter            ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hNfrag_vs_Energy                     ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    //hNfrag_vs_Energy_scaled              ->Draw("COLZ");
    //canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hNfrag_vs_Rapidity                   ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hImpactParameter_vs_Energy           ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hImpactParameter_vs_Energy_Fragment  ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hImpactParameter_vs_Rapidity         ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hRapidity_vs_Energy                  ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");
    hPx_vs_Py                            ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf");

    hPseudoRapidity_vs_A                 ->Draw("COLZ");
    canvas->Print("../output/histos/QA_dcmqgsm_lhc.pdf)");
       

    TFile *fout = new TFile(output_dirname, "recreate");
    
    hImpactParameter           ->Write();
   /* hNpart                     ->Write();
    hNspec                     ->Write();
    hNcoll                     ->Write();
    hNpart_vs_ImpactParameter  ->Write();
    hNspec_vs_ImpactParameter  ->Write();
    hNcoll_vs_ImpactParameter  ->Write();
    hNspec_vs_Npart            ->Write();
    hNspec_vs_Ncoll            ->Write();
    hNpart_vs_Ncoll            ->Write();*/
    //hEcc1_vs_B                 ->Write();
    //hEcc2_vs_B                 ->Write();
    //hEcc3_vs_B                 ->Write();
    //hEcc4_vs_B                 ->Write();
    
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
    hPseudoRapidity_nucl                 -> Write();
    hPseudoRapidity_vs_A                 ->Write();
    hPTnucl                              ->Write();
    hPT_vs_A                             ->Write();


    hEnergyA_vs_ImpactParameter_proj -> Write();
    //hEnergyA_vs_ImpactParameter_targ -> Write();
     hNspect_vs_Espect_proj -> Write();
   //hNspect_vs_Espect_targ -> Write();
    hNspect_vs_sumZ_proj -> Write();
    hNspect_vs_Zb2 -> Write();
   //hNspect_vs_sumZ_targ -> Write();
    //hEnergyProj_vs_EnergyTarg        -> Write();
    //hImpactParameter_vs_Energy_targ   -> Write();
   //hImpactParameter_vs_Energy_proj   -> Write();

    //hNspectT_vs_NspectP -> Write();
        //hNsumZT_vs_NsumZP -> Write();;
        hNnucl_vs_Nfrag_proj -> Write();
        //hNnucl_vs_Nfrag_targ -> Write();

        hEnergyE_vs_sumZ_proj -> Write();
        hEnergyE_vs_Nnucl_proj -> Write();
        hEnergyE_vs_Nimf_proj -> Write();
    
    //hImpactParameter_Nfrag_vs_Energy_scaled ->Write();
    
	fout->Close();
}
