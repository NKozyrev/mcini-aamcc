#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <string>
//#include "G4Types.hh"

#include "../include/EventInitialState.h"
#include "../include/UEvent.h"
#include "../include/URun.h"

R__LOAD_LIBRARY(libMcIniData.so)

void convertAAMCC(TString inputFileName = "particles.root", TString outputFileName = "mcini_aamcc.root", Bool_t only_not_empty = true,
	Int_t aProj = 197, Int_t zProj = 79, Double_t pProj = 12.,
	Int_t aTarg = 197, Int_t zTarg = 79, Double_t pTarg = 0.,
	Double_t bMin = 0., Double_t bMax = 20., Int_t bWeight = 0,
	Double_t phiMin = 0., Double_t phiMax = 0.,
	Double_t sigma = 0., Int_t nEvents = 0)
{
	TFile* fIn = new TFile(inputFileName);
	TTree* fTreeG = (TTree*)fIn->Get("Glauber");
	TTree* fTreeC = (TTree*)fIn->Get("Conditions");
	if (!fTreeG || !fTreeC)
	{
		std::cerr << "TTree was not found in the input file." << std::endl;
		return;
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
	Int_t Npart=0, Ncoll;
	Float_t b;
	UInt_t id;
	Double_t AonA, AonB, ZonA, ZonB, Energy;
	Float_t Ecc[10];


	fTreeG->SetBranchAddress("id", &id);
	fTreeG->SetBranchAddress("A_on_A", &MassOnSideA);
	fTreeG->SetBranchAddress("A_on_B", &MassOnSideB);
	fTreeG->SetBranchAddress("Z_on_A", &ChargeOnSideA);
	fTreeG->SetBranchAddress("Z_on_B", &ChargeOnSideB);
	fTreeG->SetBranchAddress("impact_parameter", &b);
	//fTreeG->SetBranchAddress("NpartA", &Npart);
	fTreeG->SetBranchAddress("Ncoll", &Ncoll);
	fTreeG->SetBranchAddress("pX_on_A", &pXonSideA);
	fTreeG->SetBranchAddress("pY_on_A", &pYonSideA);
	fTreeG->SetBranchAddress("pZ_on_A", &pZonSideA);
	fTreeG->SetBranchAddress("pX_on_B", &pXonSideB);
	fTreeG->SetBranchAddress("pY_on_B", &pYonSideB);
	fTreeG->SetBranchAddress("pZ_on_B", &pZonSideB);
	fTreeG->SetBranchAddress("Ecc", &Ecc);

	fTreeC->SetBranchAddress("Xsect_total", &sigma);
	fTreeC->SetBranchAddress("Mass_on_A", &AonA);
	fTreeC->SetBranchAddress("Mass_on_B", &AonB);
	fTreeC->SetBranchAddress("Charge_on_A", &ZonA);
	fTreeC->SetBranchAddress("Charge_on_B", &ZonB);
	fTreeC->SetBranchAddress("Kinetic_energy_per_nucleon_of_projectile_in_GeV", &pProj);

	fTreeC->GetEntry(0);

	aProj = AonA; zProj = ZonA; aTarg = AonB; zTarg = ZonB;

	TFile* fOut = new TFile(outputFileName, "recreate");

	TTree* iniTree = new TTree("events", "AAMCC");
	URun header("AAMCC", "final state only", aProj, zProj, pProj, aTarg, zTarg, pTarg, bMin, bMax, bWeight, phiMin, phiMax, sigma, nEvents);
	UEvent* event = new UEvent;
	EventInitialState* iniState = new EventInitialState;
	iniTree->Branch("iniState", "EventInitialState", iniState);
	iniTree->Branch("event", "UEvent", event);

	Long64_t eventCounter = 0;
	Int_t child[2] = { 0,0 };

	for (Int_t iev = 0; iev < fTreeG->GetEntries(); iev++)
	{
		fTreeG->GetEntry(iev);
		event->Clear();
		iniState->Clear();

		if (iev % 100 == 0) std::cout << "Event [" << iev << "/" << fTreeG->GetEntries() << "]" << std::endl;

		// Fill event
		event->SetEventNr(id);
		event->SetB(b);
		event->SetEcc(Ecc);
		event->SetPhi(0.);
		event->SetNes(1);
		event->SetStepNr(0);
		event->SetStepT(0);

		iniState->setNColl(Ncoll);
		iniState->setNPart(Npart);
		// Fill particle
		for (Int_t ipart = 0; ipart < (MassOnSideA->size()); ipart++)
		{
			Int_t fragment_id = 0;
			
			Int_t hundreds_mass = MassOnSideA->at(ipart)/100;
			Int_t dozens_mass = MassOnSideA->at(ipart)/10;
			Int_t hundreds_charge = ChargeOnSideA->at(ipart)/100;
			Int_t dozens_charge = ChargeOnSideA->at(ipart)/10;

			fragment_id = pow(10, 9) + pow(10, 6)*hundreds_charge + pow(10, 5)*(dozens_charge - hundreds_charge*10) + pow(10, 4)*(ChargeOnSideA->at(ipart)- dozens_charge*10) + pow(10, 3)*hundreds_mass + pow(10, 2)*(dozens_mass - hundreds_mass*10) + 10*(MassOnSideA->at(ipart)- dozens_mass*10);

			Energy = pow(pow(pXonSideA->at(ipart) / 1000, 2) + pow(pYonSideA->at(ipart) / 1000, 2) + pow(pZonSideA->at(ipart) / 1000, 2) + pow(0.9395654 * (MassOnSideA->at(ipart) - ChargeOnSideA->at(ipart)) + 0.9382721 * ChargeOnSideA->at(ipart), 2), 0.5);
			
			if (MassOnSideA->at(ipart) == 1 && ChargeOnSideA->at(ipart) == 0) event->AddParticle(ipart, 2112, 0,
				0, 0,
				0, 0, child,
				pXonSideA->at(ipart) / 1000, pYonSideA->at(ipart) / 1000, pZonSideA->at(ipart) / 1000, Energy,
				0, 0, 0, 1,
				1.);
			else if (MassOnSideA->at(ipart) == 1 && ChargeOnSideA->at(ipart) == 1) event->AddParticle(ipart, 2212, 0,
				0, 0,
				0, 0, child,
				pXonSideA->at(ipart) / 1000, pYonSideA->at(ipart) / 1000, pZonSideA->at(ipart) / 1000, Energy,
				0, 0, 0, 1,
				1.);
			else{
				event->AddParticle(ipart, fragment_id, 0,
				0, 0,
				0, 0, child,
				pXonSideA->at(ipart) / 1000, pYonSideA->at(ipart) / 1000, pZonSideA->at(ipart) / 1000, Energy,
				0, 0, 0, 1,
				1.);
			}
		}

		for (Int_t ipart = 0; ipart < (MassOnSideB->size()); ipart++)
		{
			Int_t fragment_id = 0;
			
			Int_t hundreds_mass = MassOnSideB->at(ipart)/100;
			Int_t dozens_mass = MassOnSideB->at(ipart)/10;
			Int_t hundreds_charge = ChargeOnSideB->at(ipart)/100;
			Int_t dozens_charge = ChargeOnSideB->at(ipart)/10;

			fragment_id = pow(10, 9) + pow(10, 6)*hundreds_charge + pow(10, 5)*(dozens_charge - hundreds_charge*10) + pow(10, 4)*(ChargeOnSideB->at(ipart)- dozens_charge*10) + pow(10, 3)*hundreds_mass + pow(10, 2)*(dozens_mass - hundreds_mass*10) + 10*(MassOnSideB->at(ipart)- dozens_mass*10);

			Energy = pow(pow(pXonSideB->at(ipart) / 1000, 2) + pow(pYonSideB->at(ipart) / 1000, 2) + pow(pZonSideB->at(ipart) / 1000, 2) + pow(0.9395654 * (MassOnSideB->at(ipart) - ChargeOnSideB->at(ipart)) + 0.9382721 * ChargeOnSideB->at(ipart), 2), 0.5);
			
			if (MassOnSideB->at(ipart) == 1 && ChargeOnSideB->at(ipart) == 0) event->AddParticle(ipart, 2112, 0,
				0, 0,
				0, 0, child,
				pXonSideB->at(ipart) / 1000, pYonSideB->at(ipart) / 1000, pZonSideB->at(ipart) / 1000, Energy,
				0, 0, 0, -1,
				1.);
			else if (MassOnSideB->at(ipart) == 1 && ChargeOnSideB->at(ipart) == 1) event->AddParticle(ipart, 2212, 0,
				0, 0,
				0, 0, child,
				pXonSideB->at(ipart) / 1000, pYonSideB->at(ipart) / 1000, pZonSideB->at(ipart) / 1000, Energy,
				0, 0, 0, -1,
				1.);
			else{
				event->AddParticle(ipart, fragment_id, 0,
				0, 0,
				0, 0, child,
				pXonSideB->at(ipart) / 1000, pYonSideB->at(ipart) / 1000, pZonSideB->at(ipart) / 1000, Energy,
				0, 0, 0, -1,
				1.);
			}
		}
		iniTree->Fill();
		eventCounter++;
	}
	header.SetNEvents(eventCounter);

	fOut->cd();
	header.Write();
	iniTree->Write();
	fOut->Close();
	gApplication->Terminate();
}