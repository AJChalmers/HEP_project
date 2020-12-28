//Program to read data from root file and make invarient mass hist
// Based on: https://mhance.scipp.ucsc.edu/analysisintro.php

//V2: Taken from HEP_practice program, then fitted for CERN data

# include <iostream>
# include <fstream>
# include <math.h>
# include <iomanip>
# include <cmath>
# include <stdlib.h>
# include <cstdlib>
//# include <fstream.h>
# include <string.h>
# include <string>
//# include <dos.h> //For Sleep() 

# include "TROOT.h"
# include "TFile.h"
# include "TTree.h"
# include "TBrowser.h"
# include "TH1.h"
# include "TH2.h"
# include "TH3.h"
# include "TRandom.h"

# include <TLorentzVector.h>


using namespace std;
int main(){
	
	/////Controls/////

	//int NEntries = 469384; //469384 for data.root. Click on tree and read hist to find these.
	
	char rootFilePrefix[50] = "HEPTutorial/files/"; //Prefix to root file, this is the path to it - all the way up until its name.
	char rootFileName[50] = "data.root"; //The root file name itself - this is where the tree is.
	char TreeName1[50] = "events;1";
	
	char rootOutputFileName1[50] = "InvarientMassHist.root"; //Name of the output file for nLep hist
	
	double xMin = -100;
	double xMax = 100;
	
	const int xNChannels = 500; //Number of channels for root hist on x axis
	const int yNChannels = 500; //Number of channels for root hist on y axis 
	
	bool SaveAsRootFile = true; //If true, will save canvases to root file instead of printing to screen.
	
	///////////////////////////////////////
	
	/////Varibles (Pre-stated)/////
	
	char rootFilePath1[100]; 
	
	const int nJetArray = nJets + 1;
	
	int NJet, NMuon, NElectron, Electron_Charge, NPhoton, MCleptonPDGid, NPrimaryVertices;
	
	const int NMuonArray = NMuon + 1;
	
	int Muon_Charge[NMuonArray];
	
	float Jet_Px, Jet_Py, Jet_Pz, Jet_E, Jet_btag, Muon_Px[NMuonArray], Muon_Py[NMuonArray], Muon_Pz[NMuonArray], Muon_E[NMuonArray], Muon_Iso, 
	Electron_Px, Electron_Py, Electron_Pz, Electron_E, Electron_Iso, Photon_Px, Photon_Py, Photon_Pz, Photon_E, Photon_Iso,
	MET_px, MET_py, MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz,
	MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz,
	MClepton_px, MClepton_py, MClepton_pz, MCneutrino_px, MCneutrino_py, MCneutrino_pz, EventWeight;
	
	bool Jet_ID, triggerIsoMu24;
	
	double InvarientMass;
	
	ofstream outFile1;
	
	TH1F *histM = new TH1F("histM", "Invarient Mass hist", xNChannels, 0, 150);
	TH1F *histNMuon = new TH1F("histNMuon", "NMuon hist", 10, -5, 5);
	
	TH1F *histPosPxPrime = new TH1F("histPosPxPrime", "PosPxPrime hist", xNChannels, xMin, xMax);
	TH1F *histPosPyPrime = new TH1F("histPosPyPrime", "PosPyPrime hist", xNChannels, xMin, xMax);
	TH1F *histPosPzPrime = new TH1F("histPosPzPrime", "PosPzPrime hist", xNChannels, xMin, xMax);
	TH1F *histPosEPrime = new TH1F("histPosEPrime", "PosEPrime hist", xNChannels, xMin, xMax);
	TH1F *histNegPxPrime = new TH1F("histNegPxPrime", "NegPx hist", xNChannels, xMin, xMax);
	TH1F *histNegPyPrime = new TH1F("histNegPyPrime", "NegPy hist", xNChannels, xMin, xMax);
	TH1F *histNegPzPrime = new TH1F("histNegPzPrime", "NegPz hist", xNChannels, xMin, xMax);
	TH1F *histNegEPrime = new TH1F("histNegEPrime", "NegE hist", xNChannels, xMin, xMax);
	
	//////////////////////////////////////
	
	//Creating path to root file and getting tree//
	strcpy(rootFilePath1, rootFilePrefix);
	strcat(rootFilePath1, rootFileName);
	cout<<"Root file path 1: "<<rootFilePath1<<endl;
	
	TFile *file = new TFile(rootFilePath1);
	TTree *t1 = (TTree*)file->Get(TreeName1);
	
	// t1->SetBranchAddress("NJet",&NJet);
	// t1->SetBranchAddress("Jet_Px",&Jet_Px);
	// t1->SetBranchAddress("Jet_Py",&Jet_Py);
	// t1->SetBranchAddress("Jet_Pz",&Jet_Pz);
	// t1->SetBranchAddress("Jet_E",&Jet_E);
	// t1->SetBranchAddress("Jet_btag",&Jet_btag);
	// t1->SetBranchAddress("Jet_ID",&Jet_ID);
	t1->SetBranchAddress("NMuon",&NMuon);
	t1->SetBranchAddress("Muon_Px",&Muon_Px[0]);
	t1->SetBranchAddress("Muon_Py",&Muon_Py[0]);
	t1->SetBranchAddress("Muon_Pz",&Muon_Pz[0]);
	t1->SetBranchAddress("Muon_E",&Muon_E[0]);
	t1->SetBranchAddress("Muon_Charge",&Muon_Charge[0]);
	// t1->SetBranchAddress("Muon_Iso",&Muon_Iso);
	// t1->SetBranchAddress("NElectron",&NElectron);
	// t1->SetBranchAddress("Electron_Px",&Electron_Px);
	// t1->SetBranchAddress("Electron_Py",&Electron_Py);
	// t1->SetBranchAddress("Electron_Pz",&Electron_Pz);
	// t1->SetBranchAddress("Electron_E",&Electron_E);
	// t1->SetBranchAddress("Electron_Charge",&Electron_Charge);
	// t1->SetBranchAddress("Electron_Iso",&Electron_Iso);
	// t1->SetBranchAddress("NPhoton",&NPhoton);
	// t1->SetBranchAddress("Photon_Px",&Photon_Px);
	// t1->SetBranchAddress("Photon_Py",&Photon_Py);
	// t1->SetBranchAddress("Photon_Pz",&Photon_Pz);
	// t1->SetBranchAddress("Photon_E",&Photon_E);
	// t1->SetBranchAddress("Photon_Iso",&Photon_Iso);
	// t1->SetBranchAddress("MET_px",&MET_px);
	// t1->SetBranchAddress("MET_py",&MET_py);
	// t1->SetBranchAddress("MChadronicBottom_px",&MChadronicBottom_px);
	// t1->SetBranchAddress("MChadronicBottom_py",&MChadronicBottom_py);
	// t1->SetBranchAddress("MChadronicBottom_pz",&MChadronicBottom_pz);
	// t1->SetBranchAddress("MCleptonicBottom_px",&MCleptonicBottom_px);
	// t1->SetBranchAddress("MCleptonicBottom_py",&MCleptonicBottom_py);
	// t1->SetBranchAddress("MCleptonicBottom_pz",&MCleptonicBottom_pz);
	// t1->SetBranchAddress("MChadronicWDecayQuark_px",&MChadronicWDecayQuark_px);
	// t1->SetBranchAddress("MChadronicWDecayQuark_py",&MChadronicWDecayQuark_py);
	// t1->SetBranchAddress("MChadronicWDecayQuark_pz",&MChadronicWDecayQuark_pz);
	// t1->SetBranchAddress("MChadronicWDecayQuarkBar_px",&MChadronicWDecayQuarkBar_px);
	// t1->SetBranchAddress("MChadronicWDecayQuarkBar_py",&MChadronicWDecayQuarkBar_py);
	// t1->SetBranchAddress("MChadronicWDecayQuarkBar_pz",&MChadronicWDecayQuarkBar_pz);
	// t1->SetBranchAddress("MClepton_px",&MClepton_px);
	// t1->SetBranchAddress("MClepton_py",&MClepton_py);
	// t1->SetBranchAddress("MClepton_pz",&MClepton_pz);
	// t1->SetBranchAddress("MCleptonPDGid",&MCleptonPDGid);
	// t1->SetBranchAddress("MCneutrino_px",&MCneutrino_px);
	// t1->SetBranchAddress("MCneutrino_py",&MCneutrino_py);
	// t1->SetBranchAddress("MCneutrino_pz",&MCneutrino_pz);
	// t1->SetBranchAddress("NPrimaryVertices",&NPrimaryVertices);
	// t1->SetBranchAddress("triggerIsoMu24",&triggerIsoMu24);
	// t1->SetBranchAddress("EventWeight",&EventWeight);
	
	int NEntries1 = t1->GetEntries();
	cout<<"NEntries 1: "<<NEntries1<<endl;
	
	for(int i = 0; i < NEntries1; i++){
		t1->GetEntry(i);
		TLorentzVector pos_Muon_prime;
		TLorentzVector neg_Muon_prime;
		TLorentzVector tot_Muon_prime;
		TLorentzVector pos_Muon;
		TLorentzVector neg_Muon;
		TLorentzVector tot_Muon;
		
		if(NMuon == 2){
			cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			for(int j = 0; j < NMuon; j++){
				cout<<" i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
				if(Muon_Charge[j] == 1){
					pos_Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					histPosPxPrime->Fill(Muon_Px[j]);
					histPosPyPrime->Fill(Muon_Py[j]);
					histPosPzPrime->Fill(Muon_Pz[j]);
					histPosEPrime->Fill(Muon_E[j]);
				}
				if(Muon_Charge[j] == -1){
					neg_Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					histNegPxPrime->Fill(Muon_Px[j]);
					histNegPyPrime->Fill(Muon_Py[j]);
					histNegPzPrime->Fill(Muon_Pz[j]);
					histNegEPrime->Fill(Muon_E[j]);
				}
			}
			if((Muon_Charge[0] == 1 && Muon_Charge[1] == -1) || (Muon_Charge[0] == -1 && Muon_Charge[1] == 1)){
				tot_Muon_prime = pos_Muon_prime + neg_Muon_prime;
				histM->Fill(tot_Muon_prime.M());
			}
		}
	}
	
	//////Graphics//////
	cout<<"Making histrograms..."<<endl;
	
	TCanvas *c1 = new TCanvas("c1","InvarientMass Histrogram" ,200,10,900,700);
		c1->SetFillColor(10);
		c1->SetGrid();
		c1->GetFrame()->SetFillColor(10);
		c1->GetFrame()->SetBorderSize(12);
		c1->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c1->cd();
			histM->SetXTitle("Invarient Mass");
			histM->SetYTitle("Counts");
			histM->SetStats(1);
			histM->Draw();
			if(SaveAsRootFile){
				TFile outfile1(rootOutputFileName1, "RECREATE");
				c1->Write(rootOutputFileName1);
				outfile1.Close();
			}
			

		

	return 0;
}