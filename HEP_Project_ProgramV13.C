//Program to read data from root file and make invarient mass hist
// Based on: https://mhance.scipp.ucsc.edu/analysisintro.php

//V2: Taken from HEP_practice program, then fitted for CERN data. Got Invarient mass for Z with muons

//V3: Added invarient mass for W with electrons using transverse mass

//V4: Changes made with Nick, added Iso bounds and UseAllMuonsForW bool to get only pos muons for W+. Adding more hists for mT and MET

//V5: Fixing MET_phi

//V6: Added more iso hists

//V7: Adding electron W's back. 
///Note: Add WElectrons to WMuons or keep them seperate (currently seperate)

//V8: Added bJet with Nick, not working correctly? Keeps crashing when added bJet with WMuon

//V9: Fixed V8 so it will run

//V10: Made NMuon and NJet = 2 case allowed. This is not the whole picture - Muons and Jets are not parallel arrays. Look in notebook to make V11

//V11: 

//V12: Changes from 11/10

//V13: Huge overhall of WMuon and Top loop. They are now split into two loops

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

//gDebug = 1;
using namespace std;
int main(){
	
	/////Controls/////
	
	char rootFilePrefix[50] = "HEPTutorial/files/"; //Prefix to root file, this is the path to it - all the way up until its name.
	char rootFileName[50] = "data.root"; //The root file name itself - this is where the tree is.
	char TreeName1[50] = "events;1";
	
	char rootOutputFileName1[50] = "NoBounds_InvarientMassHist.root"; //Name of an output root file 
	char rootOutputFileName2[50] = "NoBounds_PosVNegPrime.root"; //Name of an output root file
	char rootOutputFileName3[50] = "NoBounds_PosVNeg.root"; //Name of an output root file
	char rootOutputFileName4[50] = "NoBounds_MuonIso_btag.root"; //Name of an output root file
	char rootOutputFileName5[50] = "NoBounds_MET_WMuonMass.root";
	char rootOutputFileName6[50] = "NoBounds_TransverseMassHist.root";
	char rootOutputFileName7[50] = "NoBounds_InvarientMassHist_TOP.root";
	char rootOutputFileName8[50] = "NoBounds_InvarientMassHist_W.root";
	
	double MuonIsoBound_ForZ = 1000000000; //Upper bound to muon iso when looking for ZMuons 
	double MuonIsoBound_ForW = 1000000000; //Upper bound to muon iso when looking for WMuons 
	double ElectronIsoBound_ForW = 1000000000; //Upper bound to electron iso when looking for WElectrons 
	
	double RelIsoBound = 0.05; //Iso/Pt bound
	
	double Jet_btag_Bound = 1.74; //Lower bound to Jet_btag
	
	double MuonPtCut = 25;
	
	double mT_MaxCount = 1000; //Max counts for mT hists. About 1,000 for full data set and 25 for 1/44
	
	double xMin = -100;
	double xMax = 100;
	double yMin = -300;
	double yMax = 300;
	
	const int xNChannels = 1000; //Number of channels for root hist on x axis
	const int yNChannels = 500; //Number of channels for root hist on y axis 
	
	bool SaveAsRootFile = true; //If true, will save canvases to root file instead of printing to screen.
	bool FitInvarientMassHist = true; //If true, will fit Invarient mass hist for Z from muons with gaus, and W from muons with gaus
	// bool UseAllMuonsForW = false; //If true, will use all Muons for W and not just positive ones
	// bool UseAllMuonsForW = false; //If true, will use all Muons for W and not just positive ones
	
	///////////////////////////////////////
	
	/////Varibles (Pre-stated)/////
	
	char rootFilePath1[100]; 
	
	int NJet, NMuon, NElectron, NPhoton, MCleptonPDGid, NPrimaryVertices;
	
	const int NMuonArray = 20;
	const int NElectronArray = 20;
	const int NJetArray = 20;
	
	int Muon_Charge[NMuonArray], Electron_Charge[NElectronArray];
	
	float Jet_Px[NJetArray], Jet_Py[NJetArray], Jet_Pz[NJetArray], Jet_E[NJetArray], Jet_btag[NJetArray], Muon_Px[NMuonArray], Muon_Py[NMuonArray], Muon_Pz[NMuonArray], Muon_E[NMuonArray], Muon_Iso[NMuonArray], 
	Electron_Px[NElectronArray], Electron_Py[NElectronArray], Electron_Pz[NElectronArray], Electron_E[NElectronArray], Electron_Iso[NElectronArray], Photon_Px, Photon_Py, Photon_Pz, Photon_E, Photon_Iso,
	MET_px, MET_py, MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz,
	MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz,
	MClepton_px, MClepton_py, MClepton_pz, MCneutrino_px, MCneutrino_py, MCneutrino_pz, EventWeight;
	
	bool Jet_ID[NJetArray], triggerIsoMu24;
	
	double InvarientMass;
	
	ofstream outFile1;
	
	TH1F *histM_Muon = new TH1F("histM_Muon", "Invarient Mass hist: Muon", 60 , 60, 120);
	TH1F *histM_WMuon = new TH1F("histM_WMuon", "Transverse Mass hist: WMuon", xNChannels, 0, 150);
	TH1F *histM_WMuon2 = new TH1F("histM_WMuon2", "Transverse Mass hist2: WMuon", xNChannels, 0, 150);
	TH1F *histM_WMuon3 = new TH1F("histM_WMuon3", "Transverse Mass hist3: WMuon", xNChannels, 0, 150);
	
	TH1F *histM_WMuonPos = new TH1F("histM_WMuonPos", "Transverse Mass hist: WMuon Pos", xNChannels, 0, 150);
	TH1F *histM_WMuonNeg = new TH1F("histM_WMuonNeg", "Transverse Mass hist: WMuon Neg", xNChannels, 0, 150);
	
	TH1F *histM_WElectron = new TH1F("histM_WElectron", "Transverse Mass hist: WElectron", xNChannels, 0, 150);
	TH1F *histM_WElectronPos = new TH1F("histM_WElectronPos", "Transverse Mass hist: WElectron Pos", xNChannels, 0, 150);
	TH1F *histM_WElectronNeg = new TH1F("histM_WElectronNeg", "Transverse Mass hist: WElectron Neg", xNChannels, 0, 150);
	
	TH1F *histNMuon = new TH1F("histNMuon", "NMuon hist", 10, -5, 5);
	TH1F *histNElectron = new TH1F("histNElectron", "NElectron hist", 10, -5, 5);
	
	TH1F *histPosPxPrime = new TH1F("histPosPxPrime", "PosPxPrime hist", xNChannels, xMin, xMax);
	TH1F *histPosPyPrime = new TH1F("histPosPyPrime", "PosPyPrime hist", xNChannels, xMin, xMax);
	TH1F *histPosPzPrime = new TH1F("histPosPzPrime", "PosPzPrime hist", xNChannels, xMin, xMax);
	TH1F *histPosEPrime = new TH1F("histPosEPrime", "PosEPrime hist", xNChannels, xMin, xMax);
	TH1F *histNegPxPrime = new TH1F("histNegPxPrime", "NegPx hist", xNChannels, xMin, xMax);
	TH1F *histNegPyPrime = new TH1F("histNegPyPrime", "NegPy hist", xNChannels, xMin, xMax);
	TH1F *histNegPzPrime = new TH1F("histNegPzPrime", "NegPz hist", xNChannels, xMin, xMax);
	TH1F *histNegEPrime = new TH1F("histNegEPrime", "NegE hist", xNChannels, xMin, xMax);
	
	TH2F *histPosVNegPx = new TH2F("histPosVNegPx", "PosVNegPx hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH2F *histPosVNegPy = new TH2F("histPosVNegPy", "PosVNegPy hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH2F *histPosVNegPz = new TH2F("histPosVNegPz", "PosVNegPz hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH2F *histPosVNegE = new TH2F("histPosVNegE", "PosVNegE hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	
	TH2F *histPosVNegPxPrime = new TH2F("histPosVNegPxPrime", "PosVNegPxPrime hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH2F *histPosVNegPyPrime = new TH2F("histPosVNegPyPrime", "PosVNegPyPrime hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH2F *histPosVNegPzPrime = new TH2F("histPosVNegPzPrime", "PosVNegPzPrime hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH2F *histPosVNegEPrime = new TH2F("histPosVNegEPrime", "PosVNegEPrime hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	
	TH1F *histMuonIso = new TH1F("histMuonIso", "MuonIso hist", xNChannels, 0, 20);
	TH1F *histMuonPt = new TH1F("histMuonPt", "MuonPt hist", xNChannels, xMin, xMax);
	TH1F *histWMuonIso = new TH1F("histMuonIso", "MuonIso hist", xNChannels, 0, 20);
	TH1F *histWElectronIso = new TH1F("histElectronIso", "ElectronIso hist", xNChannels, 0, 20);
	
	TH2F *hist_WMuonVIso = new TH2F("hist_WMuonVIso", "WMuonVIso hist", xNChannels, 0, 20, 0, 1000);
	TH2F *hist_bJetVbtag = new TH2F("hist_bJetVbtag", "bJetVbtag hist", xNChannels, 0, 20, 0, 1000);
	 
	TH2F *histMETxy = new TH2F("histMETxy", "METxy hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH1F *histMET = new TH1F("histMET", "MET hist", xNChannels, xMin, xMax);
	TH1F *histMET_Phi = new TH1F("histMET_Phi", "MET_Phi hist", xNChannels, 0, 360);

	TH1F *histM_TopPos = new TH1F("histM_TopPos", "Invarient Mass hist: TopPos", 50, 0, 500);
	TH1F *histM_TopNeg = new TH1F("histM_TopNeg", "Invarient Mass hist: TopNeg", 50, 0, 500);
	TH1F *histM_Top = new TH1F("histM_Top", "Invarient Mass hist: Top", 50, 0, 500);
	
	TH1F *histM_WPos = new TH1F("histM_WPos", "Invarient Mass hist: WPos", xNChannels, 0, 200);
	TH1F *histM_WNeg = new TH1F("histM_WNeg", "Invarient Mass hist: WNeg", xNChannels, 0, 200);
	TH1F *histM_W = new TH1F("histM_W", "Invarient Mass hist: W", xNChannels, 0, 200);
	
	double Pi = acos(-1.);
	double Deg2Rad = Pi/180.;
	
	//////////////////////////////////////
	
	//Creating path to root file and getting tree//
	strcpy(rootFilePath1, rootFilePrefix);
	strcat(rootFilePath1, rootFileName);
	cout<<"Root file path 1: "<<rootFilePath1<<endl;
	
	TFile *file = new TFile(rootFilePath1);
	TTree *t1 = (TTree*)file->Get(TreeName1);
	
	t1->SetBranchAddress("NJet",&NJet);
	t1->SetBranchAddress("Jet_Px",&Jet_Px[0]);
	t1->SetBranchAddress("Jet_Py",&Jet_Py[0]);
	t1->SetBranchAddress("Jet_Pz",&Jet_Pz[0]);
	t1->SetBranchAddress("Jet_E",&Jet_E[0]);
	t1->SetBranchAddress("Jet_btag",&Jet_btag[0]);
	t1->SetBranchAddress("Jet_ID",&Jet_ID[0]);
	t1->SetBranchAddress("NMuon",&NMuon);
	t1->SetBranchAddress("Muon_Px",&Muon_Px[0]);
	t1->SetBranchAddress("Muon_Py",&Muon_Py[0]);
	t1->SetBranchAddress("Muon_Pz",&Muon_Pz[0]);
	t1->SetBranchAddress("Muon_E",&Muon_E[0]);
	t1->SetBranchAddress("Muon_Charge",&Muon_Charge[0]);
	t1->SetBranchAddress("Muon_Iso",&Muon_Iso[0]);
	t1->SetBranchAddress("NElectron",&NElectron);
	t1->SetBranchAddress("Electron_Px",&Electron_Px[0]);
	t1->SetBranchAddress("Electron_Py",&Electron_Py[0]);
	t1->SetBranchAddress("Electron_Pz",&Electron_Pz[0]);
	t1->SetBranchAddress("Electron_E",&Electron_E[0]);
	t1->SetBranchAddress("Electron_Charge",&Electron_Charge[0]);
	t1->SetBranchAddress("Electron_Iso",&Electron_Iso[0]);
	// t1->SetBranchAddress("NPhoton",&NPhoton);
	// t1->SetBranchAddress("Photon_Px",&Photon_Px);
	// t1->SetBranchAddress("Photon_Py",&Photon_Py);
	// t1->SetBranchAddress("Photon_Pz",&Photon_Pz);
	// t1->SetBranchAddress("Photon_E",&Photon_E);
	// t1->SetBranchAddress("Photon_Iso",&Photon_Iso);
	t1->SetBranchAddress("MET_px",&MET_px);
	t1->SetBranchAddress("MET_py",&MET_py);
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
	t1->SetBranchAddress("triggerIsoMu24",&triggerIsoMu24);
	// t1->SetBranchAddress("EventWeight",&EventWeight);
	
		int NWMuon = 0;
		int NbJet = 0;
	
		TLorentzVector pos_Muon_prime;
		TLorentzVector neg_Muon_prime;
		TLorentzVector tot_Muon_prime;
		TLorentzVector pos_Muon;
		TLorentzVector neg_Muon;
		TLorentzVector tot_Muon;
		
		TLorentzVector WMuon_prime;
		TLorentzVector WMuon_prime_Pos;
		TLorentzVector WMuon_prime_Neg;

		TLorentzVector WElectron_prime;
		TLorentzVector WElectron_prime_Pos;
		TLorentzVector WElectron_prime_Neg;
		
		TLorentzVector bJet_Pos;
		TLorentzVector bJet_Neg;
		TLorentzVector TopPos;
		TLorentzVector TopNeg;
		
		TLorentzVector MuonTemp;
		
		TLorentzVector MET;
		
		TLorentzVector Neutrino1_Pos;
		TLorentzVector Neutrino1_Neg;
		TLorentzVector Neutrino2_Pos;
		TLorentzVector Neutrino2_Neg;
		
		TLorentzVector WPos;
		TLorentzVector WNeg;
		
		TLorentzVector Muon_prime;
		TLorentzVector Muon_prime_Iso1;
		TLorentzVector Muon_prime_Iso2;
		
		TLorentzVector bJet;
		
		double Muon_Pt_pos;
		double Muon_Pt_neg;
		
		double MET_pz;
		double MET_E;
		
		double RelIso;
	
	int NEntries1 = t1->GetEntries();
	cout<<"NEntries 1: "<<NEntries1<<endl;
	
	for(int i = 0; i < NEntries1; i++){
		cout<<"*************"<<endl;
		cout<<"i: "<<i<<endl;
		t1->GetEntry(i);
		
		/*pos_Muon_prime = new TLorentzVector;
		neg_Muon_prime = new TLorentzVector;
		tot_Muon_prime = new TLorentzVector;
		pos_Muon = new TLorentzVector;
		neg_Muon = new TLorentzVector;
		tot_Muon = new TLorentzVector;
		
		WMuon_prime = new TLorentzVector;
		WMuon_prime_Pos = new TLorentzVector;
		WMuon_prime_Neg = new TLorentzVector;

		WElectron_prime = new TLorentzVector;
		WElectron_prime_Pos = new TLorentzVector;
		WElectron_prime_Neg = new TLorentzVector;
		
		bJet_prime_Pos = new TLorentzVector;
		bJet_prime_Neg = new TLorentzVector;
		TopPos = new TLorentzVector;
		TopNeg = new TLorentzVector; 
		
		MuonTemp = new TLorentzVector;
		
		MET = new TLorentzVector; */
		
		bool WPosHappened = false;
		bool WNegHappened = false;

		//histMuonIso->Fill(Muon_Iso); 
		
		
		//For the Z//
		int N_IsoMuon = 0;
		if(NMuon >= 2){
			cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			for(int j = 0; j < NMuon; j++){
				histMuonIso->Fill(Muon_Iso[j]); 
				MuonTemp.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
				RelIso = Muon_Iso[j]/MuonTemp.Pt();
				if((RelIso < RelIsoBound) && (MuonTemp.P() > 0)){
					++N_IsoMuon;
					cout<<"Z| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					//if(Muon_Charge[j] == 1){
					if (N_IsoMuon == 1){
						pos_Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
						// histPosPxPrime->Fill(Muon_Px[j]);
						// histPosPyPrime->Fill(Muon_Py[j]);
						// histPosPzPrime->Fill(Muon_Pz[j]);
						// histPosEPrime->Fill(Muon_E[j]);
						double Muon_Iso_pos = Muon_Iso[j];
						Muon_Pt_pos = pos_Muon_prime.Pt();
					}
					//if(Muon_Charge[j] == -1){
					if (N_IsoMuon == 2) {
						neg_Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
						// histNegPxPrime->Fill(Muon_Px[j]);
						// histNegPyPrime->Fill(Muon_Py[j]);
						// histNegPzPrime->Fill(Muon_Pz[j]);
						// histNegEPrime->Fill(Muon_E[j]);
						double Muon_Iso_neg = Muon_Iso[j];
						double Muon_Pt_neg = neg_Muon_prime.Pt();
					}
					histPosVNegPxPrime->Fill(pos_Muon_prime.Px(), neg_Muon_prime.Px());
					histPosVNegPyPrime->Fill(pos_Muon_prime.Py(), neg_Muon_prime.Py());
					histPosVNegPzPrime->Fill(pos_Muon_prime.Pz(), neg_Muon_prime.Pz());
					histPosVNegEPrime->Fill(pos_Muon_prime.E(), neg_Muon_prime.E());
					
				}
				
			}
			//if(((Muon_Charge[0] == 1 && Muon_Charge[1] == -1) || (Muon_Charge[0] == -1 && Muon_Charge[1] == 1)) && ((Muon_Iso[0] < MuonIsoBound_ForZ) && (Muon_Iso[1] < MuonIsoBound_ForZ)) && (MuonTemp.P() > 0)){
				tot_Muon_prime = pos_Muon_prime + neg_Muon_prime;
				histM_Muon->Fill(tot_Muon_prime.M());
				
				//histMuonIso->Fill(Muon_Iso_pos);  ///neg_Muon_prime.Pt());
				
				histMuonPt->Fill(Muon_Pt_pos);
				histMuonPt->Fill(Muon_Pt_neg);
				
				TVector3 tot_Muon_boost;
				tot_Muon_boost = tot_Muon_prime.BoostVector();
				TLorentzRotation l;
				TLorentzRotation l2;
				l.Boost(0,0,0);
				l.Boost(tot_Muon_boost);
				l2 = l.Inverse();
				neg_Muon = l2*neg_Muon_prime;
				pos_Muon = l2*pos_Muon_prime;
				tot_Muon = pos_Muon + neg_Muon;
				
				histPosVNegPx->Fill(pos_Muon.Px(), neg_Muon.Px());
				histPosVNegPy->Fill(pos_Muon.Px(), neg_Muon.Px());
				histPosVNegPz->Fill(pos_Muon.Px(), neg_Muon.Px());
				histPosVNegE->Fill(pos_Muon.Px(), neg_Muon.Px());
			//}
		}
		
		int NIsoMuon = 0;
		for(int j = 0; j < NMuon; j++){ // Muon loop: Looks at every loop
			cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			
			cout<<"W| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
			cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
			Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
			RelIso = Muon_Iso[j]/Muon_prime.Pt();
			
			if(RelIso < RelIsoBound){
				histWMuonIso->Fill(Muon_Iso[j]);
				
				if(NIsoMuon == 0){
					cout<<"Muon| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
					Muon_prime_Iso1.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					
				}
				//This one is for future implementation of the Z
				if(NIsoMuon == 1){
					cout<<"Muon| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
					Muon_prime_Iso2.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
				}
				NIsoMuon++;
			}
		}
		
		int NBJet = 0;
		bool IsSelected = false;
		if(NIsoMuon == 0){
			if(Muon_prime_Iso1.Pt() > MuonPtCut){
				NBJet = 0;
				for(int jet = 0; jet < NJet; jet++){
					if(Jet_ID[jet] && (Jet_btag[jet] > Jet_btag_Bound)){
						NBJet++;
					}
				}
				if(NBJet > 1){
					if(triggerIsoMu24){
						IsSelected = true;
						cout<<"IS SELECTED SET TO TRUE!!!!!"<<endl;
					}
				}
			}
		}
		
		if(IsSelected){
			// Jet loop for to get W and then Top
			for(int jet = 0; jet < NJet; jet++){  
				if(Jet_ID[jet] && (Jet_btag[jet] > Jet_btag_Bound)){
					//For the W and Top//
					MET.SetXYZM(MET_px, MET_py, 0., 0.);
					bJet.SetPxPyPzE(Jet_Px[jet], Jet_Py[jet], Jet_Pz[jet], Jet_E[jet]); 
					
					double mW = 80.4;
					double A = pow(Muon_prime_Iso1.E(), 2) - pow(Muon_prime_Iso1.Pz(), 2);
					double B = mW * mW / 2. + (Muon_prime_Iso1.Px() * MET_px + Muon_prime_Iso1.Py() * MET_py);
					double D = pow(Muon_prime_Iso1.E(), 2) * (pow(B, 2) - pow(MET.Pt(), 2) * A);
					
					// first solution
					if(D >= 0){
					  double MET_pz = Muon_prime_Iso1.Pz() * B / A + sqrt(D) / A;
					  double MET_E = sqrt(MET_px * MET_px + MET_py * MET_py + MET_pz * MET_pz);
					  Neutrino1_Pos.SetPxPyPzE(MET_px, MET_py, MET_pz, MET_E);
					  WPos = Neutrino1_Pos + Muon_prime_Iso1;
					  histM_WPos->Fill(WPos.M());
					  histM_W->Fill(WPos.M());
					  
					  TopPos = WPos + bJet;
					  histM_TopPos->Fill(TopPos.M());
					  histM_Top->Fill(TopPos.M());
					  
					}
					// second solution
					if(D > 0){
					  double pz =  Muon_prime_Iso1.Pz() * B / A - sqrt(D) / A;
					  double E = sqrt(MET_px * MET_px + MET_py * MET_py + MET_pz * MET_pz);
					  Neutrino2_Pos.SetPxPyPzE(MET_px, MET_py, MET_pz, MET_E);
					  WPos = Neutrino2_Pos + Muon_prime_Iso1;
					  histM_WPos->Fill(WPos.M());
					  histM_W->Fill(WPos.M());
					  
					  TopPos = WPos + bJet;
					  histM_TopPos->Fill(TopPos.M());
					  histM_Top->Fill(TopPos.M());
					}
				}
			}
		}
	}
	
	// if(FitInvarientMassHist){
		TF1 *f1 = new TF1("f1","gaus",80,100);
		TF1 *f2 = new TF1("f2","([0]*exp(-((x-[1])^2)/(2*[2]^2)))",60,85);
		
		f2->SetParameters(0, 200); //Normal
		f2->SetParameters(1, 70); //Centroid
		f2->SetParameters(2, 25); //Width
	// }
	
	//////Graphics//////
	cout<<"Making histrograms..."<<endl;
	
	TCanvas *c1 = new TCanvas("c1","InvarientMass Histrogram" ,200,10,900,700);
		c1->Divide(2,2);
		c1->SetFillColor(10);
		c1->SetGrid();
		c1->GetFrame()->SetFillColor(10);
		c1->GetFrame()->SetBorderSize(12);
		c1->Range(0,0,1,1);

					
		gStyle->SetOptStat(1);
		if(FitInvarientMassHist){
			gStyle->SetOptFit(1);
		}
			
			c1->cd(1);
			gPad->SetLogy();
			histM_Muon->SetXTitle("Invarient Mass");
			histM_Muon->SetYTitle("Counts");
			histM_Muon->SetStats(1);
			if(FitInvarientMassHist){
				histM_Muon->Fit(f1, "M");
			}
			histM_Muon->Draw();
			
			c1->cd(2);
			histM_WMuon->SetXTitle("Mt Mass");
			histM_WMuon->SetYTitle("Counts");
			histM_WMuon->SetStats(1);
			if(FitInvarientMassHist){
				histM_WMuon->Fit(f2, "M");
			}
			histM_WMuon->Draw();
			
			c1->cd(3);
			histM_WMuon2->SetXTitle("Mt 2 Mass");
			histM_WMuon2->SetYTitle("Counts");
			histM_WMuon2->SetStats(1);
			histM_WMuon2->Draw();
			
			c1->cd(4);
			histM_WMuon3->SetXTitle("Mt 3 Mass");
			histM_WMuon3->SetYTitle("Counts");
			histM_WMuon3->SetStats(1);
			histM_WMuon3->Draw();
			
			if(SaveAsRootFile){
				TFile outfile1(rootOutputFileName1, "RECREATE");
				c1->Write(rootOutputFileName1);
				outfile1.Close();
			}
		
	TCanvas *c2 = new TCanvas("c2","PosVNegPrime Histrogram" ,200,10,900,700);
		c2->Divide(2,2);
		c2->SetFillColor(10);
		c2->SetGrid();
		c2->GetFrame()->SetFillColor(10);
		c2->GetFrame()->SetBorderSize(12);
		c2->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c2->cd(1);
			histPosVNegPxPrime->SetXTitle("PosPxPrime");
			histPosVNegPxPrime->SetYTitle("NegPxPrime");
			histPosVNegPxPrime->SetStats(1);
			histPosVNegPxPrime->Draw();
			
			c2->cd(2);
			histPosVNegPyPrime->SetXTitle("PosPyPrime");
			histPosVNegPyPrime->SetYTitle("NegPxPrime");
			histPosVNegPyPrime->SetStats(1);
			histPosVNegPyPrime->Draw();
			
			c2->cd(3);
			histPosVNegPzPrime->SetXTitle("PosPzPrime");
			histPosVNegPzPrime->SetYTitle("NegPxPrime");
			histPosVNegPzPrime->SetStats(1);
			histPosVNegPzPrime->Draw();
			
			c2->cd(4);
			histPosVNegEPrime->SetXTitle("PosEPrime");
			histPosVNegEPrime->SetYTitle("NegPxPrime");
			histPosVNegEPrime->SetStats(1);
			histPosVNegEPrime->Draw();
			
			if(SaveAsRootFile){
				TFile outfile2(rootOutputFileName2, "RECREATE");
				c2->Write(rootOutputFileName2);
				outfile2.Close();
			}
	
	TCanvas *c3 = new TCanvas("c3","PosVNeg Histrogram" ,200,10,900,700);
		c3->Divide(2,2);
		c3->SetFillColor(10);
		c3->SetGrid();
		c3->GetFrame()->SetFillColor(10);
		c3->GetFrame()->SetBorderSize(12);
		c3->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c3->cd(1);
			histPosVNegPx->SetXTitle("PosPx");
			histPosVNegPx->SetYTitle("NegPx");
			histPosVNegPx->SetStats(1);
			histPosVNegPx->Draw();
			
			c3->cd(2);
			histPosVNegPy->SetXTitle("PosPy");
			histPosVNegPy->SetYTitle("NegPx");
			histPosVNegPy->SetStats(1);
			histPosVNegPy->Draw();
			
			c3->cd(3);
			histPosVNegPz->SetXTitle("PosPz");
			histPosVNegPz->SetYTitle("NegPx");
			histPosVNegPz->SetStats(1);
			histPosVNegPz->Draw();
			
			c3->cd(4);
			histPosVNegE->SetXTitle("PosE");
			histPosVNegE->SetYTitle("NegPx");
			histPosVNegE->SetStats(1);
			histPosVNegE->Draw();
			
			if(SaveAsRootFile){
				TFile outfile3(rootOutputFileName3, "RECREATE");
				c3->Write(rootOutputFileName3);
				outfile3.Close();
			}
			
	TCanvas *c4 = new TCanvas("c4","Iso and Jet_btag Histrograms" ,200,10,900,700);
		c4->Divide(2,2);
		c4->SetFillColor(10);
		c4->SetGrid();
		c4->GetFrame()->SetFillColor(10);
		c4->GetFrame()->SetBorderSize(12);
		c4->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c4->cd(1);
			histMuonIso->SetXTitle("MuonIso");
			histMuonIso->SetYTitle("Counts");
			histMuonIso->SetStats(1);
			histMuonIso->Draw();
			
			c4->cd(2);
			hist_WMuonVIso->SetXTitle("Muon_Iso");
			hist_WMuonVIso->SetYTitle("NWMuon");
			hist_WMuonVIso->SetStats(1);
			hist_WMuonVIso->Draw();
			
			c4->cd(3);
			histWMuonIso->SetXTitle("WMuonIso");
			histWMuonIso->SetYTitle("Counts");
			histWMuonIso->SetStats(1);
			histWMuonIso->Draw();
			
			c4->cd(4);
			hist_bJetVbtag->SetXTitle("Jet_btag");
			hist_bJetVbtag->SetYTitle("NbJet");
			hist_bJetVbtag->SetStats(1);
			hist_bJetVbtag->Draw();
			
			if(SaveAsRootFile){
				TFile outfile4(rootOutputFileName4, "RECREATE");
				c4->Write(rootOutputFileName4);
				outfile4.Close();
			}
			
	TCanvas *c5 = new TCanvas("c5","MET (And WMuon mass) Histrogram" ,200,10,900,700);
		c5->Divide(1,3);
		c5->SetFillColor(10);
		c5->SetGrid();
		c5->GetFrame()->SetFillColor(10);
		c5->GetFrame()->SetBorderSize(12);
		c5->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c5->cd(1);
			histMET->SetXTitle("MET");
			histMET->SetYTitle("Counts");
			histMET->SetStats(1);
			histMET->Draw();
			
			c5->cd(2);
			histMETxy->SetXTitle("MET_px");
			histMETxy->SetYTitle("MET_py");
			histMETxy->SetStats(1);
			histMETxy->Draw();
			
			c5->cd(3);
			histMET_Phi->SetXTitle("MET_Phi");
			histMET_Phi->SetYTitle("counts");
			histMET_Phi->SetStats(1);
			histMET_Phi->Draw();
			
			// c5_4->cd();
			// histWMuonM->SetXTitle("WMuonMass");
			// histWMuonM->SetYTitle("counts");
			// histWMuonM->SetStats(1);
			// histWMuonM->Draw();
			
			if(SaveAsRootFile){
				TFile outfile5(rootOutputFileName5, "RECREATE");
				c5->Write(rootOutputFileName5);
				outfile5.Close();
			}
			
	TCanvas *c6 = new TCanvas("c6","Transverse Mass Histrogram" ,200,10,900,700);
		c6->Divide(2,2);
		c6->SetFillColor(10);
		c6->SetGrid();
		c6->GetFrame()->SetFillColor(10);
		c6->GetFrame()->SetBorderSize(12);
		c6->Range(0,0,1,1);

		gStyle->SetOptStat(1);
		if(FitInvarientMassHist){
			gStyle->SetOptFit(1);
		}
			
			c6->cd(1);
			histM_WMuon2->SetXTitle("Mt Mass"); //Using WMuon2 here so all defs use Et
			histM_WMuon2->SetYTitle("Counts");
			histM_WMuon2->SetStats(1);
			histM_WMuon2->SetMaximum(mT_MaxCount);
			histM_WMuon2->Draw();
			
			c6->cd(2);
			histM_WMuonPos->SetXTitle("Mt Pos Mass");
			histM_WMuonPos->SetYTitle("Counts");
			histM_WMuonPos->SetStats(1);
			histM_WMuonPos->SetMaximum(mT_MaxCount);
			histM_WMuonPos->Draw();
			
			c6->cd(3);
			histM_WMuonNeg->SetXTitle("Mt Neg Mass");
			histM_WMuonNeg->SetYTitle("Counts");
			histM_WMuonNeg->SetStats(1);
			histM_WMuonNeg->SetMaximum(mT_MaxCount);
			histM_WMuonNeg->Draw();
			
			c6->cd(4);
			histM_WElectron->SetXTitle("Mt e- Mass");
			histM_WElectron->SetYTitle("Counts");
			histM_WElectron->SetStats(1);
			histM_WElectron->SetMaximum(mT_MaxCount);
			histM_WElectron->Draw();
			
			if(SaveAsRootFile){
				TFile outfile6(rootOutputFileName6, "RECREATE");
				c6->Write(rootOutputFileName6);
				outfile6.Close();
			}
			
	TCanvas *c7 = new TCanvas("c7","Transverse Mass Histrogram" ,200,10,900,700);
		c7->Divide(1,3);
		c7->SetFillColor(10);
		c7->SetGrid();
		c7->GetFrame()->SetFillColor(10);
		c7->GetFrame()->SetBorderSize(12);
		c7->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c7->cd(1);
			gPad->SetLogy();
			histM_Top->SetXTitle("Invarient Mass");
			histM_Top->SetYTitle("Counts");
			histM_Top->SetStats(1);
			//histM_Top->SetMaximum(mT_MaxCount);
			histM_Top->Draw();
			
			c7->cd(2);
			gPad->SetLogy();
			histM_TopPos->SetXTitle("Invarient Mass");
			histM_TopPos->SetYTitle("Counts");
			histM_TopPos->SetStats(1);
			//histM_TopPos->SetMaximum(mT_MaxCount);
			histM_TopPos->Draw();
			
			c7->cd(3);
			gPad->SetLogy();
			histM_TopNeg->SetXTitle("Invarient Mass");
			histM_TopNeg->SetYTitle("Counts");
			histM_TopNeg->SetStats(1);
			//histM_TopNeg->SetMaximum(mT_MaxCount);
			histM_TopNeg->Draw();
			
			// c7_4->cd();
			// histM_WElectron->SetXTitle("Mt e- Mass");
			// histM_WElectron->SetYTitle("Counts");
			// histM_WElectron->SetStats(1);
			// histM_WElectron->SetMaximum(mT_MaxCount);
			// histM_WElectron->Draw();
			
			if(SaveAsRootFile){
				TFile outfile7(rootOutputFileName7, "RECREATE");
				c7->Write(rootOutputFileName7);
				outfile7.Close();
			}
		
	TCanvas *c8 = new TCanvas("c8","Invariant Mass of W Histrogram" ,200,10,900,700);
	c8->Divide(1,3);
	c8->SetFillColor(10);
	c8->SetGrid();
	c8->GetFrame()->SetFillColor(10);
	c8->GetFrame()->SetBorderSize(12);
	c8->Range(0,0,1,1);

	gStyle->SetOptStat(1);
		
		c8->cd(1);
		histM_W->SetXTitle("Invarient Mass");
		histM_W->SetYTitle("Counts");
		histM_W->SetStats(1);
		//histM_W->SetMaximum(mT_MaxCount);
		histM_W->Draw();
		
		c8->cd(2);
		histM_WPos->SetXTitle("Invarient Mass");
		histM_WPos->SetYTitle("Counts");
		histM_WPos->SetStats(1);
		//histM_WPos->SetMaximum(mT_MaxCount);
		histM_WPos->Draw();
		
		c8->cd(3);
		histM_WNeg->SetXTitle("Invarient Mass");
		histM_WNeg->SetYTitle("Counts");
		histM_WNeg->SetStats(1);
		//histM_WNeg->SetMaximum(mT_MaxCount);
		histM_WNeg->Draw();
		
		// c8_4->cd();
		// histM_WElectron->SetXTitle("Mt e- Mass");
		// histM_WElectron->SetYTitle("Counts");
		// histM_WElectron->SetStats(1);
		// histM_WElectron->SetMaximum(mT_MaxCount);
		// histM_WElectron->Draw();
		
		if(SaveAsRootFile){
			TFile outfile8(rootOutputFileName8, "RECREATE");
			c8->Write(rootOutputFileName8);
			outfile8.Close();
		}
		

	return 0;
}