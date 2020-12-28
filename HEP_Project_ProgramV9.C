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

gDebug = 1;
using namespace std;
int main(){
	
	/////Controls/////
	
	char rootFilePrefix[50] = "HEPTutorial/files/"; //Prefix to root file, this is the path to it - all the way up until its name.
	char rootFileName[50] = "data.root"; //The root file name itself - this is where the tree is.
	char TreeName1[50] = "events;1";
	
	char rootOutputFileName1[50] = "InvarientMassHist.root"; //Name of an output root file 
	char rootOutputFileName2[50] = "PosVNegPrime.root"; //Name of an output root file
	char rootOutputFileName3[50] = "PosVNeg.root"; //Name of an output root file
	char rootOutputFileName4[50] = "MuonIso.root"; //Name of an output root file
	char rootOutputFileName5[50] = "MET_WMuonMass.root";
	char rootOutputFileName6[50] = "TransverseMassHist.root";
	char rootOutputFileName7[50] = "InvarientMassHist_TOP.root";
	
	double MuonIsoBound_ForZ = 1; //Upper bound to muon iso when looking for ZMuons 
	double MuonIsoBound_ForW = 1.5; //Upper bound to muon iso when looking for WMuons 
	double ElectronIsoBound_ForW = 1000; //Upper bound to electron iso when looking for WElectrons 
	
	double Jet_btag_Bound = 2; //Lower bound to Jet_btag
	
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
	
	TH1F *histM_Muon = new TH1F("histM_Muon", "Invarient Mass hist: Muon", xNChannels, 0, 150);
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
	
	TH1F *histMuonIso = new TH1F("histMuonIso", "MuonIso hist", 60, 0, 20);
	TH1F *histMuonPt = new TH1F("histMuonPt", "MuonPt hist", xNChannels, xMin, xMax);
	TH1F *histWMuonIso = new TH1F("histMuonIso", "MuonIso hist", 60, 0, 20);
	TH1F *histWElectronIso = new TH1F("histElectronIso", "ElectronIso hist", 60, 0, 20);
	 
	TH2F *histMETxy = new TH2F("histMETxy", "METxy hist", xNChannels, xMin, xMax, yNChannels, yMin, yMax);
	TH1F *histMET = new TH1F("histMET", "MET hist", xNChannels, xMin, xMax);
	TH1F *histMET_Phi = new TH1F("histMET_Phi", "MET_Phi hist", xNChannels, 0, 360);
	
	TH1F *histM_TopPos = new TH1F("histM_TopPos", "Invarient Mass hist: TopPos", xNChannels, 0, 200);
	TH1F *histM_TopNeg = new TH1F("histM_TopNeg", "Invarient Mass hist: TopNeg", xNChannels, 0, 200);
	TH1F *histM_Top = new TH1F("histM_Top", "Invarient Mass hist: Top", xNChannels, 0, 200);
	
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
	// t1->SetBranchAddress("triggerIsoMu24",&triggerIsoMu24);
	// t1->SetBranchAddress("EventWeight",&EventWeight);
	
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
		
		TLorentzVector bJet_prime_Pos;
		TLorentzVector bJet_prime_Neg;
		TLorentzVector TopPos;
		TLorentzVector TopNeg;
		
		TLorentzVector MuonTemp;
	
	int NEntries1 = t1->GetEntries();
	cout<<"NEntries 1: "<<NEntries1<<endl;
	
	//NEntries1 = NEntries1/44; //Cut for testing
	
	for(int i = 0; i < NEntries1; i++){
		cout<<"*************"<<endl;
		cout<<"i: "<<i<<endl;
		t1->GetEntry(i);
		pos_Muon_prime = new TLorentzVector;
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
		
		bool WPosHappened = false;
		bool WNegHappened = false;
		
		//For the Z//
		if(NMuon == 2){
			cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			for(int j = 0; j < NMuon; j++){
				//TLorentzVector MuonTemp
				MuonTemp.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
				if((Muon_Iso[j] < MuonIsoBound_ForZ) && (MuonTemp.P() > 0)){
					cout<<"Z| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					if(Muon_Charge[j] == 1){
						pos_Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
						// histPosPxPrime->Fill(Muon_Px[j]);
						// histPosPyPrime->Fill(Muon_Py[j]);
						// histPosPzPrime->Fill(Muon_Pz[j]);
						// histPosEPrime->Fill(Muon_E[j]);
						double Muon_Iso_pos = Muon_Iso[j];
						double Muon_Pt_pos = pos_Muon_prime.Pt();
					}
					if(Muon_Charge[j] == -1){
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
			if(((Muon_Charge[0] == 1 && Muon_Charge[1] == -1) || (Muon_Charge[0] == -1 && Muon_Charge[1] == 1)) && ((Muon_Iso[0] < MuonIsoBound_ForZ) && (Muon_Iso[1] < MuonIsoBound_ForZ)) && (MuonTemp.P() > 0)){
				tot_Muon_prime = pos_Muon_prime + neg_Muon_prime;
				histM_Muon->Fill(tot_Muon_prime.M());
				
				histMuonIso->Fill(Muon_Iso_pos);  ///neg_Muon_prime.Pt());
				
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
			}
		}
		
		//For the WMuon//
		if(NMuon == 1){
			cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			for(int j = 0; j < NMuon; j++){
				if(Muon_Iso[j] < MuonIsoBound_ForW){
					histWMuonIso->Fill(Muon_Iso[j]);
					cout<<"W| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
					WMuon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
						
					if(Muon_Charge[j] > 0){
						cout<<"W+| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
						cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
						WMuon_prime_Pos.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
						WPosHappened = true;
						
					}
					if(Muon_Charge[j] < 0){
						cout<<"W-| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
						cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
						WMuon_prime_Neg.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
						WNegHappened = true;
						
					}
				}
			}
			double MET_prime = sqrt((MET_px*MET_px) + (MET_py*MET_py));
			if(MET_px != 0){
				double MET_prime_Phi_Reference = abs(atan(MET_py/MET_px));
			}
			double MET_prime_Phi;
			if(MET_px > 0 && MET_py > 0){
				MET_prime_Phi = MET_prime_Phi_Reference;
			}
			if(MET_px > 0 && MET_py < 0){
				MET_prime_Phi = Pi + MET_prime_Phi_Reference;
			}
			if(MET_px < 0 && MET_py > 0){
				MET_prime_Phi = Pi - MET_prime_Phi_Reference;
			}
			if(MET_px < 0 && MET_py < 0){
				MET_prime_Phi = 2*Pi - MET_prime_Phi_Reference;
			}
			if(MET_px == 0 && MET_py > 0){
				MET_prime_Phi = Pi/2.;
			}
			if(MET_px == 0 && MET_py < 0){
				MET_prime_Phi = (3./2.)*Pi;
			}
			if(MET_px > 0 && MET_py == 0){
				MET_prime_Phi = 0;
			}
			if(MET_px < 0 && MET_py == 0){
				MET_prime_Phi = Pi;
			}
			cout<<MET_prime_Phi<<endl;
			histMETxy->Fill(MET_px, MET_py);
			histMET->Fill(MET_prime);
			histMET_Phi->Fill(MET_prime_Phi * (1./Deg2Rad));
			
			double WMuon_prime_Phi = WMuon_prime.Phi();
			double WMuon_prime_Pt = WMuon_prime.Pt();
			double WMuon_prime_Et = WMuon_prime.Et();
			
			double WMuon_prime_Phi_Pos = WMuon_prime_Pos.Phi();
			double WMuon_prime_Pt_Pos = WMuon_prime_Pos.Pt();
			double WMuon_prime_Et_Pos = WMuon_prime_Pos.Et();
			
			double WMuon_prime_Phi_Neg = WMuon_prime_Neg.Phi();
			double WMuon_prime_Pt_Neg = WMuon_prime_Neg.Pt();
			double WMuon_prime_Et_Neg = WMuon_prime_Neg.Et();
			
			double mT_mine = sqrt(2.*MET_prime*WMuon_prime_Pt*(1-cos(MET_prime_Phi - WMuon_prime_Phi)));
			double mT_mine2 = sqrt(2.*MET_prime*WMuon_prime_Et*(1-cos(MET_prime_Phi - WMuon_prime_Phi)));
			double mT_mine3 = sqrt(2.*(MET_prime*WMuon_prime_Et-(MET_prime*WMuon_prime_Pt*cos(MET_prime_Phi - WMuon_prime_Phi))));
			
			double mT_mine_Pos = sqrt(2.*MET_prime*WMuon_prime_Et_Pos*(1-cos(MET_prime_Phi - WMuon_prime_Phi_Pos)));
			double mT_mine_Neg = sqrt(2.*MET_prime*WMuon_prime_Et_Neg*(1-cos(MET_prime_Phi - WMuon_prime_Phi_Neg)));
			
			histM_WMuon->Fill(mT_mine);
			histM_WMuon2->Fill(mT_mine2);
			histM_WMuon3->Fill(mT_mine3);
			
			histM_WMuonPos->Fill(mT_mine_Pos);
			histM_WMuonNeg->Fill(mT_mine_Neg);
			
		}
		//For the WElectron//
		if(NElectron == 1){
			cout<<"NElectron: "<<NElectron<<endl;
			histNElectron->Fill(NElectron);
			for(int j = 0; j < NElectron; j++){
				if(Electron_Iso[j] < ElectronIsoBound_ForW){
					histWElectronIso->Fill(Electron_Iso[j]);
					cout<<"W| i: "<<i<<" j: "<<j<<" Electron_Charge: "<<Electron_Charge[j]<<endl;
					WElectron_prime.SetPxPyPzE(Electron_Px[j], Electron_Py[j], Electron_Pz[j], Electron_E[j]);
						
					if(Electron_Charge[j] > 0){
						cout<<"W+| i: "<<i<<" j: "<<j<<" Electron_Charge: "<<Electron_Charge[j]<<endl;
						WElectron_prime_Pos.SetPxPyPzE(Electron_Px[j], Electron_Py[j], Electron_Pz[j], Electron_E[j]);
						
					}
					if(Electron_Charge[j] < 0){
						cout<<"W-| i: "<<i<<" j: "<<j<<" Electron_Charge: "<<Electron_Charge[j]<<endl;
						WElectron_prime_Neg.SetPxPyPzE(Electron_Px[j], Electron_Py[j], Electron_Pz[j], Electron_E[j]);
						
					}
				}
			}
			double MET_prime = sqrt((MET_px*MET_px) + (MET_py*MET_py));
			if(MET_px != 0){
				double MET_prime_Phi_Reference = abs(atan(MET_py/MET_px));
			}
			double MET_prime_Phi;
			if(MET_px > 0 && MET_py > 0){
				MET_prime_Phi = MET_prime_Phi_Reference;
			}
			if(MET_px > 0 && MET_py < 0){
				MET_prime_Phi = Pi + MET_prime_Phi_Reference;
			}
			if(MET_px < 0 && MET_py > 0){
				MET_prime_Phi = Pi - MET_prime_Phi_Reference;
			}
			if(MET_px < 0 && MET_py < 0){
				MET_prime_Phi = 2*Pi - MET_prime_Phi_Reference;
			}
			if(MET_px == 0 && MET_py > 0){
				MET_prime_Phi = Pi/2.;
			}
			if(MET_px == 0 && MET_py < 0){
				MET_prime_Phi = (3./2.)*Pi;
			}
			if(MET_px > 0 && MET_py == 0){
				MET_prime_Phi = 0;
			}
			if(MET_px < 0 && MET_py == 0){
				MET_prime_Phi = Pi;
			}
			cout<<MET_prime_Phi<<endl;
			histMETxy->Fill(MET_px, MET_py);
			histMET->Fill(MET_prime);
			histMET_Phi->Fill(MET_prime_Phi * (1./Deg2Rad));
			
			double WElectron_prime_Phi = WElectron_prime.Phi();
			double WElectron_prime_Pt = WElectron_prime.Pt();
			double WElectron_prime_Et = WElectron_prime.Et();
			
			double WElectron_prime_Phi_Pos = WElectron_prime_Pos.Phi();
			double WElectron_prime_Pt_Pos = WElectron_prime_Pos.Pt();
			double WElectron_prime_Et_Pos = WElectron_prime_Pos.Et();
			
			double WElectron_prime_Phi_Neg = WElectron_prime_Neg.Phi();
			double WElectron_prime_Pt_Neg = WElectron_prime_Neg.Pt();
			double WElectron_prime_Et_Neg = WElectron_prime_Neg.Et();
			
			double mT_mine = sqrt(2.*MET_prime*WElectron_prime_Et*(1-cos(MET_prime_Phi - WElectron_prime_Phi)));
			double mT_mine_Pos = sqrt(2.*MET_prime*WElectron_prime_Et_Pos*(1-cos(MET_prime_Phi - WElectron_prime_Phi_Pos)));
			double mT_mine_Neg = sqrt(2.*MET_prime*WElectron_prime_Et_Neg*(1-cos(MET_prime_Phi - WElectron_prime_Phi_Neg)));
			
			histM_WElectron->Fill(mT_mine);
			histM_WElectronPos->Fill(mT_mine_Pos);
			histM_WElectronNeg->Fill(mT_mine_Neg);
		}
		if(NJet == 1){
			cout<<"NJet: "<<NJet<<endl;
			for(int j = 0; j < NJet; j++){
				if(Jet_ID[j] && (Jet_btag[j] > Jet_btag_Bound)){
					cout<<"Found good bJet"<<endl;
					if(WPosHappened){
						cout<<"b| i: "<<i<<" j: "<<j<<endl;
						bJet_prime_Pos.SetPxPyPzE(Jet_Px[j], Jet_Py[j], Jet_Pz[j], Jet_E[j]);
						cout<<"Made 4 Vect"<<endl;
						TLorentzVector Top_Pos = bJet_prime_Pos + WMuon_prime_Pos;
						cout<<"Added Vects"<<endl;
						histM_TopPos->Fill(Top_Pos.M());
						histM_Top->Fill(Top_Pos.M());
						cout<<"Filled hists"<<endl;
					}
					if(WNegHappened){
						cout<<"b-| i: "<<i<<" j: "<<j<<endl;
						bJet_prime_Neg.SetPxPyPzE(Jet_Px[j], Jet_Py[j], Jet_Pz[j], Jet_E[j]);
						cout<<"Made 4 Vect"<<endl;
						TLorentzVector Top_Neg = bJet_prime_Neg + WMuon_prime_Neg;
						cout<<"Added Vects"<<endl;
						histM_TopNeg->Fill(Top_Neg.M());
						histM_Top->Fill(Top_Neg.M());
						cout<<"Filled hists"<<endl;
					}
				}
				// if(WPosHappened){
					// cout<<"t+| i: "<<i<<endl;
					// Top_Pos = bJet_prime_Pos + WMuon_prime_Pos;
					// cout<<"Added Vects"<<endl;
					// histM_TopPos->Fill(Top_Pos.M());
					// histM_Top->Fill(Top_Pos.M());
					// cout<<"Filled hists"<<endl;
				// }
				// if(WNegHappened){
					// cout<<"t-| i: "<<i<<endl;
					// Top_Neg = bJet_prime_Neg + WMuon_prime_Neg;
					// cout<<"Added Vects"<<endl;
					// histM_TopNeg->Fill(Top_Neg.M());
					// histM_Top->Fill(Top_Neg.M());
					// cout<<"Filled hists"<<endl;
				// }
			}	
		}
		// delete pos_Muon_prime;
		// delete neg_Muon_prime;
		// delete tot_Muon_prime;
		// delete pos_Muon;
		// delete neg_Muon;
		// delete tot_Muon;
		
		// delete WMuon_prime;
		// delete WMuon_prime_Pos;
		// delete WMuon_prime_Neg;

		// delete WElectron_prime;
		// delete WElectron_prime_Pos;
		// delete WElectron_prime_Neg;
		
		// delete bJet_prime_Pos;
		// delete bJet_prime_Neg;
		// delete TopPos;
		// delete TopNeg;
		
		// delete MuonTemp;
		
		// delete WPosHappened = false;
		// delete WNegHappened = false;
	}
	
	if(FitInvarientMassHist){
		TF1 *f1 = new TF1("f1","gaus",80,100);
		TF1 *f2 = new TF1("f2","([0]*exp(-((x-[1])^2)/(2*[2]^2)))",60,85);
		
		f2->SetParameters(0, 200); //Normal
		f2->SetParameters(1, 70); //Centroid
		f2->SetParameters(2, 25); //Width
	}
	
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
			
			c1_1->cd();
			histM_Muon->SetXTitle("Invarient Mass");
			histM_Muon->SetYTitle("Counts");
			histM_Muon->SetStats(1);
			if(FitInvarientMassHist){
				histM_Muon->Fit(f1, "M");
			}
			histM_Muon->Draw();
			
			c1_2->cd();
			histM_WMuon->SetXTitle("Mt Mass");
			histM_WMuon->SetYTitle("Counts");
			histM_WMuon->SetStats(1);
			if(FitInvarientMassHist){
				histM_WMuon->Fit(f2, "M");
			}
			histM_WMuon->Draw();
			
			c1_3->cd();
			histM_WMuon2->SetXTitle("Mt 2 Mass");
			histM_WMuon2->SetYTitle("Counts");
			histM_WMuon2->SetStats(1);
			histM_WMuon2->Draw();
			
			c1_4->cd();
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
			
			c2_1->cd();
			histPosVNegPxPrime->SetXTitle("PosPxPrime");
			histPosVNegPxPrime->SetYTitle("NegPxPrime");
			histPosVNegPxPrime->SetStats(1);
			histPosVNegPxPrime->Draw();
			
			c2_2->cd();
			histPosVNegPyPrime->SetXTitle("PosPyPrime");
			histPosVNegPyPrime->SetYTitle("NegPxPrime");
			histPosVNegPyPrime->SetStats(1);
			histPosVNegPyPrime->Draw();
			
			c2_3->cd();
			histPosVNegPzPrime->SetXTitle("PosPzPrime");
			histPosVNegPzPrime->SetYTitle("NegPxPrime");
			histPosVNegPzPrime->SetStats(1);
			histPosVNegPzPrime->Draw();
			
			c2_4->cd();
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
			
			c3_1->cd();
			histPosVNegPx->SetXTitle("PosPx");
			histPosVNegPx->SetYTitle("NegPx");
			histPosVNegPx->SetStats(1);
			histPosVNegPx->Draw();
			
			c3_2->cd();
			histPosVNegPy->SetXTitle("PosPy");
			histPosVNegPy->SetYTitle("NegPx");
			histPosVNegPy->SetStats(1);
			histPosVNegPy->Draw();
			
			c3_3->cd();
			histPosVNegPz->SetXTitle("PosPz");
			histPosVNegPz->SetYTitle("NegPx");
			histPosVNegPz->SetStats(1);
			histPosVNegPz->Draw();
			
			c3_4->cd();
			histPosVNegE->SetXTitle("PosE");
			histPosVNegE->SetYTitle("NegPx");
			histPosVNegE->SetStats(1);
			histPosVNegE->Draw();
			
			if(SaveAsRootFile){
				TFile outfile3(rootOutputFileName3, "RECREATE");
				c3->Write(rootOutputFileName3);
				outfile3.Close();
			}
			
	TCanvas *c4 = new TCanvas("c4","MuonIso and Pt Histrogram" ,200,10,900,700);
		c4->Divide(2,2);
		c4->SetFillColor(10);
		c4->SetGrid();
		c4->GetFrame()->SetFillColor(10);
		c4->GetFrame()->SetBorderSize(12);
		c4->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c4_1->cd();
			histMuonIso->SetXTitle("MuonIso");
			histMuonIso->SetYTitle("Counts");
			histMuonIso->SetStats(1);
			histMuonIso->Draw();
			
			c4_2->cd();
			histMuonPt->SetXTitle("Pt");
			histMuonPt->SetYTitle("Counts");
			histMuonPt->SetStats(1);
			histMuonPt->Draw();
			
			c4_3->cd();
			histWMuonIso->SetXTitle("WMuonIso");
			histWMuonIso->SetYTitle("Counts");
			histWMuonIso->SetStats(1);
			histWMuonIso->Draw();
			
			c4_4->cd();
			histElectronIso->SetXTitle("WElectronIso");
			histElectronIso->SetYTitle("Counts");
			histElectronIso->SetStats(1);
			histElectronIso->Draw();
			
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
			
			c5_1->cd();
			histMET->SetXTitle("MET");
			histMET->SetYTitle("Counts");
			histMET->SetStats(1);
			histMET->Draw();
			
			c5_2->cd();
			histMETxy->SetXTitle("MET_px");
			histMETxy->SetYTitle("MET_py");
			histMETxy->SetStats(1);
			histMETxy->Draw();
			
			c5_3->cd();
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
			
			c6_1->cd();
			histM_WMuon2->SetXTitle("Mt Mass"); //Using WMuon2 here so all defs use Et
			histM_WMuon2->SetYTitle("Counts");
			histM_WMuon2->SetStats(1);
			histM_WMuon2->SetMaximum(mT_MaxCount);
			histM_WMuon2->Draw();
			
			c6_2->cd();
			histM_WMuonPos->SetXTitle("Mt Pos Mass");
			histM_WMuonPos->SetYTitle("Counts");
			histM_WMuonPos->SetStats(1);
			histM_WMuonPos->SetMaximum(mT_MaxCount);
			histM_WMuonPos->Draw();
			
			c6_3->cd();
			histM_WMuonNeg->SetXTitle("Mt Neg Mass");
			histM_WMuonNeg->SetYTitle("Counts");
			histM_WMuonNeg->SetStats(1);
			histM_WMuonNeg->SetMaximum(mT_MaxCount);
			histM_WMuonNeg->Draw();
			
			c6_4->cd();
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
			
			c7_1->cd();
			histM_Top->SetXTitle("Invarient Mass");
			histM_Top->SetYTitle("Counts");
			histM_Top->SetStats(1);
			//histM_Top->SetMaximum(mT_MaxCount);
			histM_Top->Draw();
			
			c7_2->cd();
			histM_TopPos->SetXTitle("Invarient Mass");
			histM_TopPos->SetYTitle("Counts");
			histM_TopPos->SetStats(1);
			//histM_TopPos->SetMaximum(mT_MaxCount);
			histM_TopPos->Draw();
			
			c7_3->cd();
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
		
		

	return 0;
}