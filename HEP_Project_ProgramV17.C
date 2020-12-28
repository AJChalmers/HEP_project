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

//V14: Changes made at the end of 11/12 class for data overlay. Now adding more to our data.

//V15: Changes made with Nick 11/15: Fixed NIsoMuon counter, needed to be at top of loop. Also fixed bug with "N_IsoMuon == 1;" not being "NIsoMuon == 1;"
//Notes from 11/15: Used NMuonIso instead of N_MuonIso and moved counter

//V16_TEMP: Changes being made to create Zs as done in the solution code. Attempt at seperating top - antitop onto seperate distributions. TEMP: Changes made during 11/17 class not present here (Should only have been some 2D IsoVPt hists)
//V16: Temp but with changes made during class 11/17.

//V17: Added 2D histogram of invariant mass of W vs invariant mass of top
///IT WORKS!!!!!

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
	char rootFileName[50] = "ttbar.root"; //The root file name itself - this is where the tree is.
	char TreeName1[50] = "events;1";
	
	char rootOutputFileName1[50] = "ttbar_InvarientMassHist.root"; //Name of an output root file 
	char rootOutputFileName2[50] = "ttbar_PosVNegPrime.root"; //Name of an output root file
	char rootOutputFileName3[50] = "ttbar_PosVNeg.root"; //Name of an output root file
	char rootOutputFileName4[50] = "ttbar_MuonIso_btag.root"; //Name of an output root file
	char rootOutputFileName5[50] = "ttbar_MET_WMuonMass.root";
	char rootOutputFileName6[50] = "ttbar_TransverseMassHist.root";
	char rootOutputFileName7[50] = "ttbar_InvarientMassHist_TOP.root";
	char rootOutputFileName8[50] = "ttbar_InvarientMassHist_W.root";
	char rootOutputFileName9[50] = "ttbar_MuonIsoVsPt.root";
	char rootOutputFileName10[50] = "ttbar_InvarientMassHist_Z.root";
	char rootOutputFileName11[50] = "ttbar_InvarientMassHist_TOP_Comparison.root";
	char rootOutputFileName12[50] = "ttbar_InvarientMassHist_WComparison.root";
	char rootOutputFileName13[50] = "ttbar_InvarientMassHist_qJet.root";
	char rootOutputFileName14[50] = "ttbar_histM_WHadVsM_TopHad.root";
	
	double MuonIsoBound_ForZ = 1000000000; //Upper bound to muon iso when looking for ZMuons 
	double MuonIsoBound_ForW = 1.5; //Upper bound to muon iso when looking for WMuons 
	double ElectronIsoBound_ForW = 1000000000; //Upper bound to electron iso when looking for WElectrons 
	
	double RelIsoBound = 0.1; //Iso/Pt bound aka MuonRelIsoCut in MyAnalysis.C
	double Jet_btag_Bound = 1.74; //Lower bound to Jet_btag
	double MuonPtCut = 25.;
	
	double mT_MaxCount = 1000; //Max counts for mT hists. About 1,000 for full data set and 25 for 1/44
	
	double Z_Mass_LowerBound = 90;  //Z_Mass bounds used for cut on muons used for top
	double Z_Mass_UpperBound = 90;
	
	double W_Mass_LowerBound = 65; //W_Mass bounds used for Hadronic branch 
	double W_Mass_UpperBound = 95;
	
	double xMin = -100;
	double xMax = 100;
	double yMin = -300;
	double yMax = 300;
	
	const int xNChannels = 1000; //Number of channels for root hist on x axis
	const int yNChannels = 500; //Number of channels for root hist on y axis 
	
	bool SaveAsRootFile = true; //If true, will save canvases to root file instead of printing to screen.
	bool FitInvarientMassHist = true; //If true, will fit Invarient mass hist for Z from muons with gaus, and W from muons with gaus
	bool DoCMRotation = false; //If true, will do rotation to CM for muons. Will make code take longer to run. Note: When false some hists will be empty!
	
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
	
	TH1F *histM_Muon = new TH1F("histM_Muon", "Invarient Mass hist: Z", 60 , 60, 120);
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
	TH1F *histM_Top_ZCut = new TH1F("histM_Top_ZCut", "Invarient Mass hist: Top (ZCut)", 50, 0, 500);
	
	
	TH1F *histM_WPos = new TH1F("histM_WPos", "Invarient Mass hist: WPos", xNChannels, 0, 200);
	TH1F *histM_WNeg = new TH1F("histM_WNeg", "Invarient Mass hist: WNeg", xNChannels, 0, 200);
	TH1F *histM_W = new TH1F("histM_W", "Invarient Mass hist: W", xNChannels, 0, 200);
	
	TH1F *histDataOverlay = new TH1F("histDataOverlay", "Data Overlay", 50, 0, 500);
	
	TH1F *histWMuonPt = new TH1F("histWMuonPt", "WMuonPt hist", xNChannels, 10, 100);
	
	TH2F *histMuonIsoVsPt = new TH2F("histMuonIsoVsPt", "MuonIsoVsPt hist", xNChannels, 0, 300, yNChannels, 0, 15);
	TH2F *histMuonIsoVsPt_AfterCut = new TH2F("histMuonIsoVsPt_AfterCut", "MuonIsoVsPt_AfterCut hist", xNChannels, 0, 300, yNChannels, 0, 15);
	
	TH1F *histM_Z = new TH1F("histM_Z", "Invarient Mass hist: Z (Solution)", 50 , 0, 500);
	
	TH1F *histM_WHad = new TH1F("histM_WHad", "Invarient Mass hist: W Hadronic", 50, 0, 500);
	TH1F *histM_TopHad = new TH1F("histM_TopHad", "Invarient Mass hist: Top Hadronic", 50, 0, 500);
	TH1F *histM_Wt = new TH1F("histM_Wt", "Transverse Mass hist: W", 50, 0, 500);
	
	TH1F *histM_qJet = new TH1F("histM_qJet", "qJet hist", 40, 0, 20);
	
	TH2F *histM_WHadVsM_TopHad = new TH2F("histM_WHadVsM_TopHad", "M_WHadVsM_TopHad hist", 100, 100, 300, 100, 0, 200);
	TH2F *histM_WHadVsM_TopHad_NoCut = new TH2F("histM_WHadVsM_TopHad_NoCut", "M_WHadVsM_TopHad_NoCut hist", 100, 100, 300, 100, 0, 200);

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
		
		TLorentzVector bJet;
		TLorentzVector bJet_Pos;
		TLorentzVector bJet_Neg;
		TLorentzVector bJet_Had;
		
		TLorentzVector Top;
		TLorentzVector TopPos;
		TLorentzVector TopNeg;
		TLorentzVector Top_ZCut;
		TLorentzVector Top_Had;
		
		TLorentzVector MuonTemp;
		
		TLorentzVector MET;
		
		TLorentzVector Neutrino1;
		TLorentzVector Neutrino2;
		TLorentzVector Neutrino1_Pos;
		TLorentzVector Neutrino1_Neg;
		TLorentzVector Neutrino2_Pos;
		TLorentzVector Neutrino2_Neg;
		
		TLorentzVector W;
		TLorentzVector WPos;
		TLorentzVector WNeg;
		TLorentzVector W_Had;
		
		TLorentzVector Muon_prime;
		TLorentzVector Muon_prime_Iso1;
		TLorentzVector Muon_prime_Iso2;
		
		TLorentzVector Z_prime;
		
		TLorentzVector qJet1;
		TLorentzVector qJet2;
		
		TLorentzVector qJet1_temp;
		TLorentzVector qJet2_temp;
		
		double Muon_Pt_pos;
		double Muon_Pt_neg;
		
		double MET_pz;
		double MET_E;
		
		double RelIso;
		
		double Z_Mass;
		double MET_prime_Phi_Reference;
		
		int N1_counter = 0;
		int jet_counter = 0;
		int NIsoMuon_counter = 0;
		int N_selectedEvents = 0;
		
		int TotalEvents = 0;
		
		int Muon1_Charge;
		int Muon2_Charge;
		
		int TerminalCounter = 0;
		int  FirstCounter = 0;
        int  SecondCounter = 0;
        int  ThirdCounter = 0;
        int  FourthCounter = 0;
	
	int NEntries1 = t1->GetEntries();
	cout<<"NEntries 1: "<<NEntries1<<endl;
	
	for(int i = 0; i < NEntries1; i++){
		
		TotalEvents++;
		if (TotalEvents % 10000 == 0){
			cout << "Next event -----> " << TotalEvents << endl;
		}
		
		//cout<<"*************"<<endl;
		//cout<<"i: "<<i<<endl;
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
		
		
		//For the Z OLD//
		int N_IsoMuon = 0;
		if(NMuon >= 2){
			//cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			for(int j = 0; j < NMuon; j++){
				histMuonIso->Fill(Muon_Iso[j]); 
				MuonTemp.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
				RelIso = Muon_Iso[j]/MuonTemp.Pt();
				if((RelIso < RelIsoBound) && (MuonTemp.P() > 0)){
					++N_IsoMuon;
					//cout<<"Z| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
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
				
				Z_Mass = tot_Muon_prime.M();
				
				//histMuonIso->Fill(Muon_Iso_pos);  ///neg_Muon_prime.Pt());
				
				histMuonPt->Fill(Muon_Pt_pos);
				histMuonPt->Fill(Muon_Pt_neg);
				
				//Uncomment for CM rotation//
				
				// TVector3 tot_Muon_boost;
				// tot_Muon_boost = tot_Muon_prime.BoostVector();
				// TLorentzRotation l;
				// TLorentzRotation l2;
				// l.Boost(0,0,0);
				// l.Boost(tot_Muon_boost);
				// l2 = l.Inverse();
				// neg_Muon = l2*neg_Muon_prime;
				// pos_Muon = l2*pos_Muon_prime;
				// tot_Muon = pos_Muon + neg_Muon;
				
				// histPosVNegPx->Fill(pos_Muon.Px(), neg_Muon.Px());
				// histPosVNegPy->Fill(pos_Muon.Px(), neg_Muon.Px());
				// histPosVNegPz->Fill(pos_Muon.Px(), neg_Muon.Px());
				// histPosVNegE->Fill(pos_Muon.Px(), neg_Muon.Px());
				
			//}
		}
		
	// W Transverse //
	//For the WMuon//
	if(NMuon == 1){
		//cout<<"NMuon: "<<NMuon<<endl;
		//histNMuon->Fill(NMuon);
		for(int j = 0; j < NMuon; j++){
			if(Muon_Iso[j] < MuonIsoBound_ForW){
				histWMuonIso->Fill(Muon_Iso[j]);
				//cout<<"W| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
				WMuon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					
				if(Muon_Charge[j] > 0){
					//cout<<"W+| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					WMuon_prime_Pos.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					
				}
				if(Muon_Charge[j] < 0){
					//cout<<"W-| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					WMuon_prime_Neg.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					
				}
			}
		}
		double MET_prime = sqrt((MET_px*MET_px) + (MET_py*MET_py));
		if(MET_px != 0){
			MET_prime_Phi_Reference = abs(atan(MET_py/MET_px));
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
		//cout<<MET_prime_Phi<<endl;
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
		
		histM_Wt->Fill(mT_mine);
		
		// histM_WMuon2->Fill(mT_mine2);
		// histM_WMuon3->Fill(mT_mine3);
		
		// histM_WMuonPos->Fill(mT_mine_Pos);
		// histM_WMuonNeg->Fill(mT_mine_Neg);
		
	}
		
		// Isolation selecton for Muons and Z construction alike solution //
		int NIsoMuon = 0;
		for(int j = 0; j < NMuon; j++){ // Muon loop: Looks at every loop
			//cout<<"NMuon: "<<NMuon<<endl;
			histNMuon->Fill(NMuon);
			
			//cout<<"W| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
			//cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
			Muon_prime.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
			RelIso = Muon_Iso[j]/Muon_prime.Pt();
			histWMuonPt->Fill(Muon_prime.Pt());
			
			if(RelIso < RelIsoBound){
				NIsoMuon++;
				histWMuonIso->Fill(Muon_Iso[j]);
				
				if(NIsoMuon == 1){
					//cout<<"Muon| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					//cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
					Muon_prime_Iso1.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					
					Muon1_Charge = Muon_Charge[j];
					
					//++N1_counter;
				}
				if(NIsoMuon == 2){
					//cout<<"Muon| i: "<<i<<" j: "<<j<<" Muon_Charge: "<<Muon_Charge[j]<<endl;
					//cout<<"Muon_Px: "<<Muon_Px[j]<<" Muon_Py: "<<Muon_Py[j]<<" Muon_Pz: "<<Muon_Pz[j]<<endl;
					Muon_prime_Iso2.SetPxPyPzE(Muon_Px[j], Muon_Py[j], Muon_Pz[j], Muon_E[j]);
					
					Muon2_Charge = Muon_Charge[j];
				}
				//cout << NIsoMuon << endl;
			}
		}
		// if((NIsoMuon > 1) && triggerIsoMu24 && (Muon_prime_Iso1.Pt() > MuonPtCut) && (NBJet < 2)){
			// Z_prime = Muon_prime_Iso1 + Muon_prime_Iso2;
			// histM_Z->Fill(Z_prime.M());
		// }	
		//cout << "number of N1s: " << N1_counter << endl;
		
		
		// Accept efficiency (jets) //
		int NBJet;
		bool IsSelected = false;
		if(NIsoMuon == 1) {
			++N1_counter;
			if(Muon_prime_Iso1.Pt() > MuonPtCut) {
				NBJet = 0;
				for(int jet = 0; jet < NJet; jet++) {
					if(Jet_ID[jet] && (Jet_btag[jet] > Jet_btag_Bound)){
						++NBJet;
						++jet_counter;
						//cout << "NBJet: " << NBJet << endl;
					}
				}
				if(NBJet > 1){ // > 1
					if(triggerIsoMu24){
						IsSelected = true;
						++N_selectedEvents; 
							
						//cout << "\n\n\nIS SELECTED SET TO TRUE!!!!!\n\n\n";
					}
				}
			}
		}
		//cout << "number of N1s: " << N1_counter << endl;
		//cout << "Number of Bjets: " << jet_counter << endl;
		//cout << "Number of selected events: " << N_selectedEvents << endl;
		
		// Z Creation (Soltution)//
		if((NIsoMuon > 1) && triggerIsoMu24 && (Muon_prime_Iso1.Pt() > MuonPtCut) ){ //&& (NBJet > 0)){
			Z_prime = Muon_prime_Iso1 + Muon_prime_Iso2;
			histM_Z->Fill(Z_prime.M());
		}	
		
		// Exercise 4
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
					  Neutrino1.SetPxPyPzE(MET_px, MET_py, MET_pz, MET_E);
					  W = Neutrino1 + Muon_prime_Iso1;
					  histM_W->Fill(W.M());
					  
					  Top = W + bJet;
					  histM_Top->Fill(Top.M());
					  
					  if((Z_Mass < Z_Mass_LowerBound) || (Z_Mass > Z_Mass_UpperBound)){
						  Top_ZCut = Top;
						  histM_Top_ZCut->Fill(Top_ZCut.M());
					  
						  if(Muon1_Charge > 0){
							  bJet_Pos = bJet;
							  WPos = W;
							  TopPos = Top;
							  histM_TopPos->Fill(TopPos.M());
							  histM_WPos->Fill(WPos.M());
						  }
						  if(Muon1_Charge < 0){
							  bJet_Neg = bJet;
							  WNeg = W;
							  TopNeg = Top;
							  histM_TopNeg->Fill(TopNeg.M());
							  histM_WNeg->Fill(WNeg.M());
						  }
						}
					}
					// second solution
					if(D > 0){
					  double MET_pz =  Muon_prime_Iso1.Pz() * B / A - sqrt(D) / A;
					  double MET_E = sqrt(MET_px * MET_px + MET_py * MET_py + MET_pz * MET_pz);
					  Neutrino2.SetPxPyPzE(MET_px, MET_py, MET_pz, MET_E);
					  W = Neutrino2 + Muon_prime_Iso1;
					  histM_W->Fill(W.M());
					  
					  Top = W + bJet;
					  histM_Top->Fill(Top.M());
					  
					   if((Z_Mass < Z_Mass_LowerBound) || (Z_Mass > Z_Mass_UpperBound)){
						  Top_ZCut = Top;
						  histM_Top_ZCut->Fill(Top_ZCut.M());
						  
						  if(Muon1_Charge > 0){
							  bJet_Pos = bJet;
							  WPos = W;
							  TopPos = Top;
							  histM_TopPos->Fill(TopPos.M());
							  histM_WPos->Fill(WPos.M());
						  }
						  if(Muon1_Charge < 0){
							  bJet_Neg = bJet;
							  WNeg = W;
							  TopNeg = Top;
							  histM_TopNeg->Fill(TopNeg.M());
							  histM_WNeg->Fill(WNeg.M());
						  }
					   }
					}
				}
			}
			
			// hadronic W
			for(int jet1 = 0; jet1 < NJet; jet1++){
				if(Jet_ID[jet1] && (Jet_btag[jet1] < Jet_btag_Bound)){
					qJet1_temp.SetPxPyPzE(Jet_Px[jet1], Jet_Py[jet1], Jet_Pz[jet1], Jet_E[jet1]);
						for(int jet2 = 0; jet2 < NJet; jet2++){
							if(Jet_ID[jet2] && (Jet_btag[jet2] < Jet_btag_Bound)){
								qJet2_temp.SetPxPyPzE(Jet_Px[jet2], Jet_Py[jet2], Jet_Pz[jet2], Jet_E[jet2]);
								if(qJet2_temp != qJet1_temp){
									qJet1.SetPxPyPzE(Jet_Px[jet1], Jet_Py[jet1], Jet_Pz[jet1], Jet_E[jet1]);
									qJet2.SetPxPyPzE(Jet_Px[jet2], Jet_Py[jet2], Jet_Pz[jet2], Jet_E[jet2]);
									W_Had = qJet1 + qJet2;
									
									//histM_qJet->Fill(qJet1.M());
									//histM_qJet->Fill(qJet2.M());
									//histM_WHad->Fill(W_Had.M());
								}
							}
						}
				}
			}
		
			// Hadronic W and Top //
			++TerminalCounter;
			for(int jet = 0; jet < NJet; jet++){
				++FirstCounter;
				if(Jet_ID[jet] && (Jet_btag[jet] > Jet_btag_Bound)){
					bJet_Had.SetPxPyPzE(Jet_Px[jet], Jet_Py[jet], Jet_Pz[jet], Jet_E[jet]);
					for(int jet1 = jet; jet1 < NJet; jet1++){
						++SecondCounter;
						qJet1_temp.SetPxPyPzE(Jet_Px[jet1], Jet_Py[jet1], Jet_Pz[jet1], Jet_E[jet1]);
						if(qJet1_temp != bJet_Had && !(Jet_ID[jet1] && (Jet_btag[jet1] > Jet_btag_Bound))){
							for(int jet2 = jet1; jet2 < NJet; jet2++){
								++ThirdCounter;
								qJet2_temp.SetPxPyPzE(Jet_Px[jet2], Jet_Py[jet2], Jet_Pz[jet2], Jet_E[jet2]);
								if(qJet2_temp != qJet1_temp && qJet2_temp != bJet_Had && !(Jet_ID[jet1] && (Jet_btag[jet1] > Jet_btag_Bound)) && !(Jet_ID[jet2] && (Jet_btag[jet2] > Jet_btag_Bound))){
									++FourthCounter;
									qJet1.SetPxPyPzE(Jet_Px[jet1], Jet_Py[jet1], Jet_Pz[jet1], Jet_E[jet1]);
									qJet2.SetPxPyPzE(Jet_Px[jet2], Jet_Py[jet2], Jet_Pz[jet2], Jet_E[jet2]);
									W_Had = qJet1 + qJet2;
									Top_Had = W_Had + bJet_Had;
									//histM_WHad->Fill(W_Had.M());
									histM_WHadVsM_TopHad_NoCut->Fill(Top_Had.M(),W_Had.M());
									
									//histM_qJet->Fill(qJet1.M());
									//histM_qJet->Fill(qJet2.M());
									//if(((qJet1.M() > 3) && (qJet1.M() < 7)) && ((qJet2.M() > 7) && (qJet2.M() < 11))){
										histM_WHad->Fill(W_Had.M());
										histM_qJet->Fill(qJet1.M());
										histM_qJet->Fill(qJet2.M());
									//}
									if((W_Had.M() > W_Mass_LowerBound) && ((W_Had.M() < W_Mass_UpperBound))){
										//histM_WHad->Fill(W_Had.M());
										Top_Had = W_Had + bJet_Had;
										histM_TopHad->Fill(Top_Had.M());
										histM_WHadVsM_TopHad->Fill(Top_Had.M(),W_Had.M());
									}	
								}
							}
						}
					}
					
				}
			}
			
				  // cout << "TerminalCounter: " << TerminalCounter << endl;
				  // cout << "FirstCounter: " << FirstCounter << endl;
				  // cout << "SecondCounter: " << SecondCounter << endl;
				  // cout << "ThirdCounter: " << ThirdCounter << endl;
				  // cout << "FourthCounter: " << FourthCounter << endl;
				  // cout << "******************************************** " <<  endl; 
		}
	}
	

		
		
		//For the WElectron//
		// if(NElectron == 1){
			// cout<<"NElectron: "<<NElectron<<endl;
			// histNElectron->Fill(NElectron);
			// for(int j = 0; j < NElectron; j++){
				// if(Electron_Iso[j] < ElectronIsoBound_ForW){
					// histWElectronIso->Fill(Electron_Iso[j]);
					// cout<<"W| i: "<<i<<" j: "<<j<<" Electron_Charge: "<<Electron_Charge[j]<<endl;
					// WElectron_prime.SetPxPyPzE(Electron_Px[j], Electron_Py[j], Electron_Pz[j], Electron_E[j]);
						
					// if(Electron_Charge[j] > 0){
						// cout<<"W+| i: "<<i<<" j: "<<j<<" Electron_Charge: "<<Electron_Charge[j]<<endl;
						// WElectron_prime_Pos.SetPxPyPzE(Electron_Px[j], Electron_Py[j], Electron_Pz[j], Electron_E[j]);
						
					// }
					// if(Electron_Charge[j] < 0){
						// cout<<"W-| i: "<<i<<" j: "<<j<<" Electron_Charge: "<<Electron_Charge[j]<<endl;
						// WElectron_prime_Neg.SetPxPyPzE(Electron_Px[j], Electron_Py[j], Electron_Pz[j], Electron_E[j]);
						
					// }
				// }
			// }
			// double MET_prime = sqrt((MET_px*MET_px) + (MET_py*MET_py));
			// if(MET_px != 0){
				// double MET_prime_Phi_Reference = abs(atan(MET_py/MET_px));
			// }
			// double MET_prime_Phi;
			// if(MET_px > 0 && MET_py > 0){
				// MET_prime_Phi = MET_prime_Phi_Reference;
			// }
			// if(MET_px > 0 && MET_py < 0){
				// MET_prime_Phi = Pi + MET_prime_Phi_Reference;
			// }
			// if(MET_px < 0 && MET_py > 0){
				// MET_prime_Phi = Pi - MET_prime_Phi_Reference;
			// }
			// if(MET_px < 0 && MET_py < 0){
				// MET_prime_Phi = 2*Pi - MET_prime_Phi_Reference;
			// }
			// if(MET_px == 0 && MET_py > 0){
				// MET_prime_Phi = Pi/2.;
			// }
			// if(MET_px == 0 && MET_py < 0){
				// MET_prime_Phi = (3./2.)*Pi;
			// }
			// if(MET_px > 0 && MET_py == 0){
				// MET_prime_Phi = 0;
			// }
			// if(MET_px < 0 && MET_py == 0){
				// MET_prime_Phi = Pi;
			// }
			// cout<<MET_prime_Phi<<endl;
			// histMETxy->Fill(MET_px, MET_py);
			// histMET->Fill(MET_prime);
			// histMET_Phi->Fill(MET_prime_Phi * (1./Deg2Rad));
			
			// double WElectron_prime_Phi = WElectron_prime.Phi();
			// double WElectron_prime_Pt = WElectron_prime.Pt();
			// double WElectron_prime_Et = WElectron_prime.Et();
			
			// double WElectron_prime_Phi_Pos = WElectron_prime_Pos.Phi();
			// double WElectron_prime_Pt_Pos = WElectron_prime_Pos.Pt();
			// double WElectron_prime_Et_Pos = WElectron_prime_Pos.Et();
			
			// double WElectron_prime_Phi_Neg = WElectron_prime_Neg.Phi();
			// double WElectron_prime_Pt_Neg = WElectron_prime_Neg.Pt();
			// double WElectron_prime_Et_Neg = WElectron_prime_Neg.Et();
			
			// double mT_mine = sqrt(2.*MET_prime*WElectron_prime_Et*(1-cos(MET_prime_Phi - WElectron_prime_Phi)));
			// double mT_mine_Pos = sqrt(2.*MET_prime*WElectron_prime_Et_Pos*(1-cos(MET_prime_Phi - WElectron_prime_Phi_Pos)));
			// double mT_mine_Neg = sqrt(2.*MET_prime*WElectron_prime_Et_Neg*(1-cos(MET_prime_Phi - WElectron_prime_Phi_Neg)));
			
			// histM_WElectron->Fill(mT_mine);
			// histM_WElectronPos->Fill(mT_mine_Pos);
			// histM_WElectronNeg->Fill(mT_mine_Neg);
		// }
		// if(NJet == 2){
			// cout<<"NJet: "<<NJet<<endl;
			// for(int j = 0; j < NJet; j++){
				// if(Jet_ID[j] && (Jet_btag[j] > Jet_btag_Bound)){
					// cout<<"Found good bJet"<<endl;
					// if(WPosHappened){
						// cout<<"b| i: "<<i<<" j: "<<j<<endl;
						// bJet_prime_Pos.SetPxPyPzE(Jet_Px[j], Jet_Py[j], Jet_Pz[j], Jet_E[j]);
						// cout<<"Made 4 Vect"<<endl;
					// }
					// if(WNegHappened){
						// cout<<"b-| i: "<<i<<" j: "<<j<<endl;
						// bJet_prime_Neg.SetPxPyPzE(Jet_Px[j], Jet_Py[j], Jet_Pz[j], Jet_E[j]);
						// cout<<"Made 4 Vect"<<endl;
					// }
				// }
				// if((WPosHappened && WNegHappened) && (Jet_ID[j] && (Jet_btag[j] > Jet_btag_Bound))){
					// NbJet = NbJet + 1;
					// hist_bJetVbtag->Fill(Jet_btag[j], NbJet);
					
					// Top_Pos = bJet_prime_Pos + WMuon_prime_Pos;
					// Top_Neg = bJet_prime_Neg + WMuon_prime_Neg;
					// cout<<"Added Vects"<<endl;
					
					// histM_TopPos->Fill(Top_Pos.M());
					// histM_TopNeg->Fill(Top_Neg.M());
					// histM_Top->Fill(Top_Pos.M());
					// histM_Top->Fill(Top_Neg.M());
					// cout<<"Filled hists"<<endl;
				// }
			// }	
		// }
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
	//}
	
	// if(FitInvarientMassHist){
		TF1 *f1 = new TF1("f1","gaus",80,100);
		TF1 *f2 = new TF1("f2","([0]*exp(-((x-[1])^2)/(2*[2]^2)))",60,85);
		
		f2->SetParameters(0, 200); //Normal
		f2->SetParameters(1, 70); //Centroid
		f2->SetParameters(2, 25); //Width
	// }
	
	//////Graphics//////
	ifstream inFile;
	inFile.open("mydata.dat",ios::in);
	const int numOfLines = 52;
	const int numArray = numOfLines + 1;
	double counts[numArray];
	double binCenter[numArray];
	double error[numArray];
	inFile >> binCenter[0] >> counts[0] >> error[0];
	for(int i = 1; i < numOfLines; i++) {
		inFile >> binCenter[i] >> counts[i] >> error[i];
		//cout << binCenter[i] << " " << counts[i] << " " << error[i] << "\n";
		histDataOverlay->Fill(binCenter[i], counts[i]);
	}
	//histDataOverlay->Print("all");
	inFile.close();
		
	
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
			
	TCanvas *c7 = new TCanvas("c7","Transverse Mass Top Histrogram" ,200,10,900,700);
		//c7->Divide(1,3);
		c7->SetFillColor(10);
		//c7->SetGrid();
		c7->GetFrame()->SetFillColor(10);
		c7->GetFrame()->SetBorderSize(12);
		c7->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			//c7->cd(1);
			gPad->SetLogy();
			histM_Top->SetXTitle("Invarient Mass");
			histM_Top->SetYTitle("Counts");
			histM_Top->SetLineColor(kRed);
			histM_Top->SetFillColor(kRed);
			histM_Top->SetStats(1);
			histM_Top->SetMaximum(mT_MaxCount);
			histM_Top->Draw();
			
			gPad->SetLogy();
			histDataOverlay->SetXTitle("Invarient Mass");
			histDataOverlay->SetYTitle("Counts");
			histDataOverlay->SetLineColor(kBlue);
			histDataOverlay->SetStats(1);
			histDataOverlay->SetMaximum(mT_MaxCount);
			
			//TGraphErrors *TGAE1 = new TGraphErrors(50, binCenter, counts, 0, error);
			histDataOverlay->SetMarkerStyle(20);
			//histDataOverlay->SetMarkerColor(kRed);
			//histDataOverlay->SetMarkerStyle(20);
			histDataOverlay->Draw("psame");
			
			
			
			/*
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
			*/
			
			if(SaveAsRootFile){
				TFile outfile7(rootOutputFileName7, "RECREATE");
				c7->Write(rootOutputFileName7);
				outfile7.Close();
			}
		
	TCanvas *c8 = new TCanvas("c8","Invariant Mass of W Histrogram" ,200,10,900,700);
		c8->Divide(2,2);
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
			
			c8->cd(4);
			histWMuonPt->SetXTitle("WMuon Pt");
			histWMuonPt->SetYTitle("Counts");
			histWMuonPt->SetStats(1);
			//histM_WNeg->SetMaximum(mT_MaxCount);
			histWMuonPt->Draw();
			
			if(SaveAsRootFile){
				TFile outfile8(rootOutputFileName8, "RECREATE");
				c8->Write(rootOutputFileName8);
				outfile8.Close();
			}
		
	TCanvas *c9 = new TCanvas("c9","Muon Iso vs Pt" ,200,10,900,700);
		c9->Divide(1,2);
		c9->SetFillColor(10);
		c9->SetGrid();
		c9->GetFrame()->SetFillColor(10);
		c9->GetFrame()->SetBorderSize(12);
		c9->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c9->cd(1);
			histMuonIsoVsPt->SetXTitle("Pt");
			histMuonIsoVsPt->SetYTitle("Muon Iso");
			histMuonIsoVsPt->SetStats(1);
			//histM_W->SetMaximum(mT_MaxCount);
			histMuonIsoVsPt->Draw();
			
			c9->cd(2);
			histMuonIsoVsPt_AfterCut->SetXTitle("Pt");
			histMuonIsoVsPt_AfterCut->SetYTitle("Muon Iso");
			histMuonIsoVsPt_AfterCut->SetStats(1);
			//histM_WPos->SetMaximum(mT_MaxCount);
			histMuonIsoVsPt_AfterCut->Draw();
			
			/*
			c9->cd(3);
			histM_WNeg->SetXTitle("Invarient Mass");
			histM_WNeg->SetYTitle("Counts");
			histM_WNeg->SetStats(1);
			//histM_WNeg->SetMaximum(mT_MaxCount);
			histM_WNeg->Draw();
			
			c9->cd(4);
			histWMuonPt->SetXTitle("WMuon Pt");
			histWMuonPt->SetYTitle("Counts");
			histWMuonPt->SetStats(1);
			//histM_WNeg->SetMaximum(mT_MaxCount);
			histWMuonPt->Draw();
			*/
			
			if(SaveAsRootFile){
				TFile outfile9(rootOutputFileName9, "RECREATE");
				c9->Write(rootOutputFileName9);
				outfile9.Close();
			}
		
		
	TCanvas *c10 = new TCanvas("c10","InvarientMass Z Histrogram" ,200,10,900,700);
		c10->Divide(1,2);
		c10->SetFillColor(10);
		c10->SetGrid();
		c10->GetFrame()->SetFillColor(10);
		c10->GetFrame()->SetBorderSize(12);
		c10->Range(0,0,1,1);

					
		gStyle->SetOptStat(1);
		if(FitInvarientMassHist){
			gStyle->SetOptFit(1);
		}
			
			c10->cd(1);
			gPad->SetLogy();
			histM_Muon->SetXTitle("Invarient Mass");
			histM_Muon->SetYTitle("Counts");
			histM_Muon->SetStats(1);
			if(FitInvarientMassHist){
				histM_Muon->Fit(f1, "M");
			}
			histM_Muon->Draw();
			
			c10->cd(2);
			gPad->SetLogy();
			histM_Z->SetXTitle("Invarient Mass (Solution)");
			histM_Z->SetYTitle("Counts");
			histM_Z->SetStats(1);
			if(FitInvarientMassHist){
				histM_Z->Fit(f1, "M");
			}
			histM_Z->Draw();
			
			// c10->cd(3);
			// histM_WMuon2->SetXTitle("Mt 2 Mass");
			// histM_WMuon2->SetYTitle("Counts");
			// histM_WMuon2->SetStats(1);
			// histM_WMuon2->Draw();
			
			// c10->cd(4);
			// histM_WMuon3->SetXTitle("Mt 3 Mass");
			// histM_WMuon3->SetYTitle("Counts");
			// histM_WMuon3->SetStats(1);
			// histM_WMuon3->Draw();
			
			if(SaveAsRootFile){
				TFile outfile10(rootOutputFileName10, "RECREATE");
				c10->Write(rootOutputFileName10);
				outfile10.Close();
			}
			
	TCanvas *c11 = new TCanvas("c11","Transverse Mass Top Comparison Histrogram" ,200,10,900,700);
		c11->Divide(2,2);
		c11->SetFillColor(10);
		c11->SetGrid();
		c11->GetFrame()->SetFillColor(10);
		c11->GetFrame()->SetBorderSize(12);
		c11->Range(0,0,1,1);

					
		gStyle->SetOptStat(1);
		if(FitInvarientMassHist){
			gStyle->SetOptFit(1);
		}
			
			c11->cd(1);
			gPad->SetLogy();
			histM_Top_ZCut->SetXTitle("Invarient Total Top Mass");
			histM_Top_ZCut->SetYTitle("Counts");
			//histM_Top_ZCut->SetLineColor(kBlue);
			histM_Top_ZCut->SetStats(1);
			histM_Top_ZCut->SetMaximum(200);
			histM_Top_ZCut->Draw();
			
			c11->cd(2);
			gPad->SetLogy();
			histM_TopPos->SetXTitle("Invarient TopPos Mass");
			histM_TopPos->SetYTitle("Counts");
			histM_TopPos->SetStats(1);
			histM_TopPos->SetMaximum(200);
			histM_TopPos->Draw();
			
			c11->cd(3);
			gPad->SetLogy();
			histM_TopNeg->SetXTitle("Invarient TopNeg Mass");
			histM_TopNeg->SetYTitle("Counts");
			histM_TopNeg->SetStats(1);
			histM_TopNeg->SetMaximum(200);
			histM_TopNeg->Draw();
			
			c11->cd(4);
			gPad->SetLogy();
			histM_TopHad->SetXTitle("Invarient TopHad Mass");
			histM_TopHad->SetYTitle("Counts");
			histM_TopHad->SetStats(1);
			histM_TopHad->SetMaximum(200);
			histM_TopHad->Draw();
			
			if(SaveAsRootFile){
				TFile outfile11(rootOutputFileName11, "RECREATE");
				c11->Write(rootOutputFileName11);
				outfile11.Close();
			}
		
	TCanvas *c12 = new TCanvas("c12","Invariant Mass of W Comparison Histrogram" ,200,10,900,700);
		c12->Divide(1,3);
		c12->SetFillColor(10);
		c12->SetGrid();
		c12->GetFrame()->SetFillColor(10);
		c12->GetFrame()->SetBorderSize(12);
		c12->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c12->cd(1);
			histM_W->SetXTitle("Invarient Mass");
			histM_W->SetYTitle("Counts");
			histM_W->SetStats(1);
			//histM_W->SetMaximum(mT_MaxCount);
			histM_W->Draw();
			
			c12->cd(2);
			gPad->SetLogy();
			histM_WHad->SetXTitle("Invarient Mass (Hadronic)");
			histM_WHad->SetYTitle("Counts");
			histM_WHad->SetStats(1);
			//histM_WHad->SetMaximum(mT_MaxCount);
			histM_WHad->Draw();
			
			c12->cd(3);
			histM_Wt->SetXTitle("Transverse Mass");
			histM_Wt->SetYTitle("Counts");
			histM_Wt->SetStats(1);
			//histM_Wt->SetMaximum(mT_MaxCount);
			histM_Wt->Draw();
			
			// c12->cd(4);
			// histWMuonPt->SetXTitle("WMuon Pt");
			// histWMuonPt->SetYTitle("Counts");
			// histWMuonPt->SetStats(1);
			// histM_WNeg->SetMaximum(mT_MaxCount);
			// histWMuonPt->Draw();
			
			if(SaveAsRootFile){
				TFile outfile12(rootOutputFileName12, "RECREATE");
				c12->Write(rootOutputFileName12);
				outfile12.Close();
			}
			
	TCanvas *c13 = new TCanvas("c13","qJet Histogram" ,200,10,900,700);
		//c13->Divide(1,3);
		c13->SetFillColor(10);
		c13->SetGrid();
		c13->GetFrame()->SetFillColor(10);
		c13->GetFrame()->SetBorderSize(12);
		c13->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c13->cd(1);
			gPad->SetLogy();
			histM_qJet->SetXTitle("Invarient Mass");
			histM_qJet->SetYTitle("Counts");
			histM_qJet->SetStats(1);
			//histM_qJet->SetMaximum(mT_MaxCount);
			histM_qJet->Draw();
			

	
			
			if(SaveAsRootFile){
				TFile outfile13(rootOutputFileName13, "RECREATE");
				c13->Write(rootOutputFileName13);
				outfile13.Close();
			}
			
	TCanvas *c14 = new TCanvas("c14","Invariant Mass of W vs Invariant Mass of Top" ,200,10,900,700);
	c14->Divide(1,2);
	c14->SetFillColor(10);
	c14->SetGrid();
	c14->GetFrame()->SetFillColor(10);
	c14->GetFrame()->SetBorderSize(12);
	c14->Range(0,0,1,1);

	gStyle->SetOptStat(1);
		
		c14->cd(1);
		histM_WHadVsM_TopHad->SetXTitle("Invariant Mass of the Top");
		histM_WHadVsM_TopHad->SetYTitle("Invariant Mass of the W");
		histM_WHadVsM_TopHad->SetStats(1);
		//histM_WHadVsM_TopHad->SetMarkerSize(200);
		histM_WHadVsM_TopHad->SetMarkerStyle(20);
		//histM_WHadVsM_TopHad->SetMaximum(mT_MaxCount);
		histM_WHadVsM_TopHad->Draw();
		
		c14->cd(2);
		histM_WHadVsM_TopHad_NoCut->SetXTitle("Invariant Mass of the Top");
		histM_WHadVsM_TopHad_NoCut->SetYTitle("Invariant Mass of the W");
		histM_WHadVsM_TopHad_NoCut->SetStats(1);
		//histM_WHadVsM_TopHad->SetMarkerSize(200);
		histM_WHadVsM_TopHad_NoCut->SetMarkerStyle(20);
		//histM_WHadVsM_TopHad->SetMaximum(mT_MaxCount);
		histM_WHadVsM_TopHad_NoCut->Draw();

		
		if(SaveAsRootFile){
			TFile outfile14(rootOutputFileName14, "RECREATE");
			c14->Write(rootOutputFileName14);
			outfile14.Close();
		}
			
	return 0;
}