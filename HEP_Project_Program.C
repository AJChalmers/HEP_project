//Program to read CERN data from trees and make hists. (So far)


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


using namespace std;
int main(){
	
	/////Controls/////

	//int NEntries = 469384; //469384 for data.root. Click on tree and read hist to find these.
	
	char rootFilePrefix[50] = "HEPTutorial/files/"; //Prefix to root file, this is the path to it - all the way up until its name.
	char rootFileName[50] = "wjets.root"; //The root file name itself - this is where the tree is.
	
	double xMin = -400;
	double xMax = 400;
	
	// double yMin = 0;
	// double yMax = 10000;
	// double zMin = 0;
	// double zMax = 10000;
	
	const int NChannels = 100; //Number of channels for root hist
	
	bool OutputToFile = false;
	
	///////////////////////////////////////
	
	/////Varibles (Pre-stated)/////
	
	char rootFilePath[100]; 
	
	int NJet, NMuon, Muon_Charge, NElectron, Electron_Charge, NPhoton, MCleptonPDGid, NPrimaryVertices;
	
	float Jet_Px, Jet_Py, Jet_Pz, Jet_E, Jet_btag, Muon_Px, Muon_Py, Muon_Pz, Muon_E, Muon_Iso, 
	Electron_Px, Electron_Py, Electron_Pz, Electron_E, Electron_Iso, Photon_Px, Photon_Py, Photon_Pz, Photon_E, Photon_Iso,
	MET_px, MET_py, MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz,
	MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz,
	MClepton_px, MClepton_py, MClepton_pz, MCneutrino_px, MCneutrino_py, MCneutrino_pz, EventWeight;
	
	bool Jet_ID, triggerIsoMu24;
	
	ofstream outFile1;
	
	//////////////////////////////////////
	
	//Creating path to root file and getting tree//
	strcpy(rootFilePath, rootFilePrefix);
	strcat(rootFilePath, rootFileName);
	cout<<"Root file path: "<<rootFilePath<<endl;
	
	TFile *file = new TFile(rootFilePath);
	TTree *t1 = (TTree*)file->Get("events");
	
	t1->SetBranchAddress("NJet",&NJet);
	t1->SetBranchAddress("Jet_Px",&Jet_Px);
	t1->SetBranchAddress("Jet_Py",&Jet_Py);
	t1->SetBranchAddress("Jet_Pz",&Jet_Pz);
	t1->SetBranchAddress("Jet_E",&Jet_E);
	t1->SetBranchAddress("Jet_btag",&Jet_btag);
	t1->SetBranchAddress("Jet_ID",&Jet_ID);
	t1->SetBranchAddress("NMuon",&NMuon);
	t1->SetBranchAddress("Muon_Px",&Muon_Px);
	t1->SetBranchAddress("Muon_Py",&Muon_Py);
	t1->SetBranchAddress("Muon_Pz",&Muon_Pz);
	t1->SetBranchAddress("Muon_E",&Muon_E);
	t1->SetBranchAddress("Muon_Charge",&Muon_Charge);
	t1->SetBranchAddress("Muon_Iso",&Muon_Iso);
	t1->SetBranchAddress("NElectron",&NElectron);
	t1->SetBranchAddress("Electron_Px",&Electron_Px);
	t1->SetBranchAddress("Electron_Py",&Electron_Py);
	t1->SetBranchAddress("Electron_Pz",&Electron_Pz);
	t1->SetBranchAddress("Electron_E",&Electron_E);
	t1->SetBranchAddress("Electron_Charge",&Electron_Charge);
	t1->SetBranchAddress("Electron_Iso",&Electron_Iso);
	t1->SetBranchAddress("NPhoton",&NPhoton);
	t1->SetBranchAddress("Photon_Px",&Photon_Px);
	t1->SetBranchAddress("Photon_Py",&Photon_Py);
	t1->SetBranchAddress("Photon_Pz",&Photon_Pz);
	t1->SetBranchAddress("Photon_E",&Photon_E);
	t1->SetBranchAddress("Photon_Iso",&Photon_Iso);
	t1->SetBranchAddress("MET_px",&MET_px);
	t1->SetBranchAddress("MET_py",&MET_py);
	t1->SetBranchAddress("MChadronicBottom_px",&MChadronicBottom_px);
	t1->SetBranchAddress("MChadronicBottom_py",&MChadronicBottom_py);
	t1->SetBranchAddress("MChadronicBottom_pz",&MChadronicBottom_pz);
	t1->SetBranchAddress("MCleptonicBottom_px",&MCleptonicBottom_px);
	t1->SetBranchAddress("MCleptonicBottom_py",&MCleptonicBottom_py);
	t1->SetBranchAddress("MCleptonicBottom_pz",&MCleptonicBottom_pz);
	t1->SetBranchAddress("MChadronicWDecayQuark_px",&MChadronicWDecayQuark_px);
	t1->SetBranchAddress("MChadronicWDecayQuark_py",&MChadronicWDecayQuark_py);
	t1->SetBranchAddress("MChadronicWDecayQuark_pz",&MChadronicWDecayQuark_pz);
	t1->SetBranchAddress("MChadronicWDecayQuarkBar_px",&MChadronicWDecayQuarkBar_px);
	t1->SetBranchAddress("MChadronicWDecayQuarkBar_py",&MChadronicWDecayQuarkBar_py);
	t1->SetBranchAddress("MChadronicWDecayQuarkBar_pz",&MChadronicWDecayQuarkBar_pz);
	t1->SetBranchAddress("MClepton_px",&MClepton_px);
	t1->SetBranchAddress("MClepton_py",&MClepton_py);
	t1->SetBranchAddress("MClepton_pz",&MClepton_pz);
	t1->SetBranchAddress("MCleptonPDGid",&MCleptonPDGid);
	t1->SetBranchAddress("MCneutrino_px",&MCneutrino_px);
	t1->SetBranchAddress("MCneutrino_py",&MCneutrino_py);
	t1->SetBranchAddress("MCneutrino_pz",&MCneutrino_pz);
	t1->SetBranchAddress("NPrimaryVertices",&NPrimaryVertices);
	t1->SetBranchAddress("triggerIsoMu24",&triggerIsoMu24);
	t1->SetBranchAddress("EventWeight",&EventWeight);
	
	int NEntries = t1->GetEntries();
	cout<<"NEntries: "<<NEntries<<endl;
	
	//Loop to get entries and fill hists//
	
	TH1F *histJetPx = new TH1F("histJetPx", "JetPx Hist", NChannels, xMin, xMax);
	TH1F *histJetPy = new TH1F("histJetPy", "JetPy Hist", NChannels, xMin, xMax);
	TH1F *histJetPz = new TH1F("histJetPz", "JetPz Hist", NChannels, xMin, xMax);
	//TH3F *histJets = new TH1F("histJet", HistTitle, NChannels, xMin, xMax, yMin, yMax, zMin, zMax);
	
	if(OutputToFile){
		outFile1.open("jetfile.dat",ios::out);
		
		for(int i = 0; i < NEntries; i++){
			t1->GetEntry(i);
			histJetPx->Fill(Jet_Px);
			histJetPy->Fill(Jet_Py);
			histJetPz->Fill(Jet_Pz);
			
			outFile1<<i<<setw(10)<<Jet_Px<<setw(10)<<Jet_Py<<setw(10)<<Jet_Pz<<endl;
		}
		
		outFile1.close();
	}
	
	if(!OutputToFile){
		for(int i = 0; i < NEntries; i++){
			t1->GetEntry(i);
			histJetPx->Fill(Jet_Px);
			histJetPy->Fill(Jet_Py);
			histJetPz->Fill(Jet_Pz);
			
		}
	}
	
	
	
	//////Graphics//////
	cout<<"Making histrograms..."<<endl;
	
	TCanvas *c1 = new TCanvas("c1","Jet_Px Histrogram" ,200,10,900,700);
		c1->SetFillColor(10);
		c1->SetGrid();
		c1->GetFrame()->SetFillColor(10);
		c1->GetFrame()->SetBorderSize(12);
		c1->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c1->cd();
			histJetPx->SetXTitle("Jet_Px");
			histJetPx->SetYTitle("Counts");
			histJetPx->SetStats(1);
			histJetPx->Draw();
			
			
	TCanvas *c2 = new TCanvas("c2","Jet_Py Histrogram" ,200,10,900,700);
		c2->SetFillColor(10);
		c2->SetGrid();
		c2->GetFrame()->SetFillColor(10);
		c2->GetFrame()->SetBorderSize(12);
		c2->Range(0,0,1,1);

		gStyle->SetOptStat(1);
			
			c2->cd();
			histJetPy->SetXTitle("Jet_Py");
			histJetPy->SetYTitle("Counts");
			histJetPy->SetStats(1);
			histJetPy->Draw();


	TCanvas *c3 = new TCanvas("c3","Jet_Pz Histrogram" ,200,10,900,700);
		c3->SetFillColor(10);
		c3->SetGrid();
		c3->GetFrame()->SetFillColor(10);
		c3->GetFrame()->SetBorderSize(12);
		c3->Range(0,0,1,1);

		gStyle->SetOptStat(1);
		
			c3->cd();
			histJetPz->SetXTitle("Jet_Pz");
			histJetPz->SetYTitle("Counts");
			histJetPz->SetStats(1);
			histJetPz->Draw();



	return 0;
}