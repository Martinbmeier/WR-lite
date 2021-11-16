//C++ CLASSES
#include <iostream>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
//LOCAL CLASSES
#include "ExoAnalysis/WR_lite/interface/eventBits2.h"
#include "ExoAnalysis/WR_lite/interface/cutFlowHistos.h"

#include "TH1D.h"
#include "TStyle.h"



cutFlowHistos::cutFlowHistos () {}
  
//void cutFlowHistos::book(TFileDirectory histoFolder) {
void cutFlowHistos::book(TFileDirectory histoFolder, int nCut) {

	gStyle->SetOptStat("omen");
	m_histoFolder=histoFolder;

	//Unmatchded lepton histos


	//Muon Histos

	m_recoMuonPt[nCut] =  {m_histoFolder.make<TH1D>("recoMuonPt","Pt for reco muons",100,0,1000)};
	m_recoMuonPt[nCut]->GetXaxis()-> SetTitle("Pt (GeV)");

	m_invMass[nCut] =  {m_histoFolder.make<TH1D>("4objectMass","invariant mass",100,0,1000)};
	m_invMass[nCut]->GetXaxis()-> SetTitle("mass (GeV/c)");

}

//General histogram filling
void cutFlowHistos::fill(double pT, double mass, int cutNumber, double weight) {
	

		m_recoMuonPt[cutNumber]->Fill(pT,weight);
		m_invMass[cutNumber]->Fill(mass,weight);



}


	
