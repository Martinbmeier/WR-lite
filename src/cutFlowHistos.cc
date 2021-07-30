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

	m_recoMuonPt[nCut] =  {m_histoFolder.make<TH1D>("recoMuonPt","Pt for reco muons",70,0,700)};
	m_recoMuonPt[nCut]->GetXaxis()-> SetTitle("Pt (GeV)");

}

//General histogram filling
void cutFlowHistos::fill(double pT, int cutNumber, double weight) {
	

		m_recoMuonPt[cutNumber]->Fill(pT, weight);



}


	
