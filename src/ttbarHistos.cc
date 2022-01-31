//C++ CLASSES
#include <iostream>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
//LOCAL CLASSES
#include "ExoAnalysis/WR_lite/interface/eventBits2.h"
#include "ExoAnalysis/WR_lite/interface/ttbarHistos.h"

#include "TH1D.h"
#include "TStyle.h"



ttbarHistos::ttbarHistos () {}
  
//void cutFlowHistos::book(TFileDirectory histoFolder) {
void ttbarHistos::book(TFileDirectory histoFolder) {

	gStyle->SetOptStat("omen");
	m_histoFolder=histoFolder;

	m_electronpT = m_histoFolder.make<TH1D>("electronpT","gen electron pT",100,0,1000);
	m_electronpT.GetXaxis()->SetTitle("pT (GeV)");

	m_electronMuonpT = m_histoFolder.make<TH2D>("electronmuonpt","electron pT vs. muon pT",50,0,800,50,0,800);
	m_electronMuonpT.GetXaxis()->SetTitle("gen muon pT (GeV)");
	m_electronMuonpT.GetYaxis()->SetTitle("genm electron pT (GeV)");

}

//General histogram filling
void ttbarHistos::fill(double mupT, double epT, double weight) {
	
		m_electronpT.Fill(epT,weight);
		m_electronMuonpT.Fill(mupT,epT,weight);

}


	
