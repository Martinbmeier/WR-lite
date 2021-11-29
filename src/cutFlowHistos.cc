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

	m_recoMuonPt[nCut] =  {m_histoFolder.make<TH1D>("recoMuonPt","Pt for reco muons",100,0,1000)};
	m_recoMuonPt[nCut]->GetXaxis()-> SetTitle("Pt (GeV)");

	m_invMassSS_mu_mu[nCut] =  {m_histoFolder.make<TH1D>("4objectMassSS_mu_mu","invariant mass (same sign, 2 muons)",80,0,8000)};
	m_invMassSS_mu_mu[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

	m_invMassOS_mu_mu[nCut] =  {m_histoFolder.make<TH1D>("4objectMassOS_mu_mu","invariant mass (opposite sign, 2 muons)",80,0,8000)};
	m_invMassOS_mu_mu[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

	m_invMassSS_e_e[nCut] =  {m_histoFolder.make<TH1D>("4objectMassSS_e_e","invariant mass (same sign, 2 electrons)",80,0,8000)};
	m_invMassSS_e_e[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

	m_invMassOS_e_e[nCut] =  {m_histoFolder.make<TH1D>("4objectMassOS_e_e","invariant mass (opposite sign, 2 electrons)",80,0,8000)};
	m_invMassOS_e_e[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

	m_invMassSS_mu_e[nCut] =  {m_histoFolder.make<TH1D>("4objectMassSS_mu_e","invariant mass (same sign, muon/electron)",80,0,8000)};
	m_invMassSS_mu_e[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

	m_invMassOS_mu_e[nCut] =  {m_histoFolder.make<TH1D>("4objectMassOS_mu_e","invariant mass (opposite sign, muon/electron)",80,0,8000)};
	m_invMassOS_mu_e[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

}

//General histogram filling
void cutFlowHistos::fill(double pT, double SSmass_mu_mu, double OSmass_mu_mu, double SSmass_e_e, double OSmass_e_e, double SSmass_mu_e, double OSmass_mu_e,  int cutNumber, double weight) {
	
		m_recoMuonPt[cutNumber]->Fill(pT,weight);
		m_invMassSS_mu_mu[cutNumber]->Fill(SSmass_mu_mu,weight);
		m_invMassOS_mu_mu[cutNumber]->Fill(OSmass_mu_mu,weight);
		m_invMassSS_e_e[cutNumber]->Fill(SSmass_e_e,weight);
		m_invMassOS_e_e[cutNumber]->Fill(OSmass_e_e,weight);
		m_invMassSS_mu_e[cutNumber]->Fill(SSmass_mu_e,weight);
		m_invMassOS_mu_e[cutNumber]->Fill(OSmass_mu_e,weight);

}


	
