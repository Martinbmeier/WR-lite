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

	m_invMassAS_mu_mu[nCut] =  {m_histoFolder.make<TH1D>("4objectMassAS_mu_mu","invariant mass (misID'ed)",80,0,8000)};
	m_invMassAS_mu_mu[nCut]->GetXaxis()-> SetTitle("mass (GeV)");

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

	m_leadmuon_pt[nCut] =  {m_histoFolder.make<TH1D>("leadmuon_pt","leading muon pT",50,0,800)};
	m_leadmuon_pt[nCut]->GetXaxis()-> SetTitle("pT (GeV)");

	m_subleadmuon_pt[nCut] =  {m_histoFolder.make<TH1D>("subleadmuon_pt","subleading muon pT",50,0,800)};
	m_subleadmuon_pt[nCut]->GetXaxis()-> SetTitle("pT (GeV)");

	m_leadsubleadmuon_pt[nCut] =  {m_histoFolder.make<TH2D>("leadsubleadmuon_pt","leading muon pT vs. subleading muon pT",50,0,800,50,0,800)};
	m_leadsubleadmuon_pt[nCut]->GetXaxis()-> SetTitle("leading pT (GeV)");
	m_leadsubleadmuon_pt[nCut]->GetYaxis()-> SetTitle("subleading pT (GeV)");

}

//General histogram filling
void cutFlowHistos::fill(double pT1, double pT2, double ASmass_mu_mu, double SSmass_mu_mu, double OSmass_mu_mu, double SSmass_e_e, double OSmass_e_e, double SSmass_mu_e, double OSmass_mu_e,  int cutNumber, double weight) {
	
		m_recoMuonPt[cutNumber]->Fill(pT1,weight);
		m_recoMuonPt[cutNumber]->Fill(pT2,weight);
		m_invMassAS_mu_mu[cutNumber]->Fill(ASmass_mu_mu,weight);
		m_invMassSS_mu_mu[cutNumber]->Fill(SSmass_mu_mu,weight);
		m_invMassOS_mu_mu[cutNumber]->Fill(OSmass_mu_mu,weight);
		m_invMassSS_e_e[cutNumber]->Fill(SSmass_e_e,weight);
		m_invMassOS_e_e[cutNumber]->Fill(OSmass_e_e,weight);
		m_invMassSS_mu_e[cutNumber]->Fill(SSmass_mu_e,weight);
		m_invMassOS_mu_e[cutNumber]->Fill(OSmass_mu_e,weight);

		m_leadmuon_pt[cutNumber]->Fill(pT1,weight);
		m_subleadmuon_pt[cutNumber]->Fill(pT2,weight);
		m_leadsubleadmuon_pt[cutNumber]->Fill(pT1,pT2,weight);

}


	
