//C++ CLASSES
#include <iostream>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
//LOCAL CLASSES
#include "ExoAnalysis/WR_lite/interface/eventBits2.h"
#include "ExoAnalysis/WR_lite/interface/eventHistos2.h"

#include "TH1D.h"
#include "TStyle.h"



eventHistos2::eventHistos2 () {}
  
//void eventHistos2::book(TFileDirectory histoFolder) {
void eventHistos2::book(TFileDirectory histoFolder, int nDirectory) {

	gStyle->SetOptStat("omen");
	m_histoFolder=histoFolder;

	//m_histoFolder.cd();

	//Gen Histos
	
	//1D
	
	m_matchStatus[nDirectory] = {m_histoFolder.make<TH1D>("matchStatus","were there successful gen-reco matches in this event?", 2,1,2)};
	m_matchStatus[nDirectory]->GetXaxis()->SetBinLabel(1,"matched");
	m_matchStatus[nDirectory]->GetXaxis()->SetBinLabel(2,"unmatched");

	m_numRecoMuons[nDirectory] = {m_histoFolder.make<TH1D>("numRecoMuons","num reco", 10,0,10)};
	m_numRecoElectrons[nDirectory] = {m_histoFolder.make<TH1D>("numRecoElectrons","num reco", 10,0,10)};

	m_numGenMuons[nDirectory] = {m_histoFolder.make<TH1D>("numGenMuons","num gen", 10,0,10)};
	m_numGenElectrons[nDirectory] = {m_histoFolder.make<TH1D>("numGenElectrons","num gen", 10,0,10)};

	m_genMinusRecoMuons[nDirectory] = {m_histoFolder.make<TH1D>("genMinusRecoMuons","gen-reco", 10,-5,5)};
	m_genMinusRecoElectrons[nDirectory] = {m_histoFolder.make<TH1D>("genMinusRecoElectrons","gen-reco", 10,-5,5)};

	m_unmatchedPt[nDirectory] = {m_histoFolder.make<TH1D>("unmatchedPt","Pt of unmatched lepton", 100,0,500)};
	m_unmatchedPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	m_unmatchedEta[nDirectory] = {m_histoFolder.make<TH1D>("unmatchedEta","Eta of unmatched lepton", 100,-3.4,3.4)};
	m_unmatchedEta[nDirectory]->GetXaxis()-> SetTitle("eta (radians)");
	m_unmatchedPhi[nDirectory] = {m_histoFolder.make<TH1D>("unmatchedPhi","Phi of unmatched lepton", 100,-3.4,3.4)};
	m_unmatchedPhi[nDirectory]->GetXaxis()-> SetTitle("phi (radians)");
	m_unmatchedDR[nDirectory] = {m_histoFolder.make<TH1D>("unmatchedDR","smallest delta R of unmatched lepton", 100,0,5)};
	gStyle->SetOptStat("omen");

	m_genRecoMuonPtRatio[nDirectory] = {m_histoFolder.make<TH1D>("genRecoMuonPtRatio" , "ratio of gen pt to reco pt for muons" , 50 , 0 , 5 )};
	gStyle->SetOptStat("omen");

	m_genMuonPt[nDirectory] = {m_histoFolder.make<TH1D>("genMuonPt","Pt for gen muons",100,0,500)};
	m_genMuonPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_genMuonPhi[nDirectory] = {m_histoFolder.make<TH1D>("genMuonPhi","Phi for gen muons",100,-3.4,3.4)};
	m_genMuonPhi[nDirectory]->GetXaxis()-> SetTitle("Phi");

	m_genMuonEta[nDirectory] = {m_histoFolder.make<TH1D>("genMuonEta","Eta for gen muons",100,-3.4,3.4)};
	m_genMuonEta[nDirectory]->GetXaxis()-> SetTitle("Eta");

	m_recoMuonPt[nDirectory] =  {m_histoFolder.make<TH1D>("recoMuonPt","Pt for reco muons",100,0,500)};
	m_recoMuonPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");

	m_recoMuonPhi[nDirectory] = {m_histoFolder.make<TH1D>("recoMuonPhi","Phi for reco muons",100,-3.4,3.4)};
	m_recoMuonPhi[nDirectory]->GetXaxis()-> SetTitle("Phi");

	m_recoMuonEta[nDirectory] = {m_histoFolder.make<TH1D>("recoMuonEta","Eta for reco muons",100,-3.4,3.4)};
	m_recoMuonEta[nDirectory]->GetXaxis()-> SetTitle("Eta");
	gStyle->SetOptStat("omen");
	
	m_WrongSignGenRecoMuonPtRatio[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignGenRecoMuonPtRatio" , "ratio of gen pt to reco pt (wrong charge matched muons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_WrongSignRecoMuonPt[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignRecoMuonPt" , "Pt for wrong sign reco muons" , 100 , 0 , 500 )};
	m_WrongSignRecoMuonPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_WrongSignRecoMuonPhi[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignRecoMuonPhi" , "phi for wrong sign muons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoMuonPhi[nDirectory]->GetXaxis()->SetTitle("phi");

	m_WrongSignRecoMuonEta[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignRecoMuonEta" , "Eta for wrong sign muons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoMuonEta[nDirectory]->GetXaxis()->SetTitle("eta");

	m_RightSignGenRecoMuonPtRatio[nDirectory] = {m_histoFolder.make<TH1D>("RightSignGenRecoMuonPtRatio" , "ratio of gen pt to reco pt (right charge matched muons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_RightSignRecoMuonPt[nDirectory] = {m_histoFolder.make<TH1D>("RightSignRecoMuonPt" , "Pt for right sign reco muons" , 100 , 0 , 500 )};
	m_RightSignRecoMuonPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	
	m_RightSignRecoMuonPhi[nDirectory] = {m_histoFolder.make<TH1D>("RightSignRecoMuonPhi" , "phi for right sign muons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoMuonPhi[nDirectory]->GetXaxis()->SetTitle("phi");

	m_RightSignRecoMuonEta[nDirectory] = {m_histoFolder.make<TH1D>("RightSignRecoMuonEta" , "eta for right sign muons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoMuonEta[nDirectory]->GetXaxis()->SetTitle("eta");

	
	m_bestMuondR[nDirectory] = {m_histoFolder.make<TH1D>("bestMuondR" , "smallest deltaR (muons)" , 50 , 0 , 0.03 )};
	m_bestMuondR[nDirectory]->GetXaxis()-> SetTitle("dR");
	gStyle->SetOptStat("omen");



	//electron histograms
	m_genRecoElectronPtRatio[nDirectory] = {m_histoFolder.make<TH1D>("genRecoElectronPtRatio" , "ratio of gen pt to reco pt for electrons" , 50 , 0 , 5 )};


	m_genElectronPt[nDirectory] = {m_histoFolder.make<TH1D>("genElectronPt","Pt for gen electrons",100,0,500)};
	m_genElectronPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_genElectronPhi[nDirectory] = {m_histoFolder.make<TH1D>("genElectronPhi","Phi for gen electrons",100,-3.4,3.4)};
	m_genElectronPhi[nDirectory]->GetXaxis()-> SetTitle("Phi");

	m_genElectronEta[nDirectory] = {m_histoFolder.make<TH1D>("genElectronEta","Eta for gen electrons",100,-3.4,3.4)};
	m_genElectronEta[nDirectory]->GetXaxis()-> SetTitle("Eta");

	m_recoElectronPt[nDirectory] = {m_histoFolder.make<TH1D>("recoElectronPt","Pt for reco electrons",100,0,500)};
	m_recoElectronPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_recoElectronPhi[nDirectory] = {m_histoFolder.make<TH1D>("recoElectronPhi","Phi for reco electrons",100,-3.4,3.4)};
	m_recoElectronPhi[nDirectory]->GetXaxis()-> SetTitle("Phi");

	m_recoElectronEta[nDirectory] = {m_histoFolder.make<TH1D>("recoElectronEta","Eta for reco electrons",100,-3.4,3.4)};
	m_recoElectronEta[nDirectory]->GetXaxis()-> SetTitle("Eta");


	m_WrongSignGenRecoElectronPtRatio[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignGenRecoElectronPtRatio" , "ratio of gen pt to reco pt (wrong charge matched electrons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_WrongSignRecoElectronPt[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignRecoElectronPt" , "Pt for wrong sign reco electrons" , 100 , 0 , 500 )};
	m_WrongSignRecoElectronPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_WrongSignRecoElectronPhi[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignRecoElectronPhi" , "phi for wrong sign electrons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoElectronPhi[nDirectory]->GetXaxis()->SetTitle("phi");

	m_WrongSignRecoElectronEta[nDirectory] = {m_histoFolder.make<TH1D>("WrongSignRecoElectronEta" , "Eta for wrong sign electrons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoElectronEta[nDirectory]->GetXaxis()->SetTitle("eta");

	m_RightSignGenRecoElectronPtRatio[nDirectory] = {m_histoFolder.make<TH1D>("RightSignGenRecoElectronPtRatio" , "ratio of gen pt to reco pt (right charge matched electrons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_RightSignRecoElectronPt[nDirectory] = {m_histoFolder.make<TH1D>("RightSignRecoElectronPt" , "Pt for right sign reco electrons" , 100 , 0 , 500 )};
	m_RightSignRecoElectronPt[nDirectory]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_RightSignRecoElectronPhi[nDirectory] = {m_histoFolder.make<TH1D>("RightSignRecoElectronPhi" , "phi for right sign electrons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoElectronPhi[nDirectory]->GetXaxis()->SetTitle("phi");

	m_RightSignRecoElectronEta[nDirectory] = {m_histoFolder.make<TH1D>("RightSignRecoElectronEta" , "Eta for right sign electrons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoElectronEta[nDirectory]->GetXaxis()->SetTitle("eta");
	
	m_bestElectrondR[nDirectory] = {m_histoFolder.make<TH1D>("bestElectrondR" , "smallest deltaR (electrons)" , 50 , 0 , 0.03 )};
	m_bestElectrondR[nDirectory]->GetXaxis()-> SetTitle("dR");
	gStyle->SetOptStat("omen");

	m_TWgenMuon[nDirectory] = {m_histoFolder.make<TH1D>("TWgenMuon" , "t->W status for muons" , 2 , 1 , 2 )};
	m_TWgenMuon[nDirectory]->GetXaxis()->SetBinLabel(1,"traced to tW");
	m_TWgenMuon[nDirectory]->GetXaxis()->SetBinLabel(2,"not traced to tW");
	m_TWgenElectron[nDirectory] = {m_histoFolder.make<TH1D>("TWgenElectron" , "t->W status for electrons" , 2 , 1 , 2 )};
	m_TWgenElectron[nDirectory]->GetXaxis()->SetBinLabel(1,"traced to tW");
	m_TWgenElectron[nDirectory]->GetXaxis()->SetBinLabel(2,"not traced to tW");

	

}

//General histogram filling
void eventHistos2::fill(eventBits2& event) {

	//m_countHisto->Fill("count", event.count);

int iSize;

if(event.twoMuons){

		m_numRecoMuons[0]->Fill(event.muonRecoCount);
		m_numGenMuons[0]->Fill(event.muonGenCount);
		m_numRecoElectrons[0]->Fill(event.electronRecoCount);
		m_numGenElectrons[0]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[0]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[0]->Fill(event.electronGenCount-event.electronRecoCount);

	if(!event.failedMatch){
		m_matchStatus[0]->Fill("matched",1);


		m_genMuonPt[0]->Fill(event.lepton1Pt);
		m_genMuonPt[0]->Fill(event.lepton2Pt);

		m_genMuonPhi[0]->Fill(event.lepton1Phi);
		m_genMuonPhi[0]->Fill(event.lepton2Phi);

		m_genMuonEta[0]->Fill(event.lepton1Eta);
		m_genMuonEta[0]->Fill(event.lepton2Eta);

		m_recoMuonPt[0]->Fill(event.Muon1Pt);
		m_recoMuonPt[0]->Fill(event.Muon2Pt);

		m_recoMuonPhi[0]->Fill(event.Muon1Phi);
		m_recoMuonPhi[0]->Fill(event.Muon2Phi);

		m_recoMuonEta[0]->Fill(event.Muon1Eta);
		m_recoMuonEta[0]->Fill(event.Muon2Eta);

		if(event.Muon1chargeMatch){
			m_RightSignGenRecoMuonPtRatio[0]->Fill(event.Muon1PtRatio); 
			m_RightSignRecoMuonPhi[0]->Fill(event.Muon1Phi);
			m_RightSignRecoMuonPt[0]->Fill(event.Muon1Pt);
			m_RightSignRecoMuonEta[0]->Fill(event.Muon1Eta);

			//std::cout<<"ptRatio: ";
			//std::cout<<event.Muon1PtRatio<<std::endl;
			//std::cout<<"event weight: ";
			//std::cout<<event.eventWeight<<std::endl;
			//std::cout<<"--------------------"<<std::endl;
		}

		if(event.Muon2chargeMatch){
			m_RightSignGenRecoMuonPtRatio[0]->Fill(event.Muon2PtRatio); 
			m_RightSignRecoMuonPhi[0]->Fill(event.Muon2Phi);
			m_RightSignRecoMuonPt[0]->Fill(event.Muon2Pt);
			m_RightSignRecoMuonEta[0]->Fill(event.Muon2Eta);
		}

		if(!event.Muon1chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[0]->Fill(event.Muon1PtRatio); 
			m_WrongSignRecoMuonPhi[0]->Fill(event.Muon1Phi);
			m_WrongSignRecoMuonPt[0]->Fill(event.Muon1Pt);
			m_WrongSignRecoMuonEta[0]->Fill(event.Muon1Eta);
		}

		if(!event.Muon2chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[0]->Fill(event.Muon2PtRatio); 
			m_WrongSignRecoMuonPhi[0]->Fill(event.Muon2Phi);
			m_WrongSignRecoMuonPt[0]->Fill(event.Muon2Pt);
			m_WrongSignRecoMuonEta[0]->Fill(event.Muon2Eta);
		}

		if(event.Muon1TW){
			m_TWgenMuon[0]->Fill("traced to tW",1);
		}
		if(event.Muon2TW){
			m_TWgenMuon[0]->Fill("traced to tW",1);
		}
		if(!event.Muon1TW){
			m_TWgenMuon[0]->Fill("not traced to tW",1);
		}
		if(!event.Muon2TW){
			m_TWgenMuon[0]->Fill("not traced to tW",1);
		}


		m_bestMuondR[0]->Fill(event.Muon1dR);
		m_bestMuondR[0]->Fill(event.Muon2dR);
		m_genRecoMuonPtRatio[0]->Fill(event.Muon1PtRatio);
		m_genRecoMuonPtRatio[0]->Fill(event.Muon2PtRatio);

	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[0]->Fill("unmatched",1);

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPhi[0]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedEta[0]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPt[0]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedDR[0]->Fill(event.unmatchedDR[iValue]);
		}
	}
} 
	
if(event.twoElectrons){

		m_numRecoMuons[1]->Fill(event.muonRecoCount);
		m_numGenMuons[1]->Fill(event.muonGenCount);
		m_numRecoElectrons[1]->Fill(event.electronRecoCount);
		m_numGenElectrons[1]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[1]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[1]->Fill(event.electronGenCount-event.electronRecoCount);

	if(!event.failedMatch){
		m_matchStatus[1]->Fill("matched",1);


		m_genElectronPt[1]->Fill(event.lepton1Pt);
		m_genElectronPt[1]->Fill(event.lepton2Pt);

		m_genElectronPhi[1]->Fill(event.lepton1Phi);
		m_genElectronPhi[1]->Fill(event.lepton2Phi);

		m_genElectronEta[1]->Fill(event.lepton1Eta);
		m_genElectronEta[1]->Fill(event.lepton2Eta);

		m_recoElectronPt[1]->Fill(event.Electron1Pt);
		m_recoElectronPt[1]->Fill(event.Electron2Pt);

		m_recoElectronPhi[1]->Fill(event.Electron1Phi);
		m_recoElectronPhi[1]->Fill(event.Electron2Phi);

		m_recoElectronEta[1]->Fill(event.Electron1Eta);
		m_recoElectronEta[1]->Fill(event.Electron2Eta);

		if(event.Electron1chargeMatch){
			m_RightSignGenRecoElectronPtRatio[1]->Fill(event.Electron1PtRatio); 
			m_RightSignRecoElectronPhi[1]->Fill(event.Electron1Phi);
			m_RightSignRecoElectronPt[1]->Fill(event.Electron1Pt);
			m_RightSignRecoElectronEta[1]->Fill(event.Electron1Eta);
		}

		if(event.Electron2chargeMatch){
			m_RightSignGenRecoElectronPtRatio[1]->Fill(event.Electron2PtRatio); 
			m_RightSignRecoElectronPhi[1]->Fill(event.Electron2Phi);
			m_RightSignRecoElectronPt[1]->Fill(event.Electron2Pt);
			m_RightSignRecoElectronEta[1]->Fill(event.Electron2Eta);
		}

		if(!event.Electron1chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[1]->Fill(event.Electron1PtRatio); 
			m_WrongSignRecoElectronPhi[1]->Fill(event.Electron1Phi);
			m_WrongSignRecoElectronPt[1]->Fill(event.Electron1Pt);
			m_WrongSignRecoElectronEta[1]->Fill(event.Electron1Eta);
		}

		if(!event.Electron2chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[1]->Fill(event.Electron2PtRatio); 
			m_WrongSignRecoElectronPhi[1]->Fill(event.Electron2Phi);
			m_WrongSignRecoElectronPt[1]->Fill(event.Electron2Pt);
			m_WrongSignRecoElectronEta[1]->Fill(event.Electron2Eta);
		}
		if(event.Electron1TW){
			m_TWgenElectron[1]->Fill("traced to tW",1);
		}
		if(event.Electron2TW){
			m_TWgenElectron[1]->Fill("traced to tW",1);
		}
		if(!event.Electron1TW){
			m_TWgenElectron[1]->Fill("not traced to tW",1);
		}
		if(!event.Electron2TW){
			m_TWgenElectron[1]->Fill("not traced to tW",1);
		}
		m_bestElectrondR[1]->Fill(event.Electron1dR);
		m_bestElectrondR[1]->Fill(event.Electron2dR);

		m_genRecoElectronPtRatio[1]->Fill(event.Electron1PtRatio);
		m_genRecoElectronPtRatio[1]->Fill(event.Electron2PtRatio);
	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[1]->Fill("unmatched",1);
		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPhi[1]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedEta[1]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPt[1]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedDR[1]->Fill(event.unmatchedDR[iValue]);
		}


	} 
}


	
if(event.muonElectron){
		m_numRecoMuons[2]->Fill(event.muonRecoCount);
		m_numGenMuons[2]->Fill(event.muonGenCount);
		m_numRecoElectrons[2]->Fill(event.electronRecoCount);
		m_numGenElectrons[2]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[2]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[2]->Fill(event.electronGenCount-event.electronRecoCount);
	if(!event.failedMatch){
		m_matchStatus[2]->Fill("matched",1);



		if(abs(event.lepton1Id)==13){
			m_genMuonPt[2]->Fill(event.lepton1Pt);
			m_genMuonPhi[2]->Fill(event.lepton1Phi);
			m_genMuonEta[2]->Fill(event.lepton1Eta);

			m_genElectronPt[2]->Fill(event.lepton2Pt);
			m_genElectronPhi[2]->Fill(event.lepton2Phi);
			m_genElectronEta[2]->Fill(event.lepton2Eta);
		}
		else{ 		
			m_genMuonPt[2]->Fill(event.lepton2Pt);
			m_genMuonPhi[2]->Fill(event.lepton2Phi);
			m_genMuonEta[2]->Fill(event.lepton2Eta);

			m_genElectronPt[2]->Fill(event.lepton1Pt);
			m_genElectronPhi[2]->Fill(event.lepton1Phi);
			m_genElectronEta[2]->Fill(event.lepton1Eta);
		}

		m_recoMuonPt[2]->Fill(event.Muon1Pt);
		m_recoMuonPhi[2]->Fill(event.Muon1Phi);
		m_recoMuonEta[2]->Fill(event.Muon1Eta);

		m_recoElectronPt[2]->Fill(event.Electron1Pt);
		m_recoElectronPhi[2]->Fill(event.Electron1Phi);
		m_recoElectronEta[2]->Fill(event.Electron1Eta);

		if(event.Muon1chargeMatch){
			m_RightSignGenRecoMuonPtRatio[2]->Fill(event.Muon1PtRatio); 
			m_RightSignRecoMuonPhi[2]->Fill(event.Muon1Phi);
			m_RightSignRecoMuonPt[2]->Fill(event.Muon1Pt);
			m_RightSignRecoMuonEta[2]->Fill(event.Muon1Eta);
		}

		if(!event.Muon1chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[2]->Fill(event.Muon1PtRatio); 
			m_WrongSignRecoMuonPhi[2]->Fill(event.Muon1Phi);
			m_WrongSignRecoMuonPt[2]->Fill(event.Muon1Pt);
			m_WrongSignRecoMuonEta[2]->Fill(event.Muon1Eta);
		}

		if(event.Electron1chargeMatch){
			m_RightSignGenRecoElectronPtRatio[2]->Fill(event.Electron1PtRatio); 
			m_RightSignRecoElectronPhi[2]->Fill(event.Electron1Phi);
			m_RightSignRecoElectronPt[2]->Fill(event.Electron1Pt);
			m_RightSignRecoElectronEta[2]->Fill(event.Electron1Eta);
		}

		if(!event.Electron1chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[2]->Fill(event.Electron1PtRatio); 
			m_WrongSignRecoElectronPhi[2]->Fill(event.Electron1Phi);
			m_WrongSignRecoElectronPt[2]->Fill(event.Electron1Pt);
			m_WrongSignRecoElectronEta[2]->Fill(event.Electron1Eta);
		}

		if(event.Muon1TW){
			m_TWgenMuon[2]->Fill("traced to tW",1);
		}
		if(event.Electron1TW){
			m_TWgenElectron[2]->Fill("traced to tW",1);
		}
		if(!event.Muon1TW){
			m_TWgenMuon[2]->Fill("not traced to tW",1);
		}
		if(!event.Electron1TW){
			m_TWgenElectron[2]->Fill("not traced to tW",1);
		}

		m_bestMuondR[2]->Fill(event.Muon1dR);
		m_bestElectrondR[2]->Fill(event.Electron1dR);
		m_genRecoMuonPtRatio[2]->Fill(event.Muon1PtRatio);
		m_genRecoElectronPtRatio[2]->Fill(event.Electron1PtRatio);
	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[2]->Fill("unmatched",1);

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPhi[2]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedEta[2]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPt[2]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedDR[2]->Fill(event.unmatchedDR[iValue]);
		}

	}
}

if(event.muonTau){
		m_numRecoMuons[3]->Fill(event.muonRecoCount);
		m_numGenMuons[3]->Fill(event.muonGenCount);
		m_numRecoElectrons[3]->Fill(event.electronRecoCount);
		m_numGenElectrons[3]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[3]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[3]->Fill(event.electronGenCount-event.electronRecoCount);
	if(!event.failedMatch){
		m_matchStatus[3]->Fill("matched",1);


		if(abs(event.lepton1Id)==13){
			m_genMuonPt[3]->Fill(event.lepton1Pt);
			m_genMuonPhi[3]->Fill(event.lepton1Phi);
			m_genMuonEta[3]->Fill(event.lepton1Eta);
		}
		else{ 		
			m_genMuonPt[3]->Fill(event.lepton2Pt);
			m_genMuonPhi[3]->Fill(event.lepton2Phi);
			m_genMuonEta[3]->Fill(event.lepton2Eta);
		}

		m_recoMuonPt[3]->Fill(event.Muon1Pt);
		m_recoMuonPhi[3]->Fill(event.Muon1Phi);
		m_recoMuonEta[3]->Fill(event.Muon1Eta);

		if(event.Muon1chargeMatch){
			m_RightSignGenRecoMuonPtRatio[3]->Fill(event.Muon1PtRatio); 
			m_RightSignRecoMuonPhi[3]->Fill(event.Muon1Phi);
			m_RightSignRecoMuonPt[3]->Fill(event.Muon1Pt);
			m_RightSignRecoMuonEta[3]->Fill(event.Muon1Eta);
		}


		if(!event.Muon1chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[3]->Fill(event.Muon1PtRatio); 
			m_WrongSignRecoMuonPhi[3]->Fill(event.Muon1Phi);
			m_WrongSignRecoMuonPt[3]->Fill(event.Muon1Pt);
			m_WrongSignRecoMuonEta[3]->Fill(event.Muon1Eta);
		}
		if(event.Muon1TW){
			m_TWgenMuon[3]->Fill("traced to tW",1);
		}
		if(!event.Muon1TW){
			m_TWgenMuon[3]->Fill("not traced to tW",1);
		}

		m_bestMuondR[3]->Fill(event.Muon1dR);
		m_genRecoMuonPtRatio[3]->Fill(event.Muon1PtRatio);

	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[3]->Fill("unmatched",1);

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPhi[3]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedEta[3]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPt[3]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedDR[3]->Fill(event.unmatchedDR[iValue]);
		}

	}
}

if(event.electronTau){
		m_numRecoMuons[4]->Fill(event.muonRecoCount);
		m_numGenMuons[4]->Fill(event.muonGenCount);
		m_numRecoElectrons[4]->Fill(event.electronRecoCount);
		m_numGenElectrons[4]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[4]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[4]->Fill(event.electronGenCount-event.electronRecoCount);
	if(!event.failedMatch){
		m_matchStatus[4]->Fill("matched",1);


		if(abs(event.lepton1Id)==11){
			m_genElectronPt[4]->Fill(event.lepton1Pt);
			m_genElectronPhi[4]->Fill(event.lepton1Phi);
			m_genElectronEta[4]->Fill(event.lepton1Eta);
		}
		else{ 		
			m_genElectronPt[4]->Fill(event.lepton2Pt);
			m_genElectronPhi[4]->Fill(event.lepton2Phi);
			m_genElectronEta[4]->Fill(event.lepton2Eta);
		}

		m_recoElectronPt[4]->Fill(event.Electron1Pt);
		m_recoElectronPhi[4]->Fill(event.Electron1Phi);
		m_recoElectronEta[4]->Fill(event.Electron1Eta);

		if(event.Electron1chargeMatch){
			m_RightSignGenRecoElectronPtRatio[4]->Fill(event.Electron1PtRatio); 
			m_RightSignRecoElectronPhi[4]->Fill(event.Electron1Phi);
			m_RightSignRecoElectronPt[4]->Fill(event.Electron1Pt);
			m_RightSignRecoElectronEta[4]->Fill(event.Electron1Eta);
		}

		if(!event.Electron1chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[4]->Fill(event.Electron1PtRatio); 
			m_WrongSignRecoElectronPhi[4]->Fill(event.Electron1Phi);
			m_WrongSignRecoElectronPt[4]->Fill(event.Electron1Pt);
			m_WrongSignRecoElectronEta[4]->Fill(event.Electron1Eta);
		}
				

		if(event.Electron1TW){
			m_TWgenElectron[4]->Fill("traced to tW",1);
		}
		if(!event.Electron1TW){
			m_TWgenElectron[4]->Fill("not traced to tW",1);
		}

		m_bestElectrondR[4]->Fill(event.Electron1dR);
		m_genRecoElectronPtRatio[4]->Fill(event.Electron1PtRatio);
	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[4]->Fill("unmatched",1);

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPhi[4]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedEta[4]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedPt[4]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedDR[4]->Fill(event.unmatchedDR[iValue]);
		}	
	}
}

//top->cd();

}


	
