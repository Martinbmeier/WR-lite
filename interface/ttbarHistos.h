#ifndef ttbarHistos_h
#define ttbarHistos_h


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

//ROOT CLASSES
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
//C++ CLASSES
#include <iostream>
//LOCAL CLASSES
#include "ttbarHistos.h"
#include "eventBits2.h"


class ttbarHistos {
	public:
	ttbarHistos();
	void book(TFileDirectory histoFolder);
	void fill(double mupT, double epT, double weight);

	private:
	
	TFileDirectory m_histoFolder;

	TH1D* m_electronpT;
	TH2D* m_electronMuonpT;

	//std::vector<TH1D*> m_recoMuonPt{3};
   
};

#endif