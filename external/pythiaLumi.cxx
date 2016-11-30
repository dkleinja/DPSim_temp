#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "Randomize.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "Pythia8/Pythia.h"

#include "DPDetectorConstruction.h"
#include "DPPrimaryGeneratorAction.h"
#include "DPRunAction.h"
#include "DPEventAction.h"
#include "DPTrackingAction.h"
#include "DPSteppingAction.h"
#include "DPVertexGenerator.h"
#include "DPSimConfig.h"

using namespace std;
using namespace Pythia8;

int main(int argc, char* argv[])
{
    //Initialize the configuration
    DPSimConfig* p_config = DPSimConfig::instance();
    p_config->init(argv[1]);

    //Initialize vertex generator
    DPVertexGenerator* p_vertexGen = DPVertexGenerator::instance();
    p_vertexGen->init();

    //Initialize Output file
    int eventID;
    int n;
    int pdg[10000];
    int parentID[10000];
    float energy[10000];
    TClonesArray* p_pos = new TClonesArray("TVector3");  p_pos->BypassStreamer(); TClonesArray& pos = *p_pos;
    TClonesArray* p_mom = new TClonesArray("TVector3");  p_mom->BypassStreamer(); TClonesArray& mom = *p_mom;

    TH1I *Hcount = new TH1I("Hcount", "p+p collisions with/without particles", 2, 0, 2);
    TH2F *Hbeam = new TH2F("Hbeam","Hbeam", 150, -1.5, 1.5, 150, -1.5, 1.5);
    
    TFile* saveFile = new TFile(p_config->outputFileName, "recreate");
    TTree* saveTree = new TTree("save", "save");

    saveTree->Branch("eventID", &eventID, "eventID/I");
    saveTree->Branch("n", &n, "n/I");
    saveTree->Branch("pdg", pdg, "pdg[n]/I");
    saveTree->Branch("parentID", parentID, "parentID[n]/I");
    saveTree->Branch("energy", energy, "energy[n]/F");
    saveTree->Branch("pos", &p_pos, 256000, 99);
    saveTree->Branch("mom", &p_mom, 256000, 99);

    //Initialize pythia for pp and pn collisions
    Pythia ppGen;
    ppGen.readFile(p_config->pythiaConfig.Data());
    ppGen.readString(Form("Random:seed = %d", p_config->seed));
    ppGen.readString("Beams:idB = 2212");
    ppGen.init();

    Pythia pnGen;
    pnGen.readFile(p_config->pythiaConfig.Data());
    pnGen.readString(Form("Random:seed = %d", p_config->seed));
    pnGen.readString("Beams:idB = 2112");
    pnGen.init();

    //simple lorentzvector
    TLorentzVector thetacheck;
    //Main loop
    int nEvents = p_config->nEvents;
    for(int i = 0; i < nEvents; ++i)
    {
        double zvtx = p_vertexGen->generateVertex();
	double xvtx = G4RandGauss::shoot(0.,0.68);
	double yvtx = G4RandGauss::shoot(0.,0.76);
	if(xvtx*xvtx/1.7/1.7*4+ yvtx*yvtx/1.9/1.9*4 > 1) continue;
	Hbeam -> Fill(xvtx, yvtx);
        double pARatio = p_vertexGen->getPARatio();

        Pythia* p_pythia = G4UniformRand() < pARatio ? &ppGen : &pnGen;
        Event& events = p_pythia->event;
        while(!p_pythia->next()) {}

        n = 0;
        eventID = i;
        int nParticles = events.size();
        for(int j = 1; j < nParticles; ++j)
        {
          thetacheck.SetPxPyPzE(events[j].px(), events[j].py(), events[j].pz(), events[j].e());
          if(events[j].status()<0)continue;
	  //if(!events[j].isCharged()) continue;
	  if(events[j].id() > 2000) continue;
	  if(abs(events[j].py()/events[j].px()) > 0.176 && abs(events[j].py()/events[j].px()) < 5.67 )continue;//#phi = #pm 10 deg at 0, 90, 180, 360
	  if(abs(thetacheck.CosTheta()) < 0.174 && (events[j].zProd()/10. + zvtx)< -0.)//#theta = #pm 10 deg
	  {

                //Fill normal particles 
                pdg[n] = events[j].id();
		parentID[n] = events[events[j].mother1()].id();
		energy[n] = events[j].e();
		new(mom[n]) TVector3(events[j].px(), events[j].py(), events[j].pz());
                new(pos[n]) TVector3(events[j].xProd()/10. + xvtx, events[j].yProd()/10. + yvtx, events[j].zProd()/10. + zvtx);

                ++n;

                /*
                //If mother particle exists, fill mother as well
                int motherID = events[j].mother1();
                if(motherID < 1 || motherID >= nParticles) continue;
                if(abs(events[motherID].id()) < 10) continue;
                pdg[n] = events[motherID].id();
                new(mom[n]) TVector3(events[motherID].px(), events[motherID].py(), events[motherID].pz());
                new(pos[n]) TVector3(events[motherID].xProd()/10., events[motherID].yProd()/10., events[motherID].zProd()/10. + zvtx);
                ++n;*/
            }
        }

        if(n == 0){
	  Hcount -> Fill(0);
	  //continue;
	}
	else{
	  Hcount -> Fill(1);
	  saveTree->Fill();
	  mom.Clear();
	  pos.Clear();
	}
    }
    
    //finalize
    saveFile->cd();
    saveTree->Write();
    Hcount->Write();
    Hbeam->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}
