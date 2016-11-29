#include "DPSensitiveDetector.h"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>

DPSensitiveDetector::DPSensitiveDetector(const G4String& name, const G4String& hitCollectionName): G4VSensitiveDetector(name), hcID(-1), theHC(NULL)
{
    collectionName.insert(hitCollectionName);
}

DPSensitiveDetector::~DPSensitiveDetector() {}

void DPSensitiveDetector::Initialize(G4HCofThisEvent* HCE)  //may have memory leaks here, refer to BZ's code
{
    theHC = new DPVHitCollection(SensitiveDetectorName, collectionName[0]);

    if(hcID < 0) hcID = GetCollectionID(0);
    HCE->AddHitsCollection(hcID, theHC);
}

G4bool DPSensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* hist)
{
    //reject neutral tracks for now
    if(aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0) return false;

    G4Track* theTrack = aStep->GetTrack();
    if(theTrack->GetMomentum()[2] < 1.*MeV) return false;

    DPVirtualHit* newHit = new DPVirtualHit();
    newHit->particleID = theTrack->GetTrackID();
    newHit->particlePDG = theTrack->GetDefinition()->GetPDGEncoding();
    newHit->edep = aStep->GetTotalEnergyDeposit();
    newHit->mom = theTrack->GetMomentum();
    newHit->pos = theTrack->GetPosition();
    newHit->detectorGroupName = theTrack->GetVolume()->GetName();
    newHit->digiHits.clear();
    theHC->insert(newHit);

#ifdef DEBUG_TR
    std::cout << __FILE__ << " " << __FUNCTION__ << " processing hits for track " << theTrack->GetTrackID()
              << " with PDG = " << theTrack->GetDefinition()->GetPDGEncoding()
              << " on " << newHit->detectorGroupName << std::endl;
#endif
    return true;
}

void DPSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{
}
