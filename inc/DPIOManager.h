#ifndef DPIOManager_H
#define DPIOManager_H

#include <iostream>
#include <vector>
#include <list>
#include <map>

#include "G4Event.hh"
#include "G4ThreeVector.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "DPSimConfig.h"
#include "DPDigitizer.h"
#include "DPTriggerAnalyzer.h"
#include "DPDummyRecon.h"
#include "DPMCRawEvent.h"
#include "DPVirtualHit.h"

class DPIOManager
{
public:
    static DPIOManager* instance();

    DPIOManager();
    ~DPIOManager();

    //!Enum type for buffer state
    enum BufferState {INTERNAL, CLEAN, FILLED, FLUSH};

    //!Enum type for save mode
    enum SaveMode {EVERYTHING, PRIMARY, HITSONLY, INACCONLY};

    //!Initialize, called at the beginning of each Run
    void initialize(int runID);

    //!reset the interal containers at the end of each event
    void reset();

    //!set the generation info
    void fillOneDimuon(double weight, const DPMCDimuon& dimuon);

    //!Fill one event, called at the end of each event
    void fillOneEvent(const G4Event*);

    //!apply all the analysis module enabled, then push the event to tree
    void finalizeEvent();

    //!Fill one track, called at the beginning of each track
    void fillOneTrack(const DPMCTrack& mcTrack, bool keep = false);
    void updateOneTrack(unsigned int trackID, G4ThreeVector pos, G4ThreeVector mom, bool keep = false);

    //!Fill the hit list from G4Event, and digitize the virtual hits called by fillOneEvent
    void fillHitsList(const G4Event* theEvent);

    //!Re-index the tracks and hits, called by fillOneEvent
    void reIndex();

    //!External control of the single buffer state
    void setBufferState(BufferState flag) { bufferState = flag; }

    //!Finalize, called at the end of each Run
    void finalize();

private:
    //!static pointer
    static DPIOManager* p_IOmamnger;

    //!Pointer to the configuration
    DPSimConfig* p_config;

    //!Pointer to the digitizer
    DPDigitizer* p_digitizer;

    //!Pointer to the trigger analyzer
    DPTriggerAnalyzer* p_triggerAna;

    //!Pointer to dummy reconstruction
    DPDummyRecon* p_dummyRecon;

    //!save mode
    SaveMode saveMode;

    //!@name Output file
    //@{
    TFile* saveFile;
    TTree* saveTree;
    TTree* configTree;
    DPMCRawEvent* rawEvent;
    //@}

    //!temporary buffer of the single event
    DPMCRawEvent* singleEvent;

    //!flag used to indicate the state of single event buffer
    BufferState bufferState;

    //!container of tracks
    //@{
    std::map<unsigned int, unsigned int> trackIDs; //maps real trackID to index in tracks vector
    std::list<DPVirtualHit> hits;
    std::vector<std::pair<DPMCTrack, bool> > tracks;
    //@}

    //!ID of the hit collection
    int sensHitColID;
};

#endif
