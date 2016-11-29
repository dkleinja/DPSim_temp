#ifndef DPSimConfig_H
#define DPSimConfig_H

#include <map>
#include <string>

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>

class DPSimConfig: public TObject
{
public:
    static DPSimConfig* instance();
    DPSimConfig();

    //!parse the config file and set the config
    void init(TString configFile);

    //!check if there are any inconsistencies
    bool sanityCheck();

    //!check file status
    bool checkFile(TString filename);

public:
    //!Random seed
    Int_t seed;

    //!Number of events to run
    Int_t nEvents;
    Int_t printFreq;

    //!MC version hash string
    TString version;

    //!beam setup
    //@{
    Int_t    bucket_size;     // How many protons per event in gun generator
    Double_t beamMomentum;
    Double_t beamCurrent;
    //@}

    //!Detector setup
    //@{
    TString geometryGDMLInput;// The GDML file for geometry
    TString detectorEffResol; // The table of channel-by-channel efficiency and resolution
    TString geometrySchema;   // The sql schema that GMC pulls the geometry information from
    TString mysqlServer;      // Address of the SQL Server, shouldn't need to modify
    Int_t   mysqlPort;        // The port number for MySQL
    TString login;            // The login for the SQL server
    TString password;         // The password for the SQL server
    //@}

    //!Trigger matrix input
    TString triggerMatrix;    // Input table file of the trigger matrix

    //!Magnetic field setup
    //@{
    TString fMagMap;          // Name of the ascii text file that contains the fmag map
    TString kMagMap;          // Name of the ascii text file that contains the kmag map
    Double_t fMagMultiplier;  // Multiplies the strength of FMAG's field
    Double_t kMagMultiplier;  // Multiplies the strength of KMAG's field
    //@}

    //!Event generation setup
    //@{
    TString generatorType;    // The type of event generator running, i.e. single/dimuon/external
    TString generatorEng;     // The generator engine, legacyDY/legacyJPsi/legacyPsip/pythia/geant
    TString externalInput;    // Input ROOT file containing the generator info
    TString pythiaConfig;     // pythia configuration file
    TString customLUT;        // custom dimuon cross section look-up table, used in custom mode
    TString testParticle;     // particle type for single particle test generator
    TString physicsList;      // physics list used to simulate detector effects
    //@}

    //!Beam line objects setup
    //@{
    bool targetInBeam;
    bool dumpInBeam;
    bool instruInBeam;
    bool airInBeam;
    //@}

    //!Beam profile parameters
    //@{
    bool beamProfile;
    double beamCenterX;
    double beamCenterY;
    //@}

    //!Event vertex biasing, i.e. generate vertex evenly instead of according to material
    bool biasVertexGen;


    //!Force pion and kaon decay by modifying particle lifetime
    bool forcePionDecay;
    bool forceKaonDecay;

    //!I/O setup
    //@{
    TString configFileName;   // name of configuration file
    TString outputFileName;   // The database name that output goes to
    TString outputMode;       // save what kind of events to disk
    //@}

    //! Optional analysis component enable flag
    bool enableDummyRecon;    // enable the dummy reconstruction

    //! phase space constrain
    //@{
    Double_t x1Min;
    Double_t x1Max;
    Double_t x2Min;
    Double_t x2Max;
    Double_t xfMin;
    Double_t xfMax;
    Double_t massMin;
    Double_t massMax;
    Double_t cosThetaMin;
    Double_t cosThetaMax;
    Double_t zOffsetMin;
    Double_t zOffsetMax;
    //@}

    //! Run-accumulated variables
    //@{
    Int_t nEventsThrown;
    Int_t nEventsPhysics;
    Int_t nEventsAccepted;
    //@}

    //! boolean flags used only internally
    //@{
    bool dimuonMode;
    bool drellyanMode;
    //@}

private: //used only for parsing
    //! general parser
    //@{
    void parseConfig(TString configFile);
    TString expandEnv(const TString& input) const;
    //@}

    //! get string/bool/float/int config entry
    //@{
    TString  pString(TString name, TString default_val = "N/A");
    bool     pBool(TString name, bool default_val);
    Double_t pDouble(TString name, Double_t default_val);
    Int_t    pInt(TString name, Int_t default_val);
    //@}

    //! config group symbol map
    std::map<TString, TString> symbols;

    //! static pointer
    static DPSimConfig* p_config;

    ClassDef(DPSimConfig, 1)
};

#endif
