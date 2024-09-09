#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

//#include "globals.hh"
#include "g4root.hh"
//#include <TFile.h>
//#include <TTree.h>


class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
MyPrimaryGenerator();
~MyPrimaryGenerator();

virtual void GeneratePrimaries(G4Event *);

private:

G4ParticleGun *fParticleGun;

//    // ROOT file and TTree for saving the data
//    TFile* rootFile;   // ROOT file pointer
//    TTree* tree;       // ROOT TTree pointer
//
//    // Variables to store in the ROOT tree
//    G4double energy;   // Photon energy (keV)
//    G4double theta;    // Theta angle (radians)
//    G4double phi;      // Phi angle (radians)
//    G4double cosTheta; // Cosine of theta
 

};

#endif
