#include "generator.hh"
#include "G4RandomTools.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{

  fParticleGun = new G4ParticleGun(1);
  G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();

  //G4ParticleDefinition * particle = particleTable->FindParticle("geantino");
  //G4ParticleDefinition * particle = particleTable->FindParticle("mu-");
  G4ParticleDefinition * particle = particleTable->FindParticle("gamma");


  //G4ThreeVector mom(0.,0.,-1.); 
  //fParticleGun->SetParticleMomentumDirection(mom);
  
  
  //fParticleGun->SetParticleMomentum(519.9*GeV);
  fParticleGun->SetParticleMomentum(519.9*keV);

  fParticleGun->SetParticleDefinition(particle);



   // // Define branches for Energy, Phi, Theta, and cos(Theta)
   // tree->Branch("Energy", &energy, "Energy/D");
   // tree->Branch("Phi", &phi, "Phi/D");
   // tree->Branch("Theta", &theta, "Theta/D");
   // tree->Branch("CosTheta", &cosTheta, "CosTheta/D");


}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;


   // // Write and close the ROOT file
   // rootFile->Write();
   // rootFile->Close();
   // delete fParticleGun;

}

void MyPrimaryGenerator::GeneratePrimaries(G4Event * anEvent)
{
  G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();

  if (particle == G4Geantino::Geantino())
  {
     G4int Z =27;
     G4int A =60;
     G4double charge = 0.*eplus;
     G4double energy = 0.*keV;

     G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z,A,energy);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(charge); 
  } 


      // Generate a random position within the 1m x 1m area
    G4double x = (G4UniformRand() - 0.5) * 7.0 * cm; // 0cm to 7cm
    G4double y = (G4UniformRand() - 0.5) * 7.0 * cm; // 0cm to 7cm
    G4double z = (G4UniformRand() - 0) * 7.0 * cm; // 0cm to 7cm

    G4ThreeVector pos(x, y, z);
    //G4ThreeVector pos(0, 0, 0.1);
    fParticleGun->SetParticlePosition(pos);


//    // Generate uniform direction for theta and phi
//    theta = std::acos(2 * G4UniformRand() - 1);  // theta in [0, pi]
//    phi = 2 * CLHEP::pi * G4UniformRand();       // phi in [0, 2*pi]
//    cosTheta = std::cos(theta);



    // Generate uniform direction for theta and phi
    G4double theta = std::acos(2 * G4UniformRand() - 1);  // theta in [0, pi]
    G4double phi = 2 * CLHEP::pi * G4UniformRand();       // phi in [0, 2*pi]
    G4double cosTheta = std::cos(theta);



    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(6, 0, 519.9);    
    analysisManager->FillNtupleDColumn(6, 1, phi);    
    analysisManager->FillNtupleDColumn(6, 2, theta);    
    analysisManager->FillNtupleDColumn(6, 3, cosTheta);    
    analysisManager->FillNtupleDColumn(6, 4, x);    
    analysisManager->FillNtupleDColumn(6, 5, y);    
    analysisManager->FillNtupleDColumn(6, 6, z);    
    analysisManager->AddNtupleRow(6);





    // Convert spherical to Cartesian coordinates
    G4double px = std::sin(theta) * std::cos(phi);
    G4double py = std::sin(theta) * std::sin(phi);
    G4double pz = std::cos(theta);

    // Set the particle momentum direction
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    

  fParticleGun->GeneratePrimaryVertex(anEvent);  

}