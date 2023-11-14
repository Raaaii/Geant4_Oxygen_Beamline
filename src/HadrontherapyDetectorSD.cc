//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyDetectorHit.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyLet.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "HadrontherapyMatrix.hh"


#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "HadrontherapySteppingAction.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UserEventAction.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "HadrontherapyRunAction.hh"
#include "G4SystemOfUnits.hh"
#include "HadrontherapyRBE.hh"
#include <G4AccumulableManager.hh>

#include <iostream> // SERENA
#include <fstream>// SERENA

// LOCK DECLARATION***********************************************
namespace {
 //Mutex to acquire access to singleton instance check/creation
 G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
 //Mutex to acquire accss to histograms creation/access
 //It is also used to control all operations related to histos
 //File writing and check analysis
 G4Mutex dataManipulationMutex = G4MUTEX_INITIALIZER;
}
// ****************************************************************

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorSD::HadrontherapyDetectorSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="HadrontherapyDetectorHitsCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
    
   
}
 
/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorSD::~HadrontherapyDetectorSD()
{}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new HadrontherapyDetectorHitsCollection(sensitiveDetectorName,
                                                             collectionName[0]);
}

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
    if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "RODetectorZDivisionPhys") return false;
    
    // Get kinetic energy
    G4Track * theTrack = aStep  ->  GetTrack();
    G4double kineticEnergy = theTrack->GetKineticEnergy();
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    //Get particle name
    G4String particleName =  particleDef -> GetParticleName();
    
    // Get particle PDG code
    G4int pdg = particleDef ->GetPDGEncoding();
    
    // Get unique track_id (in an event)
    G4int trackID = theTrack -> GetTrackID();
    
    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
    
    G4double DX = aStep -> GetStepLength();
    G4int Z = particleDef-> GetAtomicNumber();
    G4int A = particleDef-> GetAtomicMass();
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    
    // Read voxel indexes: i is the x index, k is the z index
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int k  = touchable->GetReplicaNumber(0);
    G4int i  = touchable->GetReplicaNumber(2);
    G4int j  = touchable->GetReplicaNumber(1);
    
    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    G4String namePre = volumePre->GetName();
    
    HadrontherapyMatrix* matrix = HadrontherapyMatrix::GetInstance();
    HadrontherapyLet* let = HadrontherapyLet::GetInstance();
    
    G4int* hitTrack = matrix -> GetHitTrack(i,j,k);

    // SERENA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //  ************************ Spectra *******************************
  
    if (pdg == 1000020040){
        G4AutoLock l(&dataManipulationMutex); // LOCK START
        G4double Xpos_PreStep =  aStep->GetPreStepPoint()->GetPosition().x(); // SERENA
        G4double Ek_PreStep = aStep -> GetPreStepPoint() -> GetKineticEnergy(); // SERENA
        
        if (i==0){
            std::ofstream WriteDataIn("spettro_He4_at-phantom-entrance.out", std::ios::app);
                WriteDataIn
                <<   Xpos_PreStep   <<'\t' //  5
                <<   Ek_PreStep     <<'\t' //  6
                <<   G4endl;
        }
        l.unlock(); // LOCK STOP
    
    }
    
    /*
    if (pdg != 11 && pdg != -11 && pdg != 22 && pdg !=2112){
    G4AutoLock l(&dataManipulationMutex); // LOCK START
        
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int parentID = theTrack -> GetParentID();
    G4ThreeVector init_pos = theTrack->GetVertexPosition();
    G4String process_name = "initial";
    if( theTrack->GetCreatorProcess()) process_name = theTrack->GetCreatorProcess()->GetProcessName();
    
      
    if (i==85){
              std::ofstream WriteDataIn("spettro_simulato_8_5mm.out", std::ios::app);
                 WriteDataIn
                  <<   i              << " " //  1
                  <<   eventID        << " " //  2
                  <<   trackID        << " " //  3
                  <<   parentID       << " " //  4
                  <<   pdg            << " " //  5
                  <<   Z              << " " //  6
                  <<   A              << " " //  7
                  <<   kineticEnergy  << " " //  8
                  <<   init_pos[0]    << " " //  9
                  <<   init_pos[1]    << " " //  10
                  <<   init_pos[2]    << " " //  11
                  <<   process_name          //  12
                  <<   G4endl;
      }
           
        
  
    if (i==1){
        std::ofstream WriteDataIn("spettro_simulato_0_1mm.out", std::ios::app);
           WriteDataIn
            <<   i              << " " //  1
            <<   eventID        << " " //  2
            <<   trackID        << " " //  3
            <<   parentID       << " " //  4
            <<   pdg            << " " //  5
            <<   Z              << " " //  6
            <<   A              << " " //  7
            <<   kineticEnergy  << " " //  8
            <<   init_pos[0]    << " " //  9
            <<   init_pos[1]    << " " //  10
            <<   init_pos[2]    << " " //  11
            <<   process_name          //  12
            <<   G4endl;
    
        }
    
        
        
    if (i==20){
        std::ofstream WriteDataIn("spettro_simulato_2mm.out", std::ios::app);
           WriteDataIn
            <<   i              << " " //  1
            <<   eventID        << " " //  2
            <<   trackID        << " " //  3
            <<   parentID       << " " //  4
            <<   pdg            << " " //  5
            <<   Z              << " " //  6
            <<   A              << " " //  7
            <<   kineticEnergy  << " " //  8
            <<   init_pos[0]    << " " //  9
            <<   init_pos[1]    << " " //  10
            <<   init_pos[2]    << " " //  11
            <<   process_name          //  12
            <<   G4endl;
        }
    
    if (i==30){
        std::ofstream WriteDataIn("spettro_simulato_3mm.out", std::ios::app);
           WriteDataIn
            <<   i              << " " //  1
            <<   eventID        << " " //  2
            <<   trackID        << " " //  3
            <<   parentID       << " " //  4
            <<   pdg            << " " //  5
            <<   Z              << " " //  6
            <<   A              << " " //  7
            <<   kineticEnergy  << " " //  8
            <<   init_pos[0]    << " " //  9
            <<   init_pos[1]    << " " //  10
            <<   init_pos[2]    << " " //  11
            <<   process_name          //  12
            <<   G4endl;
        }
    
    
    if (i==40){
        std::ofstream WriteDataIn("spettro_simulato_4mm.out", std::ios::app);
           WriteDataIn
            <<   i              << " " //  1
            <<   eventID        << " " //  2
            <<   trackID        << " " //  3
            <<   parentID       << " " //  4
            <<   pdg            << " " //  5
            <<   Z              << " " //  6
            <<   A              << " " //  7
            <<   kineticEnergy  << " " //  8
            <<   init_pos[0]    << " " //  9
            <<   init_pos[1]    << " " //  10
            <<   init_pos[2]    << " " //  11
            <<   process_name          //  12
            <<   G4endl;
}
    
    if (i==50){
                std::ofstream WriteDataIn("spettro_simulato_5mm.out", std::ios::app);
                   WriteDataIn
                    <<   i              << " " //  1
                    <<   eventID        << " " //  2
                    <<   trackID        << " " //  3
                    <<   parentID       << " " //  4
                    <<   pdg            << " " //  5
                    <<   Z              << " " //  6
                    <<   A              << " " //  7
                    <<   kineticEnergy  << " " //  8
                    <<   init_pos[0]    << " " //  9
                    <<   init_pos[1]    << " " //  10
                    <<   init_pos[2]    << " " //  11
                    <<   process_name          //  12
                    <<   G4endl;
        }
    
    
    if (i==60){
                std::ofstream WriteDataIn("spettro_simulato_6mm.out", std::ios::app);
                   WriteDataIn
                    <<   i              << " " //  1
                    <<   eventID        << " " //  2
                    <<   trackID        << " " //  3
                    <<   parentID       << " " //  4
                    <<   pdg            << " " //  5
                    <<   Z              << " " //  6
                    <<   A              << " " //  7
                    <<   kineticEnergy  << " " //  8
                    <<   init_pos[0]    << " " //  9
                    <<   init_pos[1]    << " " //  10
                    <<   init_pos[2]    << " " //  11
                    <<   process_name          //  12
                    <<   G4endl;
        }
    
    if (i==70){
                std::ofstream WriteDataIn("spettro_simulato_7mm.out", std::ios::app);
                   WriteDataIn
                    <<   i              << " " //  1
                    <<   eventID        << " " //  2
                    <<   trackID        << " " //  3
                    <<   parentID       << " " //  4
                    <<   pdg            << " " //  5
                    <<   Z              << " " //  6
                    <<   A              << " " //  7
                    <<   kineticEnergy  << " " //  8
                    <<   init_pos[0]    << " " //  9
                    <<   init_pos[1]    << " " //  10
                    <<   init_pos[2]    << " " //  11
                    <<   process_name          //  12
                    <<   G4endl;
        }
    
    if (i==80){
                std::ofstream WriteDataIn("spettro_simulato_8mm.out", std::ios::app);
                   WriteDataIn
                    <<   i              << " " //  1
                    <<   eventID        << " " //  2
                    <<   trackID        << " " //  3
                    <<   parentID       << " " //  4
                    <<   pdg            << " " //  5
                    <<   Z              << " " //  6
                    <<   A              << " " //  7
                    <<   kineticEnergy  << " " //  8
                    <<   init_pos[0]    << " " //  9
                    <<   init_pos[1]    << " " //  10
                    <<   init_pos[2]    << " " //  11
                    <<   process_name          //  12
                    <<   G4endl;
        }
             
             
             
    if (i==90){
                std::ofstream WriteDataIn("spettro_simulato_9mm.out", std::ios::app);
                   WriteDataIn
                    <<   i              << " " //  1
                    <<   eventID        << " " //  2
                    <<   trackID        << " " //  3
                    <<   parentID       << " " //  4
                    <<   pdg            << " " //  5
                    <<   Z              << " " //  6
                    <<   A              << " " //  7
                    <<   kineticEnergy  << " " //  8
                    <<   init_pos[0]    << " " //  9
                    <<   init_pos[1]    << " " //  10
                    <<   init_pos[2]    << " " //  11
                    <<   process_name          //  12
                    <<   G4endl;
        }

         
  l.unlock(); // LOCK STOP
  }
     */
    //  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
    
    //  ******************** let ***************************
    if (let)
    {
        if ( !(Z==0 && A==1) ) // All but not neutrons
        //if ( !(Z==0 && A==1) && Z!=9 ) // All but not neutrons and Flourine ******* SERENA
        // if (Z==1 || Z==2 || Z==3 || Z==4 || Z==5 || Z==6 || Z==7 || Z==8 || Z==9 || Z==10) // Only H, He, Li, Be, B, C, N, O, F, Ne ******* SERENA
        // if (Z==2 || Z==9) //  He & F  ******* SERENA
        //if (Z==10) // Only Ne ******* SERENA
        {
            if( energyDeposit>0. && DX >0. )// calculate only energy deposit
            {
                if (pdg !=22 && pdg !=11) // not gamma and electrons
                {
                    
                    // Get the pre-step kinetic energy
                    G4double eKinPre = aStep -> GetPreStepPoint() -> GetKineticEnergy();
                    // Get the post-step kinetic energy
                    G4double eKinPost = aStep -> GetPostStepPoint() -> GetKineticEnergy();
                    // Get the step average kinetic energy
                    G4double eKinMean = (eKinPre + eKinPost) * 0.5;
                    
                    // get the material
                    G4Material * materialStep = aStep -> GetPreStepPoint() -> GetMaterial();
                    
                    // get the secondary paticles
                    G4Step fstep = *theTrack -> GetStep();
                    // store all the secondary partilce in current step
                    const std::vector<const G4Track*> * secondary = fstep.GetSecondaryInCurrentStep();
                    
                    size_t SecondarySize = (*secondary).size();
                    G4double EnergySecondary = 0.;
                    
                    // get secondary electrons energy deposited
                    if (SecondarySize) // calculate only secondary particles
                    {
                        for (size_t numsec = 0; numsec< SecondarySize ; numsec ++)
                        {
                            //Get the PDG code of every secondaty particles in current step
                            G4int PDGSecondary=(*secondary)[numsec]->GetDefinition()->GetPDGEncoding();
                            
                            if(PDGSecondary == 11) // calculate only secondary electrons
                            {
                                // calculate the energy deposit of secondary electrons in current step
                                EnergySecondary += (*secondary)[numsec]->GetKineticEnergy();
                            }
                        }
                        
                    }
                    
                    // call the let filldatas function to calculate let
                    G4AutoLock l(&dataManipulationMutex); // LOCK START
                    let -> FillEnergySpectrum(trackID, particleDef, eKinMean, materialStep,
                                              energyDeposit,EnergySecondary,DX, i, j, k);
                    l.unlock(); // LOCK STOP
                }
            }
        }
    }
    
    
    if (matrix)
    {
        G4AutoLock l(&dataManipulationMutex);  // LOCK START
        
        // Increment Fluences & accumulate energy spectra
        // Hit voxels are marked with track_id throught hitTrack matrix
        //G4int* hitTrack = matrix -> GetHitTrack(i,j,k); // hitTrack MUST BE cleared at every eventAction!
        if ( *hitTrack != trackID )
        {
            *hitTrack = trackID;
            /*
             * Fill FLUENCE data for every single nuclide
             * Exclude e-, neutrons, gamma, ...
             */
            if ( Z >= 1) matrix -> Fill(trackID, particleDef, i, j, k, 0, true);
            
        }
        
        if(energyDeposit != 0)
        {
            /*
             *  This method will fill a dose matrix for every single nuclide.
             *  A method of the HadrontherapyMatrix class (StoreDoseFluenceAscii())
             *  is called automatically at the end of main (or via the macro command /analysis/writeDoseFile.
             *  It permits to store all dose/fluence data into a single plane ASCII file.
             */
            // if (A==1 && Z==1) // primary and sec. protons
            if ( Z>=1 )    //  exclude e-, neutrons, gamma, ...
                matrix -> Fill(trackID, particleDef, i, j, k, energyDeposit);
            /*
             * Create a hit with the information of position is in the detector
             */
            HadrontherapyDetectorHit* detectorHit = new HadrontherapyDetectorHit();
            detectorHit -> SetEdepAndPosition(i, j, k, energyDeposit);
            HitsCollection -> insert(detectorHit);
        }
        l.unlock(); // LOCK STOP
    }

    auto rbe = HadrontherapyRBE::GetInstance();
    if (rbe->IsCalculationEnabled())
    {
        if (!fRBEAccumulable)
        {
            fRBEAccumulable = dynamic_cast<HadrontherapyRBEAccumulable*>(G4AccumulableManager::Instance()->GetAccumulable("RBE"));
            if (!fRBEAccumulable)
            {
                G4Exception("HadrontherapyDetectorSD::ProcessHits", "NoAccumulable", FatalException, "Accumulable RBE not found.");
            }
        }
	if (A>0) //protect against gammas, e- , etc
	  {
	    fRBEAccumulable->Accumulate(kineticEnergy / A, energyDeposit, DX, Z, i, j, k);
	  }
    }


    return true;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

