#ifndef __PRIMARY_PARTICLE_INFORMATION_H__
#define __PRIMARY_PARTICLE_INFORMATION_H__

#include "G4VUserPrimaryParticleInformation.hh"
#include "G4ios.hh"

// Attached to every G4PrimaryParticle created from a generator particle. It
// holds the particle's index in TruthRecordGraph::m_final_state_particles, so
// that particle_* rows can be linked to the truth graph independently of the
// order in which Geant4 tracks the primaries.
class PrimaryParticleInformation : public G4VUserPrimaryParticleInformation
{
public:
    explicit PrimaryParticleInformation(int generator_idx) : fGeneratorIdx(generator_idx) {}
    int GetGeneratorIdx() const { return fGeneratorIdx; }
    void Print() const override { G4cout << "COCOA generator index " << fGeneratorIdx << G4endl; }

private:
    int fGeneratorIdx;
};

#endif // __PRIMARY_PARTICLE_INFORMATION_H__
