#ifndef COCOA_ACTION_INITIALIZATION_HH
#define COCOA_ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"
#include "Config_reader_var.hh"
#include <string>

class ActionInitialization final : public G4VUserActionInitialization
{
public:
    ActionInitialization(Geometry_definition geometry, std::string output_file,
                         bool save_truth_graph);
    void Build() const override;

private:
    Geometry_definition geometry_;
    std::string output_file_;
    bool save_truth_graph_;
};

#endif
