#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "OutputRunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

#include <utility>

ActionInitialization::ActionInitialization(Geometry_definition geometry,
                                           std::string output_file,
                                           bool save_truth_graph)
    : geometry_(std::move(geometry)),
      output_file_(std::move(output_file)),
      save_truth_graph_(save_truth_graph)
{
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new OutputRunAction(output_file_, save_truth_graph_));
    SetUserAction(new EventAction);
    SetUserAction(new TrackingAction);
    SetUserAction(new SteppingAction(geometry_));
}
