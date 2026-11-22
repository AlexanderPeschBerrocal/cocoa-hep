#include "Truth_record_data.hh"

#include "Config_reader_var.hh"
#include "Full_trajectory_info_data.hh"

#include "G4DynamicParticle.hh"
#include "G4EmProcessSubType.hh"
#include "G4HadronicProcessType.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <cmath>
#include <limits>
#include <string>

namespace
{
const char *kTrackerRegionName = "TrackerRegion";
const float kNaN = std::numeric_limits<float>::quiet_NaN();
// Daughters are matched to their nuclear-interaction vertex by position.
const double kVertexMatchTolerance = 1e-6; // mm

bool IsHadron(const G4ParticleDefinition *definition)
{
    const G4String &type = definition->GetParticleType();
    return type == "baryon" || type == "meson" || type == "nucleus";
}

bool IsScatter(int type, int subtype)
{
    return (type == fElectromagnetic && (subtype == fComptonScattering || subtype == fRayleigh)) ||
           (type == fHadronic && subtype == fHadronElastic);
}

bool IsNuclearInteraction(int type, int subtype)
{
    return type == fHadronic &&
           (subtype == fHadronInelastic || subtype == fCapture || subtype == fFission);
}
} // namespace

Truth_outgoing::Truth_outgoing()
    : calo_e(kNaN), calo_x(kNaN), calo_y(kNaN), calo_z(kNaN)
{
}

// ---------------------------------------------------------------------------
// Per-event bookkeeping and output
// ---------------------------------------------------------------------------

void Truth_record_data::clear()
{
    conv_vtx.clear();
    conv_el.clear();
    nucl_vtx.clear();
    nucl_dau.clear();

    m_current_track_id = -1;
    m_x0 = 0.0;
    m_lambda = 0.0;
    m_n_scatters = 0;

    m_primary_row.clear();
    m_ancestor_row.clear();
    m_pending_conv.clear();
    m_pending_nucl.clear();
    m_exit_watch.clear();

    m_conv_vtx_branches.clear();
    m_conv_vtx_ep_idx.clear();
    m_conv_vtx_em_idx.clear();
    m_conv_el_branches.clear();
    m_nucl_vtx_branches.clear();
    m_nucl_vtx_process_subtype.clear();
    m_nucl_vtx_n_dau_total.clear();
    m_nucl_vtx_n_dau_charged.clear();
    m_nucl_dau_branches.clear();
}

void Truth_record_data::SetEventInfo(int event_id, long long run_seed)
{
    m_event_id = event_id;
    m_run_seed = run_seed;
}

void Truth_record_data::set_event_tree_branches(TTree *outTree)
{
    outTree->Branch("event_id", &m_event_id, "event_id/I");
    outTree->Branch("run_seed", &m_run_seed, "run_seed/L");
}

void Truth_record_data::set_tree_branches(TTree *outTree)
{
    m_conv_vtx_branches.set_branches(outTree, "conv_vtx");
    outTree->Branch("conv_vtx_ep_idx", "vector<int>", &m_conv_vtx_ep_idx);
    outTree->Branch("conv_vtx_em_idx", "vector<int>", &m_conv_vtx_em_idx);
    m_conv_el_branches.set_branches(outTree, "conv_el");

    m_nucl_vtx_branches.set_branches(outTree, "nucl_vtx");
    outTree->Branch("nucl_vtx_process_subtype", "vector<int>", &m_nucl_vtx_process_subtype);
    outTree->Branch("nucl_vtx_n_dau_total", "vector<int>", &m_nucl_vtx_n_dau_total);
    outTree->Branch("nucl_vtx_n_dau_charged", "vector<int>", &m_nucl_vtx_n_dau_charged);
    m_nucl_dau_branches.set_branches(outTree, "nucl_dau");
}

void Truth_record_data::fill_var()
{
    m_conv_vtx_branches.fill(conv_vtx);
    for (const Truth_vertex &v : conv_vtx)
    {
        m_conv_vtx_ep_idx.push_back(v.ep_idx);
        m_conv_vtx_em_idx.push_back(v.em_idx);
    }
    m_conv_el_branches.fill(conv_el);

    m_nucl_vtx_branches.fill(nucl_vtx);
    for (const Truth_vertex &v : nucl_vtx)
    {
        m_nucl_vtx_process_subtype.push_back(v.process_subtype);
        m_nucl_vtx_n_dau_total.push_back(v.n_dau_total);
        m_nucl_vtx_n_dau_charged.push_back(v.n_dau_charged);
    }
    m_nucl_dau_branches.fill(nucl_dau);
}

void Truth_record_data::Vertex_branches::clear()
{
    for (auto *v : {&x, &y, &z, &in_e, &in_px, &in_py, &in_pz, &in_prod_x, &in_prod_y, &in_prod_z, &in_x0, &in_lambda})
        v->clear();
    for (auto *v : {&in_pdgid, &in_particle_idx, &in_ancestor_idx, &in_creator_type, &in_creator_subtype,
                    &in_g4_track_id, &material_idx, &volume_idx})
        v->clear();
}

void Truth_record_data::Vertex_branches::fill(const std::vector<Truth_vertex> &rows)
{
    for (const Truth_vertex &v : rows)
    {
        x.push_back(v.x);
        y.push_back(v.y);
        z.push_back(v.z);
        in_pdgid.push_back(v.in.pdgid);
        in_e.push_back(v.in.e);
        in_px.push_back(v.in.px);
        in_py.push_back(v.in.py);
        in_pz.push_back(v.in.pz);
        in_particle_idx.push_back(v.in.particle_idx);
        in_ancestor_idx.push_back(v.in.ancestor_idx);
        in_creator_type.push_back(v.in.creator_type);
        in_creator_subtype.push_back(v.in.creator_subtype);
        in_prod_x.push_back(v.in.prod_x);
        in_prod_y.push_back(v.in.prod_y);
        in_prod_z.push_back(v.in.prod_z);
        in_x0.push_back(v.in.x0);
        in_lambda.push_back(v.in.lambda);
        in_g4_track_id.push_back(v.in.g4_track_id);
        material_idx.push_back(v.material_idx);
        volume_idx.push_back(v.volume_idx);
    }
}

void Truth_record_data::Vertex_branches::set_branches(TTree *tree, const std::string &prefix)
{
    auto floats = [&](const char *name, std::vector<float> *v) {
        tree->Branch((prefix + "_" + name).c_str(), "vector<float>", v);
    };
    auto ints = [&](const char *name, std::vector<int> *v) {
        tree->Branch((prefix + "_" + name).c_str(), "vector<int>", v);
    };
    floats("x", &x);
    floats("y", &y);
    floats("z", &z);
    ints("in_pdgid", &in_pdgid);
    floats("in_e", &in_e);
    floats("in_px", &in_px);
    floats("in_py", &in_py);
    floats("in_pz", &in_pz);
    ints("in_particle_idx", &in_particle_idx);
    ints("in_ancestor_idx", &in_ancestor_idx);
    ints("in_creator_type", &in_creator_type);
    ints("in_creator_subtype", &in_creator_subtype);
    floats("in_prod_x", &in_prod_x);
    floats("in_prod_y", &in_prod_y);
    floats("in_prod_z", &in_prod_z);
    floats("in_x0", &in_x0);
    floats("in_lambda", &in_lambda);
    ints("in_g4_track_id", &in_g4_track_id);
    ints("material_idx", &material_idx);
    ints("volume_idx", &volume_idx);
}

void Truth_record_data::Outgoing_branches::clear()
{
    for (auto *v : {&q, &e, &px, &py, &pz, &calo_e, &calo_x, &calo_y, &calo_z})
        v->clear();
    for (auto *v : {&vtx_idx, &pdgid, &g4_track_id, &calo_reached})
        v->clear();
}

void Truth_record_data::Outgoing_branches::fill(const std::vector<Truth_outgoing> &rows)
{
    for (const Truth_outgoing &o : rows)
    {
        vtx_idx.push_back(o.vtx_idx);
        pdgid.push_back(o.pdgid);
        q.push_back(o.q);
        e.push_back(o.e);
        px.push_back(o.px);
        py.push_back(o.py);
        pz.push_back(o.pz);
        g4_track_id.push_back(o.g4_track_id);
        calo_reached.push_back(o.calo_reached);
        calo_e.push_back(o.calo_e);
        calo_x.push_back(o.calo_x);
        calo_y.push_back(o.calo_y);
        calo_z.push_back(o.calo_z);
    }
}

void Truth_record_data::Outgoing_branches::set_branches(TTree *tree, const std::string &prefix)
{
    auto floats = [&](const char *name, std::vector<float> *v) {
        tree->Branch((prefix + "_" + name).c_str(), "vector<float>", v);
    };
    auto ints = [&](const char *name, std::vector<int> *v) {
        tree->Branch((prefix + "_" + name).c_str(), "vector<int>", v);
    };
    ints("vtx_idx", &vtx_idx);
    ints("pdgid", &pdgid);
    floats("q", &q);
    floats("e", &e);
    floats("px", &px);
    floats("py", &py);
    floats("pz", &pz);
    ints("g4_track_id", &g4_track_id);
    ints("calo_reached", &calo_reached);
    floats("calo_e", &calo_e);
    floats("calo_x", &calo_x);
    floats("calo_y", &calo_y);
    floats("calo_z", &calo_z);
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

const G4Region *Truth_record_data::TrackerRegion()
{
    if (!m_tracker_region)
        m_tracker_region = G4RegionStore::GetInstance()->GetRegion(kTrackerRegionName, false);
    return m_tracker_region;
}

bool Truth_record_data::InTracker(const G4LogicalVolume *lv)
{
    const G4Region *region = TrackerRegion();
    return lv && region && lv->GetRegion() == region;
}

int Truth_record_data::VolumeIndex(const G4LogicalVolume *lv)
{
    if (m_volume_index.empty())
    {
        const G4LogicalVolumeStore *store = G4LogicalVolumeStore::GetInstance();
        for (size_t i = 0; i < store->size(); ++i)
            m_volume_index[(*store)[i]] = static_cast<int>(i);
    }
    auto it = m_volume_index.find(lv);
    return it != m_volume_index.end() ? it->second : -1;
}

// ---------------------------------------------------------------------------
// Tracking
// ---------------------------------------------------------------------------

void Truth_record_data::RegisterPrimary(int g4_track_id, int particle_row)
{
    m_primary_row[g4_track_id] = particle_row;
}

int Truth_record_data::BeginTrack(const G4Track *track)
{
    const int track_id = track->GetTrackID();
    const int parent_id = track->GetParentID();

    m_current_track_id = track_id;
    m_x0 = 0.0;
    m_lambda = 0.0;
    m_n_scatters = 0;

    int ancestor = -1;
    if (parent_id == 0)
    {
        const auto primary = m_primary_row.find(track_id);
        if (primary != m_primary_row.end())
            ancestor = primary->second;
    }
    else
    {
        const auto parent = m_ancestor_row.find(parent_id);
        if (parent != m_ancestor_row.end())
            ancestor = parent->second;
    }
    m_ancestor_row[track_id] = ancestor;

    if (parent_id == 0)
        return -1;
    const G4VProcess *creator = track->GetCreatorProcess();
    if (!creator)
        return -1;

    // Leg of a recorded conversion?
    const auto conv = m_pending_conv.find(parent_id);
    if (conv != m_pending_conv.end() &&
        creator->GetProcessType() == fElectromagnetic &&
        creator->GetProcessSubType() == fGammaConversion &&
        std::abs(track->GetDefinition()->GetPDGEncoding()) == 11)
    {
        const int row = static_cast<int>(conv_el.size());
        Truth_outgoing out;
        FillOutgoing(track, conv->second, out);
        conv_el.push_back(out);

        Truth_vertex &vertex = conv_vtx[conv->second];
        if (out.q > 0 && vertex.ep_idx < 0)
            vertex.ep_idx = row;
        else if (out.q < 0 && vertex.em_idx < 0)
            vertex.em_idx = row;

        m_exit_watch[track_id] = {Watched::kConvEl, row};
        return row;
    }

    // Daughter of a recorded nuclear interaction?
    const auto nucl = m_pending_nucl.find(parent_id);
    if (nucl != m_pending_nucl.end())
    {
        const G4ThreeVector &vertex_position = track->GetVertexPosition();
        for (const Pending_nuclear &pending : nucl->second)
        {
            if (pending.process != creator ||
                (pending.position - vertex_position).mag() > kVertexMatchTolerance)
                continue;
            if (track->GetKineticEnergy() >= Config_reader_var::GetInstance().truth_record.nuclear_daughter_ekin_min)
            {
                const int row = static_cast<int>(nucl_dau.size());
                Truth_outgoing out;
                FillOutgoing(track, pending.row, out);
                nucl_dau.push_back(out);
                m_exit_watch[track_id] = {Watched::kNuclDau, row};
            }
            break;
        }
    }
    return -1;
}

void Truth_record_data::Step(const G4Step *step)
{
    const G4StepPoint *pre = step->GetPreStepPoint();
    const G4VPhysicalVolume *pre_pv = pre->GetPhysicalVolume();
    if (!pre_pv || !InTracker(pre_pv->GetLogicalVolume()))
        return;

    const G4Track *track = step->GetTrack();
    if (track->GetTrackID() == m_current_track_id)
    {
        const G4Material *material = pre->GetMaterial();
        const double length = step->GetStepLength();
        m_x0 += length / material->GetRadlen();
        m_lambda += length / material->GetNuclearInterLength();
    }

    const G4StepPoint *post = step->GetPostStepPoint();
    const G4VProcess *process = post->GetProcessDefinedStep();
    if (process)
    {
        const int type = process->GetProcessType();
        const int subtype = process->GetProcessSubType();
        const G4ParticleDefinition *definition = track->GetDefinition();
        if (IsScatter(type, subtype))
            ++m_n_scatters;
        else if (type == fElectromagnetic && subtype == fGammaConversion && definition->GetPDGEncoding() == 22)
            MakeConversionVertex(step);
        else if (IsNuclearInteraction(type, subtype) && IsHadron(definition))
            MakeNuclearVertex(step, process);
    }

    // Outgoing particles leaving the tracker region, i.e. reaching the
    // calorimeter front face.
    if (!m_exit_watch.empty())
    {
        const G4VPhysicalVolume *post_pv = post->GetPhysicalVolume();
        if (!post_pv || !InTracker(post_pv->GetLogicalVolume()))
        {
            const auto watched = m_exit_watch.find(track->GetTrackID());
            if (watched != m_exit_watch.end())
            {
                Truth_outgoing &out = watched->second.first == Watched::kConvEl
                                          ? conv_el[watched->second.second]
                                          : nucl_dau[watched->second.second];
                const G4ThreeVector &position = post->GetPosition();
                out.calo_reached = 1;
                out.calo_e = post->GetTotalEnergy();
                out.calo_x = position.x();
                out.calo_y = position.y();
                out.calo_z = position.z();
                m_exit_watch.erase(watched);
            }
        }
    }
}

void Truth_record_data::EndPrimaryTrack(const G4Track *track, FullTrajectoryInfo &primary)
{
    primary.fEndPosition = track->GetPosition();
    primary.fTrackLength = track->GetTrackLength();

    const G4Step *step = track->GetStep();
    if (step)
    {
        const G4StepPoint *pre = step->GetPreStepPoint();
        const G4StepPoint *post = step->GetPostStepPoint();
        primary.fEndEnergy = pre->GetTotalEnergy();
        primary.fEndMomentum = pre->GetMomentum();

        const G4VProcess *process = post->GetProcessDefinedStep();
        primary.fEndProcessType = process ? process->GetProcessType() : -1;
        primary.fEndProcessSubType = process ? process->GetProcessSubType() : -1;

        const G4VPhysicalVolume *pre_pv = pre->GetPhysicalVolume();
        const G4LogicalVolume *lv = pre_pv ? pre_pv->GetLogicalVolume() : nullptr;
        primary.fEndInTracker = InTracker(lv) ? 1 : 0;
        primary.fEndMaterialIdx = pre->GetMaterial() ? static_cast<int>(pre->GetMaterial()->GetIndex()) : -1;
        primary.fEndVolumeIdx = VolumeIndex(lv);
    }

    if (track->GetTrackID() == m_current_track_id)
    {
        primary.fX0Tracker = m_x0;
        primary.fLambdaTracker = m_lambda;
        primary.fNScattersTracker = m_n_scatters;
    }
}

// ---------------------------------------------------------------------------
// Vertices
// ---------------------------------------------------------------------------

Truth_incoming Truth_record_data::MakeIncoming(const G4Track *track, const G4StepPoint *pre) const
{
    Truth_incoming in;
    const int track_id = track->GetTrackID();
    in.pdgid = track->GetDefinition()->GetPDGEncoding();
    in.e = pre->GetTotalEnergy();
    const G4ThreeVector momentum = pre->GetMomentum();
    in.px = momentum.x();
    in.py = momentum.y();
    in.pz = momentum.z();

    if (track->GetParentID() == 0)
    {
        const auto primary = m_primary_row.find(track_id);
        in.particle_idx = primary != m_primary_row.end() ? primary->second : -1;
    }
    const auto ancestor = m_ancestor_row.find(track_id);
    in.ancestor_idx = ancestor != m_ancestor_row.end() ? ancestor->second : -1;

    if (const G4VProcess *creator = track->GetCreatorProcess())
    {
        in.creator_type = creator->GetProcessType();
        in.creator_subtype = creator->GetProcessSubType();
    }

    const G4ThreeVector &production = track->GetVertexPosition();
    in.prod_x = production.x();
    in.prod_y = production.y();
    in.prod_z = production.z();
    if (track_id == m_current_track_id)
    {
        in.x0 = m_x0;
        in.lambda = m_lambda;
    }
    in.g4_track_id = track_id;
    return in;
}

void Truth_record_data::MakeConversionVertex(const G4Step *step)
{
    const G4Track *track = step->GetTrack();
    const G4StepPoint *pre = step->GetPreStepPoint();
    const Truth_record_var &settings = Config_reader_var::GetInstance().truth_record;
    if (track->GetParentID() != 0 &&
        (!settings.record_secondary_conversions ||
         pre->GetTotalEnergy() < settings.conversion_secondary_photon_e_min))
        return;

    Truth_vertex vertex;
    const G4ThreeVector &position = step->GetPostStepPoint()->GetPosition();
    vertex.x = position.x();
    vertex.y = position.y();
    vertex.z = position.z();
    vertex.in = MakeIncoming(track, pre);
    vertex.material_idx = static_cast<int>(pre->GetMaterial()->GetIndex());
    vertex.volume_idx = VolumeIndex(pre->GetPhysicalVolume()->GetLogicalVolume());

    m_pending_conv[track->GetTrackID()] = static_cast<int>(conv_vtx.size());
    conv_vtx.push_back(vertex);
}

void Truth_record_data::MakeNuclearVertex(const G4Step *step, const G4VProcess *process)
{
    const G4Track *track = step->GetTrack();
    const G4StepPoint *pre = step->GetPreStepPoint();
    const Truth_record_var &settings = Config_reader_var::GetInstance().truth_record;
    if (track->GetParentID() != 0 &&
        (!settings.record_secondary_nuclear_interactions ||
         pre->GetKineticEnergy() < settings.nuclear_secondary_hadron_ekin_min))
        return;

    Truth_vertex vertex;
    const G4ThreeVector &position = step->GetPostStepPoint()->GetPosition();
    vertex.x = position.x();
    vertex.y = position.y();
    vertex.z = position.z();
    vertex.in = MakeIncoming(track, pre);
    vertex.material_idx = static_cast<int>(pre->GetMaterial()->GetIndex());
    vertex.volume_idx = VolumeIndex(pre->GetPhysicalVolume()->GetLogicalVolume());
    vertex.process_subtype = process->GetProcessSubType();

    if (const std::vector<const G4Track *> *secondaries = step->GetSecondaryInCurrentStep())
    {
        for (const G4Track *secondary : *secondaries)
        {
            if (secondary->GetCreatorProcess() != process)
                continue;
            ++vertex.n_dau_total;
            if (secondary->GetDynamicParticle()->GetCharge() != 0.)
                ++vertex.n_dau_charged;
        }
    }

    const int row = static_cast<int>(nucl_vtx.size());
    m_pending_nucl[track->GetTrackID()].push_back({row, process, position});
    nucl_vtx.push_back(vertex);
}

void Truth_record_data::FillOutgoing(const G4Track *track, int vtx_idx, Truth_outgoing &out)
{
    const G4DynamicParticle *particle = track->GetDynamicParticle();
    const G4ThreeVector momentum = particle->GetMomentum();
    out.vtx_idx = vtx_idx;
    out.pdgid = track->GetDefinition()->GetPDGEncoding();
    out.q = particle->GetCharge();
    out.e = particle->GetTotalEnergy();
    out.px = momentum.x();
    out.py = momentum.y();
    out.pz = momentum.z();
    out.g4_track_id = track->GetTrackID();
}
