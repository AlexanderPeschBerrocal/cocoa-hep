#include "Run_metadata.hh"

#include "CocoaVersion.hh" // generated at build time, see cmake/CocoaGitVersion.cmake
#include "Config_reader_var.hh"
#include "GeometryMaterialWriter.hh"
#include "Truth_record_data.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Version.hh"

#include "Pythia8/Pythia.h"

#include "TMD5.h"
#include "TROOT.h"
#include "TTree.h"

#include <unistd.h>

#include <cstdio>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace
{
std::string g_start_time_utc;

std::string NowUtc()
{
    const std::time_t now = std::time(nullptr);
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    return buffer;
}

std::string HostName()
{
    char buffer[256] = {0};
    if (gethostname(buffer, sizeof(buffer) - 1) != 0)
        return "unknown";
    return buffer;
}

std::string ReadFile(const std::string &path)
{
    if (path.empty())
        return "";
    std::ifstream in(path);
    if (!in)
        return "";
    std::stringstream text;
    text << in.rdbuf();
    return text.str();
}

// MD5 over every placed volume: path, copy number, placement, solid parameters
// and material. Identical geometries give identical fingerprints.
class GeometryFingerprint
{
public:
    std::string Compute(const G4VPhysicalVolume *world)
    {
        Add(world, "");
        m_md5.Final();
        return m_md5.AsString();
    }

private:
    void Add(const G4VPhysicalVolume *pv, const std::string &parent_path)
    {
        const G4LogicalVolume *lv = pv->GetLogicalVolume();
        const std::string path = parent_path + "/" + pv->GetName() + "#" + std::to_string(pv->GetCopyNo());
        const G4ThreeVector t = pv->GetTranslation();
        const G4RotationMatrix r = pv->GetObjectRotationValue();

        char placement[256];
        std::snprintf(placement, sizeof(placement), "|%.6g %.6g %.6g|%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g|",
                      t.x(), t.y(), t.z(), r.xx(), r.xy(), r.xz(), r.yx(), r.yy(), r.yz(), r.zx(), r.zy(), r.zz());

        std::string record = path + placement + lv->GetName() + "|" +
                             (lv->GetMaterial() ? std::string(lv->GetMaterial()->GetName()) : std::string()) + "|" +
                             SolidDescription(lv->GetSolid()) + "\n";
        m_md5.Update(reinterpret_cast<const UChar_t *>(record.data()), record.size());

        for (size_t i = 0; i < lv->GetNoDaughters(); ++i)
            Add(lv->GetDaughter(i), path);
    }

    const std::string &SolidDescription(G4VSolid *solid)
    {
        auto cached = m_solids.find(solid);
        if (cached != m_solids.end())
            return cached->second;
        std::ostringstream description;
        description.precision(6);
        solid->StreamInfo(description);
        return m_solids.emplace(solid, description.str()).first->second;
    }

    TMD5 m_md5;
    std::unordered_map<const G4VSolid *, std::string> m_solids;
};
} // namespace

void Run_metadata::BeginRun()
{
    g_start_time_utc = NowUtc();
}

void Run_metadata::Write(int n_events_processed)
{
    Config_reader_var &config = Config_reader_var::GetInstance();
    Truth_record_data &truth = Truth_record_data::GetInstance();

    // ---- RunInfo ----------------------------------------------------------
    int schema_version = kSchemaVersion;
    std::string git_commit = COCOA_GIT_COMMIT;
    std::string git_describe = COCOA_GIT_DESCRIBE;
    int git_dirty = COCOA_GIT_DIRTY;
    std::string geant4_version = G4Version;
    std::string root_version = gROOT->GetVersion();
    char pythia_buffer[16];
    std::snprintf(pythia_buffer, sizeof(pythia_buffer), "%.3f", static_cast<double>(PYTHIA_VERSION));
    std::string pythia8_version = pythia_buffer;
    std::string physics_list = "FTFP_BERT";
    std::string command_line = config.command_line;
    std::string config_path = config.config_path;
    std::string config_json = config.config_json_text;
    std::string macro_path = config.macro_path_used;
    std::string macro_text = ReadFile(config.macro_path_used);
    std::string output_path = config.Output_file_path;
    Long64_t run_seed = config.run_seed;
    int n_events_requested = config.n_events_requested;
    int n_processed = n_events_processed;
    std::string type_of_running = config.Type_of_running;

    const G4VPhysicalVolume *world = G4TransportationManager::GetTransportationManager()
                                         ->GetNavigatorForTracking()
                                         ->GetWorldVolume();
    std::string geometry_fingerprint = world ? GeometryFingerprint().Compute(world) : "";

    float tracker_r_max = 0.f;
    float tracker_abs_z_max = 0.f;
    if (const G4LogicalVolume *inner = G4LogicalVolumeStore::GetInstance()->GetVolume("Inner_LV", false))
    {
        if (const auto *tubs = dynamic_cast<const G4Tubs *>(inner->GetSolid()))
        {
            tracker_r_max = tubs->GetOuterRadius() / mm;
            tracker_abs_z_max = tubs->GetZHalfLength() / mm;
        }
    }

    const Truth_record_var &record = config.truth_record;
    int truth_secondary_conversions = record.record_secondary_conversions;
    float truth_conversion_e_min = record.conversion_secondary_photon_e_min;
    int truth_secondary_nuclear = record.record_secondary_nuclear_interactions;
    float truth_nuclear_hadron_ekin_min = record.nuclear_secondary_hadron_ekin_min;
    float truth_nuclear_daughter_ekin_min = record.nuclear_daughter_ekin_min;

    float fiducial_pt_min = config.fiducial_cuts.pt_min;
    float fiducial_e_min = config.fiducial_cuts.e_min;
    float fiducial_eta_max = config.fiducial_cuts.eta_max;
    float fiducial_dr_cut = config.fiducial_cuts.dR_cut;
    float magnetic_field = config.fieldValue / tesla;
    int use_inner_detector = config.use_inner_detector;
    int use_detailed_tracker_material = config.use_detailed_tracker_material;

    std::string start_time = g_start_time_utc;
    std::string end_time = NowUtc();
    std::string host = HostName();

    TTree run_info("RunInfo", "COCOA run metadata (docs/output_schema.md)");
    run_info.Branch("schema_version", &schema_version, "schema_version/I");
    run_info.Branch("cocoa_git_commit", &git_commit);
    run_info.Branch("cocoa_git_describe", &git_describe);
    run_info.Branch("cocoa_git_dirty", &git_dirty, "cocoa_git_dirty/I");
    run_info.Branch("geant4_version", &geant4_version);
    run_info.Branch("root_version", &root_version);
    run_info.Branch("pythia8_version", &pythia8_version);
    run_info.Branch("physics_list", &physics_list);
    run_info.Branch("command_line", &command_line);
    run_info.Branch("config_path", &config_path);
    run_info.Branch("config_json", &config_json);
    run_info.Branch("macro_path", &macro_path);
    run_info.Branch("macro_text", &macro_text);
    run_info.Branch("output_path", &output_path);
    run_info.Branch("run_seed", &run_seed, "run_seed/L");
    run_info.Branch("n_events_requested", &n_events_requested, "n_events_requested/I");
    run_info.Branch("n_events_processed", &n_processed, "n_events_processed/I");
    run_info.Branch("type_of_running", &type_of_running);
    run_info.Branch("geometry_fingerprint", &geometry_fingerprint);
    run_info.Branch("tracker_region_r_max_mm", &tracker_r_max, "tracker_region_r_max_mm/F");
    run_info.Branch("tracker_region_abs_z_max_mm", &tracker_abs_z_max, "tracker_region_abs_z_max_mm/F");
    run_info.Branch("truth_record_secondary_conversions", &truth_secondary_conversions, "truth_record_secondary_conversions/I");
    run_info.Branch("truth_record_conversion_secondary_photon_e_min_mev", &truth_conversion_e_min, "truth_record_conversion_secondary_photon_e_min_mev/F");
    run_info.Branch("truth_record_secondary_nuclear", &truth_secondary_nuclear, "truth_record_secondary_nuclear/I");
    run_info.Branch("truth_record_nuclear_secondary_hadron_ekin_min_mev", &truth_nuclear_hadron_ekin_min, "truth_record_nuclear_secondary_hadron_ekin_min_mev/F");
    run_info.Branch("truth_record_nuclear_daughter_ekin_min_mev", &truth_nuclear_daughter_ekin_min, "truth_record_nuclear_daughter_ekin_min_mev/F");
    run_info.Branch("fiducial_pt_min_gev", &fiducial_pt_min, "fiducial_pt_min_gev/F");
    run_info.Branch("fiducial_e_min_gev", &fiducial_e_min, "fiducial_e_min_gev/F");
    run_info.Branch("fiducial_eta_max", &fiducial_eta_max, "fiducial_eta_max/F");
    run_info.Branch("fiducial_dr_cut", &fiducial_dr_cut, "fiducial_dr_cut/F");
    run_info.Branch("magnetic_field_tesla", &magnetic_field, "magnetic_field_tesla/F");
    run_info.Branch("use_inner_detector", &use_inner_detector, "use_inner_detector/I");
    run_info.Branch("use_detailed_tracker_material", &use_detailed_tracker_material, "use_detailed_tracker_material/I");
    run_info.Branch("start_time_utc", &start_time);
    run_info.Branch("end_time_utc", &end_time);
    run_info.Branch("host", &host);
    run_info.Fill();
    run_info.Write();

    // ---- Materials --------------------------------------------------------
    int material_idx = -1;
    std::string material_name;
    float density = 0.f, x0 = 0.f, lambda = 0.f, z_eff = 0.f, a_eff = 0.f;
    TTree materials("Materials", "Geant4 materials, target of *_material_idx");
    materials.Branch("material_idx", &material_idx, "material_idx/I");
    materials.Branch("name", &material_name);
    materials.Branch("density_g_cm3", &density, "density_g_cm3/F");
    materials.Branch("x0_mm", &x0, "x0_mm/F");
    materials.Branch("lambda_mm", &lambda, "lambda_mm/F");
    materials.Branch("z_eff", &z_eff, "z_eff/F");
    materials.Branch("a_eff_g_mol", &a_eff, "a_eff_g_mol/F");
    for (const G4Material *material : *G4Material::GetMaterialTable())
    {
        material_idx = static_cast<int>(material->GetIndex());
        material_name = material->GetName();
        density = material->GetDensity() / (g / cm3);
        x0 = material->GetRadlen() / mm;
        lambda = material->GetNuclearInterLength() / mm;
        double z = 0.0, a = 0.0;
        GeometryMaterialWriter::EffectiveZA(material, z, a);
        z_eff = z;
        a_eff = a;
        materials.Fill();
    }
    materials.Write();

    // ---- Volumes ----------------------------------------------------------
    int volume_idx = -1;
    std::string volume_name;
    int volume_material_idx = -1;
    int in_tracker = 0;
    TTree volumes("Volumes", "Geant4 logical volumes, target of *_volume_idx");
    volumes.Branch("volume_idx", &volume_idx, "volume_idx/I");
    volumes.Branch("name", &volume_name);
    volumes.Branch("material_idx", &volume_material_idx, "material_idx/I");
    volumes.Branch("in_tracker", &in_tracker, "in_tracker/I");
    const G4LogicalVolumeStore *store = G4LogicalVolumeStore::GetInstance();
    for (size_t i = 0; i < store->size(); ++i)
    {
        const G4LogicalVolume *lv = (*store)[i];
        volume_idx = truth.VolumeIndex(lv);
        volume_name = lv->GetName();
        volume_material_idx = lv->GetMaterial() ? static_cast<int>(lv->GetMaterial()->GetIndex()) : -1;
        in_tracker = truth.InTracker(lv) ? 1 : 0;
        volumes.Fill();
    }
    volumes.Write();
}
