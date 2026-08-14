#include "Config_reader_func.hh"

#include <iostream>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <set>

namespace
{
void requireNumericArray(const Json::Value &value, const std::string &name)
{
    if (!value.isArray() || value.empty())
        throw std::runtime_error(name + " must be a non-empty array");
    for (Json::ArrayIndex i = 0; i < value.size(); ++i)
    {
        if (!value[i].isNumeric())
            throw std::runtime_error(name + " must contain only numbers");
    }
}

bool isPositiveInt(const Json::Value &value)
{
    if (value.isInt())
        return value.asInt() > 0;
    if (value.isUInt())
        return value.asUInt() > 0 &&
               value.asUInt() <= static_cast<unsigned int>(std::numeric_limits<int>::max());
    return false;
}

void requirePositiveIntegerArray(const Json::Value &value, const std::string &name)
{
    if (!value.isArray() || value.empty())
        throw std::runtime_error(name + " must be a non-empty array");
    for (Json::ArrayIndex i = 0; i < value.size(); ++i)
    {
        if (!isPositiveInt(value[i]))
            throw std::runtime_error(name + " must contain only positive integers");
    }
}

void requireBoundedNumericArray(const Json::Value &value, const std::string &name,
                                double lower_bound, bool allow_lower_bound)
{
    requireNumericArray(value, name);
    for (Json::ArrayIndex i = 0; i < value.size(); ++i)
    {
        const double number = value[i].asDouble();
        if (allow_lower_bound ? number < lower_bound : number <= lower_bound)
            throw std::runtime_error(name + " contains an out-of-range value");
    }
}

void requireNonEmptyStringArray(const Json::Value &value, const std::string &name)
{
    if (!value.isArray() || value.empty())
        throw std::runtime_error(name + " must be a non-empty array");
    for (Json::ArrayIndex i = 0; i < value.size(); ++i)
    {
        if (!value[i].isString() || value[i].asString().empty())
            throw std::runtime_error(name + " must contain only non-empty strings");
    }
}

void requireSameSize(const Json::Value &left, const Json::Value &right,
                     const std::string &description)
{
    if (left.size() != right.size())
        throw std::runtime_error("Mismatched configuration lengths: " + description);
}

void requireMatchingNumeric2D(const Json::Value &left, const Json::Value &right,
                              const std::string &description)
{
    if (!left.isArray() || left.empty() || !right.isArray() ||
        left.size() != right.size())
        throw std::runtime_error("Invalid high-granularity arrays: " + description);
    for (Json::ArrayIndex i = 0; i < left.size(); ++i)
    {
        requireNumericArray(left[i], description);
        requireNumericArray(right[i], description);
        requireSameSize(left[i], right[i], description);
    }
}
}

Config_reader_func::Config_reader_func(std::string path, Config_reader_var &config_var)
{
    
    std::ifstream config_doc(path, std::ifstream::binary);
    if (!config_doc.is_open())
        throw std::runtime_error("Cannot open configuration file: " + path);

    Json::CharReaderBuilder builder;
    builder["allowComments"] = false;
    builder["collectComments"] = false;
    builder["strictRoot"] = true;
    builder["failIfExtra"] = true;
    std::string parse_errors;
    if (!Json::parseFromStream(builder, config_doc, &configs, &parse_errors))
        throw std::runtime_error("Invalid JSON in " + path + ": " + parse_errors);

    if (!configs.isObject() || !configs["Geometry_definition"].isObject())
        throw std::runtime_error("Configuration requires a Geometry_definition object");
    
    config_var.Output_file_path = configs.get("Output_file_path", "").asString();
    config_var.Type_of_running = configs.get("Type_of_running", "GeometryCheck").asString();
    if (!configs["Number_of_events"].isNull() &&
        !isPositiveInt(configs["Number_of_events"]))
        throw std::runtime_error("Number_of_events must be a positive integer");
    config_var.Number_of_events = configs.get("Number_of_events", 1).asInt();
    config_var.Macro_file_path = configs.get("Macro_file_path", "").asString();
    config_var.Save_truth_particle_graph = configs.get("Save_truth_particle_graph", false).asBool();
    config_var.Use_high_granularity = configs.get("Use_high_granularity", false).asBool();
    config_var.Skip_unuseable_tracks = configs.get("Skip_unuseable_tracks", false).asBool();
    config_var.doSuperclustering = configs.get( "Do_superclustering", false ).asBool();

    config_var.r_inn_calo = configs["Geometry_definition"]["Inner_calorimeter_layer"].asDouble();
    config_var.Layer_gap = configs["Geometry_definition"]["Layer_gap"].asDouble();
    config_var.fieldValue = configs["Geometry_definition"]["Magnetic_field"].asDouble() * tesla;
    config_var.max_eta_barrel = configs["Geometry_definition"]["Max_eta_of_barrel_region"].asDouble();
    config_var.max_eta_endcap = configs["Geometry_definition"]["Max_eta_of_endcap_region"].asDouble();
    config_var.max_phi = 2. * M_PI * rad;
    config_var.check_geometry_overlap = configs["Geometry_definition"]["Check_Geometry_overlap"].asBool();
    config_var.check_geometry_overlap_only = configs["Geometry_definition"].get( "Check_Geometry_overlap_only", false ).asBool();
    if ( config_var.check_geometry_overlap_only )
	config_var.check_geometry_overlap = true;
    config_var.use_inner_detector = configs["Geometry_definition"].get( "Use_inner_detector", true ).asBool();
    config_var.use_ID_support = configs["Geometry_definition"].get( "Use_ID_support", true ).asBool();

    if (config_var.r_inn_calo <= 0. || config_var.Layer_gap < 0. ||
        config_var.max_eta_barrel <= 0. ||
        config_var.max_eta_endcap < config_var.max_eta_barrel)
        throw std::runtime_error("Invalid calorimeter dimensions or eta coverage");

    const Json::Value &geometry_config = configs["Geometry_definition"];
    const Json::Value &granularity = geometry_config["Detector_granularity"];
    const Json::Value &pixels_ecal = granularity["Number_of_pixels_ECAL"];
    const Json::Value &pixels_hcal = granularity["Number_of_pixels_HCAL"];
    const Json::Value &width_ecal = granularity["Width_of_ECAL_layers_in_X0"];
    const Json::Value &width_hcal = granularity["Width_of_HCAL_layers_in_Lambda_int"];
    const Json::Value &noise_ecal = geometry_config["Noise_in_ECAL"];
    const Json::Value &noise_hcal = geometry_config["Noise_in_HCAL"];
    requirePositiveIntegerArray(pixels_ecal, "Number_of_pixels_ECAL");
    requirePositiveIntegerArray(pixels_hcal, "Number_of_pixels_HCAL");
    requireBoundedNumericArray(width_ecal, "Width_of_ECAL_layers_in_X0", 0., false);
    requireBoundedNumericArray(width_hcal, "Width_of_HCAL_layers_in_Lambda_int", 0., false);
    requireBoundedNumericArray(noise_ecal, "Noise_in_ECAL", 0., true);
    requireBoundedNumericArray(noise_hcal, "Noise_in_HCAL", 0., true);
    requireSameSize(pixels_ecal, width_ecal, "ECAL pixels and widths");
    requireSameSize(pixels_ecal, noise_ecal, "ECAL pixels and noise");
    requireSameSize(pixels_hcal, width_hcal, "HCAL pixels and widths");
    requireSameSize(pixels_hcal, noise_hcal, "HCAL pixels and noise");

    const double sampling_ecal = geometry_config.get("SamplingFraction_ECAL", 0.02).asDouble();
    const double sampling_hcal = geometry_config.get("SamplingFraction_HCAL", 0.02).asDouble();
    if (sampling_ecal <= 0. || sampling_ecal > 1. ||
        sampling_hcal <= 0. || sampling_hcal > 1.)
        throw std::runtime_error("Sampling fractions must be in the interval (0, 1]");

    requireSameSize(geometry_config["Material_for_ECAL"],
                    geometry_config["ECAL_material_mixing_in_volume_proportion"],
                    "ECAL materials and mixing proportions");
    requireSameSize(geometry_config["Material_for_HCAL"],
                    geometry_config["HCAL_material_mixing_in_volume_proportion"],
                    "HCAL materials and mixing proportions");
    requireNonEmptyStringArray(geometry_config["Material_for_ECAL"],
                               "Material_for_ECAL");
    requireNonEmptyStringArray(geometry_config["Material_for_HCAL"],
                               "Material_for_HCAL");
    requireBoundedNumericArray(geometry_config["ECAL_material_mixing_in_volume_proportion"],
                               "ECAL material proportions", 0., false);
    requireBoundedNumericArray(geometry_config["HCAL_material_mixing_in_volume_proportion"],
                               "HCAL material proportions", 0., false);

    if (config_var.Use_high_granularity)
    {
        const Json::Value &high = geometry_config["High_Granularity_detector"];
        if (!high.isObject())
            throw std::runtime_error("Use_high_granularity requires High_Granularity_detector");
        requireMatchingNumeric2D(high["Number_of_pixels_ECAL"],
                                 high["Width_of_ECAL_layers_in_X0"],
                                 "high-granularity ECAL pixels and widths");
        requireMatchingNumeric2D(high["Number_of_pixels_HCAL"],
                                 high["Width_of_HCAL_layers_in_Lambda_int"],
                                 "high-granularity HCAL pixels and widths");
        for (Json::ArrayIndex i = 0; i < high["Number_of_pixels_ECAL"].size(); ++i)
        {
            requirePositiveIntegerArray(high["Number_of_pixels_ECAL"][i],
                                        "high-granularity ECAL pixels");
            requireBoundedNumericArray(high["Width_of_ECAL_layers_in_X0"][i],
                                       "high-granularity ECAL widths", 0., false);
        }
        for (Json::ArrayIndex i = 0; i < high["Number_of_pixels_HCAL"].size(); ++i)
        {
            requirePositiveIntegerArray(high["Number_of_pixels_HCAL"][i],
                                        "high-granularity HCAL pixels");
            requireBoundedNumericArray(high["Width_of_HCAL_layers_in_Lambda_int"][i],
                                       "high-granularity HCAL widths", 0., false);
        }
        requireSameSize(high["Number_of_pixels_ECAL"], pixels_ecal,
                        "low- and high-granularity ECAL layers");
        requireSameSize(high["Number_of_pixels_HCAL"], pixels_hcal,
                        "low- and high-granularity HCAL layers");
    }

    if (config_var.Type_of_running == "Standard")
    {
		if (config_var.Number_of_events <= 0)
			throw std::runtime_error("Number_of_events must be positive");
        const Json::Value &jets = configs["Jet_parameters"];
        const std::set<std::string> algorithms = {
            "kt_algorithm", "cambridge_algorithm", "antikt_algorithm",
            "genkt_algorithm", "ee_kt_algorithm", "ee_genkt_algorithm"};
        const std::set<std::string> schemes = {
            "E_scheme", "pt_scheme", "pt2_scheme", "Et_scheme", "Et2_scheme",
            "BIpt_scheme", "BIpt2_scheme", "WTA_pt_scheme",
            "WTA_modp_scheme"};
        if (!jets.isObject() || algorithms.count(jets["algorithm"].asString()) == 0 ||
            schemes.count(jets["recombination_scheme"].asString()) == 0 ||
            !jets["radius"].isNumeric() || jets["radius"].asDouble() <= 0. ||
            !jets["ptmin"].isNumeric() || jets["ptmin"].asDouble() < 0.)
            throw std::runtime_error("Jet_parameters is missing or invalid");
        const std::string algorithm = jets["algorithm"].asString();
        if ((algorithm == "genkt_algorithm" || algorithm == "ee_genkt_algorithm") &&
            !jets["power"].isNumeric())
            throw std::runtime_error("Generalised-kt algorithms require a numeric Jet_parameters.power");
        if (!configs["Fiducial_cuts"].isObject())
            throw std::runtime_error("Standard runs require Fiducial_cuts");
    }

    Json::Value &layervals = configs["Graph_construction"]["max_samelayer_edges"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_samelayer_edges);
                 layervals = configs["Graph_construction"]["max_dr"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_dr);
                 layervals = configs["Graph_construction"]["max_interlayer_edges"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_interlayer_edges);
    config_var.graph_construction.max_layer_sep = configs["Graph_construction"]["max_layer_sep"].asInt();
                 layervals = configs["Graph_construction"]["max_celltrack_edges"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_celltrack_edges);
                 layervals = configs["Graph_construction"]["max_celltrack_dr"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_celltrack_dr);

    config_var.topological_clustering.sigma_threshold_for_seed_cells = configs["TopoClustering"]["sigma_threshold_for_seed_cells"].asFloat();
    config_var.topological_clustering.sigma_threshold_for_neighboring_cells = configs["TopoClustering"]["sigma_threshold_for_neighboring_cells"].asFloat();
    config_var.topological_clustering.sigma_threshold_for_last_cells = configs["TopoClustering"]["sigma_threshold_for_last_cells"].asFloat();
    config_var.topological_clustering.cluster_negative_energy_cells = configs["TopoClustering"]["cluster_negative_energy_cells"].asBool();
    config_var.topological_clustering.local_max_seed_energy = configs["TopoClustering"]["local_max_seed_energy"].asFloat();
    
    config_var.particle_flow.S_discriminant_threshold = configs["Particle_flow"].get( "S_discriminant_threshold", -1.0 ).asFloat();
    config_var.particle_flow.E_div_p_threshold = configs["Particle_flow"].get( "E_div_p_threshold", 0.1 ).asFloat();
    config_var.particle_flow.factor_sigma_E_div_p_template = configs["Particle_flow"].get( "factor_sigma_E_div_p_template", 1.25 ).asFloat();
    config_var.particle_flow.Moliere_radius = configs["Particle_flow"].get( "Moliere_radius", 0.035 ).asFloat();
    
    config_var.jet_parameter.algorithm = configs["Jet_parameters"]["algorithm"].asString();
    config_var.jet_parameter.recombination_scheme = configs["Jet_parameters"]["recombination_scheme"].asString();
    config_var.jet_parameter.ptmin = configs["Jet_parameters"]["ptmin"].asDouble();
    config_var.jet_parameter.radius = configs["Jet_parameters"]["radius"].asDouble();
    config_var.jet_parameter.power = configs["Jet_parameters"].get("power", 0.).asDouble();

    config_var.fiducial_cuts.pt_min = configs["Fiducial_cuts"].get( "pt_min_gev", 1.0 ).asFloat();
    config_var.fiducial_cuts.eta_max = configs["Fiducial_cuts"].get( "eta_max", 3.0 ).asFloat();
    config_var.fiducial_cuts.dR_cut = configs["Fiducial_cuts"].get( "dR_cut", -1.0 ).asFloat();

    if (config_var.fiducial_cuts.pt_min < 0. ||
        config_var.fiducial_cuts.eta_max <= 0.)
        throw std::runtime_error("Fiducial cuts require pt_min_gev >= 0 and eta_max > 0");

    if (config_var.Type_of_running != "Standard")
    {
        std::cout<<config_var.Type_of_running<<std::endl;
        config_var.Use_high_granularity = false;
    }

    // ==============================================================
    // *Materials
    // ==============================================================
    config_var.Material_ECAL = Material_build("ECAL");
    config_var.Material_HCAL = Material_build("HCAL");

    //* Fill low resolution
    Json::Value &characters = configs["Geometry_definition"]["Detector_granularity"]["Number_of_pixels_ECAL"];
    Fill_1D_vector(characters, config_var.low_resolution.number_of_pixels_ECAL); // Low_number_of_pixels_ECAL
    characters = configs["Geometry_definition"]["Detector_granularity"]["Number_of_pixels_HCAL"];
    Fill_1D_vector(characters, config_var.low_resolution.number_of_pixels_HCAL); //Low_resolution_number_of_pixels_HCAL
    characters = configs["Geometry_definition"]["Detector_granularity"]["Width_of_ECAL_layers_in_X0"];
    Fill_1D_vector(characters, config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0); //Low_resolution_width_of_ECAL_layers_in_X0
    characters = configs["Geometry_definition"]["Detector_granularity"]["Width_of_HCAL_layers_in_Lambda_int"];
    config_var.samplingFraction_ECAL = configs["Geometry_definition"].get("SamplingFraction_ECAL", 0.02 ).asFloat();
    config_var.samplingFraction_HCAL = configs["Geometry_definition"].get("SamplingFraction_HCAL", 0.02 ).asFloat();
    Fill_1D_vector(characters, config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int); //Low_resolution_width_of_HCAL_layers_in_Lambda_int
    characters = configs["Geometry_definition"]["Noise_in_ECAL"];
    Fill_1D_vector(characters, config_var.low_resolution.layer_noise_ECAL); //Low_layer_noise_ECAL
    characters = configs["Geometry_definition"]["Noise_in_HCAL"];
    Fill_1D_vector(characters, config_var.low_resolution.layer_noise_HCAL); //Low_layer_noise_HCAL
    
    if ( configs["Particle_flow"].isObject() ) {
	config_var.doPFlow = true;
	characters = configs["Particle_flow"]["delta_Rprime_threshold"];
	Fill_1D_vector(characters, config_var.particle_flow.delta_Rprime_threshold);
	characters = configs["Particle_flow"]["momentum_delta_Rprime_threshold"];
	Fill_1D_vector(characters, config_var.particle_flow.momentum_delta_Rprime_threshold);
    } else {
	config_var.doPFlow = false;
	config_var.particle_flow.delta_Rprime_threshold = { 2.4, 1.25, 0.8 };
	config_var.particle_flow.momentum_delta_Rprime_threshold = { 2000, 5000 };
    }
    
    config_var.low_resolution.layer_inn_radius_ECAL = config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0;
    config_var.low_resolution.layer_mid_radius_ECAL = config_var.low_resolution.layer_inn_radius_ECAL;
    config_var.low_resolution.layer_out_radius_ECAL = config_var.low_resolution.layer_mid_radius_ECAL;
    config_var.low_resolution.layer_deta_ECAL = config_var.low_resolution.layer_out_radius_ECAL;
    config_var.low_resolution.layer_dphi_ECAL = config_var.low_resolution.layer_deta_ECAL;
    config_var.low_resolution.layer_inn_radius_HCAL = config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int;
    config_var.low_resolution.layer_mid_radius_HCAL = config_var.low_resolution.layer_inn_radius_HCAL;
    config_var.low_resolution.layer_out_radius_HCAL = config_var.low_resolution.layer_mid_radius_HCAL;
    config_var.low_resolution.layer_deta_HCAL = config_var.low_resolution.layer_out_radius_HCAL;
    config_var.low_resolution.layer_dphi_HCAL = config_var.low_resolution.layer_deta_HCAL;
    
    int kNLayers = 0;
    int nLow_Layers = config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0.size();
    kNLayers += nLow_Layers;
    long double r_inn = config_var.r_inn_calo;
    for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
    {
        config_var.low_resolution.layer_noise.push_back(config_var.low_resolution.layer_noise_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_ECAL.at(ilow_layer).at(0) = 2 * config_var.max_eta_endcap / config_var.low_resolution.number_of_pixels_ECAL.at(ilow_layer).at(0); // Low_number_of_pixels_ECAL.at(ilow_layer);
        config_var.low_resolution.layer_dphi_ECAL.at(ilow_layer).at(0) = config_var.max_phi / config_var.low_resolution.number_of_pixels_ECAL.at(ilow_layer).at(0);
        config_var.low_resolution.number_of_pixels_flatten.push_back(config_var.low_resolution.number_of_pixels_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_flatten.push_back(config_var.low_resolution.layer_deta_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_dphi_flatten.push_back(config_var.low_resolution.layer_dphi_ECAL.at(ilow_layer).at(0));
        long double r_out = r_inn + config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0.at(ilow_layer).at(0) * config_var.Material_ECAL->GetRadlen();
        config_var.low_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(0) = r_inn;
        config_var.low_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(0) = r_out - 0.5 * (r_out - r_inn);
        config_var.low_resolution.layer_out_radius_ECAL.at(ilow_layer).at(0) = r_out;
        config_var.low_resolution.layer_inn_radius_flatten.push_back(config_var.low_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_mid_radius_flatten.push_back(config_var.low_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_out_radius_flatten.push_back(config_var.low_resolution.layer_out_radius_ECAL.at(ilow_layer).at(0));
        r_inn = r_out;
    }
    r_inn = r_inn + config_var.Layer_gap;
    nLow_Layers = config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int.size();
    kNLayers += nLow_Layers;
    for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
    {
        config_var.low_resolution.layer_noise.push_back(config_var.low_resolution.layer_noise_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_HCAL.at(ilow_layer).at(0) = 2 * config_var.max_eta_endcap / config_var.low_resolution.number_of_pixels_HCAL.at(ilow_layer).at(0);
        config_var.low_resolution.layer_dphi_HCAL.at(ilow_layer).at(0) = config_var.max_phi / config_var.low_resolution.number_of_pixels_HCAL.at(ilow_layer).at(0);
        config_var.low_resolution.number_of_pixels_flatten.push_back(config_var.low_resolution.number_of_pixels_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_dphi_flatten.push_back(config_var.low_resolution.layer_dphi_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_flatten.push_back(config_var.low_resolution.layer_deta_HCAL.at(ilow_layer).at(0));
        long double r_out = r_inn + config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int.at(ilow_layer).at(0) * config_var.Material_HCAL->GetNuclearInterLength();
        config_var.low_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(0) = r_inn;
        config_var.low_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(0) = r_out - 0.5 * (r_out - r_inn);
        config_var.low_resolution.layer_out_radius_HCAL.at(ilow_layer).at(0) = r_out;
        config_var.low_resolution.layer_inn_radius_flatten.push_back(config_var.low_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_mid_radius_flatten.push_back(config_var.low_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_out_radius_flatten.push_back(config_var.low_resolution.layer_out_radius_HCAL.at(ilow_layer).at(0));
        r_inn = r_out;
    }
    //* Fill low resolution end
    config_var.low_resolution.kNLayers = kNLayers;
    if (config_var.Use_high_granularity)
    {
        kNLayers = 0;
        //* Fill high resolution
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Number_of_pixels_ECAL"];
        Fill_2D_vector(characters, config_var.high_resolution.number_of_pixels_ECAL);
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Number_of_pixels_HCAL"];
        Fill_2D_vector(characters, config_var.high_resolution.number_of_pixels_HCAL);
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Width_of_ECAL_layers_in_X0"];
        Fill_2D_vector(characters, config_var.high_resolution.resolution_width_of_ECAL_layers_in_X0);
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Width_of_HCAL_layers_in_Lambda_int"];
        Fill_2D_vector(characters, config_var.high_resolution.resolution_width_of_HCAL_layers_in_Lambda_int);
        //* Fill high resolution end
        nLow_Layers = config_var.high_resolution.number_of_pixels_ECAL.size();
        config_var.high_resolution.layer_inn_radius_ECAL = config_var.high_resolution.resolution_width_of_ECAL_layers_in_X0;
        config_var.high_resolution.layer_mid_radius_ECAL = config_var.high_resolution.layer_inn_radius_ECAL;
        config_var.high_resolution.layer_out_radius_ECAL = config_var.high_resolution.layer_mid_radius_ECAL;
        config_var.high_resolution.layer_deta_ECAL = config_var.high_resolution.layer_out_radius_ECAL;
        config_var.high_resolution.layer_dphi_ECAL = config_var.high_resolution.layer_deta_ECAL;
        config_var.high_resolution.layer_inn_radius_HCAL = config_var.high_resolution.resolution_width_of_HCAL_layers_in_Lambda_int;
        config_var.high_resolution.layer_mid_radius_HCAL = config_var.high_resolution.layer_inn_radius_HCAL;
        config_var.high_resolution.layer_out_radius_HCAL = config_var.high_resolution.layer_mid_radius_HCAL;
        config_var.high_resolution.layer_deta_HCAL = config_var.high_resolution.layer_out_radius_HCAL;
        config_var.high_resolution.layer_dphi_HCAL = config_var.high_resolution.layer_deta_HCAL;
        r_inn = config_var.r_inn_calo;
        for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
        {

            int nHigh_Layers = config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).size();
            kNLayers += nHigh_Layers;
            for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
            {
                config_var.high_resolution.layer_deta_ECAL.at(ilow_layer).at(ihigh_layer) = 2 * config_var.max_eta_endcap / config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).at(ihigh_layer);
                config_var.high_resolution.layer_dphi_ECAL.at(ilow_layer).at(ihigh_layer) = config_var.max_phi / config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).at(ihigh_layer);
                config_var.high_resolution.number_of_pixels_flatten.push_back(config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_dphi_flatten.push_back(config_var.high_resolution.layer_dphi_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_deta_flatten.push_back(config_var.high_resolution.layer_deta_ECAL.at(ilow_layer).at(ihigh_layer));
                long double r_out = r_inn + config_var.high_resolution.resolution_width_of_ECAL_layers_in_X0.at(ilow_layer).at(ihigh_layer) * config_var.Material_ECAL->GetRadlen();
                config_var.high_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(ihigh_layer) = r_inn;
                config_var.high_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(ihigh_layer) = r_out - 0.5 * (r_out - r_inn);
                config_var.high_resolution.layer_out_radius_ECAL.at(ilow_layer).at(ihigh_layer) = r_out;
                config_var.high_resolution.layer_inn_radius_flatten.push_back(config_var.high_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_mid_radius_flatten.push_back(config_var.high_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_out_radius_flatten.push_back(config_var.high_resolution.layer_out_radius_ECAL.at(ilow_layer).at(ihigh_layer));
                r_inn = r_out;
            }
        }
        r_inn = r_inn + config_var.Layer_gap;
        nLow_Layers = config_var.high_resolution.number_of_pixels_HCAL.size();
        for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
        {
            int nHigh_Layers = config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).size();
            kNLayers += nHigh_Layers;
            for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
            {
                config_var.high_resolution.layer_deta_HCAL.at(ilow_layer).at(ihigh_layer) = 2 * config_var.max_eta_endcap / config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).at(ihigh_layer);
                config_var.high_resolution.layer_dphi_HCAL.at(ilow_layer).at(ihigh_layer) = config_var.max_phi / config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).at(ihigh_layer);
                long double r_out = r_inn + config_var.high_resolution.resolution_width_of_HCAL_layers_in_Lambda_int.at(ilow_layer).at(ihigh_layer) * config_var.Material_HCAL->GetNuclearInterLength(); // Material_H->GetRadlen();
                config_var.high_resolution.number_of_pixels_flatten.push_back(config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_dphi_flatten.push_back(config_var.high_resolution.layer_dphi_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_deta_flatten.push_back(config_var.high_resolution.layer_deta_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(ihigh_layer) = r_inn;
                config_var.high_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(ihigh_layer) = r_out - 0.5 * (r_out - r_inn);
                config_var.high_resolution.layer_out_radius_HCAL.at(ilow_layer).at(ihigh_layer) = r_out;
                config_var.high_resolution.layer_inn_radius_flatten.push_back(config_var.high_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_mid_radius_flatten.push_back(config_var.high_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_out_radius_flatten.push_back(config_var.high_resolution.layer_out_radius_HCAL.at(ilow_layer).at(ihigh_layer));
                r_inn = r_out;
            }
        }
    }
    
    config_var.high_resolution.kNLayers = kNLayers;

    config_var.run_hadron_test = configs.get( "run_hadron_test", false ).asBool();
    config_var.run_piZero_test = configs.get( "run_piZero_test", false ).asBool();
    config_var.run_jets_test   = configs.get( "run_jets_test", false ).asBool();
    size_t n_tests = 0;
    for ( bool test : { config_var.run_hadron_test,
		        config_var.run_piZero_test,
		config_var.run_jets_test } ) {
	if ( test)
	    ++n_tests;
    }
    if ( n_tests > 1 )
	G4cerr << "Warning: chose more than one test to perform, but tests are meant to be mutually exclusive." << G4endl;
    
}

G4Material *Config_reader_func::Material_build(std::string name)
{
    G4NistManager *nistManager = G4NistManager::Instance();
    long double a, z, density_liq, fracMass;
    std::vector<std::string> sMaterial_vector;
    std::vector<G4Material *> Material_vector;
    G4Material *Material;

    Json::Value &characters = configs["Geometry_definition"]["Material_for_" + name];
    Fill_1D_vector(characters, sMaterial_vector);
    int icustom_element = 0;
    for (unsigned long imat_ecal = 0; imat_ecal < sMaterial_vector.size(); imat_ecal++)
    {
        G4Material *element = nistManager->FindOrBuildMaterial(sMaterial_vector[imat_ecal]);
        if (!element)
        {
            const Json::Value &custom =
                configs["Geometry_definition"]["Characteristic_of_custom_material_for_" + name];
            if (!custom.isArray() || icustom_element >= static_cast<int>(custom.size()) ||
                !custom[icustom_element].isArray() || custom[icustom_element].size() != 3 ||
                !custom[icustom_element][0].isNumeric() ||
                !custom[icustom_element][1].isNumeric() ||
                !custom[icustom_element][2].isNumeric() ||
                custom[icustom_element][0].asDouble() <= 0. ||
                custom[icustom_element][1].asDouble() <= 0. ||
                custom[icustom_element][2].asDouble() <= 0.)
                throw std::runtime_error("Custom material " + sMaterial_vector[imat_ecal] +
                                         " requires positive [Z, A, density] values");
            element = new G4Material(sMaterial_vector[imat_ecal], 
                                     z = custom[icustom_element][0].asDouble(),
                                     a = custom[icustom_element][1].asDouble() * g / mole,
                                     density_liq = custom[icustom_element][2].asDouble() * g / cm3);
            icustom_element++;
        }
        Material_vector.push_back(element);
    }
    long double numerator = 0;
    long double denominator = 0;
    for (unsigned long irho_ecal = 0; irho_ecal < configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"].size(); irho_ecal++)
    {
        numerator += Material_vector[irho_ecal]->GetDensity() * configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"][(int)irho_ecal].asDouble();
        denominator += configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"][(int)irho_ecal].asDouble();
    }
    Material = new G4Material("Material_" + name, numerator / denominator, sMaterial_vector.size());
    for (unsigned long imat_ecal = 0; imat_ecal < Material_vector.size(); imat_ecal++)
    {
        Material->AddMaterial(Material_vector[imat_ecal],
                              fracMass = (Material_vector[imat_ecal]->GetDensity() * configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"][(int)imat_ecal].asDouble()) / (numerator)*100 * perCent);
    }
    return Material;
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<std::vector<long double>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<long double> vec_elem;
        vec_elem.push_back(Json_list[i].asDouble());
        vec.push_back(vec_elem);
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<float> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        vec.push_back(Json_list[i].asDouble());
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<int> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        vec.push_back(Json_list[i].asInt());
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<std::vector<int>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<int> vec_elem;
        vec_elem.push_back(Json_list[i].asUInt());
        vec.push_back(vec_elem);
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<std::string> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        vec.push_back(Json_list[i].asString());
    }
}
void Config_reader_func::Fill_2D_vector(const Json::Value &Json_list, std::vector<std::vector<int>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<int> vec_elem;
        vec_elem.clear();
        for (unsigned int j = 0; j < Json_list[i].size(); j++)
        {
            vec_elem.push_back(Json_list[i][j].asUInt());
        }
        vec.push_back(vec_elem);
    }
}
void Config_reader_func::Fill_2D_vector(const Json::Value &Json_list, std::vector<std::vector<long double>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<long double> vec_elem;
        vec_elem.clear();
        for (unsigned int j = 0; j < Json_list[i].size(); j++)
        {
            vec_elem.push_back(Json_list[i][j].asDouble());
        }
        vec.push_back(vec_elem);
    }
}
