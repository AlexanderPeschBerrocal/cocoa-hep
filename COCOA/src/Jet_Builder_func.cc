#include "Jet_Builder_func.hh"

void Jet_Builder_func::build_jets(std::vector<fastjet::PseudoJet> &input_jets,
				  Jet_Builder_data &jet_data,
				  Jet_parameters jet_par)
{
	const fastjet::JetAlgorithm jet_algorithm = algorithm(jet_par.algorithm);
	const fastjet::RecombinationScheme scheme =
		recombination_scheme(jet_par.recombination_scheme);
	fastjet::JetDefinition jet_def;
	if (jet_algorithm == fastjet::ee_kt_algorithm)
	{
		jet_def = fastjet::JetDefinition(jet_algorithm, scheme, fastjet::Best);
	}
	else if (jet_algorithm == fastjet::genkt_algorithm ||
			 jet_algorithm == fastjet::ee_genkt_algorithm)
	{
		jet_def = fastjet::JetDefinition(jet_algorithm, jet_par.radius,
			jet_par.power, scheme, fastjet::Best);
	}
	else
	{
		jet_def = fastjet::JetDefinition(jet_algorithm, jet_par.radius,
			scheme, fastjet::Best);
	}
	cs = new fastjet::ClusterSequence(input_jets, jet_def);
	jet_data.jets = sorted_by_pt(cs->inclusive_jets(jet_par.ptmin));
	jet_data.set_n_constituents( input_jets.size() );
}
void Jet_Builder_func::reset()
{
	delete cs;
}
fastjet::JetAlgorithm Jet_Builder_func::algorithm(std::string algo)
{
	fastjet::JetAlgorithm Algorithm = fastjet::kt_algorithm;
	if (algo == "kt_algorithm")
		Algorithm = fastjet::kt_algorithm;
	else if (algo == "cambridge_algorithm")
		Algorithm = fastjet::cambridge_algorithm;
	else if (algo == "antikt_algorithm")
		Algorithm = fastjet::antikt_algorithm;
	else if (algo == "genkt_algorithm")
		Algorithm = fastjet::genkt_algorithm;
	else if (algo == "ee_kt_algorithm")
		Algorithm = fastjet::ee_kt_algorithm;
	else if (algo == "ee_genkt_algorithm")
		Algorithm = fastjet::ee_genkt_algorithm;
	return Algorithm;
}
fastjet::RecombinationScheme Jet_Builder_func::recombination_scheme(std::string reco)
{
	fastjet::RecombinationScheme Scheme =  fastjet::E_scheme;
	if (reco == "E_scheme")
		Scheme = fastjet::E_scheme;
	else if (reco == "pt_scheme")
		Scheme = fastjet::pt_scheme;
	else if (reco == "pt2_scheme")
		Scheme = fastjet::pt2_scheme;
	else if (reco == "Et_scheme")
		Scheme = fastjet::Et_scheme;
	else if (reco == "Et2_scheme")
		Scheme = fastjet::Et2_scheme;
	else if (reco == "BIpt_scheme")
		Scheme = fastjet::BIpt_scheme;
	else if (reco == "BIpt2_scheme")
		Scheme = fastjet::BIpt2_scheme;
	else if (reco == "WTA_pt_scheme")
		Scheme = fastjet::WTA_pt_scheme;
	else if (reco == "WTA_modp_scheme")
		Scheme = fastjet::WTA_modp_scheme;
	return Scheme;
}
