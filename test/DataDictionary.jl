# ----------------------------------------------------------------------------------- #
# Copyright (c) 2018 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2018-01-06T13:57:00.029
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	# What is the system dimension? -
	(number_of_species,number_of_reactions) = size(stoichiometric_matrix)

	# Setup default flux bounds array -
	default_bounds_array = [
		0	187.20000000000002	;	# 1 A --> B
		0	187.20000000000002	;	# 2 B --> C
		0	187.20000000000002	;	# 3 A_e --> A
		0	187.20000000000002	;	# 4 C_e --> C
		0	187.20000000000002	;	# 5 C --> C_e
		0	187.20000000000002	;	# 6 G1+RNAP --> G1_RNAP
		0	187.20000000000002	;	# 7 G1_RNAP --> G1+RNAP+mRNA_E1
		0	187.20000000000002	;	# 8 G2+RNAP --> G2_RNAP
		0	187.20000000000002	;	# 9 G2_RNAP --> G2+RNAP+mRNA_E2
		0	187.20000000000002	;	# 10 G3+RNAP --> G3_RNAP
		0	187.20000000000002	;	# 11 G3_RNAP --> G3+RNAP+mRNA_E3
		0	187.20000000000002	;	# 12 G4+RNAP --> G4_RNAP
		0	187.20000000000002	;	# 13 G4_RNAP --> G4+RNAP+mRNA_E4
		0	187.20000000000002	;	# 14 mRNA_E1+RIBOSOME --> mRNA_E1_RIBOSOME
		0	187.20000000000002	;	# 15 mRNA_E1_RIBOSOME --> E1+mRNA_E1+RIBOSOME
		0	187.20000000000002	;	# 16 mRNA_E2+RIBOSOME --> mRNA_E2_RIBOSOME
		0	187.20000000000002	;	# 17 mRNA_E2_RIBOSOME --> E2+mRNA_E2+RIBOSOME
		0	187.20000000000002	;	# 18 mRNA_E3+RIBOSOME --> mRNA_E3_RIBOSOME
		0	187.20000000000002	;	# 19 mRNA_E3_RIBOSOME --> E3+mRNA_E3+RIBOSOME
		0	187.20000000000002	;	# 20 mRNA_E4+RIBOSOME --> mRNA_E4_RIBOSOME
		0	187.20000000000002	;	# 21 mRNA_E4_RIBOSOME --> E4+mRNA_E4+RIBOSOME
		0	187.20000000000002	;	# 22 E1 --> []
		0	187.20000000000002	;	# 23 E2 --> []
		0	187.20000000000002	;	# 24 E3 --> []
		0	187.20000000000002	;	# 25 E4 --> []
		0	187.20000000000002	;	# 26 mRNA_E1 --> []
		0	187.20000000000002	;	# 27 mRNA_E2 --> []
		0	187.20000000000002	;	# 28 mRNA_E3 --> []
		0	187.20000000000002	;	# 29 mRNA_E4 --> []
		0	187.20000000000002	;	# 30 RIBOSOME --> []
		0	187.20000000000002	;	# 31 [] --> RIBOSOME
		0	187.20000000000002	;	# 32 RNAP --> []
		0	187.20000000000002	;	# 33 [] --> RNAP
		0	187.20000000000002	;	# 34 [] --> G1
		0	187.20000000000002	;	# 35 [] --> G2
		0	187.20000000000002	;	# 36 [] --> G3
		0	187.20000000000002	;	# 37 [] --> G4
		0	187.20000000000002	;	# 38 X_e --> []
		0	187.20000000000002	;	# 39 [] --> X_e
	];


	# Setup species abundance array -
	species_abundance_array = [
		0.0	;	# 1 A
		0.0	;	# 2 B
		0.0	;	# 3 C
		0.0	;	# 4 E1
		0.0	;	# 5 E2
		0.0	;	# 6 E3
		0.0	;	# 7 E4
		0.0	;	# 8 G1
		0.0	;	# 9 G1_RNAP
		0.0	;	# 10 G2
		0.0	;	# 11 G2_RNAP
		0.0	;	# 12 G3
		0.0	;	# 13 G3_RNAP
		0.0	;	# 14 G4
		0.0	;	# 15 G4_RNAP
		0.0	;	# 16 RIBOSOME
		0.0	;	# 17 RNAP
		0.0	;	# 18 mRNA_E1
		0.0	;	# 19 mRNA_E1_RIBOSOME
		0.0	;	# 20 mRNA_E2
		0.0	;	# 21 mRNA_E2_RIBOSOME
		0.0	;	# 22 mRNA_E3
		0.0	;	# 23 mRNA_E3_RIBOSOME
		0.0	;	# 24 mRNA_E4
		0.0	;	# 25 mRNA_E4_RIBOSOME
		0.0	;	# 26 A_e
		0.0	;	# 27 C_e
		0.0	;	# 28 X_e
	];

	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1 A
		0.0	0.0	;	# 2 B
		0.0	0.0	;	# 3 C
		0.0	0.0	;	# 4 E1
		0.0	0.0	;	# 5 E2
		0.0	0.0	;	# 6 E3
		0.0	0.0	;	# 7 E4
		0.0	0.0	;	# 8 G1
		0.0	0.0	;	# 9 G1_RNAP
		0.0	0.0	;	# 10 G2
		0.0	0.0	;	# 11 G2_RNAP
		0.0	0.0	;	# 12 G3
		0.0	0.0	;	# 13 G3_RNAP
		0.0	0.0	;	# 14 G4
		0.0	0.0	;	# 15 G4_RNAP
		0.0	0.0	;	# 16 RIBOSOME
		0.0	0.0	;	# 17 RNAP
		0.0	0.0	;	# 18 mRNA_E1
		0.0	0.0	;	# 19 mRNA_E1_RIBOSOME
		0.0	0.0	;	# 20 mRNA_E2
		0.0	0.0	;	# 21 mRNA_E2_RIBOSOME
		0.0	0.0	;	# 22 mRNA_E3
		0.0	0.0	;	# 23 mRNA_E3_RIBOSOME
		0.0	0.0	;	# 24 mRNA_E4
		0.0	0.0	;	# 25 mRNA_E4_RIBOSOME
		-1.0	1.0	;	# 26 A_e
		-1.0	1.0	;	# 27 C_e
		-1.0	1.0	;	# 28 X_e
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# Setup the objective coefficient array -
	objective_coefficient_array = [
		0.0	;	# 1 R1_E1::A --> B
		0.0	;	# 2 R2_E2::B --> C
		0.0	;	# 3 Ex1_E3::A_e --> A
		0.0	;	# 4 Ex2_E4::C_e --> C
		0.0	;	# 5 Ex2_E4_r::C --> C_e
		0.0	;	# 6 T1_open::G1+RNAP --> G1_RNAP
		0.0	;	# 7 T1::G1_RNAP --> G1+RNAP+mRNA_E1
		0.0	;	# 8 T2_open::G2+RNAP --> G2_RNAP
		0.0	;	# 9 T2::G2_RNAP --> G2+RNAP+mRNA_E2
		0.0	;	# 10 T3_open::G3+RNAP --> G3_RNAP
		0.0	;	# 11 T3::G3_RNAP --> G3+RNAP+mRNA_E3
		0.0	;	# 12 T4_open::G4+RNAP --> G4_RNAP
		0.0	;	# 13 T4::G4_RNAP --> G4+RNAP+mRNA_E4
		0.0	;	# 14 X1_open::mRNA_E1+RIBOSOME --> mRNA_E1_RIBOSOME
		0.0	;	# 15 X1::mRNA_E1_RIBOSOME --> E1+mRNA_E1+RIBOSOME
		0.0	;	# 16 X2_open::mRNA_E2+RIBOSOME --> mRNA_E2_RIBOSOME
		0.0	;	# 17 X2::mRNA_E2_RIBOSOME --> E2+mRNA_E2+RIBOSOME
		0.0	;	# 18 X3_open::mRNA_E3+RIBOSOME --> mRNA_E3_RIBOSOME
		0.0	;	# 19 X3::mRNA_E3_RIBOSOME --> E3+mRNA_E3+RIBOSOME
		0.0	;	# 20 X4_open::mRNA_E4+RIBOSOME --> mRNA_E4_RIBOSOME
		0.0	;	# 21 X4::mRNA_E4_RIBOSOME --> E4+mRNA_E4+RIBOSOME
		0.0	;	# 22 E1_degradation::E1 --> []
		0.0	;	# 23 E2_degradation::E2 --> []
		0.0	;	# 24 E3_degradation::E3 --> []
		0.0	;	# 25 E4_degradation::E4 --> []
		0.0	;	# 26 mRNA_E1_degradation::mRNA_E1 --> []
		0.0	;	# 27 mRNA_E2_degradation::mRNA_E2 --> []
		0.0	;	# 28 mRNA_E3_degradation::mRNA_E3 --> []
		0.0	;	# 29 mRNA_E4_degradation::mRNA_E4 --> []
		0.0	;	# 30 RIBOSOME_SYNTHESIS_DEGRADATION::RIBOSOME --> []
		0.0	;	# 31 RIBOSOME_SYNTHESIS_DEGRADATION_reverse::[] --> RIBOSOME
		0.0	;	# 32 RNAP_SYNTHESIS_DEGRADATION::RNAP --> []
		0.0	;	# 33 RNAP_SYNTHESIS_DEGRADATION_reverse::[] --> RNAP
		0.0	;	# 34 G1_SYNTHESIS::[] --> G1
		0.0	;	# 35 G2_SYNTHESIS::[] --> G2
		0.0	;	# 36 G3_SYNTHESIS::[] --> G3
		0.0	;	# 37 G4_SYNTHESIS::[] --> G4
		0.0	;	# 38 CELLGROWTH::X_e --> []
		0.0	;	# 39 CELLGROWTH_reverse::[] --> X_e
	];

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"R1_E1::A --> B"
		"R2_E2::B --> C"
		"Ex1_E3::A_e --> A"
		"Ex2_E4::C_e --> C"
		"Ex2_E4_r::C --> C_e"
		"T1_open::G1+RNAP --> G1_RNAP"
		"T1::G1_RNAP --> G1+RNAP+mRNA_E1"
		"T2_open::G2+RNAP --> G2_RNAP"
		"T2::G2_RNAP --> G2+RNAP+mRNA_E2"
		"T3_open::G3+RNAP --> G3_RNAP"
		"T3::G3_RNAP --> G3+RNAP+mRNA_E3"
		"T4_open::G4+RNAP --> G4_RNAP"
		"T4::G4_RNAP --> G4+RNAP+mRNA_E4"
		"X1_open::mRNA_E1+RIBOSOME --> mRNA_E1_RIBOSOME"
		"X1::mRNA_E1_RIBOSOME --> E1+mRNA_E1+RIBOSOME"
		"X2_open::mRNA_E2+RIBOSOME --> mRNA_E2_RIBOSOME"
		"X2::mRNA_E2_RIBOSOME --> E2+mRNA_E2+RIBOSOME"
		"X3_open::mRNA_E3+RIBOSOME --> mRNA_E3_RIBOSOME"
		"X3::mRNA_E3_RIBOSOME --> E3+mRNA_E3+RIBOSOME"
		"X4_open::mRNA_E4+RIBOSOME --> mRNA_E4_RIBOSOME"
		"X4::mRNA_E4_RIBOSOME --> E4+mRNA_E4+RIBOSOME"
		"E1_degradation::E1 --> []"
		"E2_degradation::E2 --> []"
		"E3_degradation::E3 --> []"
		"E4_degradation::E4 --> []"
		"mRNA_E1_degradation::mRNA_E1 --> []"
		"mRNA_E2_degradation::mRNA_E2 --> []"
		"mRNA_E3_degradation::mRNA_E3 --> []"
		"mRNA_E4_degradation::mRNA_E4 --> []"
		"RIBOSOME_SYNTHESIS_DEGRADATION::RIBOSOME --> []"
		"RIBOSOME_SYNTHESIS_DEGRADATION_reverse::[] --> RIBOSOME"
		"RNAP_SYNTHESIS_DEGRADATION::RNAP --> []"
		"RNAP_SYNTHESIS_DEGRADATION_reverse::[] --> RNAP"
		"G1_SYNTHESIS::[] --> G1"
		"G2_SYNTHESIS::[] --> G2"
		"G3_SYNTHESIS::[] --> G3"
		"G4_SYNTHESIS::[] --> G4"
		"CELLGROWTH::X_e --> []"
		"CELLGROWTH_reverse::[] --> X_e"
	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"A"
		"B"
		"C"
		"E1"
		"E2"
		"E3"
		"E4"
		"G1"
		"G1_RNAP"
		"G2"
		"G2_RNAP"
		"G3"
		"G3_RNAP"
		"G4"
		"G4_RNAP"
		"RIBOSOME"
		"RNAP"
		"mRNA_E1"
		"mRNA_E1_RIBOSOME"
		"mRNA_E2"
		"mRNA_E2_RIBOSOME"
		"mRNA_E3"
		"mRNA_E3_RIBOSOME"
		"mRNA_E4"
		"mRNA_E4_RIBOSOME"
		"A_e"
		"C_e"
		"X_e"
	];


	# Metabolic kcat array (units:hr^-1) -
	metabolic_rate_constant_array = [
		46800.0	# 1	A --> B
		46800.0	# 2	B --> C
		46800.0	# 3	A_e --> A
		46800.0	# 4	C_e --> C
		46800.0	# 5	C --> C_e
	];

	# Metabolic saturation constant array (units mM) -
	number_of_metabolic_rates = length(metabolic_rate_constant_array)
	metabolic_saturation_constant_array = 0.130*ones(number_of_metabolic_rates*number_of_species)

	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                             # mum
	number_of_rnapII = 4600            	            # copies/cells
	number_of_ribosome = 50000         	            # copies/cells
	mRNA_half_life_TF = 0.083                       # hrs
	protein_half_life = 70                          # hrs
	doubling_time_cell = 0.33                       # hrs
	max_translation_rate = 16.5                     # aa/sec
	max_transcription_rate = 60.0                   # nt/sec
	transcription_initiation_time_contstant = 400  # sec
	average_transcript_length = 1200   	            # nt
	average_protein_length = 400       	            # aa
	fraction_nucleus = 0.0             	            # dimensionless
	av_number = 6.02e23                             # number/mol
	avg_gene_number = 2                             # number of copies of a gene
	polysome_number = 4					            # number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                         # nM
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/V)*1e9                     # nM

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(e,0.5)                           # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(e,0.5)                        # hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

	# kcat for transcription initiation -
	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
	kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(e,2)                          # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                        # nM

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/V)*1e9                                 # nM
	saturation_translation = 150000*(1/av_number)*(1/V)*1e9                                 # nM
	# -------------------------------------------------------------------------------------------#

	# Alias the txtl parameters -
	txtl_parameter_dictionary = Dict{AbstractString,Float64}()
	txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	txtl_parameter_dictionary["death_rate_constant"] = death_rate_constant
	txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
	txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
	txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation
	txtl_parameter_dictionary["average_transcript_length"] = average_transcript_length
	txtl_parameter_dictionary["average_protein_length"] = average_protein_length
	txtl_parameter_dictionary["kcat_transcription_initiation"] = kcat_transcription_initiation
	txtl_parameter_dictionary["kcat_translation_initiation"] = kcat_translation_initiation

	# txtl_parameter_dictionary["gene_coding_length_array"] = gene_coding_length_array
	# txtl_parameter_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	# txtl_parameter_dictionary["protein_coding_length_array"] = protein_coding_length_array


	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_abundance_array"] = species_abundance_array
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	data_dictionary["number_of_species"] = number_of_species
	data_dictionary["number_of_reactions"] = number_of_reactions
	data_dictionary["metabolic_saturation_constant_array"] = metabolic_saturation_constant_array
	data_dictionary["metabolic_rate_constant_array"] = metabolic_rate_constant_array
	data_dictionary["txtl_parameter_dictionary"] = txtl_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
