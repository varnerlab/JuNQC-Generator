# Alias the txtl parameters -
txtl_parameter_dictionary = Dict{AbstractString,Any}()
txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
txtl_parameter_dictionary["degrdation_constant_infrastructure"] = degrdation_constant_infrastructure  # hr^-1
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
