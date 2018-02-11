# ------------------------------------------------------------------------------------------#
# constants (from bionumbers)       units
# ------------------------------------------------------------------------------------------#
cell_diameter = 1.1                             # mum
mass_of_single_cell = 2.8e-13                   # g
number_of_rnapII = 4600            	            # copies/cells
number_of_ribosome = 50000         	            # copies/cells
mRNA_half_life_TF = 0.083                       # hrs
protein_half_life = 70                          # hrs
infrastructure_half_life = 300					# hrs
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
rnapII_concentration = number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*1e6       # mumol/gdw
ribosome_concentration = number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*1e6   # mumol/ddw

# degrdation rate constants -
degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(e,0.5)                           # hr^-1
degradation_constant_protein = -(1/protein_half_life)*log(e,0.5)                        # hr^-1
degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(e,0.5)			# hr^-1

# kcats for transcription and translation -
kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

# kcat for transcription initiation -
kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

# Maximum specific growth rate -
maximum_specific_growth_rate = (1/doubling_time_cell)*log(e,2)                          # hr^-1

# What is the average gene concentration -
avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/V)*1e9              # nmol/gdw

# How fast do my cells die?
death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

# Saturation constants for translation and trascription -
saturation_transcription = 4600*(1/av_number)*(1/mass_of_single_cell)*1e9               # nmol/gdw
saturation_translation = 150000*(1/av_number)*(1/mass_of_single_cell)*1e6               # nmol/gdw
# -------------------------------------------------------------------------------------------#
