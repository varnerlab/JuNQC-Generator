# Script to convert kcat data from:
# Adadi R, Volkmer B, Milo R, Heinemann M, Shlomi T (2012)
# Prediction of Microbial Growth Rate versus Biomass Yield by a Metabolic Network with Kinetic Parameters.
# PLoS Comput Biol 8(7): e1002575. doi:10.1371/journal.pcbi.1002575

# path to csv data file -
path_to_data_file = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/julia_work/simple_HL60_metabolic_model/EnzymeData.csv"
path_to_json_file = "../../config/Enzymes.json"

# read the data file -
data_array = readdlm(path_to_data_file,',')

# how many rows and cols do we have?
(number_of_rows,number_of_cols) = size(data_array)

@show size(data_array)

# Initialize the buffer -
buffer = ""
buffer *= "{\n"
buffer *= "\t\"enzyme_data\":[\n"

# go through each row -
for row_index = 1:number_of_rows

    # grab each chunck of data -
    reaction_abbreviation_string = data_array[row_index,1]
    reaction_name_string = data_array[row_index,2]
    reaction_chemical_string = data_array[row_index,3]
    median_turnover = data_array[row_index,4]
    source_string = data_array[row_index,5]
    reaction_ec_number = data_array[row_index,6]

    buffer *= "\t\t{\n"
    buffer *= "\t\t\t\"ec_number\":\"$(reaction_ec_number)\",\n"
    buffer *= "\t\t\t\"kcat_value\":\"$(median_turnover)\",\n"
    buffer *= "\t\t\t\"reaction_abbreviation\":\"$(reaction_abbreviation_string)\"\n"

    # close -
    if (row_index<number_of_rows)
        buffer *= "\t\t},\n"
    else
        buffer *= "\t\t}\n"
    end
end

# close -
buffer *= "\t]\n"
buffer *= "}\n"

# write -
outfile = open(path_to_json_file, "w")
write(outfile,buffer);
close(outfile);
