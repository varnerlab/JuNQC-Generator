# script to convert a VFF extended format file to JSON -

using JSON
include("../Types.jl")
include("../Macros.jl")
include("../Parser.jl")
include("../Common.jl")
include("../Problem.jl")


# helper methods -
function encode_species_object(species_object::SpeciesObject)

    # initialize -
    buffer = ""

    # grab the data -
    species_type::Symbol = species_object.species_type
    species_bound_type::Symbol = species_object.species_bound_type
    species_symbol::AbstractString = species_object.species_symbol
    species_compartment::Symbol = species_object.species_compartment

    # encode -
    buffer *= "\t\t\t\"species_symbol\":\"$(species_symbol)\",\n"
    buffer *= "\t\t\t\"species_bound_type\":\"$(species_bound_type)\",\n"
    buffer *= "\t\t\t\"species_compartment\":\"$(species_compartment)\"\n"

    # return -
    return buffer
end

function encode_list_of_species(list_of_species)

    # init -
    buffer = ""

    # fill the buffer -
    number_of_species = length(list_of_species)
    species_counter = 1
    for species_object in list_of_species

        # get data from the species object -
        species_symbol = species_object.species_symbol
        stcoeff = species_object.stoichiometric_coefficient

        # encode -
        buffer *= "\t\t\t\t{\n"
        buffer *= "\t\t\t\t\t\"symbol\":\"$(species_symbol)\",\n"
        buffer *= "\t\t\t\t\t\"stoichiometry\":\"$(stcoeff)\"\n"
        buffer *= "\t\t\t\t}"

        # add a comma?
        if (species_counter<number_of_species)
            buffer *= ",\n"
        else
            buffer *= "\n"
        end

        # update the counter -
        species_counter = species_counter + 1
    end

    # return -
    return buffer
end

function main()

    # path to the vff file -
    path_to_vff_model_file = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/julia_work/core_ecoli_model_dipp_v1/Network.net"   # Path to the vff that you want to convert
    path_to_json_model_file = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/julia_work/core_ecoli_model_dipp_v1/Network.json"  #

    # Ok, create a collection of metabolic statements -
    metabolic_statement_vector::Array{VFFSentence} = parse_vff_metabolic_statements(path_to_vff_model_file)

    # Load the JSON configuration file -
    config_dict = JSON.parsefile("../../config/Configuration.json")

    # Load the enzyme kinetics JSON file -
    enzyme_kinetics_dictionary = JSON.parsefile("../../config/Enzymes.json")

    # Generate the problem object -
    problem_object = generate_problem_object(metabolic_statement_vector,config_dict)
    problem_object.configuration_dictionary = config_dict
    problem_object.enzyme_kinetics_dictionary = enzyme_kinetics_dictionary

    # get list of species and reactions -
    list_of_species::Array{SpeciesObject} = problem_object.list_of_species
    list_of_reactions::Array{ReactionObject} = problem_object.list_of_reactions

    # what is the size of the system?
    number_of_species = length(list_of_species)
    number_of_reactions = length(list_of_reactions)


    # Initialize the buffer -
    buffer = ""
    buffer *= "{\n"
    buffer *= "\t\"metadata\":{\n"
    buffer *= "\t},\n"
    buffer *= "\n"

    # Encode list of species -
    buffer *= "\t\"list_of_species\":[\n"

        # go through species objects, grab the data and then re-encode in JSON
        species_counter = 1
        for species_object in list_of_species

            buffer *= "\t\t{\n"
            buffer *= encode_species_object(species_object)
            buffer *= "\t\t}"

            # check: are we done?
            if (species_counter<number_of_species)
                buffer *=",\n"
            else
                buffer *= "\n"
            end

            # update the counter, and go around again -
            species_counter = species_counter + 1
        end

    buffer *= "\t],\n"

    # Encode list of reactions -
    buffer *= "\t\"list_of_reactions\":[\n"

        # go through each reaction object, grab the data, and then re-encode in a JSON structure
        reaction_counter = 1
        for reaction_object in list_of_reactions

            # Get the data for this reaction -
            enyzme_generation_flag = reaction_object.enyzme_generation_flag
            reaction_type = reaction_object.reaction_type_symbol
            reaction_index = reaction_object.reaction_index
            is_reaction_reversible = reaction_object.is_reaction_reversible
            reaction_name = reaction_object.reaction_name

            list_of_reactants = reaction_object.list_of_reactants::Array{SpeciesObject}
            list_of_products = reaction_object.list_of_products::Array{SpeciesObject}

            # holds the type flag -
            reaction_type_flag = reaction_object.reaction_type_flag
            catalyst_lexeme = reaction_object.catalyst_lexeme
            catalyst_ec_number = reaction_object.catalyst_ec_number

            # encode -
            buffer *= "\t\t{\n"
            buffer *= "\t\t\t\"reaction_name\":\"$(reaction_name)\",\n"
            buffer *= "\t\t\t\"reaction_type_flag\":\"$(reaction_type_flag)\",\n"
            buffer *= "\t\t\t\"reaction_type\":\"$(reaction_type)\",\n"
            buffer *= "\t\t\t\"catalyst_lexeme\":\"$(catalyst_lexeme)\",\n"
            buffer *= "\t\t\t\"catalyst_ec_number\":\"$(catalyst_ec_number)\",\n"

            # do the list of species -
            # reactants -
            buffer *="\t\t\t\"list_of_reactants\":[\n"

            # get reactant block -
            reactant_block = encode_list_of_species(list_of_reactants)
            buffer *= "$(reactant_block)"
            buffer *="\t\t\t],\n"

            # products -
            buffer *="\t\t\t\"list_of_products\":[\n"

            # get reactant block -
            product_block = encode_list_of_species(list_of_products)
            buffer *= "$(product_block)"
            buffer *="\t\t\t]\n"
            buffer *="\t\t}"

            # check: are we done?
            if (reaction_counter<number_of_reactions)
                buffer *=",\n"
            else
                buffer *= "\n"
            end

            # update the counter, and go around again -
            reaction_counter = reaction_counter + 1
        end

    buffer *= "\t]\n"
    buffer *= "}\n"

    # write -
    outfile = open(path_to_json_model_file, "w")
    write(outfile,buffer);
    close(outfile);
end

main()
