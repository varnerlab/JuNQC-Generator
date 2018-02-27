function generate_problem_object(model_dictionary::Dict{String,Any},configuration_dictionary::Dict{String,Any})

    # Initilize an empty problem object -
    problem_object::ProblemObject = ProblemObject()

    # initialize -
    species_array::Array{SpeciesObject} = SpeciesObject[]
    reaction_array::Array{ReactionObject} = ReactionObject[]

    # get the list of species -
    list_of_species = model_dictionary["list_of_species"]
    for species_dictionary in list_of_species

        # Get the data -
        species_symbol = species_dictionary["species_symbol"]
        species_bound_type = species_dictionary["species_bound_type"]

        # Build new species object -
        species_object = SpeciesObject()
        species_object.species_symbol = species_symbol
        species_object.species_bound_type = Symbol(species_bound_type)

        # cache -
        push!(species_array,species_object)
    end

    # Build the reaction object from the JSON reaction dictionary -
    list_of_reactions = model_dictionary["list_of_reactions"]
    for reaction_dictionary in list_of_reactions

        # Get the reaction data -
        reaction_name = reaction_dictionary["reaction_name"]
        reaction_type_flag = reaction_dictionary["reaction_type_flag"]
        catalyst_lexeme = reaction_dictionary["catalyst_lexeme"]
        catalyst_ec_number = reaction_dictionary["catalyst_ec_number"]

        # build a new reaction object -
        reaction_object = ReactionObject()
        reaction_object.reaction_name = reaction_name
        reaction_object.reaction_type_flag = parse(Int,reaction_type_flag)

        # Hack to work w/code generator -
        if (reaction_type_flag == 0)
            reaction_object.reaction_type_symbol = :metabolic
        else
            reaction_object.reaction_type_symbol = :unknown
        end


        reaction_object.catalyst_lexeme = catalyst_lexeme
        reaction_object.catalyst_ec_number = catalyst_ec_number

        # Add the list of reactants -
        list_of_reactants_json = reaction_dictionary["list_of_reactants"]
        list_of_reactants_array = SpeciesObject[]
        for species_dictionary in list_of_reactants_json

            species_symbol = species_dictionary["symbol"]
            stcoeff = parse(Float64,species_dictionary["stoichiometry"])

            species_object = SpeciesObject()
            species_object.species_symbol = species_symbol
            species_object.stoichiometric_coefficient = stcoeff
            push!(list_of_reactants_array,species_object)
        end
        reaction_object.list_of_reactants = list_of_reactants_array

        # Add the list of products -
        list_of_products_json = reaction_dictionary["list_of_products"]
        list_of_products_array = SpeciesObject[]
        for species_dictionary in list_of_products_json

            species_symbol = species_dictionary["symbol"]
            stcoeff = parse(Float64,species_dictionary["stoichiometry"])

            species_object = SpeciesObject()
            species_object.species_symbol = species_symbol
            species_object.stoichiometric_coefficient = stcoeff
            push!(list_of_products_array,species_object)
        end
        reaction_object.list_of_products = list_of_products_array

        # cache -
        push!(reaction_array,reaction_object)
    end

    # Partition species -
    partition!(species_array)

    # set data on problem_object -
    problem_object.list_of_species = species_array
    problem_object.list_of_reactions = reaction_array

    # return#the problem_object -
    return problem_object
end


function generate_problem_object(metabolic_statement_vector::Array{VFFSentence},configuration_dictionary::Dict{String,Any})

  # Initilize an empty problem object -
  problem_object::ProblemObject = ProblemObject()

  # construct the array of species -
  species_array::Array{SpeciesObject} = build_species_list(metabolic_statement_vector,configuration_dictionary)

  # Partition species -
  partition!(species_array)

  # construct the array of reactions -
  reaction_array::Array{ReactionObject} = build_reaction_list(metabolic_statement_vector,configuration_dictionary)

  # partition the reactions -
  partition!(reaction_array)

  # set data on problem_object -
  problem_object.list_of_species = species_array
  problem_object.list_of_reactions = reaction_array

  # return#the problem_object -
  return problem_object
end

function build_species_list!(reaction_clause::AbstractString,list_of_species::Array{SpeciesObject},configuration_dictionary::Dict{String,Any})

  # We need to set the unbalanced species (type :unbalanced)
  unbalanced_species_suffix = configuration_dictionary["unbalanced_species_suffix"]["symbol"]

  if (contains(reaction_clause,"+") == true)

    # split around the +, and recursivley call me ..
    tmp_split_array = split(reaction_clause,"+")
    for fragment in tmp_split_array
      build_species_list!(fragment,list_of_species,configuration_dictionary)
    end

  else

    species_object::SpeciesObject = SpeciesObject()

    # ok, no +, but maybe still a stoichiometric_coefficient -
    if (contains(reaction_clause,"*") == true)

      symbol = strip(split(reaction_clause,"*")[end])
      coefficient = split(reaction_clause,"*")[1]

      if (symbol != "[]")

        # Build the species object -
        species_object.species_index = 0.0
        species_object.species_type = :metabolite

        # how many chars is the


        if (contains(symbol[end-1:end],unbalanced_species_suffix) == true)
          species_object.species_bound_type = :unbalanced
        else
          species_object.species_bound_type = :balanced
        end

        species_object.species_symbol = symbol
        species_object.stoichiometric_coefficient = parse(Float64,coefficient)
        species_object.species_compartment = :reactor

        # add to list -
        push!(list_of_species,species_object)
      end
    else

      symbol = strip(reaction_clause)
      coefficient = 1.0

      if (symbol != "[]")

        # Build the species object -
        species_object.species_index = 0.0
        species_object.species_type = :metabolite

        # R we unbalanced?
        species_object.species_bound_type = balanced_or_unbalanced(symbol,unbalanced_species_suffix)
        species_object.species_symbol = symbol
        species_object.stoichiometric_coefficient = coefficient
        species_object.species_compartment = :reactor

        # add to list -
        push!(list_of_species,species_object)
      end
    end
  end
end

function build_species_list(statement_vector::Array{VFFSentence},configuration_dictionary::Dict{String,Any})

  species_set::Set{AbstractString} = Set{AbstractString}()
  for vff_sentence in statement_vector

      # grab the reactant and prodict strings -
      reactant_string = vff_sentence.sentence_reactant_clause
      product_string = vff_sentence.sentence_product_clause

      # build the sets -
      reactant_set = build_species_set_from_clause(reactant_string)
      product_set = build_species_set_from_clause(product_string)

      # add these to the species_set -
      for item in reactant_set
        push!(species_set,item)
      end

      for item in product_set
        push!(species_set,item)
      end
  end

  # ok, I have a set of symbols = it should be unique, sort and then create an array of SpeciesObjets -
  tmp_species_array = AbstractString[]
  for species_symbol in species_set
    push!(tmp_species_array,species_symbol)
  end

  # sort this in place (alphbetical ..)
  sort!(tmp_species_array)

  # We need to set the unbalanced species (type :unbalanced)
  unbalanced_species_suffix = configuration_dictionary["unbalanced_species_suffix"]["symbol"]

  # Finally ... build the list of SpeciesObjects -
  list_of_species = SpeciesObject[]
  for (index,species_symbol) in enumerate(tmp_species_array)

    # species_symbol::AbstractString
    # species_index::Int
    # stoichiometric_coefficient::Float64
    # species_compartment::Symbol

    species_object = SpeciesObject()
    species_object.species_symbol = species_symbol
    species_object.species_index = index
    species_object.stoichiometric_coefficient = 0.0
    species_object.species_type = :metabolite
    species_object.species_compartment = :reactor

    # check - do we have an unbalanced species?
    species_object.species_bound_type = balanced_or_unbalanced(species_symbol,unbalanced_species_suffix)

    # push -
    push!(list_of_species,species_object)
  end

  # return my sorted list of species objects -
  return list_of_species
end

function build_reaction_list(statement_vector::Array{VFFSentence},configuration_dictionary::Dict{String,Any})

  # Initialize an empty array -
  reaction_array::Array{ReactionObject} = ReactionObject[]

  # Iterate through senetences, build reaction list -
  for vff_sentence in statement_vector

    # build an empty reaction object -
    reaction_object::ReactionObject = ReactionObject()

    # grab the reactant and product strings -
    reactant_string = vff_sentence.sentence_reactant_clause
    product_string = vff_sentence.sentence_product_clause
    enyzme_generation_flag = vff_sentence.sentence_type_flag
    reaction_type_flag = vff_sentence.sentence_type_flag
    catalyst_lexeme = vff_sentence.catalyst_lexeme
    catalyst_ec_number = vff_sentence.sentence_ec_number

    # recatants -
    list_of_reactants::Array{SpeciesObject} = SpeciesObject[]
    build_species_list!(reactant_string,list_of_reactants,configuration_dictionary)

    # products -
    list_of_products::Array{SpeciesObject} = SpeciesObject[]
    build_species_list!(product_string,list_of_products,configuration_dictionary)

    # populate -
    reaction_object.is_reaction_reversible = false
    reaction_object.enyzme_generation_flag = enyzme_generation_flag
    reaction_object.list_of_reactants = list_of_reactants
    reaction_object.list_of_products = list_of_products
    reaction_object.reaction_name = vff_sentence.sentence_name
    reaction_object.reaction_type_flag = reaction_type_flag
    reaction_object.catalyst_lexeme = catalyst_lexeme
    reaction_object.catalyst_ec_number = catalyst_ec_number

    # store -
    push!(reaction_array,reaction_object)
  end

  # return -
  return reaction_array
end

function build_species_set_from_clause(reaction_clause::AbstractString)

  species_symbol_set::Set{AbstractString} = Set{AbstractString}()

  # does this reaction clause have a +?
  if (contains(reaction_clause,"+") == true)

    # split around the + -
    tmp_split_array = split(reaction_clause,"+")
    for fragment in tmp_split_array

      if (contains(fragment,"*") == true)
        value = strip(split(fragment,"*")[end])
        if (value != "[]")
          push!(species_symbol_set,value)
        end
      else
        if (strip(fragment) != "[]")
          push!(species_symbol_set,strip(fragment))
        end
      end
    end
  else

    # ok, no +, but maybe still a stoichiometric_coefficient -
    if (contains(reaction_clause,"*") == true)

      value = strip(split(reaction_clause,"*")[end])
      if (value != "[]")
        push!(species_symbol_set,value)
      end
    else

      # no + -or- *, so just a bare species -
      if (strip(reaction_clause) != "[]")
        push!(species_symbol_set,strip(reaction_clause))
      end
    end
  end

  return species_symbol_set
end

function balanced_or_unbalanced(species_symbol,unbalanced_species_suffix)

  species_bound_type = :balanced
  if (length(species_symbol)>1 && contains(species_symbol[end-1:end],unbalanced_species_suffix) == true)
    species_bound_type = :unbalanced
  end

  return species_bound_type
end
