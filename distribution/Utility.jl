function show_flux_profile(flux_array::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what fluxes are > epsilon?
  idx_cutoff = find(flux_array.>epsilon)

  # what is the list of reaction strings?
  list_of_reaction_strings = data_dictionary["list_of_reaction_strings"]

  # create a list of reactions?
  list_of_flux_records = String[]
  for flux_index in idx_cutoff

    # key,value -
    key = list_of_reaction_strings[flux_index]
    value = flux_array[flux_index]
    record = "$(flux_index),$(key),$(value)"
    push!(list_of_flux_records,record)

  end

  return list_of_flux_records
end

function generate_mode_file_buffer(uptake_archive::Array{Float64,2},data_dictionary::Dict{AbstractString,Any})

  # Grab the list of metabolites -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]
  stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]

  # Write the modes mapping file -
  buffer = ""
  (number_of_species,number_of_modes) = size(uptake_archive)
  for mode_index in 1:number_of_modes

    # record -
    buffer *= "M$(mode_index)"

    # ok, for this mode, find the pivot index (if we have multiple, choose the first)
    idx_non_zero = find(uptake_archive[:,mode_index].<0)
    if (length(idx_non_zero)>1)
      idx_non_zero = idx_non_zero[1]
    end

    @show idx_non_zero

    idx_pivot = find(stoichiometric_matrix[idx_non_zero,:].<0)[1]
    buffer *= ",$(idx_pivot)"

    # grab the flux -
    reaction_string = generate_net_reaction_string(uptake_archive[:,mode_index],data_dictionary)

    buffer *=",$(reaction_string)"

    # what species are consumed by this mode?
    idx_reactants = find(uptake_archive[:,mode_index].<0.0)
    for reactant_index in idx_reactants
      metabolite_symbol = list_of_metabolite_symbols[reactant_index]
      buffer *=",$(metabolite_symbol)"
    end

    buffer *= "\n"
  end

  return buffer
end

function generate_net_reaction_string(uptake_array::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

  # get list of metabolite symbols -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]

  # which elememts are positive (products)?
  idx_product_array = find(uptake_array.>0)

  # which elements are negative (reactants?)
  idx_reactant_array = find(uptake_array.<0)

  # build the string ...
  net_reaction_buffer = ""
  for idx_reactant in idx_reactant_array

    metabolite_symbol = list_of_metabolite_symbols[idx_reactant]
    st_coeff = abs(uptake_array[idx_reactant])

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # add the arrow -
  net_reaction_buffer *= " --> "

  # write the trailing stuff -
  for idx_product in idx_product_array

    metabolite_symbol = list_of_metabolite_symbols[idx_product]
    st_coeff = abs(uptake_array[idx_product])

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # return -
  return net_reaction_buffer
end
