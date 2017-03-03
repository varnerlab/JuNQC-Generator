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

function generate_net_reaction(index_reactions::Array{Int64,1},data_dictionary::Dict{AbstractString,Any})

  # get the stoichiometric matrix -
  stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]

  # get list of metabolite symbols -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]

  # ok, combine the cols in index_reactions -
  coefficient_array = sum(stoichiometric_matrix[:,index_reactions],2)

  # which elememts are positive (products)?
  idx_product_array = find(coefficient_array.>0)

  # which elements are negative (reactants?)
  idx_reactant_array = find(coefficient_array.<0)

  # build the string ...
  net_reaction_buffer = ""
  for idx_reactant in idx_reactant_array

    metabolite_symbol = list_of_metabolite_symbols[idx_reactant]
    st_coeff = abs(coefficient_array[idx_reactant])

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
    st_coeff = abs(coefficient_array[idx_product])

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
