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

  # ok, combine the cols in index_reactions -
  

end
