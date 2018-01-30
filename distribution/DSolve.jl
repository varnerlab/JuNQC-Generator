include("Include.jl")

# helper methods ...
function archive_solution_array(state_archive,state_array)

    # add an addional row to the archive w/the new solution -
    (nrows,ncols) = size(state_archive)

    # initialize new archive -
    new_archive_array = zeros(nrows+1,ncols)

    # copy old data into new archive -
    for row_index = 1:nrows
        for col_index = 1:ncols
            new_archive_array[row_index,col_index] = state_archive[row_index,col_index]
        end
    end

    # add new row -
    for col_index = 1:ncols
        new_archive_array[end,col_index] = state_array[col_index]
    end

    # return -
    return new_archive_array
end

# update the flux bounds -
function update_flux_bounds_array(time,time_step_size,state_array,data_dictionary)
    return user_update_flux_bounds(time,time_step_size,state_array,data_dictionary)
end

# update the species bounds -
function update_species_bounds_array(time,time_step_size,state_array,data_dictionary)
    return user_update_species_bounds(time,time_step_size,state_array,data_dictionary)
end

# update the state array -
function time_stepping_routine(time,time_step_size,state_array,data_dictionary)
    return user_time_stepping_routine(time,time_step_size,state_array,data_dictionary)
end

# update the data dictionary -
function update_model_data_dictionary(time_start,time_stop,time_step_size,data_dictionary)
    return user_update_data_dictionary(time_start,time_stop,time_step_size,data_dictionary)
end

# main execution method -
function main(time_start::Float64,time_stop::Float64,time_step_size::Float64,data_dictionary::Dict{AbstractString,Any})

    # initialize -
    is_ok_to_continue = true
    current_time = time_start
    time_archive = Float64[]

    # how many states do we have => initialize state archive
    number_of_species = data_dictionary["number_of_species"]
    state_archive = zeros(1,number_of_species)

    # how many fluxes do we have?
    number_of_reactions = data_dictionary["number_of_reactions"]
    flux_archive = zeros(1,number_of_reactions)

    # get the initial condition array -
    initial_condition_array = data_dictionary["species_abundance_array"]

    # setup current state -
    current_state = initial_condition_array

    # archive the state -
    for col_index = 1:number_of_species
        state_archive[end,col_index] = current_state[col_index]
    end

    # archive the time -
    push!(time_archive,time_start)

    # main time step loop -
    while (is_ok_to_continue)

        # update the flux bounds function -
        flux_bounds_array = update_flux_bounds_array(time,time_step_size,current_state,data_dictionary)

        # update the species bounds function -
        species_bounds_array = update_species_bounds_array(time,time_step_size,current_state,data_dictionary)

        # update the species and flux bounds arrays -
        data_dictionary["default_flux_bounds_array"] = flux_bounds_array
        data_dictionary["species_bounds_array"] = species_bounds_array

        # update the current state -
        (current_state,flux_array) = time_stepping_routine(time,time_step_size,current_state,data_dictionary)

        # update time -
        current_time = current_time + time_step_size
        if (current_time>time_stop)
            is_ok_to_continue = false
        else

            # archive data, go around again -

            # archive time -
            push!(time_archive,current_time)

            # archive state -
            state_archive = archive_solution_array(state_archive,current_state)

            # archive the flux -
            flux_archive = archive_solution_array(flux_archive,flux_array)
        end
    end

    return (time_archive,state_archive,flux_archive)
end


#  ============ setup call to main ============================================ #

# setup the simulation time scale -
time_start = 0.0        # hr
time_stop = 1.0         # hr
time_step_size = 0.1    # hr

# load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# update the data dictionary w/model specific data -
data_dictionary = update_model_data_dictionary(time_start,time_stop,time_step_size,data_dictionary)

# call main -
(time_archive,state_archive,flux_archive) = main(time_start,time_stop,time_step_size,data_dictionary)

# output method -
# ...
#  ============ setup call to main ============================================ #
