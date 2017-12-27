include("Include.jl")

# helper methods ...

# update the flux bounds -
function update_flux_bounds_array(time,time_step_size,state_array,data_dictionary)
end

# update the species bounds -
function update_species_bounds_array(time,time_step_size,state_array,data_dictionary)
end

# update the species array -
function time_stepping_routine(time,time_step_size,state_array,data_dictionary)
end

# main execution method -
function main(time_start::Float64,time_stop::Float64,time_step_size::Float64,data_dictionary::Dict{String,Any})

    # initialize -
    is_ok_to_continue = true
    current_time = time_start
    time_archive = Float64[]

    # how many states do we have => initialize state archive
    number_of_species = data_dictionary["number_of_species"]
    state_archive = zeros(1,number_of_species)

    # get the initial condition array -
    initial_condition_array = data_dictionary["initial_condition_array"]

    # setup current state -
    current_state = initial_condition_array

    # main time step loop -
    while (is_ok_to_continue)

        # update the flux bounds function -
        flux_bounds_array = update_flux_bounds_array(time,time_step_size,current_state,data_dictionary)

        # update the species bounds function -
        species_bounds_array = update_species_bounds_array(time,time_step_size,current_state,data_dictionary)

        # update the current state -
        current_state = time_stepping_routine(time,time_step_size,current_state,data_dictionary)

        # update time -
        current_time = current_time + time_step_size
        if (current_time>time_stop)
            is_ok_to_continue = false
        else

            # archive data, go around again -
            # archive time -
            push!(time_archive,current_time)

            # archive state -
            state_archive = [state_archive ; current_state]
        end
    end

    return (time_archive,state_archive[2:end,:])
end


#  ============ setup call to main ============================================ #

# setup the simulation time scale -
time_start = 0.0        # hr
time_step = 1.0         # hr
time_step_size = 0.1    # hr

# load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# call main -
main()
