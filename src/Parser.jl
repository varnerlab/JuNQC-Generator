function parse_vff_measured_species_statements(path_to_model_file::AbstractString)

  # initialize -
  measured_species_vector = SpeciesObject[]
  tmp_array::Array{AbstractString} = AbstractString[]
  desired_handler_symbol::Symbol = :measured_species_handler
  current_handler_symbol::Symbol = :metabolic_reaction_handler

  try

    # Open the model file, and read each line into a vector -
    open(path_to_model_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    counter = 1
    for sentence in tmp_array

      # Check the current pragma -
      if (contains(sentence,"#pragma") == true)

        # Extract the current handler -
        current_handler_symbol = extract_vff_handler_symbol(sentence)
      end

      # Depending upon the handler_symbol, we will process diff types of records -
      if (current_handler_symbol == desired_handler_symbol && contains(sentence,"#pragma") == false)

        measured_species_object = vff_measured_species_object_factory(sentence,current_handler_symbol)
        measured_species_object.species_index = counter
        push!(measured_species_vector,measured_species_object)

        # Update the counter -
        counter = counter + 1
      end
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return measured_species_vector
end


function parse_vff_metabolic_statements(path_to_model_file::AbstractString)

  # We are going to load the sentences in the file into a vector
  # if not a valid model file, then throw an error -
  sentence_vector = VFFSentence[]
  expanded_sentence_vector::Array{VFFSentence} = VFFSentence[]
  tmp_array::Array{AbstractString} = AbstractString[]
  desired_handler_symbol::Symbol = :metabolism_handler
  current_handler_symbol::Symbol = :metabolism_handler

  try

    # Open the model file, and read each line into a vector -
    open(path_to_model_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n") == 0:-1 && line != "")
            push!(tmp_array,chomp(line))
          end
      end
    end

    for sentence in tmp_array

      # convert substring to string -
      local_sentence = convert(String,sentence)

      # Check the current pragma -
      if (contains(local_sentence,"#pragma") == true)

        # Extract the current handler -
        current_handler_symbol = extract_vff_handler_symbol(local_sentence)
      end

        # Depending upon the handler_symbol, we will process diff types of records -
        if (current_handler_symbol == desired_handler_symbol && contains(local_sentence,"#pragma") == false)

            local_vff_sentence_array = vff_metabolic_sentence_factory(local_sentence,current_handler_symbol)
            for local_vff_sentence in collect(local_vff_sentence_array)
              push!(expanded_sentence_vector,local_vff_sentence)
            end
        elseif (current_handler_symbol == :transcription_handler && contains(local_sentence,"#pragma") == false)

            local_vff_sentence_array = vff_transcription_sentence_factory(local_sentence,current_handler_symbol)
            for local_vff_sentence in collect(local_vff_sentence_array)
              push!(expanded_sentence_vector,local_vff_sentence)
            end

        elseif (current_handler_symbol == :translation_handler && contains(local_sentence,"#pragma") == false)

            local_vff_sentence_array = vff_translation_sentence_factory(local_sentence,current_handler_symbol)
            for local_vff_sentence in collect(local_vff_sentence_array)
              push!(expanded_sentence_vector,local_vff_sentence)
            end

        elseif (current_handler_symbol == :degradation_handler && contains(local_sentence,"#pragma") == false)

            local_vff_sentence_array = vff_degradation_sentence_factory(local_sentence,current_handler_symbol)
            for local_vff_sentence in collect(local_vff_sentence_array)
              push!(expanded_sentence_vector,local_vff_sentence)
            end

        elseif (current_handler_symbol == :infrastructure_handler && contains(local_sentence,"#pragma") == false)

            local_vff_sentence_array = vff_infrastructure_sentence_factory(local_sentence,current_handler_symbol)
            for local_vff_sentence in collect(local_vff_sentence_array)
              push!(expanded_sentence_vector,local_vff_sentence)
            end
        end
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return expanded_sentence_vector
end


# ========================================================================================= #
# Helper functions -
# ========================================================================================= #
function extract_vff_handler_symbol(sentence::String)

    # Default - reaction handler -
    handler_symbol = :metabolic_reaction_handler

    # split the sentence -
    split_array = split(sentence,"::")

    # grab handler -
    handler_string = split_array[2]
    if (handler_string == "metabolism_handler")
        handler_symbol = :metabolism_handler
    elseif (handler_string == "transcription_handler")
        handler_symbol = :transcription_handler
    elseif (handler_string == "translation_handler")
        handler_symbol = :translation_handler
    elseif (handler_string == "degradation_handler")
        handler_symbol = :degradation_handler
    elseif (handler_string == "infrastructure_handler")
        handler_symbol = :infrastructure_handler
    end

    # return the handler -
    return handler_symbol
end

function vff_transcription_sentence_factory(sentence::String,handler_symbol::Symbol)
    return vff_default_sentence_factory(sentence,handler_symbol)
end

function vff_translation_sentence_factory(sentence::String,handler_symbol::Symbol)
    return vff_default_sentence_factory(sentence,handler_symbol)
end

function vff_degradation_sentence_factory(sentence::String,handler_symbol::Symbol)
    return vff_default_sentence_factory(sentence,handler_symbol)
end

function vff_infrastructure_sentence_factory(sentence::String,handler_symbol::Symbol)
    return vff_default_sentence_factory(sentence,handler_symbol)
end

function vff_default_sentence_factory(sentence::String,handler_symbol::Symbol)

    # initialize -
    sentence_vector = VFFSentence[]

    # Ok, so now we have the array for sentences -
    vff_sentence = VFFSentence()
    vff_sentence.original_sentence = sentence

    # split the sentence -
    split_array = split(sentence,",")
    vff_sentence.sentence_name = split_array[1]
    vff_sentence.sentence_type_flag = parse(Int,convert(String,split_array[2]))
    vff_sentence.catalyst_lexeme = convert(String,split_array[3])
    vff_sentence.sentence_reactant_clause = split_array[4]
    vff_sentence.sentence_product_clause = split_array[5]
    vff_sentence.sentence_reverse_bound = parse(Float64,split_array[6])
    vff_sentence.sentence_forward_bound = parse(Float64,split_array[7])
    vff_sentence.sentence_handler = handler_symbol
    vff_sentence.sentence_delimiter = ','
    vff_sentence.sentence_ec_number = "[]"

    # add sentence to sentence_vector -
    push!(sentence_vector,vff_sentence)

    # Check - is this reversible?
    if (vff_sentence.sentence_reverse_bound == -Inf)

      # ok, so we have a reversible reaction -
      # first change lower bound to 0 -
      vff_sentence.sentence_reverse_bound = 0.0

      # create a new copy of sentence -
      vff_sentence_copy = deepcopy(vff_sentence)
      vff_sentence_copy.sentence_name = (vff_sentence_copy.sentence_name)*"_reverse"
      vff_sentence_copy.sentence_reactant_clause = vff_sentence.sentence_product_clause
      vff_sentence_copy.sentence_product_clause = vff_sentence.sentence_reactant_clause
      vff_sentence_copy.sentence_handler = vff_sentence.sentence_handler
      vff_sentence_copy.sentence_type_flag = vff_sentence.sentence_type_flag
      vff_sentence_copy.sentence_ec_number = "[]"

      # add sentence and sentence copy to the expanded_sentence_vector -
      push!(sentence_vector,vff_sentence_copy)
    end

    return sentence_vector
end

function vff_metabolic_sentence_factory(sentence::String,handler_symbol::Symbol)

    # Initialize -
    sentence_vector = VFFSentence[]

    # Ok, so now we have the array for sentences -
    vff_sentence = VFFSentence()
    vff_sentence.original_sentence = sentence

    # split the sentence -
    split_array = split(sentence,",")
    vff_sentence.sentence_name = split_array[1]
    vff_sentence.sentence_type_flag = parse(Int,convert(String,split_array[2]))
    vff_sentence.catalyst_lexeme = convert(String,split_array[3])
    vff_sentence.sentence_reactant_clause = split_array[4]
    vff_sentence.sentence_product_clause = split_array[5]
    vff_sentence.sentence_reverse_bound = parse(Float64,split_array[6])
    vff_sentence.sentence_forward_bound = parse(Float64,split_array[7])
    vff_sentence.sentence_handler = handler_symbol
    vff_sentence.sentence_delimiter = ','
    vff_sentence.sentence_ec_number = split_array[8]

    # add sentence to sentence_vector -
    push!(sentence_vector,vff_sentence)

    # Check - is this reversible?
    if (vff_sentence.sentence_reverse_bound == -Inf)

        # ok, so we have a reversible reaction -
        # first change lower bound to 0 -
        vff_sentence.sentence_reverse_bound = 0.0

        # create a new copy of sentence -
        vff_sentence_copy = deepcopy(vff_sentence)
        vff_sentence_copy.sentence_name = (vff_sentence_copy.sentence_name)*"_reverse"
        vff_sentence_copy.sentence_reactant_clause = vff_sentence.sentence_product_clause
        vff_sentence_copy.sentence_product_clause = vff_sentence.sentence_reactant_clause
        vff_sentence_copy.sentence_handler = vff_sentence.sentence_handler
        vff_sentence_copy.sentence_type_flag = vff_sentence.sentence_type_flag
        vff_sentence_copy.sentence_ec_number = vff_sentence.sentence_ec_number

        # add sentence and sentence copy to the expanded_sentence_vector -
        push!(sentence_vector,vff_sentence_copy)
    end

    return sentence_vector
end
