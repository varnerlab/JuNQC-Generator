
type VFFSentence

  original_sentence::AbstractString
  sentence_name::AbstractString
  sentence_type_flag::Int
  sentence_reactant_clause::AbstractString
  sentence_product_clause::AbstractString
  sentence_reverse_bound::Float64
  sentence_forward_bound::Float64
  sentence_delimiter::Char
  sentence_handler::Symbol
  sentence_ec_number::AbstractString

  # holds the symbol for the catalyst -
  catalyst_lexeme::String

  function VFFSentence()
    this = new()
  end
end


type ProgramComponent

  component_type::Symbol
  filename::AbstractString
  buffer::AbstractString
  matrix::Array{Float64,2}

  function ProgramComponent()
    this = new()
  end

end

type SpeciesObject

  species_type::Symbol
  species_bound_type::Symbol
  species_symbol::AbstractString
  species_index::Int
  stoichiometric_coefficient::Float64
  species_compartment::Symbol

  function SpeciesObject()
    this = new()
  end
end

type ReactionObject

  enyzme_generation_flag::Int
  reaction_type_symbol::Symbol
  reaction_index::Int
  is_reaction_reversible::Bool
  reaction_name::AbstractString
  list_of_reactants::Array{SpeciesObject}
  list_of_products::Array{SpeciesObject}

  # holds the type flag -
  reaction_type_flag::Int
  catalyst_lexeme::String
  catalyst_ec_number::String

  function ReactionObject()
    this = new()
  end
end


type ProblemObject

    enzyme_kinetics_dictionary::Dict{AbstractString,Any}
    configuration_dictionary::Dict{AbstractString,Any}
    list_of_species::Array{SpeciesObject}
    list_of_reactions::Array{ReactionObject}

  function ProblemObject()
    this = new()
  end
end
