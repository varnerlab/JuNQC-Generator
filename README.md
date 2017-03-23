## Network Quality Control in Julia (JuNQC)
JuNQC-Generator is a code generation system for static stoichiometric models (SSM) written in the [Julia](http://julialang.org) programming language.
JuNQC-Generator transforms a simple comma/space delimited flat-file into fully commented metabolic model code in the [Julia](http://julialang.org)  programming language.
The generated code uses the GLPK solver to solve the metabolic flux balance analysis program.

### Installation and Requirements
You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/JuNQC-Generator.git

or

	$ git clone https://github.com/varnerlab/JuNQC-Generator.git

To execute a code generation job, Julia must be installed on your machine along with the  [ArgParse](https://github.com/carlobaldassi/ArgParse.jl) and [JSON](https://github.com/JuliaIO/JSON.jl) Julia packages.
Julia can be downloaded/installed on any platform.
The required [Julia](http://julialang.org) packages can be installed by executing the commands:

	julia> Pkg.add("ArgParse")

and

	julia> Pkg.add("JSON")

in the Julia REPL. Lastly, the generated code uses the Julia plugin for the [GLPK](https://github.com/JuliaOpt/GLPK.jl) linear programming solver. To install the GLPK program issue the command:

  	julia> Pkg.add("GLPK")

in the Julia REPL.

### How do I generate FBA code?
To generate flux balance analysis code, issue the command ``make_julia_model.jl`` from the command line:

	$ julia make_julia_model.jl -m <input path> -o <output path>

The ``make_julia_model.jl`` command takes two command line arguments:

Argument | Required | Default | Description
--- | --- | --- | ---
-m | Yes	| none | Path to model input file
-o | No	| current directory | Path where files are written

#### What is the structure of model input file?
JuNQC converts model files in the Varner Flat File (VFF) format to program code. VFF is a simple comma delimited flat file with a record structure:
	
	Reaction_name (unique),Reactants,Products,reverse_string,forward_string
	
For example, the first few steps in glycolysis are written as:

	// ----------------------------------------------------------------------------------- //
	// metabolic reactions -
	#pragma::metabolic_reaction_handler 
	//
	// Record:
	// name (unique),reactants,products,reverse,forward
	//
	// Notes:
	// [] = system (not tracked)
	// ----------------------------------------------------------------------------------- //

	
	// Phosphoglucose isomerase (pgi)
	R_pgi,M_g6p_c,M_f6p_c,-inf,inf

	// Phosphofructokinase (pfk)
	R_pfk,M_atp_c+M_f6p_c,M_adp_c+M_fdp_c+M_h_c,0,inf

	// Fructose-1,6-bisphosphate aldolase (fdp)
	R_fdp,M_fdp_c+M_h2o_c,M_f6p_c+M_pi_c,0,inf
	
	// R_FUMt2_2::Fumarate transport via proton symport (2 H)
	R_FUMt2_2,M_fum_e+2.0*M_h_e,M_fum_c+2.0*M_h_c,0,inf
	
Comments are written with C++ style `//` and can contain any text e.g., references, notes etc. Chemical species must not contain spaces, or special characters. Stoichiometric coefficients other than `1.0` are added to a species using the `*` symbol, e.g., `2.0*M_h_e`. Lastly, while the header comment is optional, the:
	
	#pragma::metabolic_reaction_handler

parser directive must appear at the top of the metabolic reaction section of the file. 

#### What files get generated and what do they do?
JuNQC generates files to setup and solve a flux balance analysis problem:

file | description 
--- | ---
Atom.txt | Contains molecular formula library (in CSV format) for metabolites. `M x 6` where each row encodes a metabolite symbol e.g., `M_g6p_c` and each describes the formula in terms of `C,H,N,O,P` and `S`.  
DataDictionary.jl | Encodes the species (`M x 2`) and reaction (`R x 2`) bounds arrays and the objective array (`R x 1`). Data is stored in a [Julia dictionary](http://docs.julialang.org/en/stable/stdlib/collections/?highlight=dict#Base.Dict) type and can be accessed through the appropriate key.
Debug.txt | List of reactions and species used to generate the model code
FluxDriver.jl | Julia interface with the [GLPK](https://github.com/JuliaOpt/GLPK.jl) solver. Users should `NEVER, UNDER ANY CIRCUMSTANCES, EVER` edit this file. 
Include.jl | Encodes all the include statements for the project. Should be included at the top of top-level driver scripts. 
Network.dat | Stoichiometric array for the model
Solve.jl | Default top-level driver implementation. 
Utility.jl | Encodes utility functions for interacting with estimated flux profiles, and reading/writing auxiliary files such as the Atom matrix. 
