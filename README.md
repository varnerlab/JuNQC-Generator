## Network Quality Control in Julia (JuNQC)
JuNQC-Generator is a code generation system for static stoichiometric models (SSM) written in the [Julia](http://julialang.org) programming language.
JuNQC-Generator transforms a simple comma/space delimited flat-file into fully commented metabolic model code in the [Julia](http://julialang.org)  programming language.
The generated code uses the GLPK solver to solve the metabolic flux balance analysis program.

### Installation and Requirements
You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/JuNQC-Generator.git

or

	$ git clone https://github.com/varnerlab/JuNQC-Generator.git

To execute a code generation job, Julia must be installed on your machine along with the ``ArgParse`` and ``JSON`` Julia packages.
Julia can be downloaded/installed on any platform.
The required [Julia](http://julialang.org) packages can be installed by executing the commands:

	julia> Pkg.add("ArgParse")

and

	julia> Pkg.add("JSON")

in the Julia REPL.  
Lastly, the generated code uses the Julia plugin for the GLPK linear programming solver.
To install the GLPK program issue the command:

  julia> Pkg.add("GLPK")

in the Julia REPL.

### Problem generation
To generate a SSM, issue the command ``make_julia_model.jl`` from the command line:

	$ julia make_julia_model.jl -m <input path> -o <output path>

The ``make_julia_model.jl`` command takes two command line arguments:

Argument | Required | Default | Description
--- | --- | --- | ---
-m | Yes	| none | Path to model input file
-o | No	| current directory | Path where files are written
