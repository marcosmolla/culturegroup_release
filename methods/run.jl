#!/usr/bin/env julia

### Load libraries
using Distributions
using StatsBase
using RCall
R"library(igraph)"

# Test whether this is run on the cluster, or locally
global cluster = any(collect(keys(ENV)).=="SLURM_ARRAY_JOB_ID");

# IF CLUSTER (Used when running on an HPC or HTC cluster)
if cluster
	println("You are using SLRUM, your SLURM job ID is ", ENV["SLURM_JOB_ID"])
  global paths="/home/uma/culturegroup/XXXX"
  # Initialise Methods
  include(string(paths, "/methods.jl"));
  global queue = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]);
  # Initialise Parameters
  include(string(paths, "/main.jl"));
  # Run Simulation
  @time runit();
else # (Used when running locally)
  # Initialise Methods
  # include("/Users/marco/Documents/Programming/julia/marcosProjects/culturegroup/methodsJL/methods.jl")
  include("/Users/marco/Documents/Programming/julia/marcosProjects/culturegroup/methods/methods.jl");
  # Run Simulation
  for q in 1 #:1120
    global queue = q
	# Initialise Parameters
    @time include("/Users/marco/Documents/Programming/julia/marcosProjects/culturegroup/methods/main.jl");
	# Run Simulation
	@time runit();
	println(" with ", queue,".")
  end
end
