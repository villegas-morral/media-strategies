""" This is the entry file, to be included, for
using this not as a package but a simple project """

using Revise
using LinearAlgebra


LinearAlgebra.BLAS.set_num_threads(parse(Int,get(ENV, "SLURM_JOB_CPUS_PER_NODE", "20")))
print("Running on ", string(Threads.nthreads())," threads with ",string(LinearAlgebra.BLAS.get_num_threads())," BLAS threads.")


includet("setting-general.jl")
includet("solve-pde.jl")
includet("solve-abm.jl")
includet("solve-abm-ensemble.jl")
includet("plots.jl")
includet("controlstrategies.jl")
includet("setting-control.jl")
includet("tests.jl")
includet("plots-paper.jl")
