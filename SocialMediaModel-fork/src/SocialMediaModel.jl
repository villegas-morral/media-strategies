module SocialMediaModel
    using TerminalLoggers
    function __init__()
        Base.CoreLogging.global_logger(TerminalLogger()) # progress bars for PDE solve and ABM loop
        Plots.default(fontfamily="Computer Modern") # LaTeX-style plots (titles are LaTeX strings in plots.jl)
    end
    include("setting.jl")
    include("pde.jl")
    include("abm.jl")
    include("pde-ensemble.jl")
    include("plots.jl")
    include("controlstrategies.jl")
    include("strategyexperiments.jl")
    include("tests.jl")
    include("plots-paper.jl")

    export  ABMsolve,  ABMsolveplot, PDEsolve, PDEsolveplot, PDEplotsnapshots, ABMplotsnapshots, PDEgifsingle,
     ABMplotfollowernumbers, ABMgifsingle, plotensemblesnapshots, PDEsolveensemble, ABMsolveensemble, runensembles,
     parameters, parametersstronginf, ABMconstruct, PDEconstruct, PDEconstructcoarse, PDEconstructmeso, parameterscontrol,
     PDEsolvefixedtargets, PDEsolveopt, GLmed, GLinf, infrightcorner, infcounteraction, medcounteraction,
     testabm, testpde, testensemble, testcontrol, testoptimization, runtests, paperfigures
end
