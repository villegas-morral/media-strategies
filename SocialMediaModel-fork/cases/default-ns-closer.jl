# No switching, influencers and media closer: at (±0.6,±0.6), T=2.
# Run from the SocialMediaModel-fork folder (plots go to figures/img/): include("cases/default-ns-closer.jl")

using SocialMediaModel


q = parameters(b=2., c=4.)

ABMsolveplot(NT=200, ts=[1 11 51 101 151 201], init="closer", q=q, name="default-ns-closer")
PDEsolveplot(tmax=2.0, ts=[0. 0.1 0.4 1.0 1.5 2.0], init="closer", q=q, name="default-ns-closer")
