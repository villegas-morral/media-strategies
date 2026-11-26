# Default without switching: 4 influencers and 2 media at (±1,±1), T=2.
# Run from the SocialMediaModel-fork folder (plots go to img-default-ns/): include("cases/default-ns.jl")

using SocialMediaModel


q = parameters(b=2., c=4.)

ABMsolveplot(NT=200, ts=[1 11 51 101 151 201], init="4inf", q=q, name="default-ns")
PDEsolveplot(tmax=2.0, ts=[0. 0.1 0.4 1.0 1.5 2.0], init="4inf", q=q, name="default-ns")
