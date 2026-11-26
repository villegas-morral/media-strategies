# Default (with switching) up to T=10. The PDE takes about 6 minutes.
# Run from the SocialMediaModel-fork folder (plots go to img-default-t10/): include("cases/default-t10.jl")

using SocialMediaModel


q = parameters(b=2., c=4., switching=true)

ABMsolveplot(NT=1000, ts=[1 11 51 101 501 1001], init="4inf", q=q, name="default-t10")
PDEsolveplot(tmax=10.0, ts=[0. 0.1 0.4 1.0 5.0 10.0], init="4inf", q=q, name="default-t10")
