function paperfigures()

    # single ABM realization
    ABMsolveplot();

    # for appendix
    ABMsolveplot(init="4infstrongb",q=parameters(a=1,b=3,c=1))
    ABMsolveplot(init="4infstronga",q=parameters(a=3,b=1,c=1))
    # further good experiments for appendix
    # ABMsolveplot(init="4infweakgamma2.5",q=parameters(frictionI=2.5))
    # ABMsolveplot(init="4infweakgamma1",q=parameters(frictionI=1))
    # ABMsolveplot(init="4infweakgamma2.5noise",q=parameters(sigmahat=1,frictionI=2.5),seed=3)

    # ensemble
    # saved by runensembles(...; savedata=true)
    @load string("data/pde-ensemble/pde_ensemble_4inf.jld2") us ys zs p q
    @load string("data/pde-ensemble/abm_ensemble_4inf.jld2") us2 ys2 zs2 p2 q2
    plotensemblesnapshots(us, ys, zs, (p,q), us2, ys2, zs2, (p2,q), 2.; name="4inf",save=true)

    # optimizations
    infrightcorner();
    infcounteraction();
    medcounteraction();
end
