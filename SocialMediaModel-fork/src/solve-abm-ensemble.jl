function ABMsolveensemble(NT=200, N=10; savepoints = 4, init="4inf", p = ABMconstruct(), q= parameters(),sigma=0.1,save=true,savedata=false)
    (; X, Y, domain, dx, dy) = p
    (; n, L) = q
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    N_x = size(x_arr,1)
    N_y = size(y_arr,1)
    ys = zeros(2, 2, savepoints, N)
    zs = zeros(2, L, savepoints, N)
    us = zeros(N_x, N_y, 2, L, savepoints, N)
    savetimes = Int.(round.(LinRange(1, NT, savepoints)))
    states = [-1 1]

    for k=1:N
        xs, stateinfs, infs, meds, state, _ = ABMsolve(NT;  p=p, q=q, init=init,chosenseed=k)

        for m in 1:savepoints
            t = savetimes[m]
            x = xs[t]
            stateinf = stateinfs[t]
            inf = infs[t]
            media = meds[t]
            for i in 1:2
                for j in 1:L
                    xi = findall(x -> x == j, stateinf)
                    xm = findall(x -> x == states[i], state)
                    choice = intersect(xi, xm)
                    us[:,:,i,j,m,k] = [(1/n)*sumgaussian([X[i,j], Y[i,j]], x[choice,:],sigma=sigma) for i in 1:size(X,1), j in 1:size(X,2)]
                    zs[:,j,m,k] = inf[j,:]
                end
                ys[:,i,m,k] = media[i,:]
            end
        end

    end
    if savedata
        us2 = us;
        ys2 = ys;
        zs2 = zs;
        p2 = p;
        q2 = q;
        mkpath("data/pde-ensemble")
        @save string("data/pde-ensemble/abm_ensemble_",init,".jld2") us2 ys2 zs2 p2 q2
    end
    return us, ys, zs, (p,q)
end


function PDEsolveensemble(tmax=0.1, N=10; savepoints = 4, alg=nothing, init="4inf", p = PDEconstruct(), q= parameters(),save=true,savedata=false)
    (; N_x, N_y) = p
    L = q.L

    ys = zeros(2, 2, savepoints, N)
    zs = zeros(2, L, savepoints, N)
    us = zeros(N_x, N_y, 2, L, savepoints,  N)
    savetimes = LinRange(0, tmax, savepoints)

    # Threads.@threads
    local P
    for i=1:N
        sol, P = PDEsolve(tmax; alg=alg, init=init, p=p, q=q)

        for j in 1:savepoints
            u,y,z = sol2uyz(sol, savetimes[j])
            us[:,:,:,:,j,i] = u
            ys[:,:,j,i] = y
            zs[:,:,j,i] = z
        end

    end
    if savedata
        (p,q) = P
        mkpath("data/pde-ensemble")
        @save string("data/pde-ensemble/pde_ensemble_",init,".jld2") us ys zs p q
    end
    return us, ys, zs, P
end


function runensembles(N; NT=200, tmax=2., savepoints = 5, q=parameters(),sigma=0.1,save=true,savedata=false, init="4inf", name=init)
    # pde
    us, ys, zs, (p,q)= PDEsolveensemble(tmax, 1;savepoints=savepoints, init=init, q=q,save=save,savedata=savedata)
    # abm
    us2, ys2, zs2, (p2,q2) = ABMsolveensemble(NT,N; savepoints=savepoints, init=init, q=q,sigma=sigma,save=save,savedata=savedata)
    plotensemblesnapshots(us, ys, zs, (p,q), us2, ys2, zs2, (p2,q), tmax; name=name,save=save)
    return us, ys, zs, (p,q), us2, ys2, zs2, (p2,q)
end
