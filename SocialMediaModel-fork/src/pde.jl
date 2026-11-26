using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings


"""Grid matrices (built once by PDEconstruct in setting.jl)"""

# interaction kernel on the grid: w(x,x') = exp(-|x'-x|) and k(x,x') = w(x,x') (x'-x)
function generate_K(grid_points)
    N, d = size(grid_points)
    @assert d == 2
    k = zeros(N, N, 2)
    w = zeros(N, N)

    @inbounds for j in 1:N, i in 1:N
        x = grid_points[j,1] - grid_points[i,1] # first component of x' -x
        y = grid_points[j,2] - grid_points[i,2] # 2nd component of x' -x
        w[i,j] = ww = exp(-sqrt(x^2 + y^2))
        k[i,j, 1] = ww * x
        k[i,j, 2] = ww * y
    end

    return k, w
end

# second-difference matrix for the diffusion term (Neumann at the walls: first and last diagonal entries -1)
function second_derivative((N_x,N_y), (dx, dy))
    # here matrix should have different shape for Ny not equal Nx
    M = Tridiagonal(ones(N_x-1), fill(-2., N_x), ones(N_x-1))
    # boundary conditions change in time to ensure flux balance, they are added when solving pde
    M[1,1] = -1.
    M[end,end] = -1.
    M .= (1/dx^2)*M

    return M
end

# centred first-difference matrix for the transport term (one-sided at the walls)
function centered_difference((N_x,N_y), (dx, dy))
    # centered first difference for force, doesnt work for different x, y grids
    C = 1/(2*dx)*Tridiagonal(-ones(N_x-1), zeros(N_x), ones(N_x-1))
    # at the boundary a one-sided difference scheme is used
    C[1,1:2] = 1/(dx)* [-1,1]
    C[end,end-1:end] = 1/(dx)* [-1,1]

    return C
end


"""Set initial conditions"""

# DEFAULT and closer

# uniform density on [-2,2]^2, each quadrant follows the influencer in it, split equally over both media;
# media at ±(pos,pos), influencers at (±pos,±pos): pos = 1 (4inf) or 0.6 (closer)
function uniforminit((p,q); pos=1.0)
    (; N_x , N_y,dV, domain, dx) = p
    (; L, M) = q
    @assert M == 2
    @assert L == 4

    rho_0 = zeros(N_x, N_y, 2, 4)
    mid_y = Int(round(N_y/2))
    mid_x = Int(round(N_x/2))
    start_x = Int(round((domain[1,2] - 2)/dx + 1))
    end_x = N_x - start_x +1
    rho_0[start_x:mid_x, start_x:mid_y,:,4] .= 1
    rho_0[start_x:mid_x, mid_y+1:end_x,:,2] .= 1
    rho_0[mid_x+1:end_x, start_x:mid_y,:,3] .= 1
    rho_0[mid_x+1:end_x, mid_y+1:end_x,:,1] .= 1

    u0 = rho_0/(sum(rho_0)*dV)
    y2_0 = [pos,pos]
    y1_0 = [-pos,-pos]
    y0 = [y1_0  y2_0]
    z4_0 = [-pos,-pos]
    z2_0 = [-pos,pos]
    z3_0 = [pos,-pos]
    z1_0 = [pos,pos ]

    z0 = [z1_0  z2_0 z3_0 z4_0]
    controlled_inf = zeros(L)
    controlled_med = zeros(2)
    return ArrayPartition(u0,y0,z0), controlled_inf, controlled_med
end

# 2 INFLUENCERS

# uniform density, above the line y=-x follows influencer 1 at (1,-1), below influencer 2 at (-1,1)
function uniforminit_2inf((p,q))
    (; N_x , N_y, dV, domain, dx, grid_points) = p
    (; L, M) = q
    @assert M == 2
    @assert L == 2

    rho_0 = zeros(N_x, N_y, 2, L)

    # above y=-x → influencer 1, below → influencer 2
    for idx in 1:N_x*N_y
        gx, gy = grid_points[idx, 1], grid_points[idx, 2]
        if gy > -gx
            rho_0[mod1(idx, N_x), div(idx-1, N_x)+1, :, 1] .= 1
        else
            rho_0[mod1(idx, N_x), div(idx-1, N_x)+1, :, 2] .= 1
        end
    end

    u0 = rho_0 / (sum(rho_0) * dV)
    y0 = [-1. 1.; -1. 1.]
    z0 = [1. -1.; -1. 1.]

    controlled_inf = zeros(L)
    controlled_med = zeros(M)
    return ArrayPartition(u0, y0, z0), controlled_inf, controlled_med
end

# RANDOM (unused)

# n random individuals, each smoothed into a gaussian; influencers at the mean of their followers
function randominit((p,q))
    (; N_x , N_y,dV, grid_points,domain) = p
    (; L, n, M) = q
    @assert M == 2
    @assert L == 4

    random_pos = rand(n,2).*(domain[1,2]-domain[1,1]) .+domain[1,1] # n uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, L)
    counts = zeros(L,2)
    z0 = zeros(2,L)

    function add_agent(agenti, infi, state,  rho_0, counts, z0)
        rho_0[:,:,state,infi] += reshape([gaussian(grid_points[j,:], random_pos[agenti,:]) for j in 1:N_x*N_y], N_x, N_y)
        counts[infi, state] += 1
        z0[:,infi] += random_pos[agenti,:]
        return rho_0, counts, z0
    end

    for i in 1:n
        state = rand([1,2])
        if random_pos[i,2] > 0
            if random_pos[i,1] > 0
                rho_0, counts, z0 = add_agent(i, 1, state,  rho_0, counts, z0)
            else
                rho_0, counts, z0 = add_agent(i, 2, state,  rho_0, counts, z0)
            end
        else
            if random_pos[i,1] > 0
                rho_0, counts, z0 = add_agent(i, 3, state,  rho_0, counts, z0)
            else
                rho_0, counts, z0 = add_agent(i, 4, state,  rho_0, counts, z0)
            end
        end
    end

    rho_0 = rho_0/(sum(rho_0)*dV)
    u0 = rho_0
    y2_0 = [1.,1.]
    y1_0 = [-1.,-1.]
    y0 = [y1_0  y2_0]
    z0 = z0./dropdims(sum(counts, dims=2), dims=2)'
    controlled_inf = zeros(L)
    controlled_med = zeros(2)
    return ArrayPartition(u0,y0,z0),  controlled_inf,controlled_med
end


"""Right-hand side of the PDE"""

# interaction force at every grid point: (k * rho) / (w * rho)
function agent_force(rho, K_matrix, W_matrix,  dV)
    force = zeros(size(rho)..., 2)
    @views for d in 1:2
        f = vec(force[:,:,d])
        f .= dV * K_matrix[:,:,d] * vec(rho)
    end

    norms = dV * W_matrix * vec(rho)

    return force./reshape(norms, size(rho))
end

# recommender function r, bounded below by 0.1
function r(x)
    if x > 0.1
        return x
    else
        return 0.1
    end
end

# switching rates from influencer j to j2: eta * exp(-distance to j2) * r(...)
function gamma(grid_points, rho, z, eta,dV)
    Nx, Ny, M, L = size(rho)
    rate = zeros((Nx*Ny, M, L,L))
    m = dV*sum(rho, dims=(1,2))[1,1,:,:]
    dists = pairwise(Euclidean(), grid_points', z)

    for i in 1:M, j in 1:L, j2 in 1:L
        j == j2 && continue
        @views @. rate[:,i,j,j2] = eta * exp.(-dists[:, j2]) * r((m[i,j2]-m[mod1(i+1,Int(M)),j2])/(m[i,j2]+m[mod1(i+1,Int(M)),j2]))
    end

    return reshape(rate,(Nx, Ny, M, L,L))
end

# pull towards a medium or influencer at pos: pos - x at every grid point
function follower_force(pos, grid_points, N_x, N_y)
    pointwise_int = pos' .- grid_points
    force = reshape(pointwise_int,N_x,N_y, 2)
    return force
end

# in place Y = a*(A*B + B*A'), used for the diffusion term
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end

# dρ/dt for the 8 densities (diffusion - transport + switching) and dy/dt, dz/dt
# for media and influencers (relaxation towards their followers, or controlled motion)
function f(duyz,uyz,(p,q),t)
    yield()
    (; a, b, c, sigma, eta, switching, L, frictionM, frictionI, controlled_inf, controlled_med,controlspeed1,controltarget1, controlspeed2,controltarget2) = q
    (; grid_points, N_x, N_y,N, K_matrix, W_matrix, dx,dy, dV, C,  MassM) = p

    D = sigma^2 * 0.5
    u, y, z = uyz.x # density, media, influencers
    du, ydot, zdot = duyz.x

    rhosum = sum(u, dims=(3,4))[:,:,1,1]
    rhosum_j = sum(u, dims=3)
    rhosum_i = sum(u, dims=4)
    m_i = dV*sum(u,dims = (1,2,4))
    m_j = dV*sum(u,dims=(1,2,3))
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    rate_matrix = gamma(grid_points, u, z, eta, dV)
    reac = zeros(N_x, N_y)
    dif = similar(reac)
    dive = similar(reac)
    rhoforce = zeros(N_x, N_y, 2)
    for i in 1:2
        for j in 1:L
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            yi = @view y[:,i]
            dyi = @view ydot[:,i]
            zj = @view z[:,j]
            dzj = @view zdot[:,j]

            # divergence term
            rhoforce .= a .* Fagent .+ b .* follower_force(yi, grid_points, N_x, N_y) .+ c .* follower_force(zj, grid_points, N_x, N_y)
            rhoforce .= rho .* rhoforce
            @views mul!(dive, C, rhoforce[:,:,1]) # dive = C* rhoforce[:,:,1] (matrix-matrix product)
            @views mul!(dive, rhoforce[:,:,2], C', 1, 1) # dive = rhoforce[:,:,2] * C' + dive

            # reaction term
            reac .= 0
            for j2=1:L
                if j2 != j
                    @. @views reac += -rate_matrix[:,:,i,j,j2] * rho + rate_matrix[:,:,i,j2,j] * u[:,:,i,j2]
                end
            end

            # diffusion term
            a_AB_BAT!(dif, D, MassM, rho) # inplace dif = D*(M*rho + rho*M')

            # balance fluxes at boundary (part of boundary conditions)
            dif[1,:] += -1/dx * (rhoforce[1,:,1])
            dif[end,:] += 1/dx * (rhoforce[end,:,1])
            dif[:,1] += -1/dy * (rhoforce[:,1,2])
            dif[:,end] += 1/dy * (rhoforce[:,end,2])

            drho .=  dif .- dive
            if switching
                drho .+= reac
            end

            # note: missing boundary conditions for influencer and media movement
            if controlled_inf[j] == 0
                mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
                dzj .= 1/(frictionI) * (mean_rhoj' - zj)

            elseif controlled_inf[j] == 1 # first controlled influencer, movement defined by specified target and speed
                if norm(controltarget1' - zj) > 0
                    dzj .= controlspeed1* (controltarget1' - zj)./ norm(controltarget1' - zj)
                else
                    dzj .= 0
                end
            elseif controlled_inf[j] == 3 # second controlled influencer
                if norm(controltarget2' - zj) > 0
                    dzj .= controlspeed2* (controltarget2' - zj)./ norm(controltarget2' - zj)
                else
                    dzj .= 0
                end
            end
            if controlled_med[i] == 1 # controlled medium
                if norm(controltarget1' - yi) > 0
                    dyi .= controlspeed1* (controltarget1' - yi)./ norm(controltarget1' - yi)
                else
                    dyi .= 0
                end
            else
                mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
                dyi .= 1/(frictionM) * (mean_rhoi' - yi)
            end
        end
    end
end


"""Solve the PDE"""

# choose the initial condition from init and add the (zero) control flags to q
function constructinitial(init,(p,q))
    if init == "4inf"
        uyz0,  controlled_inf,controlled_med = uniforminit((p,q))
    elseif init == "closer"
        uyz0, controlled_inf,controlled_med = uniforminit((p,q); pos=0.6)
    elseif init == "2inf"
        uyz0, controlled_inf,controlled_med = uniforminit_2inf((p,q))
    elseif init == "uniform"
        uyz0, controlled_inf,controlled_med = uniforminit((p,q))
    elseif init == "random"
        uyz0, controlled_inf,controlled_med = randominit((p,q))
    end
    q = (; q..., controlled_inf,controlled_med)
    return uyz0, q
end

# build the initial condition, then solve the ODE system in time with a stiff solver
function PDEsolve(tmax=0.1; alg=nothing, init="4inf", p = PDEconstruct(), q= parameters())

    uyz0, q = constructinitial(init,(p,q))

    # Solve the ODE
    prob = ODEProblem(f,uyz0,(0.0,tmax),(p,q))
    @time sol = DifferentialEquations.solve(prob, alg, alg_hints = [:stiff], save_start=true, progress=true, progress_steps=1)

    return sol, (p,q)
end

# solution converter for plotting
function sol2uyz(sol, t)
    u = sol(t).x[1]
    y = sol(t).x[2]
    z = sol(t).x[3]
    return u,y,z
end

# solve, then plot snapshots at the times ts
function PDEsolveplot(; tmax=2.0, ts = [0. 0.1 0.4 1.0 1.5 2.0], alg=nothing, init="4inf", p = PDEconstruct(), q= parameters(),save=true, name=init)
    sol,(p,q) = PDEsolve(tmax; alg=alg, init=init, p=p, q=q)
    PDEplotsnapshots(sol, (p,q), ts;  name=name,save=save)

    return sol, (p,q)
end
