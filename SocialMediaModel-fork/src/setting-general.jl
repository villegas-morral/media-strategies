function parameters(;
        L = 4,  # number of influencers
        M = 2, # number of media, warning: in some places the code assumes that M=2 and is not generalized
        n = 250, # number of agents
        eta = 15.0, # rate constant for changing influencer
        switching = false, # whether individuals switch influencer (at rate eta)
        a = 1., # interaction strength between agents
        b = 3., # interaction strength between agents and media
        c = 3., # interaction strength between agents and influencers
        sigma = 0.5, # noise on individual agents
        sigmahat = 0., # noise on influencers, warning: in some places in the code it is assumed that sigmahat=0
        sigmatilde = 0., # noise on media, warning: in some places in the code it is assumed that sigmatilde=0
        frictionI = 10., # friction for influencers
        frictionM = 100.,  # friction for media
        controlspeed1 = sqrt(1.5^2 +1.5^2)/5, # control definition of first influencer or medium
        controltarget1 = [1.5 1.5],
        controlspeed2 = sqrt(1.5^2 +1.5^2)/5, # control definition of another influencer
        controltarget2 = [1.5 1.5]
    )

    q = (; n, L, M, frictionM, frictionI, a, b, c, eta, switching, sigma, sigmahat, sigmatilde, controlspeed1,controltarget1, controlspeed2,controltarget2)
    return q
end


function parametersstronginf()
    q = parameters(;eta=1.)
    return q
end


function parameterscontrol(;ntarg=3, speedbound = 2, Tmax = 5,tequil = 5.,dtmin =0.001, start="influencer", maximize="follower",alpha=0.05,mindist=5.)
    r = (; ntarg, speedbound, Tmax,tequil,dtmin,start,maximize,alpha,mindist)
    # options for start: "zero", "mean", "influencer"
    # options for maximize:  "follower", "counter_follower"
    return r
end
