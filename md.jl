using Plots

# forces
"""Toy model, FORCE due to harmonic potential
"""
function ho( x, x_ref, k )
    return -k * ( x - x_ref )
end

# velocity verlet
"""Velocity verlet function to update the position.
It uses the Trotter decomposition version of the
algorithm.

x0 = position at time t
v0 = velocity at time t
F  = function of the forces in play
Fargs = tuple with x_ref, k force constant
m  = 1 (default)

It returns:
xt
vt
"""
function verlet_step( x0, v0, F, Fargs; m = 1, dt=0.01 )

    v_dt2 = v0 + 0.5 * F(x0, Fargs...)/m * dt
    x_dt  = x0 + v_dt2 * dt
    a_dt  = F(x_dt, Fargs...) / m
    v_dt  = v_dt2 + 0.5 * a_dt * dt

    return x_dt, v_dt
end

function main()
    # particle position
    xt = 0.
    # particle velocity
    vt = 5 #arb units
    # Force variable
    x_ref = 0.
    k = 5

    positions = Float64[]
    velocities = Float64[]
    # Iteration
    for ts in 1:1000
        xt, vt = verlet_step( xt, vt, ho, (x_ref, k) )
        append!( positions, xt )
        append!( velocities, vt )
    end

    return positions, velocities
end

xx, vv = main()
plot(xx)


# dumping the traj
