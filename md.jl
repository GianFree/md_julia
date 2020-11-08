# forces
"""Toy model, harmonic potential
"""
function ho( x, x_ref, k )
    return -k/2 * ( x - x_ref )^2
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
function verlet_step( x0, v0, F, Fargs; m = 1 )

    v_dt2 = v0 + 0.5 * F(x0, Fargs...)/m * dt
    x_dt  = x0 + v_dt2 * dt
    a_dt  = F(x_dt, Fargs...) / m
    v_dt  = v_dt2 + 0.5 * a_dt * dt

    #=
    v(t) -> v(t+dt/2) = v(t) + 1/2. a(t) * dt
    x(t) -> x(t+dt) = x(t) + v(t+dt/2)*dt
    a(t+dt) = F(t+dt)/m
    v(t+dt) = v(t+dt/2)+1/2 a(t+dt)dt
    =#
    return x_dt, v_dt
end


# particle position
xt = 0.
# particle velocity
vt = 5. #arb units
# Force variable
x_ref = 0
k = 4

positions = []
velocities = []
# Iteration
for ts in 1:10
    xt, vt = verlet_step(xt, vt, ho, (x_ref, k))
    append!( positions, xt )
    append!( velocities, vt )
# dumping the traj
