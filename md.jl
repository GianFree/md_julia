using Plots
using LinearAlgebra

# forces
"""Toy model, FORCE due to harmonic potential
"""
function f_ho( x, x_ref, k )
    return -k .* ( x - x_ref )
end

"""Return """
function f_lj( x; ε=1, σ=1 )
    nparticles = size(x)[1]
    f_ij = zeros((nparticles, nparticles, 3))
    for i in 1:nparticles
        for j in i+1:nparticles
            Rij = x[i, :] - x[j, :]
            rij = norm(Rij)
            f_ij[i,j,:] = 48 * ε / rij^2 * ( (σ/rij)^12 - 0.5 * (σ/rij)^6 ) * Rij
        end
    end
    return f_ij
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

    v_dt2 = v0 + 0.5 * F.(x0, Fargs...)/m * dt
    x_dt  = x0 + v_dt2 * dt
    a_dt  = F.(x_dt, Fargs...) / m
    v_dt  = v_dt2 + 0.5 * a_dt * dt

    return x_dt, v_dt
end

function force_matrix( x, Fargs)
    matrix_force = ones((2,2,3))
    #return a matrix NxNx3
    return matrix_force
end

function force( x, f_matrix )
    return [ 1 0.0 0; -1 0. 0. ]
end

function update_step( x0, v0; m = 1, dt=0.01 )
    #x0 is a Nx3 matrix
    #v0 is a Nx3 matrix
    #Fargs is a tuple

    #F = force()

    # forza su particella 1:
    # risultante= np.sum(force[i, :])

    v_dt2 = v0 + 0.5 * force(x0)/m * dt
    x_dt  = x0 + v_dt2 * dt
    a_dt  = force(x_dt) / m
    v_dt  = v_dt2 + 0.5 * a_dt * dt

    return x_dt, v_dt
end

function main()
    # particle position
    xt    = [0.0 0.0 0.0; 0.0 0.0 0.0]
    # particle velocity
    vt    = [0.0 0.0 0.0; 0.0 0.0 0.0] #arb units
    # Force variable
    x_ref = [0.0, 0.0, 0.0]
    k     = [10, 10.0, 0.0]
    nparticles = 2

    positions  = xt
    velocities = vt

    println(f_lj([1 1 1; 3 3 3]))

    # Iteration
    io = open("traj.xyz", "w")
    for ts in 1:20000
        if (ts % 100 == 0)
            println(io, "$nparticles")
            println(io, "SCIAO GARI")
        end

        xt, vt = update_step( xt, vt )
        #for p in 1:nparticles
        #xt[p,:], vt[p,:] = verlet_step( xt[p,:], vt[p,:], f_ho, (x_ref, k) )
        #positions = cat(positions, xt, dims=2)
        #velocities = cat(velocities, vt, dims=2)
        if (ts % 100 == 0)
            for p in 1:nparticles
                println(io, "AR $(xt[p,1]) $(xt[p,2]) $(xt[p,3])")
            end
            #write(io, "AR $(xt[p,1]) $(xt[p,2]) $(xt[p,3])\n")
        end
    end

    close(io)

    #println( ho.( xt, x_ref, k) )
    return transpose(positions), transpose(velocities)
end

xx,vv = main()

#plot(xx[1],)
plot(xx[:,1], xx[:,2], )
# plot(xx[2])

#
# lj(r) = 4 * (-1/r^6 + 1/r^12)
# r = 0.95:0.01:5
# plot(r, lj.(r))

#println( ho.( [1,2,3], [0,0,0], 3  ) )
# dumping the traj
