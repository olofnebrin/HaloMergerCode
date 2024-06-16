""" This program, based on extended Press-Schechter theory, generates a merger
    tree for a halo of mass M_final at redshift z_final. These parameters, as
    well as the mass resolution can be changed below. The output is a textfile
    containing the masses and redshifts of all progenitors.

    ALGORITHM BY: H. Parkinson et al. (2008), Generating Dark Matter Halo
                  Merger Trees, MNRAS, vol. 383, 2, pp. 557-564.

    CODE AUTHOR:  Olof Nebrin
    DATE:         July, 2019 (updated Dec, 2021)  """

module MergerTree

# EXPORTED FUNCTIONS:

export M_main                    # This is our main function for generating merger trees, so
                                 # this will be exported for use in other files.

export FindStarFormingHalos      # For finding atomic-cooling halos.


# IMPORTED PACKAGES:

using QuadGK                     # To evaluate a needed integral
using Dierckx                    # For interpolation of the needed integral (to speed up code)
using Printf                     # For printing info
using DelimitedFiles             # For saving data in text files

z_reion    = 6.0
v_reion    = 30.0

###############################################################
#######                                                 #######
#######               C O N S T A N T S                 #######
#######                                                 #######
###############################################################

# Physical & astronomical constants (in cgs units)

G      = 6.673e-8     # Gravitational constant
Msun   = 1.989e33     # Solar mass
pc     = 3.086e18     # Parsec
mpc    = 1e6*pc       # Mpc
yr     = 3.1556926e7  # 1 year in seconds

# Cosmological parameters:

h      = 0.674        # Normalized Hubble constant
H0     = h*100e5/mpc  # Hubble constant (in s^-1)
Ω_M0   = 0.315        # Present matter density parameter
Ω_Λ0   = 0.685        # Present dark energy density parameter
Ω_B0   = 0.049        # Present baryon density parameter
δ_crit = 1.686        # Critical linear overdensity for collapse
σ8_obs = 0.811        # Observed RMS over density at R = 8.0/h Mpc

# Recovery time-scale of gas (to get the threshold halo mass for
# Pop II star formation)

t_recov = 100e6*yr

# Parameters related to the conditional halo mass function:

G0 = 0.57
γ1 = 0.38
γ2 = - 0.01

# Parameters related to the redshift step size in the algorithm
# (they should be << 1):

ϵ1 = 0.1
ϵ2 = 0.1

###############################################################
#######                                                 #######
#######               F U N C T I O N S                 #######
#######                                                 #######
###############################################################

function σ(M)

    """ RMS mass fluctuation as a function of halo mass, evaluated
        at the present (z = 0). This is a convenient yet fairly
        accurate analytical approximation. """

    M_eq = 2.4e17*((Ω_M0*(h^2)/0.14)^(-1/2))
    m    = 8*M/M_eq
    N    = 0.0845*σ8_obs

    σ = N*sqrt( (36/(1 + 3*m)) - log(m/(1 + m))^3 )

end

function Δ_vir(z)

    """ The virial overdensity ρ_vir/ρ_M at redshift z.
        Uses the accurate fit found in Eq. (A23) in:

        M. Tegmark et al. (2006), Dimensionless constants,
        cosmology, and other dark matters, Phys. Rev. D, vol. 73, 2. """

    x     = (Ω_Λ0/Ω_M0)*(1/(1 + z)^3)

    Δ_vir = 18*(pi^2) + 52.8*(x^0.7) + 16*x

end

function R_vir(M,z)

    """ The virial radius of a halo of mass M
        at redshift z (cm). """

    # The mean matter density:

    H0  = 3.24e-18*h
    ρ_M = Ω_M0*(3*(H0^2)/(8*pi*G))*((1+z)^3)

    # The virial radius:

    R_vir = (3*M*Msun/( 4*pi*Δ_vir(z)*ρ_M ))^(1/3)

end

function v_vir(M,z)

    """ The virial velocity (in km/s) of of a halo of mass M
        at redshift z. """

    # Put it all together:

    v_vir = ( (G*M*Msun/R_vir(M,z))^0.5 )/1e5

end

function α(M)

    """ The function α(M) = - dlnσ/dlnM. """

    # Step size:

    δM = 1e-6*M

    # The derivative dσ/dM:

    dσ_dM = (1/(2*δM))*( σ(M + δM) - σ(M - δM) )

    # Put it all together:

    α = - (M/σ(M))*dσ_dM

end

# Some quantities derived from the above functions:

α2(M2)   = α(M2)
αh(M2)   = α(M2/2)
α1(q,M2) = α(q*M2)

σ2(M2)   = σ(M2)
σh(M2)   = σ(M2/2)
σ1(q,M2) = σ(q*M2)

V(q,M2)     = (σ1(q,M2)^2)/( (σ1(q,M2)^2 - σ2(M2)^2)^(3/2) )

β(q_res,M2) = log( V(q_res,M2)/V(1/2, M2) )/log( 2*q_res )
B(q_res,M2) = V(q_res,M2)/(q_res^β(q_res,M2))

μ(M2)       = αh(M2)
η(q_res,M2) = β(q_res,M2) - 1 - γ1*μ(M2)

# Other needed functions:

function GrowthFactor(z)

    """ The unnormalized growth factor. """

    Ω_tot = Ω_Λ0 + Ω_M0*((1+z)^3)
    Ω_M   = Ω_M0*((1+z)^3)/Ω_tot
    Ω_Λ   = Ω_Λ0/Ω_tot

    FirstFactor  = 5*Ω_M/(2*(1+z))
    SecondFactor = Ω_M^(4/7) - Ω_Λ + (1 + 0.5*Ω_M)*(1 + Ω_Λ/70)

    GrowthFactor = FirstFactor/SecondFactor

end

# Normalize the growth factor:

D(z) = GrowthFactor(z)/GrowthFactor(0.0)

function δ(z)

    """ Simply the function δ_crit/D(z). """

    δ = δ_crit/D(z)

end

function dδ_dz(z)

    """ The derivative of δ(z). """

    # Step size:

    δz = 1e-5

    # The derivative:

    dδ_dz = ( δ(z + δz) - δ(z - δz) )/(2*δz)

end

function S(q,z,M2,q_res)

    """ The function S(q) needed for dN/dq. """

    # Convenient variables:

    F1 = q^(η(q_res,M2)-1)
    F2 = ( δ(z)/σ2(M2) )^γ2
    F3 = ( σh(M2)/σ2(M2) )^γ1

    # Put it all together:

    S = sqrt(2/pi)*B(q_res,M2)*αh(M2)*F1*(G0/2^(μ(M2)*γ1))*F2*F3*dδ_dz(z)

end

function R(q, M2, q_res)

    """ The function R(q), also needed for dN/dq. """

    R = ( ( α1(q,M2)/αh(M2) )*( V(q,M2)/(B(q_res,M2)*(q^β(q_res,M2))) )
        *( (((2*q)^μ(M2))*σ1(q,M2)/σh(M2))^γ1 ) )

end

function Δz(z, M2, q_res)

    """ The redshift step size for the merger
        tree algorithm. """

    # First we will need the integral
    # of S(q) from q = q_res to q = 1/2.
    # This integral can be evaluated analytically
    # to yield:

    Integral = (S(1,z,M2,q_res)/η(q_res,M2))*(
                0.5^η(q_res,M2) - q_res^η(q_res,M2) )

    # We also need the following:
    dz = sqrt(2)*sqrt( σh(M2)^2 - σ2(M2)^2 )/dδ_dz(z)

    # The minimum redshift interval:

    Δz = min( ϵ1*dz, ϵ2/Integral )

end

function N_upper(z, M2, q_res)

    """ The upper limit to the expected number of
        resolved fragments produced in a redshift
        step Δz.  """

    N_upper = (S(1,z,M2,q_res)/η(q_res,M2))*(
               0.5^η(q_res,M2) - q_res^η(q_res,M2) )*Δz(z,M2,q_res)

end

function J_exact(u_res)

    """ Needed integral for the calculation
        of F. """

    Integrand(u) = (1 + 1 ./u.^2).^(γ1/2)

    Integral = quadgk( Integrand, 0, u_res )[1]

end

# Using J_exact for J leads to unnecessarily long computation time.
# We can speed up the algorithm by a factor of ~ 4 by making an
# interpolating approximation of J_exact:

println("")
println(" C R E A T I N G  I N T E R P O L A T I O N  F U N C T I O N ...")
println("")

x     = range( 1e-5, 300.0, length = 600 )
data  = Float64[ J_exact(el) for el in x ]
J     = Spline1D(x, data)

println(" D O N E !")
println("")

function F(z, M2, q_res)

    """ The fraction of mass in progenitors to M2
        below the resolution limit. """

    M_res = q_res*M2
    σres  = σ(M_res)
    u_res = σ2(M2)/sqrt( σres^2 - σ2(M2)^2 )
    δ2    = δ(z)

    F = sqrt(2/pi)*J(u_res)*(G0/σ2(M2))*((δ2/σ2(M2))^γ2)*dδ_dz(z)*Δz(z,M2,q_res)

end

function M_main(M0, z0, z_max, M_res, DataFileName1, DataFileName2,
                print_info = true, save = true)

    """ The algorithm to find the evolution of the
        main progenitor, given that it has a mass M0
        at redshift z0. The main branch is followed
        back to a maximum redshift z_max.

        RETURNS: M_main, z_main, v_vir, M_sub, z_sub

        where M_main[i] is the mass of the main branch
        at redshift z_main[i], and v_vir[i] is the corresponding
        virial velocity at that redshift. Furthermore,
        M_sub[i] is the mass of an accreted subhalo at redshift
        z_sub[i].

        Furthermore, the data is saved in data files with names
        DataFile1 and DataFile2 containing the
        following data:

        DataFileName1: M_main, z_main, v_vir
        DataFileName2: M_sub,  z_sub """

    # First we initialize the arrays we need:

    M_main = Float64[M0]
    z_main = Float64[z0]
    V_vir  = Float64[v_vir(M0,z0)]
    M_sub  = Float64[]
    z_sub  = Float64[]

    ###############################################################
    #######                                                 #######
    #######               A L G O R I T H M                 #######
    #######                                                 #######
    ###############################################################

    if print_info == true

        println("")
        println(" S T A R T I N G  I T E R A T I O N ... ")
        println("")

        println("Main progenitor mass M (M_⊙)   Redshift (z)")
        println("")

    end

    while z_main[end] < z_max && M_main[end] > 2*M_res

        # Set the mass of the main progenitor:

        M2 = M_main[end]

        if print_info == true

            # Print the values of the main progenitor mass
            # and the corresponding redshift (with 3
            # significant figures):

            info = @sprintf( "%.2E %.2E", M2, z_main[end] )
            println(info)

        end

        # Determine q_res:

        q_res = M_res/M2

        # Draw the first of three random numbers, uniformally
        # between 0 and 1:

        r1 = rand()

        if r1 > N_upper(z_main[end], M2, q_res)

            # In this case no split occurs, but the halo
            # mass M2 is reduced to M2(1 - F):

            M_new  = M2*( 1 - F(z_main[end], M2, q_res) )

            # Save (i.e. append) the new mass to M_main list,
            # and similarly the virial velocity:

            push!(M_main, M_new)

            # We also calculate and add the next redshift:

            z_new  = z_main[end] + Δz(z_main[end], M2, q_res)
            push!( z_main, z_new )

            # Finally, also save the virial velocity at z_new:

            push!(V_vir, v_vir(M_new,z_new))


        else

            # In this case we generate a second random variable:

            r2 = rand()

            # We transform it to give a value for q:

            q = (q_res^η(q_res,M2) +
                 ( 2^(-η(q_res,M2)) - q_res^η(q_res,M2) )*r2)^( 1/η(q_res,M2) )

            # We draw yet another random value:

            r3 = rand()

            if r3 < R(q, M2, q_res)

                # In this case we generate we get
                # two progenitors of mass q*M2 and
                # M2*(1 - F - q):

                Mdraw1 = M2*q
                Mdraw2 = M2*( 1 - F(z_main[end], M2, q_res) - q )

                # The larger of the two will be
                # assigned as the next main progenitor,
                # with the other being labelled as a
                # subhalo:

                Mlarger  = max(Mdraw1, Mdraw2)
                Msmaller = min(Mdraw1, Mdraw2)

                push!( M_main, Mlarger )
                push!( M_sub,  Msmaller )

                # Also calculate and save the needed
                # redshifts:

                push!( z_sub, z_main[end] )

                z_new  = z_main[end] + Δz(z_main[end], M2, q_res)
                push!( z_main, z_new )

                # Finally, also save the virial velocity at z_new:

                push!(V_vir, v_vir(Mlarger,z_new))

            else

                # Otherwise there is no split,
                # and we simply reduce the progenitor mass:

                M_new  = M2*( 1 - F(z_main[end], M2, q_res) )

                # Save (i.e. append) the new mass to M_main list:

                push!(M_main, M_new)

                # We also calculate and add the next redshift:

                z_new  = z_main[end] + Δz(z_main[end], M2, q_res)
                push!( z_main, z_new )

                # Finally, also save the virial velocity at z_new:

                push!(V_vir, v_vir(M_new,z_new))

            end

        end

    end

    if print_info == true

        println("")
        println(" I T E R A T I O N  C O M P L E T E ")
        println("")

    end

    if save == true

        # Next we save the data in text files. First we get the
        # directory of the current script (i.e. MergerTree.jl):

        CurrentDir = string( @__DIR__ )

        # We then get the following file paths:

        FilePath1  = string( CurrentDir, "\\", DataFileName1 )
        FilePath2  = string( CurrentDir, "\\", DataFileName2 )

        writedlm( FilePath1, [M_main z_main V_vir] )
        writedlm( FilePath2, [M_sub z_sub] )

    end

    # When done with the iteration, we can return all the
    # quantities of interest:

    M_main, z_main, V_vir, M_sub, z_sub

end


#####################################################################
#######                                                       #######
#######      F I N D  S T A R - F O R M I N G  H A L O S      #######
#######                                                       #######
#####################################################################

# Here the algorithm is developed for finding the number
# of star-forming halos accreted onto the host halo, as well
# the redshift when they reach the cooling threshold.

# First we use the M_main routine to find all the subhalos
# of interest#

function v_thresh(z)

    """ The cooling threshold expressed
        in terms of the virial velocity. For
        redshifts z < z_reion, it is set to v_reion. """

    if z > z_reion

        # Atomic-cooling threshold virial velocity prior
        # to reionization:

        Tvir   = 9946.0*(((1+z)/10)^(-0.052))
        v_cool = 13.5*sqrt(Tvir/1e4)

    else

        # After reionization:

        v_cool = v_reion
    end

    return v_cool

end

function d_host(z, z_acc, Mhost)

    """ This function computes the distance (in kpc) between an accreted
        subhalo and its host (of mass Mhost(z_acc)) at earlier
        redshifts z > z_acc (i.e. prior to accretion at redshift z_acc).
        Uses a fit to the spherical collapse model. """

    # Redshift-time relation (ignoring normalization because it
    # cancels in the ratio t/t_vir):

    t     = (1/(1+z))^(3/2)
    tacc  = (1/(1+z_acc))^(3/2)

    # Relation between t_vir and accretion time t_acc:

    tvir  = tacc/0.909

    # Distance from the host halo:

    stuff  = 1 - 0.237*(sin(pi*t/tvir)^(3/5))
    d_host = 0.5*((12*sin(pi*t/tvir))^(2/3))*stuff*R_vir(Mhost,z_acc)

    return d_host/(1e3*pc)

end

function M_th(z, J_21)

    """ Cooling threshold at redshift z given a local Lyman-Werner background
        J_21. """

    # Convenient definition:

    Z = (1+z)/10.0

    # H2-cooling threshold in the absence of a Lyman-Werner background
    # (derived in the thesis):

    M_th_H2_noLW  = 1.4e6*exp(0.015*(Z^0.2707))*(Z^(-1.094))/(log(1
                    + 1.8*(Z^(-2.25)))^0.2703)

    # H2-cooling threshold in a Lyman-Werner background (Latif & Khochfar, 2019),
    # also including a uniform background J_21_bg:

    J_21_bg       = 10*exp(-(z-6.0)/3.3)

    J_21_tot      = J_21 + J_21_bg

    M_th_LW       = exp(19.38 - 4.3*exp(-0.16*log(J_21_tot)))

    # Resulting H2-cooling threshold:

    M_th_H2       = max( M_th_H2_noLW, M_th_LW )

    # Atomic-cooling threshold due to Lyman-α cooling
    # (derived in the thesis):

    Tvir4    = 0.9946*(Z^(-0.052))
    M_th_Lyα = 5.11e7*(Tvir4^(3/2))*(Z^(-3/2))

    # Putting it all together:

    M_th     = min( M_th_Lyα, M_th_H2 )

    return M_th

end

function FindStarFormingHalos( M0, z0, z_max, M_res, DataFileNameHalos,
         savedata = true )

    """ This routine finds all 1st-order subhalos around a host of mass
        M0 at redshift z0 that have reached a given threshold for
        star formation (v_vir = v_thresh), and when their main
        progenitor reached this threshold. """

    # First we use the M_main routine to find all the subhalos
    # of interest:

    HostAndSubHalos = M_main(M0, z0, z_max, M_res,
                      "Mmain_zmain_vvir.txt", "Msub_zsub.txt", true, true)

    SubHaloMasses   = HostAndSubHalos[4]
    z_sub           = HostAndSubHalos[5]

    # Also make an interpolation for the mass of the host progenitor as a function
    # of redshift. This will be needed to compute the Lyman-Werner
    # feedback from the host progenitor.

    HostProgenitorMassData  = HostAndSubHalos[1]
    HostRedshiftData        = HostAndSubHalos[2]

    HostProgenitorMass      = Spline1D(HostRedshiftData, HostProgenitorMassData)

    # Next we use M_main again, but for each of the candidate
    # subhalos of interest. We do this in order to get their
    # mass assembly history, and to find the redshift when
    # they reach the star-forming threshold threshold.

    z_PopII     = Float64[]  # Stores redshift when gas is recovered and Pop II stars form
    M_threshold = Float64[]  # Stores halo mass at z_PopII
    M_acc       = Float64[]  # Stores the halo mass when accreted onto the host
    z_acc       = Float64[]  # Stores accretion redshift
    v_peak      = Float64[]  # Stores the peak virial velocity

    for j in 1:length( SubHaloMasses )

        Msubhalo       = SubHaloMasses[j]
        zsubhalo_acc   = z_sub[j]

        SubHaloHistory = M_main(Msubhalo, zsubhalo_acc , z_max, M_res,
                         "HostMassRedshiftVirialVelocity.txt", "placeholder.txt", false, true)
        MHistory       = SubHaloHistory[1]
        zHistory       = SubHaloHistory[2]
        vHistory       = SubHaloHistory[3]

        # Compute the host progenitor mass at each redshift point
        # of the subhalo. This is needed for calculating Lyman-Werner
        # feedback:

        HostMass       = HostProgenitorMass.( zHistory )

        # Host halo mass at the accretion redshift (needed for the
        # distance from the host as a function of redshift):

        HostMassAcc    = HostProgenitorMass(zsubhalo_acc )

        # Calculate the Lyman-Werner background at the position
        # of the subhalo at a given redshift in its history:

        ZHistory     = (1 .+ zHistory)/10.0
        rkpc         = d_host.(zHistory, zsubhalo_acc, HostMassAcc)

        J_21_gal     = 400.0*((HostMass/1e10).^1.58).*(ZHistory.^2.20).*(rkpc.^(-2))

        # We want to find the redshift z_cool where the subhalo
        # reaches the cooling threshold M = M_th.
        # To find z_cool we can find the redshift z for which
        # abs( M(z) - M_th(z) )/M_th(z) is minimized:

        M_thresh_z     = M_th.( zHistory, J_21_gal )
        difference     = abs.( MHistory - M_thresh_z )
        index_min      = argmin( difference )

        # Only save the redshift z_cool and the corresponding halo
        # mass if it indeed is close to the cooling threshold:

        upperbound_index = min( index_min + 3, length(MHistory) )
        lowerbound_index = max( index_min - 3, 1 )

        if MHistory[upperbound_index] < M_th( zHistory[index_min], J_21_gal[index_min] )
            if MHistory[lowerbound_index] > M_th( zHistory[index_min], J_21_gal[index_min] )

                # Redshift when the gas first cool and form Pop III stars:

                z_Cool   = zHistory[index_min]

                # Next we find the redshift when the gas is recovered:

                t_Cool   = (2/(3*H0*sqrt(Ω_M0)))*1/((1+z_Cool)^(3/2))

                z_recov  = (2/(3*H0*sqrt(Ω_M0)*(t_Cool + t_recov)))^(2/3) - 1

                if z_recov > z0

                    # Find the redshift closest to z_recov and the corresponding halo mass:

                    index_recov_min = argmin( abs.(zHistory.-z_recov) )

                    z_Recov         = zHistory[index_recov_min]
                    M_Recov         = MHistory[index_recov_min]

                    push!( z_PopII, z_Recov )
                    push!( M_threshold, M_Recov )
                    push!( M_acc, Msubhalo )
                    push!( z_acc, zsubhalo_acc  )

                    println("Threshold halo mass:")
                    println(M_Recov)

                    # We also save the peak virial velocity:

                    push!( v_peak, findmax(vHistory)[1] )

                end

            end
        end

    end

    # When finished, print the number of star-forming halos
    # and save the data in textfiles:

    if savedata == true
        println("")
        println("NUMBER OF STAR-FORMING SUBHALOS: ", length(z_PopII) )
        println("")

        CurrentDir = string( @__DIR__ )
        FilePath   = string( CurrentDir, "\\", DataFileNameHalos )
        writedlm( FilePath, [z_PopII z_acc M_threshold M_acc v_peak] )
    end

    # Return the data:

    z_PopII, M_threshold

end

end
