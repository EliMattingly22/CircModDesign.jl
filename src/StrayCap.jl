using Optim
using MPI_Tools
##This uses equations presented in "Stray capacitances of single-layer air-core inductors for high-frequency applications"
# By"Grandi, G., Kazimierczuk, M. K.,Massarini, A. Reggiani, U."

function F(p,t,r,ϵr)
    return (p/(2*r))./((1+t/r)^(1-1/ϵr))
end
function Ct(D,r,p,t,ϵr)
    #D - solenoid diameter in meters
    #r - wire radius
    #p - pitch in meters
    #t - wire insulation thickness
    #ϵr - insulation relative permittivity

    ϵ0 = 8.85e-12

    FVal = F(p,t,r,ϵr)
    CtVal = pi^2*D*ϵ0 / (log(FVal+sqrt(FVal^2 - (1+t/r)^(2/ϵr))))
end


## Using derivation in appendix
function CC(r,t,ϵr)
    return 2*π*ϵ₀*ϵr / (log(1+t/r))
end

function CG(p,r,t)
    (π * ϵ₀) ./ log( (p/(2*r))/(1+t/r) + sqrt( ( (p/(2*r))/(1+t/r) )^2 - 1) )
end

"""
This function calculates the equivalent turn-turn capacitcance of a single layer solenoid
    #D - solenoid diameter in meters
    #r - wire radius
    #p - pitch in meters
    #t - wire insulation thickness
    #ϵr - insulation relative permittivity
"""
function Ceq(D,r,p,t,ϵr)
    CC(r,t,ϵr)*CG(p,r,t) / (CC(r,t,ϵr)+2*CG(p,r,t)) * π*D
end

# Ceq(.2,.15e-3,1.33e-3,.5e-3,4)/10


## Medhurst
function CMedh(L,D)
    ϵ₀ = 8.85*1e-12
    return (4*ϵ₀*L/π + 8e-12*D + 27e-12*D*sqrt(D/L)) #Farads
end

function CMedh_Knight(L,D)
    ϵ₀ = 8.85*1e-12
    return 4*ϵ₀*L/π*(1+ 0.8249* (D/L)+ 2.329*(D/L)^(3/2))#Farads
end

## Equation 5.3, Cₗ-DAE in Knight, p.52

function CL_DAE_Knight(
    L, #Coil length (axial)
    D, #Coil diamter
    N, #N turns
    ϵᵣₓ, #external permittivity
    ϵᵣᵢ) #Internal permittivity
    kc(L,D) =0.717439(D/L) +0.933048(D/L)^(3/2) +0.106(D/L)^2

    ϵ₀ = 8.85*1e-12
    p = L/N

    Ψ = atan(p/(π*D))

    C_L = (4*ϵ₀*ϵᵣₓ)/π *L * (1+kc(L,D)*(1+ϵᵣᵢ/ϵᵣₓ)/2)*1/(cos(Ψ)^2)
    return C_L
end
