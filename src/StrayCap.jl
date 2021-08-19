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

function Ceq(D,r,p,t,ϵr)
    CC(r,t,ϵr)*CG(p,r,t) / (CC(r,t,ϵr)+2*CG(p,r,t)) * π*D
end
