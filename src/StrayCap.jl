using Optim

##This uses equations presented in "Stray capacitances of single-layer air-core inductors for high-frequency applications"
# By"Grandi, G., Kazimierczuk, M. K.,Massarini, A. Reggiani, U."

function F(p,t,r,ϵr)
    return (p/(2*r))./((1-t/r)^(1-1/ϵr))
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
