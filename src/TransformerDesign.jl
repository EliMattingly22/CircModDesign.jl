using Plots
using Optim

function CalcFlux(Ae,Le,μr,N,V,ω)
    μ0 = 4*pi*1e-7
    MagnetizingInductance = N^2*μr*μ0*Ae /Le
    Impedance = ω*MagnetizingInductance
    I = V/Impedance
    Flux = N*I*μr*μ0*Ae /Le

    B = Flux/Ae
    return B
end
function CalcPermeability(AL,Ae,Le)
    #This function takes in the nominal inductance (inductance at one turn),effective length and area
    #And calculates the permeabily
    # L = N^2 μ_r*μ_o A_e /L_e ⟶ N=1, A_L = μA_e/L_e
    μ0 = 4*pi*1e-7
    return AL*Le/Ae/μ0
end

function BMax(V,f,N,A;k=4.44)
    return V ./ (k.*N.*A.*f)
end

function MaxTurns(CoreID,WireD)
    InnerWindingCircum = pi*(CoreID-WireD)
    return InnerWindingCircum/WireD
end

function RunMagnetizingFlux(WireD,CoreID,Ae,Le,AL,NRatio,VPk,ω)
    μ0 = 4*pi*1e-7
    μr = CalcPermeability(AL,Ae,Le)
    NTot  = MaxTurns(CoreID,WireD)
    NPrimary = NTot* (1/(1+NRatio))
    BCore = CalcFlux(Ae,Le,μr,NPrimary,VPk,ω)
    println("
    μr              = $(round(μr;sigdigits=3))
    Primary Turns   = $(round(NPrimary;sigdigits=3))
    Total Turns     = $(round(NTot;sigdigits=3))
    BCore(mT)       = $(round(BCore*1e3;sigdigits=3))
    ")
end


function I2(L1,L2,M,R1,R2)
    TransformerDiff = M*s / ((R1+L1*s)*(R2+L2*s) - M^2*s^2)
end

Vin(freq) = s/(s^2 + (2*pi*freq)^2)
Y = Vin(25e3)*I2(1e4,1e4,1e4,1,1)
data = simulate(Y,.001,nb_time_points=10000)
