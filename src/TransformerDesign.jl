using Optim
using Elliptic
using Struve
using SpecialFunctions

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

function BMax(V,f,N,A;k=(2*π/√2))
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


"""
Function to determine the leakage inductnace of a toroidal transformer
Ref. `Calculation of self and mutual impedances between sections of tranformer windings`, Wilcox et al. IEE Proceedings, 1989


"""
function CalcLeakageInductance_Toroid()

end
"""
Function to determine the leakage inductnace of a toroidal transformer
Ref. `Calculation of self and mutual impedances between sections of tranformer windings`, Wilcox et al. IEE Proceedings, 1989


"""
function CalcLₖₘ(Nₖ,Nₘ,r,a,z)
    μ₀ = 4*π*1e-7
    κ =κFunk(a,r,z)
    Lₖₘ = μ₀*Nₖ*Nₘ*√(r*a)*2/κ*(  (1 - κ^2/2) * Elliptic.K(κ) - Elliptic.E(κ) )

end



function Z₁₍ₖₘ₎(Nₖ,Nₘ,b,μᵣ,ω,λ,ρ)
    
    μ₁ = 4*π*1e-7 #for consistency with eqn. using 1 subscript
    μₓ =4*π*1e-7*μᵣ
    m = √(im*ω*μₓ/ρ)
    Z = 1*im*ω*Nₖ*Nₘ *π*b^2/λ * (2*μₓ*besseli(1,m*b) / (m*b*besseli(0,m*b)) - μ₁)
    

end

"""
Equation 12,Ref. `Calculation of self and mutual impedances between sections of tranformer windings`, Wilcox et al. IEE Proceedings, 1989
"""
function BigP₁(x,y,β)
    (p₁(x) - p₁(y)) ./ (β.^2)

end


function p₁(α)
    π*α/2 * (bessely(0,α)*Struve.L(0, α)+bessely(1,α)*Struve.L(1, α))
end


"""
Equation 10,Ref. `Calculation of self and mutual impedances between sections of tranformer windings`, Wilcox et al. IEE Proceedings, 1989
"""
function κFunk(a,r,z)
    √(4*a*r ./ (z^2 + (a+r)^2))
end

