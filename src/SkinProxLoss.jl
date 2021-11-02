

"""
This function calculates the field components given two wire center position coordinates
    It assumes you want to fild the fields AT C₂ FROM C₁ as well as the (fully ANTIsymmetric) fields from the wire at C₂ on itself given some wire diameter WireD and current I₂ (assume 1 Amp)

    This will be some constant term, which will be ⟂ to the vector between the two(C₂-C₁), and have some gradient that decays in the direction of C₂-C₁. 

This is a 2D function and returns the field 
"""
function WireDist2Field(C₁::Vector,C₂::Vector,WireD=0.001,I₁=1,I₂=1)
    μ₀ = 4π*1e-7
    r = WireD/2
    DistVec = C₂ .- C₁
    Dist = √(sum(DistVec.^2))
    # D̂ = Dist ./ sum(Dist.^2) #the unit vector between the two coordinates
    # D′ = D̂ * [0 -1;1 0] #the unit vector rotated 90°

    MagH₂₁_Near = WireNearField_B(Dist-r,r;I=I₁)/μ₀ # The field from 1 on 2 is ⟂ to the vector between them
    MagH₂₁_Far = WireNearField_B(Dist+r,r;I=I₁)/μ₀ # The field from 1 on 2 is ⟂ to the vector between them
    GradH₂₁ = DerivField(Dist,r;I=I₁)/μ₀

    MagH₂₂ = WireNearField_B(r,r,I=I₂)/μ₀

    ## There are assumed to be four boundary conditions on wire 2, K₁,₂,₃,₄ where 1 is the point farthest from wire 1, 4 is the closest, and 2 is the point between 1,4 in the direction of D′ and 4 is the opposite of 2


    
    HK₁ =      MagH₂₂ + (MagH₂₁_Far)
    HK₃ = -1 * MagH₂₂ + (MagH₂₁_Near)

    HK₂ = -1 * MagH₂₂ 
    HK₄ =      MagH₂₂ 
    


    return HK₁,HK₂,HK₃,HK₄
end


"""
This function calculates the H field distribution in a Conductor (1D) where the width of the wire is W(meters),The skin depth is δ, the external field is ANTIsymmetric (skin effect) and directed in the Z direction (current is in Y) with a magnitude of K. X is the dimension that is ⟂ to I,H, and extends from -W/2, to W/2 

    Boundary conditions are ±K A/m at x=±W/2
"""
function HAntiSym(x,δ,W = 1., K=1.)
    if (x)>(W/2)
        return -K
    elseif (x)<(-1*W/2)
        return K
    else
        return  K*(exp.(-1*(1+im)*x/δ) - exp.((1+im)*x/δ)) / (exp.((1+im)*W/(2*δ)) - exp.(-1*(1+im)*W/(2*δ)))
    end

end

"""
This resultant current distribution that corresponds to HAntiSym(...). J = -dH/dx

With anti-symmetric fields, the current on opposite sides of the conductor is in phase(non-zero net current)
"""
function JAntiSym(x,δ,W = 1., K=1.)
    if abs(x)>(W/2)
        return 0
    else
        return  K*((1+im)/δ)*(exp.(-1*(1+im)*x/δ) + exp.((1+im)*x/δ)) / 
                      (exp.((1+im)*W/(2*δ)) - exp.(-1*(1+im)*W/(2*δ)))
    end
 
 end

 function HSym(x,δ,W = 1, K=1)
    if (x)>(W/2)
        return K
    elseif (x)<(-1*W/2)
        return K
    else
        return  K*(exp.(-1*(1+im)*x/δ) + exp.((1+im)*x/δ)) / (exp.((1+im)*W/(2*δ)) + exp.(-1*(1+im)*W/(2*δ)))
    end

end

"""
This resultant current distribution that corresponds to HSym(...). J = -dH/dx
With symmetric fields (even, as in proximity effect), the current on opposite sides of the conductor is in phase(zero net current)
"""
function JSym(x,δ,W = 1, K=1)
    if abs(x)>(W/2)
        return 0
    else
        return  K*((1+im)/δ)*(exp.(-1*(1+im)*x/δ) - exp.((1+im)*x/δ)) / 
                      (exp.((1+im)*W/(2*δ)) - exp.(-1*(1+im)*W/(2*δ)))
    end
 
 end

 
 function SkinDepth(f,ρ = 1.67e-8)
    μ₀ = 4π*1e-7
    SkinDepth = sqrt.(ρ ./ (π * f * μ₀))
 end

 

 function plot1DMagDiffusion(K₋,K₊,W,δ,n=1000)
    
    HTot, JTot = MagDiffusion1D(K₋,K₊,W,δ,1000)
    δX = W/n
    x = (-1*W/2):δX:(W/2)
    figure()
    pygui(true)

    subplot(221)
    plot(x,abs.(abs.(HTot)))
    title("|H|")
    subplot(223)
    plot(x,rad2deg.(angle.(HTot)))
    title("∠H")
    subplot(222)
    plot(x,abs.(JTot))
    plot(x,real.(JTot))
    
    subplot(224)
    plot(x,rad2deg.(angle.(JTot)))

    println(sum(JTot)/length(JTot)*pi/4*.001^2)
    return HTot, JTot
 end


 """
Calculates the proximity effect losses in one wire at C₁ due to many others in an Nx2 array, Cₐ
 """
 function MultiWireProximity(C₁::Vector,Cₐ::Array,f,WireD=0.001,I₁=1,Iₐ=1)


    
    Nₐ = length(Cₐ[:,1])
    Losses = 0.
    Losses += OneWireLoss(f,WireD,I₁)
    for i in 1:Nₐ
        Losses += sum(TwoWireLoss(Cₐ[i,:],C₁,f,WireD,1,0))*2
    end
    println(Losses)

    return Losses
 end


 """
Calculates the losses IN WIRE 2 FROM WIRE 1. By default only proximity losses are considered so the default is to have zero current in wire 2

 """
 function TwoWireLoss(C₁::Vector,C₂::Vector,f,WireD=0.001,I₁=1,I₂=0)

    δ = SkinDepth(f)
    HK₁,HK₂,HK₃,HK₄ = WireDist2Field(C₁,C₂,WireD,I₁,I₂)
    println("$([HK₁,HK₂,HK₃,HK₄])")
    HTot_LR, JTot_LR = MagDiffusion1D(HK₁,HK₃,WireD,δ,100000)
    HTot_UD, JTot_UD = MagDiffusion1D(HK₂,HK₄,WireD,δ,100000)
    J_LR² =  abs.(JTot_LR).^2
    J_UD² =  abs.(JTot_UD).^2
    # JTot_LR = WeightRadial!(JTot_LR)
    # JTot_UD = WeightRadial!(JTot_UD)
    σ=58e6
    LossLR = 1/2*(sum(J_LR²)/length(JTot_LR)*pi/4*WireD^2)/σ
    LossUD = 1/2*(sum(J_UD²)/length(JTot_UD)*pi/4*WireD^2)/σ
    # println(LossLR)
    # println(LossUD)
    # θ = LinRange(-π,π,  100)
    # TotLoss =√(sum( (sin.(θ)*LossUD).^2 .+ (cos.(θ)*LossLR).^2))/length(θ)
    return LossLR, LossUD

 end

 function OneWireLoss(f,WireD=0.001,I=1)

    Area = pi/4*WireD^2
    μ₀ = 4π*1e-7
    r = WireD/2
    MagH = WireNearField_B(r,r,I=I)/μ₀
    δ = SkinDepth(f)
    HTot, JTot = MagDiffusion1D(-MagH,MagH,WireD,δ,100000)
    JTot = 2 .* JTot ##account for X and Y
    σ=58e6
    # JTot = WeightRadial!(JTot)
    # J² =  WeightRadial!(abs.(JTot).^2)
    J²RMS =  WeightRadial!(abs.(JTot ./ (√2)).^2)
    print("Net current:")
    println(round( (sum(JTot)/length(JTot)*Area);sigdigits=3))
    Losses = (sum(J²RMS)/length(JTot)*Area)/σ 
    
    #println(Losses)

    return Losses
 end


function MagDiffusion1D(HK₋,HK₊,W,δ,n=1000)
    δX = W/n
    x = (-1*W/2):δX:(W/2)
    
     

    meanK = (HK₊+HK₋)
    ΔK = (HK₋-HK₊)/2
   
    HTot = HAntiSym.(x,δ,W, ΔK) .+ HSym.(x,δ,W, meanK)
    JTot = JAntiSym.(x,δ,W, ΔK) .+ JSym.(x,δ,W, meanK)

    return HTot, JTot
end


 
function SeriesWireBundleLosses(Cᵦ::Array, f,I,WireD)

    Nₐ = length(Cᵦ[:,1])
    Losses = 0.
    
    for i in 1:Nₐ
        Cₐ = Cᵦ[setdiff(1:end, i), :]
        C₁ = Cᵦ[i,:]

        Losses += MultiWireProximity(C₁,Cₐ,f,WireD,I,I)
    end
    return Losses
end

"""
This function takes in a vector and weights it such that the values are weighted by distance from the middle (for example when integrating a radial system)
"""
function WeightRadial!(V::Vector)
    Weights = abs.(collect(LinRange(-2,2,length(V))))
    Weights = Weights ./(sum(Weights))*length(V)
    V = V.*Weights
    return V
end


"""
This function takes in a vector and weights it such that the values are weighted by distance from the middle except with more weight towards the center. 
"""
function WeightCenter(V::Vector)
    Weights = abs.(abs.(collect(LinRange(-1,1,length(V)))) .- 2)
    # Weights = Weights ./(sum(Weights))*length(V)
    V = V.*Weights
    return V
end
