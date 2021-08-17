
using PyPlot
using Optim

function PlotNoise_VsR(en,in;YLims = [1e-10,1e-8])
    k_b = 1.38e-23
    R_L_List = 1:1:1e5
    noise = zeros(length(R_L_List),1)
    Thnoise = zeros(length(R_L_List),1)
    NoiseZ = en/in

    for i in 1:length(R_L_List)
        noise[i] = Noise(R_L_List[i],en,in)
        Thnoise[i] = ThermalNoise(R_L_List[i])
    end
    display(PyPlot.plot(R_L_List,noise[:],label="Total Noise"))
    PyPlot.plot(R_L_List,Thnoise[:],"r:",label="Thermal Noise Only")
    PyPlot.plot(NoiseZ,Noise(NoiseZ,en,in),"k^",label="Noise Matched Location",ms=10)
    ylim(YLims[1],YLims[2])
    legend()
    xlabel("Termination Resistance(Ω)")
    ylabel("Input Noise (V/√Hz)")
    xscale("log")
    yscale("log")
end

function PlotNoiseFigure_VsR(en,in;YLims = [1 10])
    k_b = 1.38e-23
    R_L_List = 1:1:1e5
    noise = zeros(length(R_L_List),1)
    NF = zeros(length(R_L_List),1)
    Thnoise = zeros(length(R_L_List),1)
    NoiseZ = en/in

    for i in 1:length(R_L_List)
        noise[i] = Noise(R_L_List[i],en,in)
        Thnoise[i] = ThermalNoise(R_L_List[i])

    end
    NF = 1 .+ noise[:] ./ Thnoise[:]

    display(PyPlot.plot(R_L_List,NF,label="Total Noise"))

    ylim(YLims[1],YLims[2])
    legend()
    xlabel("Termination Resistance(Ω)")
    ylabel("Noise Figure (NF)")

end

function DetermineNoise(RVals::Vector,MeasNoise::Vector)
    ## This function takes in a vector of measured noise (in V/√Hz) and the
    ## List of resistors that were used for termination, and estimates the preamp's en,in
    Error(NoiseParams) = 1e6*sum((MeasNoise - Noise(RVals,NoiseParams[1],NoiseParams[2])).^2)
    OptimResults = optimize(Error,[0.0,0.0],g_tol=1e-13)
    OptimResults.minimizer[1]
    PlotNoise_VsR(OptimResults.minimizer[1],OptimResults.minimizer[2])
    PyPlot.plot(RVals,MeasNoise,"ro",label = "Measured",ms=5)
    legend()
    return (OptimResults.minimizer)
end

function Noise(R_L,en,in;T=300)
    k_b = 1.38e-23
     sqrt.(4*k_b*T*R_L) .+ en .+ in*R_L
 end

function ThermalNoise(R_L)
    k_b = 1.38e-23
     sqrt.(4*k_b*300*R_L)
end
