using Optim
using PyPlot
using FFTW
using ACME
using ElectricalEngineering
include("ToroidOptimizer.jl")


"""
This function takes in  the following parameters:

    LDrive, drive coil inductance (H)

    RDrive, Drive coil resistance (Ω)

    TargetZ, Amplifier/target impedance to match to  (Ω)

    DriveFreq; drive frequency (Hz)

    CDrive = 1e6, Any series capacitance with drive coil

    NumDriveElements = 1, The number of series copies of the L,C,R for example if there are 6 pancakes of LDrive, RDrive, Cdrive, then this would be 6

    WireDiam = 2e-3, Wire diameter for toroids

    WireFillFac = 0.75, If litz wire, this is the ratio of copper to air

    Example (Copy/paste):
    RDrive = 400e-3
    LDrive = 400e-6
    Amp_Z = 4
    DriveFreq = 25e3

    DesignDriveFilter(LDrive,RDrive,Amp_Z,DriveFreq)




"""
function DesignDriveFilter(
    LDrive,
    RDrive,
    TargetZ,
    DriveFreq;
    CDrive = 1e6,
    NumDriveElements = 1,
    WireDiam = 2e-3,
    WireFillFac = 0.75,
)
    ## This function takes in the drive coil impedances and designs a filter to match it to an amplifier
    #LDrive is the drive coil inductance in Henries
    #RDrive is the drive coil resistance (series) in Ohms
    #TargetZ is the ideal source impedance
    #DriveFreq is the drive frequency in hertz
    #CDrive is the SERIES capacitance, it defaults to 1e6 Farads (zero impedance at AC)
    #NumDriveElements is if there are multiple series repeated drive coil elements e.g. 2 drive assembles

    ZeroVal = 1e-12
    ωDr = 2 * pi * DriveFreq
    ZDrive =
        RDrive * NumDriveElements +
        im * ωDr * LDrive * NumDriveElements +
        NumDriveElements ./ (im * 2 * pi * DriveFreq * CDrive) #Equivalent complex impednace of the load

    Reactance_Load =
        ωDr * LDrive * NumDriveElements - NumDriveElements ./ (ωDr * CDrive)
    Reactance_Load = round(Reactance_Load; sigdigits = 3)
    println("The reactance of the load is: $Reactance_Load Ω")

    if (Reactance_Load > 0)

        println("Load is inductive")
        EquivSerL = Reactance_Load / ωDr
        println(
            "Load appears to be a inductor with a value of $(round(EquivSerL*1e6;sigdigits=3))μH and a ESR of $(real(ZDrive))Ω ",
        )

        matchRatio = real(TargetZ) / real(ZDrive)
        Q = sqrt(matchRatio - 1)
        Xs = Q * real(ZDrive) #Target reactance

        X_SerCap = Reactance_Load - Xs #Reactance of series capacitor
        SerCap = 1 / (ωDr * X_SerCap) #Series Capacitor value
        CParAct = Q / (ωDr * TargetZ)
        LTee_2 = ZeroVal
        LTee_2_ESR = ZeroVal
        LTee_1 = ZeroVal
        LTee_1_ESR = ZeroVal


    elseif (Reactance_Load < 0)
        println("Load is capacitive")
        EquivSerC = abs.(1 / (Reactance_Load * ωDr))
        println(
            "Load appears to be a capacitor with a value of $(round(EquivSerC*1e6;sigdigits=3))μF and a ESR of $(real(ZDrive))Ω ",
        )
        LSer = abs.(Reactance_Load) / (ωDr)
        matchRatio = real(TargetZ) / real(ZDrive)
        Q = sqrt(matchRatio - 1)
        Xs = Q * real(ZDrive)
        LSer2 = Xs / (ωDr)
        LTee_2 = LSer2 + LSer
        LTee_1 = 0
        LTee_1_ESR = 0
        CParAct = findResPair((1 + Q^(-2)) * LSer2, DriveFreq)
        LTee_2_Geom =
            ToroidOptimizer(WireDiam, LTee_2; CuFillFactor = WireFillFac)
        LTee_2_ESR = LTee_2_Geom.DCore.Resistance
        println("Q = $Q")
        SerCap = ZeroVal
    elseif (Reactance_Load == 0)
        matchRatio = real(TargetZ) / real(ZDrive)
        Q = sqrt(matchRatio - 1)
        Xs = Q * real(ZDrive)
        LSer2 = Xs / (ωDr)
        LTee_2 = LSer2
        LTee_1 = ZeroVal
        LTee_1_ESR = ZeroVal
        CParAct = findResPair((1 + Q^(-2)) * LSer2, DriveFreq)
        LTee_2_Geom =
            ToroidOptimizer(WireDiam, LTee_2; CuFillFactor = WireFillFac)
        LTee_2_ESR = LTee_2_Geom.DCore.Resistance
        SerCap = ZeroVal
    end
    println("L Tee 1 =  $(round(LTee_1*1e6;sigdigits=3))μH ")
    println("L Tee 2 =  $(round(LTee_2*1e6;sigdigits=3))μH ")
    println("CParAct =  $(round(CParAct*1e6;sigdigits=3))μF ")
    println("SerCap =  $(round(SerCap*1e6;sigdigits=3))μF ")


    (LFilt, CFilt) = Butterworth_2(TargetZ, DriveFreq)
    LFiltMatch_C = findResPair(LFilt, DriveFreq)
    CFiltMatch_L = findResPair(CFilt, DriveFreq)
    println(
        "LFilt =  $(round(LFilt*1e6;sigdigits=3))μH matched with: LFiltMatch_C =  $(round(LFiltMatch_C*1e6;sigdigits=3))μF ",
    )
    println(
        "CFilt =  $(round(CFilt*1e6;sigdigits=3))μF matched with: CFiltMatch_L = $(round(CFiltMatch_L*1e6;sigdigits=3))μH ",
    )
    circ = @circuit begin
        j_in = voltagesource()
        rs = resistor(TargetZ)
        ESR = resistor(RDrive)
        ESL = inductor(NumDriveElements * LDrive)
        ESC = capacitor(CDrive / NumDriveElements)
        MatchCSer = capacitor(SerCap)
        LTee_2_mod = inductor(LTee_2)
        LTee_2_ESRmod = resistor(LTee_2_ESR)
        LTee_1_mod = inductor(LTee_1)
        LTee_1_ESRmod = resistor(LTee_1_ESR)
        CParAct_mod = capacitor(CParAct)
        LFilt_mod = inductor(LFilt)
        CFilt_mod = capacitor(CFilt)
        CFiltMatch_L_mod = inductor(CFiltMatch_L)
        LFiltMatch_C_mod = capacitor(LFiltMatch_C)
        i_out = currentprobe()
        j_in[+] ⟷ rs[1] #\longleftrightarrow
        j_in[-] ⟷ gnd
        rs[2] ⟷ LFilt_mod[1]
        LFilt_mod[2] ⟷ LFiltMatch_C_mod[1]
        LFiltMatch_C_mod[2] ⟷ CFilt_mod[1]
        LFiltMatch_C_mod[2] ⟷ CFiltMatch_L_mod[1]
        CFiltMatch_L_mod[2] ⟷ gnd
        CFilt_mod[2] ⟷ gnd
        LFiltMatch_C_mod[2] ⟷ LTee_1_mod[1]
        LTee_1_mod[2] ⟷ LTee_1_ESRmod[1]
        LTee_1_ESRmod[2] ⟷ CParAct_mod[1]
        CParAct_mod[2] ⟷ gnd
        LTee_1_ESRmod[2] ⟷ LTee_2_mod[1]
        LTee_2_mod[2] ⟷ LTee_2_ESRmod[1]
        LTee_2_ESRmod[2] ⟷ ESC[1]
        ESC[2] ⟷ MatchCSer[1]
        MatchCSer[2] ⟷ ESL[1]
        ESL[2] ⟷ ESR[1]
        ESR[2] ⟷ i_out[+]
        i_out[-] ⟷ gnd
    end
    println("Starting to model")
    fs_model = 1e6
    model = DiscreteModel(circ, 1 / fs_model)
    FreqList = 1000:100:100e3
    NumPeriods = 100
    df = FreqList[2] - FreqList[1]
    Mags = zeros(length(FreqList), 1)
    WindowHanning(N) =
        0.5 .- 0.5 .* cos.(2 .* pi .* collect(0:(N-1)) ./ (N - 1))
    println("Starting to run sim")
    for F in FreqList
        y = run!(
            model,
            [
                2 .* sin(2π * F / fs_model * n) for c = 1:1,
                n = 1:NumPeriods*fs_model/F
            ],
        )
        y = y[:] .* WindowHanning(length(y))[:]
        FFT_y = abs.(fft(y)) / length(y) .* 4

        Mags[Int((F - FreqList[1]) / df)+1] = maximum(FFT_y)
    end

    semilogy(FreqList, Mags[:])
    xlabel("Frequency")
    ylabel("Drive Current")

    MaxCurrent = maximum(Mags[:])
    println("Max amp per volt is $(round(MaxCurrent,sigdigits=3))")
end

function Par(Z1, Z2)
    return 1 / (1 / Z1 + 1 / Z2)
end

"""
Converts from a polar coordinate (Mag, phase) to cartesian
if you enter a tuple or array it assumes first is mag.

"""
function Polar2Cart(Mag, Phase)
    Real = real.(Mag * exp(-im * Phase / 360 * 2 * pi))
    Imag = imag.(Mag * exp(-im * Phase / 360 * 2 * pi))
    return (Real, Imag)
end

function Polar2Cart(MagPhaseTup::Union{Tuple, Array})
    Polar2Cart(MagPhaseTup[1], MagPhaseTup[2])
end

function Z_Cap(C, f)
    return (1 ./ (im * 2 * pi .* f .* C))
end
function Z_Ind(L, f)
    return im * 2 * pi .* f .* L
end


"""
This function takes in some reactance element (E) and a frequency, and finds the pair


"""
function findResPair(E, f)
    ##Finds the resonant pair for a reactive element, E
    ## ωL = 1/(ωC) → L = 1/(ω^2C),and vice versa
    return 1 ./ ((2 * pi * f) .^ 2 * E)
end

function Butterworth_2(Z, f)
    #Sqrt(L/C) = Z ; f = 1/(2π*sqrt(LC))
    #L = Z^2*C
    #f = 1/(2π*Z C) → C = 1/(2*pi*f*Z)
    #Butterworth 2nd order = [LX = sqrt(2)Z, CX = sqrt(2)Z]
    C = 1 / (2 * pi * f) / Z * sqrt(2)
    L = 1 / (2 * pi * f) * Z * sqrt(2)
    return (L, C)
end




function PlotPhasor(
    Impedance::Tuple,
    PrevImpedance;
    Color = nothing,
    NumSteps = 1,
    Admittance = false,
    HeadWidth = 0.01,
)

    ## This function takes in an impedance, and plots it using 'phasor()'
    #optionally, it takes in a previous impedance which can be added in series or parallel
    #It then plots the line using PrevImpedance as a base
    #If Admittance is true, it means the inputs are Admittance not impedance, so parallel and series are switched
    NewImp = 0 + 0im
    dx = 0
    dy = 0
    for I = 1:NumSteps
        if I > 1
            #    Color = "gray"
        end
        if lowercase(Impedance[2]) == "series"
            if !Admittance
                NewImp = PrevImpedance + Impedance[1] / NumSteps
            else
                NewImp = Par(PrevImpedance, Impedance[1] * NumSteps)
            end
        elseif lowercase(Impedance[2]) == "parallel"
            if !Admittance
                NewImp = Par(PrevImpedance, Impedance[1] * NumSteps)
            else
                NewImp = PrevImpedance + Impedance[1] / NumSteps
            end
        else
            error("must be series or parallel")
        end

        RealParts = [real(NewImp), real(PrevImpedance)]
        ImagParts = [imag(NewImp), imag(PrevImpedance)]

        dx = -(RealParts[2] - RealParts[1])
        dy = -(ImagParts[2] - ImagParts[1])
        if Color != nothing



            PyPlot.plot(RealParts, ImagParts, color = "gray")

            if I == 1
                phasor(PrevImpedance, color = Color)
            end
        else
            PyPlot.plot(RealParts, ImagParts)
            if I == 1
                phasor(PrevImpedance)
            end
        end
        PrevImpedance = NewImp
    end
    # arrow(real(NewImp)-dx,imag(NewImp)-dy,dx,dy,head_width=HeadWidth,fc="gray",ec="k",alpha=1.,width = 0.00001,head_length= 2*HeadWidth)

    xlabel("Real")
    ylabel("Imag.")
    return NewImp
end


"""
This function takes in an array of tuples structures such as:

    [(<Complex impedance #1> , <"series" or "Parallel">),
    (<Complex impedance #2> , <"series" or "Parallel">)
    .
    .
    .
    (<Complex impedance #n> , <"series" or "Parallel">)]

    e.g.

    ZList = [(1+1*im, "series"),
             (1+5*im, "series"),
             (5+1*im, "parallel"),]
    PlotImpedanceTransformList(ZList)
"""
function PlotImpedanceTransformList(ImpList; InitialImp = nothing, ArrHeadWidth =1)
    ColorList = ["r", "g", "b", "magenta", "teal"]
    L = zeros(length(ImpList))
    for I = 1:length(ImpList)
        L[I] = abs(ImpList[I][1])
    end
    LongestVector = maximum(L)
    HeadWidth = LongestVector / ArrHeadWidth
    ##Impedance Plotting
    if InitialImp != nothing
        PrevImp = InitialImp
        StartVal = 1
    else
        PrevImp = ImpList[1][1]
        StartVal = 2
    end
    for I = StartVal:length(ImpList)
        PrevImp = PlotPhasor(
            ImpList[I],
            PrevImp;
            NumSteps = 500,
            Color = ColorList[I],
            HeadWidth = HeadWidth,
        )
    end
    phasor(PrevImp, color = ColorList[length(ImpList)+1])
    title("Impedance")
    ##Admittance Plotting
    figure()

    L = zeros(length(ImpList))
    for I = 1:length(ImpList)
        L[I] = abs(1 / ImpList[I][1])
    end
    LongestVector = maximum(L)
    HeadWidth = LongestVector / ArrHeadWidth
    if InitialImp != nothing
        PrevImp = 1 / InitialImp
        StartVal = 1
    else
        PrevImp = 1 / ImpList[1][1]
        StartVal = 2
    end
    for I = StartVal:length(ImpList)
        PrevImp = PlotPhasor(
            (1 / ImpList[I][1], ImpList[I][2]),
            PrevImp;
            NumSteps = 500,
            Color = ColorList[I],
            Admittance = true,
            HeadWidth = HeadWidth,
        )
    end
    phasor(PrevImp, color = ColorList[length(ImpList)+1])
    title("Admittance")
end
