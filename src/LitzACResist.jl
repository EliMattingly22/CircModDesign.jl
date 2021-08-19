using PyPlot

function LitzWire_ACRes_CalcFunk(
        StrandDiam,
        CoilTurns,
        CoilLength,
        CoilDiam,
        NumLayers = 1,
        NumStrands = 1,
        BundleDiam = StrandDiam;
        f = 1:1000:1e6,
)

        if BundleDiam < (NumStrands * StrandDiam)
                error("Too many strands--will not fit if they are side-by-side")
        end

        ρ = 1.7e-8 #Ohm-meter
        WireL = pi * CoilDiam * CoilTurns * NumLayers
        CoilProps = Dict("WireL" => WireL)
        CoilProps["StrandDiam"] = StrandDiam
        CoilProps["BundleDiam"] = BundleDiam
        CoilProps["NumStrands"] = NumStrands

        Rdc = 4 * ρ * WireL ./ (pi * StrandDiam .^ 2) ./ NumStrands
        CoilProps["Rdc"] = Rdc

        SkinDepth = sqrt.(ρ ./ (pi * f * 4 * pi * 1e-7))

        G = (StrandDiam ./ (4 * SkinDepth)) .^ 4
        F = G / 3

        Pitch = 1 / (CoilTurns ./ CoilLength)
        if Pitch < BundleDiam
                error("Pitch cannot be finer than a single bundle of wire")
        end
        k = 2 - 1.4 ./ NumStrands
        u = 1 ./ (0.18 .* exp.(-0.68 .* CoilLength ./ CoilDiam) .+ 0.11)
        Rac_Straight =
                Rdc .* (
                        1 .+ F .+
                        k .* (NumStrands .* StrandDiam ./ BundleDiam) .^ 2 .* G
                )
        Rac_Coiled =
                Rdc .* (
                        1 .+ F .+
                        k .* (NumStrands .* StrandDiam ./ BundleDiam) .^ 2 .*
                        G .+ u .* (NumStrands .* StrandDiam ./ Pitch) .^ 2 .* G
                )

        figure(1)
        semilogy(f, Rac_Coiled)

        plot([f[1], f[end]], [Rdc, Rdc])

        return CoilProps, Rac_Straight, Rac_Coiled, SkinDepth
end
