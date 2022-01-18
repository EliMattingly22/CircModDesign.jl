# using PyPlot

# using Interpolations

"""
Equations utilized are from Terman, 1943 Radio Engineers' Handbook

StrandDiam - Strand diameter (m) if one conductor is used, this is the conductor thickness

CoilTurns

CoilLength -Axial length of the coil (m)

CoilDiam - Diameter of the full solenoid (mean for elliptical ones)

NumStrands (default= 1) - Number of strands

BundleDiam (default= StrandDiam)

;
f = 1:1000:1e6, Test freuencies

PlotOn=true
"""
function LitzWire_ACRes(
        StrandDiam,
        CoilTurns,
        CoilLength,
        CoilDiam,
        NumStrands = 1,
        BundleDiam = StrandDiam;
        f = 1:1000:1e6,
        PlotOn=false
)

        if BundleDiam < (NumStrands * StrandDiam)
                error("Too many strands--will not fit if they are side-by-side")
        end

        ρ = 1.7e-8 #Ohm-meter
        WireL = pi * CoilDiam * CoilTurns
        CoilProps = Dict("WireL" => WireL)
        CoilProps["StrandDiam"] = StrandDiam
        CoilProps["BundleDiam"] = BundleDiam
        CoilProps["NumStrands"] = NumStrands
        StrandArea = π/4 *StrandDiam^ 2
        Rdc = ρ * WireL ./ (StrandArea) ./ NumStrands
        CoilProps["Rdc"] = Rdc

        u = Table20Terman(CoilTurns)

        μ₀ = 4*π*1e-7
        SkinDepth = sqrt.(ρ ./ (π * f * μ₀))
        Pitch = CoilLength / CoilTurns
        if Pitch < BundleDiam
                error("Pitch cannot be finer than a single bundle of wire")
        end
        Tmp = [[Table18Terman(StrandDiam,ρ,ff;units="meters")[1],
        Table18Terman(StrandDiam,ρ,ff;units="meters")[2],
        Table18Terman(StrandDiam,ρ,ff;units="meters")[3]] for ff in f]
        Tmp = hcat(Tmp...)
        X =Tmp[1,:]
        G = Tmp[2,:]
        H = Tmp[3,:]
        k = Table23Terman(NumStrands)
        Rac =[Rdc .* ( H[i] .+ (k + u .* (StrandDiam ./ Pitch) .^ 2) .*
                (StrandDiam/BundleDiam)^2 .* NumStrands^2 * G[i])
                for i in 1:length(f)]



        if PlotOn==true
                figure(1)
                PyPlot.plot(f, Rac,label="AC Resistance","k")

                PyPlot.plot([f[1], f[end]], [Rdc, Rdc] , "b--", label = "DC")
                xlabel("Freq. (Hz)")
                ylabel("Resistance (Ω)")
                legend()
        end

        return CoilProps, Rac, SkinDepth
end



function Table18Terman(X::AbstractFloat)
Table = [0.5   1       0.00097
         0.6   1.001   0.00202
         0.8   1.002   0.00632
         1.0   1.005   0.01519
         1.1   1.008   0.02196
         1.3   1.015   0.04127
         1.5   1.026   0.0691
         1.7   1.042   0.1055
         1.9   1.064   0.1489
         2.2   1.111   0.2214
         2.5   1.175   0.2949
         2.8   1.256   0.3412
         3.0   1.318   0.4049
         3.3   1.420   0.4626
         3.5   1.492   0.4987
         3.8   1.603   0.5503
         4.0   1.679   0.5842
         4.2   1.752   0.618
         4.5   1.863   0.669
         4.8   1.971   0.720
         5.0   2.043   0.755
         5.4   2.184   0.826
         6.0   2.394   0.932
         6.6   2.603   1.038
         7.0   2.743   1.109
         7.6   2.954   1.216
         8.0   3.094   1.287
         8.6   3.306   1.393
         9.0   3.446   1.464
         10.0  3.799   1.641
         13.0  4.856   2.171
         16.0  5.915   2.702
         20.0  7.328   3.409
         23.0  8.388   3.940
         25.0  9.094   4.294
         30.0  10.86   5.177
         40.0  14.40   6.946
         50.0  17.93   8.713
         60.0  21.46   10.48
         70.0  25.00   12.25
         80.0  28.54   14.02
         90.0  32.07   15.78
         100.  35.61   17.55]

         if X<0.5
                 G=X^4/64
                 H = 1.0
         elseif X>100
                 G = (sqrt(2)*X+1)/4
                 H = (sqrt(2)*X-1)/8
         else
                 interp_linearG = LinearInterpolation(Table[:,1], Table[:,3])
                 interp_linearH = LinearInterpolation(Table[:,1], Table[:,2])
                 G =interp_linearG(X)
                 H =interp_linearH(X)
         end


         return X,G,H
end

function Table18Terman(d,ρ,f;units="meters")
        if units=="meters"
                d = d*100 ##Units are cm
                ρ = ρ*1e11 ## units are abΩ⋅cm
        elseif (units=="centimeters")|(units=="cm")
        else
                error("Units not recognized")
        end
        X = π*d*sqrt(2*f/ρ)
        return Table18Terman(X)
end

"""
This function interpolates Table 20 in Terman Radio Engineers' Handbook
        It applies for solenoids with relatively low density (pitch > 2x wire diameter)
        It takes in number of turns, and returns `u`
"""
function Table20Terman(N)
Table = [2 1
         4 1.8
         6 2.16
         8 2.37
         10 2.51
         12 2.61
         16 2.74
         24 2.91
         32 3.0]
         if N>30
                 return 3.29
         else
                 interp_linear = LinearInterpolation(Table[:,1], Table[:,2])
                 return interp_linear(N)
         end
 end

function Table23Terman(n)
        Table = [3   1.55
                 9   1.84
                 27  1.92]
                 if n>27
                         return 2
                 elseif n<3
                         return 0
                 else
                         interp_linear = LinearInterpolation(Table[:,1], Table[:,2])
                         return interp_linear(n)
                 end

end
