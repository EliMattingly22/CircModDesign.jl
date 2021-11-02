module CircModDesign
using LinearAlgebra
using PyPlot
using MPI_Tools
using Interpolations
using Optim
using Elliptic
using Struve
using SpecialFunctions


using FFTW
using ACME
using ElectricalEngineering

include("Filter_Designer.jl")
include("NoiseModeling.jl")
include("ToroidOptimizer.jl")
include("TransformerDesign.jl")
include("StrayCap.jl")
include("LitzACResist.jl")
include("InductanceCalc.jl")
include("MakeSolenoid.jl")
include("RunAnalysis.jl")
include("SPICE2Matrix.jl")
include("CoupledMagnetics.jl")
include("ThermalModeling.jl")
include("Toroid_Inductance.jl")
include("SkinProxLoss.jl")


export  PlotImpedanceTransformList,
        BiotSav,
        MakeEllip,
        FieldOnAxis_Circ,
        DesignDriveFilter,
        Par,
        Polar2Cart,
        Z_Cap,
        Z_Ind,
        findResPair,
        Butterworth_2,
        CalcInductance,
        SimpleInduct,
        LitzWire_ACRes,
        MakeEllipSolenoid,
        EvalInduct_Biot,
        EvalField_Centered,
        Noise,
        ThermalNoise,
        CL_DAE_Knight,
        ToroidOptimizer,
        Length2Resist,
        FieldMapPointPath,
        CompareD_Circ_Rect_Toroid,
        GetElementVal,
        UpdateElementVal!,
        UpdateElementESR!,
        SPICE_DF2Matrix_ω,
        PipeFlow,
        OneWireLoss,
        TwoWireLoss,
        MultiWireProximity,
        SeriesWireBundleLosses,
        plot1DMagDiffusion




# Write your package code here.

end
