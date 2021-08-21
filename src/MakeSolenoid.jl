using LinearAlgebra
function MakeEllipSolenoid(
    N, #Number of turns
    r₁, #Radius 1 of ellipse
    r₂,#Radius 2 of ellipse
    L; #Length of solenoid in meters
    StrandDiam = 0.3e-3,
    InsDiam = 1.e-3,
    ϵr_ins = 4,
    NPts_Coil = 100,
    x = -.3:0.1:0.3, #FOV Vectors
    y = x,
    z = x,
)

##
MeasLayers = 55#Number of concentric layers to calc field
PtsPerLayer =  6
A_All = zeros(MeasLayers*PtsPerLayer,3,N)
Weights = zeros(MeasLayers*PtsPerLayer)
NPtsPerSlice = MeasLayers*PtsPerLayer
for wind_num in 1:N
    A_All[:,1,wind_num],A_All[:,2,wind_num],A_All[:,3,wind_num]  =
     MakeEllipTestPoints(r₁,r₂;Center = [0,0,(wind_num - 1) / N * L],
                                            NPts=PtsPerLayer,
                                            Layers = MeasLayers)

end
SliceArea = π*r₁*r₂
Weights = sum(A_All[:,1:2,1].^2,dims=2) #A matrix that weights each value based on radius
Weights = Weights / sum(Weights)*NPtsPerSlice  #normalizing

ΦAll = zeros(N,N,3)

for Wᵢ in 1:N
    for Aᵢ in 1:N
        DeconstructA =
                    [([A_All[i,1,Aᵢ],A_All[i,2,Aᵢ],A_All[i,3,Aᵢ]]) for i in 1:length(Weights)]
        Field3DSlice = [
            BiotSav(
                MakeEllip(
                            r₁, #Radius 1 of ellipse
                            r₂;#Radius 2 of ellipse
                            Center = [0, 0, (Wᵢ - 1) / N * L], #The solenoid's axis is in Z
                            NPts = NPts_Coil, #How discretized the windings are
                            ),
                    TestLoc;
                    MinThreshold = 0.001,
                    )
                    for TestLoc in DeconstructA
                    ]

     Bweighted = vcat( [(Field3DSlice[i].*Weights[i]) for i in 1:length(Weights)]...)
      # Bweighted = vcat( [(Field3DSlice[i]) for i in 1:length(Weights)]...)

     ΦAll[Wᵢ,Aᵢ,:] =sum(Bweighted,dims=1).*SliceArea ./NPtsPerSlice
    # ΦAll[Wᵢ,Aᵢ,:] =sum(Bweighted,dims=1) ./NPtsPerSlice
    end
    println("$Wᵢ of $N windings calculated")
end

LMat = LfromΦMat(ΦAll)
LTotal = sum(LMat)*1e6

##

    X, Y, Z = MPI_Tools.meshgrid(x, y, z)
    MeanRad = (r₁ + r₂) / 2
    ZCoord = [(wind_num - 1) / N * L for wind_num = 1:N]
    RCoords = MeanRad .* ones(length(ZCoord))
    Coords = hcat(RCoords[:], ZCoord[:])
    L_Coil = CalcInductance(Coords)
    C_Coil =
        Ceq(2 * MeanRad, StrandDiam / 2, L / N, InsDiam / 2, ϵr_ins) / N * 3
    #D - solenoid diameter in meters
    #r - wire radius
    #p - pitch in meters
    #t - wire insulation thickness
    #ϵr - insulation relative permittivity
    SRF = 1 / (2 * π * sqrt(L_Coil * C_Coil))

    CoilProps, Rac_Straight, Rac_Coiled, SkinDepth =
        LitzWire_ACRes_CalcFunk(StrandDiam, N, L, MeanRad * 2; f = [100e3 1e6])



    return L_Coil, C_Coil, SRF, DefField, X, Y, Z
end


Field3D =DefField
IndNum = 2
Xs = X[IndNum, :, :]
Ys = Y[IndNum, :, :]
Zs = Z[IndNum, :, :]
FieldSlice = sum(Field3D, dims = 4)[IndNum, :, :]

pcolormesh( FieldSlice)
X,Y = MPI_Tools.meshgrid(x,y)


function LfromΦMat(ΦMat)

    N = size(ΦMat)[1]
    LMat = zeros(N,N)
    for i in 1:N
        for j in 1:N
            MagΦᵢ = sqrt(sum(ΦMat[i,i,:].^2))
            LMat[i,j] = sum(ΦMat[i,j,:].*ΦMat[i,i,:]) / MagΦᵢ
        end
    end

    return LMat
end
