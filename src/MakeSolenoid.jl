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
Xtp,Ytp,Ztp = MakeEllipTestPoints(r₁,r₂;Center = [0,0,0],NPts=6,Layers = 3)
for wind_num = 2:N
    Xtp1,Ytp1,Ztp1 = MakeEllipTestPoints(r₁,r₂;Center = [0,0,(wind_num - 1) / N * L],NPts=6,Layers = 3)
    Xtp = vcat(Xtp,Xtp1)
    Ytp= vcat(Ytp,Ytp1)
    Ztp = vcat(Ztp,Ztp1)
end


    DefField = [
        BiotSav(
            MakeEllip(
                r₁, #Radius 1 of ellipse
                r₂;#Radius 2 of ellipse
                Center = [0, 0, (wind_num - 1) / N * L], #The solenoid's axis is in Z
                NPts = NPts_Coil, #How discretized the windings are
            ),
            [i, j, k];
            MinThreshold = 0.05,
        )[3] for i in Xtp, j in Ytp, k in Ztp, wind_num = 1:N
    ]
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



# IndNum = 2
# Xs = X[IndNum, :, :]
# Ys = Y[IndNum, :, :]
# Zs = Z[IndNum, :, :]
# FieldSlice = sum(Field3D, dims = 4)[IndNum, :, :]
#
# pcolormesh(Ys, Zs, FieldSlice)
# X,Y = MPI_Tools.meshgrid(x,y)
