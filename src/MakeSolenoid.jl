function MakeEllipSolenoid(N,r₁,r₂,L;
    NPts_Coil = 100,
    x = -.3:.1:.3,
    y = x,
    z = x
)

    DefField = [BiotSav(MakeEllip(1,1;Center = [0,0,(wind_num-1)/N*L],NPts=NPts_Coil),[i,j,k];MinThreshold= 0.05)[3]
                    for i in x,
                        j in y,
                        k in z,
                        wind_num in 1:N]
X,Y,Z = MPI_Tools.meshgrid(x,y,z)
return DefField,X,Y,Z
end
                    
IndNum = 2
Xs = X[IndNum,:,:]
Ys = Y[IndNum,:,:]
Zs = Z[IndNum,:,:]
FieldSlice = sum(Field3D,dims=4)[IndNum,:,:]

pcolormesh(Ys,Zs,FieldSlice)
        # X,Y = MPI_Tools.meshgrid(x,y)
