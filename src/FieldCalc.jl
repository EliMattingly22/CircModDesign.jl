## Biot savart law
# using LinearAlgebra
# using PyPlot
# using MPI_Tools
# using Interpolations

VecDist(X::Array) = [ X[2,1]-X[1,1], X[2,2]-X[1,2] , X[2,3]-X[1,3]]
#X is the 2 coordinates to take distance between in format of [x₁,y₁,z₁;x₂,y₂,z₂]

"""
This function takes in:
    "Point Path" which is an Nx3 set of coordinates that define a current loop
    "r" which is a 1x3 vector that is a point in space to evaluate the 3D B field at

    kwargs:
        "current" is a scalar multiplier. Default is 1 Amp 
        "MinThreshold" is a minimum distance that the "r" point must be from any wire
            The function tends to inf near wires

"""
function BiotSav(PointPath,r;Current=1 ,MinThreshold = 0.01)
    PointPath = vcat(PointPath,PointPath[1,:]')
    NPts = length(PointPath[:,1])
    # dB = zeros(NPts-1,3)
    dB = [0. , 0., 0.]
    minDist = minimum(sqrt.(sum((PointPath .- repeat(transpose(r[:]),NPts,1)).^2,dims=2)))
    if minDist>= MinThreshold
        for I in 2:NPts
            dL = VecDist(PointPath[I-1:I,:])
            MeanPoint = sum(PointPath[I-1:I,:],dims=1)[:] ./ 2
            Rprime = r .- MeanPoint

            RDist = sqrt(sum(Rprime.^2))
            R̂ = Rprime/RDist

                if Current isa AbstractArray
                    # dB[I-1,:] = μ₀/(4*π) .* Current[I].*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2
                    dB .+= μ₀/(4*π) .* Current[I].*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2
                
                else
                    # dB[I-1,:] = μ₀/(4*π) .* Current.*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2
                    dB .+= μ₀/(4*π) .* Current.*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2
                end

        end
    # else
    #    plot(r[1],r[2],"r*")
    end
    sum(dB,dims=1)
end

function BiotSav(PointPath,dL,r,L;Current=1.0 ,MinThreshold = 0.01)
    # PointPath = vcat(PointPath,PointPath[1,:]')
    # NPts = length(PointPath[:,1])
    # dB = zeros(NPts-1,3)
    dB = [0. , 0., 0.]
    # minDist = minimum(sqrt.(sum((PointPath .- repeat(transpose(r[:]),NPts,1)).^2,dims=2)))
    # if minDist>= MinThreshold
        for I in 2:L
          
            # MeanPoint = 
            Rprime = r .- PointPath[I,:]

            RDist = (sum(Rprime.^2))
            R̂ = Rprime/sqrt(RDist)

                # if Current isa AbstractArray
                #     # dB[I-1,:] = μ₀/(4*π) .* Current[I].*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2
                #     dB .+= 1e-7 .* Current[I].*LinearAlgebra.cross(dL[I-1,:],R̂[:]) ./ (RDist)^2
                
                # else
                    # dB[I-1,:] = μ₀/(4*π) .* Current.*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2
                    dB .+= 1e-7 .* Current.*LinearAlgebra.cross(dL[I-1,:],R̂[:]) ./ (RDist)
                # end

        end
    # else
    #    plot(r[1],r[2],"r*")
    # end
    dB
end

"""
Makes a point path of an ellipse with r₁,r₂ as two radii
normal direction is assumed to be in Z
kwargs: Center is the center of path [X,Y,Z]
NPts is the number of points along path


"""
function MakeEllip(r₁,r₂;Center = [0,0,0],NPts=100)
    X,Y,Z = ([r₁ .* cos.(I/NPts *2*π) for I in 1:NPts].+Center[1],
             [r₂ .* sin.(I/NPts *2*π) for I in 1:NPts].+Center[2],
             [0 for I in 1:NPts].+Center[3])
    
    Coords = hcat(X,Y,Z)

end


"""
Simple analytical formula for the field on axis of a Circle
    R - radius in meters
    z - axial offset in meters
    
    kwargs:
        I = current in coil

"""
function FieldOnAxis_Circ(R,z;I=1)
    Bz = μ₀*2*pi*R^2*I ./ (4*π* (z^2 + R^2)^(3/2))
end

# x = -.3:.01:.3
# y = x
# DefField = [BiotSav(MakeEllip(1,1;Center = [0,0,0],NPts=1000),[i,j,k];MinThreshold= 0.05)[3]
#             for i in x,
#                 j in y,
#                 k=0]
# X,Y = MPI_Tools.meshgrid(x,y)
# PyPlot.surf(X[:],Y[:],DefField[:],cmap="jet")

function MakeEllipTestPoints(r₁,r₂;Center = [0,0,0],NPts=6,Layers = 3)
    LL=1
    X,Y,Z = ([r₁*LL/(Layers+1) .* cos.(I/NPts *2*π) for I in 1:NPts].+Center[1],
             [r₂ *LL/(Layers+1) .* sin.(I/NPts *2*π) for I in 1:NPts].+Center[2],
             [0 for I in 1:NPts].+Center[3])
             Coords = hcat(X,Y,Z)
    for LL in 2:Layers
        X₁,Y₁,Z₁ = ([r₁*LL/(Layers+1) .* cos.(I/NPts *2*π) for I in 1:NPts].+Center[1],
                 [r₂ *LL/(Layers+1) .* sin.(I/NPts *2*π) for I in 1:NPts].+Center[2],
                 [0 for I in 1:NPts].+Center[3])
                 Coords = hcat(X,Y,Z)
        X = vcat(X,X₁)
        Y = vcat(Y,Y₁)
        Z = vcat(Z,Z₁)

     end
     return X,Y,Z
 end


"""
This function takes in a pointpath [X,Y,Z] array and maps the interior

Inputs are
PointPath [X,Y,Z] , N x 3 array
NumLayers - number of interior test points

"""
 function FieldMapPointPath(PointPath, NumLayers; WeightRadius = false,InvWeights = false, PlotAxes = nothing)

     xTP = PointPath[:,1]
     x   = PointPath[:,1]
     yTP = PointPath[:,2]
     y   = PointPath[:,2]
     x₁ = minimum(xTP)
     x₂ = maximum(xTP)
     y₁ = minimum(yTP)
     y₂ = maximum(yTP)
     xmid = (x₁+ x₂)/2
     ymid = (y₁+ y₂)/2



     for i in 1:(NumLayers-2)
         xTP = vcat(xTP,(x .- xmid) .* (i/NumLayers) .+ xmid)
         yTP = vcat(yTP,(y .- ymid) .* (i/NumLayers) .+ ymid)
     end
     TPList = hcat(xTP,yTP,zeros(size(yTP)))
     TP_Arr =[([TPList[i,1],TPList[i,2],TPList[i,3]]) for i in 1:length(xTP)]
     if WeightRadius
         Weights = PointPath[:,1]
     else
          Weights = ones(size( PointPath[:,1]))
      end
      if InvWeights
          Weights = 1 ./Weights
      end
      push!(Weights,Weights[1])
     BMag = [BiotSav(PointPath,TP_Arr[i];Current = Weights,MinThreshold =1e-8)[3] for i in 1:length(TP_Arr)]

     BMag = abs.(BMag)

     BPlot = findall(xx-> xx!=0,BMag)
     if PlotAxes!== nothing
        PyPlot.axes(PlotAxes)
     end

     surf(xTP[BPlot],yTP[BPlot],BMag[BPlot],cmap="jet")


 end

"""
Makes a point path for the BiotSav function
    Input is an Nx3 array of coordinates where it is structured as:
    [x₁, y₁, z₁
    ..
    ..
    xₙ, yₙ, zₙ]

Coords = [1 0 0
          2 0 0
          2 1 0
          1 1 0]
"""
function MakeRectPointPath(Coords;NElement = 50, PlotOn=false)
    Coords = vcat(Coords, reshape(Coords[1,:],1,3))
    xVec = Coords[:,1]
    yVec = Coords[:,2]
    zVec = Coords[:,3]
    xVecUp = [Coords[1,1]]
    yVecUp = [Coords[1,2]]
    zVecUp = [Coords[1,3]]
    for i in 1:(length(xVec)-1)
        UpSamp = LinRange(xVec[i],xVec[i+1],NElement)[:]
        xVecUp = vcat(xVecUp,UpSamp,Coords[i+1,1])
        UpSampy = LinRange(yVec[i],yVec[i+1],NElement)[:]
        yVecUp = vcat(yVecUp,UpSampy,Coords[i+1,2])
        UpSampz = LinRange(zVec[i],zVec[i+1],NElement)[:]
        zVecUp = vcat(zVecUp,UpSampz,Coords[i+1,3])
    end
    PointPath = hcat(xVecUp,yVecUp,zVecUp)
    if PlotOn
        scatter3D(xVecUp,yVecUp,zVecUp)
    end
    return PointPath

end


function WireNearField_B(r,Wire_R;I=1)
    μ₀ = 4*π*1e-7
    if r<Wire_R
        B = μ₀*r*I / (2*π*Wire_R^2)
    elseif r==Wire_R
        B = μ₀*I / (2*π*Wire_R)
    else
        B = μ₀*I / (2*π*r)
    end
    return B
end

function DerivField(r,Wire_R=0.001;I=1)
    μ₀ = 4*π*1e-7
    if r<Wire_R
        B = 0 ##Update later.
    elseif r==Wire_R
        B=0
    else
        B = -1*μ₀*I / (2*π*r^2)
    end
    return B
end


function BiotSavSpeedTest(NPathPoints, NTestPoints)
    PointPath =  MakeEllip(1. ,1. ;NPts = NPathPoints)
    TestPoints,W = PP_2_TestPoints(PointPath)
    PointPath = vcat(PointPath,PointPath[1,:]')
    dL = vcat([CircModDesign.VecDist(PointPath[I-1:I,:]) for I in 2:length(PointPath[:,1])]'...)
    
    L = length(PointPath[:,1])    
    dBMat = zeros(size(TestPoints))
   @time begin
    for I in 1:length(TestPoints[:,1])
    dBMat[I,:] = BiotSav(PointPath,dL,TestPoints[I,:],L)
    end
   end
   dBMat
end