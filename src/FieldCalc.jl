## Biot savart law
using LinearAlgebra
using PyPlot
using MPI_Tools
using Interpolations

VecDist(X::Array) = [ X[2,1]-X[1,1], X[2,2]-X[1,2] , X[2,3]-X[1,3]]
#X is the 2 coordinates to take distance between in format of [x₁,y₁,z₁;x₂,y₂,z₂]


function BiotSav(PointPath,r;Current=1,MinThreshold = 0.01)
    PointPath = vcat(PointPath,PointPath[1,:]')
    NPts = length(PointPath[:,1])
    dB = zeros(NPts-1,3)

    minDist = minimum(sqrt.(sum((PointPath .- repeat(transpose(r[:]),NPts,1)).^2,dims=2)))
    if minDist>= MinThreshold
        for I in 2:NPts
            dL = VecDist(PointPath[I-1:I,:])
            MeanPoint = sum(PointPath[I-1:I,:],dims=1)[:] ./ 2
            Rprime = r .- MeanPoint

            RDist = sqrt(sum(Rprime.^2))
            R̂ = Rprime/RDist

                dB[I-1,:] = μ₀/(4*π) .* Current.*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^2

        end
    else
    #    plot(r[1],r[2],"r*")
    end
    sum(dB,dims=1)
end


function MakeEllip(r₁,r₂;Center = [0,0,0],NPts=100)
    X,Y,Z = ([r₁ .* cos.(I/NPts *2*π) for I in 1:NPts].+Center[1],
             [r₂ .* sin.(I/NPts *2*π) for I in 1:NPts].+Center[2],
             [0 for I in 1:NPts].+Center[3])
             Coords = hcat(X,Y,Z)

end

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
 function FieldMapPointPath(PointPath, NumLayers; WeightRadius = false,InvWeights = false)

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
         Weights = xTP
     else
          Weights = ones(size(xTP))
      end
      if InvWeights
          Weights = 1 ./Weights
      end
     BMag = [BiotSav(PointPath,TP_Arr[i];Current = Weights[i],MinThreshold =1e-8)[3] for i in 1:length(TP_Arr)]


     BPlot = findall(xx-> xx!=0,BMag)
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
        xVecUp = vcat(xVecUp,UpSamp,Coords[i,1])
        UpSampy = LinRange(yVec[i],yVec[i+1],NElement)[:]
        yVecUp = vcat(yVecUp,UpSampy,Coords[i,2])
        UpSampz = LinRange(zVec[i],zVec[i+1],NElement)[:]
        zVecUp = vcat(zVecUp,UpSampz,Coords[i,3])
    end
PointPath = hcat(xVecUp,yVecUp,zVecUp)
if PlotOn
    scatter3D(xVecUp,yVecUp,zVecUp)
end
    return PointPath

end
