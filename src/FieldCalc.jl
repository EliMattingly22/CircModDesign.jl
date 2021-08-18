## Biot savart law
μ₀ = 4*π*1e-7

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

                dB[I-1,:] = μ₀/(4*π) .* Current.*LinearAlgebra.cross(dL,R̂[:]) ./ (RDist)^3

        end
    else
        plot(r[1],r[2],"r*")
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

x = -1:.1:1
y = x
DefField = [BiotSav(PointPath,[i,j,k];MinThreshold= 0.05)[3]
            for i in x,
                j in y,
                k=0]
X,Y = MPI_Tools.meshgrid(x,y)
surf(DefField,cmap="jet")
