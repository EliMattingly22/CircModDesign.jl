# using Elliptic
"""
This function calculates the inductance for an N array of coaxial wire loops

The input is an N x 2 array of coordinates [R₁ Θ₁;...;Rₙ Θₙ]
The output is the overall inductance including any mutual inductances

For strange arrangements of loops, this may be useful, but for simple solenoids, other methods are desireable
Empirical formulas (``SimpleInduct'') will be more accurate and less cumbersome in those cases.
"""
function CalcInductance(Coords)
    #Thus function takes in two vectors listing the centroid [in polar coords]
    #of series-wound strands of wire to form a coil. It also takes in the index
    #(of the vector) that different coil "segments" are wound from. An example
    #being if you wanted to break up a Helmholtz pair into the two halves to
    #look at the self-inductances & mutual inductances of each independently
    #you would mark the middle index value in Segment Index List

    for i = 1:(length(Coords[:, 1])-1)
        println(i)
        for j = (i+1):length(Coords[:, 1])
            if Coords[j, :] == Coords[i, :]
                Coords[j, :] = Coords[j, :] .+ 1e-8
            end
        end
    end

    R_List = Coords[:, 1]
    Z_List = Coords[:, 2]
    L = length(R_List)

    μ₀ = 4 * π * 1e-7

    RMesh, AMesh = MPI_Tools.meshgrid(R_List,dims=2)

    KVal = similar(RMesh)
    EValue = similar(RMesh)

    KVal, EValue = (
        [ellipke.(fFunk(R_List, Z_List) .^ 2)[i, j][1] for i = 1:L, j = 1:L],
        [ellipke.(fFunk(R_List, Z_List) .^ 2)[i, j][2] for i = 1:L, j = 1:L],
    )


    M =
        μ₀ .* sqrt.(RMesh .* AMesh) .* (
            (2 .* (KVal .- EValue)) ./ fFunk(R_List, Z_List) .-
            fFunk(R_List, Z_List) .* KVal
        ) #https://thompsonrd.com/OSEE-inductance.pdf
    TotalL = sum(M)


    return TotalL
end

function fFunk(r,z)

    RMesh, AMesh = MPI_Tools.meshgrid(r,dims=2)
    ZMesh_A, ZMesh_B = MPI_Tools.meshgrid(z,dims=2)
    ZDelta = ZMesh_A .- ZMesh_B
    out = sqrt.( (4 .* RMesh.*AMesh ) ./ (ZDelta.^2 .+ (RMesh .+ AMesh).^2) )
    replace!(x -> x>=1 ? 1-eps() : x,out)


    return out
end

"""
SIMPLE INDUCTANCE FORMULAS FOR RADIO COILS*
BY HAROLD A. WHEELER
This is an empirical formula for inductance
    Takes in turns(N), radius (a), and length (b)
    the formulation was designed for inches, but by default
    this function assumes meters.
"""
function SimpleInduct(N,a,b;units = "meters")
    if units=="meters"
        a = a/25.4e-3
        b = b/25.4e-3
    elseif units == "inches"
    else
    error("unknown units")
    end

    L = (a^2*N^2 / (9*a + 10*b))*1e-6
end
