

abstract type Toroid end

struct DCore<:Toroid

    FlatHeight     ::Float64
    MaxHeight      ::Float64 #Height from CENTER LINE, so the actual core is 2 times as tall overall. This is the half-height effectively
    RadiusAtPeak   ::Float64 #This is how far from the center axis the peak of the geometry is, not the radius of curvature at the apex
    ID             ::Float64 #Thner diameter of the core
    OD             ::Float64 #Thner diameter of the core
    Turns          ::Int64
    Layers         ::Int64
    WireLength     ::Float64
    Resistance     ::Float64
end

struct Circle<:Toroid

    CoreRadius     ::Float64
    CenterRadius   ::Float64
    ID             ::Float64 #Thner diameter of the core
    OD             ::Float64 #Thner diameter of the core
    Turns          ::Int64
    Layers         ::Int64
    WireLength     ::Float64
    Resistance     ::Float64
end
struct GeneralParams<:Toroid
    SingleLayerInductance::Float64
    WireDiameter::Float64
    TargetInductance::Float64

end

mutable struct Geom
    Circ::DCore
    DCore::Circle
    General::GeneralParams
end

function ToroidOptimizer(Dia,LTarget;NumLayers=2,CoreMu=1,Alpha=2,CuFillFactor = 1)
## Summary: To use this function enter the wire diameter and the target
#inductance (in meters and Henries). The calculations assume "NumLayers" fully
#wound layers as it is a good balance between heat dissipation and
#geometric efficiency (2 is a good option, unless heat is an issue). The output is a struct containing the geometric
#information for a toroid with a D shaped core, and a circular core.


## Constants
μ_o = 4*pi*10^(-7)
μ = μ_o*CoreMu
L0 = μ*Dia/(2*pi)
L = LTarget/NumLayers^2 #Use the assumption that there will be 'NumLayers' fully wound layers...
#If there are two fully wound layers, it can be calculated at just one
#layer with 1/4 the inductance. Then double it when winding
LperL0 = L/L0 #An intermediate dimensionless parameter that comes from the papers which these formulas are derived.
SingleLayerInductance = L
WireDiameter = Dia
TargetInductance = LTarget
General = GeneralParams(SingleLayerInductance,WireDiameter,TargetInductance)

## D shaped core
#Alpha = 2; #User can change this to an integer between 2 and 5 (inclusive) varies the aspect ratio of the D
#A high alpha value such as 5 will be a taller/bigger core with fewer
#turns. It is theoretically more efficient but will probably have more
#leakage flux.

# AlphaMat Key. The columns represent [alpha e h p s t]

AlphaMat = [2 0.26 0.645 3.6  0.72 1.77
    3 0.85 1.5   8.0  2.74 4.42
    4 1.6  2.4   12.8 5.76 8.09
    5 2.45 3.4   17.9 9.61 12.6]; #Table 1 From  "D-Shaped toroidal cage inductors" P.N. Murgatroyd & D. Belahrache,1989

T = AlphaMat[Alpha-1,6] #Picking a value in the table for clarity
S = AlphaMat[Alpha-1,5]; #Also Selecting value for clarity
P = AlphaMat[Alpha-1,4]; # As above

KFunk_D(k) =(2*pi)^.5*(S/(P^(3/2)))*k.^(3/2)+0.25*k-LperL0 #From "The Optimal Form for Coreless Inductor", P. Murgatroyd IEEE TMI, 25 No. 3 1989
K =mybisect(KFunk_D, 0,1e5,50); #solving for ideal K with bisection method
WireLength_D = K*Dia; #Equation 6 in "Economic designs for..."
Turns_D = (2*pi*K/P)^(1/2); # Also from P. Murgatroyd, section 4, Eq. 20 in "Economic designs for single-layer toroidal inductors"
# TurnsD = sqrt(2*pi*L/(T*μ*B)); #Turns needed in a D shaped core toroid
# ID_DToroid = 2*B;
# B = WireLength_D/(Turns_D*P); #B=inner radius, so B = total length/(number of turns*Perimeter of each turn)
B = Dia/2+Turns_D*Dia/(2*pi); #Making it so that the wires all tough on the inner edge of the core
## For the following reference Fig 4 in "D-Shaped toroidal cage inductors" P.N. Murgatroyd & D. Belahrache,1989
FlatHeight = B*AlphaMat[Alpha-1,2];
MaxHeight = B*AlphaMat[Alpha-1,3]; #Height from CENTER LINE, so the actual core is 2 times as tall overall. This is the half-height effectively
RadiusAtPeak = Alpha^.5*B; #This is how far from the center axis the peak of the geometry is, not the radius of curvature at the apex
ID = B*2;#Thner diameter of the core
OD = 2*Alpha*B;
Turns = round(Turns_D*NumLayers);
Layers = NumLayers;
WireLength = WireLength_D*NumLayers;
D_WireR = Length2Resist(WireLength,Dia;FillFac =CuFillFactor )
DCoreGeom = DCore(FlatHeight,MaxHeight,RadiusAtPeak,ID,OD,round(Turns),Layers,WireLength,D_WireR)

## Two-layer circular core


KFunk(k) = 0.2722*k.^(3/2)+0.25*k-LperL0; #From "The Optimal Form for Coreless Inductor", P. Murgatroyd IEEE TMI, 25 No. 3 1989
K =mybisect(KFunk, 0,1e5,50); #solving for ideal K with bisection method
WireLength = K*Dia;
Turns = 0.8165*K^(1/2); # Also from P. Murgatroyd, section 2 in "Economic designs for single-layer toroidal inductors"
R = WireLength/(2*pi*Turns);# Equation 2 in "Economic designs for single-layer toroidal inductors",P. Murgatroyd
T = Dia/(sin(pi/Turns)*2)+R;# Equation 3 in "Economic designs for single-layer toroidal inductors"

Turns = Turns*NumLayers
WireLength = WireLength*NumLayers
Layers=NumLayers
ID = 2*(T-R)
OD = 2*(T+R)
CoreRadius = R
CenterRadius = T
Circ_WireR = Length2Resist(WireLength,Dia;FillFac =CuFillFactor )
CircleGeom = Circle(CoreRadius,CenterRadius,ID,OD,round(Turns),Layers,WireLength,Circ_WireR )


## Back-checking the work with two approximations.

Approx_EnclosedArea_DCore = DCoreGeom.MaxHeight*2*(DCoreGeom.OD-DCoreGeom.ID)/2*.75;#This ROUGH approximation is assuming the core fills 75# of the circumscribed rectangle around the core (roughly a semi-cirle)
SanityCheckInductance = μ*DCoreGeom.Turns^2*Approx_EnclosedArea_DCore/(pi*(DCoreGeom.OD+DCoreGeom.ID)/2);



MeanPathLength = (DCoreGeom.OD+DCoreGeom.ID)*pi

CoreFlux(MeanPathLength,Approx_EnclosedArea_DCore,DCoreGeom.Turns;CoreMu=CoreMu,Current=1e-3,AirGapSize = 0)


Laa = μ*CircleGeom.Turns^2*R^2/(2*T);#Another formula for inductance of circular
# toroid from hyperphysics
println("Target inductance is: $(round(LTarget*1e6)) μH")
println("Sanity check inductance is: $(round(SanityCheckInductance*1e6)) μH")
println("
FlatHeight     :$(round(DCoreGeom.FlatHeight*1000;sigdigits=2)) mm
MaxHeight      :$(round(DCoreGeom.MaxHeight*1000;sigdigits=2)) mm
RadiusAtPeak   :$(round(DCoreGeom.RadiusAtPeak*1000;sigdigits=2)) mm
ID             :$(round(DCoreGeom.ID*1000;sigdigits=2)) mm
OD             :$(round(DCoreGeom.OD*1000;sigdigits=2)) mm
Turns          :$(DCoreGeom.Turns) Turns
Layers         :$(DCoreGeom.Layers) Turns
WireLength     :$(round(DCoreGeom.WireLength;sigdigits=2)) meters
WireResist     :$(round(DCoreGeom.Resistance;sigdigits=2)) Ω")
CoilGeom = Geom(DCoreGeom,CircleGeom,General)
return CoilGeom
end

function CoreFlux(Le,Ae,Turns;CoreMu=1,Current=1,AirGapSize = 0)
    ## This assumes a simple geometry where the gap area equals Ferrite area
    μ_o = 4*pi*1e-7
    μ = μ_o*CoreMu
    Φ = μ_o*Turns*Current/((Le-AirGapSize)/(CoreMu*Ae)+AirGapSize/Ae)
    L = Turns*Φ
    Flux = Φ/Ae
    println("The total B field within the core given $(Current*1000)mA
     is $(round(Flux*1000;sigdigits=2))mT")
     println("~~~~~~~~~~~~")
    return Flux

end

function Length2Resist(L,Diam;ρ=1.68e-8,FillFac = 1)
    if Diam>.1
        println("Diameter must be in METERS, it was automatically multiplied by 1e-3 in case you entered it in mm")
        Diam = Diam*1e-3
    end
    A = pi/4 .* Diam.^2
    return L*ρ*FillFac/A
end
