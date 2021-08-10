include("NoiseModeling.jl")

function NF2NoiseDensity(R,NF_db)
    #NF = 1+en/Thermal Noise assuming minimal current noise

    NF = 10^(NF_db/20)
    en = (NF-1)*ThermalNoise(R)
end
JFET_NF = 1
JFET_NF_RVal = 1
RLoad = 20
NumJFETs = 20
TransferAdmittance = 35e-3#Siemens.
TotalTransconductance = NumJFETs*TransferAdmittance
TotalGain = RLoad*TotalTransconductance
ExpectedNoise_JFets = NF2NoiseDensity(R,NF_db)
ThermalNoise_RLoad_InputRef = ThermalNoise(RLoad)/TotalGain
