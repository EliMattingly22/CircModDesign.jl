function determineK(L₁₁,L₂₂,Lₘ)
    K = Lₘ/√(L₁₁*L₂₂)
end

function determineCancelledCenter(L₁₁,L₂₂)
    Lₘ = L₁₁
    Lₐ = Lₘ
    Lᵦ = L₂₂-Lₘ
    Lc = L₁₁-Lₘ
    K = Lₘ/√(L₁₁*L₂₂)
    return K,Lₐ,Lᵦ,Lc
end



