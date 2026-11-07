function copied(model::ExactPiModel, b::Int, p, h::Real, sequence::Symbol)
    br = model.branches[b]
    len = br.len isa Function ? br.len(p) : br.len
    z, y = br.zy(p, model.f1 * h, sequence)
    γ, zc = sqrt(z * y), sqrt(z / y)
    zser = zc * sinh(γ * len)
    yhalf = tanh(γ * len / 2) / zc
    return (br_r = real(zser) / br.zbase, br_x = imag(zser) / br.zbase,
        g_fr = real(yhalf) * br.zbase, b_fr = imag(yhalf) * br.zbase,
        g_to = real(yhalf) * br.zbase, b_to = imag(yhalf) * br.zbase)
end
