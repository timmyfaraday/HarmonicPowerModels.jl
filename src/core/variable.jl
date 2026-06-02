################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
################################################################################

# fairness principle
""
function variable_fairness_principle(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    # maximum efficiency
    if pm.data["principle"] == "maximum efficiency"
        # no additional variables
    end

    # absolute equality
    if pm.data["principle"] == "absolute equality"
        # no additional variables 
    end

    # maximin
    if pm.data["principle"] == "maximin"
        cmh = _PMs.var(pm, nw)[:cmh] = JuMP.@variable(pm.model, base_name="$(nw)_cmh",
                start = 0.0
        )

        if bounded
            JuMP.set_lower_bound(cmh, 0.0)
            # JuMP.set_upper_bound(cmd[d], c_rating)                            # @Hakan: dit is ook bij :cmd, van waar komt deze c rating, the fundamental component - wat als die er niet is?
        end

        report && (_PMs.sol(pm, nw, :fairness)[:cmh] = cmh)
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["principle"] == "Kalai-Smorodinsky bargaining"
        fh =  _PMs.var(pm, nw)[:fh] = JuMP.@variable(pm.model, base_name="$(nw)_fh",
                start = 0.0
        )

        if bounded
            JuMP.set_lower_bound(fh, 0.0)
            JuMP.set_upper_bound(fh, 1.0)
        end

        report && (_PMs.sol(pm, nw, :fairness)[:fh] = fh)
    end
end

# bus
""
function variable_bus_voltage_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    vr = _PMs.var(pm, nw)[:vr] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_vr",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "vr_start", 1.0)
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vr[i], -bus["vmax"])
            JuMP.set_upper_bound(vr[i],  bus["vmax"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vr, _PMs.ids(pm, nw, :bus), vr)
end

""
function variable_bus_voltage_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    vi = _PMs.var(pm, nw)[:vi] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_vi",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "vi_start")
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vi[i], -bus["vmax"])
            JuMP.set_upper_bound(vi[i],  bus["vmax"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vi, _PMs.ids(pm, nw, :bus), vi)
end

# branch
""
function variable_branch_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cr = _PMs.var(pm, nw)[:cr] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_cr",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "cr_start")
    )

    if bounded
        for (l,i,j) in _PMs.ref(pm, nw, :arcs)
            branch = _PMs.ref(pm, nw, :branch, l)
            JuMP.set_lower_bound(cr[(l,i,j)], -branch["c_rating"])
            JuMP.set_upper_bound(cr[(l,i,j)],  branch["c_rating"])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), cr)
end
""
function variable_branch_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    ci = _PMs.var(pm, nw)[:ci] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_ci",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "ci_start")
    )

    if bounded
        for (l,i,j) in _PMs.ref(pm, nw, :arcs)
            branch = _PMs.ref(pm, nw, :branch, l)
            JuMP.set_lower_bound(ci[(l,i,j)], -branch["c_rating"])
            JuMP.set_upper_bound(ci[(l,i,j)],  branch["c_rating"])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), ci)
end
""
function variable_branch_series_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csr = _PMs.var(pm, nw)[:csr] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_csr",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "csr_start", 0.0)
    )

    if bounded
        for (b, branch) in _PMs.ref(pm, nw, :branch)
            JuMP.set_lower_bound(csr[b], -branch["c_rating"])
            JuMP.set_upper_bound(csr[b],  branch["c_rating"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :branch, :csr_fr, _PMs.ids(pm, nw, :branch), csr)
end
""
function variable_branch_series_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csi = _PMs.var(pm, nw)[:csi] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_csi",
        start=_PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "csi_start", 0.0)
    )

    if bounded
        for (b, branch) in _PMs.ref(pm, nw, :branch)
            JuMP.set_lower_bound(csi[b], -branch["c_rating"])
            JuMP.set_upper_bound(csi[b],  branch["c_rating"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :branch, :csi_fr, _PMs.ids(pm, nw, :branch), csi)
end

# filter
""
function variable_filter_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    crf = _PMs.var(pm, nw)[:crf] = JuMP.@variable(pm.model,
            [f in _PMs.ids(pm, nw, :filter)], base_name="$(nw)_crf",
            start=_PMs.comp_start_value(_PMs.ref(pm, nw, :filter, f), "crf_start", 0.0)
    )
    
    if bounded
        for (f, filter) in _PMs.ref(pm, nw, :filter)
            JuMP.set_lower_bound(crf[f], -filter["c_rating"])
            JuMP.set_upper_bound(crf[f],  filter["c_rating"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :filter, :crf, _PMs.ids(pm, nw, :filter), crf)
end
""
function variable_filter_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cif = _PMs.var(pm, nw)[:cif] = JuMP.@variable(pm.model,
            [f in _PMs.ids(pm, nw, :filter)], base_name="$(nw)_cif",
            start=_PMs.comp_start_value(_PMs.ref(pm, nw, :filter, f), "cif_start", 0.0)
    )
    
    if bounded
        for (f, filter) in _PMs.ref(pm, nw, :filter)
            JuMP.set_lower_bound(cif[f], -filter["c_rating"])
            JuMP.set_upper_bound(cif[f],  filter["c_rating"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :filter, :cif, _PMs.ids(pm, nw, :filter), cif)
end

# xfmr 
""
function variable_xfmr_voltage_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    vrx = _PMs.var(pm, nw)[:vrx] = JuMP.@variable(pm.model, 
            [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_vrx",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "vrx_start", 1.0)
    )

    if bounded
        for (x, i, j) in _PMs.ref(pm, nw, :xfmr_arcs)
            vrx_min = - _PMs.ref(pm, nw, :bus, i)["vmax"] * _PMs.ref(pm, nw, :bus, i)["ihdmax"]
            vrx_max =   _PMs.ref(pm, nw, :bus, i)["vmax"] * _PMs.ref(pm, nw, :bus, i)["ihdmax"]

            JuMP.set_lower_bound(vrx[(x, i, j)], vrx_min)
            JuMP.set_upper_bound(vrx[(x, i, j)], vrx_max)
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :vrx_fr, :vrx_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), vrx)
end
""
function variable_xfmr_voltage_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    vix = _PMs.var(pm, nw)[:vix] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_vix",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "vix_start", 0.0)
    )

    if bounded
        for (x, i, j) in _PMs.ref(pm, nw, :xfmr_arcs)
            vix_min = - _PMs.ref(pm, nw, :bus, i)["vmax"] * _PMs.ref(pm, nw, :bus, i)["ihdmax"]
            vix_max =   _PMs.ref(pm, nw, :bus, i)["vmax"] * _PMs.ref(pm, nw, :bus, i)["ihdmax"]

            JuMP.set_lower_bound(vix[(x, i, j)], vix_min)
            JuMP.set_upper_bound(vix[(x, i, j)], vix_max)
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :vix_fr, :vix_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), vix)
end
""
function variable_xfmr_voltage_excitation_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, epsilon::Float64=1E-6)
    erx = _PMs.var(pm, nw)[:erx] = JuMP.@variable(pm.model,
            [x in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_erx",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "erx_start", 0.0)
    )

    if bounded
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            if haskey(xfmr, "erx_min")
                JuMP.set_lower_bound(erx[x], xfmr["erx_min"] + epsilon)
                JuMP.set_upper_bound(erx[x], xfmr["erx_max"] - epsilon)
            else
                JuMP.set_lower_bound(erx[x], -_PMs.ref(pm, nw, :bus, xfmr["f_bus"])["vmax"])
                JuMP.set_upper_bound(erx[x],  _PMs.ref(pm, nw, :bus, xfmr["f_bus"])["vmax"])
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :erx, _PMs.ids(pm, nw, :xfmr), erx)
end
function variable_xfmr_voltage_excitation_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, epsilon::Float64=1E-6)
    eix = _PMs.var(pm, nw)[:eix] = JuMP.@variable(pm.model,
            [x in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_eix",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "eix_start", 0.0)
    )

    if bounded
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            if haskey(xfmr, "erx_min")
                JuMP.set_lower_bound(eix[x], xfmr["eix_min"] + epsilon)
                JuMP.set_upper_bound(eix[x], xfmr["eix_max"] - epsilon)
             else
                JuMP.set_lower_bound(eix[x], -_PMs.ref(pm, nw, :bus, xfmr["f_bus"])["vmax"])
                JuMP.set_upper_bound(eix[x],  _PMs.ref(pm, nw, :bus, xfmr["f_bus"])["vmax"])
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :eix, _PMs.ids(pm, nw, :xfmr), eix)
end

""
function variable_xfmr_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    crx = _PMs.var(pm, nw)[:crx] = JuMP.@variable(pm.model,
            [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_crx",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "crx_start", 0.0)
    )

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            xfmr = _PMs.ref(pm, nw, :xfmr, x)
            JuMP.set_lower_bound(crx[(x,i,j)], -xfmr["c_rating"])
            JuMP.set_upper_bound(crx[(x,i,j)],  xfmr["c_rating"])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :crx_fr, :crx_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), crx)
end
""
function variable_xfmr_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cix = _PMs.var(pm, nw)[:cix] = JuMP.@variable(pm.model,
            [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_cix",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "cix_start", 0.0)
    )

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            xfmr = _PMs.ref(pm, nw, :xfmr, x)
            JuMP.set_lower_bound(cix[(x,i,j)], -xfmr["c_rating"])
            JuMP.set_upper_bound(cix[(x,i,j)],  xfmr["c_rating"])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :cix_fr, :cix_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), cix)
end
""
function variable_xfmr_current_series_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csrx = _PMs.var(pm, nw)[:csrx] = JuMP.@variable(pm.model,
            [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_csrx",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "csrx_start", 0.0)
    )

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            xfmr = _PMs.ref(pm, nw, :xfmr, x)
            JuMP.set_lower_bound(csrx[(x,i,j)], -xfmr["c_rating"])
            JuMP.set_upper_bound(csrx[(x,i,j)],  xfmr["c_rating"])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :csrx_fr, :csrx_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), csrx)
end
""
function variable_xfmr_current_series_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csix = _PMs.var(pm, nw)[:csix] = JuMP.@variable(pm.model,
            [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_csix",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "csix_start", 0.0)
    )

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            xfmr = _PMs.ref(pm, nw, :xfmr, x)
            JuMP.set_lower_bound(csix[(x,i,j)], -xfmr["c_rating"])
            JuMP.set_upper_bound(csix[(x,i,j)],  xfmr["c_rating"])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :csix_fr, :csix_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), csix)
end
""
function variable_xfmr_current_magnetizing_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cmrx = _PMs.var(pm, nw)[:cmrx] = JuMP.@variable(pm.model,
            [x in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_cmrx",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "cmrx_start", 0.0)
    )

    if bounded
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            JuMP.set_lower_bound(cmrx[x], 0)
            JuMP.set_upper_bound(cmrx[x], xfmr["c_rating"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :cmrx, _PMs.ids(pm, nw, :xfmr), cmrx)
end
""
function variable_xfmr_current_magnetizing_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cmix = _PMs.var(pm, nw)[:cmix] = JuMP.@variable(pm.model,
            [x in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_cmix",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "cmix_start", 0.0)
    )

    if bounded
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            JuMP.set_lower_bound(cmix[x], 0)
            JuMP.set_upper_bound(cmix[x], xfmr["c_rating"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :cmix, _PMs.ids(pm, nw, :xfmr), cmix)
end

# generator 
""
function variable_gen_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    crg = _PMs.var(pm, nw)[:crg] = JuMP.@variable(pm.model,
            [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_crg",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :gen, g), "crg_start", 0.0)
    )

    if bounded
        for (g, gen) in _PMs.ref(pm, nw, :gen)
            c_rating = gen["c_rating"]
            JuMP.set_lower_bound(crg[g], -c_rating)
            JuMP.set_upper_bound(crg[g],  c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :crg, _PMs.ids(pm, nw, :gen), crg)
end
""
function variable_gen_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cig = _PMs.var(pm, nw)[:cig] = JuMP.@variable(pm.model,
            [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_cig",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :gen, g), "cig_start", 0.0)
    )

    if bounded
        for (g, gen) in _PMs.ref(pm, nw, :gen)
            c_rating = gen["c_rating"]
            JuMP.set_lower_bound(cig[g], -c_rating)
            JuMP.set_upper_bound(cig[g],  c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :cig, _PMs.ids(pm, nw, :gen), cig)
end

# load 
""
function variable_load_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    crd = _PMs.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
            [d in _PMs.ids(pm, nw, :load)], base_name="$(nw)_crd",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, d), "crd_start", 0.0)
    )

    if bounded
        for (d, load) in _PMs.ref(pm, nw, :load)
            c_rating = load["c_rating"]
            JuMP.set_lower_bound(crd[d], -c_rating)
            JuMP.set_upper_bound(crd[d],  c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :load, :crd, _PMs.ids(pm, nw, :load), crd)
end
""
function variable_load_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cid = _PMs.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
            [d in _PMs.ids(pm, nw, :load)], base_name="$(nw)_cid",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, d), "cid_start", 0.0)
    )

    if bounded
        for (d, load) in _PMs.ref(pm, nw, :load)
            c_rating = load["c_rating"]
            JuMP.set_lower_bound(cid[d], -c_rating)
            JuMP.set_upper_bound(cid[d],  c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :load, :cid, _PMs.ids(pm, nw, :load), cid)
end
""
function variable_load_current_magnitude(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cmd = _PMs.var(pm, nw)[:cmd] = JuMP.@variable(pm.model,
            [d in _PMs.ids(pm, nw, :load)], base_name="$(nw)_cmd",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, d), "cmd_start", 0.0)
    )

    if bounded
        for (d, load) in _PMs.ref(pm, nw, :load)
            c_rating = load["c_rating"]
            JuMP.set_lower_bound(cmd[d], 0.0)
            JuMP.set_upper_bound(cmd[d], c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :load, :cmd, _PMs.ids(pm, nw, :load), cmd)
end

################################################################################
# sHHC_QCQP — PCE-expanded variables                                           #
# K = 0:δ  (e.g., 0:2 for δ = 2).  Indexed as (component_id, k).             #
# Only the k=0 coefficient is bounded; higher-order coefficients are           #
# unconstrained in sign (they are polynomial chaos expansion coefficients,     #
# not physical magnitudes).                                                    #
################################################################################

# ── bus voltage PCE ────────────────────────────────────────────────────────────

"""
    variable_bus_voltage_pce(pm::sHHC_QCQP; nw, K, bounded, report)

PCE coefficients of the real and imaginary bus voltage for harmonic network `nw`.

Creates `vr_pce[i,k]` and `vi_pce[i,k]` for each bus `i` and PCE index `k ∈ K`.

Bounds (when `bounded=true`): only the k=0 moment coefficient is box-constrained
to `[-vmax, vmax]`; higher-order coefficients are unbounded.
"""
function variable_bus_voltage_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                   K=0:2, bounded::Bool=true, report::Bool=true)
    vr_pce = _PMs.var(pm, nw)[:vr_pce] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus), k in K], base_name="$(nw)_vr_pce",
        start = 0.0
    )
    vi_pce = _PMs.var(pm, nw)[:vi_pce] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus), k in K], base_name="$(nw)_vi_pce",
        start = 0.0
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            vmax = bus["vmax"]
            JuMP.set_lower_bound(vr_pce[i, 0], -vmax)
            JuMP.set_upper_bound(vr_pce[i, 0],  vmax)
            JuMP.set_lower_bound(vi_pce[i, 0], -vmax)
            JuMP.set_upper_bound(vi_pce[i, 0],  vmax)
        end
    end

    if report
        for i in _PMs.ids(pm, nw, :bus)
            _PMs.sol(pm, nw, :bus, i)[:vr_pce] = Dict(k => vr_pce[i, k] for k in K)
            _PMs.sol(pm, nw, :bus, i)[:vi_pce] = Dict(k => vi_pce[i, k] for k in K)
        end
    end
end

# ── branch current PCE ─────────────────────────────────────────────────────────

"""
    variable_branch_current_pce(pm::sHHC_QCQP; nw, K, bounded, report)

PCE coefficients of the real and imaginary branch current for harmonic network `nw`.

Creates `cr_pce[(l,i,j),k]` and `ci_pce[(l,i,j),k]` for each arc `(l,i,j)`
and PCE index `k ∈ K`.

Bounds (when `bounded=true`): only k=0 is constrained to `[-c_rating, c_rating]`.
"""
function variable_branch_current_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                      K=0:2, bounded::Bool=true, report::Bool=true)
    cr_pce = _PMs.var(pm, nw)[:cr_pce] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs), k in K], base_name="$(nw)_cr_pce",
        start = 0.0
    )
    ci_pce = _PMs.var(pm, nw)[:ci_pce] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs), k in K], base_name="$(nw)_ci_pce",
        start = 0.0
    )

    if bounded
        for (l,i,j) in _PMs.ref(pm, nw, :arcs)
            c_rating = _PMs.ref(pm, nw, :branch, l)["c_rating"]
            JuMP.set_lower_bound(cr_pce[(l,i,j), 0], -c_rating)
            JuMP.set_upper_bound(cr_pce[(l,i,j), 0],  c_rating)
            JuMP.set_lower_bound(ci_pce[(l,i,j), 0], -c_rating)
            JuMP.set_upper_bound(ci_pce[(l,i,j), 0],  c_rating)
        end
    end

    if report
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            _PMs.sol(pm, nw, :branch, l)[:cr_pce_fr] = Dict(k => cr_pce[(l,i,j), k] for k in K)
            _PMs.sol(pm, nw, :branch, l)[:ci_pce_fr] = Dict(k => ci_pce[(l,i,j), k] for k in K)
        end
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_to)
            _PMs.sol(pm, nw, :branch, l)[:cr_pce_to] = Dict(k => cr_pce[(l,i,j), k] for k in K)
            _PMs.sol(pm, nw, :branch, l)[:ci_pce_to] = Dict(k => ci_pce[(l,i,j), k] for k in K)
        end
    end
end

# ── load (harmonic unit) current PCE ──────────────────────────────────────────

"""
    variable_load_current_pce(pm::sHHC_QCQP; nw, K, bounded, report)

PCE coefficients of the real and imaginary harmonic unit injection current for
network `nw`.

Creates `crd_pce[d,k]` and `cid_pce[d,k]` for each load (harmonic unit) `d`
and PCE index `k ∈ K`.

Bounds (when `bounded=true`): only k=0 is constrained to `[-c_rating, c_rating]`.
"""
function variable_load_current_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                    K=0:2, bounded::Bool=true, report::Bool=true)
    crd_pce = _PMs.var(pm, nw)[:crd_pce] = JuMP.@variable(pm.model,
        [d in _PMs.ids(pm, nw, :load), k in K], base_name="$(nw)_crd_pce",
        start = 0.0
    )
    cid_pce = _PMs.var(pm, nw)[:cid_pce] = JuMP.@variable(pm.model,
        [d in _PMs.ids(pm, nw, :load), k in K], base_name="$(nw)_cid_pce",
        start = 0.0
    )

    if bounded
        for (d, load) in _PMs.ref(pm, nw, :load)
            c_rating = load["c_rating"]
            JuMP.set_lower_bound(crd_pce[d, 0], -c_rating)
            JuMP.set_upper_bound(crd_pce[d, 0],  c_rating)
            JuMP.set_lower_bound(cid_pce[d, 0], -c_rating)
            JuMP.set_upper_bound(cid_pce[d, 0],  c_rating)
        end
    end

    if report
        for d in _PMs.ids(pm, nw, :load)
            _PMs.sol(pm, nw, :load, d)[:crd_pce] = Dict(k => crd_pce[d, k] for k in K)
            _PMs.sol(pm, nw, :load, d)[:cid_pce] = Dict(k => cid_pce[d, k] for k in K)
        end
    end
end

# ── lifted auxiliary variables ─────────────────────────────────────────────────

"""
    variable_lifted_voltage_pce(pm::sHHC_QCQP; nw, K, report)

Lifted auxiliary variables `w_pce[i,k]` for each bus `i` and PCE index `k ∈ K`.

Semantics: `w_pce[i,k]` is the k-th PCE coefficient of the squared voltage
magnitude `|V_{h,i}|² = Vr_{h,i}² + Vi_{h,i}²`. These variables are linked
to the voltage PCE coefficients via:
  - k=0 : SOC inequality  `w_pce[i,0] ≥ vr_pce[i,·]² + vi_pce[i,·]²` (lifted)
  - k≥1 : quadratic equality constraint from the multiplication tensor

No explicit bounds are imposed; the lifting constraints enforce `w_pce[i,0] ≥ 0`
and the Cantelli chance constraints use `w_pce` to linearise the variance term.
"""
function variable_lifted_voltage_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                      K=0:2, report::Bool=true)
    w_pce = _PMs.var(pm, nw)[:w_pce] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus), k in K], base_name="$(nw)_w_pce",
        start = 0.0
    )

    if report
        for i in _PMs.ids(pm, nw, :bus)
            _PMs.sol(pm, nw, :bus, i)[:w_pce] = Dict(k => w_pce[i, k] for k in K)
        end
    end
end

# ── Transformer PCE variables ─────────────────────────────────────────────────

"""
    variable_xfmr_voltage_pce(pm::sHHC_QCQP; nw, K, bounded, report)

PCE coefficients of the transformer winding and core excitation voltages.

Creates:
- `vrx_pce[(x,i,j), k]`, `vix_pce[(x,i,j), k]` — real/imaginary winding arc voltages
- `erx_pce[x, k]`, `eix_pce[x, k]` — real/imaginary core excitation voltages

Bounds (k=0 only): winding arcs ±vmax·ihdmax; excitation ±vmax (from f_bus).
Higher-order PCE coefficients are unbounded in sign.
"""
function variable_xfmr_voltage_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                    K=0:2, bounded::Bool=true, report::Bool=true)
    # ── winding arc voltages ──────────────────────────────────────────────────
    vrx_pce = _PMs.var(pm, nw)[:vrx_pce] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs), k in K],
        base_name="$(nw)_vrx_pce", start=0.0
    )
    vix_pce = _PMs.var(pm, nw)[:vix_pce] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs), k in K],
        base_name="$(nw)_vix_pce", start=0.0
    )

    if bounded
        for (x, i, j) in _PMs.ref(pm, nw, :xfmr_arcs)
            bus   = _PMs.ref(pm, nw, :bus, i)
            vlim  = bus["vmax"] * get(bus, "ihdmax", 1.0)
            JuMP.set_lower_bound(vrx_pce[(x,i,j), 0], -vlim)
            JuMP.set_upper_bound(vrx_pce[(x,i,j), 0],  vlim)
            JuMP.set_lower_bound(vix_pce[(x,i,j), 0], -vlim)
            JuMP.set_upper_bound(vix_pce[(x,i,j), 0],  vlim)
        end
    end

    # ── core excitation voltages ──────────────────────────────────────────────
    erx_pce = _PMs.var(pm, nw)[:erx_pce] = JuMP.@variable(pm.model,
        [x in _PMs.ids(pm, nw, :xfmr), k in K],
        base_name="$(nw)_erx_pce", start=0.0
    )
    eix_pce = _PMs.var(pm, nw)[:eix_pce] = JuMP.@variable(pm.model,
        [x in _PMs.ids(pm, nw, :xfmr), k in K],
        base_name="$(nw)_eix_pce", start=0.0
    )

    if bounded
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            vlim = _PMs.ref(pm, nw, :bus, xfmr["f_bus"])["vmax"]
            JuMP.set_lower_bound(erx_pce[x, 0], -vlim)
            JuMP.set_upper_bound(erx_pce[x, 0],  vlim)
            JuMP.set_lower_bound(eix_pce[x, 0], -vlim)
            JuMP.set_upper_bound(eix_pce[x, 0],  vlim)
        end
    end

    if report
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
            f_idx = (x, f_bus, t_bus); t_idx = (x, t_bus, f_bus)
            sol_x = _PMs.sol(pm, nw, :xfmr, x)
            sol_x[:vrx_pce_fr] = Dict(k => vrx_pce[f_idx, k] for k in K)
            sol_x[:vix_pce_fr] = Dict(k => vix_pce[f_idx, k] for k in K)
            sol_x[:vrx_pce_to] = Dict(k => vrx_pce[t_idx, k] for k in K)
            sol_x[:vix_pce_to] = Dict(k => vix_pce[t_idx, k] for k in K)
            sol_x[:erx_pce]    = Dict(k => erx_pce[x, k] for k in K)
            sol_x[:eix_pce]    = Dict(k => eix_pce[x, k] for k in K)
        end
    end
end

"""
    variable_xfmr_current_pce(pm::sHHC_QCQP; nw, K, bounded, report)

PCE coefficients of the transformer winding currents (total, series, magnetising).

Creates:
- `crx_pce[(x,i,j), k]`, `cix_pce[(x,i,j), k]` — total winding arc currents
- `csrx_pce[(x,i,j), k]`, `csix_pce[(x,i,j), k]` — series winding arc currents
- `cmrx_pce[x, k]`, `cmix_pce[x, k]` — magnetising currents (set to zero by
  `constraint_xfmr_core_magnetization_pce`; variables kept for interface uniformity)

Bounds (k=0 only): ±c_rating for all arc current variables; k≥1 are unbounded.
"""
function variable_xfmr_current_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                    K=0:2, bounded::Bool=true, report::Bool=true)
    # ── total winding currents ────────────────────────────────────────────────
    crx_pce = _PMs.var(pm, nw)[:crx_pce] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs), k in K],
        base_name="$(nw)_crx_pce", start=0.0
    )
    cix_pce = _PMs.var(pm, nw)[:cix_pce] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs), k in K],
        base_name="$(nw)_cix_pce", start=0.0
    )

    # ── series winding currents ───────────────────────────────────────────────
    csrx_pce = _PMs.var(pm, nw)[:csrx_pce] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs), k in K],
        base_name="$(nw)_csrx_pce", start=0.0
    )
    csix_pce = _PMs.var(pm, nw)[:csix_pce] = JuMP.@variable(pm.model,
        [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs), k in K],
        base_name="$(nw)_csix_pce", start=0.0
    )

    # ── magnetising currents ──────────────────────────────────────────────────
    cmrx_pce = _PMs.var(pm, nw)[:cmrx_pce] = JuMP.@variable(pm.model,
        [x in _PMs.ids(pm, nw, :xfmr), k in K],
        base_name="$(nw)_cmrx_pce", start=0.0
    )
    cmix_pce = _PMs.var(pm, nw)[:cmix_pce] = JuMP.@variable(pm.model,
        [x in _PMs.ids(pm, nw, :xfmr), k in K],
        base_name="$(nw)_cmix_pce", start=0.0
    )

    if bounded
        for (x, i, j) in _PMs.ref(pm, nw, :xfmr_arcs)
            c_rating = _PMs.ref(pm, nw, :xfmr, x)["c_rating"]
            for var in (crx_pce, cix_pce, csrx_pce, csix_pce)
                JuMP.set_lower_bound(var[(x,i,j), 0], -c_rating)
                JuMP.set_upper_bound(var[(x,i,j), 0],  c_rating)
            end
        end
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            c_rating = xfmr["c_rating"]
            JuMP.set_lower_bound(cmrx_pce[x, 0], -c_rating)
            JuMP.set_upper_bound(cmrx_pce[x, 0],  c_rating)
            JuMP.set_lower_bound(cmix_pce[x, 0], -c_rating)
            JuMP.set_upper_bound(cmix_pce[x, 0],  c_rating)
        end
    end

    if report
        for (x, xfmr) in _PMs.ref(pm, nw, :xfmr)
            f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
            f_idx = (x, f_bus, t_bus); t_idx = (x, t_bus, f_bus)
            sol_x = _PMs.sol(pm, nw, :xfmr, x)
            sol_x[:crx_pce_fr]  = Dict(k => crx_pce[f_idx, k] for k in K)
            sol_x[:cix_pce_fr]  = Dict(k => cix_pce[f_idx, k] for k in K)
            sol_x[:csrx_pce_fr] = Dict(k => csrx_pce[f_idx, k] for k in K)
            sol_x[:csix_pce_fr] = Dict(k => csix_pce[f_idx, k] for k in K)
            sol_x[:crx_pce_to]  = Dict(k => crx_pce[t_idx, k] for k in K)
            sol_x[:cix_pce_to]  = Dict(k => cix_pce[t_idx, k] for k in K)
            sol_x[:csrx_pce_to] = Dict(k => csrx_pce[t_idx, k] for k in K)
            sol_x[:csix_pce_to] = Dict(k => csix_pce[t_idx, k] for k in K)
            sol_x[:cmrx_pce]    = Dict(k => cmrx_pce[x, k] for k in K)
            sol_x[:cmix_pce]    = Dict(k => cmix_pce[x, k] for k in K)
        end
    end
end

"""
    variable_lifted_current_pce(pm::sHHC_QCQP; nw, K, report)

Lifted auxiliary variables `j_pce[(l,i,j),k]` for each arc `(l,i,j)` and PCE
index `k ∈ K`.

Semantics: `j_pce[(l,i,j),k]` is the k-th PCE coefficient of the squared branch
current magnitude `|I_{l,h}|² = Cr_{l,h}² + Ci_{l,h}²`. These variables are
linked to the branch current PCE coefficients via SOC (k=0) and quadratic
equality constraints (k≥1). No explicit bounds are imposed.
"""
function variable_lifted_current_pce(pm::sHHC_QCQP; nw::Int=fundamental(pm),
                                      K=0:2, report::Bool=true)
    j_pce = _PMs.var(pm, nw)[:j_pce] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs), k in K], base_name="$(nw)_j_pce",
        start = 0.0
    )

    if report
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            _PMs.sol(pm, nw, :branch, l)[:j_pce_fr] = Dict(k => j_pce[(l,i,j), k] for k in K)
        end
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_to)
            _PMs.sol(pm, nw, :branch, l)[:j_pce_to] = Dict(k => j_pce[(l,i,j), k] for k in K)
        end
    end
end