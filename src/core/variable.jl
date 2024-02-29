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
################################################################################

# bus
""
function variable_bus_voltage_real(pm::_PMs.AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
function variable_bus_voltage_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
function variable_branch_current_real(pm::_PMs.AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
function variable_branch_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
function variable_branch_series_current_real(pm::_PMs.AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
function variable_branch_series_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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

    ## bounds are needed

    report && _PMs.sol_component_value(pm, nw, :xfmr, :cmrx, _PMs.ids(pm, nw, :xfmr), cmrx)
end
""
function variable_xfmr_current_magnetizing_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cmix = _PMs.var(pm, nw)[:cmix] = JuMP.@variable(pm.model,
            [x in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_cmix",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, x), "cmix_start", 0.0)
    )

    ## bounds are needed

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