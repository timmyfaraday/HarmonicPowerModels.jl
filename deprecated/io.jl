function plot_harmonic_current_injections(result, h, path)
    res = result["solution"]["nw"]["$h"]
    _EE.figure()
    _EE.rc("text", usetex = true)
    for (l, load) in res["load"]
        hc = load["crd"] + load["cid"]im
        label = join(["\$I_{L","$l","}\$"])
        _EE.phasor(hc, label = label, labeltsep = abs(hc)/100, headwidth = 1, headlength = 5)
    end
    plot_path = joinpath(path, join(["current_phasors_h_", "$h",".pdf"]))
    _EE.savefig(plot_path)
    _EE.close("all")
end

function plot_harmonic_voltages(result, h, path)
    res = result["solution"]["nw"]["$h"]
    _EE.figure()
    _EE.rc("text", usetex = true)
    for (b, bus) in res["bus"]
        hv = bus["vr"] + bus["vi"]im
        label = join(["\$U_{","$b","}\$"])
        _EE.phasor(hv, label = label, labeltsep = abs(hv)/100, headwidth = 2, headlength = 5)
    end
    plot_path = joinpath(path, join(["voltage_phasors_h_", "$h",".pdf"]))
    _EE.savefig(plot_path)
    _EE.close("all")
end