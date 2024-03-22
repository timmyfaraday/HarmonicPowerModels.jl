function csv_export(result, filename)
    csvfile = open(filename,"w")
    write(csvfile, "Harmonic,", "Node,", "vr,", "vi,", "er,", "ei," ,"vr_fr,", "vr_to,", "vi_fr,", "vi_to,", "type", "\n")
    for h in sort(parse.(Int, keys(result["solution"]["nw"])))
        for b in sort(parse.(Int, keys(result["solution"]["nw"]["$h"]["bus"])))
            vr = result["solution"]["nw"]["$h"]["bus"]["$b"]["vr"]
            vi = result["solution"]["nw"]["$h"]["bus"]["$b"]["vi"]
            er = 0
            ei = 0
            vr_fr = 0
            vr_to = 0
            vi_fr = 0
            vi_to = 0
            line = (h, b, vr, vi, ei, er, vr_fr, vr_to, vi_fr, vr_to, "bus")
            write(csvfile, join(line, ","), "\n")
        end
        buses = length(sort(parse.(Int, keys(result["solution"]["nw"]["$h"]["bus"]))))
        for x in sort(parse.(Int, keys(result["solution"]["nw"]["$h"]["xfmr"])))
            b = buses + (x - 1)
            vr = 0
            vi = 0
            er = result["solution"]["nw"]["$h"]["xfmr"]["$x"]["erx"]
            ei = result["solution"]["nw"]["$h"]["xfmr"]["$x"]["eix"]
            vr_fr = result["solution"]["nw"]["$h"]["xfmr"]["$x"]["vrx_fr"]
            vr_to = result["solution"]["nw"]["$h"]["xfmr"]["$x"]["vrx_to"]
            vi_fr = result["solution"]["nw"]["$h"]["xfmr"]["$x"]["vix_fr"]
            vi_to = result["solution"]["nw"]["$h"]["xfmr"]["$x"]["vix_to"]
            line = (h, b, vr, vi, ei, er, vr_fr, vr_to, vi_fr, vr_to, "xfmr")
            write(csvfile, join(line, ","), "\n")
        end
    end
    close(csvfile)
end