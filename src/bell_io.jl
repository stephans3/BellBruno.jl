function write_bell_poly(folder_path :: String, J_comp,n;max_digit=2)
    num = string(n, pad=max_digit)
    for (i, elem) in enumerate(J_comp)
        num_i = string(i, pad=max_digit)
        path_bp = folder_path*string("bp-",num,"-",num_i,"-f.txt")
        open(path_bp, "w") do io
            writedlm(io, elem)
        end
    end
end


function read_bell_poly(;path_to_folder="bell_results/", bp=[[],[ones(Int8,1)]] )
    
    if endswith(path_to_folder, "/") == false
        path_to_folder = path_to_folder *"/"
    end

    flist = readdir(path_to_folder)
    bp_i1 = findfirst(isequal(1),startswith.(flist,"bp"))
    bp_i2 = findlast(isequal(1),startswith.(flist,"bp"))
    bp_files = flist[bp_i1:bp_i2]

    bp_outer = bp
    bp_inner = []

    bp_idx = -1;
    for (i, elem) in enumerate(bp_files)
        path = path_to_folder*elem
        data = readdlm(path, '\t', Int8)
        
        bp_idx_now = parse(Int, split(elem,"-")[2])
        if i==1
            bp_inner = [data]
            bp_idx = bp_idx_now
        elseif bp_idx != bp_idx_now 
            bp_outer = vcat(bp_outer, [bp_inner])
            bp_inner = [data]
            bp_idx = bp_idx_now
        else
            bp_inner = vcat(bp_inner, [data])
        end
    end
    bp_outer = vcat(bp_outer, [bp_inner])
    return bp_outer;
end
