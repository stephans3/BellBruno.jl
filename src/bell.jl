function hasColumn(vec1::Vector{T}, vec2::Vector{T}) where T <:Real
   
    if vec1 == vec2
        return true
    end

    return false
end

function hasColumn(mat::Matrix{T}, vec::Vector{T}) where T <:Real
    for jdx in axes(mat, 2)
        if mat[:,jdx] == vec
            return true
        end
    end

    return false
end


"""
    bell_poly_step(n, bp)

Computes all Bell polynomials up to order of N.

Returns the recent partial Bell polynomial as .
"""
function bell_poly_step(n, bp)
    J1 = vcat(zeros(Int8,n-1),Int8(1))
    Jn = vcat(Int8(n),zeros(Int8,n-1))

    J_comp = [J1]
    
    for k=2:n-1

        J_int = Vector{Int8}(undef,n) # Vector{Int8}(undef,8) # zeros(Int8,N,0)

        for i=1:(n-k+1)
            J_prev = bp[n-i+1][k-1]
            
            J_next = similar(J1)
            m=length(J_prev[1,:])

            if m==1
                L = vcat(zeros(Int8,i-1),Int8(1),zeros(Int8,n-i))
                J_next = vcat(J_prev,zeros(Int8,i)) + L
            else
                L = vcat(zeros(Int8,i-1,m),ones(Int8,1,m),zeros(Int8,n-i,m))
                J_next = vcat(J_prev,zeros(Int8,i,m)) + L
            end
            
            if i==1
                J_int = J_next
            else
                for m1 = 1 : m
                    if hasColumn(J_int, J_next[:,m1]) == false
                        J_int = hcat(J_int, J_next[:,m1])
                    end
                end
            end
        end
        J_comp = vcat(J_comp, [J_int])
    end
    J_comp = vcat(J_comp, [Jn])
    return J_comp
end




"""
bell_poly(Nmax :: Int64; 
            bp=[[],[ones(Int8,1)]], 
            save_on_disk=false,     
            path_to_folder="bell_results/", 
            print_iteration=false)

Computes all Bell polynomials up to order of Nmax.

Returns a Vector of N Bell polynomials (saved as vectors and matrices).
"""
function bell_poly(Nmax :: Int64; 
                    bp=[[],[ones(Int8,1)]], 
                    save_on_disk=false, 
                    path_to_folder="bell_results/", 
                    print_iteration=false)

        max_digit = length(digits(Nmax))

        if save_on_disk == true 
            
            if endswith(path_to_folder, "/") == false
                path_to_folder = path_to_folder *"/"
            end

            if isdir(path_to_folder) == false
                mkdir(path_to_folder)    # Make folder
            end

            println("Files are saved in folder: ", path_to_folder, "\n")

        end

        n_start = length(bp) 
    
        for n=n_start:Nmax
            J_comp = bell_poly_step(n, bp)
            
            if save_on_disk == true
                write_bell_poly(path_to_folder,J_comp,n,max_digit=max_digit)
            end

            bp = vcat(bp, [J_comp])

            if print_iteration == true
                day_now = Dates.Date(Dates.now())
                time_now = Dates.Time(Dates.now())
                iter_now = string(n, pad=max_digit)
                println("Iteration: n = ",iter_now, " Day: ", day_now, " Time: ", time_now)
            end

        end
    
    return bp
end




function bell_poly_original(N::Int64)
    
    Minit = vcat(Int8(1),zeros(Int8,N-1))
    Mbell = [[zeros(Int8,N)], [Minit]]
    
    
    for n=2:N
        J1 = vcat(zeros(Int8,n-1),Int8(1),zeros(Int8,N-n))
        Jn = vcat(Int8(n),zeros(Int8,N-1))
    
        J_comp = [J1]
        
        for k=2:n-1
    
            J_int = Vector{Int8}(undef,8) # zeros(Int8,N,0)
    
            for i=1:(n-k+1)
                J_prev = Mbell[n-i+1][k-1]
                m=1
                if typeof(J_prev) != typeof(Mbell[1][1])
                    m = size(J_prev)[2]
                    L = vcat(zeros(Int8,i-1,m),ones(Int8,1,m),zeros(Int8,N-i,m))
                else
                    L = vcat(zeros(Int8,i-1),Int8(1),zeros(Int8,N-i))
                end
         
                J_next = J_prev + L
    
                if i==1
                    J_int = J_next
                else
                    for m1 = 1 : m
                        if hasColumn(J_int, J_next[:,m1]) == false
                            J_int = hcat(J_int, J_next[:,m1])
                        end
                    end
                end
            end
            J_comp = vcat(J_comp, [J_int])
        end
    
        J_comp = vcat(J_comp, [Jn])
        Mbell = vcat(Mbell, [J_comp])
        
    end
    
    return Mbell
end

function _bell_coeff(mat :: Matrix{Int8})
    res = Int64[];

    for col in eachcol(mat)
        res = vcat(res, _bell_coeff(col[:]))
    end

    return res;
end

function _bell_coeff_original(vec :: Vector{Int8})
    res = 1;

    N = length(vec)
    k = sum(vec)
    n = collect(1:N)' * vec

    for idx in 1:(n-k+1)
        elem = vec[idx]
        res = res * factorial(elem)*factorial(idx)^elem
    end

    return round(Int64,factorial(n)/res);
end

function _bell_coeff(vec :: Vector{Int8})
    res = BigInt(1);

    N = length(vec)
    #k = sum(vec)
    n = collect(1:N)' * vec

    for (idx, elem) in enumerate(vec)
        fac_idx = factorial(big(idx));
        fac_elem = factorial(big(elem));
        res = res * fac_elem * fac_idx^elem
    end

    fac_n = factorial(big(n))    
    return round(BigInt,fac_n/res)
end


"""
    bell_coeff(bp)

Computes and returns the coefficients of all Bell polynomials saved in `bp`.
"""
function bell_coeff(bp)
    bc = [];

    for (_, e1) in enumerate(bp)
        bc_out = [];
        for (_, e2) in enumerate(e1)
            bc_in = _bell_coeff(e2)
            bc_out = vcat(bc_out, [bc_in])
        end
    
        bc = vcat(bc, [bc_out])
    end

    return bc
end

"""
    reduce_bell_poly(bp; max_order = 2 :: Int64)

Removes all monomials from a Bell polynomial with an higher order than `max_order`.

Returns the reduced vector of Bell polynomials.
"""
function reduce_bell_poly(bp; max_order = 2 :: Int64)

    bp_red = []
   
    for (idx, elem) in enumerate(bp)
        del_elems = Int64[];
        elem_cat = hcat(elem...)
        for jdx in axes(elem_cat, 2)
            if sum(elem_cat[max_order+1:end,jdx]) > 0
                del_elems = vcat(del_elems, jdx)
            end
        end
        mat_reduced = elem_cat[1:max_order, 1:end .∉ [del_elems]]
        bp_red = vcat(bp_red, [mat_reduced]) 
        
    end

    return bp_red
end

"""
    note_bell_poly(bp, index)

Returns the Bell polynomial of `index` as a string. 
"""
function note_bell_poly(bp, index)
    p = bp
    n=index;
    bn_poly = hcat(p[n+1]...)
    cn = _bell_coeff(bn_poly)
    term_text = "B"*string(n);
    bn_poly_nz = findall(!iszero, bn_poly);
    col_old = 0;
    for (idx, elem) in enumerate(bn_poly_nz)
        row = elem[1];
        col = elem[2];
        c = cn[col]
        b_exp = bn_poly[row,col]
        if col_old < col
            if idx == 1
                term_text = term_text * " = " * string(c)
            else    
                term_text = term_text * " + " * string(c)
            end
            col_old = col;
        end
        
        if b_exp > 1
            term_text = term_text * " x" * string(row) * "^" * string(b_exp); 
        else
            term_text = term_text * " x" * string(row);
        end
    end
    return term_text
end
