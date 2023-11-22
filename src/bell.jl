
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

        J_int = Vector{Int8}(undef,n) 

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

            J_int = (i==1 ? J_next : hcat(J_int, J_next))
        end
        J_int = unique(J_int,dims=2)
        J_comp = vcat(J_comp, [J_int])
    end
    J_comp = vcat(J_comp, [Jn])
    return J_comp
end




"""
bell_poly(N:: Int64; 
            bp=[[],[ones(Int8,1)]], 
            save_on_disk=false,     
            path_to_folder="bell_results/", 
            print_iteration=false)

Computes all Bell polynomials up to order of N.

Returns a Vector of N Bell polynomials (saved as vectors and matrices).
"""
function bell_poly(N:: Int64; 
                    bp=[[],[ones(Int8,1)]], 
                    save_on_disk=false, 
                    path_to_folder="bell_results/", 
                    print_iteration=false)

        max_digit = length(digits(N))

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
    
        for n=n_start:N
            J_comp = bell_poly_step(n, bp)
            
            if save_on_disk == true
                write_bell_poly(path_to_folder,J_comp,n,max_digit=max_digit)
                # write_bell_poly(path_to_folder,bp[n],n,max_digit=max_digit)
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


function _bell_coeff(mat :: Matrix{<: Integer},n,k)
    res = Int64[];
    for col in eachcol(mat)
        res = vcat(res, _bell_coeff(col[:],n,k))
    end
    return res;
end

function _bell_coeff(v :: Vector{<: Integer},n,k)
    a1 = mapreduce(i-> factorial(big(v[i]))*factorial(big(i))^v[i],*, 1:n-k)
    a2 = factorial(big(n-1))*inv(a1)
    return round(BigInt, a2)
end

"""
    bell_coeff(bp)

Computes and returns the coefficients of all Bell polynomials saved in `bp`.
"""
function bell_coeff(bp)
    bc = [];

    for (n1, e1) in enumerate(bp)
        bc_out = [];
        for (n2, e2) in enumerate(e1)
            bc_in = _bell_coeff(e2,n1,n2)
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
