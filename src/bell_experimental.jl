
"""
    bell_poly_step_combi(bp, bc, n)

Computes all Bell polynomials and coefficients for order n.
"""
function bell_poly_step_combi(bp, bc, n)
    J1 = vcat(zeros(Int8,n-1),Int8(1))
    Jn = vcat(Int8(n),zeros(Int8,n-1))

    J_comp = [J1]
    a_comp = [1]
    for k=2:n-1

        J_int = Vector{Int8}(undef,n) 
        a_int = []
        for i=1:(n-k+1)
            J_prev = bp[n-i+1][k-1]
            a_prev = bc[n-i+1][k-1]
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
            
            a_next = a_prev * binomial(big(n-1), big(i-1))
            a_int = vcat(a_int,a_next)
        end
        J_int_uni = unique(J_int,dims=2)
        J_comp = vcat(J_comp, [J_int_uni])
        
        id_bp = map(i-> findall( all(J_int .== J_int_uni[:,i], dims=1)[1,:]), 1:size(J_int_uni)[2] )
        a_int = map(i-> sum(a_int[i]), id_bp)

        a_comp = vcat(a_comp,[a_int])
    end
    J_comp = vcat(J_comp, [Jn])
    a_comp = vcat(a_comp, 1)

    return J_comp, a_comp
end




"""
    bell_poly_coeff_combi(N :: Int64; 
                    bp=[[],[Int8[1]]], 
                    bc=[[],[1]],
                    save_on_disk=false, 
                    path_to_folder="bell_results/", 
                    print_iteration=false)

Computes all Bell polynomials and coefficients up to order of N.

Returns a Vector of N Bell polynomials (saved as vectors and matrices) and coeffiencts.
"""
function bell_poly_coeff_combi(N:: Int64; 
                    bp=[[],[Int8[1]]], 
                    bc=[[],[1]],
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
            J_comp, a_comp = bell_poly_step_combi(bp, bc, n)
            
            if save_on_disk == true
                write_bell_poly(path_to_folder,bp[n],n,max_digit=max_digit)
            end

            bp = vcat(bp, [J_comp])
            bc = vcat(bc, [a_comp])

            if print_iteration == true
                day_now = Dates.Date(Dates.now())
                time_now = Dates.Time(Dates.now())
                iter_now = string(n, pad=max_digit)
                println("Iteration: n = ",iter_now, " Day: ", day_now, " Time: ", time_now)
            end

        end
    
    return bp,bc
end



"""
    bell_poly_coeff_step(bp, n)

Computes all Bell polynomials and coefficients for order n.
"""
function bell_poly_coeff_step(bp, n)
    J1 = vcat(zeros(Int8,n-1),Int8(1))
    Jn = vcat(Int8(n),zeros(Int8,n-1))

    J_comp = [J1]
    a_comp = [1]
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
        a_int = _bell_coeff(J_int,n+1)
        a_comp = vcat(a_comp,[a_int])
    end
    J_comp = vcat(J_comp, [Jn])
    a_comp = vcat(a_comp, 1)

    return J_comp, a_comp
end




"""
bell_poly_coeff(N :: Int64; 
                    bp=[[],[Int8[1]]], 
                    bc=[[],[1]],
                    save_on_disk=false, 
                    path_to_folder="bell_results/", 
                    print_iteration=false)

Computes all Bell polynomials and coefficients up to order of N.

Returns a Vector of N Bell polynomials (saved as vectors and matrices) and coeffiencts.
"""
function bell_poly_coeff(N:: Int64; 
                    bp=[[],[Int8[1]]], 
                    bc=[[],[1]],
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
            J_comp, a_comp = bell_poly_step(bp, bc, n)
            
            if save_on_disk == true
                write_bell_poly(path_to_folder,bp[n],n,max_digit=max_digit)
            end

            bp = vcat(bp, [J_comp])
            bc = vcat(bc, [a_comp])

            if print_iteration == true
                day_now = Dates.Date(Dates.now())
                time_now = Dates.Time(Dates.now())
                iter_now = string(n, pad=max_digit)
                println("Iteration: n = ",iter_now, " Day: ", day_now, " Time: ", time_now)
            end

        end
    
    return bp,bc
end
