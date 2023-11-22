# Computing the iterations of bell_poly

function sum_poly(x)
    val = 0
    for k=2:x-1
        for i=1:x-k+1
            val = val+1;
        end
    end
    return val
end

N_bell = 40;
cnt_poly = zeros(Int64,40);
@. cnt_poly = sum_poly(1:N_bell)

M = vcat(map(x-> [x^2 x 1], 10:10:30)...)
y = cnt_poly[10:10:30]
c = round.(inv(M)*y,digits=1)

f(x) = c[1]*x^2 + c[2]*x + c[3]

f.(10:5:40)