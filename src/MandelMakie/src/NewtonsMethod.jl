using Symbolics

@variables z c

function symbolic_to_function(f)
    txt = string(build_function(f, c))
    code = replace(txt, "NaNMath." => "")
    lines = split(code, '\n')
    lines[1] = "function (c::Vector{ComplexF64})"
    lines[5] = "    return @. " * lstrip(lines[5])
    new_code = """"""
    for line in lines
        if startswith(lstrip(line), "#")
            continue
        end
        new_code = new_code * line * "\n"
    end
    return eval(Meta.parse(new_code))
end
p = z^3 + c * z^2 + z
a = -1 / 3 * c - 1 / 3 * sqrt(c^2 - 3)
b = -1 / 3 * c + 1 / 3 * sqrt(c^2 - 3)

# pfunc = symbolic_to_function(p)
function pfunc(z, c)
    return @. z^3 + c * z^2 + z
end
function afunc(c)
    return @. -1 / 3 * c - 1 / 3 * sqrt(c^2 - 3)
end
function bfunc(c)
    return @. -1 / 3 * c + 1 / 3 * sqrt(c^2 - 3)
end

function pp_z_func(z, c)
    return @. 3z^2 + 2z * c + 1
end
function pp_c_func(z) # derivative with respect to c
    return @. z^2
end
function apfunc(c) # derivative with respect to c
    return @. -0.3333333333333333 + (-0.3333333333333333c) / sqrt(-3 + c^2)
end
function bpfunc(c) # derivative with respect to c
    return @. -0.3333333333333333 + (0.3333333333333333c) / sqrt(-3 + c^2)
end

function _iterate(z, c, n)
    # evaluates composition
    ret = z
    for _ in 1:n
        ret = pfunc(ret, c)
    end
    return ret
end

function comp_prime_c(zv, cv, n)
    # derivative of composition n times
    f_c = 1

    z = zv
    for _ in 1:n
        f_c = @. f_c * pp_z_func(z, cv) + pp_c_func(z)
        z = pfunc(z, cv)
    end

    return f_c
end

function comp_prime_z(zv, cv, n)
    # derivative of composition n times
    ret = 1

    z = zv
    for _ in 1:n
        ret = @. ret * pp_z_func(z, cv) # Multiply by derivative of f at current z
        z = pfunc(z, cv) # Update z to f(z)
    end

    return ret
end

function value(cv, n, m)
    return _iterate(afunc(cv), cv, n) .- _iterate(bfunc(cv), cv, m)
end

function derivative(cv, n, m)
    rhs = comp_prime_z(afunc(cv), cv, n) .* apfunc(cv) .+ comp_prime_c(afunc(cv), cv, n)
    lhs = @. comp_prime_z(bfunc(cv), cv, m) * bpfunc(cv) + comp_prime_c(bfunc(cv), cv, m)
    return @. rhs - lhs
end

function solve(n, m, x₀s; abstol = 1e-8, maxiters = 50)
    x₀s = collect(x₀s)
    xₙs = copy(x₀s)
    converged = falses(length(x₀s))

    solutions = similar(xₙs) # Preallocate space for solutions
    for _ in 1:maxiters
        non_converged_indices = findall(.!converged)
        xₙ₊₁s =
            xₙs[non_converged_indices] .-
            value(xₙs[non_converged_indices], n, m) ./
            derivative(xₙs[non_converged_indices], n, m)

        new_converged_indices =
            non_converged_indices[abs.(xₙ₊₁s .- xₙs[non_converged_indices]).<abstol]
        converged[new_converged_indices] .= true
        solutions[new_converged_indices] .= xₙs[new_converged_indices]
        xₙs[non_converged_indices] .= xₙ₊₁s
        if all(converged)
            break
        end
    end

    # I don't know if this is necessary
    confirmed_indices = findall(abs.(value(solutions, n, m)) .< abstol)
    confirmed_solutions = solutions[confirmed_indices]
    # confirmed_solutions = solutions

    filtered_solutions = ComplexF64[]
    for solution in confirmed_solutions
        if isempty(filtered_solutions) ||
           all(abs(solution - fs) > 2e-8 for fs in filtered_solutions)
            push!(filtered_solutions, solution)
        end
    end

    return filtered_solutions
end

function search_rectangle(n, m, z1::ComplexF64, z2::ComplexF64)
    # Extract real and imaginary parts of the corners
    x1, y1 = real(z1), imag(z1)
    x2, y2 = real(z2), imag(z2)

    # Ensure x1 is the leftmost point and y1 is the bottommost point
    x_min, x_max = min(x1, x2), max(x1, x2)
    y_min, y_max = min(y1, y2), max(y1, y2)

    # Generate a grid of points within the rectangle
    num_points_per_side = ceil(Int, sqrt(100)) # Adjust this to change the total number of points
    real_points = range(x_min, stop = x_max, length = num_points_per_side)
    imag_points = range(y_min, stop = y_max, length = num_points_per_side)
    grid = [Complex(x, y) for x in real_points for y in imag_points]
    # grid = [(z1 + z2) / 2]

    # Use the grid points as initial guesses for the Newton's method solver
    solutions = solve(n, m, grid)

    return solutions
end

# println(value(ComplexF64(1), 1, 1))
# 0.0 + 0.6285393610547088im -- expected

# ex1 = c^2
# ex2 = 2.0im
# eq = Equation(ex1, ex2)
# solver = NewtonSolver(eq, c)
# -0.4855697078611445 - 1.0345273182409185im
# -0.4823432393234249 - 1.0313008497031988im
# result = search_rectangle(2, 2, Complex(-0 - 0im), Complex(2.0 + 2.0im))
