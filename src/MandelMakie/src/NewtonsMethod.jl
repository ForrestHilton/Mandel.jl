using Symbolics

@variables z c
f = z^3 + c * z^2 + z
function symbolic_iterate(ex, n)
    for _ in 1:n
        ex = substitute(f, Dict([z => ex]))
    end
    return ex
end

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

struct NewtonSolver
    eq::Any
    x::Any
    f::Any
    f′::Any
    ffunc::Any
    fpfunc::Any
    abstol::Float64
    maxiters::Int

    function NewtonSolver(eq, x; abstol = 1e-8, maxiters = 50)
        f_ns = eq.lhs - eq.rhs
        ffunc = symbolic_to_function(f_ns)
        f′ = Symbolics.derivative(f_ns, x)
        fpfunc = symbolic_to_function(f′)

        new(eq, x, f_ns, f′, ffunc, fpfunc, abstol, maxiters)
    end
end

function solve(solver, x₀s)
    x₀s = collect(x₀s)
    xₙs = copy(x₀s)
    converged = falses(length(x₀s))

    solutions = similar(xₙs) # Preallocate space for solutions
    for _ in 1:solver.maxiters
        non_converged_indices = findall(.!converged)
        xₙ₊₁s =
            xₙs[non_converged_indices] .-
            Base.invokelatest(solver.ffunc, xₙs[non_converged_indices]) ./
            Base.invokelatest(solver.fpfunc, xₙs[non_converged_indices])

        new_converged_indices =
            non_converged_indices[abs.(xₙ₊₁s .- xₙs[non_converged_indices]).<solver.abstol]
        converged[new_converged_indices] .= true
        solutions[new_converged_indices] .= xₙs[new_converged_indices]
        xₙs[non_converged_indices] .= xₙ₊₁s
        if all(converged)
            break
        end
    end

    # I don't know if this is necessary
    # confirmed_indices =
    #     findall(abs.(Base.invokelatest(solver.ffunc, solutions)) .< solver.abstol)
    #
    # confirmed_solutions = solutions[confirmed_indices]
    filtered_solutions = ComplexF64[]
    confirmed_solutions = solutions
    for solution in confirmed_solutions
        if isempty(filtered_solutions) ||
           all(abs(solution - fs) > 2e-8 for fs in filtered_solutions)
            push!(filtered_solutions, solution)
        end
    end

    return filtered_solutions
end

function search_rectangle(solver, z1::ComplexF64, z2::ComplexF64)
    # Extract real and imaginary parts of the corners
    println(solver.eq, z1, z2)
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

    # Use the grid points as initial guesses for the Newton's method solver
    solutions = solve(solver, grid)
    println(solutions)

    return solutions
end

# ex1 = c^2
# ex2 = 2.0im
# eq = Equation(ex1, ex2)
# solver = NewtonSolver(eq, c)
#
# result = search_rectangle(solver, Complex(-1.2 - 1.2im), Complex(1.2 + 1.2im))
nothing
