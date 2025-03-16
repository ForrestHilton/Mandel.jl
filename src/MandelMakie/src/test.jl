using Symbolics
@variables z c

f = z^3 + c * z^2 + z
cpts = [-1 / 3 * c - 1 / 3 * sqrt(c^2 - 3), -1 / 3 * c + 1 / 3 * sqrt(c^2 - 3)]

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
    println(new_code)
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

        converged_indices =
            non_converged_indices[abs.(xₙ₊₁s .- xₙs[non_converged_indices]).<solver.abstol]
        converged[converged_indices] .= true
        solutions[converged_indices] .= xₙs[converged_indices]
        xₙs[.!converged] .= xₙ₊₁s[.!converged]
        if all(converged)
            break
        end
    end

    # I don't know if this is necessary
    confirmed_indices =
        findall(abs.(Base.invokelatest(solver.ffunc, solutions)) .< solver.abstol)

    confirmed_solutions = solutions[confirmed_indices]
    filtered_solutions = ComplexF64[]
    for solution in confirmed_solutions
        if isempty(filtered_solutions) ||
           all(abs(solution - fs) > epsilon for fs in filtered_solutions)
            push!(filtered_solutions, solution)
        end
    end

    return filtered_solutions
end

ex1 = c^2
ex2 = 2.0im
eq = Equation(ex1, ex2)
solver = NewtonSolver(eq, c)

result = solve(solver, [Complex(1.2), Complex(-1.2)])
println(result)
