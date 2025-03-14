using Symbolics
@variables z c

f = f = z^3 + c * z^2 + z
cpts = [-1 / 3 * c - 1 / 3 * sqrt(c^2 - 3), -1 / 3 * c + 1 / 3 * sqrt(c^2 - 3)]

n = 2
m = 2
function symbolic_iterate(ex, n)
    for _ in 1:n
        ex = substitute(f, Dict([z => ex]))
    end
    return ex
end
ex1 = symbolic_iterate(cpts[1], n)
ex2 = symbolic_iterate(cpts[2], m)

eq = ex1 ~ ex2

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
    # println(new_code)
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
        f = eq.lhs - eq.rhs
        ffunc = symbolic_to_function(f)
        f′ = Symbolics.derivative(f, x)
        fpfunc = symbolic_to_function(f′)

        new(eq, x, f, f′, ffunc, fpfunc, abstol, maxiters)
    end
end

function solve(solver, x₀s)
    x₀s = collect(x₀s)
    xₙs = copy(x₀s)

    for _ in 1:solver.maxiters
        xₙ₊₁s =
            Complex(1) .-
            Base.invokelatest(solver.ffunc, xₙs) ./ Base.invokelatest(solver.fpfunc, xₙs)

        converged = abs.(xₙ₊₁s .- xₙs) .< solver.abstol
        if all(converged)
            return xₙ₊₁s
        else
            xₙs[.!converged] .= xₙ₊₁s[.!converged]
        end
    end
    error("Newton's method failed to converge for some initial guesses")
end

solver = NewtonSolver(eq, c)
initial_guesses = [Complex(0.5 + (i / 25 * 1im)) for i in -50:50]
println(solve(solver, initial_guesses[1:2]))
results = @time solve(solver, initial_guesses)

println.(results)

# print(solve(ex1-ex2,c))
# Sym{PythonCall.Core.Py}[-1.73205080756888, 1.73205080756888, -1.85790551759459*I, 1.8579055175
# 9459*I, -1.76126765566016 - 0.501155962471181*I, -1.76126765566016 + 0.501155962471181*I, 1.76
# 126765566016 - 0.501155962471181*I, 1.76126765566016 + 0.501155962471181*I]
