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

function solve_newton(eq, x, x₀; abstol = 1e-8, maxiters = 50)
    # symbolic expressions for f(x) and f′(x)
    f = eq.lhs - eq.rhs # want to find root of f(x)
    f′ = Symbolics.derivative(f, x)

    xₙ = x₀ # numerical value of the initial guess
    for i in 1:maxiters
        # calculate new guess by numerically evaluating symbolic expression at previous guess
        xₙ₊₁ = substitute(x - f / f′, x => xₙ)
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁ # converged
        else
            xₙ = xₙ₊₁
        end
    end
    error("Newton's method failed to converge")
end

function print_many()
    for i in -50:50
        println(solve_newton(eq, c, Complex(0.5 + (i / 25 * 1im))))
    end
end
print_many()

# print(solve(ex1-ex2,c))
# Sym{PythonCall.Core.Py}[-1.73205080756888, 1.73205080756888, -1.85790551759459*I, 1.8579055175
# 9459*I, -1.76126765566016 - 0.501155962471181*I, -1.76126765566016 + 0.501155962471181*I, 1.76
# 126765566016 - 0.501155962471181*I, 1.76126765566016 + 0.501155962471181*I]
