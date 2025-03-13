using Symbolics
@variables c

f = c^3

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
    global new_code = new_code * line * "\n"
end

println(new_code)
# comparison = """
# function (c::Vector{ComplexF64})
#     return @. (^)(c, 3)
# end
# """
# println(length(comparison))
# println(length(new_code))
# println(new_code == comparison)

# ffunc = eval(Meta.parse(comparison))
# println(ffunc([Complex(1.0)]))

ffunc2 = eval(Meta.parse(new_code))
println(ffunc2([Complex(1.0), Complex(2.0)]))
