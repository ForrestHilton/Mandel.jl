import Pkg
Pkg.instantiate()

function complex_sqrt(z)
  # modified branch cut
  r = sqrt(real(z)^2 + imag(z)^2)
  θ = atan(imag(z), real(z))
  if θ < 0
    θ += 2 * pi
  end
  θ = θ / 2
  return sqrt(r) * (cos(θ) + im * sin(θ))
end

f(z::ComplexF64, c::ComplexF64) = z^3 / 3 - c / 2 * z^2 + 3 / 4 * c + 3 * complex_sqrt(c^2 / 16 - 1 / 3)
crit(c) = c
include("./src/MandelMakie.jl")
MandelMakie.Viewer(f; crit=crit, mandel_diameter=1.0)
c = -1.8368098237963426 - 0.3721564996080768im
println(1 / 3, c / 2, 0, 3 / 4 * c, 3 * complex_sqrt(c^2 / 16 - 1 / 3))
