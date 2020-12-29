import Pkg; Pkg.activate(joinpath(@__DIR__,".."))
using Test
using OptiViz
using Makie
using StaticArrays, LinearAlgebra
using ForwardDiff 

include("../src/newton.jl")

f(x) = ((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2) / 1000
gradf(x) = ForwardDiff.gradient(f, x)
hessf(x) = ForwardDiff.hessian(f, x)


# newton(f, gradf, hessf, x0)

N = 30
function xy_data(x, y)
    r = sqrt(x^2 + y^2)
    r == 0.0 ? 1f0 : (sin(r)/r)
end
lspace = range(-3, stop = 3, length = N)
z = Float64[f([x,y]) for x in lspace, y in lspace]
r = range(0, stop = 3, length = N)
r = lspace
surface!(scene,
    r, r, z,
    colormap = :Spectral,
)

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

scene = Scene()
lim = FRect3D((0,0,0.), (50,50,1))
surface!(scene,xs, ys, zs, scale=(1,1,1))

x = Node(2.0)
y = Node(1.0)
p = @lift Point3($x,$y,f([$x, $y]))
scatter!(scene, p, markersize=100)
x[] = 1 
y[] = 1