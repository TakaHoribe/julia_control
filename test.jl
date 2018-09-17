using Pkg
# Pkg.add("IJulia")
# Pkg.add("PyCall")
# Pkg.add("PyPlot")
# Pkg.add("Gadfly")
Pkg.add("Plots")

using Plots
using LinearAlgebra
# set_default_plot_size(9inch, 9inch/golden);
t = 0:0.1:10
x = sin.(t)
y = rand(10)
plot(x, y, linewidth=2, title="My Plot")
scatter(x, y)



a = 1
b = 2
println(a+b)
c=4


function MATX(x, y)
  return x * y
end


a = 3
b = 2
c = MATX(a,b)
println("c = ", c)


A = [1 2; 3 4]
B = [5 6; 7 8]
C = A * B
D = MATX(A, B)
eigvals(A)
eigvecs(A)
X = rand(2,3)
I = X * pinv(X)

println(C)
println(D)
