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
y = rand(101)
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


# ===== find nearest point =====
function findNearestPoint(path_x, path_y, self_x, self_y)
    dist_squared = (path_x .- self_x).^2 + (path_y .- self_y).^2
    (min_dist_squared, min_index) = findmin(dist_squared)
    return (sqrt(min_dist_squared), min_index)
end


# ===== trajectory design =====
x = 0:1:10
y = sin.(x)
scatter(x, y)
# dy/dx = cos(x)からthetaを求める
dydx = cos.(x) # = tan(theta)
theta = atan.(dydx)
scatter(x, theta)


# ===== robot model =====
#
# -- continuous --
# dot x     = v * cos (theta)
# dot y     = v * sin (theta)
# dot theta = w
#
# X = [x, y, theta]'
# U = [v, w]'.
#
#  -- discrete --
# x_k+1     = x_k + v_k * cos(theta_k) * dt
# y_k+1     = y_k + v_k * sin(theta_k) * dt
# theta_k+1 = theta_k + w_k * dt
#
# dF(X,U)/dX = [1, 0, -v*sin(theta)*dt;
#               0, 1, v*cos(theta)*dt;
#               0, 0, 1]
#
# dF(X,U)/dU = [cos(theta)*dt, 0;
#               sin(theta)*dt, 0;
#               0,            dt]



# ==== test ====
self_x = 3
self_y = 3

(min_dist, i) = findNearestPoint(x, y, self_x, self_y)
