# Pkg.add("Pkg") # not work in v0.6 ?
# using Pkg # not work in v0.6 ?
# Pkg.add("IJulia")
# Pkg.add("PyCall")
# Pkg.add("PyPlot")
# Pkg.add("Gadfly")
Pkg.add("BenchmarkTools") # ベンチマーク
Pkg.add("Plots")
Pkg.add("Clp") # 線形計画法
Pkg.add("JuMP") # 最適化パーサー: supported only up to v0.6

using Plots
# using LinearAlgebra # used as default in v0.6, needed only for v1.0
using BenchmarkTools
# set_default_plot_size(9inch, 9inch/golden);
t = 0:0.1:10
x = sin.(t)
y = rand(101)
plot(x, y, linewidth=2, title="My Plot")
plot!(x, y.-1)
# scatter(x, y)




A = [1. 2.; 3. 4.]
B = [5 6; 7 8]
C = A * B
eigvals(A)
eigvecs(A)
X = rand(2,3)
Imat = X * pinv(X)

println(C)
println(D)


# ===== find nearest point =====
function findNearestPoint(path_x, path_y, self_x, self_y)
    dist_squared = (path_x .- self_x).^2 + (path_y .- self_y).^2
    (min_dist_squared, min_index) = findmin(dist_squared)
    return (sqrt(min_dist_squared), min_index)
end


# ===== find nearest point =====
function findNearestPointQ(path_x, path_y, self_x, self_y)
    dist_squared = zeros(length(path_x))
    @simd for i=1:length(path_x)
        @inbounds dist_squared[i] = (path_x[i] - self_x)^2 + (path_y[i] - self_y)^2
    end
    (min_dist_squared, min_index) = findmin(dist_squared)
    return (sqrt(min_dist_squared), min_index)
end

# ===== set params =====
dim_x = 3
dim_u = 2

# # ===== trajectory design 1 =====
# x = 0:0.1:10
# y = sin.(x)
# scatter(x, y)
# # dy/dx = cos(x)からthetaを求める
# dydx = cos.(x) # = tan(theta)
# theta = atan.(dydx)
# scatter(x, theta)

# ===== trajectory design 2 =====
dt = 0.1
t = 0:dt:10
T = 5.0
traj_v = 0.2 + 0.1 * sin.(2.0 * pi * t / T)
traj_w = 0.5 * cos.(2.0 * pi * t / T)
fig1 = plot(t, traj_v, linewidth=2, label="v(t)")
plot!(t, traj_w, linewidth=2, label="w(t)")
traj_x = zeros(length(t))
traj_y = zeros(length(t))
traj_theta = zeros(length(t))
traj_x[1] = 0.
traj_y[1] = 0.
traj_theta[1] = 0.
for i = 2:length(t)
    traj_x[i] = traj_x[i-1] + traj_v[i-1]*cos(traj_theta[i-1])*dt
    traj_y[i] = traj_y[i-1] + traj_v[i-1]*sin(traj_theta[i-1])*dt
    traj_theta[i] = traj_theta[i-1] + traj_w[i-1]*dt
end
traj_X = zeros(length(t), 3)
traj_X[:,1] = traj_x
traj_X[:,2] = traj_y
traj_X[:,3] = traj_theta
traj_U = zeros(length(t), 2)
traj_U[:,1] = traj_v
traj_U[:,2] = traj_w
# fig2 = scatter(x,y)

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
# dF(X,U)/dX = [1, 0, -v_k*sin(theta_k)*dt;
#               0, 1,  v_k*cos(theta_k)*dt;
#               0, 0, 1]
#
# dF(X,U)/dU = [cos(theta)*dt, 0;
#               sin(theta)*dt, 0;
#               0,            dt]

# ===== linearize A & B around trajectory =====
linAs = zeros(length(t), dim_x, dim_x)
linBs = zeros(length(t), dim_x, dim_u)
linAs[:,1,1] = 1.
linAs[:,2,2] = 1.
linAs[:,3,3] = 1.
linAs[:,1,3] = -traj_v.*cos.(traj_theta)*dt
linAs[:,2,3] = traj_v.*sin.(traj_theta)*dt
linBs[:,1,1] = cos.(traj_theta)*dt
linBs[:,2,1] = sin.(traj_theta)*dt
linBs[:,3,2] = dt


# ===== calc nearest point, set Xref & Uref =====
N = 3
self_x = 0.
self_y = 0.
self_theta = 0.
self_X = [self_x, self_y, self_theta]
dist, idx = findNearestPointQ(traj_x, traj_y, self_x, self_y)
pred_num = min(N, length(traj_x)-idx-N-1)
Xref = zeros(0)
Uref = zeros(0)
for i = 0:pred_num-1
    append!(Xref, traj_X[idx+i+1,:])
    append!(Uref, traj_U[idx+i,:])
end



# ===== set Aex, Bex, Q, R =====
Aex = zeros(0)
Bex = zeros(0)
Aprev = eye(dim_x)
for i=0:pred_num-1
    Atmp = linAs[idx+i, :,:] * Aprev
    Aex = vcat(Aex, Atmp)
    Aprev = Atmp

    Btmp = zeros(dim_x, N*dim_u)
    for j=0:i
        tmp = linBs[idx+j,:,:]
        for k=j+1:i
            tmp = linAs[idx+k,:,:] * tmp
        end
        Btmp[:, j*dim_u+1:(j+1)*dim_u] = tmp
    end
    Bex = vcat(Bex, Btmp)
end

# ===== set Q & R for MPC =====
Qex = eye(dim_x*pred_num)
Rex = eye(dim_u*pred_num)
mat1 = inv(Bex' * Qex * Bex + Rex)
mat2 = (self_X'*Aex'*Qex*Bex - Xref'*Qex*Bex - Uref'*Rex)'
Umpc = -inv(mat1) * mat2


# ===== cost function J =====

Xmpc = Aex * self_X  + Bex * Umpc
J = (Xmpc - Xref)'*Qex*(Xmpc - Xref) + (Umpc-Uref)'*Rex*(Umpc-Uref)

Xmpc2 = Aex * self_X  + Bex * Uref
J = (Xmpc2 - Xref)'*Qex*(Xmpc2 - Xref) + (Uref-Uref)'*Rex*(Uref-Uref)


# ==== test ====
self_x = 3
self_y = 3

a = @benchmark findNearestPoint(x, y, self_x, self_y)
b = @benchmark findNearestPointQ(x, y, self_x, self_y)

println(a)
println(b)


for i=0:-1:-1
    println(i)
end
