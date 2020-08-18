using Plots
using CPUTime
#--------------------------------#
#define function #
#--------------------------------#
function partion_grid(nx, x1, x2)
        h = (x2 - x1)/(nx-1)
        y = Array{Float64}(undef, nx)
        for j = 1:nx
            y[j] = (j-1) * h
        end
    return y
end

function BCL(k,s,r,t,opt)
    if opt == "put"
        return K * exp(r * t)
    elseif opt == "call"
        return 0
    end
end
function BCR(k,s,r,t,opt)
    if opt == "put"
        return 0
    elseif opt == "call"
        return s-k * exp(r * t)
    end
end
function IC(k,s,opt)
    if opt == "put"
        return max(k - s, 0)
    elseif opt == "call"
        return max(s - k, 0)
    end
end

#function stepping(m,n,dt,P,S,deltaS,r,i,j,opt)
#    rhs = 0.5*var*S[i]^2*(P[j,i+1]-2*P[j,i]+P[j,i-1])/(deltaS^2) + r*S[i]*(P[j,i+1]-P[j,i-1])/2/deltaS - r*P[j,i]
#    P[j+1,i] = rhs * dt + P[j,i]
#    return P[j+1,i]
#end

#--------------------------------#
#Initialize field #
#--------------------------------#
opt="call"
CPUtic()
K = 100  # strike
S0 = 120  # Initial Stock Price
r = 0.1  # Risk Free Interest Rate
sig = 0.16  # Volatility
Smax = 150  # Largest value of the underlying asset
Tmax = 0.25   # Time to Expiration in years
var = sig * sig
N = 2000  # Number of time steps
M = 200  # Number of stock steps
Time = partion_grid(N, 0, Tmax)
Time = Time[end:-1:1]
Stock = partion_grid(M, 0, Smax)
deltaS = Stock[2] - Stock[1]
deltaT = -Time[2] + Time[1] # delta T should be poitive if the BS is inversed for tau=T-t
RKstep = [0.25, 0.3333, 0.5, 1]


S = Stock
P = zeros(N, M)
for ii = 1:M
    P[1, ii] = IC(K,ii * deltaS, opt)
end
for jj = 1:N
    P[jj, M] = BCR(K,Smax,r,jj * deltaT, opt)
    P[jj, 1] = BCL(K,Smax,r,jj * deltaT, opt)
end

for j = 1:N-1 for k = 1:4 for i = 2:M-1
    coeff = RKstep[k]
    dt = coeff*deltaT
    rhs = 0.5*var*S[i]^2*(P[j,i+1]-2*P[j,i]+P[j,i-1])/(deltaS^2) + r*S[i]*(P[j,i+1]-P[j,i-1])/2/deltaS - r*P[j,i]
    P[j+1,i] = rhs * dt + P[j,i]
end end end

plot(P[N, 1:M])
CPUtoc()
