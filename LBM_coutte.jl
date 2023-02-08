using Statistics
using Plots
using BenchmarkTools

# Main 

function lbm(t)
        rho = sum(fprop,dims = 3)
        
        # Momentum Components
        
        global u = sum(fprop[:,:,[1 5 8]], dims = 4) - sum(fprop[:,:,[3 6 7]], dims = 4)
        global v = sum(fprop[:,:,[2 5 6]], dims = 4) - sum(fprop[:,:,[4 7 8]], dims = 4)
         
        for k=1:NPOP
            # Compute equilibrium distribution (linear equilibrium with incompressible model)
            #feq(:,:,k)=w(k)*(rho + 3*(u*cx(k)+v*cy(k)));
            global feq[:,:,k]=w[k].*rho.*(1 .+ 3*(u*cx[k]+v*cy[k]) + 9*((u*cx[k]+v*cy[k]).^2)/2 - 3*(u.^2+v.^2)/2  )
        end
        
        #Collision step
        fnew = (1-omega)*fprop + omega*feq
        
        for k=1:NPOP
            for j=1:NY
                for i=1:NX
                    # Streaming step (Periodic streaming of whole domain)
                    newx=1+mod(i-1+cx[k]+NX,NX)
                    newy=1+mod(j-1+cy[k]+NY,NY)
                    global fprop[newx,newy,k]=fnew[i,j,k]
                end
            end
        end
        
        # Boundary condition (bounce-back)
        # Top wall (moving with tangential velocity u_max)
        fprop[:,NY,4]=fnew[:,NY,2]
        fprop[:,NY,7]=fnew[:,NY,5].-(1/6)*u_max
        fprop[:,NY,8]=fnew[:,NY,6].+(1/6)*u_max
        
        # Bottom wall (rest)
        fprop[:,1,2]=fnew[:,1,4];
        fprop[:,1,5]=fnew[:,1,7];
        fprop[:,1,6]=fnew[:,1,8];
        
        for i=1:NX
            error[i]=(sqrt(sum((u[i,:]-u_analy).^2)))./sqrt(sum(u_analy.^2))
        end
        err[1,t] = 1/NX*sum(error) 
end

# Variables

scale = 1
NX = 10*scale;
NY = 10*scale;
NSTEPS = 10000
tau = 0.9
omega = 1/tau
u_max = 1/scale
nu = (2*tau-1)/6
Re = NY*u_max/nu

NPOP = 9                                       # number of velocities
w  = [1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9] # weights
cx = [1 0 -1  0 1 -1 -1  1 0]                  # velocities, x components
cy = [0 1  0 -1 1  1 -1 -1 0]

x = 1-0.5:NX-0.5
y = 1-0.5:NY-0.5

u_analy=u_max/NY*y

feq = zeros(NX,NY,NPOP)
fnew = zeros(NX,NY,NPOP)
fprop = zeros(NX,NY,NPOP)

u = zeros(NX,NY)
v = zeros(NX,NY)

for k = 1:NPOP
    feq[:,:,k] .= w[k]
    fnew[:,:,k] .= w[k]
    fprop[:,:,k] .= w[k]
end

# convergence parameters
tol = 1e-6      # tolerance to steady state convergence
teval = 100      # time step to evaluate convergence
u_old = zeros(NX,NY)

error = zeros(NX,1)
err = zeros(1,NSTEPS)

@time for t = 1:NSTEPS
    u_old = u
    lbm(t) 
    conv = abs(mean(u[:])/mean(u_old[:])-1)
    global step = t
    
    # Check Convergence
            
    if conv <= tol
        break
    end
        
end

err = err[1,(1:step)]
time = zeros(1,step)
for i = 1:step
    time[1,i] = i;
end

println("Converges in " ,string(Int(time[end])), " time steps")

# Error vs Time step Plot
plot(time[1,:],err, title="Error vs Time Steps Graph", legend = false)

Plots.scatter(u[1,:], y/NY, title="Velocity vs Normalized Height y/NY ", label="LBM")
plot!(u_analy,y/NY, label="Analytical")
xlabel!("Normalized Velocity")
ylabel!("y/NY")

u


