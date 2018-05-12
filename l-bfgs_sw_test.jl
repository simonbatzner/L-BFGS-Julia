function lbfgs!(F, x0,maxIt,m,τgrad=1e-5)
#INPUT
# F: function to be optimized
# x0: initial guess
# maxIt: maximum Iteration
# m: last m input differences and gradient differences are stored
# τgrad: tolerance for norm of the slope


#OUTPUT
#x1: optimized variable
#f1: function value at x1
#k iteration number

    k=0
    n=length(x0)
    Sm=zeros(n,m) #S_k=x_k+1-x_k
    Ym=zeros(n,m) #Y_k=g_k+1-g_k
    f0,g0=F(x0)
    #use the simplest line search to find step size
    α, f1, g1=strongwolfe(F,-g0,x0,f0,g0)
    x1 = x0 - α.*g0
    k=1

    while true
        if k>maxIt
            break
        end
        gnorm=norm(g0)
        if gnorm < τgrad
            break
        end
        s0=x1-x0
        y0=g1-g0
        #println("y0=$y0")
        H0=s0'*y0/(y0'*y0) #hessian diagonal satisfying secant condition

        #update Sm and Ym
        if k<=m
            Sm[:,k]=s0
            Ym[:,k]=y0
            p=-approxInvHess(g1,Sm[:,1:k],Ym[:,1:k],H0)
        # only keep m entries in Sm and Ym so purge the old ones
        elseif (k>m)
            Sm[:,1:(m-1)]=Sm[:,2:m]
            Ym[:,1:(m-1)]=Sm[:,2:m]
            Sm[:,m]=s0
            Ym[:,m]=y0
            p=-approxInvHess(g1,Sm,Ym,H0)
        end
        # new direction=p, find new step size
        α, fs, gs=strongwolfe(F,p,x1,f1,g1)
        #update for next iteration
        x0=x1
        g0=g1
        x1=x1+α.*p
        f1=fs
        g1=gs
        k=k+1
        println("It=$k,x=$x1")
    end
    k=k-1
    return x1, f1, k
end

function strongwolfe(F,d,x0,fx0,gx0,maxIt=5)
# strong wolfe
    α_m=20
    α_p=0
    c1=1e-4
    c2=0.9
    α_x=1
    gx0=copy(gx0'*d)
    fxp=copy(fx0)
    gxp=copy(gx0)
    i=1
    α_s=0
    fs=copy(fx0)
    gs=copy(gx0)
    while true
        xx=x0+α_x*d
        fxx,gxx=F(xx)
        fs=copy(fxx)
        gs=copy(gxx)
        gxx=copy(gxx'*d)

        if (fxx>(fx0+c1*α_x*gx0)[1]) || (i>1) & (fxx>=fxp)
            α_s,fs,gs=zoom(F,x0,d,α_p,α_x,fx0,gx0)
            return α_s,fs,gs
        end
        if abs(gxx)<=-c2*(gx0)
            α_s=copy(α_x)
            return α_s,fs,gs
        end
        if gxx>=0
        #if abs.(gxx)[1]>=0 && abs.(gxx)[2]>=0
            α_s,fs,gs=zoom(F,x0,d,α_x,α_p,fx0,gx0)
            return α_s,fs,gs
        end
        α_p=copy(α_x)
        fxp=copy(fxx)
        gxp=copy(gxx)

        if i>maxIt
            α_s=α_x
            return α_s,fs,gs

        end
        r=0.8
        #r=0.8
        α_x=α_x+(α_m-α_x)*r
        i=i+1

    end
    return α_s,fs,gs
end

function zoom(F,x0,d,α_l,α_h,fx0,gx0,maxIt=5)
    c1=1e-4
    c2=0.9
    i=0
    α_s=0
    fs=copy(fx0)
    gs=copy(gx0)
    while true
        α_x=0.5*(α_l+α_h)
        α_s=copy(α_x)
        xx=x0+α_x*d
        fxx,gxx=F(xx)
        fs=copy(fxx)
        gs=copy(gxx)
        gxx=gxx'*d
        xl=x0+α_l*d
        fxl,gxl=F(xl)
        if (fxx>(fx0+c1*α_x*gx0)[1]) || fxx>=fxl
            α_h=copy(α_x)
        else
            if abs(gxx)[1]<=-c2*(gx0)
                α_s=copy(α_x)
                return α_s,fs,gs
            end
            if gxx*(α_h-α_l)[1]>=0
                α_h=copy(α_l)
            end
            α_l=copy(α_x)
        end
        i=i+1
        if i>maxIt
            α_s=copy(α_x)
            return α_s,fs,gs
        end
    end
    return α_s,fs,gs
end

function approxInvHess(g,S,Y,H0)
    #INPUT

    #g: gradient nx1 vector
    #S: nxk matrixs storing S[i]=x[i+1]-x[i]
    #Y: nxk matrixs storing Y[i]=g[i+1]-g[i]
    #H0: initial hessian diagnol scalar

    #OUTPUT
    # p:  the approximate inverse hessian multiplied by the gradient g
    #     which is the new direction
    #notation follows:
    #https://en.wikipedia.org/wiki/Limited-memory_BFGS

    n,k=size(S)
    rho=zeros(k)
    for i=1:k
        rho[i]=1/(Y[:,i]'*S[:,i])
    end

    q=zeros(n,k+1)
    r=zeros(n,1)
    α=zeros(k,1)
    β=zeros(k,1)

    q[:,k+1]=g

    for i=k:-1:1
        α[i]=rho[i]*S[:,i]'*q[:,i+1]
        q[:,i]=q[:,i+1]-α[i]*Y[:,i]
    end

    z=H0*q[:,1]


    for i=1:k
        β[i]=rho[i]*Y[:,i]'*z
        z=z+S[:,i]*(α[i]-β[i])
    end

    p=z

    return p
end

x0=[1,3]
include("testfns.jl")

x1, f1, k=lbfgs!(himmelblau,x0,100,2)
