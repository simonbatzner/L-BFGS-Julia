function Engval(x::Vector)
    #if n=1000, fmin=1108
    #if n=5000, fmin=5548

    F=0.0
    f=zeros(n-1)
    g=zeros(n)
    for i=2:n
        f=((x[i-1])^2+(x[i])^2)^2-4*(x[i-1])+3

        F +=f
    end
    g[1]=x[1]*((x[1])^2+(x[2])^2)-1
    g[n]=x[n]*((x[n-1])^2+(x[n])^2)
    for i=2:n-1
        g[i]=x[i]*((x[i-1])^2+2*(x[i])^2+(x[i+1])^2)-1
    end
    return F, g
end
