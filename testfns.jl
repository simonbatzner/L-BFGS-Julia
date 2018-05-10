function rosenbrock(x::Vector)
    # valley shaped
    #   (1,1)
    F=(1-x[1])^2+100*(x[2]-x[1]^2)^2
    d=[-2*(1-x[1])-400*(x[2]-x[1]^2)*x[1],200 * (x[2] - x[1]^2)]
   return F, d
end

function himmelblau(x::Vector)
    # many local minima
    #   (3,2)
    #   (-2.805118,3.131312)
    #   (-3.779310, -3.283186)
    #   (3.584428,-1.848126)
    F=(x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
    d=[4*x[1]*(x[1]^2+x[2]-11)+2*(x[1]+x[2]^2-7),2*(x[1]^2+x[2]-11)+4*x[2]*(x[1]+x[2]^2-7)]
    return F, d
end

function booth(x::Vector)
    # plate shaped
    #   (1,3)
    F=(x[1]+2*x[2]-7)^2+(2*x[1]+x[2]-5)^2
    d=[8*x[2]+10*x[1]-34,10*x[2]+8*x[1]-38]
    return F, d
end

function bohachevsky1(x::Vector)
    # bowl shaped
    #   (0,0)
    F= x[1]^2+2*x[2]^2-0.3*cos(3*pi*x[1])-0.4*cos(4*pi*x[2])+0.7
    d=[2*x[1]+0.3*3*pi*sin(3*pi*x[1]),4*x[2]+0.4*4*pi*sin(4*pi*x[2])]
    return F, d
end

function easom(x::Vector)
    # steep ridge and drop
    #   (pi,pi)
    F=-cos(x[1])*cos(x[2])*exp(-(x[1]-pi)^2-(x[2]-pi)^2)
    d=[exp(-(x[1]-pi)^2-(x[2]-pi)^2)*cos(x[2])*(sin(x[1])+2*(x[1]-pi)*cos(x[1])),2*(x[2]-pi)*exp(-(x[2]-pi)^2-(x[1]-pi)^2)*cos(x[2])*cos(x[1])+exp(-(x[2]-pi)^2-(x[1]-pi)^2)*cos(x[1])*sin(x[2])]
    return F, d
end
