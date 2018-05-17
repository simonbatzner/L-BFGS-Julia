function rosenbrock(x::Vector)
    # valley shaped
    #   (1,1)
    #   x_i ∈ [-2.048, 2.048]
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
    #   x_i ∈ [-6, 6]
    F=(x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
    d=[4*x[1]*(x[1]^2+x[2]-11)+2*(x[1]+x[2]^2-7),2*(x[1]^2+x[2]-11)+4*x[2]*(x[1]+x[2]^2-7)]
    return F, d
end

function booth(x::Vector)
    # plate shaped
    #   (1,3)
    #   x_i ∈ [-10, 10]
    F=(x[1]+2*x[2]-7)^2+(2*x[1]+x[2]-5)^2
    d=[8*x[2]+10*x[1]-34,10*x[2]+8*x[1]-38]
    return F, d
end

