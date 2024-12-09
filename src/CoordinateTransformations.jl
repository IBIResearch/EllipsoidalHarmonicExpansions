function getEllipsoidalTDesign(L,a)
    if 2L <= 2
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2,9,1u"m"))
    elseif 2L == 4
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,21,1u"m"))
    elseif 2L == 6
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,35,1u"m"))
    elseif 2L == 8
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,45,1u"m"))
    elseif 2L == 10
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,62,1u"m"))
    elseif 2L == 12
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,86,1u"m"))
    elseif 2L == 14
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,114,1u"m"))
    elseif 2L == 16
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(2L,146,1u"m"))
    elseif 2L > 16 && 2L <= 18 
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(17,156,1u"m"))
    else
        sphericalTDesign = map(x -> uconvert.(Unitful.m, x), MPIFiles.loadTDesign(12,86,1u"m"))
    end
    return getCoordinatesOnEllipsoid.(sphericalTDesign,[a])
end 

function getCoordinatesOnSphere(ϑ,φ)
    return (sin(ϑ)*cos(φ),sin(ϑ)*sin(φ),cos(ϑ))
end

function getCoordinatesOnEllipsoid(y,a)
    ρ=a[1]
    h2 = (a[2]^2-a[3]^2,a[1]^2-a[3]^2,a[1]^2-a[2]^2)
    return (y[3]*ρ,y[1]*sqrt(ρ^2-h2[3]),y[2]*sqrt(ρ^2-h2[2]))
end

function getEllipsoidalCoords(x,a)
    x = (y->ustrip.(y)).(x)
    c2 = a[1]^2+a[2]^2+a[3]^2-x[1]^2-x[2]^2-x[3]^2
    c1 = a[1]^2*a[2]^2+a[1]^2*a[3]^2+a[2]^2*a[3]^2-(a[2]^2+a[3]^2)*x[1]^2-(a[1]^2+a[3]^2)*x[2]^2-(a[1]^2+a[2]^2)*x[3]^2
    c0 = a[1]^2*a[2]^2*a[3]^2-a[2]^2*a[3]^2*x[1]^2-a[1]^2*a[3]^2*x[2]^2-a[1]^2*a[2]^2*x[3]^2
    p = (c2^2-3*c1)/9
    q = (9*c1*c2-27*c0-2*c2^3)/54
    abs(q) < abs(p^(3/2)) ? ω = acos(q/p^(3/2)) : ω = 0 
    #p,q,ω
    s = (2*sqrt(p)*cos(ω/3)-c2/3, 2*sqrt(p)*cos((ω-2*pi)/3)-c2/3, 2*sqrt(p)*cos((ω-4*pi)/3)-c2/3)
    ρ = sqrt(round(a[1]^2+s[1],digits=10))
    μ = sqrt(round(a[1]^2+s[2],digits=10))
    ν = sign(x[1])*sqrt(round(a[1]^2+s[3],digits=10))
    signs = sign.(x)
    return (ρ,μ,ν,signs)
end

function getCartesianCoords(y,a)
    (ρ,μ,ν,signs) = y
    (hy,hx,hz) = (sqrt(round(a[1]^2-a[3]^2,digits=16)),sqrt(round(a[2]^2-a[3]^2,digits=16)),sqrt(round(a[1]^2-a[2]^2,digits=16)))
    x1=signs[1]*abs((ρ*μ*ν)/(hy*hz))
    x2=(sqrt(round(ρ^2-hz^2,digits=16))*sqrt(round(μ^2-hz^2,digits=16)) * signs[2]*sqrt(round(hz^2-ν^2,digits=16)))/(hx*hz)
    x3=(sqrt(round(ρ^2-hy^2,digits=16))*signs[3]*sqrt(round(hy^2-μ^2,digits=16))*sqrt(round(hy^2-ν^2,digits=16)))/(hy*hx)
    return (x1,x2,x3)
end

function getellipsepoints(cx, cy, rx, ry, θ)
    t = range(0, 2*pi, length=100)
    ellipse_x_r = @. rx * cos(t)
    ellipse_y_r = @. ry * sin(t)
    R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
    r_ellipse = [ellipse_x_r ellipse_y_r] * R
    x = @. cx + r_ellipse[:,1]
    y = @. cy + r_ellipse[:,2]
    (x,y)
end

function inside_ellipsoid(x,y,z;a=(1.0,1.0,1.0))
	if x^2/a[1]^2 + y^2/a[2]^2 + z^2/a[3]^2 <= 1
        return true
    else
        return false
    end
end
