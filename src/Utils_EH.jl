using Polynomials

# PolynomialMatrices @Neveritt (PolynomialMatrices Package not up to date anymore)
include("PolyMatrixDet.jl")

#############################################################################################################
# Ellipsoidal Harmonic Weights and Normalization
#############################################################################################################

function getGammaNormalization(m,n,Lame,ellipsoidalTDesign)
    x = map((x->getEllipsoidalCoords(x,a)),ellipsoidalTDesign)
    return 4*pi/length(x) * sum([getSmn(x[i],m,n,Lame)^2 for i in eachindex(x)])
end

function getWeightsA(m,n,Lame,ellipsoidalTDesign,F,a)
    x = map((x->getEllipsoidalCoords(x,a)),ellipsoidalTDesign)
    Integral = real(4*pi/length(x) * sum([F[i]*getSmn(x[i],m,n,Lame) for i in eachindex(x)]))
    Normalization = real(1/(getGammaNormalization(m,n,Lame,ellipsoidalTDesign)*Lame[m,n+1](a[1])))
    return Integral*Normalization
end

#############################################################################################################
# Build Up Inner Ellipsoidal Harmonics and Surface Ellipsoidal Harmonics
#############################################################################################################

# Surface EH

function getSmn(coord,m,n,Lame)
    μ,ν,signs = coord[2:4]
    if iseven(n)
        if m <= n/2+1 # Kn: Vorzeichen schon korrekt
            return Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= n+1 #Ln: sign(y) == sign(sqrt(h2[3]-μ^2)
            return signs[2]*Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= n+1+n/2 #Mn: sign(z) == sign(sqrt(h2[2]-ν^2)
            return Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        else #Nn: sign(y) == sign(sqrt(h2[3]-μ^2) && sign(z) == sign(sqrt(h2[2]-ν^2)
            return signs[2]*Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        end
    else
        if m <= (n+1)/2 # Kn: Vorzeichen schon korrekt
            return Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= n+1 #Ln: sign(y) == sign(sqrt(h2[3]-μ^2)
            return signs[2]*Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= (n+1)/2+n+1 #Mn: sign(z) == sign(sqrt(h2[2]-ν^2)
            return Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        else #Nn: sign(y) == sign(sqrt(h2[3]-μ^2) && sign(z) == sign(sqrt(h2[2]-ν^2)
            return signs[2]*Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        end    
    end 
end

# inner EH

function getEmn(coord,m,n,Lame)
    ρ,μ,ν,signs = coord
    if iseven(n)
        if m <= n/2+1 # Kn: Vorzeichen schon korrekt
            return Lame[m,n+1](ρ)*Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= n+1 #Ln: sign(y) == sign(sqrt(h2[3]-μ^2)
            return Lame[m,n+1](ρ)*signs[2]*Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= n+1+n/2 #Mn: sign(z) == sign(sqrt(h2[2]-ν^2)
            return Lame[m,n+1](ρ)*Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        else #Nn: sign(y) == sign(sqrt(h2[3]-μ^2) && sign(z) == sign(sqrt(h2[2]-ν^2)
            return Lame[m,n+1](ρ)*signs[2]*Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        end
    else
        if m <= (n+1)/2 # Kn: Vorzeichen schon korrekt
            return Lame[m,n+1](ρ)*Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= n+1 #Ln: sign(y) == sign(sqrt(h2[3]-μ^2)
            return Lame[m,n+1](ρ)*signs[2]*Lame[m,n+1](μ)*Lame[m,n+1](ν)
        elseif m <= (n+1)/2+n+1 #Mn: sign(z) == sign(sqrt(h2[2]-ν^2)
            return Lame[m,n+1](ρ)*Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        else #Nn: sign(y) == sign(sqrt(h2[3]-μ^2) && sign(z) == sign(sqrt(h2[2]-ν^2)
            return Lame[m,n+1](ρ)*signs[2]*Lame[m,n+1](μ)*signs[3]*Lame[m,n+1](ν)
        end    
    end 
end

# Lame Function Build Up

function evalLame(x,Lame)
    map(f -> f(x), Lame)
end

function getLameFunctions(N,a)
    hcat([append!(getLameFunctionsN(n,a),[(x -> nothing) for i=1:(2N-2n)]) for n=0:N]...)
end

function getLameFunctionsN(n,a)
    # calculate Hs and r
    h2 = (a[2]^2-a[3]^2,a[1]^2-a[3]^2,a[1]^2-a[2]^2)
    
    if n==0
        return vcat( getLameKn(n,h2)
                    )
    elseif n==1
        return vcat( getLameKn(n,h2),
                     getLameLnMn(n,h2)
                    )   
    else
        return vcat(getLameKn(n,h2),
                    getLameLnMn(n,h2),
                    getLameNn(n,h2)
                    )
    end
end

#############################################################################################################
# Low-Level Lame Functions
#############################################################################################################

function getLameKn(n,h2)
    r = iseven(n) ? Int(n/2) : Int((n-1)/2);

    K=Array{Polynomial{Float64}}(undef,r+1,r+1)

    K.=Polynomial(0,:p)
    for i=1:r+1
        K[i,i] = Polynomial([(h2[3]+h2[2])*(n-2i+2)^2,-(h2[3]+h2[2])], :p)
    end
    for j=1:r
        K[j,j+1] = Polynomial([2j*(2n-2j+1)],:p)
        K[j+1,j] = Polynomial([-(h2[3]*h2[2])*(n-2j+1)*(n-2j+2)],:p)
    end

    K=PolyMatrix(K)

    p=roots(det(K))

    # Faktoren α
    α=ones(Complex{Float64},r+1,r+1)

    for i=1:r+1, k=1:r
        if k==1
            α[i,k+1] = ( (h2[3]+h2[2])*(p[i]-(n-2(k-1))^2)*α[i,k] ) / ((2(k-1)+2)*(2n-1-2(k-1)))
        else
            α[i,k+1] = ( (h2[3]+h2[2])*(p[i]-(n-2(k-1))^2)*α[i,k] + h2[3]*h2[2]*(n-2(k-1)+2)*(n-2(k-1)+1)*α[i,k-1] ) / ((2(k-1)+2)*(2n-1-2(k-1)))
        end
    end

    Kn = Array{Polynomial{Complex{Float64}}}(undef,r+1)
    tmp=zeros(Complex{Float64},r+1,n+1)
    for i=1:r+1
        for k=0:r
            tmp[i,1+n-2k]=α[i,k+1]
        end    
        Kn[i] = Polynomial(tmp[i,:])
    end
    v = Function[]
    for i=1:r+1
        push!(v,(x -> Kn[i](x)))
    end
    return v
end

#-------------------------------------------------------------------------------

function getLameLn(n,h2,r)
    L=Array{Polynomial{Float64}}(undef,r,r)

    L.=Polynomial(0,:p)
    for i=1:r
        L[i,i] = Polynomial([(h2[3]+h2[2])*(n-2i+1)^2+(2n-4i+3)*h2[2],-(h2[3]+h2[2])], :p)
    end
    for j=1:r-1
        L[j,j+1] = Polynomial([2j*(2n-2j+1)],:p)
        L[j+1,j] = Polynomial([-(h2[3]*h2[2])*(n-2j+1)*(n-2j)],:p)
    end
    
    L=PolyMatrix(L)

    p=roots(det(L))

    # Faktoren  β
    β=ones(Complex{Float64},r,r)

    for i=1:r, k=1:r-1
        if k==1
            β[i,2] = ((h2[3]+h2[2])*(p[i]-(n-1)^2)-(2n-1)*h2[2]) / (2*(2n-1))
        else
            β[i,k+1] =  ( ((h2[3]+h2[2])*(p[i]-(n-2(k-1)-1)^2)-(2n-4(k-1)-1)*h2[2])*β[i,k] 
                            + h2[3]*h2[2]*(n-2(k-1)+1)*(n-2(k-1))*β[i,k-1] ) / 
                        ((2((k-1)+1))*(2n-1-2(k-1)))
        end
    end

    Ln = Array{Polynomial{Complex{Float64}}}(undef,r) #eig Ln-1
    tmp=zeros(Complex{Float64},r,n)
    for i=1:r
        for k=0:r-1
            tmp[i,1+n-1-2k]=β[i,k+1]
        end    
        Ln[i] = Polynomial(tmp[i,:])
    end

    return Ln
end


function getLameMn(n,h2,r)
    R=Array{Polynomial{Float64}}(undef,r,r)

    R.=Polynomial(0,:p)
    for i=1:r
        R[i,i] = Polynomial([(h2[2]+h2[3])*(n-2i+1)^2+(2n-4i+3)*h2[3],-(h2[2]+h2[3])], :p)
    end
    for j=1:r-1
        R[j,j+1] = Polynomial([2j*(2n-2j+1)],:p)
        R[j+1,j] = Polynomial([-(h2[2]*h2[3])*(n-2j+1)*(n-2j)],:p)
    end

    R=PolyMatrix(R)

    p=roots(det(R))

    # Faktoren  β
    β=ones(Complex{Float64},r,r)

    for i=1:r, k=1:r-1
        if k==1
            β[i,2] = ((h2[2]+h2[3])*(p[i]-(n-1)^2)-(2n-1)*h2[3]) / (2*(2n-1))
        else
            β[i,k+1] = ( ((h2[2]+h2[3])*(p[i]-(n-2(k-1)-1)^2)-(2n-4(k-1)-1)*h2[3])*β[i,k] + h2[2]*h2[3]*(n-2(k-1)+1)*(n-2(k-1))*β[i,k-1] ) / ((2((k-1)+1))*(2n-1-2(k-1)))
        end
    end

    Rn = Array{Polynomial{Complex{Float64}}}(undef,r) #eig Rn-1
    tmp=zeros(Complex{Float64},r,n)
    for i=1:r
        for k=0:r-1
            tmp[i,1+n-1-2k]=β[i,k+1]
        end    
        Rn[i] = Polynomial(tmp[i,:])
    end

    return Rn
end

function getLameLnMn(n,h2)
    r = iseven(n) ? Int(n/2) : Int((n+1)/2);
    
    if r == 0
        @error "No Ln/Mn for n = 0"
    else
        Ln=getLameLn(n,h2,r)
        Mn=getLameMn(n,h2,r)
        return vcat( [(x ->  sqrt(abs(x^2-h2[3]))*Ln[i](x)) for i=1:r],
                     [(x ->  sqrt(abs(x^2-h2[2]))*Mn[i](x)) for i=1:r]
                )
    end
end


#-------------------------------------------------------------------------------

function getLameSn(n,h2,r)
    S=Array{Polynomial{Float64}}(undef,r,r)

    S.=Polynomial(0,:p)
    for i=1:r
        S[i,i] = Polynomial([(h2[3]+h2[2])*(n-2i+1)^2,-(h2[3]+h2[2])], :p)
    end
    for j=1:r-1
        S[j,j+1] = Polynomial([2j*(2n-2j+1)],:p)
        S[j+1,j] = Polynomial([-(h2[3]*h2[2])*(n-2j-1)*(n-2j)],:p)
    end

    S=PolyMatrix(S)

    p=roots(det(S))

    # Faktoren  β
    γ=ones(Complex{Float64},r,r)

    for i=1:r, k=1:r-1
        if k==1
            γ[i,2] = ( (h2[3]+h2[2])*(p[i]-(n-1)^2)*γ[i,k] ) / (2*(2n-1))
        else
            γ[i,k+1] = ( (h2[3]+h2[2])*(p[i]-(n-2(k-1)-1)^2)*γ[i,k] + h2[3]*h2[2]*(n-2(k-1)-1)*(n-2(k-1))*γ[i,k-1] ) / ((2((k-1)+1))*(2n-1-2(k-1)))
        end
    end

    Sn = Array{Polynomial{Complex{Float64}}}(undef,r) #eig Sn-2
    tmp=zeros(Complex{Float64},r,n-1)
    for i=1:r
        for k=0:r-1
            tmp[i,1+n-2-2k]=γ[i,k+1]
        end    
        Sn[i] = Polynomial(tmp[i,:])
    end

    return Sn
end

function getLameNn(n,h2)
    r = iseven(n) ? Int(n/2) : Int((n-1)/2);
    if r == 0
        @error "No Nn for n ∈ {0,1}"
    else
        Sn = getLameSn(n,h2,r)
        return [(x -> ( sqrt(abs(x^2-h2[3])) * sqrt(abs(x^2-h2[2])) ) * Sn[i](x)) for i=1:r]
    end
end

#############################################################################################################
# Coordinate Transformations and Helper Functions
#############################################################################################################

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
