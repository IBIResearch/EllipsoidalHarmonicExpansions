using Polynomials, PolynomialMatrices, LinearAlgebra

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

#-------------------------------------------------------------------------------

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

