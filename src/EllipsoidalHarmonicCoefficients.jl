using MPIFiles, Polynomials, PolynomialMatrices, LinearAlgebra, MPISphericalHarmonics, Unitful, MPIMagneticFields

include("LameFunctions.jl")
include("CoordinateTransformations.jl")

#############################################################################################################
# Coefficients
#############################################################################################################

mutable struct EllipsoidalHarmonicCoefficients{T<:Real}
    c::Matrix{T} # coefficients
    radius::Tuple{T,T,T} # half-axis of ellipsoid
end

mutable struct EllipsoidalMagneticFieldCoefficients 
    coeffs::Array{EllipsoidalHarmonicCoefficients,2} # coefficients
    radius::Tuple{Float64,Float64,Float64} # radii of measured ellipsoid
    center::Array{Float64,2} # maybe we need it later
    ffp::Union{Array{Float64,2},Nothing} # field-free-point (if available, depends on expansion point)
end

EllipsoidalMagneticFieldCoefficients(coeffs::Array{EllipsoidalHarmonicCoefficients,2},radius::Tuple{Float64,Float64,Float64}) =
    EllipsoidalMagneticFieldCoefficients(coeffs, radius, [0.0,0.0,0.0])

EllipsoidalMagneticFieldCoefficients(
        coeffs::Array{EllipsoidalHarmonicCoefficients,2},
        radius::Tuple{Float64,Float64,Float64},
        center::VecOrMat{Float64}) = 
        EllipsoidalMagneticFieldCoefficients(coeffs, radius, center, nothing)

EllipsoidalMagneticFieldCoefficients(
        coeffs::Array{EllipsoidalHarmonicCoefficients,2},
        radius::Tuple{Float64,Float64,Float64},
        center::Vector{Float64},
        ffp::Union{Matrix{Float64},Nothing}) = 
        EllipsoidalMagneticFieldCoefficients(
                coeffs,
                radius,
                hcat([center for p = 1:size(coeffs, 2)]...),
                ffp)

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

# get ellipsoidal harmonic coefficients for a given field, reference ellipsoid (half-axis) and degree L

function magneticField(
    L::Int,
    field_data::Union{AbstractArray{T,2},AbstractArray{T,3}},
    a::Tuple{Float64,Float64,Float64}) where {T<:Real}

    # get ellipsoidal tDesign positions [m] and removing the unit
    ellipsoidalTDesign = getEllipsoidalTDesign(L,a)
    Lame=getLameFunctions(L,a)
    A = zeros(L+1,2L+1,3)
    for n=0:L
        for m=1:2*n+1
            A[n+1,m,1] = Real(getWeightsA(m,n,Lame,ellipsoidalTDesign,field_data[1,:],a))
            A[n+1,m,2] = Real(getWeightsA(m,n,Lame,ellipsoidalTDesign,field_data[2,:],a))
            A[n+1,m,3] = Real(getWeightsA(m,n,Lame,ellipsoidalTDesign,field_data[3,:],a))
        end
    end
    # until now only one Patch
    c=Array{EllipsoidalHarmonicCoefficients,2}(undef,3,1)    
    c[1,1]=EllipsoidalHarmonicCoefficients(A[:,:,1],a)
    c[2,1]=EllipsoidalHarmonicCoefficients(A[:,:,2],a)
    c[3,1]=EllipsoidalHarmonicCoefficients(A[:,:,3],a)
    
    return c
end

## return EllipsoidalMagneticFieldCoefficients for given field, reference ellipsoid (half-axis) and degree L

function EllipsoidalMagneticFieldCoefficients(
    L::Int,
    field_data::Union{AbstractArray{T,2},AbstractArray{T,3}},
    a::Tuple{Float64,Float64,Float64}) where {T<:Real}
    
    return EllipsoidalMagneticFieldCoefficients(magneticField(L,field_data,a),a)
end

#############################################################################################################
# Fields
#############################################################################################################

Base.@kwdef mutable struct EllipsoidalHarmonicsDefinedField <: AbstractMagneticField
    coeffs::EllipsoidalMagneticFieldCoefficients
    Lame::Matrix{Function}
    patch::Integer = 1
end

function EllipsoidalHarmonicsDefinedField(coeffs::EllipsoidalMagneticFieldCoefficients)
    return EllipsoidalHarmonicsDefinedField(coeffs=coeffs,Lame=getLameFunctions(Int(size(coeffs.coeffs[1].c,1)-1),coeffs.radius))
end

MPIMagneticFields.FieldStyle(::EllipsoidalHarmonicsDefinedField) = OtherField()
MPIMagneticFields.FieldDefinitionStyle(::EllipsoidalHarmonicsDefinedField) =
    SphericalHarmonicsDataBasedFieldDefinition()
MPIMagneticFields.FieldTimeDependencyStyle(::EllipsoidalHarmonicsDefinedField) =
    TimeConstant()

function evaluateField(
    field::EllipsoidalHarmonicsDefinedField,
    grid::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}})

    a = field.coeffs.radius

    VF=Array{Tuple{Float64,Float64,Float64}}(undef,length.(grid)...)

    for x=1:length(grid[1]),y=1:length(grid[2]),z=1:length(grid[3])
        coord=getEllipsoidalCoords((grid[1][x], grid[2][y], grid[3][z]),a)
        #@info "$x,$y,$z inside Ellipsoid"
        fx,fy,fz = 0.0,0.0,0.0
        for n=0:size(field.Lame,2)-1
            for m=1:2*n+1
                E=Real(getEmn(coord,m,n,field.Lame))
                fx+=coeffs.coeffs[1,1].c[n+1,m]*E
                fy+=coeffs.coeffs[2,1].c[n+1,m]*E
                fz+=coeffs.coeffs[3,1].c[n+1,m]*E
            end
        end
        VF[x,y,z]=(fx,fy,fz)
    end
    return reshape(VF,length.(grid)...)
end

function evaluateField(
    field::EllipsoidalHarmonicsDefinedField,
    r::Vector{Float64})

    coord=getEllipsoidalCoords(r,field.coeffs.radius)
    
    fx,fy,fz = 0.0,0.0,0.0
    for n=0:size(field.Lame,2)-1
        for m=1:2*n+1
            E=Real(getEmn(coord,m,n,field.Lame))
            fx+=field.coeffs.coeffs[1,1].c[n+1,m]*E
            fy+=field.coeffs.coeffs[2,1].c[n+1,m]*E
            fz+=field.coeffs.coeffs[3,1].c[n+1,m]*E
        end
    end
    return [fx,fy,fz]
end

# get field values at position r as [field_x field_y field_z]
MPIMagneticFields.value_(field::EllipsoidalHarmonicsDefinedField, r) =
    evaluateField(field,r)
