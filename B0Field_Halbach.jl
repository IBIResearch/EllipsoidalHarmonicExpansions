using Pkg
Pkg.activate(".")
Pkg.instantiate()

using LinearAlgebra, Interpolations, DataFrames, CSV, Unitful, Plots

include("src/EllipsoidalHarmonicCoefficients.jl")

###############################################################################
# pre definitions #
###############################################################################

gridSizeX = 61
gridSizeY = 27
gridSizeZ = 27

shape = [gridSizeX,gridSizeY,gridSizeZ]
FoV 		= [0.6,0.26,0.26]
center 		= [0,0,0]

grid = [range(center[i]-FoV[i]/2.0,center[i]-FoV[i]/2.0 + FoV[i],shape[i]) for i=1:3]

###############################################################################
# load Field Data and interpolate to grid
###############################################################################

# load simulated field
H_ = Array(CSV.read("data/b0FieldHalbach.csv",DataFrame;header=false))

# set up for interpolation
itpH_X = Interpolations.scale(Interpolations.interpolate(reshape(H_[:,1],shape...), BSpline(Cubic(Line(OnCell())))),grid...)
itpH_Y = Interpolations.scale(Interpolations.interpolate(reshape(H_[:,2],shape...), BSpline(Cubic(Line(OnCell())))),grid...)
itpH_Z = Interpolations.scale(Interpolations.interpolate(reshape(H_[:,3],shape...), BSpline(Cubic(Line(OnCell())))),grid...)

###############################################################################
# Build up Reference Ellipsoid and load T-Design
###############################################################################

# define polynomial degree for t-design
L=7; N=2*L

# define reference ellipsoid
a = (0.26,0.125,0.1249)

# build t-design on reference ellipsoid and extract nodes
ellipsoidalDesign = getEllipsoidalTDesign(L,a)
points = (y->ustrip.(y)).(ellipsoidalDesign)

# interpolate field data to the t-design
points_X = (x->itpH_X(x...)).(points)
points_Y = (x->itpH_Y(x...)).(points)
points_Z = (x->itpH_Z(x...)).(points)
pointValues = hcat(points_X,points_Y,points_Z)'

###############################################################################
# Calculate Field using Ellipsoidal Harmonics
###############################################################################

# build ellipsoidal field coefficients
ellipsoidal_coeffs = EllipsoidalMagneticFieldCoefficients(L,pointValues,a)

# build ellipsoidal field
efield = EllipsoidalHarmonicsDefinedField(ellipsoidal_coeffs)

# evaluate ellipsoidal field at grid points
VF = [efield[x,y,z] for x in grid[1], y in grid[2], z in grid[3]]

###############################################################################
# Plotting
###############################################################################

function arrow0!(x, y, u, v; as=0.07, lw=1, lc=:black, la=1)
    nuv = sqrt(u^2 + v^2)
    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
    v5 = v4 - 2*(v4'*v2)*v2
    v4, v5 = as*nuv*v4, as*nuv*v5
    plot!([x,x+u], [y,y+v], lw=lw, lc=lc, la=la)
    plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lw=lw, lc=lc, la=la)
    plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lw=lw, lc=lc, la=la)
end

gr(legend=false, dpi=600)
sc_ar = 0.12
c=(0.045,0.06)
c_diff = (-20.0,20.0)

H=[(itpH_X(x,y,z),itpH_Y(x,y,z),itpH_Z(x,y,z)) for x in grid[1], y in grid[2], z in grid[3]]
normH = norm.(H)
normH_EH = norm.(VF)

# XY 

Nxy = hcat([[grid[1][i],grid[2][j]] for i=1:gridSizeX for j=1:gridSizeY]...)
H_X = [itpH_X(x,y,0) for x in grid[1], y in grid[2]]
H_Y = [itpH_Y(x,y,0) for x in grid[1], y in grid[2]]

Plots.heatmap(grid[1],grid[2],normH[:,:,14]',aspect_ratio=:equal,clims=c,cmap=:bam,xlabel = "X",ylabel = "Y", title="SimField");
arrow0!.(Nxy[1,:],Nxy[2,:],sc_ar .* vec(permutedims(H_X,[2,1])),sc_ar .* vec(permutedims(H_Y,[2,1])); as=0.3, lw=0.5,lc=:black, la=1);
P_XY = Plots.current();

Plots.heatmap(grid[1],grid[2],normH_EH[:,:,14]',aspect_ratio=:equal,clims=c,cmap=:bam,xlabel = "X",ylabel = "Y", title="EH_Field");
arrow0!.(Nxy[1,:],Nxy[2,:],sc_ar .* vec(permutedims((h->h[1]).(VF)[:,:,14],[2,1])),sc_ar .* vec(permutedims((h->h[2]).(VF)[:,:,14],[2,1])); as=0.3, lw=0.5,lc=:black, la=1);
plot!(getellipsepoints(0, 0, a[1], a[2], 0), lw=1, lc=:white, la=0.7);
P_XY_EH = Plots.current();

XY_diff = [norm(H[i] .- VF[i]) ./ norm(H[i]) for i in CartesianIndices(H)][:,:,14]' .* 100
Plots.heatmap(grid[1],grid[2],XY_diff,aspect_ratio=:equal,clims=c_diff,cmap=:vik,xlabel = "X",ylabel = "Y", title="EH Diff in %");
plot!(getellipsepoints(0, 0, a[1], a[2], 0), lw=1, lc=:black, la=0.7);
P_XY_diff = Plots.current();

# YZ 

Nyz = hcat([[grid[2][i],grid[3][j]] for i=1:gridSizeY for j=1:gridSizeZ]...)
H_Z = [itpH_Z(0,y,z) for y in grid[2], z in grid[3]]
H_Y = [itpH_Y(0,y,z) for y in grid[2], z in grid[3]]

Plots.heatmap(grid[2],grid[3],normH[26,:,:]',aspect_ratio=:equal,clims=c,cmap=:bam,xlabel = "Y",ylabel = "Z", title="SimField");
arrow0!.(Nyz[1,:],Nyz[2,:],sc_ar .* vec(permutedims(H_Y,[2,1])),sc_ar .* vec(permutedims(H_Z,[2,1])); as=0.3, lw=0.5,lc=:black, la=1);
P_YZ = Plots.current();

Plots.heatmap(grid[2],grid[3],normH_EH[26,:,:]',aspect_ratio=:equal,clims=c,cmap=:bam,xlabel = "Y",ylabel = "Z", title="EH_Field");
arrow0!.(Nyz[1,:],Nyz[2,:],sc_ar .* vec(permutedims((h->h[2]).(VF)[26,:,:],[2,1])),sc_ar .* vec(permutedims((h->h[3]).(VF)[26,:,:],[2,1])); as=0.3, lw=0.5,lc=:black, la=1);
plot!(getellipsepoints(0, 0, a[2], a[3], 0), lw=1, lc=:white, la=0.7);
P_YZ_EH = Plots.current();

YZ_diff = (normH[26,:,:]' .- normH_EH[26,:,:]') ./ normH[26,:,:]' .* 100
Plots.heatmap(grid[2],grid[3],YZ_diff,aspect_ratio=:equal,clims=c_diff,cmap=:vik,xlabel = "Y",ylabel = "Z", title="EH Diff in %");
plot!(getellipsepoints(0, 0, a[2], a[3], 0), lw=1, lc=:black, la=0.7);
P_YZ_diff = Plots.current();

# XZ 

Nxz = hcat([[grid[1][i],grid[3][j]] for i=1:gridSizeX for j=1:gridSizeZ]...)
H_Z = [itpH_Z(x,0,z) for x in grid[1], z in grid[3]]
H_X = [itpH_X(x,0,z) for x in grid[1], z in grid[3]]

Plots.heatmap(grid[1],grid[3],normH[:,14,:]',aspect_ratio=:equal,clims=c,cmap=:bam,xlabel = "X",ylabel = "Z", title="SimField");
arrow0!.(Nxz[1,:],Nxz[2,:],sc_ar .* vec(permutedims(H_X,[2,1])),sc_ar .* vec(permutedims(H_Z,[2,1])); as=0.3, lw=0.5,lc=:black, la=1);
P_XZ = Plots.current();

Plots.heatmap(grid[1],grid[3],normH_EH[:,14,:]',aspect_ratio=:equal,clims=c,cmap=:bam,xlabel = "X",ylabel = "Z", title="EH_Field");
arrow0!.(Nxz[1,:],Nxz[2,:],sc_ar .* vec(permutedims((h->h[1]).(VF)[:,14,:],[2,1])),sc_ar .* vec(permutedims((h->h[3]).(VF)[:,14,:],[2,1])); as=0.3, lw=0.5,lc=:black, la=1);
plot!(getellipsepoints(0, 0, a[1], a[3], 0), lw=1, lc=:white, la=0.7);
P_XZ_EH = Plots.current();

XZ_diff = (normH[:,14,:]' .- normH_EH[:,14,:]') ./ normH[:,14,:]' .* 100
Plots.heatmap(grid[1],grid[3],XZ_diff,aspect_ratio=:equal,clims=c_diff,cmap=:vik,xlabel = "X",ylabel = "Z", title="EH Diff in %");
plot!(getellipsepoints(0, 0, a[1], a[3], 0), lw=1, lc=:black, la=0.7);
P_XZ_diff = Plots.current();

h2 = scatter([0,0], [0,1], zcolor=[c[1],c[2]], clims=c.*1000,
                 xlims=(1,1.1), label="", framestyle=:none, c=:bam, cbar=:true);
h3 = scatter([0,0], [0,1], zcolor=[c_diff[1],c_diff[2]], clims=c_diff,
                 xlims=(1,1.1), label="", framestyle=:none, c=:vik, cbar=:true);                 

l = @layout [a{0.04w} Plots.grid(3, 3) b{0.04w}];
plot(h2,P_XY,P_YZ,P_XZ,P_XY_EH,P_YZ_EH,P_XZ_EH,P_XY_diff,P_YZ_diff,P_XZ_diff,h3,layout=l, link=:all, size=(1000,700))

savefig("results/B0Field.png")


