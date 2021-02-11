using Images
using ImageFiltering, ColorVectorSpace
push!(LOAD_PATH, "../VisualizingDelaunay.jl/")
using VisualizingDelaunay


img = load("hyrodium.jpg")

d = 1000
img = imresize(img,(d,d))

get_r(c) = c.r
get_g(c) = c.g
get_b(c) = c.b

colors_r = get_r.(img)
colors_g = get_g.(img)
colors_b = get_b.(img)

r = mapwindow(ImageFiltering.median!, colors_r, (5, 5))
g = mapwindow(ImageFiltering.median!, colors_g, (5, 5))
b = mapwindow(ImageFiltering.median!, colors_b, (5, 5))

[RGB(r[i,j], g[i,j], b[i,j]) for i in 1:d, j in 1:d]


function clp(c)
    return clamp(c,Gray(0.0),Gray(1.0))
end

img_LoG_r = imfilter(gray_r, Kernel.LoG(2))*50

img_LoG_r .+ Gray(0.5)

img_LoG_r.^2


kernel_x = [0 0 0;-1 0 1;0 0 0]
img_x = clp.(imfilter(img_LoG_r, kernel_x))

kernel_y = [0 -1 0;0 0 0;0 1 0]
img_y = clp.(imfilter(img_LoG_r, kernel_y))

E = -Float64.(img_x.^2 + img_y.^2)
Gray.(-E)

sortedindex = sortperm(E[:])

function energy(i::Integer)
    return E[i]
end

function randnext(::Integer)
    return sortedindex[rand(1:200000)]
end

# mcmc part
function mcmcnext(X, randnext, energy::Function; β=1.0)
    X′ = randnext(X)
    E = energy(X)
    E′ = energy(X′)
    exp(-β*E′)/exp(-β*E) < rand() && return X
    return X′
end

function iterate_mcmc(X_init, randnext, energy::Function, N::Int; β=1.0)
    X_tmp = X_init
    Xs = Array{typeof(X_init),1}(undef, N)
    for i in 1:N
        X_tmp = mcmcnext(X_tmp, randnext, energy, β=β)
        Xs[i] = X_tmp
    end
    return Xs
end

xx = iterate_mcmc(1, randnext, energy, 100000, β=8.0)

x = [Gray(0.0) for i in 1:d, j in 1:d]
for I in unique(xx)
    x[I...] = 1.0
end
x

xxx = unique(xx)

N = 3000

reduced = Array{Int}(undef,N)

reduced[1] = rand(xxx)
for i in 2:N
    reduced[i] = rand(xxx)
end

x = [Gray(0.0) for i in 1:d, j in 1:d]
for I in reduced
    x[I...] = 1.0
end
x



x = [Gray(0.0) for i in 1:d, j in 1:d]
for I in reduced
    x[mod(I,d), I÷d] = 1.0
end
x


ps = [VisualizingDelaunay.Point2D(I÷d, mod(I,d)) for I in reduced]

VisualizingDelaunay.main(ps)

VisualizingDelaunay.Mesh2D()


tess = DelaunayTessellation(1)


push!(tess, VoronoiDelaunay.Point(1.5, 1.5))

@less DelaunayTessellation()




E
E[sortedindex[1]]
E[sortedindex[2]]
E[sortedindex[3]]
E[sortedindex[4]]
E[sortedindex[40]]

E[sortedindex[100]]

E[sortedindex[200]]

E[sortedindex[4000]]


length(sortedindex)


using VoronoiDelaunay
tess = VoronoiDelaunay.from_image(img, 25000)

img


import VoronoiDelaunay: from_image, voronoiedges, getplotxy
import Images: imread
import Gadfly: set_default_plot_size, plot, Geom, Scale, cm, draw, SVG, inch

img = imread("julia.png")
# placing points in places that represent the image.
tess = from_image(img, 25000)

# making the plot
set_default_plot_size(10cm,10cm)
x, y = getplotxy(voronoiedges(tess))
p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.25, maxvalue=1.75))

# save as SVG
draw(SVG("voroimage.svg", 8inch, 4inch), p)

import VoronoiDelaunay: from_image, voronoiedges, getplotxy
import Images: imread
import Gadfly: set_default_plot_size, plot, Geom, Scale, cm, draw, SVG, inch





img = load("julia.png")
# placing points in places that represent the image.
tess = from_image(img, 25000)

# making the plot
set_default_plot_size(30cm,15cm)
# x, y = getplotxy(voronoiedges(tess))
# p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.25, maxvalue=1.75))




# save as SVG
draw(SVG("voroimage.svg", 8inch, 4inch), p)



