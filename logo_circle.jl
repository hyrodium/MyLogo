using Images
using Luxor
using Statistics

using ImageFiltering, ColorVectorSpace

## median image
img = load("hyrodium.jpg")
d = 1000
img = imresize(img,(d,d))

get_r(c) = c.r
get_g(c) = c.g
get_b(c) = c.b

colors_r = get_r.(img)
colors_g = get_g.(img)
colors_b = get_b.(img)

r = mapwindow(ImageFiltering.median!, colors_r, (11, 11))
g = mapwindow(ImageFiltering.median!, colors_g, (11, 11))
b = mapwindow(ImageFiltering.median!, colors_b, (11, 11))

img = [RGB(r[i,j], g[i,j], b[i,j]) for i in 1:d, j in 1:d]
img_interpolated = ImageTransformations.box_extrapolation(img)
save("1000.png",img)

f(z::Union{Complex,Real})::Complex=exp(im*pi/4)*(exp(z)-1)/(exp(z)+1)
f′(z::Union{Complex,Real})::Complex=exp(im*pi/4)*2exp(z)/(exp(z)+1)^2

lxrpt(w::Complex)=Luxor.Point(500*real(w),-500*imag(w))
imgpt(w::Complex)=(999/2-999/2*imag(w)+1,999/2*real(w)+999/2+1)

struct Cell
    m::Int  # fineness
    i::Int  # i-index
    j::Int  # j-index
end

function sampling_z(c::Cell)
    Δ = π/2c.m
    i, j = c.i, c.j
    zs = [z_r+im*z_i for z_r in range(Δ*i,Δ*(i+1),length=27), z_i in range(Δ*j,Δ*(j+1),length=27)]
    return zs
end

function sampling_c(c::Cell)
    zs = sampling_z(c)
    cs = [img_interpolated[i...] for i in imgpt.(f.(zs))]
    return cs
end

function meancolor(c::Cell)
    cs = sampling_c(c)
    return mean(cs)
end

function isflat(c::Cell)
    # if c.m>31
    #     return true
    # end
    if sizeof(c) < 0.02
        return true
    end
    gs = Gray.(sampling_c(c))
    max_gray = maximum(Float64.(gs))
    min_gray = minimum(Float64.(gs))
    diff = (max_gray-min_gray)/(min_gray+0.2)
    return diff<0.5
end

function sizeof(c::Cell)
    Δ = π/2c.m
    i, j = c.i, c.j
    z = Δ*(i+1/2) + im*Δ*(j+1/2)
    return abs(Δ*f′(z))
end

"""
control points of Bézier curve from given function and range
"""
function BézPts(f,t0,t1)
    f(t0),
    3*(f(2*t0/3+t1/3)-(8*f(t0)+f(t1))/27)-3*(f(t0/3+2*t1/3)-(f(t0)+8*f(t1))/27)/2,
    -3*(f(2*t0/3+t1/3)-(8*f(t0)+f(t1))/27)/2+3*(f(t0/3+2*t1/3)-(f(t0)+8*f(t1))/27),
    f(t1)
end

function drawcell(c::Cell)
    Δ = π/2c.m
    i, j = c.i, c.j
    real_min = Δ*i
    real_max = Δ*(i+1)
    imag_min = Δ*j
    imag_max = Δ*(j+1)

    # setcolor
    Luxor.setcolor(meancolor(c))

    # polygon
    # zs1 = [z_r+im*imag_min for z_r in range(real_min,real_max,length=25)]
    # zs2 = [real_max+im*z_i for z_i in range(imag_min,imag_max,length=25)]
    # zs3 = [z_r+im*imag_max for z_r in range(real_max,real_min,length=25)]
    # zs4 = [real_min+im*z_i for z_i in range(imag_max,imag_min,length=25)]
    # zs = vcat(zs1,zs2,zs3,zs4)
    # ps = lxrpt.(f.(zs))
    # Luxor.poly(ps,:fill)

    # bezier
    path = BezierPath([
        BezierPathSegment(lxrpt.(BézPts(z_r->f(z_r+im*imag_min),real_min,real_max))...),
        BezierPathSegment(lxrpt.(BézPts(z_i->f(real_max+im*z_i),imag_min,imag_max))...),
        BezierPathSegment(lxrpt.(BézPts(z_r->f(z_r+im*imag_max),real_max,real_min))...),
        BezierPathSegment(lxrpt.(BézPts(z_i->f(real_min+im*z_i),imag_max,imag_min))...)
    ])
    drawbezierpath(path,:fill,close=true)
    Luxor.setcolor(RGB(0.1,0.05,0.05))
    drawbezierpath(path,:stroke,close=true)

    # circle
    # p1 = lxrpt(f((real_min+real_max)/2+im*imag_min))
    # p2 = lxrpt(f(real_min+im*(imag_min+imag_max)/2))
    # p3 = lxrpt(f((real_min+real_max)/2+im*imag_max))
    # circle(p1,p2,p3,:fill)
end

function drawcaps()
    Luxor.setcolor(RGB(0.1,0.05,0.05))
    r = (π/2)*3

    p1 = lxrpt(f(r))
    p2 = lxrpt(f(r+im))
    p3 = lxrpt(f(r-im))
    q1 = lxrpt(f(-r))
    q2 = lxrpt(f(-r+im))
    q3 = lxrpt(f(-r-im))

    circle(p1,p2,p3,:clip)
    circle(O,500,:fill)
    circle(O,500,:stroke)
    clipreset()
    circle(q1,q2,q3,:clip)
    circle(O,500,:fill)
    circle(O,500,:stroke)
    clipreset()
end

function draw(cs;filename="MyLogo.png",size=1000)
    Drawing(size, size, filename)
    setline(1.2)
    w = _img.width
    h = _img.height
    origin()
    for c in cs
        drawcell(c)
    end
    sethue(RGB(1,0,0))
    drawcaps()
    finish()
end

function splitcell(c::Cell)
    m = c.m
    i = c.i
    j = c.j
    m′ = 2m
    return [Cell(m′,2i,2j),Cell(m′,2i+1,2j),Cell(m′,2i,2j+1),Cell(m′,2i+1,2j+1)]
end

cs_tmp = vec([Cell(1,i,j) for i in -3:2, j in -1:0])
cs = Cell[]

while !isempty(cs_tmp)
    c = pop!(cs_tmp)
    if isflat(c)
        push!(cs,c)
    else
        append!(cs_tmp,splitcell(c))
    end
end

draw(cs,filename="MyLogo.png",size=1024)
