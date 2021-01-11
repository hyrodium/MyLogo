using LinearAlgebra, Luxor, Colors

struct Point2D
    x::Float64
    y::Float64
end

struct Vertex
    i1::Int
end

struct Edge
    i1::Int
    i2::Int
    function Edge(i1,i2)
        return new(min(i1,i2), max(i1,i2))
    end
end

struct Face
    i1::Int
    i2::Int
    i3::Int
    function Face(i1,i2,i3)
        return new(sort([i1,i2,i3])...)
    end
end

Base.:(==)(v1::Vertex, v2::Vertex) = v1.i1 == v2.i1
Base.:(==)(e1::Edge, e2::Edge) = Set((e1.i1, e1.i2)) == Set((e2.i1, e2.i2))
Base.:(==)(f1::Face, f2::Face) = Set((f1.i1, f1.i2, f1.i3)) == Set((f2.i1, f2.i2, f2.i3))

Luxor.Point(p::Point2D) = Point(p.x,p.y)

Base.:-(p1::Point2D, p2::Point2D) = Point2D(p1.x-p2.x, p1.y-p2.y)

function _homogeneouscoordinate(p::Point2D, t::Face, ps)
    p1 = ps[t.i1]
    p2 = ps[t.i2]
    p3 = ps[t.i3]

    v = [p.x, p.y, 1.0]
    A = [p1.x p2.x p3.x
         p1.y p2.y p3.y
         1.0  1.0  1.0]
    return A\v
end

function _pin(p::Point2D, t::Face, ps)
    hc = _homogeneouscoordinate(p,t, ps)
    return minimum(hc) > 0
end

function draw(ps::Vector{Point2D}, vs::Vector{Vertex} ,es::Vector{Edge} ,fs::Vector{Face}; withcircle=true, name=nothing)
    @svg begin
        d = 1000
        if isnothing(name)
            Drawing(d, d)
        else
            Drawing(d, d, name)
        end
        origin()
        translate(-500, -500)
        background(RGB(1,1,1))

        for f in fs
            p1 = Point(ps[f.i1])
            p2 = Point(ps[f.i2])
            p3 = Point(ps[f.i3])
            sethue(RGB(0,1,0))
            poly([p1,p2,p3,p1],:fill)
        end

        for f in fs
            p1 = Point(ps[f.i1])
            p2 = Point(ps[f.i2])
            p3 = Point(ps[f.i3])
            setline(20)
            sethue(RGB(1,1,1))
            poly([p1,p2,p3,p1],:stroke)
        end

        if withcircle
            setline(3)
            sethue(RGB(1,0,1))
            for f in fs
                p1 = Point(ps[f.i1])
                p2 = Point(ps[f.i2])
                p3 = Point(ps[f.i3])
                circle(p1, p2, p3, :stroke)
            end
            # f = fs[1]
            # p1 = Point(0,0)
            # p2 = Point(d,0)
            # p3 = Point(ps[f.i3])
            # circle(p1, p2, p3, :stroke)
        end

        setline(3)
        sethue(RGB(0,1,1))
        for e in es
            line(Point(ps[e.i1]), Point(ps[e.i2]), :stroke)
        end

        setcolor(RGB(1,0,0))
        for v in vs
            circle(Point(ps[v.i1]), 10,:fill)
        end

        finish()
        preview()
    end
end

LinearAlgebra.norm(p::Point2D) = √(p.x^2+p.y^2)
LinearAlgebra.dot(p1::Point2D, p2::Point2D) = p1.x*p2.x + p1.y*p2.y

## ここから試し① とりあえず全結合

# N = 10
# ps = [Point2D(rand(1:d), rand(1:d)) for i in 1:N]
# vs = [Vertex(i) for i in 1:N]
# es = reshape([Edge(i,j) for i in 1:N, j in 1:N],N^2)
# fs = Face[]
# draw(ps, vs, es, fs)

## ここから試し②

# N = 10
# ps = [Point2D(rand(1:d), rand(1:d)) for i in 1:N]
# vs = [Vertex(i) for i in 1:N]
# es = Edge[]
# fs = Face[]

# index = 2
# mindistance = norm(ps[2]-ps[1])
# for i in 3:N
#     l = norm(ps[i]-ps[1])
#     if l < mindistance
#         index = i
#         mindistance = l
#     end
# end

# push!(es,Edge(1,index))

# draw(ps, vs, es, fs)


function edges(t::Face)
    es_tmp = Edge[]
    I = [t.i1,t.i2,t.i3]
    for e in es
        J = [e.i1,e.i2]
        if J ⊆ I
            push!(es_tmp, e)
        end
    end
    return es_tmp
end

function faces(s::Edge)
    ts_tmp = Face[]
    for f in fs
        if s ∈ edges(f)
            push!(ts_tmp, f)
        end
    end
    return ts_tmp
end

function _vertices(e::Edge)
    return [Vertex(e.i1), Vertex(e.i2)]
end

function _vertices(f::Face)
    return [Vertex(f.i1), Vertex(f.i2), Vertex(f.i3)]
end

function _isvalid(e::Edge, fs, ps)
    f_touch = _faces(e, fs)
    if length(f_touch) == 1
        return true
    elseif length(f_touch) == 2
        v2, v4 = _vertices(e)
        vs_tmp = unique(vcat(_vertices.(f_touch)...))
        v1, v3 = setdiff(vs_tmp, [v2, v4])
        p1, p2, p3, p4 = ps[v1.i1], ps[v2.i1], ps[v3.i1], ps[v4.i1]
        α = acos(dot(p2-p1, p4-p1)/norm(p2-p1)/norm(p4-p1))
        β = acos(dot(p2-p3, p4-p3)/norm(p2-p3)/norm(p4-p3))
        return α+β < π
    else
        error("eeerrr")
    end
end

function _flip!(e::Edge, es, fs, es_vague)
    f_touch = _faces(e, fs)
    v2, v4 = _vertices(e)
    vs_tmp = unique(vcat(_vertices.(f_touch)...))
    v1, v3 = setdiff(vs_tmp, [v2, v4])

    i1 = v1.i1
    i2 = v2.i1
    i3 = v3.i1
    i4 = v4.i1

    setdiff!(es, [e])
    setdiff!(es_vague, [e])
    push!(es, Edge(i1,i3))
    setdiff!(fs, f_touch)
    push!(fs, Face(i1,i2,i3), Face(i1,i3,i4))
    push!(es_vague, Edge(i1, i2), Edge(i2, i3), Edge(i3, i4), Edge(i1, i4))
    return
end

## ここから試し③
begin
    N = 8
    d = 1000
    ps = [Point2D(0,0), Point2D(1000,0), Point2D(0,1000), Point2D(1000,1000)]
    vs = [Vertex(i) for i in 1:N]
    es = [Edge(1,2),Edge(1,3),Edge(2,4),Edge(3,4)]

    ps = push!(ps,Point2D(rand(1:d), rand(1:d)))
    # ps = push!(ps,ppp[5])
    es = append!(es, [Edge(5,1),Edge(5,2),Edge(5,3),Edge(5,4)])
    fs = [Face(1,2,5),Face(1,3,5),Face(2,4,5),Face(3,4,5)]

    for i in 6:N
        new_point = Point2D(rand(1:d), rand(1:d))
        # new_point = ppp[i]
        ps = push!(ps,new_point)

        f_l = length(fs)
        for f_i in 1:f_l
            f = fs[f_i]
            if ps[end] in f
                i1 = f.i1
                i2 = f.i2
                i3 = f.i3

                global es_vague = edges(f)
                push!(es,Edge(i1,i),Edge(i2,i),Edge(i3,i))
                push!(fs,Face(i1,i2,i),Face(i2,i3,i),Face(i1,i3,i))
                deleteat!(fs,f_i)
                break
            end
        end
        while true
            if isempty(es_vague)
                println("finished!")
                break
            end

            l = length(es_vague)
            e = es_vague[end]
            if _isvalid(e)
                pop!(es_vague)
            else
                flip!(e)
            end
        end
    end

    draw(ps, vs, es, fs)
end

draw(ps, vs, es, fs, withcircle=false, name="D.png")
draw(ps, vs, es, fs, name="C.png")

length(vs) - length(es) + length(fs)

function _edges(f::Face, es::Vector{Edge})
    es_tmp = Edge[]
    i_face = [f.i1,f.i2,f.i3]
    for e in es
        i_edge = [e.i1,e.i2]
        if i_edge ⊆ i_face
            push!(es_tmp, e)
        end
    end
    return es_tmp
end

function _faces(e::Edge, fs::Vector{Face})
    fs_tmp = Face[]
    i_edge = [e.i1, e.i2]
    for f in fs
        i_face = [f.i1,f.i2,f.i3]
        if i_edge ⊆ i_face
            push!(fs_tmp, f)
        end
    end
    return fs_tmp
end


function _vertices(e::Edge)
    return [Vertex(e.i1), Vertex(e.i2)]
end

function _vertices(f::Face)
    return [Vertex(f.i1), Vertex(f.i2), Vertex(f.i3)]
end

## adfa

function adfa()
    N = 8
    d = 1000
    _ps_ = [Point2D(0,0), Point2D(1000,0), Point2D(0,1000), Point2D(1000,1000)]
    vs = [Vertex(i) for i in 1:N]
    _es = [Edge(1,2),Edge(1,3),Edge(2,4),Edge(3,4)]

    push!(_ps_,Point2D(rand(1:d), rand(1:d)))
    append!(_es, [Edge(5,1),Edge(5,2),Edge(5,3),Edge(5,4)])
    _fs = [Face(1,2,5),Face(1,3,5),Face(2,4,5),Face(3,4,5)]

    for i in 6:N
        new_point = Point2D(rand(1:d), rand(1:d))
        _ps_ = push!(_ps_, new_point)

        f_l = length(_fs)
        local _es_vague
        for f_i in 1:f_l
            f = _fs[f_i]
            if _pin(_ps_[end], f, _ps_)
                i1 = f.i1
                i2 = f.i2
                i3 = f.i3

                _es_vague = _edges(f, _es)
                push!(_es,Edge(i1,i),Edge(i2,i),Edge(i3,i))
                push!(_fs,Face(i1,i2,i),Face(i2,i3,i),Face(i1,i3,i))
                deleteat!(_fs,f_i)
                break
            end
        end
        while true
            if isempty(_es_vague)
                break
            end

            l = length(_es_vague)
            e = _es_vague[end]
            if _isvalid(e, _fs, _ps_)
                pop!(_es_vague)
            else
                _flip!(e, _es, _fs, _es_vague)
            end
        end
    end

    draw(_ps_, vs, _es, _fs)
end

adfa()

