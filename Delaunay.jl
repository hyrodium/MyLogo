using LinearAlgebra, Luxor, Colors, Random

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

LinearAlgebra.norm(p::Point2D) = √(p.x^2+p.y^2)
LinearAlgebra.dot(p1::Point2D, p2::Point2D) = p1.x*p2.x + p1.y*p2.y

function _homogeneouscoordinate(p::Point2D, f::Face, ps)
    p1 = ps[f.i1]
    p2 = ps[f.i2]
    p3 = ps[f.i3]

    v = [p.x, p.y, 1.0]
    A = [p1.x p2.x p3.x
         p1.y p2.y p3.y
         1.0  1.0  1.0]
    return A\v
end

function _pinf(p::Point2D, f::Face, ps)
    hc = _homogeneouscoordinate(p,f, ps)
    return minimum(hc) > 0
end

function _pinf(p::Point2D, f::Face, ps)
    p1 = ps[f.i1]
    p2 = ps[f.i2]
    p3 = ps[f.i3]

    S = (p1.x-p2.x)*(p1.y-p3.y)-(p1.y-p2.y)*(p1.x-p3.x)

    k12 = S * ((p1.x-p2.x)*(p1.y-p.y)-(p1.y-p2.y)*(p1.x-p.x))
    k23 = S * ((p2.x-p3.x)*(p2.y-p.y)-(p2.y-p3.y)*(p2.x-p.x))
    k31 = S * ((p3.x-p1.x)*(p3.y-p.y)-(p3.y-p1.y)*(p3.x-p.x))
    println(S)
    println([k12, k23, k31])
    S * ((p1.x-p2.x)*(p1.y-p.y)-(p1.y-p2.y)*(p1.x-p.x)) < 0 && return false
    S * ((p2.x-p3.x)*(p2.y-p.y)-(p2.y-p3.y)*(p2.x-p.x)) < 0 && return false
    S * ((p3.x-p1.x)*(p3.y-p.y)-(p3.y-p1.y)*(p3.x-p.x)) < 0 && return false
    return true
end

function _p_in_f_score(p::Point2D, f::Face, ps)
    hc = _homogeneouscoordinate(p,f, ps)
    return minimum(hc)
end


function indexedcolor(i)
    m = MersenneTwister(i)
    RGB(rand(m),rand(m),rand(m))
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

        for i in 1:length(fs)
            f = fs[i]
            p1 = Point(ps[f.i1])
            p2 = Point(ps[f.i2])
            p3 = Point(ps[f.i3])
            sethue(indexedcolor(i))
            setlinejoin("round")
            poly([p1,p2,p3,p1],:fill)
        end

        for i in 1:length(fs)
            f = fs[i]
            p1 = Point(ps[f.i1])
            p2 = Point(ps[f.i2])
            p3 = Point(ps[f.i3])
            setline(20)
            sethue(RGB(1,1,1))
            poly([p1,p2,p3,p1],:stroke)
        end

        if withcircle
            setline(3)
            for i in 1:length(fs)
                f = fs[i]
                p1 = Point(ps[f.i1])
                p2 = Point(ps[f.i2])
                p3 = Point(ps[f.i3])
                sethue(indexedcolor(i))
                circle(p1, p2, p3, :stroke)
            end
        end

        setline(3)
        sethue(RGB(0.5,0.5,0.5))
        for e in es
            line(Point(ps[e.i1]), Point(ps[e.i2]), :stroke)
        end

        setcolor(RGB(0,0,0))
        for v in vs
            circle(Point(ps[v.i1]), 10,:fill)
        end

        finish()
        preview()
    end
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
        α = acos(clamp(dot(p2-p1, p4-p1)/norm(p2-p1)/norm(p4-p1), -1, 1))
        β = acos(clamp(dot(p2-p3, p4-p3)/norm(p2-p3)/norm(p4-p3), -1, 1))
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

function _face_include_p(p,ps,fs)
    tmp_score = -1.0
    tmp_f = fs[1]
    for f in fs
        if _pinf(p, f, ps)
            return f
        end
    end
    println("best score: $(tmp_score)")
    println("best face: $(tmp_f)")
    error("no triangle")
end

function _face_include_p(p,ps,fs)
    tmp_score = -1.0
    tmp_f = fs[1]
    for f in fs
        score = _p_in_f_score(p, f, ps)
        if score > 0
            return(f)
        elseif score > tmp_score
            tmp_score = score
            tmp_f = f
        end
    end
    println("best score: $(tmp_score)")
    println("best face: $(tmp_f)")
    error("no triangle")
end

function adfa()
    N = 250
    d = 1000
    _ps = [Point2D(0,0), Point2D(1000,0), Point2D(0,1000), Point2D(1000,1000)]
    _vs = [Vertex(i) for i in 1:N]
    _es = [Edge(1,2),Edge(1,3),Edge(2,4),Edge(3,4)]

    push!(_ps,Point2D(rand(1:d), rand(1:d)))
    append!(_es, [Edge(5,1),Edge(5,2),Edge(5,3),Edge(5,4)])
    _fs = [Face(1,2,5),Face(1,3,5),Face(2,4,5),Face(3,4,5)]

    for i in 6:N
        new_point = Point2D(rand(1:d), rand(1:d))
        _ps = push!(_ps, new_point)

        f_l = length(_fs)
        _f = _face_include_p(_ps[end], _ps, _fs)

        i1 = _f.i1
        i2 = _f.i2
        i3 = _f.i3
        _es_vague = _edges(_f, _es)
        push!(_es,Edge(i1,i),Edge(i2,i),Edge(i3,i))
        push!(_fs,Face(i1,i2,i),Face(i2,i3,i),Face(i1,i3,i))
        setdiff!(_fs,[_f])

        while true
            if isempty(_es_vague)
                break
            end

            l = length(_es_vague)
            e = _es_vague[end]
            if _isvalid(e, _fs, _ps)
                pop!(_es_vague)
            else
                _flip!(e, _es, _fs, _es_vague)
            end
        end
    end

    draw(_ps, _vs, _es, _fs, name="save.png")
end

adfa()



exit()




