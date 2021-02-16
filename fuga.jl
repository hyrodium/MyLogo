using LinearAlgebra, Luxor, Colors


@svg begin
    d = 1000
    Drawing(d, d, "logo.svg")
    origin()
    translate(-500, -500)
    N = 10

    ps = [Point(rand(1:d), rand(1:d)) for i in 1:N]
    background(RGB(1,1,1))

    setcolor(RGB(1,0,0))
    for i in 1:N
        circle(ps[i],20,:fill)
    end

    for i in 1:N, j in 1:N
        norm(ps[i]-ps[j])
    end

    finish()
    preview()
end

3

