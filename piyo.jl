using VoronoiDelaunay
tess = DelaunayTessellation()
push!(tess, VoronoiDelaunay.Point(1.5, 1.5))


using VoronoiDelaunay
tess = from_image(img, 25000)

width = max_coord - min_coord
a = Point2D[Point(min_coord + rand() * width, min_coord + rand() * width) for i in 1:100]
push!(tess, a)


Point2D(1,2)

a = [Point2D(1,2),Point2D(1,3), Point2D(2,2)]
DelaunayTessellation(a)

tess = DelaunayTessellation2D(3)

push!(tess, a)


3

