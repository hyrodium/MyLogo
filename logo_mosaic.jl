using Images
using ImageFiltering, ColorVectorSpace
using Statistics

img = load("hyrodium.jpg")

N = 9

d = 2^N
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

th = 0.1

img

norm(img.-mean(img))

a = maximum(Gray.(img.-mean(img)))
b = minimum(Gray.(img.-mean(img)))

function isflat(img; threshold=0.1)
    mean_color = mean(img)
    gray = Gray.(img.-mean_color)
    return maximum(gray)-minimum(gray)<threshold
end

isflat(img)
isflat(img[1:d÷2,1:d÷2])

function quadrisection(img)
    d,_ = size(img)
    return (img[1:d÷2,1:d÷2], img[d÷2+1:d,1:d÷2], img[1:d÷2,d÷2+1:d], img[d÷2+1:d,d÷2+1:d])
end

quadrisection(img)[4]

isflat(quadrisection(img)[1], threshold=0.955)

isflat(img, threshold=0.955)


function hoge(img; threshold=0.1)
    img_new = similar(img)
    d,_ = size(img)
    if isflat(img, threshold=threshold)
        img_new .= mean(img)
    else
        img1, img2, img3, img4 = quadrisection(img)
        img_new[1:d÷2,1:d÷2] = hoge(img1,threshold=threshold)
        img_new[d÷2+1:d,1:d÷2] = hoge(img2,threshold=threshold)
        img_new[1:d÷2,d÷2+1:d] = hoge(img3,threshold=threshold)
        img_new[d÷2+1:d,d÷2+1:d] = hoge(img4,threshold=threshold)
    end
    return img_new
end

for i in 1:100
    im = similar(img)
    save("similar$(lpad(i,3,'0')).png",im)
end

for i in 1:20
    t = 0.8*(1.001-i/20)
    _img = hoge(img,threshold=t)
    save("iimag$(lpad(i,3,'0')).png",_img)
end
