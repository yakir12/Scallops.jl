using Luxor, QRCode

qrc = qrcode("https://github.com/yakir12/Scallops.jl")

@pdf begin
    tiles = Tiler(250, 250, size(qrc)..., margin=0)
    setline(0.5)
    squares = first.(tiles)
    @layer begin
        for (n, pt) in enumerate(box(O, 140, 140, vertices=true))
            sethue([Luxor.julia_purple, Luxor.julia_blue,  Luxor.julia_red,  Luxor.julia_green][n])
            setopacity(0.75)
            circle(pt, 50, :fill)
        end
    end
    for x in eachindex(qrc)
        if qrc[x]
            sethue("black")
            box(squares[x], tiles.tilewidth, tiles.tileheight, :fill)
        end
    end
end 250 250 "/home/yakir/documents/formalia/Presentations/ICIV2019/betterposter-latex-template/img/qrcode"
