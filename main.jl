using DataFrames
using CSV
using PyPlot

dfcounts = DataFrame("d" => [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.75,1.0],"C" =>
[901,652,443,339,283,281,240,220,180,154], "σ" => [30.0,25.5,21.0,18.4,16.8,16.8,15.5,14.8,13.4,
12.4])

CSV.write("gammacounts.csv",dfcounts)

dfcountswithx = transform(dfcounts,:d => (x -> @. 1/x^2) => :x)
dfcountswithy = transform(dfcountswithx,:C => (x -> x) => :y)

function fitcountstoline(data)
    Δ = sum(@. 1/data.σ^2 )
    Δ = Δ * sum(@. data.x^2 / data.σ^2)
    Δ = Δ - (sum(@. data.x / data.σ^2)) ^ 2

    oneoverdelta = (1 / Δ)
    a = sum(@. data.x^2/data.σ^2) * sum(@. data.y / data.σ^2)
    a = a - sum(@. data.x / data.σ^2) * sum(@. data.x * data.y / data.σ^2)
    a = a * oneoverdelta    

    b = sum(@. 1/data.σ^2) * sum(@. data.x * data.y / data.σ^2)
    b = b - sum(@. data.x /data.σ^2) * sum(@. data.y / data.σ^2)
    b = oneoverdelta * b

    w = @. 1 / data.σ^2
    σₐ² = oneoverdelta * sum(@. w * data.x^2) 
    σb² = oneoverdelta * sum(w)

    σₐ = sqrt(σₐ²)
    σb = sqrt(σb²)
    println("line fit with coefficients a:($a) and b:($b) uncertainy a:($σₐ) and b:($σb)")
    return a,σₐ,b,σb
end

a,σₐ,b,σb = fitcountstoline(dfcountswithy)

line(a,b,x) = a * x + b
f = figure()

plot(dfcountswithy.x,dfcountswithy.y,color="blue",linewidth=2.0,linestyle="--")

xlabel(L"distance in m")
ylabel(L"counts")

title("counts vs distance from source")

savefig("gammactsplot.png")
savefig("gammactsplot.pdf")

close(fig)