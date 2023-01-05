using DataFrames
using CSV
using PyPlot

dfcounts = DataFrame("d" => [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.75,1.0],"C" =>
[901,652,443,339,283,281,240,220,180,154], "σ" => [30.0,25.5,21.0,18.4,16.8,16.8,15.5,14.8,13.4,
12.4])

CSV.write("gammacounts.csv",dfcounts)

dfcountswithx = transform(dfcounts,:d => (x -> @. 1/x^2) => :x)
dfcountswithy = transform(dfcountswithx,:C => (x -> x) => :y)

struct Linefit
    a
    σₐ
    b
    σb
    function Linefit(a,siga,b,sigb)
        println("a:$a,siga:$siga,b:$b,sigb:$sigb")
        new(a,siga,b,sigb)
    end
end
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
    #println("line fit with coefficients a:($a) and b:($b) uncertainy a:($σₐ) and b:($σb)")
    return Linefit(a,σₐ,b,σb)
end

Gamma_Line_Fit = fitcountstoline(dfcountswithy)


line(a,b,x) = @. b * x + a
function line(lf::Linefit,x)
    @. lf.b * x + lf.a 
end

x = dfcountswithy.x 
y = dfcountswithy.y

function plotfittedgammacounts(x,y,model,filename::String)
    f = figure()

    minx = minimum(x)
    maxx = maximum(x)

    r = range(minx,maxx,length(x) * 2)
    modelx = collect(r)
    modely = model(Gamma_Line_Fit,modelx)

    plot(x,y,color="blue",linewidth=2.0,linestyle="--",
    marker="o", label=L"gamma counts")
    plot(modelx,modely, color="red",linewidth=2.0,linestyle="-",
    marker="x", label=L"fitted counts")

    xlabel("distance in m")
    ylabel("counts")

    legend(loc="upper right",fontsize="x-large")

    title("counts vs distance from source")

    savefig(filename * ".png")

    close(f)
end

plotfittedgammacounts(x,y,line,"gammacounts2")