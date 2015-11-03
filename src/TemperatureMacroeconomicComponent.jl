using Mimi

@defcomp temperaturemacroeconomic begin
    regions = Index()

    baseline = Variable(index=[regions])
    changeypcgrowth = Variable(index=[time,regions])

    macrocoeffs::Vector{Float64} = Parameter()
    temp90 = Parameter(index=[regions])
    regtmp = Parameter(index=[time,regions])
end

function timestep(s::temperaturemacroeconomic, t::Int)
    v = s.Variables
    p = s.Parameters
    d = s.Dimensions

    if t == 1
        for r in d.regions
            v.changeypcgrowth[t, r] = 1
            v.baseline[r] = sum(p.macrocoeffs .* [p.regtmp[1, r] + p.temp90[r], (p.regtmp[1, r] + p.temp90[r])^2])
        end
    else
        for r in d.regions
            v.changeypcgrowth[t, r] = exp(sum(p.macrocoeffs .* [p.regtmp[t, r] + p.temp90[r], (p.regtmp[t, r] + p.temp90[r])^2]) - v.baseline[r])
        end
    end
end
