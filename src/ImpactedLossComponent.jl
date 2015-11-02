using Mimi

@defcomp impactedloss begin
    regions = Index()

    loss = Variable(index=[time,regions])

    income = Parameter(index=[time, regions])
    impactedincome = Parameter(index=[time, regions])
end

function timestep(s::impactedloss, t::Int)
    v = s.Variables
    p = s.Parameters
    d = s.Dimensions

    for r in d.regions
        v.loss[t, r] = -(p.impactedincome[t, r] - p.income[t, r]) * 1e9
    end
end
