include("helper.jl")

include("SocioEconomicComponent.jl")
include("ImpactedSocioEconomicComponent.jl")
include("PopulationComponent.jl")
include("EmissionsComponent.jl")
include("GeographyComponent.jl")
include("ScenarioUncertaintyComponent.jl")
include("ClimateCO2CycleComponent.jl")
include("ClimateCH4CycleComponent.jl")
include("ClimateN2OCycleComponent.jl")
include("ClimateSF6CycleComponent.jl")
include("ClimateForcingComponent.jl")
include("ClimateDynamicsComponent.jl")
include("BioDiversityComponent.jl")
include("ClimateRegionalComponent.jl")
include("OceanComponent.jl")
include("TemperatureMacroeconomicComponent.jl")
include("ImpactedLossComponent.jl")

function constructfund(;nsteps=1050)
    m = Model()

    setindex(m, :time, collect(1950:1950+nsteps))
    setindex(m, :regions, ["USA", "CAN", "WEU", "JPK", "ANZ", "EEU", "FSU", "MDE", "CAM", "LAM", "SAS", "SEA", "CHI", "MAF", "SSA", "SIS"])

    # ---------------------------------------------
    # Create components
    # ---------------------------------------------
    addcomponent(m, scenariouncertainty, :scenariouncertainty)
    addcomponent(m, population, :population)
    addcomponent(m, geography, :geography)
    addcomponent(m, socioeconomic, :socioeconomic)
    addcomponent(m, impactedsocioeconomic, :impactedsocioeconomic)
    addcomponent(m, emissions, :emissions)
    addcomponent(m, climateco2cycle, :climateco2cycle)
    addcomponent(m, climatech4cycle, :climatech4cycle)
    addcomponent(m, climaten2ocycle, :climaten2ocycle)
    addcomponent(m, climatesf6cycle, :climatesf6cycle)
    addcomponent(m, climateforcing, :climateforcing)
    addcomponent(m, climatedynamics, :climatedynamics)
    addcomponent(m, biodiversity, :biodiversity)
    addcomponent(m, climateregional, :climateregional)
    addcomponent(m, ocean, :ocean)
    addcomponent(m, temperaturemacroeconomic, :temperaturemacroeconomic)
    addcomponent(m, impactedloss, :impactedloss)

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    setparameter(m, :geography, :landloss, zeros(nsteps+1, 16))

    connectparameter(m, :population, :pgrowth, :scenariouncertainty, :pgrowth)
    setparameter(m, :population, :enter, zeros(nsteps+1, 16))
    setparameter(m, :population, :leave, zeros(nsteps+1, 16))
    setparameter(m, :population, :dead, zeros(nsteps+1, 16))

    connectparameter(m, :socioeconomic, :area, :geography, :area)
    connectparameter(m, :socioeconomic, :globalpopulation, :population, :globalpopulation)
    connectparameter(m, :socioeconomic, :populationin1, :population, :populationin1)
    connectparameter(m, :socioeconomic, :population, :population, :population)
    connectparameter(m, :socioeconomic, :pgrowth, :scenariouncertainty, :pgrowth)
    connectparameter(m, :socioeconomic, :ypcgrowth, :scenariouncertainty, :ypcgrowth)
    setparameter(m, :socioeconomic, :eloss, zeros(nsteps+1, 16))
    setparameter(m, :socioeconomic, :sloss, zeros(nsteps+1, 16))
    connectparameter(m, :socioeconomic, :mitigationcost, :emissions, :mitigationcost)

    connectparameter(m, :impactedsocioeconomic, :area, :geography, :area)
    connectparameter(m, :impactedsocioeconomic, :globalpopulation, :population, :globalpopulation)
    connectparameter(m, :impactedsocioeconomic, :populationin1, :population, :populationin1)
    connectparameter(m, :impactedsocioeconomic, :population, :population, :population)
    connectparameter(m, :impactedsocioeconomic, :pgrowth, :scenariouncertainty, :pgrowth)
    connectparameter(m, :impactedsocioeconomic, :ypcgrowth, :scenariouncertainty, :ypcgrowth)
    connectparameter(m, :impactedsocioeconomic, :changeypcgrowth, :temperaturemacroeconomic, :changeypcgrowth)
    connectparameter(m, :impactedsocioeconomic, :mitigationcost, :emissions, :mitigationcost)

    connectparameter(m, :emissions, :income, :impactedsocioeconomic, :income)
    connectparameter(m, :emissions, :population, :population, :population)
    connectparameter(m, :emissions, :forestemm, :scenariouncertainty, :forestemm)
    connectparameter(m, :emissions, :aeei, :scenariouncertainty, :aeei)
    connectparameter(m, :emissions, :acei, :scenariouncertainty, :acei)
    connectparameter(m, :emissions, :ypcgrowth, :scenariouncertainty, :ypcgrowth)

    connectparameter(m, :climateco2cycle, :mco2, :emissions, :mco2)

    connectparameter(m, :climatech4cycle, :globch4, :emissions, :globch4)

    connectparameter(m, :climaten2ocycle, :globn2o, :emissions, :globn2o)
    connectparameter(m, :climateco2cycle, :temp, :climatedynamics, :temp)

    connectparameter(m, :climatesf6cycle, :globsf6, :emissions, :globsf6)

    connectparameter(m, :climateforcing, :acco2, :climateco2cycle, :acco2)
    connectparameter(m, :climateforcing, :acch4, :climatech4cycle, :acch4)
    connectparameter(m, :climateforcing, :acn2o, :climaten2ocycle, :acn2o)
    connectparameter(m, :climateforcing, :acsf6, :climatesf6cycle, :acsf6)

    connectparameter(m, :climatedynamics, :radforc, :climateforcing, :radforc)

    connectparameter(m, :climateregional, :inputtemp, :climatedynamics, :temp)

    connectparameter(m, :biodiversity, :temp, :climatedynamics, :temp)

    connectparameter(m, :ocean, :temp, :climatedynamics, :temp)

    connectparameter(m, :temperaturemacroeconomic, :regtmp, :climateregional, :regtmp)
    macroparams = readdlm("../nasdata/macroeconomic.csv", ',')
    setparameter(m, :temperaturemacroeconomic, :macrocoeffs, MvNormal(squeeze(macroparams[1,:], 1), macroparams[2:3,:]))

    connectparameter(m, :impactedloss, :income, :socioeconomic, :income)
    connectparameter(m, :impactedloss, :impactedincome, :impactedsocioeconomic, :income)

    return m
end

function getfund(;nsteps=1050, datadir="../data", params=nothing)
    # ---------------------------------------------
    # Load parameters
    # ---------------------------------------------
    if params==nothing
        parameters = loadparameters(datadir)
    else
        parameters = params
    end

    # ---------------------------------------------
    # Construct model
    # ---------------------------------------------
    m = constructfund(nsteps=nsteps)

    # ---------------------------------------------
    # Load remaining parameters from file
    # ---------------------------------------------
    setleftoverparameters(m, parameters)


    setbestguess(m)
    # ---------------------------------------------
    # Return model
    # ---------------------------------------------

    return m
end
