module HarmonicPowerModels

    # import pkgs
    import JuMP
    import PowerModels
    import InfrastructureModels
    import Memento

    # pkg constants 
    const _PMs = PowerModels
    const _IMs = InfrastructureModels

    function __init__()
        global _LOGGER = Memento.getlogger(PowerModels)
    end
end
