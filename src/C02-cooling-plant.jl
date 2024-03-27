### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 33181610-dcb2-11ee-004b-63c218681aa3
begin
    @info "Loading packages..."

    import PlutoUI
    import Luxor

    using CairoMakie
    using DocStringExtensions
    using Polynomials
    using Psychro
    using Roots
    using Statistics
    using SteamTables
    using Unitful
end

# ╔═╡ b7a018ef-3fb7-42b7-8765-375d9437b0b7
md"""
# Crushing process balance

$(PlutoUI.TableOfContents())
"""

# ╔═╡ f03d3f68-27fa-4255-8cb3-c6ae4ebdb303
md"""
## Study elements

### Reference water cooling

Using current reference data[^1] regarding water cooling circuit flow rate and measured inlet and outlet temperatures, the abstracted heat

```math
\dot{Q} = ρ\dot{V}\Delta{H}=ρ\dot{V}\left\{H(T_{out})-H(T_{in})\right\}
```

Using this simple analysis we identify the cooling capacity of the system when operated under water cooling. In what follows, data from previous acquisitions will be ignored since no water flow rate was actually available and estimations were purely based on *known* production values.
"""

# ╔═╡ c50940b2-bb86-4640-b333-f71874d1ab6f
md"""
### Experimental balances

Performed under air cooling conditions [^2].
"""

# ╔═╡ f5499bbb-9a5b-43ab-adab-d15ebbe7311f
@bind reset1 PlutoUI.Button("Reset values")

# ╔═╡ a065dad4-f32f-4b2c-8be1-dbb9d4c9aabc
@bind reset2 PlutoUI.Button("Reset values")

# ╔═╡ 27ba243d-b2b7-4952-b490-b44a9ef6f1c4
md"""
### Footnotes

[^1]: Reference data from Mar 05th 2024.
[^2]: Reference data from Feb 15th 2024 13h50.
[^T]: Tuning argument.
[^M]: Measured argument.
[^C]: Controls argument.
"""

# ╔═╡ 8682da47-4589-40c6-8bb0-b723abe27bbb
md"""
## Documentation
"""

# ╔═╡ 3d16ba7f-abd5-448a-a740-01c9318c2adb
begin
@info "Constants"

"Ideal gas constant [J/(mol.K)]."
const RGAS::Float64 = 8.314_462_618_153_24

"Reference atmospheric pressure [Pa]."
const PREF::Float64 = 101325.0

"Normal atmospheric temperature [K]."
const TREF::Float64 = 273.15

"Air mean molecular mass [kg/mol]."
const M_AIR::Float64 = 0.0289647

"Normal state concentration [mol/m³]. "
const C_AIR::Float64 = PREF / (RGAS * TREF)

"Coefficients for air enthalpy polynomial [J/kg]."
const H_AIR::Vector{Float64} = [
    -2.6257123774377e+05,
     9.8274248481342e+02,
     4.9125599795629e-02
]
end;

# ╔═╡ 4e0a432c-6de6-4a6d-bf1d-bff569eb0f57
begin
@info "Abstract types"

abstract type AbstractMaterial end
abstract type AbstractLiquidMaterial <: AbstractMaterial end
abstract type AbstractSolidMaterial <: AbstractMaterial end
abstract type AbstractGasMaterial <: AbstractMaterial end

issolid(m::AbstractMaterial)  = m isa AbstractSolidMaterial
isliquid(m::AbstractMaterial) = m isa AbstractLiquidMaterial
isgas(m::AbstractMaterial)    = m isa AbstractGasMaterial
end;

# ╔═╡ 08377a1d-1f7c-4d8f-86aa-d46c926e873d
begin
@info "Concrete types (materials)"

"Liquid water material."
struct Water <: AbstractLiquidMaterial
end

"Solid clinker material."
struct Clinker <: AbstractSolidMaterial
    ρ::Float64
    h::Polynomial

    function Clinker(; ρ = 900.0, h = [0, 850.0])
        return new(ρ, Polynomial(h, :T))
    end
end

"Gas air material."
struct Air <: AbstractGasMaterial
    M̄::Float64
    h::Polynomial

    function Air(; h = H_AIR)
        return new(M_AIR, Polynomial(h, :T))
    end
end
end;

# ╔═╡ d61e2dd1-eea2-45b7-9eca-4374d8e540ce
begin
@info "Library implementation"

##############################
# Concrete types (simulation)
##############################

"Array of materials to include in a stream."
struct StreamPipeline
    materials::Vector{AbstractMaterial}
end

""" Represents a material stream.

Attributes
----------
$(FIELDS)
"""
struct MaterialStream
    "Material mass flow rate [kg/s]."
    ṁ::Float64

    "Stream temperature [K]."
    T::Float64

    "Stream pressure [Pa]."
    P::Float64

    "Components mass fractions [-]."
    Y::Vector{Float64}

    "Materials pipeline associated to `Y`."
    pipeline::StreamPipeline
end

""" Represents an energy stream.

Attributes
----------
$(FIELDS)
"""
struct EnergyStream
    "Energy flow provided by stream [W]."
    ḣ::Float64
end

""" Represents a pipeline with heat transfer.

Models
------
1. `:TARGET_EXIT_TEMP` evaluates the heat transfer lost to environment \
  provided a target final stream temperature given by keyword argument \
  `temp_out`. Product temperature is updated through an `EnergyStream` \
  built with energy exchange computed through `exchanged_heat`, so that \
  numerical value can be slightly different from target value.
1. `:USING_GLOBAL_HTC` makes use of a global heat transfer coefficient \
  to evaluate heat flux across the pipe.

To-do's
-------
- Implement heat transfer losses through a convective heat transfer
  coefficient (HTC) computed from a suitable Nusselt number, for use
  of pipeline in *simulation* mode.

Attributes
----------
$(FIELDS)
"""
struct TransportPipeline
    "The output material stream at the end of pipeline."
    product::MaterialStream

    "The heat exchanged in pipeline [W]."
    power::EnergyStream

    function TransportPipeline(;
            product,
            model,
            verbose = true,
            kwargs...
        )
        ##########
        # INITIAL
        ##########

        Δq = 0.0
        power = EnergyStream(Δq)

        ##########
        # MODEL
        ##########

        if model == :TARGET_EXIT_TEMP
            # Compute enthalpy change with environment.
            Δq = exchanged_heat(product, kwargs[:temp_out])

            # Stream of energy to correct system temperature.
            power = EnergyStream(Δq)

            # Correct energy in both streams.
            product += power
        end

        if model == :USING_GLOBAL_HTC
            # Compute enthalpy change with environment.
            T∞ = kwargs[:temp_env]
            T₂ = kwargs[:temp_out]
            U = kwargs[:glob_htc]
            Δq = U * (T∞ - 0.5 * (T₂ + product.T))

            # Stream of energy to correct system temperature.
            power = EnergyStream(Δq)

            # Correct energy in both streams.
            product += power
        end

        ##########
        # POST
        ##########

        verbose && begin
            rounder(v) = round(v; digits = 1)
            p = rounder(ustrip(uconvert(u"kW", Δq * u"W")))
            T = rounder(ustrip(uconvert(u"°C", product.T * u"K")))

            @info """
            TransportPipeline with model $(model)

            Heat loss to environment..........: $(p) kW
            Product stream final temperature..: $(T) °C
            """
        end

        ##########
        # NEW
        ##########

        return new(product, power)
    end
end

""" Represents a solids separator with efficiency η.

To-do's
-------
- Add inverse model to automatically tune efficiency η.

Attributes
----------
$(FIELDS)
"""
struct SolidsSeparator
    "Solids separation efficiency [-]."
    η::Float64

    "The stream to be separated into solids and others."
    source::MaterialStream

    "The output solids stream."
    solids::MaterialStream

    "The output remaining stream."
    others::MaterialStream

    function SolidsSeparator(source; η = 1.0)
        # Retrieve elements kept constant.
        T, P, pipe = source.T, source.P, source.pipeline

        # Retrieve array of solids mass fractions, zeroeing other materials.
        Ys0 = map((m, Y)->issolid(m) ? Y : 0, pipe.materials, source.Y)

        # The mass flow of split solids stream corresponds to the mass flow
        # of solids multiplied by the separation efficiency of separator.
        m_sol = source.ṁ * sum(Ys0) * η

        # The mass flow of secondary stream (recirculating solids and other
        # phases) is the total nass flow minus separated solids.
        m_oth = source.ṁ - m_sol

        # Mass flow rates of each species in original flow.
        mk0 = source.ṁ * source.Y

        # For solids we multiply their total flow rate by the renormalized
        # mass fractions of individual species, balance for other stream.
        mk1 = m_sol * Ys0 ./ sum(Ys0)
        mk2 = mk0 .- mk1

        # Mass fractions are recomputed for each stream.
        Ys1 = mk1 / m_sol
        Ys2 = mk2 / m_oth

        # Create new streamsp
        solids = MaterialStream(m_sol, T, P, Ys1, pipe)
        others = MaterialStream(m_oth, T, P, Ys2, pipe)

        return new(η, source, solids, others)
    end
end

""" Represents a crushing device with cooling system.

Models
------
1. `:TARGET_COOLANT_TEMP` evaluates the heat transfer lost to coolant \
  provided a target final stream temperature given by keyword argument \
  `temp_out`. Product temperature is updated through an `EnergyStream` \
  built with energy exchange computed through `exchanged_heat`, so that \
  numerical value can be slightly different from target value.
1. `:USING_GLOBAL_HTC` makes use of a global heat transfer coefficient \
  to evaluate heat flux across the cooling stream.

Attributes
----------
$(FIELDS)
"""
struct CooledCrushingMill
    "The input meal applied to crushing process."
    rawmeal::MaterialStream

    "The output material stream at the end of product pipeline."
    product::MaterialStream

    "The output material stream at the end of cooling pipeline."
    coolant::MaterialStream

    "The power applied to the crushing process [W]"
    power::EnergyStream

    "The heat exchanged in between product and cooling pipelines [W]."
    loss::EnergyStream

    "Global heat transfer coefficient [W/K]."
    globalhtc::Union{Nothing, Float64}

    function CooledCrushingMill(;
            product,
            coolant,
            power,
            model,
            verbose = true,
            kwargs...
        )
        ##########
        # INITIAL
        ##########

        Δq = 0.0
        loss = EnergyStream(Δq)
        meal = product
        ghtc = nothing
        
        ##########
        # MODEL
        ##########

        if model == :TARGET_COOLANT_TEMP
            # Compute enthalpy change in cooling stream.
            Δq = exchanged_heat(coolant, kwargs[:temp_out])

            # Stream of energy to correct system temperature.
            loss = EnergyStream(Δq)

            # Correct energy in both streams.
            product += power - loss
            coolant += loss
        end

        if model == :USING_GLOBAL_HTC
            # Compute enthalpy change with coolant.
            T∞ = 0.5*(kwargs[:temp_out] + coolant.T)
            T₂ = 0.5*(kwargs[:temp_cru] + product.T)
            ghtc = kwargs[:glob_htc]

            # TODO add losses to environment through shell like in:
            # 10.45*(T_ext - 0.5*(kwargs[:temp_env] + coolant.T))
            
            # The value must be the absolute intake by the coolant
            # thus a minus sign in front of it.
            Δq = -ghtc * (T∞ - T₂)

            # Stream of energy to correct system temperature.
            loss = EnergyStream(Δq)

            # Correct energy in both streams.
            product += power - loss
            coolant += loss
        end

        ##########
        # POST
        ##########

        verbose && begin
            rounder(v) = round(v; digits = 1)
            Q = rounder(ustrip(uconvert(u"kW", power.ḣ * u"W")))
            p = rounder(ustrip(uconvert(u"kW", Δq * u"W")))
            T = rounder(ustrip(uconvert(u"°C", product.T * u"K")))

            @info """
            CooledCrushingMill with model $(model)

            Power applied to product stream...: $(Q) kW
            Heat extracted by cooling system..: $(p) kW
            Product stream final temperature..: $(T) °C
            """
        end

        ##########
        # NEW
        ##########

        return new(meal, product, coolant, power, loss, ghtc)
    end
end

##############################
# Model options managers
##############################

"Manage use of `TransportPipeline` with different models."
function transport_pipe(product, temp_out, temp_env, glob_htc)
    verbose = false

    function target_pipeline()
        return TransportPipeline(; model = :TARGET_EXIT_TEMP, verbose,
            product, temp_out)
    end

    function simul_pipeline()
        temp_out = isnothing(temp_out) ? product.T : temp_out

        pipe = TransportPipeline(; model = :USING_GLOBAL_HTC, verbose,
            product, temp_out, glob_htc, temp_env)

        temp_out = pipe.product.T
        return pipe, temp_out
    end

    pipe = isnothing(glob_htc) ? target_pipeline() : let
        pipe, temp_out = simul_pipeline()
        pipe
    end

    return pipe, temp_out
end

"Manage use of `CooledCrushingMill` with different models."
function cooled_crushing(; product, coolant, power, temp_out, temp_cru,
                           glob_htc, α = 1.0e-04)
    verbose = false

    function target_pipeline()
        return CooledCrushingMill(; model = :TARGET_COOLANT_TEMP, verbose,
            product, coolant, power, temp_out)
    end

    function simul_pipeline()
        temp_out = isnothing(temp_out) ? coolant.T : temp_out
        temp_cru = isnothing(temp_cru) ? product.T : temp_cru
        
        pipe = CooledCrushingMill(; model = :USING_GLOBAL_HTC, verbose, 
            product, coolant, power, temp_out, glob_htc, temp_cru)

        temp_out = α * pipe.coolant.T + (1-α) * temp_out
        temp_cru = α * pipe.product.T + (1-α) * temp_cru
        
        # temp_out = pipe.coolant.T
        # temp_cru = pipe.product.T
        
        return pipe, temp_out, temp_cru
    end

    pipe = isnothing(glob_htc) ? target_pipeline() : let
        pipe, temp_out, temp_cru = simul_pipeline()
        pipe
    end

    return pipe, temp_out, temp_cru
end

##############################
# System model
##############################

""" Process `CooledCrusherModel` inputs in SI units.

Parameters are the same as the attributes. The next are optional key-words.
Operating pressure `P_env` defaults to atmospheric pressure. If `T_in_cool`
is not provided, it is set to `T_env`. To compute cooling efficiency it is
recommended to provide the measured `T_out_cool` to get a precise balance.
Other values `T_in_sep` and `T_out_rec` are used for better balance.

Attributes
----------
$(FIELDS)
"""
struct CooledCrusherInputs
    "Environment temperature."
    T_env::Unitful.Quantity{Float64}

    "Operating pressure."
    P_env::Unitful.Quantity{Float64}

    "Efficiency of solids separator."
    ηseparator::Float64

    "Power applied to crusher."
    power_crusher::Unitful.Quantity{Float64}

    "Mass flow rate of cooling fluid."
    ṁ_cooler::Unitful.Quantity{Float64}

    "Mass flow rate of clinker."
    ṁ_clinker::Unitful.Quantity{Float64}

    "Mass flow rate of crushing auxiliary air."
    ṁ_cru_air::Unitful.Quantity{Float64}

    "Mass flow rate of separator auxiliary air."
    ṁ_sep_air::Unitful.Quantity{Float64}

    "Mass flow rate of parasite air in system."
    ṁ_par_air::Unitful.Quantity{Float64}

    "Mass flow rate of air at system exit."
    ṁ_tot_air::Unitful.Quantity{Float64}

    "Inlet temperature of cooling fluid."
    T_in_cool::Unitful.Quantity{Float64}

    "Outlet temperature of cooling fluid (forced mode)."
    T_out_cool::Union{Nothing, Unitful.Quantity{Float64}}

    "Inlet temperature of product in separator (forced mode)."
    T_in_sep::Union{Nothing, Unitful.Quantity{Float64}}

    "Outlet temperature of product in recirculation (forced mode)."
    T_out_rec::Union{Nothing, Unitful.Quantity{Float64}}

    "Crushing outlet temperature for simulation initial guess."
    T_out_crush::Union{Nothing, Unitful.Quantity{Float64}}

    function CooledCrusherInputs(; T_env, ηseparator, power_crusher, ṁ_cooler,
            ṁ_clinker, ṁ_cru_air, ṁ_sep_air, ṁ_par_air, ṁ_tot_air, kwargs...)

        kw = Dict(kwargs)

        kelvin_or_na(v) = isnothing(v) ? v : uconvert(u"K", v)

        return new(
            uconvert(u"K", T_env),
            uconvert(u"Pa", get(kw, :P_env, 1.0u"atm")),
            0.01ηseparator,
            uconvert(u"W", power_crusher),
            uconvert(u"kg/s", ṁ_cooler),
            uconvert(u"kg/s", ṁ_clinker),
            uconvert(u"kg/s", ṁ_cru_air),
            uconvert(u"kg/s", ṁ_sep_air),
            uconvert(u"kg/s", ṁ_par_air),
            uconvert(u"kg/s", ṁ_tot_air),
            uconvert(u"K", get(kw, :T_in_cool, T_env)),
            kelvin_or_na(get(kw, :T_out_cool, nothing)),
            kelvin_or_na(get(kw, :T_in_sep, nothing)),
            kelvin_or_na(get(kw, :T_out_rec, nothing)),
            kelvin_or_na(get(kw, :T_out_crush, nothing))
        )
    end
end

""" Solver parameters for `CooledCrusherModel`.

All attributes are optional and can be provided through key-words.

Attributes
----------
$(FIELDS)
"""
struct CooledCrusherSolverPars
    "Maximum number of solver iterations. Default is `100`."
    max_iter::Int64

    "Absolute tolerance for temperature. Default is `1.0e-10`."
    T_tol::Float64

    "Absolute tolerance for mass flow. Default is `1.0e-10`."
    ṁ_tol::Float64

    "Verbosity control flag. Default is `true`."
    verbose::Bool

    function CooledCrusherSolverPars(; kwargs...)
        kw = Dict(kwargs)
        return new(
            get(kw, :max_iter, 100),
            get(kw, :T_tol,    1.0e-10),
            get(kw, :ṁ_tol,    1.0e-10),
            get(kw, :verbose,  true)
        )
    end
end

""" Creates unit operations for `CooledCrusherModel`.

Specific key-word arguments of this structure include the materials specification
with the following *pure fluid* defaults. If the user needs specific material
properties, the material must implement all the method interfaces provided for the
default materials below.

- `cooler  = (Air(),     [1.0])`
- `clinker = (Clinker(), [1.0, 0.0])`
- `air     = (Air(),     [0.0, 1.0])`

Attributes
----------
$(FIELDS)
"""
struct CooledCrusherUnits
    "Cooling fluid materials pipeline."
    pipe_cool::StreamPipeline

    "Production fluid materials pipeline."
    pipe_prod::StreamPipeline

    "Energy stream applied to crushing mill."
    milling_power::EnergyStream

    "Crusher cooling stream."
    cooling_stream::MaterialStream

    "Clinker main feed stream."
    clinker_stream::MaterialStream

    "Auxiliary crusher air stream."
    crusher_air_stream::MaterialStream

    "Auxiliary separator air stream."
    separator_air_stream::MaterialStream

    "Parasite air stream mixed with clinker."
    parasite_air_stream::MaterialStream

    "Product recirculation stream."
    recirc_stream::MaterialStream

    "Balance air as complement to parasite air."
    balance_air_stream::MaterialStream

    "Premix of clinker and parasite air."
    meal_stream::MaterialStream

    "Crusher operation model."
    crusher::CooledCrushingMill

    "Separator operation model."
    separator::SolidsSeparator

    "Cyclone operation model."
    cyclone::SolidsSeparator

    "Pipeline between crusher and separator."
    transport_sep::TransportPipeline

    "Product recirculation pipeline."
    transport_rec::TransportPipeline

    function CooledCrusherUnits(inputs, solver; kwargs...)
        kw = Dict(kwargs)

        # Strip physical units.
        T_env         = ustrip(inputs.T_env)
        P_env         = ustrip(inputs.P_env)
        T_in_cool     = ustrip(inputs.T_in_cool)
        power_crusher = ustrip(inputs.power_crusher)
        ṁ_cooler      = ustrip(inputs.ṁ_cooler)
        ṁ_clinker     = ustrip(inputs.ṁ_clinker)
        ṁ_cru_air     = ustrip(inputs.ṁ_cru_air)
        ṁ_sep_air     = ustrip(inputs.ṁ_sep_air)
        ṁ_par_air     = ustrip(inputs.ṁ_par_air)
        ṁ_tot_air     = ustrip(inputs.ṁ_tot_air)

        val_or_na(qty) = isnothing(qty) ? qty : ustrip(qty)

        T_out_cool    = val_or_na(inputs.T_out_cool)
        T_in_sep      = val_or_na(inputs.T_in_sep)
        T_out_rec     = val_or_na(inputs.T_out_rec)
        T_out_crush   = val_or_na(inputs.T_out_crush)
            
        # Get material streams.
        cooler, Y_cool     = get(kw, :cooler,  (Air(),     [1.0]))
        clinker, Y_clinker = get(kw, :clinker, (Clinker(), [1.0, 0.0]))
        air, Y_air         = get(kw, :air,     (Air(),     [0.0, 1.0]))

        pipe_cool = StreamPipeline([cooler])
        pipe_prod = StreamPipeline([clinker, air])

        # Aliases
        T, P = T_env, P_env

        # Applied power at mill.
        milling_power = EnergyStream(power_crusher)

        # Cooling material stream.
        cooling_stream = MaterialStream(
            ṁ_cooler, T_in_cool, P_env, Y_cool, pipe_cool)

        # Clinker material stream.
        clinker_stream = MaterialStream(
            ṁ_clinker, T_env, P_env, Y_clinker, pipe_prod)

        # Milling air stream.
        crusher_air_stream = MaterialStream(
            ṁ_cru_air, T_env, P_env, Y_air, pipe_prod)

        # Separator air stream.
        separator_air_stream = MaterialStream(
            ṁ_sep_air, T_env, P_env, Y_air, pipe_prod)

        # Air leaks in mill.
        parasite_air_stream = MaterialStream(
            ṁ_par_air, T_env, P_env, Y_air, pipe_prod)

        # Dummy recirculation for initialization.
        recirc_stream = MaterialStream(
            0.0, T_env, P_env, Y_clinker, pipe_prod)

        # XXX Experimental code.
        htc_sep = get(kw, :htc_pipe_sep, nothing)
        htc_rec = get(kw, :htc_pipe_rec, nothing)
        htc_cru = get(kw, :htc_cooling,  nothing)
        # temp_cru = T_out_crush
        
        # Premix meal that is not iterated upon.
        meal_stream = clinker_stream
        meal_stream += crusher_air_stream
        meal_stream += parasite_air_stream

        # Run iterative procedure till steady-state.
        itercount = 0
        ṁnow = 1.0e+09
        Tnow = 1.0e+09

        crusher = nothing
        separator = nothing
        transport_sep = nothing
        transport_rec = nothing

        # Tentative
        if !isnothing(htc_cru)
            _, T_out_cool, T_out_crush = cooled_crushing(
                product  = meal_stream + recirc_stream,
                coolant  = cooling_stream,
                power    = milling_power,
                temp_out = T_out_cool,
                temp_cru = T_out_crush, 
                glob_htc = nothing
            )
        end
        
        while itercount <= solver.max_iter
            # Mix meal and recirculation
            product = meal_stream + recirc_stream
            
            # Add crushing energy and cool down system.
            # TODO T_out_crush can be actually *computed*!
            crusher, T_out_cool, T_out_crush = cooled_crushing(
                product  = product,
                coolant  = cooling_stream,
                power    = milling_power,
                temp_out = T_out_cool,
                temp_cru = T_out_crush, 
                glob_htc = htc_cru
            )

            # Loose some heat in vertical pipeline.
            rets = transport_pipe(crusher.product, T_in_sep, T_env, htc_sep)
            transport_sep, T_in_sep = rets

            # Add air stream to crushed product.
            toseparator2 = transport_sep + separator_air_stream

            # Recover streams from separator
            separator = SolidsSeparator(toseparator2; η = inputs.ηseparator)

            # Loose some heat in recirculation pipeline.
            rets = transport_pipe(separator.solids, T_out_rec, T_env, htc_rec)
            transport_rec, T_out_rec = rets

            # Compute property changes.
            ΔT = abs(Tnow - transport_rec.product.T)
            Δṁ = abs(ṁnow - transport_rec.product.ṁ)

            # Check convergence.
            if ΔT < solver.T_tol && Δṁ < solver.ṁ_tol
                solver.verbose && begin
                    @info "CooledCrusherModel: $(itercount) iterations"
                end
                break
            end

            # Prepare next iteration.
            recirc_stream = transport_rec.product
            Tnow = recirc_stream.T
            ṁnow = recirc_stream.ṁ
            itercount += 1
        end

        cyclone = SolidsSeparator(separator.others; η = 1.0)

        # Total air leaving the system (after cyclone)
        balance_air_stream = MaterialStream(
            ṁ_tot_air - cyclone.others.ṁ, T, P, Y_air, pipe_prod)

        return new(
            pipe_cool,
            pipe_prod,
            milling_power,
            cooling_stream,
            clinker_stream,
            crusher_air_stream,
            separator_air_stream,
            parasite_air_stream,
            recirc_stream,
            balance_air_stream,
            meal_stream,
            crusher,
            separator,
            cyclone,
            transport_sep,
            transport_rec
        )
    end
end

""" Implements an air cooled crusher circuit.

Constructor parameters are described in `CooledCrusherInputs` together with
additional key-word arguments, which are extended in `CooledCrusherSolverPars`
and `CooledCrusherUnits`.

Attributes
----------
$(FIELDS)
"""
struct CooledCrusherModel
    "Processed model inputs with units."
    inputs::CooledCrusherInputs

    "Solver operations to manage operation units."
    solver::CooledCrusherSolverPars

    "*De facto* process integration model."
    unitops::CooledCrusherUnits

    function CooledCrusherModel(; T_env, ṁ_cooler, ṁ_clinker, ṁ_cru_air,
        ṁ_sep_air, ṁ_par_air, ṁ_tot_air, power_crusher, ηseparator, kwargs...)

        inputs = CooledCrusherInputs(; T_env, ηseparator, power_crusher, ṁ_cooler,
            ṁ_clinker, ṁ_cru_air, ṁ_sep_air, ṁ_par_air, ṁ_tot_air, kwargs...)

        solver = CooledCrusherSolverPars(; kwargs...)

        unitops = CooledCrusherUnits(inputs, solver; kwargs...)

        return new(inputs, solver, unitops)
    end
end

##############################
# density()
##############################

# density(mat::AbstractMaterial, T, P) = error("Not implemented")

density(mat::Clinker, T, P) = mat.ρ

density(mat::Water, T, P) = 1.0 / SpecificV(P, T)

density(mat::Air, T, P) = (P * mat.M̄) / (RGAS * T)

@doc "Evaluates the density of material [kg/m³]." density

##############################
# enthalpy()
##############################

# TODO consider using clamp()
    
# enthalpy(mat::AbstractMaterial, pars...) = error("Not implemented")

enthalpy(mat::Clinker, T, P) = mat.h(T)

enthalpy(mat::Water, T, P) = 4182.0T

enthalpy(mat::Air, T, P) = mat.h(T)

function enthalpy(pipe::StreamPipeline, T, P, Y)
    return sum(Y .* enthalpy.(pipe.materials, T, P))
end

function enthalpy(stream::MaterialStream, T, P)
    return enthalpy(stream.pipeline, T, P, stream.Y)
end

function enthalpy(stream::MaterialStream; kwargs...)
    kw = Dict(kwargs)
    T = get(kw, :T, stream.T)
    P = get(kw, :P, stream.P)
    Y = get(kw, :Y, stream.Y)
    return enthalpy(stream.pipeline, T, P, Y)
end

@doc "Evaluates the enthalpy of material [J/kg]." enthalpy

##############################
# enthalpyflowrate()
##############################

enthalpyflowrate(s::MaterialStream) = s.ṁ * enthalpy(s)

enthalpyflowrate(e::EnergyStream) = e.ḣ

@doc "Enthalpy flow rate of given stream [W]." enthalpyflowrate

##############################
# General
##############################

"Heat exchanged with stream to match outlet temperature."
function exchanged_heat(s::MaterialStream, T_out)
    # The rate of heat leaving the system [kg/s * J/kg = W].
    ḣ_out = s.ṁ * enthalpy(s; T = T_out)

    # The change of rate across the system [W].
    return ḣ_out - enthalpyflowrate(s)
end

"Compute normal flow rate based on measurements"
function normal_flow_rate(; T, ⌀, U)
    return 3600 * (TREF / T) * U * π * ⌀^2 / 4
end

"Convert [Nm³/h] to [kg/h]."
nm3h_to_kg_h(q) = C_AIR * M_AIR  * q

"Convert [kg/h] to [Nm³/h]."
kg_h_to_nm3h(m) = m / (C_AIR * M_AIR)

end;

# ╔═╡ 6bddbc2a-3775-4c45-a303-0c6f430d0092
begin
@info "Overloads from Base"

function Base.:+(s₁::MaterialStream, s₂::MaterialStream;
                 verbose = false, message = "")
    # Can only mix streams using same material pipeline.
    @assert s₁.pipeline == s₂.pipeline

    # Retrieve flow rates.
    ṁ₁ = s₁.ṁ
    ṁ₂ = s₂.ṁ

    # Total mass flow is the sum of stream flows.
    ṁ = ṁ₁ + ṁ₂

    # Mass weighted average pressure.
    P = (ṁ₁ * s₁.P + ṁ₂ * s₂.P) / ṁ

    # Compute species total mass flow and divide by total flow
    # rate to get resulting stream mass fractions.
    Y = (ṁ₁ * s₁.Y + ṁ₂ * s₂.Y) / ṁ

    # Energy flow is the sum of individual stream flows.
    ĥ = enthalpyflowrate(s₁) + enthalpyflowrate(s₂)

    function f(t)
        # Create a stream with other conditions fixed.
        sn = MaterialStream(ṁ, t, P, Y, s₁.pipeline)

        # Check if with temperature `t` it matches `ĥ`.
        return enthalpyflowrate(sn) - ĥ
    end

    # Find new temperature starting from the top.
    T = find_zero(f, max(s₁.T, s₂.T))

    # Create resulting stream.
    sₒ = MaterialStream(ṁ, T, P, Y, s₁.pipeline)

    verbose && begin
        rounder(v) = round(v; digits = 1)
        T1 = rounder(ustrip(uconvert(u"°C", s₁.T * u"K")))
        T2 = rounder(ustrip(uconvert(u"°C", s₂.T * u"K")))
        To = rounder(ustrip(uconvert(u"°C", sₒ.T * u"K")))

        # TODO make more informative!
        @info """
        MaterialStream addition (+) $(message)

        First stream temperature..........: $(T1) °C
        Second stream temperature.........: $(T2) °C
        Resulting stream temperature......: $(To) °C
        """
    end

    return sₒ
end

function Base.:+(s::MaterialStream, e::EnergyStream)
    # Energy flow is the sum of individual stream flows.
    ĥ = enthalpyflowrate(s) + enthalpyflowrate(e)

    function f(t)
        # Create a stream with other conditions fixed.
        sn = MaterialStream(s.ṁ, t, s.P, s.Y, s.pipeline)

        # Check if with temperature `t` it matches `ĥ`.
        return enthalpyflowrate(sn) - ĥ
    end

    # # Find new temperature starting from current temperature.
    T = find_zero(f, TREF)

    # Create resulting stream.
    return MaterialStream(s.ṁ, T, s.P, s.Y, s.pipeline)
end

function Base.:/(s::MaterialStream, n::Number)
     MaterialStream(s.ṁ / n, s.T, s.P, s.Y, s.pipeline)
end
    
function Base.:-(e₁::EnergyStream, e₂::EnergyStream)
    return EnergyStream(e₁.ḣ - e₂.ḣ)
end

function Base.:+(e₁::EnergyStream, e₂::EnergyStream)
    return EnergyStream(e₁.ḣ + e₂.ḣ)
end

function Base.:+(p::TransportPipeline, s::MaterialStream)
    return p.product + s
end

end

# ╔═╡ 5569ce89-9a36-4dd1-acf8-f7a1a4fbf556
begin
    struct Solid <: AbstractSolidMaterial
        ρ::Float64
        h::Polynomial

        Solid() = new(900, Polynomial([0, 1000], :T))
    end

    enthalpy(s::Solid, T, P) = s.h(T)
    density(s::Solid, T, P) = s.ρ

    solid = Solid()
    enthalpy(solid, 1000, 1)
    
    prod = StreamPipeline([Solid()])
    cool = StreamPipeline([Air()])

    meal  = MaterialStream(1.0, 300, PREF, [1.0], prod)
    wind  = MaterialStream(4.0, 300, PREF, [1.0], cool)
    recir = MaterialStream(0.0, 300, PREF, [1.0], prod)
    power = EnergyStream(500_000)
    crusher = nothing

    max_iter = 10
    T_cru = zeros(max_iter)
    T_rec = zeros(max_iter)
    T_coo = zeros(max_iter)
    loss = zeros(max_iter)
    
    for iter in 1:max_iter
        global recir, crusher
        
        the_meal = meal + recir

        crusher = CooledCrushingMill(;
                verbose  = false,
                product  = the_meal,
                coolant  = wind,
                power    = power,
                model    = :TARGET_COOLANT_TEMP,
                temp_out = 350
        )

        recir = crusher.product / 2
        
        T_cru[iter] = the_meal.T
        T_rec[iter] = recir.T
        T_coo[iter] = crusher.coolant.T
        loss[iter] = crusher.loss.ḣ / 1000
    end

    with_theme(theme_light()) do
        f = Figure()
        ax = Axis(f[1, 1])
        lines!(ax, T_cru; color = :green, label = "Meal")
        lines!(ax, T_rec; color = :red,   label = "Recycle")
        lines!(ax, T_coo; color = :blue,  label = "Coolant")
        lines!(ax, loss; color = :orange, label = "Loss")
        axislegend(ax)
        f
    end
end

# ╔═╡ 4f04674f-ec72-41ef-8395-8fa7bd249a94
let
    M = [
        21.2   1.2   13.0  56.0;
        18.8   1.1   10.0  55.0;
        19.6   1.4   17.0  56.0;
        18.6   1.5   10.0  55.0;
        19.0   1.2   11.0  53.0;
        20.0   1.4   19.0  55.0;
        20.0   1.4   18.0  55.0;
        20.0   1.4   17.0  54.0;
    ]

    h(T) = SpecificH(PREF * u"Pa", T)
    ρ(T) = density(Water(), T, PREF * u"Pa")

    # h(T) =  4182.0u"J/(kg*K)" * T
    # ρ(T) = 996.0u"kg/m^3"

    V̇ = uconvert.(u"m^3/s", M[:, 2] * u"m^3/hr")
    T1 = uconvert.(u"K", M[:, 1] * u"°C")
    T2 = uconvert.(u"K", M[:, 4] * u"°C")

    ρw = map(ρ, T1)
    ΔH = map(h, T2) - map(h, T1)

    Q̇ = mean(uconvert.(u"kW", @. ρw * V̇ * ΔH))
    V̇ = mean(M[:, 2])

    @info """
    Reference cooling power.... $(round(ustrip(Q̇); digits = 1)) kW
    Mean water flow rate....... $(round(V̇; digits = 1)) m³/h
    """
end

# ╔═╡ 8e57eefe-2115-40d9-bb6e-eadadbdbe762
begin
    @info "Reference measured data"

    const ϕ             = 31.5
    const ηseparator    = 47.65

    const T_env         = 5.0
    const T_out_cool    = 93.0
    const T_in_sep      = 73.0
    const T_out_rec     = 39.0

    const q̇_tot_air     = 3600.0
    const q̇_cru_air     = 1881
    const q̇_sep_air     = 431.0

    const ṁ_clinker     = 820.0
    const power_crusher = 107.0

    # Flow rate at cooling system outlet
    const Q̇cool = normal_flow_rate(; T = T_out_cool + TREF, ⌀ = 0.05, U = 13.6)

    # Separator air inlet.
    const Q̇seps = normal_flow_rate(; T = T_env + TREF, ⌀ = 0.16, U = 10.0)

    # Leak flow rates at balls loading and clinker charger.
    const Q̇leak = let
        Q̇7 = normal_flow_rate(; T = T_env + TREF, ⌀ = 0.20, U = 1.6)
        Q̇9 = normal_flow_rate(; T = T_env + TREF, ⌀ = 0.18, U = 2.5)
        Q̇7 + Q̇9
    end
end;

# ╔═╡ 3aec0b9b-5aa2-4f47-baf3-ac0593c5fe6d
begin
    @info "Postprocessing implementation"

    "Default slider with displayed value."
    slider(rng, def) = PlutoUI.Slider(rng, default=def, show_value=true)

    "Standardized report for `CooledCrusherModel`."
    function report(model::CooledCrusherModel; show_tree = true)
        fmt(x) = round(x; digits = 1)

        cooling = model.unitops.cooling_stream
        crusher = model.unitops.crusher
        separator = model.unitops.separator
        trans_rec = model.unitops.transport_rec
        # T_out_rec = model.inputs.T_out_rec

        @info """
        CooledCrusherModel

        Cooling circuit:
        - inlet temperature...... $(fmt(cooling.T - TREF)) °C
        - outlet temperature..... $(fmt(crusher.coolant.T - TREF)) °C
        - energy intake.......... $(fmt(crusher.loss.ḣ / 1000)) kW

        Crusher circuit:
        - inlet temperature...... $(fmt(crusher.rawmeal.T - TREF)) °C
        - outlet temperature..... $(fmt(crusher.product.T - TREF)) °C

        Recirculation circuit:
        - flow rate.............. $(fmt(3600separator.solids.ṁ)) kg/h
        - initial temperature.... $(fmt(separator.solids.T - TREF)) °C
        - final temperature...... $(fmt(trans_rec.product.T - TREF)) °C

        """

        if show_tree
            @info model.inputs
            @info model.unitops
        end
    end

    "Graphical display of crusher balance results."
    function get_results_diagram(model; kwargs...)
        function inlet_main(p0, p1, p2, c)
            Luxor.sethue(c)
            Luxor.arrow(p0, p1)
            Luxor.move(p0)
            Luxor.line(p1)
            Luxor.line(p2)
            Luxor.strokepath()
        end

        kwargs_dict = Dict(kwargs)

        # For @svg
        height = get(kwargs_dict, :height, 300)
        width  = get(kwargs_dict, :width, 700)
        saveas = get(kwargs_dict, :saveas, "crusher.svg")

        # Display control
        showcrusher   = get(kwargs_dict, :showcrusher, true)
        showseparator = get(kwargs_dict, :showseparator, true)

        Luxor.@svg let
            colorbkg = "#EEEEEE"
            colorair = "#0055FF"
            colorsol = "#00AA44"
            colormix = "#FF552F"
            colorrfr = "#0099FF"

            Luxor.background(colorbkg)

            let # Leak air (7+9).
                p0 = Luxor.Point(-265, -50)
                p1 = Luxor.Point(-250, -50)
                p2 = Luxor.Point(-200, 0)
                inlet_main(p0, p1, p2, colorair)
            end

            let # Clinker inlet.
                p0 = Luxor.Point(-265, 50)
                p1 = Luxor.Point(-250, 50)
                p2 = Luxor.Point(-200, 0)
                inlet_main(p0, p1, p2, colorsol)
            end

            let # Crushing pipeline.
                p0 = Luxor.Point(-200, 0)
                p1 = Luxor.Point(125, 0)
                p2 = Luxor.Point(125, -100)
                p3 = Luxor.Point(200, -100)
                Luxor.sethue(colormix)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.line(p2)
                Luxor.line(p3)
                Luxor.strokepath()
            end

            let # Crushing air inlet.
                p0 = Luxor.Point(-150, -50)
                p1 = Luxor.Point(-150, 0)
                Luxor.sethue(colorair)
                Luxor.arrow(p0, p1)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.strokepath()
            end

            showcrusher && let # Crusher.
                Luxor.move(-50, 30)
                Luxor.line(Luxor.Point(100, 30))
                Luxor.line(Luxor.Point(100, -30))
                Luxor.line(Luxor.Point(-50, -30))
                Luxor.closepath()
                Luxor.sethue("orange"); Luxor.fillpreserve()
                Luxor.sethue("black"); Luxor.strokepath()
            end

            let # Cooling system.
                arrowheadlength = 10

                p0 = Luxor.Point(75, 60)
                p1 = Luxor.Point(75, 0)
                p2 = Luxor.Point(-25, 0)
                p3 = Luxor.Point(-25, 60)

                pm = Luxor.Point(75, 40)
                pn = Luxor.Point(-25, 40 + arrowheadlength)

                Luxor.sethue(colorrfr)
                Luxor.arrow(p0, pm; arrowheadlength)
                Luxor.arrow(p2, pn; arrowheadlength)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.line(p2)
                Luxor.line(p3)
                Luxor.strokepath()
            end

            let # Separator air.
                p0 = Luxor.Point(100, -70)
                p1 = Luxor.Point(125, -70)
                Luxor.sethue(colorair)
                Luxor.arrow(p0, p1)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.strokepath()
            end

            let # Recirculation pipe.
                p0 = Luxor.Point(125, -100)
                p1 = Luxor.Point(-100, -100)
                p2 = Luxor.Point(-100, 0)
                Luxor.sethue(colorsol)
                Luxor.arrow(p1, p2)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.line(p2)
                Luxor.strokepath()
            end

            showseparator && let # Separator.
                Luxor.move(125, -80)
                Luxor.line(Luxor.Point(110, -114))
                Luxor.line(Luxor.Point(140, -114))
                Luxor.closepath()
                Luxor.sethue("gray"); Luxor.fillpreserve()
                Luxor.sethue("black"); Luxor.strokepath()
            end

            let # Packing products.
                p0 = Luxor.Point(200, -100)
                p1 = Luxor.Point(200, 0)
                Luxor.sethue(colorsol)
                Luxor.arrow(p0, p1)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.strokepath()

                p0 = Luxor.Point(200, -100)
                p1 = Luxor.Point(250, -100)
                Luxor.sethue(colorair)
                Luxor.arrow(p0, p1)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.strokepath()


                p0 = Luxor.Point(230, 0)
                p1 = Luxor.Point(230, -100)
                Luxor.sethue(colorair)
                Luxor.arrow(p0, p1)
                Luxor.move(p0)
                Luxor.line(p1)
                Luxor.strokepath()
            end

            showseparator && let # Packing.
                Luxor.move(200, -80)
                Luxor.line(Luxor.Point(185, -114))
                Luxor.line(Luxor.Point(215, -114))
                Luxor.closepath()
                Luxor.sethue("gray"); Luxor.fillpreserve()
                Luxor.sethue("black"); Luxor.strokepath()
            end

            let # Joining points.
                radius = 2
                Luxor.sethue("black")
                Luxor.circle(Luxor.Point(-200, 0), radius; action = :fill)
                Luxor.circle(Luxor.Point(-150, 0), radius; action = :fill)
                Luxor.circle(Luxor.Point(-100, 0), radius; action = :fill)
                Luxor.circle(Luxor.Point(125, -70), radius; action = :fill)
                Luxor.circle(Luxor.Point(230, -100), radius; action = :fill)
            end

            let # Add annotations
                inp = model.inputs
                ops = model.unitops

                crusher   = ops.crusher
                separator = ops.separator
                cyclone   = ops.cyclone

                trans_sep = ops.transport_sep
                trans_rec = ops.transport_rec

                liq_cool = ops.pipe_cool.materials[1] isa AbstractLiquidMaterial

                # halign = :center
                valign = :middle

                round1(x) = round(x; digits = 1)

                celsius(T) = round1(ustrip(T) - TREF)

                crush_power = round1(crusher.power.ḣ / 1000)
                cooling_power = round1(crusher.loss.ḣ / 1000)

                ṁ_clinker = round1(3600ops.clinker_stream.ṁ)
                ṁ_recircs = round1(3600separator.solids.ṁ)
                ṁ_product = round1(3600cyclone.solids.ṁ)

                q_par_air = round1(kg_h_to_nm3h(3600ops.parasite_air_stream.ṁ))
                q_cru_air = round1(kg_h_to_nm3h(3600ops.crusher_air_stream.ṁ))
                q_sep_air = round1(kg_h_to_nm3h(3600ops.separator_air_stream.ṁ))
                q_tot_air = round1(kg_h_to_nm3h(3600ustrip(inp.ṁ_tot_air)))
                q_oth_air = round1(kg_h_to_nm3h(3600ops.balance_air_stream.ṁ))

                if liq_cool
                    q_cooling = "$(round1(3600ops.cooling_stream.ṁ)) kg/h (0)"
                else
                    q_cooling = "$(round1(kg_h_to_nm3h(3600ops.cooling_stream.ṁ))) Nm³/h (0)"
                end

                T_env     = celsius(inp.T_env)
                T_crush1  = celsius(crusher.rawmeal.T)
                T_crush2  = celsius(crusher.product.T)
                T_coolant = celsius(crusher.coolant.T)
                T_recircs = celsius(separator.solids.T)
                T_in_sep  = celsius(trans_sep.product.T)
                T_out_rec = celsius(trans_rec.product.T)

                let # Controls
                    Luxor.sethue("black")

                    Luxor.text("Environment at $(T_env) °C",
                         Luxor.Point(-340, -140); valign, halign = :left)

                    Luxor.text("Crushing @ $(crush_power) kW",
                         Luxor.Point(25, -20); valign, halign = :center)

                    Luxor.text("Cooling @ $(-1cooling_power) kW",
                         Luxor.Point(25, 10); valign, halign = :center)

                    Luxor.text(q_cooling,
                         Luxor.Point(75, 65);  valign, halign = :left, angle = π/2)

                    Luxor.text("$(q_tot_air) Nm³/h (5)",
                         Luxor.Point(255, -100); valign, halign = :left)

                    Luxor.text("$(q_sep_air) Nm³/h (4)",
                         Luxor.Point(96,  -70);  valign, halign = :right)

                    Luxor.text("$(q_cru_air) Nm³/h (3)",
                         Luxor.Point(-150, -60);  valign, halign = :center)

                    Luxor.text("$(q_par_air) Nm³/h (2)",
                         Luxor.Point(-269, -50); valign, halign = :right)

                    Luxor.text("$(ṁ_clinker) kg/h (1)",
                         Luxor.Point(-269, 50); valign, halign = :right)

                    Luxor.text("$(q_oth_air) Nm³/h (6)",
                         Luxor.Point(240, -70); valign, halign = :left, angle = π/2)

                    Luxor.text("$(ṁ_product) kg/h",
                         Luxor.Point(210, -70); valign, halign = :left, angle = π/2)
                end

                let # Measurements
                    Luxor.sethue("#FF2299")

                    Luxor.text("$(T_coolant) °C",
                         Luxor.Point(-25, 65); valign, halign = :left, angle = π/2)

                    Luxor.text("$(T_out_rec) °C",
                         Luxor.Point(-90, -40); valign, halign = :left, angle = π/2)

                    Luxor.text("$(T_in_sep) °C",
                         Luxor.Point(135, -69); valign, halign = :left, angle = π/2)
                end

                let # Fitted
                    Luxor.sethue("#9922FF")

                    Luxor.text("$(ṁ_recircs) kg/h @ $(T_recircs) °C",
                         Luxor.Point(25, -110); valign, halign = :center)
                end

                let # Main result
                    Luxor.sethue("#FF0000")
                    Luxor.fontsize(16)

                    Luxor.text("$(T_crush1) °C",
                         Luxor.Point(-55, 15); valign, halign = :right)

                    Luxor.text("$(T_crush2) °C",
                         Luxor.Point(105, 15); valign, halign = :left)
                end
            end

        end width height saveas
    end
end;

# ╔═╡ 922da43b-46a6-4335-957d-9d2637d35679
refmodel, reffig = let
    # Contribution of leaks *through* crusher.
    Q̇avai = q̇_tot_air - q̇_cru_air - q̇_sep_air
    ϕwarn = 100Q̇leak / Q̇avai

    # Select how to compute leak.
    q̇_par_air = (q̇_tot_air - q̇_cru_air - q̇_sep_air) * ϕ / 100

    model = CooledCrusherModel(;
        ηseparator    = ηseparator,
        T_env         = T_env * u"°C",
        ṁ_cooler      = nm3h_to_kg_h(Q̇cool) * u"kg/hr",
        ṁ_clinker     = ṁ_clinker * u"kg/hr",
        ṁ_cru_air     = nm3h_to_kg_h(q̇_cru_air) * u"kg/hr",
        ṁ_sep_air     = nm3h_to_kg_h(q̇_sep_air) * u"kg/hr",
        ṁ_par_air     = nm3h_to_kg_h(q̇_par_air) * u"kg/hr",
        ṁ_tot_air     = nm3h_to_kg_h(q̇_tot_air) * u"kg/hr",
        power_crusher = power_crusher * u"kW",
        T_in_cool     = T_env * u"°C",
        T_out_cool    = T_out_cool * u"°C",
        T_in_sep      = T_in_sep * u"°C",
        T_out_rec     = T_out_rec * u"°C",
        verbose       = false
    )

    report(model; show_tree = false)

    fig = get_results_diagram(model)

    model, fig
end;

# ╔═╡ 13c18387-5dc4-4d60-97bb-04377e724431
reffig

# ╔═╡ 1f6c5710-9b8c-43b9-bba4-8248aedf63f8
refmodel

# ╔═╡ 13f9fbe4-4f5c-4059-819e-28f8ae023c0a
let
    function get_pipe_global_htc(pipe, T∞, T₁)
        return 2pipe.power.ḣ / (2T∞ - pipe.product.T - T₁)
    end

    model = refmodel
    ops = model.unitops
    T∞ = ustrip(model.inputs.T_env)

    U_sep = let
        pipe = ops.transport_sep
        T₁ = ops.crusher.product.T
        get_pipe_global_htc(pipe, T∞, T₁)
    end

    U_rec = let
        pipe = ops.transport_rec
        T₁ = ops.separator.solids.T
        get_pipe_global_htc(pipe, T∞, T₁)
    end

    U_cru =  let
        pipe = ops.crusher
        T₁ = pipe.rawmeal.T
        T∞ = 0.5*(ops.cooling_stream.T + pipe.coolant.T)
        # @info 2T∞
        # @info pipe.rawmeal.T + pipe.product.T
        # T∞ = pipe.coolant.T
        -get_pipe_global_htc(pipe, T∞, T₁)
    end
    
    @info """
    Coefficient estimation for model use in simulation mode:

    - From crusher to separator.... $(round(U_sep; digits = 1)) W/K
    - Along recirculation stream... $(round(U_rec; digits = 1)) W/K
    - Cooling global coefficient... $(round(U_cru; digits = 1)) W/K
    """
end

# ╔═╡ 2f0f105d-ba0c-467f-87f7-664316976be5
let
reset1

md"""
### Air cooling simulator

| Quantity              | Value                                          | Unit |
|----------------------:|:-----------------------------------------------|:----:|
Use leak % below [^T]   | $(@bind useϕ1 PlutoUI.CheckBox(default=false)) |
Leak percentage [^T]    | $(@bind ϕleaks1    slider(0:0.5:100, 31.5))    | [%]
Separator eff. [^T]     | $(@bind ηsep1      slider(45:0.05:65, 47.65))  | [%]
Environment temp. [^M]  | $(@bind T_env1     slider(-2.0:0.5:45, 5.0))   | [°C]
Cooler feed rate [^C]   | $(@bind q̇_cooler1  slider(50:5:1000, 75))      | [Nm³/h]
Milling power [^C]      | $(@bind power1     slider(90:1:120, 107))      | [kW]
Clinker feed rate [^C]  | $(@bind ṁ_clinker1 slider(450:10:900, 820))    | [kg/h]
Crusher air flow [^C]   | $(@bind q̇_cru_air1 slider(1600:10:2500, 1881)) | [Nm³/h]
Separator air flow [^C] | $(@bind q̇_sep_air1 slider(300:10:800, 431))    | [Nm³/h]
Total air flow [^C]     | $(@bind q̇_tot_air1 slider(2500:50:4000, 3600)) | [Nm³/h]
"""
end

# ╔═╡ c5c0fcac-c4ab-4ec0-baf6-4ec57f38833c
let
    # Contribution of leaks *through* crusher.
    Q̇avai = q̇_tot_air1 - q̇_cru_air1 - q̇_sep_air1
    ϕwarn = 100Q̇leak / Q̇avai

    # Select how to compute leak.
    ϕ = useϕ1 ? ϕleaks1 : ϕwarn
    q̇_par_air1 = Q̇avai * ϕ/100

    @warn md"""
    | Quantity      | Measured              | Using |
    |:--------------|:---------------------:|:-----:|
    | Recommended ϕ | $(round(ϕwarn, digits=1))% | $(round(ϕ, digits=1))%
    | Parasite air  | $(round(Q̇leak)) Nm³/h | $(round(q̇_par_air1))
    | Separator air | $(round(Q̇seps)) Nm³/h | $(round(q̇_sep_air1)) Nm³/h
    """

    model = CooledCrusherModel(;
        ηseparator    = ηsep1,
        T_env         = T_env1 * u"°C",

        ṁ_clinker     = ṁ_clinker1 * u"kg/hr",
        ṁ_cooler      = nm3h_to_kg_h(q̇_cooler1) * u"kg/hr",
        ṁ_cru_air     = nm3h_to_kg_h(q̇_cru_air1) * u"kg/hr",
        ṁ_sep_air     = nm3h_to_kg_h(q̇_sep_air1) * u"kg/hr",
        ṁ_par_air     = nm3h_to_kg_h(q̇_par_air1) * u"kg/hr",
        ṁ_tot_air     = nm3h_to_kg_h(q̇_tot_air1) * u"kg/hr",
        power_crusher = power1 * u"kW",

        T_in_cool     = T_env1 * u"°C",
        T_out_cool    = 0 * u"°C",
        T_out_crush   = 200.0 * u"°C",
            
        verbose       = true,
        max_iter      = 10000,
    
        htc_pipe_sep  = 362.1,
        htc_pipe_rec  = 97.9,
        # htc_cooling   = 23679.7,
        htc_cooling   = 20000.0,

    )

    report(model; show_tree = true)

    get_results_diagram(model)
end

# ╔═╡ 23eb1f32-8d85-4f2c-a2fb-ebed13797703
let
reset2

md"""
### Water cooling simulator

| Quantity              | Value                                          | Unit |
|----------------------:|:-----------------------------------------------|:----:|
Use leak % below [^T]   | $(@bind useϕ2 PlutoUI.CheckBox(default=false)) |
Leak percentage [^T]    | $(@bind ϕleaks2    slider(0.0:0.5:100.0, 31.5))| [%]
Separator eff. [^T]     | $(@bind ηsep2      slider(45:0.05:65, 50.8))   | [%]
Environment temp. [^M]  | $(@bind T_env2     slider(-2.0:0.5:45.0, 5.0)) | [°C]
Cooler feed rate [^C]   | $(@bind ṁ_cooler2  slider(0.5:0.1:2.0, 1.3))   | [m³/h]
Cooler start temp. [^M] | $(@bind T_ic2      slider(5:1:40, 20))         | [°C]
Cooler end temp. [^M]   | $(@bind T_oc2      slider(30:1:70, 55))        | [°C]
Milling power [^C]      | $(@bind power2     slider(90:1:120, 105))      | [kW]
Clinker feed rate [^C]  | $(@bind ṁ_clinker2 slider(450:10:900, 550))    | [kg/h]
Crusher air flow [^C]   | $(@bind q̇_cru_air2 slider(1600:10:2500, 2232)) | [Nm³/h]
Separator air flow [^C] | $(@bind q̇_sep_air2 slider(300:10:800, 431))    | [Nm³/h]
Total air flow [^C]     | $(@bind q̇_tot_air2 slider(2500:50:4000, 3500)) | [Nm³/h]
"""
end

# ╔═╡ 0610abd3-0bd4-429d-9544-ba76a9b66dd6
let
    # Contribution of leaks *through* crusher.
    Q̇avai = q̇_tot_air2 - q̇_cru_air2 - q̇_sep_air2
    ϕwarn = 100Q̇leak / Q̇avai

    # Select how to compute leak.
    ϕ = useϕ2 ? ϕleaks2 : ϕwarn
    q̇_par_air2 = Q̇avai * ϕ/100

    model = CooledCrusherModel(;
        ηseparator    = ηsep2,
        T_env         = T_env2 * u"°C",

        ṁ_clinker     = ṁ_clinker2 * u"kg/hr",
        ṁ_cooler      = 1000ṁ_cooler2 * u"kg/hr",
        ṁ_cru_air     = nm3h_to_kg_h(q̇_cru_air2) * u"kg/hr",
        ṁ_sep_air     = nm3h_to_kg_h(q̇_sep_air2) * u"kg/hr",
        ṁ_par_air     = nm3h_to_kg_h(q̇_par_air2) * u"kg/hr",
        ṁ_tot_air     = nm3h_to_kg_h(q̇_tot_air2) * u"kg/hr",
        power_crusher = power2 * u"kW",

        T_in_cool     = T_ic2 * u"°C",
        T_out_cool    = T_oc2 * u"°C",

        verbose       = true,
        cooler        = (Water() , [1.0]),
        htc_pipe_sep  = 362.1,
        htc_pipe_rec  = 97.9,
    )

    report(model; show_tree = true)

    get_results_diagram(model)
end

# ╔═╡ 6980e81b-b5e0-4179-9cba-9a35aa1113bb
md"""
## Verifications

### Water cooling power

!!! warning "Do the same as for air!"

    Implement own properties (i.e. polynomials) in water to be self-contained.

    ```julia
    # NOTE: inputs in kJ/kg = f(MPa, K)
    enthalpy(mat::Water, T, P) = 1.0e+03SpecificH(1.0e-06P, T)
    ```
"""

# ╔═╡ 6bc03d9b-c794-434e-8453-b01316c503af
let
    water = Water()

    V̇ = 2.0u"m^3/hr"

    T = [20u"°C"; 40u"°C"]

    ρ = density(water, T[1], PREF * u"Pa")

    H = map(x->SpecificH(PREF * u"Pa", x), T)

    Q̇ = uconvert(u"kW", ρ * V̇ * diff(H)[1])

    @info "Reference cooling power of $(Q̇)"
end

# ╔═╡ 088a0e4d-ba74-423d-8747-5d36b1689ed1
md"""
### Air properties

!!! warning "IMPORTANT"

     The polynomial fit in the following cell was manually copied and stored as a constant. If any updates are provided to the air enthalpy computation method, this needs to be restored. This choice was made to make it easier to export developed code to an external module.

"""

# ╔═╡ 41752cc8-0409-43cc-9325-4b8a4d75eaea
let
    "Create air enthalpy function for given pressure and dew point."
    function get_air_enthalpy_function(; P = 1.0u"atm", D = 0.0u"°C")
        return x->Psychro.enthalpy(Psychro.MoistAir, x, Psychro.DewPoint, D, P)
    end

    # XXX: the behaviour in `LinRange` is different from what happens after
    # one calls `ustrip` for a single value. Here outputs are already in K.
    T = LinRange(0.0u"°C", 200.0u"°C", 1000+1)
    H = get_air_enthalpy_function().(T)

    cut = 100
    hf = fit(ustrip(T), ustrip(H), 2; var = :T)
    Hf = hf.(ustrip(T)[begin:cut:end])

    f = Figure(size = (700, 400))
    ax = Axis(f[1, 1]; title = string(hf))
    lines!(ax, ustrip(T), ustrip(H) ./ 1000; color = :black)
    scatter!(ax, ustrip(T)[begin:cut:end], Hf ./ 1000; color = :red)
    ax.xlabel = "Temperature [K]"
    ax.ylabel = "Enthalpy [kJ/kg]"
    f
end

# ╔═╡ b23d604a-0f34-42df-91c7-19a525343646
md"""
### Additional tests
"""

# ╔═╡ bf8e1963-5852-4ca1-bfb1-79000f6e2ea0
let
    @info "Tests implementation and run"

    "Mass flow rate per species in stream."
    mass_flow_spec(s) = s.ṁ * s.Y

    "Verify specific heat approximation from enthalpy."
    function test_specific_heat(; T = [20u"°C"; 40u"°C"], atol = 0.005)
        H = map(x->SpecificH(101325.0u"Pa", x), T)
        c = mean(x->SpecificCP(101325.0u"Pa", x), T)
        return c - (H[2] - H[1]) / (T[2] - T[1]) < atol * unit(c)
    end

    "Test enthalpy evaluations."
    function test_enthalpy(; materials = [Water(), Clinker(), Air()])
        function the_test(m, T, P)
            pipe = StreamPipeline([m])
            stream = MaterialStream(1.0, T, P, [1.0], pipe)

            h1 = enthalpy(m, T, P)
            h2 = enthalpy(pipe, T, P, [1.0])
            h3 = enthalpy(stream)

            result = h1 ≈ h2 ≈ h3
            !result && @warn "Failed testing $(typeof(m))"

            return result
        end

        return all(the_test.(materials, 293.15, 101325.0))
    end

    "Test mass balance in solids separator."
    function test_solids_separator()
        pipe = StreamPipeline([Clinker(), Air()])
        source = MaterialStream(1.0, TREF, PREF, [0.2, 0.8], pipe)
        splits = SolidsSeparator(source; η = 0.5)

        m0 = mass_flow_spec(source)
        m1 = mass_flow_spec(splits.solids)
        m2 = mass_flow_spec(splits.others)

        return sum(m1 + m2 - m0) ≈ 0.0
    end

    @assert test_specific_heat()
    @assert test_enthalpy()
    @assert test_solids_separator()
end;

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
Psychro = "9516f557-4a54-5a79-b954-c272e753c77a"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
SteamTables = "43dc94dd-f011-5c5d-8ab2-5073432dc0ba"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
CairoMakie = "~0.11.9"
DocStringExtensions = "~0.9.3"
Luxor = "~3.8.0"
PlutoUI = "~0.7.58"
Polynomials = "~4.0.6"
Psychro = "~0.3.0"
Roots = "~2.1.2"
SteamTables = "~1.4.1"
Unitful = "~1.19.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "9dc21d347d4e0dafa61db1b16e09c7ab6e3076e1"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "222ee9e50b98f51b5d78feb93dd928880df35f06"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.3.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "c0d491ef0b135fd7d63cbc6404286bc633329425"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.36"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "e2a9873379849ce2ac9f9fa34b0e37bde5d5fe0a"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c5aeb516a84459e0318a02507d2261edad97eb75"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "588e0d680ad1d7201d4c6a804dcb1cd9cba79fbb"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "6dc1bbdd6a133adf4aa751d12dbc2c6ae59f873d"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.11.9"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "9b1ca1aa6ce3f71b3d1840c538a8210a043625eb"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "0f4b5d62a88d8f59003e43c25a8a90de9eb76317"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.18"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "d4e9dc4c6106b8d44e40cd4faf8261a678552c7c"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.8.12"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "c5c28c245101bd59154f649e19b038d15901b5dc"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.2"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "73d1214fec245096717847c62d389a5d2ac86504"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.22.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Format]]
git-tree-sha1 = "f3cf88025f6d03c194d73f5d13fee9004a108329"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.6"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "055626e1a35f6771fe99060e835b72ca61a52621"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "d4f85701f569584f2cff7ba67a137d03f0cfb7d0"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.3"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "5694b56ccf9d15addedc35e9a4ba9c317721b788"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.10"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "af13a277efd8a6e716d79ef635d5342ccb75be61"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.10.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "RoundingEmulator"]
git-tree-sha1 = "2d6d22fe481eff6e337808cc0880c567d7324f9a"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.8"
weakdeps = ["DiffRules", "ForwardDiff", "RecipesBase"]

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "fee018a29b60733876eb557804b5b109dd3dd8a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.8"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "6355fb9a4d22d867318db186fd09b09b35bd2ed7"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.6.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5edc369a598dfde567269dc6add5812cfa10cd5"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Primes", "SimplePolynomials"]
git-tree-sha1 = "d76cec8007ec123c2b681269d40f94b053473fcf"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.2.7"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "Juno", "LaTeXStrings", "PrecompileTools", "Random", "Requires", "Rsvg"]
git-tree-sha1 = "aa3eb624552373a6204c19b00e95ce62ea932d32"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "3.8.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "27af6be179c711fb916a597b6644fbb5b80becc0"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.20.8"

[[deps.MakieCore]]
deps = ["Observables", "REPL"]
git-tree-sha1 = "248b7a4be0f92b497f7a331aed02c1e9a878f46b"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.7.3"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "679c1aec6934d322783bd15db4d18f898653be4f"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.27.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "96ca8a313eb6437db5ffe946c457a401bbb8ce1d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mods]]
git-tree-sha1 = "924f962b524a71eef7a21dae1e6853817f9b658f"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "2.2.4"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "302fd161eb1c439e4115b51ae456da4e9984f130"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "MathOptInterface", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d024bfb56144d947d4fafcd9cb5cafbe3410b133"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4745216e94f71cb768d58330b059c9b76f32cb66"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.14+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "eb3f9df2457819bf0a9019bd93cc451697a0751e"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.20"

[[deps.PikaParser]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "d6ff87de27ff3082131f31a714d25ab6d0a88abf"
uuid = "3bbf5609-3e7b-44cd-8549-7c69f321e792"
version = "0.6.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a9c7a523d5ed375be3983db190f6a5874ae9286d"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.6"
weakdeps = ["ChainRulesCore", "FFTW", "MakieCore", "MutableArithmetics"]

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "763a8ceb07833dd51bb9e3bbca372de32c0605ad"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.0"

[[deps.Psychro]]
deps = ["Unitful"]
git-tree-sha1 = "6b988fb38c5f69810bc976a02cf8c259095d915b"
uuid = "9516f557-4a54-5a79-b954-c272e753c77a"
version = "0.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "f39da63aa6d2d88e0c1bd20ed6a3ff9ea7171ada"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.8"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["Accessors", "ChainRulesCore", "CommonSolve", "Printf"]
git-tree-sha1 = "754acd3031a9f2eaf6632ba4850b1c01fe4460c1"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.1.2"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "f65caa24a622f985cc341de81d3f9744435d0d0f"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.6"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "e182b9e5afb194142d4668536345a365ea19363a"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.2"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "7063828369cafa93f3187b3d0159f05582011405"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.17"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["Compat", "PikaParser", "SHA", "Tables", "TupleTools"]
git-tree-sha1 = "10dc702932fe05a0e09b8e5955f00794ea1e8b12"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "1.1.8"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteamTables]]
deps = ["ForwardDiff", "PrecompileTools", "Roots", "Unitful"]
git-tree-sha1 = "a03b9ec4ee2d895b9483544d4c61c21089fd2733"
uuid = "43dc94dd-f011-5c5d-8ab2-5073432dc0ba"
version = "1.4.1"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "34cc045dd0aaa59b8bbe86c644679bc57f1d5bd0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.8"

[[deps.TranscodingStreams]]
git-tree-sha1 = "3caa21522e7efac1ba21834a03734c57b4611c7e"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.4"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "41d61b1c545b06279871ef1a4b5fcb2cac2191cd"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.5.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "07e470dabc5a6a4254ffebc29a1b3fc01464e105"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "37195dcb94a5970397ad425b95a9a26d0befce3a"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.0+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "86e7731be08b12fa5e741f719603ae740e16b666"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.10+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1ea2ebe8ffa31f9c324e8c1d6e86b4165b84a024"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─b7a018ef-3fb7-42b7-8765-375d9437b0b7
# ╟─5569ce89-9a36-4dd1-acf8-f7a1a4fbf556
# ╟─f03d3f68-27fa-4255-8cb3-c6ae4ebdb303
# ╟─4f04674f-ec72-41ef-8395-8fa7bd249a94
# ╟─c50940b2-bb86-4640-b333-f71874d1ab6f
# ╟─8e57eefe-2115-40d9-bb6e-eadadbdbe762
# ╟─13c18387-5dc4-4d60-97bb-04377e724431
# ╟─922da43b-46a6-4335-957d-9d2637d35679
# ╟─1f6c5710-9b8c-43b9-bba4-8248aedf63f8
# ╟─13f9fbe4-4f5c-4059-819e-28f8ae023c0a
# ╟─2f0f105d-ba0c-467f-87f7-664316976be5
# ╟─f5499bbb-9a5b-43ab-adab-d15ebbe7311f
# ╟─c5c0fcac-c4ab-4ec0-baf6-4ec57f38833c
# ╟─23eb1f32-8d85-4f2c-a2fb-ebed13797703
# ╟─a065dad4-f32f-4b2c-8be1-dbb9d4c9aabc
# ╟─0610abd3-0bd4-429d-9544-ba76a9b66dd6
# ╟─27ba243d-b2b7-4952-b490-b44a9ef6f1c4
# ╟─8682da47-4589-40c6-8bb0-b723abe27bbb
# ╟─33181610-dcb2-11ee-004b-63c218681aa3
# ╟─3d16ba7f-abd5-448a-a740-01c9318c2adb
# ╟─4e0a432c-6de6-4a6d-bf1d-bff569eb0f57
# ╟─08377a1d-1f7c-4d8f-86aa-d46c926e873d
# ╟─d61e2dd1-eea2-45b7-9eca-4374d8e540ce
# ╟─6bddbc2a-3775-4c45-a303-0c6f430d0092
# ╟─3aec0b9b-5aa2-4f47-baf3-ac0593c5fe6d
# ╟─6980e81b-b5e0-4179-9cba-9a35aa1113bb
# ╟─6bc03d9b-c794-434e-8453-b01316c503af
# ╟─088a0e4d-ba74-423d-8747-5d36b1689ed1
# ╟─41752cc8-0409-43cc-9325-4b8a4d75eaea
# ╟─b23d604a-0f34-42df-91c7-19a525343646
# ╟─bf8e1963-5852-4ca1-bfb1-79000f6e2ea0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
