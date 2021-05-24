#=
- Issues
    Choice of initial parameters
    Choice of price when submitting limit orders
=#

using Random, Distributions, CSV
include("CoinTossXUtilities.jl")
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Data -----#
priceHistory = CSV.File("Data/BRICs-EOD-Compact-2005-2015.xlsx") |> DataFrame
#---------------------------------------------------------------------------------------------------

#----- Structures -----#
abstract type Agent end
mutable struct Chartist <: Agent
    d::Int64 # Time lag
    T::Float64 # Position entry threshold
    τ::Float64 # Position exit threshold
    c::Float64 # Positive constant
    xₜ₋₁::Float64 # Previous position
    xₜ::Float64 # Next position
    tradingGateway::Client # Object to submit orders
end
mutable struct Fundamentalist <: Agent
    vₜ::Float64 # Current perceived value
    T::Float64 # Position entry threshold
    τ::Float64 # Position exit threshold
    c::Float64 # Positive constant
    xₜ₋₁::Float64 # Previous position
    xₜ::Float64 # Next position
    tradingGateway::Client # Object to submit orders
end
struct Parameters # Immutable
    N::Int64 # Number of agents of each type
    a::Float64 # Scale parameter for capital assignment
    dₘᵢₙ::Int64 # Minimum time delay for trend followers
    dₘₐₓ::Int64 # Maximum time delay for trend followers
    vₘᵢₙ::Float64 # Minimum offset for perceived value
    vₘₐₓ::Float64 # Maximum offset for perceived value
    Tₘᵢₙ::Float64 # Minimum threshold for entering positions
    Tₘₐₓ::Float64 # Maximum threshold for entering positions
    τₘᵢₙ::Float64 # Minimum threshold for exiting positions
    τₘₐₓ::Float64 # Maximum threshold for exiting positions
    μₙ::Float64 # Mean of the noise process in vₜ
    σₙ::Float64 # Standard deviation of the noise process in vₜ
end
#---------------------------------------------------------------------------------------------------

#----- State-dependent threshold strategies -----#
function StateDependentThresholdStrategy(agent::Agent, pₜ::Float; pₜ₋::Float)
    mₜ = typeof(agent) == Fundamentalist ? pₜ - agent.vₜ : pₜ₋ - pₜ # Mispricing
    if agent.xₜ₋₁ < 0 # Currently in a short position
        if mₜ > agent.τ # Exit threshold not met
            agent.xₜ = agent.xₜ₋₁ # Keep position and don't trade
        elseif - agent.T < mₜ < agent.τ # Exit threshold met
            agent.xₜ = 0 # Exit short position but don't enter long position
        end
    elseif agent.xₜ₋₁ > 0 # Currently in a long position
        if mₜ < - agent.τ # Exit threshold not met
            agent.xₜ = agent.xₜ₋₁ # Keep position and don't trade
        elseif - agent.τ < mₜ < agent.T # Exit threshold met
            agent.xₜ = 0 # Exit long position but don't enter short position
        end
    else # Currently in neutral position
        if - agent.T < mₜ < agent.T # Entry threshold not met
            agent.xₜ = 0 # Keep position and don't trade
        end # Else statement not necessary since the position is already set by UpdatePosition!
    end # If none of these conditions are met, then the entry thresholds have been met and the agent's position is the one given in UpdatePosition!
    ωₜ = agent.xₜ - agent.xₜ₋₁
    return ωₜ
end
#---------------------------------------------------------------------------------------------------

#----- Initialize agents -----#
function InitializeAgents(priceHistory::Vector{Int64}, parameters::Parameters)
    chartists = Vector{Chartist}(); fundamentalists = Vector{Fundamentalist}()
    for i in 1:parameters.N
        # Chartist
        d = rand(parameters.dₘᵢₙ:parameters.dₘₐₓ)
        T = rand(Uniform(parameters.Tₘᵢₙ, parameters.Tₘₐₓ))
        τ = rand(Uniform(max(-T, parameters.τₘᵢₙ), min(T, parameters.τₘₐₓ))) # Ensure -T < τ < T
        c = parameters.a * (T - τ)
        xₜ = c * sign(priceHistory[end] - priceHistory[end - d]) # History[end] is the price at the previous time step
        push!(chartists, Chartist(d, T, τ, c, 0, xₜ)) # Previous position is set to be neutral
        # Fundamentalist
        vₜ = rand(Unifrom(parameters.vₘᵢₙ, parameters.vₘₐₓ)) # The agent's current value perception
        T = rand(Uniform(parameters.Tₘᵢₙ, parameters.Tₘₐₓ))
        τ = rand(Uniform(max(-T, parameters.τₘᵢₙ), min(T, parameters.τₘₐₓ))) # Ensure -T <= τ < T
        c = parameters.a * (T - τ)
        xₜ = c * sign(vₜ - priceHistory[end]) # History[end] is the price at the previous time step
        push!(fundamentalists, Fundamentalist(vₜ, T, τ, c, 0, xₜ)) # Previous position is set to be neutral
    end
    return chartists, fundamentalists
end
#---------------------------------------------------------------------------------------------------

#----- Update positions -----#
function UpdatePosition!(agent::T, paramaters::Parameters, pₜ::Int64; pₜ₋::Int64) where T <: Agent
    agent.xₜ₋₁ = agent.xₜ # Update previous position
    if isa(agent, Chartist) # Chartist
        agent.xₜ = agent.c * (pₜ - pₜ₋) # Update current position based on the updated price
    else # Fundamentalist
        η = rand(Normal(paramaters.μₙ, paramaters.σₙ))
        agent.vₜ += η # Update value perception to be used in the evaluation of the position to take
        agent.xₜ = agent.c * (agent.vₜ - pₜ) # Update current position. Position is based on the agent's new value perception and the updated price (the agent revises their value perception after the price is updated)
    end
end
#---------------------------------------------------------------------------------------------------

#----- Simulation -----#
function Simulate(priceHistory::Vector{Int64}, horizon::Int64, parameters::Parameters; seed = 1)
    Random.seed!(seed)
    StartJVM()
    client = Login(1, 1)
    try # This ensures that the client gets logged out whether an error occurs or not
        chartists, fundamentalists = InitializeAgents(priceHistory, parameters) # Initialize agents
        Juno.progress() do id # Progress bar
            for t in 1:horizon # Rolling closing auctions
                for (chartist, fundamentalist) in zip(chartists, fundamentalists) # Iterate through both sets of agents
                    # Update positions at the start of each day for orders to be sent. Price history has been updated for the previous orders. Agents base their positions on this updated price
                    UpdatePosition!(chartist, pₜ; pₜ₋ = priceHistory[end - agent.d])
                    UpdatePosition!(fundamentalist, pₜ)
                    # Apply state dependent thresholds to agents' updated positions to determine if they should trade or not based on transaction costs
                    ωₜᶜ = StateDependentThresholdStrategy(chartist, pₜ; pₜ₋ = priceHistory[end - agent.d]) # Chartist - after positions have been revised positions are either entered, exited or kept the same
                    if ωₜᶜ != 0
                        SubmitOrder(chartist.tradingGateway, Order("1", ωₜᶜ > 0 ? "Buy" : "Sell", "Market", abs(round(Int, ωₜᶜ))))
                    end
                    ωₜᶠ = StateDependentThresholdStrategy(agent, pₜ) # Fundamentalist - after positions have been revised positions are either entered, exited or kept the same
                    if ωₜᶠ != 0
                        SubmitOrder(fundamentalist.tradingGateway, Order("1", ωₜᶠ > 0 ? "Buy" : "Sell", "Market", abs(round(Int, ωₜᶠ))))
                    end
                end
                # Receive the closing price
                pₜ = StaticPriceReference(client) # Update price once all orders have been sent
                push!(priceHistory, pₜ) # Update the history of prices
                @info "Trading" progress=(t / horizon) _id=id # Update progress
            end
        end
    finally
        Logout(client)
    end
    return priceHistory
end
#---------------------------------------------------------------------------------------------------
