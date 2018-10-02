module AuditoryBistabilityLE
using ShammaModel
using DataFrames
using Requires
using Statistics
using LinearAlgebra
using SparseArrays
using SpecialFunctions
using TOML

import ShammaModel: Δt, Δf, times, freqs, scales, rates

export adaptmi, drift, scale_weighting, ncomponents, nunits, CoherenceModel,
    fusion_ratio, object_SNR, mask, scene_object_ratio,
    object_SNR2, ab_match, mean_spect, mean_spect2, AdaptMI

using ProgressMeter
next!(x::Progress) = ProgressMeter.next!(x)
next!(x::Nothing) = nothing

include("adaptmi.jl")
include("cortmi.jl")
include("cohere.jl")
include("nmf.jl")
include("simple_tracking.jl")
include("prior_tracking.jl")
include("multi_prior_tracking.jl")
include("stim.jl")
include("peaks.jl")
include("lengths.jl")
include("biapply.jl")
include("bimodel.jl")
include("compress.jl")

# include(joinpath(@__DIR__,"rplots.jl"))
# include(joinpath(@__DIR__,"rplots.jl"))
function __init__()
  @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" include("rplots.jl")
end

end
