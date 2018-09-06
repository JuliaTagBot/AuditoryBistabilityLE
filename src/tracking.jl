using RMath
using Distributions
using StatsFuns
using Combinatorics

@with_kw struct CoherenceTrack
end

const max_sources = 4

mutable struct LogWeighted{T}
  val::T
  log_weight::Float64
end
sample(x::LogWeighted{T},x,args...) where T =
  reweight!.(w -> x.log_weight + w,sample(x.val,x,args...))
log_weight(x::LogWeighted{T}) = x.log_weight
reweight!(x::LogWeighted{T},log_weight::Float64) = x.log_weight = log_weight
reweight!(fn::Function,x) = reweight(x,fn(log_weight(x)))

function normalize_weights(x)
  s = logsumexp(log_weight.(x))
  reweight.(w -> w - s,x)
end

function broadcast(fn::typeof(log_weight),
                   nested::AbstractArray{<:AbstractArray})
  weights = Array{Float64}(sum(length.(x)))
  i = 0
  for xs in nested
    for x in xs
      weights[i+=1] = log_weight(x)
    end
  end
  weights
end
reweight(fn::Function,x::AbstractArray{<:LogWeighted}) = reweight.(fn,x)

resample(N) = xs -> resample(xs,N)
function resample(xs,N)
  weights = rand(Multinomial(N,exp.(log_weight.(xs))))
  filter!(!isinf ∘ log_weight,reweight.(xs,log.(weights)))
end

function resample(nested::AbstractArray{<:AbstractArray},N)
  weights = rand(Multinomial(N,exp.(log_weight.(xs))))
  samples = Array{eltype(eltype(nested))}(N)
  i = 0
  j = 0
  for xs in nested
    for x in xs
      i+=1
      if weights[i] > 0
        samples[j+=1] = reweight(x,log(weights[i]))
      end
    end
  end
  samples
end

# the plan, do the *SIMPLEST* form of the particle filter
# (jus track normals and alpha, nothing more)

mutable struct MultiNormalStats{T}
  x::Vector{T}
  x2::Matrix{T}
  n::Int
  x2_offset::Int
end
function MultiNormalStats(mean,cov,n::Int,x2_offset=length(mean))
  MultiNormalStats(mean*n,cov .+ mean.*mean',n,x2_offset)
end
function MultiNormalStats(data::AbstractArray{T,3} where T,
                          n=size(data,1),x2_offset=size(data,2))
  x = similar(data,size(data,2))
  x2 = similar(data,size(data,2,2))
  for k in size(data,3)
    d = data[:,:,k]
    x .+= vec(mean(d,1))*n
    x2 .+= n .* d'd ./ size(data,1)
  end
  MultiNormalStats(x,x2,n,x2_offset)
end

function MultiNormalStats(data::AbstractArray{T,2} where T,
                          n=size(data,1),x2_offset=size(data,2))
  x = vec(mean(data,1))*n
  x2 = n .* data'data ./ size(data,1)
  MultiNormalStats(x,x2,n,x2_offset)
end

function diff(a::MultiNormalStats{T},b::MultiNormalStats{T}) where T
  total = 0.0
  total += sum(abs,a.x .- b.x)
  total += sum(abs,a.x2 .- b.x2)
  total += sum(abs,a.n - b.n)

  total / (length(a.x)^2 + length(a.x) + 1)
end

function diff(as::AbstractArray{MultiNormalStats{T}},
              bs::AbstractArray{MultiNormalStats{T}})
  ks,hs = length(as) <= length(bs) ? as,bs : bs,as
  max_kn = maximum(k.n for k in ks)

  function diff_helper(ks,hs)
    total = diff.(ks,hs[1:length(k)])
    total += sum(h.n / max_kn for h in hs)

    total / length(h)
  end

  minimum(diff_helper(ks,hs[permute])
          for permute in perumutations(1:length(hs)))
end

function update!(stats::MultiNormalStats{T},x::AbstractVector{T}) where T
  stats.x .+= x
  stats.x2 .+= x.*x'
  stats.n += 1

  stats
end

function logpdf(stats::MultiNormalStats{T},x::AbstractVector{T}) where T
  d = length(x)
  Λ = (stats.x2 .- (stats.x.*stats.x')) ./ stats.n
  μ = stats.x ./ stats.n
  n = stats.n
  nu = stats.n + stats.x2_offset - d + 1

  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu))))
end
pdf(stats::MultiNormalStats{T},x::AbstractVector{T}) where T =
  exp(logpdf(stats,x))

function logpdf_mvt(v,μ,Σ,x)
  d = length(μ)
  C = lgamma((v+1)/2) - (lgamma(v/2)+log(v*π)^(d/2)) - 0.5logabsdet(Σ)
  diff = μ-x

  C*log(1+1/v*(diff'/Σ)*diff)*-(v+d)/2
end

mutable struct BinomialStats
  counts::Int
  n::Int
end
function update!(stats::BinomialStats,x::Bool)
  stats.counts += x
  stats.n += 1

  stats
end
pdf(stats::BinomialStats,x::Bool) =
  x ? stats.counts/stats.n : 1-(stats.counts/stats.n)
logpdf(stats::BinomialStats,x::Bool) = log(pdf(stats,x))

struct Particle{T}
  sources::Vector{MultiNormalStats{T}}
  freqs::Vector{BinomialStats}
  indices::Vector{Int}

  source_prior::MultiNormalStats{T}
  freq_prior::Float64
end
function Particle(source_prior,freq_prior)
  Particle(MultiNormalStats{T}[],BinomialStats[],Int[],source_prior,freq_prior)
end

function update!(p::Particle{T},sources,indices::Vector{Int})
  if indices[end] > length(p.freqs)
    push!(p.freqs,0)
    push!(p.sources,copy(p.source_prior))
  end

  for c in indices
    update!(p.sources[c],sources[c])
  end

  for i in 1:length(p.freqs)
    update!(p.freqs[i],i ∈ indices)
  end
  p.indices = indices

  p
end

function logpdf(p::Particle{T},sources,indices::Vector{Int})
  logsum = 0.0
  freq_total = 0.0

  oldindices = indices
  if indices[end] > length(p.freqs)
    oldindices = indices[1:end-1]
    freq_total += log(p.alpha)
    logsum += logpdf(p.source_prior,sources[:,indices[end]])
  end

  for c in oldindices
    logsum += logpdf(p.sources[c],sources[c])
  end

  for i in 1:length(p.freqs)
    freq_total += pdf(p.freqs[i],i ∈ oldindices)
  end

  logsum += log(freq_total / (length(p.freqs) + p.alpha))

  logsum
end
pdf(p::Particle{T},sources,indices::Vector{Int}) =
  log(pdf(p,sources,indices))

const max_sources = 6

function orderings(n,N)
  @assert n <= N <= max_sources
  [p for subset in combinations(1:N,n) for p in permutations(subset)]
end

function sample(p::Particle,sources::AbstractMatrix)
  orders = orderings(size(sources,2),1:min(length(p.freqs)+1,max_sources))
  samples = Array{LogWeighted{typeof(p)}}(length(orders))

  for (i,(columns,indices)) in enumerate(orders)
    newp = deepcopy(p)
    cols = [sources[:,c] for c in columns]
    update!(newp,cols,indices)
    samples[i] = LogWeighted(newp,logpdf(newp,cols,indices))
  end
  samples
end

function update(particles::AbstractArray{LogWeighted{Particle{T}}},
  sources::AbstractMatrix{T}) where T

  map(p -> sample(p,sources),particles) |>
    normalize_weights |>
    resample(N)
end

function remove_similar!(particles::AbstractArray{LogWeighted{Particle{T}}},
  sources::AbstractMatrix{T}) where T

  i = 0
  for p in particles
    if !any(nearly_the_same(q,p,sources) for q in particles[1:i])
      particles[i+=1] = p
    end
  end
  @views particles[1:i]
end

function nearly_the_same(p::LogWeighted{Particle{T}},
  q::LogWeighted{Particle{T}},sources::AbstractMatrix{T},
  pdftol=0.05,source_diff_tol=1e-4) where T

  if abs(logpdf(p.val,sources,p.indices) -
    logpdf(q.val,soruces,q.indices)) > log(1+pdftol)
    return false
  end

  diff(p.sources,q.sources) > source_diff_tol
end

@with_kw struct ParticleTracking
  N::Int = 1000
  diversify_every::Int = 10
  grouping_prior::Float64 = 1
  data_prior = nothing
end
Tracking(::Val{:particles};params...) = ParticleTracking(;params...)

function track(C::Coherence,params::ParticleTracking)
  time = Axis{:time}
  component = Axis{:component}
  K = ncomponents(C)
  if params.data_prior == nothing
    data_prior = MultiNormalStats(reshape(C,ntimes(C),:,ncomponents(C)),1)
  else
    data_prior = params.data_prior
  end

  particles = fill(LogWeighted(Particle(data_prior,params.grouping_prior),0.0),
                   params.N)
  assignments = Vector{Vector{Int}}(ntimes(C))

  @showprogress "Tracking Sources..." for t in eachindex(times(C))
    particles = update(particles,reshape(C[time(t)],:,ncomponents(C)))
    if diversify_every % t == 0
      particles = diversify!(particles)
    end

    assignments[t] = particles[indmax(log_weight.(particles))].val.indices
  end
  # mmm.... how do we handle labels...., for each
  # component we could ask how likely it is each component
  # is
end

# TODO: how to report the state of the particl filter...

#=
struct OldParticle{T}
  ϵ::Float64

  z::Vector{T}
  ζ::Vector{Float64}
  v::Float64

  x::Matrix{T}
end
function OldParticle(x::TrackParticle{T})
  OldParticle(p.ϵ,p.z,p.ζ,p.v,p.x)
end

struct TrackParticle{T}
  ϵ::Float64

  z::Vector{Bool}
  ζ::Vector{Float64}
  v::Vector{T}

  x::Matrix{T}
  Σ::Array{Cholesky{T}}

  old::OldParticle{T}
end

function TrackParticle(sources::MVector{SufficientStats},
                       counts::MVector{Float64},
                       active::MVector{Bool},
                       order::MVector{Int},
                       old::TrackParticle{T})
  TrackParticle(sources,counts,active,order,old.n + 1.0,
                old.rate_μ,old.rate_σ,old.scale_μ,old.scale_σ)
end

update(stats::SufficientStats,x::Array{T}) where T = update(stas,vec(x))
update(stats::S,x::Vector{T}) where {T,S} =
  error("No update defined for type $S")
update(stats::MvNormalKnownCovStats{T},x::Vector{T}) where T =
  MvNormalKnownCovStats(stats.invΣ,stats.sx + x,stats.tw + 1.0)

function update(stats::MvNormalStats{T},x::Vector{T}) where T
  s = stats.s + x
  tw = stats.tw + 1.0
  m = s/tw

  z = (s-m)
  s2 = z*z'

  MvNormalStats(s,m,s2,tw)
end

possible_active_sets(n) = []

# TODO: generate list of possible compositions of proposals,
# expressly exclude compositions that leave any components
# that are below some threshold.
function sample(s::TrackState{T},p::TrackPraticle{T},
                proposals::Matrix{T}) where T
  samples = TrackParticle{T}[]

  for r in 1:s.particle_N
    K = length(p.z)

    ϵ = if rand(Bool)
      rand(Normal(0.0,s.σ))
    else
      rand(Normal(p.ϵ,s.ϵ_jump))
    end

    v = if rand(Bool)
      rand(InverseGamma(s.V_α,s.V_θ),K)
    else
      # likelihood-ish update
      abs.(rand(Normal(p.v,s.V_jump),K))
    end
    ζ = rand.(Beta.(p.ζ .* v))
    z = [ζ[k] > rand() for k in 1:K]

    for mapping in possible_mappings(proposals,z)
      x = squeeze(sum(proposals[:,mapping],3),3)
      Σ = similar(p.Σ)
      for i in size(x,2)
        Σ[i] = LinAlg.lowrankupdate(p.Σ[i],x[:,i])
      end

      push!(samples,TrackParticle(ϵ,z,ζ,v,x,Σ,p))
    end
  end

  samples
end

function particle_logpdf(s::TrackState{T},p::TrackParticle{T})
  lp = 0.0

  lp += logsumexp(logpdf(Normal(0.0,s.σ),p.ϵ),
                  logpdf(Normal(p.old.ϵ,s.ϵ_jump),p.ϵ))
  lp += logsumexp([logpdf.(InverseGamma(s.V_α,s.V_θ),p.v);
                   logpdf.(Normal(p.old.v,s.V_jump),p.v)])
  lp += sum(logpdf.(Beta.(p.old.ζ .* p.v),p.ζ))
  lp += sum(p.z ? log(p.ζ[k]) : log(1-p.ζ[k]) for k in 1:k)

  lp
end

struct TrackState{T,M}
  method::M
  particles::Array{LogWeighted{TrackParticle{T}}}

  # TODO: fill in parameters
end

function logpdf_active(active::Array{Int},counts::AbstractVector{Float64},
                       n::Float64)
  logsumexp(logpdf(i ∈ active ? log(counts/n) : log((1-counts)/n)
                   for i in 1:length(counts)))
end

function logpdf(s::TrackState{T},p::TrackParticle{T},x::Array{T}) where T
  # TODO: implement

  lp = 0.0
  # the observation is a noisy view of a sum of sources
  lp += logpdf(vec(x),MvNormal(squeeze(sum(p.x[:,p.z],2),2),I*p.ϵ))
  lp += logpdf(p.ϵ,Normal(0.0,s.σ))

  # the sources are some variation of their prior values
  lp += logsumexp(
    logpdf.(MvNormal.(p.old.x[:,p.z],PDMats.(p.Σ[p.z])),p.x[:,p.z]))

  # which sources are active is some variation of their prior values
  lp += sum(p.z[k] ? log(p.ζ[k]) : log(1-p.ζ[k]) for k in eachindex(p.z)))
  lp += sum(logpdf.(p.ζ,Beta.(p.old.ζ .* p.v)))
  lp += logpdf(p.v,InverseGamma(s.V_α,s.V_θ))

  lp
end

# basic steps
# generate proposals by each method

# methods:
# draws from the prior
# transformations of old particles
#    max(0,x-y), max(0,y-x)
# draws from the observed data (i.e. no split)
#

# TODO:

# resample proposals, but without removing
# general classes of methods, i.e. maintain
# a range of filters on scales and rates,
# maintain some particles with multiple sources
# maintain some particles from each method???

# for resampling we need to account for the chance of drawing each type of
# particle, each method should assign a probability to each of its own
# proposals, and to any of the other methods proposals
# and we take a weighted sum these probabilities accorording
# to the weight for each method.


# TODO: we can have multiple, parallel particle filters, each with its own
# method

function sample(state::TrackState,y)
  # generate proposals for the components by various methods
  particles = map(weight(state,y),proposals(state,y))
  TrackState(particles |> normalize_weights |> resample(N),state)
end

weight(state::TrackState,y) = p -> weight(state,p,y)
weight(state::TrackState,particle::TrackParticle,y) =
  LogWeighted(particle,logpdf(state,particle,y) - proposal_logpdf(state,particle))

function sample(states::Array{TrackState{T}},y) where T
  all_particles = Array{Array{LogWeighted{TrackParticle{T}}}}(length(states))
  for (i,state) in enumerate(states)
    all_particles[i] = map(weight(state,y),proposals(state,y))
    map_weights!(all_particles[i],w -> w - log(length(all_particles[i])))
  end

  log_total = map(all_particles) do particles
    logsumexp(log_weight.(particles) .- log(length(particles)))
  end |> logsumexp

  for particles in all_particles
    map_weights!(particles,w -> w - log_total)
  end

  map(enumerate(states)) do (i,state)
    N = state.particles
    TrackState(all_particles[i] |> resample(N),state)
  end
end

struct FactorizingSampler end

struct NMFObservation
  nmf::NMF.Result
  data
end
components(x::NMFObservation) = x.nmf.W
# TODO: normalize components to ensure sum is accuate
data(x::NMFObservation) = mean(x.data,1)

function proposals(state::TrackState{T,FactorizingSampler} where T)
  [proposal for proposal in sample(state,particle,y.data,components(result))
            for particle in state.particles]
end

function (cohere::CoherenceModel{CoherenceTrack})(x)
  # TODO: apply the particle filter above to each observed step of NMF
  # NOTE: this is a
end
