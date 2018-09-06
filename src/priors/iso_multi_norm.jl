export isonorm

mutable struct IsoMultiNormStats{T} <: Stats{T}
  μ::Vector{T}
  S::Vector{T}
  n::Float64
  x2_offset::Float64
end
Statistics.std(x::IsoMultiNormStats) = x.n > 0 ? x.S ./ x.n : Inf
Statistics.zero(x::IsoMultiNormStats{T},d) where T =
  IsoMultiNormStats(zero(x.μ),zero(x.S),0.0,0.0)

function isonorm(prior::Coherence,n::Number,
                 x2_offset::Number=coherence_dim(C);thresh=1e-3)
  error("Outdated prior")
  # data = reshape(mean(prior,4),size(prior,1),:)
  # μ = squeeze(mean(data,1),1)
  # S = squeeze(sum(data.^2,1),1) .* n./size(data,1) .+ n.*thresh.^2
  # IsoMultiNormStats{eltype(data)}(μ,S,n,x2_offset)
end

# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
function update!(stats::IsoMultiNormStats{T},x::AbstractVector{T},
                 w=1.0) where T
  stats.n += w

  μ₀ = copy(stats.μ)
  stats.μ .+= (x .- stats.μ).*(w ./ stats.n)
  stats.S .+= (x .- stats.μ).*(x .- μ₀).*w

  stats
end

function downdate!(stats::IsoMultiNormStats{T},x::AbstractVector{T},
                   w=1.0) where T
  if w ≈ stats.n
    stats.μ .= 0
    stats.S .= 0
    stats.n = 0
  elseif w < stats.n
    μₙ = copy(stats.μ)
    stats.μ .= (stats.μ .- x.*(w ./ stats.n)) ./ (1 .- (w ./ stats.n))
    stats.S .-= (x .- stats.μ).*(x .- μₙ).*w

    stats.n -= w
  else
    error("Sum of weights would become negative. This happens when trying ",
          "to remove more samples then were added using `update!`.")
  end

  stats
end

function mult!(stats::IsoMultiNormStats,c)
  stats.n *= c

  stats
end

function Base.:(+)(a::IsoMultiNormStats{T},b::IsoMultiNormStats{T}) where T
  n = a.n + b.n
  IsoMultiNormStats(a.μ.*(a.n./n) .+ b.μ.*(b.n./n), a.S.+b.S,
                        n,a.x2_offset+b.x2_offset)
end

function logpdf(stats::IsoMultiNormStats,x::AbstractVector)
  d = length(x)
  σ = stats.S ./ stats.n
  nu = stats.n + stats.x2_offset - d + 1
  Λ = Diagonal(σ .* ((stats.n + 1)/(stats.n * nu)))

  if nu < 1
    if length(x) != length(stats.μ)
      error("Dimension mistmatch between data ($(length(x))) and ",
            "the model ($(length(stats.μ))).")
    else
      error("Insufficient data for predictive distribution, ",
            "try collecting more data or increasing x2_offset")
    end
  end

  logpdf_mvt(nu,stats.μ,Λ,x)
end

struct ConstIsoPrior{T} <: Stats{T}
  S::T
  N::Float64
  x2_offset::Float64
end

function isonorm(prior::Number,N::Number,dims,x2_offset::Number=prod(dims))
  ConstIsoPrior{typeof(prior)}(prior*N,N,x2_offset)
end

function Base.zero(prior::ConstIsoPrior,d)
  IsoMultiNormStats(fill(zero(prior.S),d),fill(zero(prior.S),d),0.0,0.0)
end

function Base.:(+)(a::ConstIsoPrior{T},b::IsoMultiNormStats{T}) where T
  d = length(b.μ)
  n = b.n + a.N
  IsoMultiNormStats(b.μ.*(b.n/n),a.S.+b.S,n,a.x2_offset+b.x2_offset)
end

function logpdf(stats::ConstIsoPrior,x::AbstractVector)
  d = length(x)
  σ = stats.S / stats.N
  nu = max(1,floor(Int,stats.N + stats.x2_offset - d + 1))

  logpdf_mvt(nu,0,I*(σ * ((stats.N + 1)/(stats.N * nu))),x)
end
