# WIP: conditional interface so non-condition priors can be used in cases where
# we condition observations on prior ones
history(x::Stats) = 0
cond_update!(stats::Stats{T},x::AbstractMatrix{T},w=1.0) where T =
  update!(stats,x[:,1],w)
cond_logpdf_thresh(stats::Stats{T},x::AbstractMatrix{T},thresh) where T =
  logpdf_thresh(stats,x[:,1],thresh)

struct VelocityNormalStats{T,Fn} <: Stats{T}
  μ::Vector{T}
  S::SparseMatrixCSC{T,Int} # correlation of [x_t-w,x_t-w+1,...x_t]
  Sindices::Vector{CartesianIndex{2}}
  w::Int
  n1::Float64
  n2::Float64
  x2_offset::Float64
end
Base.zero(x::VelocityNormalStats) =
  VelocityNormalStats(zero(x.μ),zero(x.S),w,0.0,0.0,0.0,
                      x.distance,x.distance_thresh)
history(x::VelocityNormalStats) = x.w

function cond_update!(stats::VelocityNormalStats{T},x::AbstractMatrix{T},w=1.0) where
  T

  δ = vec(x) .- stats.μ
  stats.n1 += w
  stats.μ .+ δ * (w/stats.n2)

  δ2 = x .- stats.μ
  stats.n2 += w^2
  w2n = (w2/stats.n2)
  for ci in stats.Sindices; stats.S[ci] += δ[ci[1]]*δ[ci[2]] * w2n; end

  stats
end

function mult!(stats::VelocityNormalStats,c)
  stats.n1 *= c
  stats.n2 *= c^2

  stats
end

function Base.:(+)(a::VelocityNormalStats{T},b::VelocityNormalStats{T}) where T
  @assert a.Sindices == b.Sindices "Incompatible velocity structures."
  N1 = a.n1 + b.n1
  aw = a.n1/N1
  bw = b.n2/N1
  μ = a.μ.*aw .+ b.μ.*bw

  N2 = a.n2 + b.n2
  aw = a.n2/N2
  bw = b.n2/N2
  S = a.S.*aw + b.S.*bw
  MultiNormalStats(μ,S,N1,N2,a.x2_offset + b.x2_offset)
end

function cond_logpdf_thresh(stats::VelocityNormalStats{T},
                            x::AbstractMatrix{T},thresh) where T
  ii = find(diagonal(stats.S) ./ stats.n .> thresh)

  d = length(ii)
  Λ = stats.S[ii,ii] ./ stats.n2
  μ = stats.μ[ii]
  n = stats.n1
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n2 + stats.x2_offset - d + 1))

  jj = setdiff(1:length(x),ii)
  logjoint = logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x[ii]) +
    sum(logpdf.(Normal(0,thresh/2),x[jj] .- stats.μ[jj]))

  hi = find(ii .> size(x,1))
  hj = find(jj .> size(x,1))
  loghistory = logpdf_mvt(nu,μ[hi],Λ[hi] .* ((n + 1)/(n * nu)),x[ii][hi]) +
    sum(logpdf.(Normal(0,thresh/2),x[jj][hj] .- stats.μ[hj]))

  logjoint - loghistory
end


