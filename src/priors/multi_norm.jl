
mutable struct MultiNormalStats{T} <: Stats{T}
  μ::Vector{T}
  S::Matrix{T}
  n1::Float64
  n2::Float64
  x2_offset::Float64
end
Base.zero(x::MultiNormalStats{T}) where T =
  MultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0)
function MultiNormalStats(data::AbstractMatrix,n=size(data,1),
                          x2_offset=size(data,2))
  μ = mean(data,1)
  S = (data.-μ)'(data.-μ) .* n./size(data,1)
  MultiNormalStats{eltype(data)}(μ,S,n,n,x2_offset)
end
MultiNormalStats(::Type{T},n) where T =
  MultiNormalStats{T}(zeros(T,n),zeros(T,n,n),0,0)

function update!(stats::MultiNormalStats{T},x::AbstractVector{T},w=1.0) where T
  δ = x .- stats.μ
  stats.n1 += w
  stats.μ .+= δ * (w/stats.n2)

  δ2 = x .- stats.μ
  stats.n2 += w^2
  stats.S .+= δ.*δ2' * (w2/stats.n2)

  stats
end

function mult!(stats::MultiNormalStats,c)
  stats.n1 *= c
  stats.n2 *= c^2

  stats
end

function Base.:(+)(a::MultiNormalStats{T},b::MultiNormalStats{T}) where T
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

function logpdf(stats::MultiNormalStats,x::AbstractVector)
  d = length(x)
  Λ = stats.S ./ stats.n2
  μ = stats.μ
  n = stats.n1
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n2 + stats.x2_offset - d + 1))

  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x)
end
pdf(stats::MultiNormalStats,x::AbstractVector) = exp(logpdf(stats,x))

logabsdet_helper(x::AbstractArray,n) = logabsdet(x)[1]
logabsdet_helper(x::Diagonal,n) = logdet(x)
# logabsdet_helper(x::AbstractSparseMatrix,n) = sum(log ∘ abs,diag(lufact(x)[:U]))
logabsdet_helper(x::AbstractSparseMatrix,n) = log(abs(det(x)))
logabsdet_helper(x::Symmetric{<:Any,<:AbstractSparseMatrix},n) = 
  log(abs(det(x)))
logabsdet_helper(x::UniformScaling,n) = n*log(abs(x.λ))
# logabsdet_helper(x::SymmetricToeplitz) = log(abs(prod(diag(chol(x)))))

function logpdf_mvt(v,μ,Σ,x)
  d = length(x)
  C = (lgamma((v+d)/2)) - (lgamma(v/2)+(d/2)log(v*π)) -
    (0.5 * (logabsdet_helper(Σ,d)))
  diff = abs.(μ.-x)
  # sometime Σ is not perfectly positive definite....
  normdiff = max.(0,(diff'/Σ)*diff)
  C + -(v+d)/2*log(1+normdiff/v)
end

function logpdf_thresh(stats::MultiNormalStats,x::AbstractVector,thresh)
  ii = find(diagonal(stats.S) ./ stats.n .> thresh)

  d = length(ii)
  Λ = stats.S[ii,ii] ./ stats.n2
  μ = stats.μ[ii]
  n = stats.n1
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n2 + stats.x2_offset - d + 1))

  jj = setdiff(1:length(x),ii)
  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x[ii]) +
    sum(logpdf.(Normal(0,thresh),x[jj] .- stats.μ[jj]))
end
