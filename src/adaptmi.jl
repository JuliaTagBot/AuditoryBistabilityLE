using Parameters
using Unitful: s, ms, ustrip

# question: should we somehow limit the unit response???
# e.g. with:
sig(x) = 1/(1+exp(-10(x-0.5)))

@with_kw struct AdaptMI{S,I}
  c_x::Float64 = 1.0
  τ_x::typeof(1.0s) = 300ms

  shape_y::S = identity

  c_a::Float64 = 5
  τ_a::typeof(1.0s) = 1.5s

  c_m::Float64 = 10
  τ_m::typeof(1.0s) = 50ms
  W_m::I = inhibit_uniform

  τ_σ::typeof(1.0s) = 100ms
  c_σ::Float64 = 0.3
end

# efficient representation of W*y where W = (ones(n,n) - I)/(n-1) for all n
function inhibit_uniform(y)
  x = (sum(y) .- vec(y)) ./ max(1,length(y)-1)
  reshape(x,size(y)...)
end

########################################
# the generic interface a type must implement for adaptmi to work on it

# empty_timeslice(x) = a single empty time slice of the object x, where x
# consists of multiple time slices
function empty_timeslice
end

# time_indices(x) = the number of time indices for x
function time_indices
end

# after set_timeslice!(y,t,ya_t), time slice t of y === ya_t
function set_timeslice!
end

# get_timeslice(y,t) is the t^th time slice of y
function get_timeslice
end

# approx(f,x,xs...) ≈ f(x,xs...)

# x can be an 'important' value, meaning it is used to determine how the
# remaining values are approximated
function approx
end

# get a time slice in the same format used for values
# within the function called by approx
function approx_empty_timeslice
end

# get a value with a similar format to data, but differing in that it has the
# same format as used within the function called by approx
function approx_similar
end

##############################
# default implementation
# - we assume an array like object and treat the first dimension
#   as the time dimension
# - there is no loss due to approximation

approx(f,args...) = f(args...)
approx_similar(x) = similar(x)
approx_zeros(x,dims...) = fill!(similar(x,dims...),zero(eltype(x)))

########################################
# complex number implementation
# - we 'approximate' values by their magnitude
#   so that m and a remain real, the resulting output
#   can be multiplied by the phase of the original signal
#   to find the final, complex output
function approx(f,x::AbstractArray{<:Complex},args...)
  y = f(abs.(x),args...)
  withangle(angle.(x),y)
end
function withangle(angle,y::Tuple)
  # @show angle
  # @show y[1]
  (y[1].*exp.(angle.*im),y[2:end]...)
end
function withangle(angle,y)
  y.*exp.(angle.*im)
end
approx_similar(x::AbstractArray{<:Complex{T}}) where T =
  similar(x,T)
approx_zeros(y::AbstractArray{<:Complex{T}}) where T =
  zeros(eltype(y),size(y)...)
approx_zeros(y::AbstractArray{<:Complex{T}},dims...) where T =
  zeros(eltype(y),dims...)

################################################################################
# generic adaptation and mutual-inhibition operation
adaptmi(x;progressbar=true,kw...) = adaptmi(x,AdaptMI(;kw...),progressbar)
function adaptmi(x,params::AdaptMI,progressbar=true)
  @assert :time ∈ axisnames(x)
  time = Axis{:time}

  y = similar(x)

  shape_y = params.shape_y
  τ_x = params.τ_x; c_x = params.c_x
  τ_a = params.τ_a; c_a = params.c_a
  τ_m = params.τ_m; c_m = params.c_m; W_m = params.W_m

  a = approx_similar(y) # a = adaptation
  m = approx_similar(y) # m = mutual inhibition
  n = approx_similar(y) # n = normalization
  x_out = approx_similar(y)

  x_t = approx_zeros(y[time(1)])
  n_t = copy(x_t)
  y_t = copy(x_t)
  a_t = copy(x_t)
  m_t = copy(x_t)

  dt_x = eltype(a_t)(Δt(x) / τ_x)
  dt_a = eltype(a_t)(Δt(x) / τ_a)
  dt_m = eltype(a_t)(Δt(x) / τ_m)

  progress = progressbar ? Progress(desc="Adapt/Inhibit: ",ntimes(y)) : nothing
  for ti in Base.axes(times(y),1)
    t = time(ti)

    y[t],m[t],a[t],x_out[t],n[t] = approx(x[t]) do raw_x_t
      @. begin # across all indices...
        # smooth input
        x_t += (c_x*raw_x_t - x_t)*dt_x

        # apply adaptation and inhibition
        y_t = shape_y((1 - c_a*a_t)*x_t - c_m*m_t)

        # update adaptation and inhibition
        a_t += (y_t - a_t)*dt_a
        m_t += ($W_m(y_t) - m_t)*dt_m
      end
      y_t,m_t,a_t,x_t,n_t
    end

    next!(progress)
  end
  y,a,m,x_out,n
end

################################################################################
# generic drifting noise function

drift(x,along_axes...;progressbar=true,kw...) =
  drift(x,AdaptMI(;kw...),along_axes,progressbar)
function drift(x,params::AdaptMI,along_axes=typeof.(AxisArayys.axes(x)),progressbar=true)
  τ_σ, c_σ = params.τ_σ, log(1+params.c_σ)
  time = Axis{:time}

  y = similar(x)
  temp = zero(x[time(1),(ax(1) for ax in along_axes)...])
  σ_t = approx_zeros(temp)
  dims = size(σ_t)
  progress = progressbar ? Progress(desc="Drift: ",ntimes(y)) : nothing
  for t in Base.axes(times(y),1)
    ya_t = x[time(t)]
    ya_t,σ_t = approx(ya_t,σ_t) do ya_t,σ_t
      σ_t .+= -σ_t.*(Δt(x)/τ_σ) .+ randn(dims).*(c_σ*sqrt(2Δt(x)/τ_σ))
      ya_t .*= exp.(σ_t)

      ya_t,σ_t
    end

    y[time(t)] = ya_t

    next!(progress)
  end
  y
end
