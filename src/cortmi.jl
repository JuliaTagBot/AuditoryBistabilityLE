export scale_weighting, freq_weighting, weighting

function weighting(x,condition,params)
  condition == :scales ? scale_weighting(x,params[:W_m_σ],params[:W_m_c]) :
  condition == :freqs ? freq_weighting(x,params[:W_m_σ],params[:W_m_c]) :
  condition == :track ? track_weighting(x,params[:W_m_σ_t],params[:W_m_σ_ϕ],
                                        params[:W_m_c]) :
  error("No condition named $condition.")
end

function scale_weight(cort,σ,c)
  s = log.(ustrip.(uconvert.(cycoct,scales(cort))))
  W = @. c*(1 - exp(-(s - s')^2 / (σ*log(2))^2))

  W
end

function scale_weighting(cort,σ,c)
  W = scale_weight(cort,σ,c)
  function helper(x::AbstractArray{T,2}) where T
    reshape(W*vec(sum(x,2)),:,1)
  end

  function helper(x::AbstractArray{T,3}) where T
    reshape(W*vec(sum(x,(1,3))),1,:,1)
  end

  function helper(x::AbstractArray{T,1}) where T
    W*x
  end

  helper
end

function track_weighting(tracks,σ_t,σ_p,c)
  tcs, priors = zip(axisvalues(AxisArrays.axes(tracks,Axis{:params}()))[1]...)
  p = collect(priors)
  t = collect(log.(ustrip.(uconvert.(s,tcs))))
  # @show p
  # @show t
  W = @. (c*(1 - exp(-(p - p')^2 / (σ_p*log(2))^2))*
          (1 - exp(-(t - t')^2 / (σ_t*log(2))^2)))
  # @show size(W)

  x -> W*x
end

function freq_weighting(spect,σ,c)
  f = log.(ustrip.(uconvert.(Hz,freqs(spect))))
  W = @. c*(1 - exp(-(f - f')^2 / (σ*log(2))^2))

  x -> W*x
end
