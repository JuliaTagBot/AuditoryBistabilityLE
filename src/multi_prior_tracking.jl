export map_components

@with_kw struct MultiPriorTracking{C} <: Tracking
  cohere::C
  time_constants::Array{typeof(1.0s)}
  time_constant_bias::Vector{Float64}
  source_priors::AxisArray
  source_prior_bias::Vector{Float64}
  freq_prior
  max_sources::Int = 4
  min_norm::Float64 = Inf
  normalize::Bool = false
end
function Tracking(C,::Val{:multi_prior};time_constants_s=[4],
                  time_constants=time_constants_s*s,
                  time_constant_bias=zeros(length(time_constants)),
                  source_prior_sds=nothing,source_prior_N=nothing,
                  source_prior_bias=zeros(length(source_prior_sds)),
                  freq_ridge=0.0, scale_ridge=0.0,
                  source_priors=nothing,
                  freq_prior=nothing,
                  ridge_threshold=0.05,
                  freq_prior_N = 2, freq_prior_bias = 0,
                  normalize=false,min_norm=Inf,
                  params...)
  checknorm(normalize,min_norm)
  if source_priors == nothing
    @assert(source_prior_sds != nothing,
            "Missing keyword argument `source_prior_sds`.")
    @assert(source_prior_N != nothing,
            "Missing keyword argument `source_prior_N`.")
    source_priors = if iszero(freq_ridge) && iszero(scale_ridge)
      AxisArray([isonorm(sd,source_prior_N, (size(C,2),size(C,3)))
                 for sd in source_prior_sds], Axis{:prior}(source_prior_sds))
    else
      AxisArray([ridgenorm(sd,source_prior_N, (size(C,2),size(C,3)),
                           freq=freq_ridge,scale=scale_ridge,
                           threshold=ridge_threshold)
                 for sd in source_prior_sds], Axis{:prior}(source_prior_sds))
    end
  end

  if freq_prior == nothing
    freq_prior = freqprior(freq_prior_bias,freq_prior_N)
  end
  MultiPriorTracking(;cohere=ShammaModel.Params(C),
                     source_priors=source_priors,
                     source_prior_bias=source_prior_bias,
                     freq_prior=freq_prior,
                     time_constants=time_constants,
                     time_constant_bias=time_constant_bias,
                     normalize=normalize,min_norm=min_norm,
                     params...)
end

function expand_params(params::MultiPriorTracking)
  AxisArray([(PriorTracking(params.cohere,tc,prior,params.freq_prior,
                            params.max_sources,params.min_norm), tb + pb)
             for (tc,tb) in zip(params.time_constants,
                                params.time_constant_bias)
             for (prior,pb) in zip(params.source_priors,
                                   params.source_prior_bias)],
            Axis{:params}([(tc,prior) for tc in params.time_constants
                           for prior in axisvalues(params.source_priors)[1]]))
end

function nitr(C::Coherence,params::MultiPriorTracking)
  ntimes(C) * length(params.time_constants) * length(params.source_priors)
end

function track(C::Coherence,params::MultiPriorTracking,progressbar=true,
               progress=track_progress(progressbar,nitr(C,params),"multi-prior"))

  C_ = prepare_coherence(C,params.min_norm)
  all_params = expand_params(params)
  S = Array{SourceTracking}(undef,size(all_params,1))
  lp = Array{Array{Float64}}(undef,size(all_params,1))

  #=@Threads.threads=# for (i,(p,bias)) in collect(enumerate(all_params))
    S[i], lp[i] = track(C_,p,true,progress)
    lp[i] .+= bias
  end

  (AxisArray(S, AxisArrays.axes(all_params,1)),
   AxisArray(hcat(lp...), AxisArrays.axes(C,1), AxisArrays.axes(all_params,1)))
end

function map_components(fn,tracks::AxisArray{<:SourceTracking},
                        tracks_lp::AxisArray{<:Float64};
                        window=500ms,step=250ms)
  windows = windowing(tracks[1],length=window,step=step)

  result = map(enumerate(windows)) do (i,ixs)
    best_track = argmax(dropdims(mean(Array(tracks_lp[ixs,:]),dims=1),dims=1))
    # TODO: figure out why this isn't working and why it's
    # returning a cartesian index now
    fn(tracks[best_track][Axis{:time}(ixs)])
  end

  AxisArray(result,AxisArrays.axes(windows,Axis{:time}))
end

function mask(sp::ShammaModel.AuditorySpectrogram,
              tracks::AxisArray{<:SourceTracking},
              tracks_lp::AxisArray{<:Float64},
              settings;progressbar=false,kwds...)
  scales = settings["scales"]["values"] .* cycoct
  freql,freqh = settings["rates"]["freq_limits"] .* Hz

  cr = cortical(sp,scales=scales,progressbar=progressbar)
  cr = cr[:,:,freql .. freqh]

  mask(cr,tracks,tracks_lp;progressbar=progressbar,kwds...)
end

function mask(cr::ShammaModel.Cortical,
              tracks::AxisArray{<:SourceTracking},
              tracks_lp::AxisArray{<:Float64},
              order=1;window=500ms,step=250ms,progressbar=false)

  @assert axisdim(cr,Axis{:time}) == 1
  @assert axisdim(cr,Axis{:scale}) == 2
  @assert axisdim(cr,Axis{:freq}) == 3
  @assert size(cr)[2:3] == size(tracks[1])[1:2] "Dimension mismatch"

  windows = windowing(tracks[1],length=window,step=step)

  progress = progressbar ? Progress(length(windows),"Masking: ") : nothing
  mask_helper(cr,tracks,tracks_lp,order,windows,progress)
end

function mask_helper(cr,tracks,tracks_lp,order,windows,progress)
  y = zero(AxisArray(cr))
  norm = similar(y,real(eltype(cr)))
  norm .= zero(real(eltype(cr)))

  cohere_windows =
    collect(windowing(cr,tracks[1].params.cohere))

  for (i,ixs) = enumerate(windows)
    best_track = argmax(dropdims(mean(Array(tracks_lp[ixs,:]),dims=1),dims=1))
    components = view(tracks[best_track],:,:,:,ixs)
    sorting = sortperm(component_means(components),rev=true)
    c = sorting[order]

    for (ti,t) in enumerate(ixs)
      resh = reshape(view(components,:,:,c,ti),1,size(y,2),size(y,3))
      y[Axis{:time}(cohere_windows[t])] .+= resh
      norm[Axis{:time}(cohere_windows[t])] .+= 1
    end

    next!(progress)
  end
  y ./= max.(1,norm)
  y ./= maximum(abs,y)
  y .*= cr

  cortical(y,ShammaModel.Params(cr))
end


