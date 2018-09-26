export bistable_model, as_namedtuple

function bistable_model(stim_count::Int,params,settings;kwds...)
  # stimulus generation
  stim = ab(params[:Δt]/2,params[:Δt]/2,1,stim_count,params[:f],params[:Δf]) |>
         normpower |> amplify(-10dB)
  bistable_model(stim,params,settings;kwds...)
end

read_settings(x) = as_namedtuple(x)
read_settings(x::String) = as_namedtuple(TOML.parsefile(x))
as_namedtuple(xs) = xs
function as_namedtuple(xs::Dict{<:AbstractString,<:Any})
  kt = ((Symbol(x) for x in keys(xs))...,)
  vt = (values(xs)...,)
  NamedTuple{(kt...,)}(as_namedtuple.(vt))
end

function bistable_model(stim::AbstractArray,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)
  settings = read_settings(settings)

  # auditory spectrogram
  spect = audiospect(stim,progressbar=progressbar; settings.freqs.analyze...)
  spectat = apply_bistable(spect,:freqs,params,progressbar=progressbar,
                           intermediate_results=intermediate_results;
                           settings.freqs.bistable...)
  specta = spectat.result

  # cortical scales
  cs = cortical(specta, progressbar=progressbar; settings.scales.analyze...)
  csclean = cortical(spect, progressbar=progressbar; settings.scales.analyze...)
  csat = apply_bistable(cs,:scales,params,progressbar=progressbar,
                        intermediate_results=intermediate_results;
                        settings.scales.bistable...)

  # cortical rates
  csa = csat.result
  crs = cortical(csa, progressbar=progressbar; settings.rates...)

  # temporal coherence (simultaneous grouping)
  C = cohere(crs, method=:nmf, progressbar=progressbar;
             settings.nmf...)

  # source tracking (sequential grouping)
  tracks,track_lp,groupings = track(C, method=:multi_prior, 
                                   progressbar=progressbar;
                                   settings.track.analyze...)

  track_lp_at = apply_bistable(track_lp,:track,params,
                               progressbar=progressbar,
                               intermediate_results=intermediate_results;
                               settings.track.bistable...)
  track_lp = track_lp_at.result

  # decision making
  sratio = component_ratio(
    tracks,track_lp,
    window=settings.percept_lengths.window_ms.*ms,
    step=settings.percept_lengths.delta_ms.*ms,
    intermediate_results=intermediate_results,
  )

  startHz, stopHz = settings.rates.freq_limits_Hz.*Hz
  bratio, mask = component_bandwidth_ratio(
    csclean[:,:,startHz .. stopHz],
    tracks,track_lp,
    window=settings.percept_lengths.window_ms.*ms,
    step=settings.percept_lengths.delta_ms.*ms,
    threshold=settings.percept_lengths.bandwidth_threshold
  )

  threshold = settings.percept_lengths.threshold
  # the counts are not perfect at this point but they are used to diagnose
  # serious problems in a simulation, subsequent analysis will examine `sratio`
  # and `bratio` across various thresholds
  counts = percept_lengths(AxisArray(sratio .< threshold,
                                     AxisArrays.axes(sratio,Axis{:time})),
                           settings.percept_lengths.min_length_ms.*ms)

  if intermediate_results
    (primary_source=mask,
     percepts=(counts=counts,sratio=sratio,bratio=bratio),
     sources=merge((tracks=tracks,groupings=groupings),track_lp_at),
     cohere=C,
     cortical=merge((clean=csclean[:,:,startHz .. stopHz],),csat),
     spect=spectat)
  else
    (percepts=(counts=counts,sratio=sratio,bratio=bratio),)
  end
end

