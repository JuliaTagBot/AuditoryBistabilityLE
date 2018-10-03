export bistable_model, as_namedtuple

read_settings(x) = as_namedtuple(x)
read_settings(x::String) = as_namedtuple(TOML.parsefile(x))
as_namedtuple(xs) = xs
function as_namedtuple(xs::Dict{<:AbstractString,<:Any})
  kt = ((Symbol(x) for x in keys(xs))...,)
  vt = (values(xs)...,)
  NamedTuple{(kt...,)}(as_namedtuple.(vt))
end

const cache = Dict{Vector{<:Number},AbstractMatrix}()
function bistable_model(stim_count::Int,params,settings;kwds...)
  settings = read_settings(settings)

  # cache stimulus generation
  spect = get!(cache, [stim_count,params[:Δt],params[:Δf],params[:f]]) do 
    stim = ab(params[:Δt]/2,params[:Δt]/2,1,stim_count,params[:f],params[:Δf]) |>
      normpower |> amplify(-10dB)
    audiospect(stim,progressbar=false; settings.freqs.analyze...)
  end

  bistable_model(spect,params,settings;kwds...)
end

function bistable_model(stim::AbstractVector,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)
  settings = read_settings(settings)
  spect = audiospect(stim,progressbar=progressbar; settings.freqs.analyze...)
  bistable_model(spect,params,settings;progressbar=progressbar,
                 intermediate_results=intermediate_results)
end

function bistable_model(spect::AbstractMatrix,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)
  settings = read_settings(settings)

  # auditory spectrogram
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
  tracks,track_lp = track(C, method=:multi_prior, progressbar=progressbar;
                          settings.track.analyze...)

  track_lp_at = apply_bistable!(track_lp,:track,params,
                               progressbar=progressbar,
                               intermediate_results=intermediate_results;
                               settings.track.bistable...)
  track_lp = track_lp_at.result

  # compute the mask for the primary source
  startHz, stopHz = settings.rates.freq_limits_Hz.*Hz
  crmask = mask(csclean[:,:,startHz .. stopHz],tracks,track_lp;settings.mask...)
  spmask = audiospect(crmask,progressbar=progressbar)

  # compute the ratio of the mask's and the scene's bandwidth
  ratio = bandwidth_ratio(spmask, spect; settings.bandwidth_ratio...)

  if intermediate_results
    (percepts=(ratio=ratio,counts=percept_lengths(ratio,settings)),
     primary_source=spmask,
     sources=merge((tracks=tracks,),track_lp_at),
     cohere=C,
     cortical=csat,
     spect=spectat,
     input=spect)
  else
    (percepts=(ratio=ratio,counts=percept_lengths(ratio,settings)),
     primary_source=spmask)
  end
end

