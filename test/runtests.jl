using Test
using ShammaModel
using AuditoryBistabilityLE
using SampledSignals
using TOML
using LinearAlgebra
using Statistics

N = 50
function len_val(x,sthresh,bthresh)
  bratio = interpolate_times(x.percepts.bratio, to=x.percepts.sratio)
  sratio = x.percepts.sratio
  thresh = similar(bratio,Bool)
  thresh .= (sratio .< sthresh) .| (bratio .< bthresh)
  percept_lengths(thresh)
end
mean_log(len,val) = any(val) ? mean(log.(len[val])) : -Inf
pstream((len,val)) = mean_log(len,val) - mean_log(len,.!val)

@testset "Stremaing Bistability" begin
  @testset "Basic Streaming" begin
    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :none,
     )

    bistable_model(10, params, "test_settings.toml", interactive=true,
                   progressbar=false)
    # uncomment to profile
    # @time bistable_model(100, params, "test_settings.toml", interactive=true,
    #                      progressbar=true)

    # TODO: test with the always-fused stimulus from Elhilali et al 2009
    # (not at immediate concern, but should be present for the final
    # mod)
    params[:Δf] = 12
    df12 = bistable_model(50, params, "test_settings.toml", 
                          interactive=true, progressbar=true)

    params[:Δf] = 0.5
    df05 = bistable_model(50, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    @test pstream(len_val(df12,2.1,0.6)) > 0.9
    @test pstream(len_val(df05,2.1,0.6)) < -0.9
  end

  @testset "Track-level Bistability" begin

    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :track,
      :τ_x        => 500ms, :c_x       => 3.0,
      :W_m_σ      => 15.0, #5.0
      :W_m_σ_t    => 8.0,   :W_m_σ_ϕ   => 6.0,
      :W_m_c      => 6.0,
      :τ_m        => 350ms, :c_m       => 100,
      :τ_a        => 3s,    :c_a       => 6,
      :τ_σ        => 500ms, :c_σ       => 0.2
     )

    params[:Δf] = 12
    df12 = bistable_model(100, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(100, params, "test_settings.toml", interactive=true,
                                  progressbar=true)

    params[:Δf] = 0.5
    df05 = bistable_model(100, params, "test_settings.toml", interactive=true,
                                  progressbar=true)


    # bare minimum tests that should pass if there's stimulus
    # selective bistability
    @test pstream(len_val(df12,2.1,0.6)) > 0.9
    @test pstream(len_val(df05,2.1,0.6)) < -0.9
    @test pstream(len_val(df12,2.1,0.6)) > pstream(len_val(df3,2.1,0.6))
    @test pstream(len_val(df3,2.1,0.6)) > pstream(len_val(df05,2.1,0.6))
    @test length(len_val(df12,2.1,0.6)[1]) < length(len_val(df3,2.1,0.6)[1])
    @test length(len_val(df05,2.1,0.6)[1]) < length(len_val(df3,2.1,0.6)[1])
  end

  @testset "Scale-level Bistability" begin
    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :scales,
      :τ_x        => 500ms, :c_x       => 3.0,
      :W_m_σ      => 15.0, #5.0
      :W_m_σ_t    => 10.0,   :W_m_σ_ϕ   => 10.0,
      :W_m_c      => 6.0,
      :τ_m        => 350ms, :c_m       => 100,
      :τ_a        => 3s,    :c_a       => 6,
      :τ_σ        => 500ms, :c_σ       => 0.2
     )

    params[:Δf] = 12
    df12 = bistable_model(100, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(100, params, "test_settings.toml", interactive=true,
                         progressbar=true)

    params[:Δf] = 0.5
    df05 = bistable_model(100, params, "test_settings.toml", interactive=true,
                          progressbar=true)


    # bare minimum tests that should pass if there's stimulus
    # selective bistability
    @test pstream(len_val(df12,1.6,0.6)) > 0.9
    @test pstream(len_val(df05,1.6,0.6)) < -0.9
    @test pstream(len_val(df12,1.6,0.6)) > pstream(len_val(df3,1.6,0.6))
    @test pstream(len_val(df3,1.6,0.6)) > pstream(len_val(df05,1.6,0.6))
    @test length(len_val(df12,2.1,0.6)[1]) < length(len_val(df3,2.1,0.6)[1])
    @test length(len_val(df05,2.1,0.6)[1]) < length(len_val(df3,2.1,0.6)[1])
  end


  @testset "Freq-level Bistability" begin
    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :freqs,
      :τ_x        => 500ms, :c_x       => 3.0,
      :W_m_σ      => 5.0,
      :W_m_σ_t    => 10.0,   :W_m_σ_ϕ   => 10.0,
      :W_m_c      => 6.0,
      :τ_m        => 350ms, :c_m       => 100,
      :τ_a        => 3s,    :c_a       => 6,
      :τ_σ        => 500ms, :c_σ       => 0.2
     )

    params[:Δf] = 12
    df12 = bistable_model(50, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(50, params, "test_settings.toml", interactive=true,
                         progressbar=true)

    params[:Δf] = 0.5
    df05 = bistable_model(50, params, "test_settings.toml", interactive=true,
                          progressbar=true)


    # bare minimum tests that should pass if there's stimulus
    # selective bistability
    @test pstream(len_val(df12,1.6,0.6)) > 0.9
    @test pstream(len_val(df05,1.6,0.6)) < -0.9
    @test pstream(len_val(df12,2.1,0.6)) > pstream(len_val(df3,2.1,0.6))
    @test pstream(len_val(df3,2.1,0.6)) > pstream(len_val(df05,2.1,0.6))
    @test length(len_val(df12,2.1,0.6)[1]) < length(len_val(df3,2.1,0.6)[1])
    @test length(len_val(df05,2.1,0.6)[1]) < length(len_val(df3,2.1,0.6)[1])
  end
end

