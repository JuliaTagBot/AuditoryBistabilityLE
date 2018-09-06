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
      :f_c_a => 0, :f_c_m => 0, :f_c_σ => 0,
      :s_c_a => 0, :s_c_m => 0, :s_c_σ => 0,
      :t_W_m_σ      => 15.0, #5.0
      :t_W_m_σ_t    => 8.0,   :t_W_m_σ_ϕ   => 6.0,
      :t_W_m_c      => 6.0,
      :t_τ_m        => 350ms, :t_c_m       => 100,
      :t_τ_a        => 3s,    :t_c_a       => 6,
      :t_τ_σ        => 500ms, :t_c_σ       => 0.2
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


    # bare minimum tests that should pass if there's stimulus selective
    # bistability the results are somewhat random, so not all tests will pass
    # on all runs, but the should normally *mostly* pass.
    # once the parameters have been fit well, the goal will be to select
    # parameter values that consistently pass
    sthresh, bthresh = 2.1, 0.25
    @test pstream(len_val(df12,sthresh,bthresh)) > 0.9
    @test pstream(len_val(df05,sthresh,bthresh)) < -0.9
    @test pstream(len_val(df12,sthresh,bthresh)) > pstream(len_val(df3,sthresh,bthresh))
    @test pstream(len_val(df3,sthresh,bthresh)) > pstream(len_val(df05,sthresh,bthresh))
    @test length(len_val(df12,sthresh,bthresh)[1]) < length(len_val(df3,sthresh,bthresh)[1])
    @test length(len_val(df05,sthresh,bthresh)[1]) < length(len_val(df3,sthresh,bthresh)[1])

  end

  @testset "Scale-level Bistability" begin

    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :track,
      :τ_x        => 500ms, :c_x       => 3.0,
      :f_c_a => 0, :f_c_m => 0, :f_c_σ => 0,
      :s_W_m_σ      => 15.0,
      :s_W_m_c      => 6.0,
      :s_τ_m        => 350ms, :s_c_m       => 65,
      :s_τ_a        => 3s,    :s_c_a       => 10,
      :s_τ_σ        => 500ms, :s_c_σ       => 0.2,
      :t_c_a => 0, :t_c_m => 0, :t_c_σ => 0
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
    sthresh, bthresh = 1.6, 0.25
    @test pstream(len_val(df12,sthresh,bthresh)) > 0.9
    @test pstream(len_val(df05,sthresh,bthresh)) < -0.9
    @test pstream(len_val(df12,sthresh,bthresh)) > pstream(len_val(df3,sthresh,bthresh))
    @test pstream(len_val(df3,sthresh,bthresh)) > pstream(len_val(df05,sthresh,bthresh))
    @test length(len_val(df12,sthresh,bthresh)[1]) < length(len_val(df3,sthresh,0.25)[1])
    @test length(len_val(df05,sthresh,bthresh)[1]) < length(len_val(df3,sthresh,bthresh)[1])

  end


  @testset "Freq-level Bistability" begin

    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :track,
      :τ_x        => 500ms, :c_x       => 3.0,
      :f_W_m_σ      => 5.6,
      :f_W_m_c      => 6.0,
      :f_τ_m        => 350ms, :f_c_m       => 100,
      :f_τ_a        => 3s,    :f_c_a       => 6,
      :f_τ_σ        => 500ms, :f_c_σ       => 0.2,
      :s_c_a => 0, :s_c_m => 0, :s_c_σ => 0,
      :t_c_a => 0, :t_c_m => 0, :t_c_σ => 0
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
    sthresh, bthresh = 2.1, 0.25
    @test pstream(len_val(df12,sthresh,bthresh)) > 0.9
    @test pstream(len_val(df05,sthresh,bthresh)) < -0.9
    @test pstream(len_val(df12,sthresh,bthresh)) > pstream(len_val(df3,sthresh,bthresh))
    @test pstream(len_val(df3,sthresh,bthresh)) > pstream(len_val(df05,sthresh,bthresh))
    @test length(len_val(df12,sthresh,bthresh)[1]) < length(len_val(df3,sthresh,bthresh)[1])
    @test length(len_val(df05,sthresh,bthresh)[1]) < length(len_val(df3,sthresh,bthresh)[1])

  end
end

