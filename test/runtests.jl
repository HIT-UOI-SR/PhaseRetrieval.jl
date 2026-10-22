using PhaseRetrieval
using FFTW
using Distances
using OpticalPropagation
using Test
using Unitful
using LinearAlgebra
using Random

@testset "PhaseRetrieval.jl" begin
    @test begin # see Gerchberg's original paper
        c1 = [Int(abs(2t)<=1)*exp(im*15*pi*t^2) for t in -1:(1/256):(255/256)]
        c2 = fft(c1)
        r1, r2 = gerchbergsaxton(abs.(c1), abs.(c2), maxiter=130, init=c1)
        msd(c1, r1) < 0.01 && msd(c2, r2) < 0.01
    end

    include("opticalpropagation.jl")
end

# Independent reference expressions: do not use internal projection helpers.
amplitude(a, x) = a .* exp.(im .* angle.(x))
phase_error(x, truth) = norm(x .* exp(-im * angle(dot(truth, x))) - truth) / norm(truth)

include("gerchbergsaxton.jl")
include("yanggu.jl")
include("projectionmethods.jl")
include("sbmir.jl")
include("validation.jl")
