using PhaseRetrieval
using FFTW
using Distances
using Test

@testset "PhaseRetrieval.jl" begin
    @test begin # see Gerchberg's original paper
        c1 = [Int(abs(2t)<=1)*exp(im*15*pi*t^2) for t in -1:(1/256):(255/256)]
        c2 = fft(c1)
        r1, r2 = gerchbergsaxton(abs.(c1), abs.(c2), maxiter=130, init=c)
        msd(c1, r1) < 0.01 && msd(c2, r2) < 0.01
    end
end
