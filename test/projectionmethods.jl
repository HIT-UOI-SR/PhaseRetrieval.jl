@testset "HIO and RAAR" begin
    rng = MersenneTwister(72)
    init = randn(rng, ComplexF64, 4, 5)
    aout = abs.(fft(randn(rng, ComplexF64, 4, 5)))
    support = falses(4,5)
    support[2:3, 2:4] .= true
    saved = deepcopy((aout, support, init))
    pm = ifft(amplitude(aout, fft(init)))
    beta = 0.73
    expected = ifelse.(support, pm, init .- beta .* pm)
    x, y = hio(aout, support; init=init, maxiter=1, beta=beta, tol=0)
    @test x ≈ expected
    @test y ≈ fft(x)
    accepted = support .& (real.(pm) .>= 0)
    expected = ifelse.(accepted, real.(pm), init .- beta .* pm)
    x, _ = hio(aout, support; init=init, maxiter=1, beta=beta, nonnegative=true, tol=0)
    @test x ≈ expected

    for positive in (false, true), beta in (0.0, 0.73, 1.0)
        reflected = 2pm - init
        projected = ifelse.(support, positive ? max.(real.(reflected), 0) : reflected, 0)
        expected = beta/2 .* (2projected - reflected + init) + (1-beta)*pm
        x, y = raar(aout, support; init=init, maxiter=1, beta=beta, nonnegative=positive, tol=0)
        @test x ≈ expected
        @test y ≈ fft(x)
    end
    @test (aout, support, init) == saved

    for solver in (hio, raar)
        x, y = solver(aout, support; init=init, maxiter=0)
        @test x == init
        @test y ≈ fft(init)
        @test x !== init
        for positive in (false, true)
            truth = ifelse.(support, positive ? complex.(abs.(init)) : init, 0.0im)
            x, y = solver(abs.(fft(truth)), support; init=truth, nonnegative=positive, maxiter=3, tol=0)
            @test x ≈ truth
            @test y ≈ fft(truth)

            # A feasible initial object must stop before applying a feedback
            # update, for both complex and real nonnegative constraints.
            calls = Ref(0)
            counted = x -> (calls[] += 1; fft(x))
            x, y = solver(abs.(fft(truth)), support; init=truth,
                nonnegative=positive, maxiter=9, tol=1e-20,
                projector=(counted, ifft))
            @test calls[] == 1
            @test x == truth
            @test y ≈ fft(truth)
        end
        # Amplitude feasibility alone cannot trigger stopping outside support.
        calls = Ref(0)
        forward = x -> (calls[] += 1; fft(x))
        solver(abs.(fft(init)), support; init=init, maxiter=1, tol=1e-12, projector=(forward,ifft))
        @test calls[] == 2
    end

    # Oversampled 2D nonnegative object, initialized without its true phase.
    support = falses(24,24)
    support[9:15, 8:16] .= true
    truth = zeros(24,24)
    truth[support] .= 0.2 .+ rand(rng, count(support))
    aout = abs.(fft(truth))
    initial = rand(rng, 24,24) .* support
    initial_error = msd(abs.(fft(initial)), aout)
    for solver in (hio, raar)
        x, _ = solver(aout, support; init=initial, nonnegative=true, beta=0.95, maxiter=600, tol=0)
        estimate = ifelse.(support, max.(real.(x), 0), 0)
        @test msd(abs.(fft(estimate)), aout) < initial_error * 1e-4
    end
end
