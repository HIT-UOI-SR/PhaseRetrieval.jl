@testset "Yang–Gu" begin
    rng = MersenneTwister(46)
    # A complex rectangular nonunitary transform catches inverse/adjoint mistakes.
    G = randn(rng, ComplexF64, 9, 5)
    truth = randn(rng, ComplexF64, 5)
    ain, aout = abs.(truth), abs.(G * truth)
    init = randn(rng, ComplexF64, 5)
    snapshots = deepcopy((ain, aout, G, init))
    expected = amplitude(ain, init)
    v = amplitude(aout, G * expected)
    W = G' * G
    offdiag = W - Diagonal(diag(W))
    for _ in 1:3
        expected = amplitude(ain, G' * v - offdiag * expected)
    end
    x, y = yanggu(ain, aout, G; init=init, maxiter=1, inneriter=3, tol=0, innertol=0)
    @test x ≈ expected
    @test y ≈ G * x
    @test (ain, aout, G, init) == snapshots
    xop, yop = yanggu(ain, aout; projector=(x -> G*x, y -> G'*y),
        normaldiag=real.(diag(W)), init=init, maxiter=1, inneriter=3, tol=0, innertol=0)
    @test xop ≈ x
    @test yop ≈ y

    # Scaled-unitary reduction to GS, including FFT's non-unit normalization.
    truth2 = randn(rng, ComplexF64, 5, 4)
    init2 = randn(rng, ComplexF64, size(truth2))
    a, b = abs.(truth2), abs.(fft(truth2))
    xgs, ygs = gerchbergsaxton(a, b; init=init2, maxiter=6, tol=0)
    xyg, yyg = yanggu(a, b; init=init2, maxiter=6, inneriter=4, tol=0, innertol=0)
    @test xyg ≈ xgs
    @test yyg ≈ ygs

    # A feasible field is fixed even for a nonunitary transform; ordinary
    # adjoint-plus-modulus GS would not have this property.
    x, y = yanggu(ain, aout, G; init=truth, maxiter=4, tol=0)
    @test x ≈ truth
    @test y ≈ G * truth
    @test norm(amplitude(ain, G' * (G * truth)) - truth) > 0.1

    # Convergence from unknown phase with an overdetermined, mildly nonunitary G.
    U = Matrix(qr(randn(rng, ComplexF64, 30, 5)).Q)[:, 1:5]
    H = U * (Matrix{ComplexF64}(I, 5, 5) + 0.05randn(rng, ComplexF64, 5, 5))
    x, y = yanggu(ain, abs.(H * truth), H; init=ones(5), maxiter=400, tol=1e-16)
    @test msd(abs.(y), abs.(H * truth)) < 1e-10
    @test phase_error(x, truth) < 1e-4

    # Singular transforms and completely unobserved coordinates are valid.
    H[:, 3] .= 0
    x, y = yanggu(ain, abs.(H*truth), H; init=truth, maxiter=3, tol=0)
    @test x ≈ truth
    @test all(isfinite, x)
    x, y = yanggu(ain, zeros(3), zeros(3,5); init=init, maxiter=2, tol=0)
    @test x ≈ amplitude(ain, init)
    @test iszero(y)
end
