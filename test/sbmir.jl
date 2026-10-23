@testset "Fresnel propagation" begin
    rng = MersenneTwister(91)
    for shape in ((9,), (8,), (7,8), (5,4,3))
        x = randn(rng, ComplexF64, shape)
        saved = copy(x)
        p = (wavelength=0.03, pixelsize=ntuple(j -> 0.4 + j/10, length(shape)))
        y = fresnelpropagate(x, 2.0; p...)
        @test fresnelpropagate(y, -2.0; p...) ≈ x
        @test fresnelpropagate(x, 0.0; p...) ≈ x
        @test norm(y) ≈ norm(x)
        @test fresnelpropagate(y, 3.0; p...) ≈ fresnelpropagate(x, 5.0; p...)
        @test x == saved
    end
    # An analytic tilted plane wave catches sign, frequency ordering and units.
    for (n,m) in ((7,8), (8,7)), (kx,ky) in ((1,-2), (-3,2), (0,0))
        dx, dy, wavelength, z = 0.2, 0.3, 0.01, 1.7
        wave = [cis(2pi*(kx*i/n + ky*j/m)) for i in 0:n-1, j in 0:m-1]
        factor = cis(-pi*wavelength*z*((kx/(n*dx))^2 + (ky/(m*dy))^2))
        @test fresnelpropagate(wave, z; wavelength=wavelength, pixelsize=(dx,dy)) ≈ wave .* factor
    end
end

@testset "SBMIR" begin
    rng = MersenneTwister(104)
    n = 8
    U = Matrix(qr(randn(rng, ComplexF64, n,n)).Q)
    V = Matrix(qr(randn(rng, ComplexF64, n,n)).Q)
    truth = randn(rng, ComplexF64, n)
    intensities = abs2.([truth U*truth V*U*truth])
    intensities = [intensities[:,j] for j in 1:3]
    init = randn(rng, ComplexF64, n)
    saved = deepcopy((intensities, init))
    propagators = [(x -> U*x, x -> U'*x), (x -> V*x, x -> V'*x)]
    a = sqrt.(hcat(intensities...))
    expected = amplitude(a[:,1], init)
    expected = amplitude(a[:,2], U * expected)
    expected = amplitude(a[:,3], V * expected)
    expected = amplitude(a[:,2], V' * expected)
    expected = amplitude(a[:,1], U' * expected)
    fields = sbmir(intensities; propagators=propagators, init=init, maxiter=1, tol=0)
    @test fields[1] ≈ expected
    @test fields[2] ≈ U * expected
    @test fields[3] ≈ V * U * expected
    @test (intensities, init) == saved
    fields0 = sbmir(intensities; propagators=propagators, init=init, maxiter=0)
    @test abs2.(fields0[1]) ≈ intensities[1]
    fields = sbmir(intensities; propagators=propagators, init=truth, maxiter=3, tol=0)
    @test fields[1] ≈ truth

    # Two-plane SBMIR has exactly the GS update with square-rooted intensities.
    fields = sbmir(intensities[1:2]; propagators=propagators[1:1], init=init, maxiter=4, tol=0)
    gs, _ = gerchbergsaxton(a[:,1], a[:,2]; projector=propagators[1], init=init, maxiter=4, tol=0)
    @test fields[1] ≈ gs

    # Unevenly spaced physical measurement planes and recovery from flat phase.
    truth = (0.6 .+ rand(rng, 8,9)) .* exp.(0.7im .* randn(rng, 8,9))
    positions = [0.0, 0.8, 2.1, 3.9, 6.2]
    p = (wavelength=0.08, pixelsize=(0.2,0.25))
    measurements = [abs2.(fresnelpropagate(truth, z; p...)) for z in positions]
    fields = sbmir(measurements, positions; p..., init=ones(size(truth)), maxiter=250, tol=1e-14)
    @test maximum(msd(abs2.(f), intensity) for (f,intensity) in zip(fields,measurements)) < 1e-8
    @test phase_error(fields[1], truth) < 1e-3
    for j in 2:length(positions)
        @test fields[j] ≈ fresnelpropagate(fields[1], positions[j]; p...)
    end
    reverse_fields = sbmir(reverse(measurements), reverse(positions); p...,
        init=fresnelpropagate(truth, positions[end]; p...), maxiter=2, tol=0)
    @test reverse_fields[end] ≈ truth

    # Every plane participates in the stopping test; imposed amplitudes alone
    # would incorrectly report zero error for inconsistent measurements.
    calls = Ref(0)
    counted = x -> (calls[] += 1; copy(x))
    sbmir([ones(2), fill(2.,2), ones(2)]; propagators=[(counted,copy),(counted,copy)],
          init=ones(2), maxiter=3, tol=1e-10)
    @test calls[] == 14 # two initial predictions + (two forward + two predictions)*3

    # Generic operators may connect arrays of different dimensionality.
    fields = sbmir([ones(4), ones(2,2)];
        propagators=[(x -> reshape(x,2,2), vec)], init=ones(4), maxiter=1)
    @test size.(fields) == [(4,), (2,2)]
end
