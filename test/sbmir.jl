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
    widths = size(truth) .* p.pixelsize
    measurements = [abs2.(fresnel2(truth, z, p.wavelength, widths...)) for z in positions]
    fields = sbmir(measurements, positions; p..., init=ones(size(truth)), maxiter=250, tol=1e-14)
    @test maximum(msd(abs2.(f), intensity) for (f,intensity) in zip(fields,measurements)) < 1e-8
    @test phase_error(fields[1], truth) < 1e-3
    for j in 2:length(positions)
        @test fields[j] ≈ fresnel2(fields[1], positions[j], p.wavelength, widths...)
    end
    reverse_fields = sbmir(reverse(measurements), reverse(positions); p...,
        init=fresnel2(truth, positions[end], p.wavelength, widths...), maxiter=2, tol=0)
    @test reverse_fields[end] ≈ truth

    # The matrix convenience interface must agree with OpticalPropagation's
    # physical-grid interface, including its complex carrier phase.
    light = MonoLightField2D(truth; wavelength=p.wavelength * u"m",
        size=widths .* u"m")
    pairs = [(fresnel2(d * u"m"), fresnel2(-d * u"m")) for d in diff(positions)]
    lightmeasurements = [abs2.(fresnel2(z * u"m")(light)) for z in positions]
    wrapped = sbmir(lightmeasurements; propagators=pairs,
        init=MonoLightField2D(light; data=ones(size(light))), maxiter=3, tol=0)
    matrices = sbmir(measurements, positions; p..., init=ones(size(truth)), maxiter=3, tol=0)
    @test all(parent(w) ≈ m for (w, m) in zip(wrapped, matrices))

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
