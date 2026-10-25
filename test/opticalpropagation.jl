@testset "OpticalPropagation 0.4" begin
    # An offset, rectangular, odd/even grid exercises the physical-grid metadata.
    data = [exp(-3 * (x^2 + y^2)) * cis(0.8x^2 - 0.6y)
            for x in range(-1, 1; length=9), y in range(-1, 1; length=8)]
    truth = MonoLightField2D(data; wavelength=532u"nm",
                            size=(0.8u"mm", 0.6u"mm"),
                            center=(0.2u"mm", -0.1u"mm"))
    distance = 20u"mm"
    checkgrid = function (result, reference)
        @test result isa MonoLightField2D
        @test size(result) == size(reference)
        @test physicalsize(result) == physicalsize(reference)
        @test gridcenter(result) == gridcenter(reference)
        @test samplepitch(result) == samplepitch(reference)
        @test result.wavelength == reference.wavelength
    end

    for method in (fresnel2, angularspectrum)
        @testset "$method" begin
            forward, backward = method(distance), method(-distance)
            measured = forward(truth)

            # A solution must remain a fixed point of the GS iterations.
            r1, r2 = gerchbergsaxton(abs.(truth), abs.(measured);
                                    projector=(forward, backward), init=copy(truth),
                                    maxiter=3, tol=0)
            @test r1 ≈ truth
            @test r2 ≈ measured
            for result in (r1, r2)
                @test result isa MonoLightField2D
                @test physicalsize(result) == physicalsize(truth)
                @test gridcenter(result) == gridcenter(truth)
                @test samplepitch(result) == samplepitch(truth)
                @test result.wavelength == truth.wavelength
            end

            # Retrieve from an incorrect initial phase, using ordinary matrices.
            matrixforward = x -> parent(forward(MonoLightField2D(truth; data=x)))
            matrixbackward = x -> parent(backward(MonoLightField2D(truth; data=x)))
            ain, aout = abs.(parent(truth)), abs.(parent(measured))
            init = complex.(ain)
            initial_error = msd(abs.(matrixforward(init)), aout)
            r1, r2 = gerchbergsaxton(ain, aout;
                                    projector=(matrixforward, matrixbackward),
                                    init=copy(init), maxiter=30, tol=0)
            @test abs.(r1) ≈ ain
            @test r2 ≈ matrixforward(r1)
            @test msd(abs.(matrixforward(r1)), aout) < initial_error

            # Wrapped measurements retain their physical grid, including when
            # the solver generates its own initial phase.
            ainfield, aoutfield = abs.(truth), abs.(measured)
            lightinit = MonoLightField2D(truth; data=init)
            saved = copy(lightinit)
            gs, predicted = gerchbergsaxton(ainfield, aoutfield;
                projector=(forward, backward), init=lightinit, maxiter=3, tol=0)
            yg, ygpredicted = yanggu(ainfield, aoutfield;
                projector=(forward, backward), normaldiag=1,
                init=lightinit, maxiter=3, inneriter=2, innertol=0, tol=0)
            @test yg ≈ gs
            @test ygpredicted ≈ predicted
            @test lightinit == saved
            for result in (gs, predicted, yg, ygpredicted)
                checkgrid(result, truth)
            end
            for solver in (gerchbergsaxton, yanggu)
                extra = solver === yanggu ? (; normaldiag=1) : (;)
                x, y = solver(ainfield, aoutfield;
                    projector=(forward, backward), maxiter=1, tol=0, extra...)
                checkgrid(x, truth)
                checkgrid(y, truth)
                @test abs.(parent(x)) ≈ ain
                @test y ≈ forward(x)
            end

            # Support and positivity predicates must remain Boolean arrays;
            # broadcasting predicates on a light field yields complex storage.
            support = falses(size(truth))
            support[3:7, 3:6] .= true
            object = MonoLightField2D(truth; data=ain .* support)
            detector = abs.(forward(object))
            for solver in (hio, raar), positive in (false, true)
                x, y = solver(detector, support; projector=(forward, backward),
                    init=lightinit, maxiter=2, tol=0, nonnegative=positive)
                refx, refy = solver(real.(parent(detector)), support;
                    projector=(matrixforward, matrixbackward), init=init,
                    maxiter=2, tol=0, nonnegative=positive)
                @test parent(x) ≈ refx
                @test parent(y) ≈ refy
                checkgrid(x, object)
                checkgrid(y, object)
                x, y = solver(detector, support; projector=(forward, backward),
                    maxiter=1, tol=0, nonnegative=positive)
                checkgrid(x, object)
                @test y ≈ forward(x)
            end

            # Adjacent distances need not be equal or positive. Compare the
            # full multi-plane iteration to the matrix-only propagation API.
            steps = (7u"mm", -2u"mm")
            propagators = [(method(d), method(-d)) for d in steps]
            lights = [truth, method(steps[1])(truth), method(sum(steps))(truth)]
            intensities = [abs2.(light) for light in lights]
            matrixpairs = [(x -> parent(f(MonoLightField2D(truth; data=x))),
                            x -> parent(b(MonoLightField2D(truth; data=x))))
                           for (f, b) in propagators]
            fields = sbmir(intensities; propagators, init=lightinit, maxiter=3, tol=0)
            reference = sbmir([real.(parent(i)) for i in intensities];
                propagators=matrixpairs, init=init, maxiter=3, tol=0)
            for (field, ref) in zip(fields, reference)
                checkgrid(field, truth)
                @test parent(field) ≈ ref
            end
            @test fields[2] ≈ propagators[1][1](fields[1])
            @test fields[3] ≈ propagators[2][1](fields[2])
            randomfields = sbmir(intensities; propagators, maxiter=1, tol=0)
            for field in randomfields
                checkgrid(field, truth)
            end
            @test lightinit == saved
        end
    end

    @testset "Measurement validation and physical grids" begin
        for data in (fill(1 + 0.1im, size(truth)), fill(-1, size(truth)), fill(NaN, size(truth)))
            bad = MonoLightField2D(truth; data)
            @test_throws ArgumentError gerchbergsaxton(bad, abs.(truth))
            @test_throws ArgumentError yanggu(abs.(truth), bad)
            @test_throws ArgumentError hio(bad, trues(size(truth)))
            @test_throws ArgumentError raar(bad, trues(size(truth)))
            @test_throws ArgumentError sbmir([abs2.(truth), bad]; propagators=[(identity, identity)])
        end
        wronggrid = MonoLightField2D(truth; center=(0u"mm", 0u"mm"))
        @test_throws DimensionMismatch gerchbergsaxton(abs.(truth), abs.(truth);
            init=wronggrid, projector=(identity, identity), maxiter=0)
        @test_throws DimensionMismatch hio(abs.(truth), trues(3, 4);
            projector=(identity, identity))
    end

    @testset "Fresnel 0.4 phase and FFT grids" begin
        wavelength = ustrip(u"m", truth.wavelength)
        widths = ustrip.(u"m", physicalsize(truth))
        for d in (0.0, 0.02, -0.02)
            result = fresnelpropagate(parent(truth), d;
                wavelength, pixelsize=widths ./ size(truth))
            # OpticalPropagation includes the constant carrier; this package's
            # fixed-grid helper explicitly omits it.
            reference = fresnel2(parent(truth), d, wavelength, widths...) .* cispi(-2d/wavelength)
            @test result ≈ reference
        end
    end

    @testset "Yang–Gu with independently sampled Rayleigh–Sommerfeld output" begin
        input = MonoLightField2D(data[1:3, 1:2]; wavelength=truth.wavelength,
            size=(0.1u"mm", 0.08u"mm"), center=gridcenter(truth))
        forward = rayleighsommerfeld(20u"mm";
            output_size=(0.12u"mm", 0.09u"mm"), output_samples=(4, 3))
        measured = forward(input)
        basis = Matrix{ComplexF64}(I, length(input), length(input))
        G = hcat([vec(parent(forward(MonoLightField2D(input;
            data=reshape(basis[:, j], size(input)))))) for j in axes(basis, 2)]...)
        # Back-propagation by -distance is not the Hermitian adjoint of this
        # finite-window, resampled operator. Use the actual discrete adjoint.
        backward = y -> MonoLightField2D(input; data=reshape(G' * vec(y), size(input)))
        diagonal = reshape(vec(sum(abs2, G; dims=1)), size(input))
        init = MonoLightField2D(input; data=abs.(parent(input)))
        field, predicted = yanggu(abs.(input), abs.(measured);
            projector=(forward, backward), normaldiag=diagonal, init,
            maxiter=3, inneriter=2, innertol=0, tol=0)
        ref, refpredicted = yanggu(abs.(parent(input)), abs.(parent(measured)), G;
            init=parent(init), maxiter=3, inneriter=2, innertol=0, tol=0)
        @test parent(field) ≈ ref
        @test parent(predicted) ≈ refpredicted
        @test predicted ≈ forward(field)
        checkgrid(field, input)
        checkgrid(predicted, measured)
    end
end
