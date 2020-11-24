"""
    gerchbergsaxton(Ain, Aout, kw...)

reconstruct the complex signals from their amplitude values based on the Gerchberg-Saxton algorithm.

## References

1. R. W. Gerchberg and W. O. Saxton, "A practical algorithm for the determination of the phase from image and diffraction plane pictures," Optik 35, 237 (1972)
"""
function gerchbergsaxton(Ain, Aout; 
        maxiter=100,
        projector=(fft, ifft),
        errorfunc=msd,
        tol=eps(),
        init=randn(ComplexF64, size(Ain))
    )
    Cin = init
    @. Cin = Ain * exp(im * angle(Cin))
    Cout = similar(Aout)
    for i in 1:maxiter
        Cout = projector[1](Cin)
        @. Cout = Aout * exp(im * angle(Cout))
        Cin = projector[2](Cout)
        if errorfunc(abs.(Cin), Ain) < tol
            break
        end
        @. Cin = Ain * exp(im * angle(Cin))
    end
    (Cin, Cout)
end
