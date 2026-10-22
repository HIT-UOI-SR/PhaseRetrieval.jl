function _fresnel_transfer(sz, distance, wavelength, pixelsize)
    distance isa Real && isfinite(distance) ||
        throw(ArgumentError("distance must be finite and real"))
    wavelength isa Real && isfinite(wavelength) && wavelength > 0 ||
        throw(ArgumentError("wavelength must be finite and positive"))
    spacing = pixelsize isa Real ? ntuple(_ -> pixelsize, length(sz)) : Tuple(pixelsize)
    length(spacing) == length(sz) || throw(DimensionMismatch("one pixel spacing is required per dimension"))
    all(d -> d isa Real && isfinite(d) && d > 0, spacing) ||
        throw(ArgumentError("pixel spacings must be finite and positive"))
    # FFT ordering: 0, 1, ..., floor((N-1)/2), -floor(N/2), ..., -1.
    T = promote_type(typeof(float(distance)), typeof(float(wavelength)),
                     map(d -> typeof(float(d)), spacing)...)
    f2 = zeros(T, sz)
    for dim in 1:length(sz)
        n = sz[dim]
        freq = [T(k <= (n-1) ÷ 2 ? k : k-n) / (n * spacing[dim]) for k in 0:n-1]
        shape = ntuple(j -> j == dim ? n : 1, length(sz))
        f2 .+= reshape(freq .^ 2, shape)
    end
    exp.((-im * T(pi) * wavelength * distance) .* f2)
end

"""
    fresnelpropagate(field, distance; wavelength, pixelsize)

Paraxial Fresnel propagation on a uniformly sampled periodic grid using
`ifft(fft(field) .* exp.(-im*pi*wavelength*distance .* sum(frequency.^2)))`.
`pixelsize` is a positive scalar or a tuple in array-dimension order. Distances,
wavelength and pixel spacings must use the same length unit. Negative distance
is the inverse propagation. Supports one or more dimensions and unshifted FFTs.

The spatially constant carrier phase `exp(im*2pi*distance/wavelength)` is omitted.
Sampling, adequate padding, and validity of the paraxial approximation are the
caller's responsibility. The field is not modified.
"""
function fresnelpropagate(field, distance; wavelength, pixelsize)
    field isa AbstractArray{<:Number} || throw(ArgumentError("field must be a numeric array"))
    isempty(field) && throw(ArgumentError("field must not be empty"))
    Base.require_one_based_indexing(field)
    all(isfinite, field) || throw(ArgumentError("field must contain finite values"))
    transfer = _fresnel_transfer(size(field), distance, wavelength, pixelsize)
    ifft(fft(field) .* transfer)
end
