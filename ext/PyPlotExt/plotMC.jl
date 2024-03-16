function plotProperty!(ax, y; x = nothing, color = defaultColor, kwargs...)
	xAxis = x === nothing ? range(1, length(y)) : x
	yAxis = getindex.(y, 1)
	Δy = getindex.(y, 2)
	#ax.plot(xAxis, y, color=color, zorder=defaultZOrder; kwargs...)

	ax.errorbar(xAxis, yAxis, Δy, color = color, zorder = defaultZOrder, ecolor = errorColor; kwargs...)
	ax.scatter(xAxis, yAxis, marker = defaultMarker, s = defaultMarkerSize, color = color, zorder = defaultZOrder + 1; kwargs...)
end

function getSubplots(N::Int64; tile = false)
	if tile
		nearestSquare = floor(Int64, sqrt(N))
		additionalRows = ceil(Int64, (N - nearestSquare^2) / nearestSquare)
		return (nearestSquare + additionalRows, nearestSquare)
	else
		return (N, 1)
	end
end

function extendLimits(l, u; c = 0.05)
	δ = (u - l) * c
	return (l - δ, u + δ)
end

function setTicks(ax, l::T1, u::T2) where {T1, T2}
	ax.set_yticks(round.(range(l, u, length = 5), digits = 1))
end

function setTicks(ax, l::Rational, u::Rational)
	ax.set_yticks((l, u), labels = [i.den == 1 ? "$(i.num)" : "$(i.num)/$(i.den)" for i in (l, u)])
end

function plotTransitionLines(axis, temperature, observables, properties)
	CvIndex = findfirst(x -> x == :specificHeat, properties)
	Cv = getindex.(observables[:, CvIndex], 1)
	transitionIndex = argmax(Cv)
	Tc = temperature[transitionIndex]
	for ax in axis
		ax.plot([Tc, Tc], [-10, 10], c = transitionColor)
	end
end

function plotObservables!(axis, temperature::Vector{Float64}, observables::Matrix{T}, properties::Vector{Symbol};
	color = defaultColor, transition = true, markLines = true, kwargs...) where {T}
	plt.minorticks_off()
	plt.xlim(0.0, maximum(temperature))
	plt.xlabel("temperature")

	transition && (:specificHeat in properties) && plotTransitionLines(axis, temperature, observables, properties)

	for (i, (ax, prop)) in enumerate(zip(axis, properties))
		y = getindex.(observables[:, i], 1)

		ax.set_ylabel(get(labelsDict, prop, prop), wrap = true)
		ax.set_ylim(extendLimits(getLimits(y, prop)...))
		setTicks(ax, getTicks(y, prop)...)
		markLine = getMarkerLine(y, prop)
		if markLines && !(markLine === nothing)
			if length(markLine) == 2
				ax.plot(ax.get_xlim(), markLine, c = transitionColor, zorder = defaultZOrder)
			else
				ax.plot(ax.get_xlim(), (markLine, markLine), c = transitionColor, zorder = defaultZOrder)
			end
		end
		plotProperty!(ax, observables[:, i], x = temperature, color = color)
	end
end

function plotObservables(temperature::Vector{Float64}, observables::Matrix{T}, properties::Vector{Symbol};
	color = defaultColor, transition = true, markLines = true, kwargs...) where {T}
	axisSize = getSubplots(length(properties))
	fig, axis = subplots(axisSize..., sharex = true)
	axis = permutedims(axis)
	#width, height = setSize(subplots=axisSize, ratio=1.0 / sqrt(8))
	#fig.set_figwidth(width)
	#fig.set_figheight(height)

	#If plotting only single property, then axs needs to be converted to an array for zip.
	if length(properties) == 1 && setup
		axis = [axis]
	end
	plotObservables!(axis, temperature, observables, properties, color = color, transition = transition, markLines = markLines; kwargs...)
	fig.align_ylabels()
	return fig, axis
end

function plotMC(fileLocation, properties::Vector{Symbol}; saveLocation = nothing, kwargs...)
	fn = getFileNames(fileLocation)
	temperatures, β, observables = loadObservables(fn, properties)
	fig, axis = plotObservables(temperatures, observables, properties; kwargs...)

	!(saveLocation === nothing) && savefig("$(fileLocation)/obs.pdf", format = "pdf")
	return fig, axis
end
