"""
Set defaults for plots to be imported in two-column LaTeX
"""
const columnWidth = 246.0 #LaTeX column width in points
PyPlot.svg(true) #For SVG plots in IJulia
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 10
rcParams["axes.labelsize"] = 6
rcParams["legend.fontsize"] = 6
rcParams["xtick.labelsize"] = 8
rcParams["ytick.labelsize"] = 8
rcParams["text.usetex"] = true
rcParams["figure.dpi"] = 100
rcParams["font.family"] = "serif"

#Importing the path effects in matplotlib seperately
pathEffects = pyimport("matplotlib.patheffects")

"""
  function setSize(; width=columnWidth, fraction=1.0, subplots=(1, 1), ratio=(sqrt(5) - 1.0) / 2.0)
Resize the figure to look better in LaTeX
"""
function setSize(; width=columnWidth, fraction=1.0, subplots=(1, 1), ratio=(sqrt(5) - 1.0) / 2.0)
    widthPt = width * fraction
    inchesPerPoint = 1 / 72.27
    width = widthPt * inchesPerPoint
    height = width * ratio * first(subplots) / last(subplots)
    return width, height
end

function redefineColorScheme()
    #Paul Tol colours. see https://personal.sron.nl/~pault/
    colors = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
    rcParams["axes.prop_cycle"] = matplotlib.cycler(color=colors)
    return colours
end

redefineColorScheme()
