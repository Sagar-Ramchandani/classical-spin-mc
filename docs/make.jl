using Documenter
using ClassicalSpinMC

makedocs(
    sitename = "ClassicalSpinMC",
    format = Documenter.HTML(repolink = "git@gitsrv-focal.thp.uni-koeln.de:sramchan/classical-spin-mc.git"),
    modules = [ClassicalSpinMC],
    repo = "git@gitsrv-focal.thp.uni-koeln.de:sramchan/classical-spin-mc.git"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
