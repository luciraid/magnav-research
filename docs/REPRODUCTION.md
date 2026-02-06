# MagNav Reproduction Notes

## Environment
- OS: Ubuntu 24.04 LTS
- Julia: 1.10.0
- MagNav.jl: installed via Julia Pkg
- MATLAB: Installed (not required for this run)

## Execution
- Example executed: `examples/simple_magnav.jl`
- Command sequence:
  ```julia
  using MagNav
  pkgroot = dirname(dirname(pathof(MagNav)))
  include(joinpath(pkgroot, "examples", "simple_magnav.jl"))
