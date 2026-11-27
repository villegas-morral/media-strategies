# Instructions for Claude on the workstation

Run these steps in order. Do not change any source code (`SocialMediaModel-fork/src/`, `cases/`) and do not touch `figures/`. If something fails, stop and report; don't improvise fixes.

## Context

- `SocialMediaModel-fork/` is a Julia package. The PDE solver (`src/pde.jl`, "balance fluxes at boundary") used to add the wall correction with an extra factor `D = σ²/2` (`D/dx`, `D/dy`). It now uses `1/dx`, `1/dy`. Details: `code-bug.md`.
- `figures/` holds all figures made with the **fixed** code, in LaTeX style. That includes the old-vs-fixed comparison (`figures/img-default-ns-2inf_comparison/`).
- **Your job:** make `figures-bug/`, the same set of case figures but with the PDE computed with the **old (buggy)** term. Same LaTeX style. No comparison folder.
- Plots are saved relative to the working directory: `img-<name>/` inside `SocialMediaModel-fork/`. So run Julia from `SocialMediaModel-fork/`.
- `include` inside a script resolves relative paths against the script's own folder, so use absolute paths (the script below does).

## 1. Setup

If the repo is already cloned:

```bash
git fetch && git checkout boundary-fix && git pull
```

Julia **1.8.5** only (the Manifest fails on newer versions). If `juliaup` is missing, install it first with `curl -fsSL https://install.julialang.org | sh` (then open a new shell):

```bash
juliaup add 1.8.5
```

If the workstation has no display (SSH, server), GR must write to files instead of opening windows. Set this in the shell that runs Julia:

```bash
export GKSwstype=100
```

```bash
cd SocialMediaModel-fork && julia +1.8.5 --project=. -e 'using Pkg; Pkg.instantiate()'
```

Check that `SocialMediaModel-fork/img/` exists and is empty, and that there are no `SocialMediaModel-fork/img-*` folders.

## 2. Run the bug versions (about 30 min)

Write this script to a file **outside the repo** (e.g. `/tmp/figures_bug.jl`), then run it from `SocialMediaModel-fork/` in the background with output to a log:

```bash
julia +1.8.5 --project=. /tmp/figures_bug.jl > /tmp/figures_bug.log 2>&1
```

```julia
using SocialMediaModel
const SMM = SocialMediaModel
R = pwd()                                   # must be .../SocialMediaModel-fork
@assert basename(R) == "SocialMediaModel-fork"

# old (buggy) boundary term: put the factor D back, in a temporary copy of pde.jl
src = read(joinpath(R, "src/pde.jl"), String)
old = replace(src, "+= -1/dx *" => "+= -D/dx *", "+= 1/dx *" => "+= D/dx *",
                   "+= -1/dy *" => "+= -D/dy *", "+= 1/dy *" => "+= D/dy *")
@assert count("D/dx *", old) == 2 && count("D/dy *", old) == 2
tmp = joinpath(tempdir(), "pde_old.jl"); write(tmp, old)
Base.include(SMM, tmp)

step(f, l) = try f(); println("OK ", l); flush(stdout) catch e; println("ERROR in ", l, ": ", sprint(showerror, e)); flush(stdout) end

# expected total mass at the end of each run with the OLD code (to confirm the old term is active)
expected = Dict("default" => 1.0113, "default-ns" => 1.0112, "default-ns-closer" => 1.0182,
                "default-ns-2inf" => 1.0418, "default-t10" => 1.0113)
for c in ["default", "default-ns", "default-ns-closer", "default-ns-2inf", "default-t10"]
    step(c) do
        sol, P = include(joinpath(R, "cases", c * ".jl"))   # each case script ends with PDEsolveplot
        m = sum(sol(sol.t[end]).x[1]) * P[1].dV
        println("MASS ", c, " = ", round(m, digits=4), "  (expected ", expected[c], ")"); flush(stdout)
    end
end
step("default gif") do
    xs, si, infs, meds, st, P = ABMsolve(200; init="4inf", q=parameters(b=2., c=4., switching=true))
    ABMgifsingle(xs, si, st, infs, meds, P; name="default")
end
step("ensemble default") do
    runensembles(100; q=parameters(b=2., c=4., switching=true), name="default")
end
println("ALL DONE")
```

When it finishes, check the log:
- Every step says `OK`.
- Each `MASS` line matches its expected value to about ±0.001. If it doesn't, the old term was not active: stop and report.

## 3. Check and organise

1. The run creates these folders in `SocialMediaModel-fork/`:
   - `img-default/` with `default_abm.pdf`, `default_abm-props.pdf`, `default_abm.gif`, `default_pde.pdf`, `default_pde-ensemble.pdf`;
   - `img-default-ns/`, `img-default-ns-closer/`, `img-default-ns-2inf/` and `img-default-t10/`, each with `<name>_abm.pdf`, `<name>_abm-props.pdf`, `<name>_pde.pdf`.
2. Look at `img-default-ns-2inf/default-ns-2inf_pde.pdf`. The panel at t = 0.1 should show artefacts near the walls; that's the bug, and it's expected.
3. Move every `SocialMediaModel-fork/img-<name>/` folder to `figures-bug/img-<name>/` at the repo root. Afterwards `SocialMediaModel-fork/` should have only the empty `img/`.
4. The ABM figures should be identical to the ones in `figures/` (the ABM doesn't use the PDE boundary term). Don't replace anything in `figures/`.

## 4. Commit and push

```bash
git add figures-bug && git commit -m "Add figures-bug: case figures with the old PDE boundary term (LaTeX style)" && git push
```

End the commit message with the attribution line your session gives you, if any.

Then add a short **Log** section at the end of this file: date, run time, the `MASS` lines, and anything unexpected. Commit and push that too.

## Log

**2026-11-26**, Windows 11 workstation, Julia 1.8.5 (juliaup). Run time 1676 s (about 28 min). Every step `OK`, then `ALL DONE`.

```
MASS default = 1.0113  (expected 1.0113)
MASS default-ns = 1.0112  (expected 1.0112)
MASS default-ns-closer = 1.0182  (expected 1.0182)
MASS default-ns-2inf = 1.0418  (expected 1.0418)
MASS default-t10 = 1.0113  (expected 1.0113)
```

- All expected files were made and moved to `figures-bug/`. `default-ns-2inf_pde.pdf` shows the wall artefacts at t = 0.1.
- The ABM PDFs render pixel for pixel the same as the ones in `figures/` (the files differ only in metadata). `default_abm.gif` is identical byte for byte.

Unexpected:

- `Pkg.instantiate()` failed while building PyCall: the pinned Conda.jl 1.7.0 downloads the Miniforge installer into the folder it installs to, and the current Miniforge (26.9.4) refuses a folder that isn't empty. Fixed without touching the Manifest: `ENV["PYTHON"]` points to the system Python 3.12.10, with `setuptools` installed there (PyCall 1.94.1 imports `distutils`, which Python 3.12 removed), then `Pkg.build("PyCall")`.
- PyPlot 2.11 doesn't load with matplotlib 3.11 (`get_cmap` was removed). Not needed: the package plots with GR only.
- `SocialMediaModel-fork/img/` doesn't exist (git doesn't keep empty folders). Harmless: the code makes it with `mkpath` when needed.
- The first push failed because this machine had no GitHub credentials; it worked after `gh auth login` and `gh auth setup-git`.
