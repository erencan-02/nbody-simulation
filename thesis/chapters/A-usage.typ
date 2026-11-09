#import "../lib/utils.typ": *

= User Guide <app:usage>

== Installation

```bash
pip install -e ".[viewer,experiments,dev]"   # from the repository root
python -m pytest                               # run the 103 tests
```

The package requires Python $>=$ 3.10, NumPy and Matplotlib. pygame is needed only for the viewer, SciPy only for the experiments and tests.

== Command-line interface

```bash
python -m nbody list                                   # scenarios and integrators
python -m nbody simulate kepler -i yoshida4 --dt 0.01 --t-end 100 --plot
python -m nbody simulate kepler -p e=0.9 mass_ratio=0.01 --out run.npz
python -m nbody compare figure_eight -i leapfrog rk4 yoshida4 yoshida6a
python -m nbody view solar_system -i leapfrog --steps-per-frame 4
python -m nbody view scenarios/three_body_legacy.json
```

In the viewer: `SPACE` pause, `T` trails, `P` forecast, `V` velocity vectors, `C` follow the centre of mass, `F` fit view, `+`/`-`/mouse wheel zoom, arrow keys pan, `.`/`,` faster/slower, `1`--`9` switch integrator, `R` reset, `S` screenshot, `H` help.

== Python API

```python
from nbody import load_scenario, simulate
from nbody.analysis import conservation_report

sc = load_scenario("kepler", e=0.6)
traj = simulate(sc.model(), sc.initial_state(), "yoshida4", t_end=100.0, dt=0.01)
print(conservation_report(traj))   # conservation errors
traj.save("kepler.npz")
```

== Scenario file format

```json
{
  "name": "binary",
  "description": "Equal-mass binary on a circular orbit",
  "G": 1.0,
  "softening": 0.0,
  "units": {"length": "", "time": "", "mass": ""},
  "recommended": {"integrator": "leapfrog", "dt": 0.01, "t_end": 50.0},
  "bodies": [
    {"name": "A", "mass": 0.5, "color": "#0083CC",
     "position": [-0.5, 0], "velocity": [0, -0.5]},
    {"name": "B", "mass": 0.5, "color": "#E6001A",
     "position": [0.5, 0], "velocity": [0, 0.5]}
  ]
}
```

All keys except `bodies` are optional, so the files of the original course project (which contain only `G` and `bodies`) remain valid. Positions and velocities may be two- or three-dimensional.

== Reproducing the results

```bash
python scripts/fetch_fonts.py                # template fonts
python -m experiments.run_all                # all experiments (~6 min)
typst compile thesis/main.typ --font-path thesis/assets/fonts
```
