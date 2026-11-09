"""Command-line interface.

Examples::

    python -m nbody list
    python -m nbody simulate kepler --integrator yoshida4 --dt 0.01 --t-end 100 --plot
    python -m nbody compare kepler -i euler leapfrog rk4 yoshida4 --dt 0.02 --t-end 200
    python -m nbody view solar_system --integrator leapfrog --steps-per-frame 4
    python -m nbody view scenarios/two_body_legacy.json
"""

from __future__ import annotations

import argparse
import json
import sys

from .analysis import conservation_report
from .integrators import available_integrators, create_integrator
from .scenarios import available_scenarios, load_scenario, scenario_description
from .simulation import Simulation


def _parse_params(items: list[str] | None) -> dict:
    """Parse ``key=value`` scenario parameters (values are JSON-decoded if possible)."""
    params = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"scenario parameter '{item}' must have the form key=value")
        key, value = item.split("=", 1)
        try:
            params[key] = json.loads(value)
        except json.JSONDecodeError:
            params[key] = value
    return params


def _integrator(name: str, rtol: float | None, atol: float | None):
    kwargs = {}
    if rtol is not None:
        kwargs["rtol"] = rtol
    if atol is not None:
        kwargs["atol"] = atol
    return create_integrator(name, **kwargs)


def cmd_list(_args) -> int:
    print("Scenarios:")
    for name in available_scenarios():
        print(f"  {name:<20s} {scenario_description(name)}")
    print("\nIntegrators (increasing sophistication):")
    for name in available_integrators():
        p = create_integrator(name).properties()
        print(
            f"  {name:<18s} order {p['order']}  "
            f"{'symplectic' if p['symplectic'] else 'non-symplectic':<15s} "
            f"{'reversible' if p['time_reversible'] else '':<11s}"
            f"{p['evals_per_step']} evals/step{'  (adaptive)' if p['adaptive'] else ''}"
        )
    return 0


def _resolve(args, sc):
    dt = args.dt if args.dt is not None else sc.recommended.get("dt")
    t_end = args.t_end if args.t_end is not None else sc.recommended.get("t_end")
    if dt is None or t_end is None:
        raise SystemExit("this scenario has no recommended dt/t_end; pass --dt and --t-end")
    return float(dt), float(t_end)


def cmd_simulate(args) -> int:
    sc = load_scenario(args.scenario, **_parse_params(args.param))
    dt, t_end = _resolve(args, sc)
    name = args.integrator or sc.recommended.get("integrator", "leapfrog")
    sim = Simulation(sc.model(), sc.initial_state(), _integrator(name, args.rtol, args.atol))
    traj = sim.run(t_end, dt, record_every=args.record_every)
    traj.meta["scenario"] = sc.name
    rep = conservation_report(traj)
    print(f"scenario    {sc.name} (N={sc.n_bodies})")
    print(f"integrator  {sim.integrator.label}")
    print(f"steps       {traj.meta['n_steps']}   force evaluations {traj.meta['n_force_evals']}")
    print(f"wall time   {traj.meta['wall_time']:.3f} s")
    for key, value in rep.as_dict().items():
        print(f"  {key:<22s} {value: .3e}")
    if args.out:
        traj.save(args.out)
        print(f"trajectory written to {args.out}")
    if args.plot:
        import matplotlib.pyplot as plt

        from .viz.plotting import quicklook

        quicklook(traj, sc.names, sc.colors, title=f"{sc.name} - {sim.integrator.label}")
        plt.show()
    return 0


def cmd_compare(args) -> int:
    sc = load_scenario(args.scenario, **_parse_params(args.param))
    dt, t_end = _resolve(args, sc)
    header = f"{'integrator':<26s}{'evals':>10s}{'max|dE/E0|':>13s}{'drift/t':>12s}{'max|dL/L0|':>13s}{'time[s]':>9s}"
    print(f"{sc.name}: t_end={t_end:g}, dt={dt:g}\n{header}\n{'-' * len(header)}")
    for name in args.integrators:
        sim = Simulation(sc.model(), sc.initial_state(), _integrator(name, args.rtol, args.atol))
        try:
            traj = sim.run(t_end, dt, record_every=args.record_every)
        except (RuntimeError, FloatingPointError) as exc:
            print(f"{sim.integrator.label:<26s} failed: {exc}")
            continue
        r = conservation_report(traj)
        print(
            f"{sim.integrator.label:<26s}{traj.meta['n_force_evals']:>10d}"
            f"{r.energy_max:>13.2e}{r.energy_drift_rate:>12.2e}"
            f"{r.angular_momentum_max:>13.2e}{traj.meta['wall_time']:>9.2f}"
        )
    return 0


def cmd_view(args) -> int:
    from .viz.viewer import Viewer

    sc = load_scenario(args.scenario, **_parse_params(args.param))
    name = args.integrator or sc.recommended.get("integrator", "leapfrog")
    if create_integrator(name).adaptive:
        print("note: adaptive integrators use --dt as the maximum step in the viewer")
    viewer = Viewer(
        sc,
        integrator=name,
        dt=args.dt,
        steps_per_frame=args.steps_per_frame,
        trail_length=args.trail_length,
        prediction_horizon=args.horizon,
    )
    viewer.run(max_frames=args.frames, screenshot=args.screenshot)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="nbody", description="Gravitational N-body simulator")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="list scenarios and integrators").set_defaults(func=cmd_list)

    def common(sp, multi=False):
        sp.add_argument("scenario", help="builtin scenario name or path to a JSON file")
        if multi:
            sp.add_argument("-i", "--integrators", nargs="+", default=available_integrators())
        else:
            sp.add_argument("-i", "--integrator", default=None)
        sp.add_argument("--dt", type=float, default=None, help="(initial) step size")
        sp.add_argument("--t-end", type=float, default=None, help="final time")
        sp.add_argument("--record-every", type=int, default=10)
        sp.add_argument("--rtol", type=float, default=None, help="adaptive methods only")
        sp.add_argument("--atol", type=float, default=None, help="adaptive methods only")
        sp.add_argument("-p", "--param", nargs="*", help="scenario parameters key=value")

    s = sub.add_parser("simulate", help="run a simulation and report conservation metrics")
    common(s)
    s.add_argument("--out", help="write the trajectory to an .npz file")
    s.add_argument("--plot", action="store_true", help="show orbits and energy error")
    s.set_defaults(func=cmd_simulate)

    c = sub.add_parser("compare", help="compare integrators on one scenario")
    common(c, multi=True)
    c.set_defaults(func=cmd_compare)

    v = sub.add_parser("view", help="interactive real-time viewer with trajectory overlays")
    v.add_argument("scenario")
    v.add_argument("-i", "--integrator", default=None)
    v.add_argument("--dt", type=float, default=None)
    v.add_argument("--steps-per-frame", type=int, default=1)
    v.add_argument("--trail-length", type=int, default=400)
    v.add_argument("--horizon", type=float, default=None, help="prediction horizon")
    v.add_argument("--frames", type=int, default=None, help="stop after N frames")
    v.add_argument("--screenshot", default=None, help="save the last frame as PNG")
    v.add_argument("-p", "--param", nargs="*", help="scenario parameters key=value")
    v.set_defaults(func=cmd_view)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
