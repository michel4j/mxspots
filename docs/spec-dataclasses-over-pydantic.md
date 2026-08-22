# Spec: Standard Library Dataclasses for Data Modeling in mxspots

> **Triage Label**: `ready-for-agent`

## Problem Statement

`mxspots` requires structured data models to represent detected spots, diffraction frame quality scoring metrics, and lattice indexing results. Using `pydantic` introduces an external third-party dependency, increasing installation overhead, potential breaking changes across major versions, and deployment complexity in lightweight environments.

## Solution

Replace `pydantic` with Python standard library `@dataclass` structures across all `mxspots` modules (`mxspots.findspots`, `mxspots.score`, `mxspots.index`). Utilize built-in `dataclasses` features (such as `__post_init__` for validation and `dataclasses.asdict` for serialization) to maintain strong type safety and clean interfaces while eliminating external data modeling dependencies.

## User Stories

1. As a crystallographer using `mxspots.findspots`, I want spot detection outputs to be represented as standard library dataclass instances containing `(x, y, d-spacing, intensity)`, so that I can inspect detected spots without requiring external library runtime dependencies.
2. As an automated beamline pipeline developer, I want `mxspots.score` metrics (such as spot count, signal-to-noise ratio, and resolution limit) to be exposed via standard dataclass objects, so that I can reliably parse diffraction frame quality indicators.
3. As an MX data analysis engineer, I want `mxspots.index` results (including unit cell parameters and percentage indexed) returned as dataclasses, so that lattice indexing outcomes are typed and predictable.
4. As a developer installing `mxspots`, I want `pyproject.toml` to exclude `pydantic`, so that installation is fast, minimal, and dependency conflicts are minimized.
5. As a developer building Python integrations, I want dataclass fields to feature precise Python standard library type annotations, so that IDE auto-completion and static type checkers function seamlessly.
6. As a CLI user running `mxspots.findspots`, `mxspots.score`, or `mxspots.index`, I want structured dataclass outputs to convert cleanly to JSON or dictionaries for stdout reporting, so that downstream shell tools can consume the results.
7. As a data pipeline developer, I want dataclass instantiation to perform input validation (such as ensuring positive intensity values or valid $d$-spacing thresholds) via `__post_init__`, so that invalid spot finding parameters are rejected early.
8. As a maintainer testing `mxspots`, I want data models to be lightweight and zero-dependency, so that unit tests execute quickly without third-party model overhead.
9. As a developer extending `mxspots`, I want immutable dataclass configurations where applicable, so that frame analysis settings cannot be accidentally mutated during Spot Engine operations.

## Implementation Decisions

- **Standard Library Data Models**: Use Python `dataclasses` (`@dataclass`) for all core domain structures across `mxspots.findspots`, `mxspots.score`, and `mxspots.index`.
- **Validation Strategy**: Use `__post_init__` methods on dataclasses to enforce domain constraints (such as non-negative spot intensity, positive $d$-spacing, valid frame dimensions, and percentage boundaries for percentage indexed).
- **Serialization Strategy**: Use `dataclasses.asdict` combined with standard library `json` serialization for CLI JSON outputs and logging.
- **Type Annotations**: Use standard library `typing` primitives (`float`, `int`, `list`, `tuple`, `Optional`) for all dataclass fields.
- **Dependency Elimination**: Ensure `pydantic` is completely removed from project configuration and runtime dependencies.
- **Architectural Alignment**: Fully respects ADR 0005 ([`0005-dataclasses-over-pydantic.md`](file:///home/michel/Projects/mxspots/docs/adr/0005-dataclasses-over-pydantic.md)) and domain terminology in [`CONTEXT.md`](file:///home/michel/Projects/mxspots/CONTEXT.md).

## Testing Decisions

- **Testing Seam**: Primary testing seam is at the public Python module interface level (`mxspots.findspots`, `mxspots.score`, `mxspots.index` entry points) returning standard library `@dataclass` instances.
- **Good Test Philosophy**: Tests must verify external behavior through public module interfaces and assert on the attributes and serialization behavior of returned dataclasses. Avoid testing private internal mechanisms or language built-in behavior.
- **Modules Tested**: Public interface modules (`mxspots.findspots`, `mxspots.score`, `mxspots.index`) and their exported data structures.
- **Prior Art**: Standard Python unit tests using `pytest` asserting on object field attributes, validation error handling on `__post_init__`, and dictionary/JSON conversion outputs.

## Out of Scope

- Modifying the underlying C Spot Engine or `ctypes` bindings (`libmxspots`).
- Adding external serialization or schema validation libraries (e.g. `marshmallow`, `attrs`, `msgspec`).
- Changing CLI command names or argument flag conventions (`mxspots.findspots`, `mxspots.score`, `mxspots.index`).

## Further Notes

- Respects single-context domain terms defined in `CONTEXT.md` and architectural decisions in `docs/adr/0005-dataclasses-over-pydantic.md`.
- Note on Issue Tracker Publication: The GitHub CLI (`gh`) is currently unauthenticated in this environment (`gh auth login` required). This specification has been saved to [`docs/spec-dataclasses-over-pydantic.md`](file:///home/michel/Projects/mxspots/docs/spec-dataclasses-over-pydantic.md) and formatted with the `ready-for-agent` triage label for immediate issue publication once authentication is configured.
