# MagNav Research

Purpose
- Host original research, analyses, and documentation that extend or investigate results from the `MagNav` reference repository.

Background & relationship to `MagNav`
- This repository contains only original analysis artifacts and documentation derived from the `MagNav` reference codebase. The original `MagNav` repository remains the authoritative location for production code and must not be modified by this project.

Repository status
- Initialized and migrated core analysis artifacts into `docs/` and runtime logs into `logs/`.
- No production source code from the reference repo was copied.
- Current: awaiting external flight data before running further experiments.

Repository layout
- `docs/` — analysis reports and experiment notes (migrated).
- `logs/` — runtime logs and shell history (migrated).
- `data/` — (data not included) place flight datasets here following `data/README.md`.
- `src/` — analysis helpers and experimental scripts (keep original code out of this repo unless authored here).
- `notebooks/` — exploratory notebooks and reproducible analysis.

Getting started
1. Clone:

	git clone https://github.com/luciraid/magnav-research.git

2. Inspect `docs/` for migrated reports and `logs/` for run artifacts.

Data
- Flight/experimental datasets are not bundled here. See `data/README.md` for instructions to add data and expected formats.

Contributing
- See `CONTRIBUTING.md` for contribution guidelines and `CODE_OF_CONDUCT.md`.

License
- This repository is available under the MIT license. See `LICENSE` for details.

Contact
- Repository owner: https://github.com/luciraid

Notes for maintainers
- Preserve provenance: when copying artifacts from `MagNav`, keep original timestamps and add a short note in the corresponding document stating origin.

--
_Created as a research workspace derived from `MagNav`; production code remains in the original repository._
