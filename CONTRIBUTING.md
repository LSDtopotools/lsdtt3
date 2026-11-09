# Contributing to lsdtt3

## Pull requests: collaborators only, for now

**Until the code base is stable, pull requests are accepted only from
repository collaborators.** The repository is set so that only collaborators
can open them. When the code base counts as stable is decided by Simon M. Mudd
(SMM); this page will say so when that changes.

The reason is practical. lsdtt3 is large and its parts depend on each other:
CLI parameter names, output formats, the Python API and the QGIS plugin are
all still moving, and a change in one place often needs matching changes
elsewhere. Reviewing outside changes against a moving target costs more than
it gains at this stage.

## What you can do now

- **Report a bug or ask for a feature** by opening an
  [issue](https://github.com/LSDtopotools/lsdtt3/issues). For a bug, include
  the tool, the output of `<tool> -v` (it reports the lsdtt3, GDAL and PROJ
  versions), the parameter file, and what you expected to happen.
- **Fork the repository and change it for your own use.** lsdtt3 is released
  under the GNU Affero General Public Licence v3 or later (see `LICENSE`).
- **Become a collaborator** if you want to work on the code with us: get in
  touch with Simon M. Mudd.

## If you are a collaborator

Work on a branch and open a pull request against `main`; never commit to
`main` directly. The workflow is in
[Git from the command line](docs/developer_notes/git_notes.md), and the
conventions code must follow are in
[Developer onboarding](docs/developer_notes/developer_onboarding.md),
[the architecture document](docs/developer_notes/lsdtt3_architecture.md) and
[parameter naming conventions](docs/developer_notes/cli_parameter_conventions.md).
