# Publishing lsdtt3 to PyPI

How the wheels get to PyPI, so that `pip install lsdtt3` works. Every
`vX.Y.Z` tag publishes them from GitHub Actions; this page covers the one-time
setup, what the release workflow does, how to redo or add an upload by hand,
and the by-hand route with `twine`.

For how the wheels themselves are built, see
[`building_python_wheels.md`](building_python_wheels.md). For the local
all-artefacts script, see `tools/README.md`.

---

## What a user gets from `pip install lsdtt3`

Worth being explicit, because it is more than a Python API. A released wheel
contains:

| Component | Path in the wheel |
|---|---|
| The C++ runtime | `lsdtt3/liblsdtt3_runtime.so` |
| The bindings extension | `lsdtt3/_core.abi3.so` (`_core.pyd` on Windows) |
| Every command-line tool (FFTW-dependent `lsdtt-spectral-analysis` only when FFTW is found) | `lsdtt3/bin/lsdtt-*` |
| GDAL and PROJ support data | `lsdtt3/share/gdal`, `lsdtt3/share/proj` |
| Cosmogenic parameters | `lsdtt3/share/lsdtt3/crn_parameters.nc` |
| Vendored dependencies | `lsdtt3.libs/` (added by `auditwheel`) |

`[project.scripts]` declares a console script per CLI, so after
`pip install lsdtt3` the `lsdtt-*` tools are on `PATH` as well as the Python
API being importable. **There is no separate backend install**, and no system
GDAL is required: `auditwheel repair` (Linux), `delocate` (macOS) and
`delvewheel` (Windows) copy GDAL, PROJ, netCDF and FFTW into the wheel and
rewrite the library names to match.

!!! danger "Only ever upload a repaired wheel"

    A wheel built outside the manylinux image links the build machine's
    libraries and has no `lsdtt3.libs/`. It will import on the machine that
    built it and fail everywhere else. `tools/make_all_packages.sh` prefixes
    such a wheel `LOCAL-ONLY-` for this reason. **Never upload a
    `LOCAL-ONLY-` wheel**, and never upload anything with a `linux_x86_64`
    platform tag — PyPI requires `manylinux*` and will reject it, but it is
    better not to try.

---

## One-time setup

PyPI and TestPyPI are separate sites with separate accounts. On each, the
uploads use **trusted publishing** (OpenID Connect): GitHub proves to the index
which repository, workflow and environment a job is running in, and no token
is stored anywhere.

1. **Trusted publishers.** On <https://pypi.org> (account → *Publishing*, or
   the project's *Manage → Publishing* once it exists) add, for owner
   `LSDtopotools`, repository `lsdtt3`, project `lsdtt3`:

    | Workflow | Environment | Used by |
    |---|---|---|
    | `release.yml` | `pypi` | every `vX.Y.Z` tag |
    | `python-wheels.yml` | `pypi` | manual runs (`.postN` rebuilds, a redone upload) |

    Do the same on <https://test.pypi.org> with environment `testpypi`. Before
    the first upload the project does not exist, so add them as *pending*
    publishers. A publisher with the environment left blank accepts any
    environment; that works too.

    PyPI checks the workflow that **started** the run. `release.yml` calls
    `python-wheels.yml` as a reusable workflow, and an upload from inside a
    called workflow is not accepted, so the tag uploads are jobs of
    `release.yml` itself, and `python-wheels.yml` uploads only on a manual
    run. That is why both files are registered.

2. **GitHub environments** `testpypi` and `pypi` (repository *Settings →
   Environments*). GitHub creates an environment the first time a job uses
   it, but without rules. A **required reviewer** on `pypi` makes every PyPI
   upload wait for approval, so the TestPyPI files can be installed and tried
   first. Recommended at least for the first releases.

3. **Storage.** A PyPI project may hold 10 GB by default. A release is four
   wheels of about 35–65 MB, about 200 MB, so the default holds roughly fifty
   releases. A higher limit can be requested at
   <https://github.com/pypi/support> (*project limit request*).

The metadata on the project page comes from `python_bindings/pyproject.toml`
(`[project]`: description, authors, keywords, classifiers, URLs) and
`python_bindings/README.md`, which is the page's text, so its links must be
absolute.

---

## What a tag does

`release.yml`, on a `vX.Y.Z` tag, after the backends, the wheels and the
GitHub release:

1. `python-wheels` builds the four wheels (Linux x86_64, macOS ARM64 and
   Intel, Windows x64), tests each one on Python 3.12, 3.13 and 3.14, and
   attaches them to the release;
2. `publish-release` publishes the GitHub release as Latest;
3. `publish-testpypi` downloads the wheels from the run, checks them, and
   uploads them to TestPyPI;
4. `publish-pypi` checks and uploads the same files to PyPI (after approval,
   if the `pypi` environment requires a reviewer).

Steps 3 and 4 run **only if every wheel built**: a version reaches PyPI with
all four platforms or not at all, because a version cannot be added to
cleanly afterwards. If a platform failed, fix it and publish with a manual run
(below).

Before uploading, each job runs
`python_bindings/tools/check_dist_wheels.py`, which refuses anything but
exactly one `cp312-abi3` wheel per platform, all carrying the tag's version
and a repaired platform tag, each under PyPI's 100 MB file limit, and then
`twine check --strict`. Uploads use `skip-existing`, so re-running a job that
uploaded part-way is harmless.

To check a release from TestPyPI before approving PyPI:

```bash
python -m venv /tmp/tp && /tmp/tp/bin/pip install \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    lsdtt3==X.Y.Z
/tmp/tp/bin/python -c "import lsdtt3; print(lsdtt3.__version__)"
/tmp/tp/bin/lsdtt-raster-preprocessing -v
```

The `--extra-index-url` matters: numpy and geopandas are not on TestPyPI.

### Uploading by hand from Actions

*Actions → Build Python wheels → Run workflow*, with the tag, platforms
**all**, and `publish_to` **testpypi** or **pypi**. It rebuilds and retests
the wheels from the tag and, if all four built, uploads them (job
`publish-pypi`). Use it:

- to publish a `vX.Y.Z.postN` wheel-only rebuild (the `.postN` section
  below). A `.postN` is always rebuilt for **all** platforms, even if only
  one was broken, because a manual upload needs the full set;
- to upload a tag whose tag run did not reach PyPI, for example because one
  platform failed and was fixed.

A failed upload job of a tag run can also simply be re-run from that run's
page: the wheels are still the run's artefacts (kept for the repository's
artefact retention period).

---

## Publishing by hand with twine

For experiments on TestPyPI without touching CI. **Do not publish a release
to PyPI this way**: it builds Linux wheels only, and a version should reach
PyPI with all four platforms at once (see [What a tag does](#what-a-tag-does)).

The catch is that **you still need manylinux wheels**, which means building
inside the manylinux image. That needs a container engine, and this project
does not use Docker — install podman:

```bash
sudo apt install podman
```

### Build

```bash
# every supported Python, repaired into manylinux wheels
CIBW_CONTAINER_ENGINE=podman \
CIBW_MANYLINUX_X86_64_IMAGE=manylinux_2_28 \
CIBW_ENVIRONMENT='CMAKE_ARGS="-DLSDTT_PYTHON_RELEASE_WHEEL=ON"' \
    pipx run cibuildwheel python_bindings --output-dir wheelhouse
```

`tools/make_all_packages.sh --only wheel` does the same thing once podman is
present, and collects the result into `dist/`.

This produces **Linux wheels only**. macOS and Windows wheels cannot be
cross-built; if you want those you need the machines, or the release workflow.

### Check before uploading

```bash
python -m twine check wheelhouse/*.whl
```

`twine check` validates the metadata and the long description render. It does
**not** check that the wheel works — do that separately:

```bash
ls wheelhouse/                                   # expect manylinux_2_28_x86_64
unzip -l wheelhouse/*.whl | grep -c lsdtt3.libs  # expect > 0: repair ran
```

If the filename says `linux_x86_64` rather than `manylinux_...`, the repair
step did not run. Stop: that wheel is not publishable.

### Upload

```bash
python -m twine upload wheelhouse/*.whl
```

`twine` prompts for the username (`__token__`) and the token. To avoid
re-entering it, put it in `~/.pypirc`:

```ini
[pypi]
username = __token__
password = pypi-AgEIcHlwaS5vcmc...
```

`chmod 600 ~/.pypirc`. **Never commit it** — it is a credential with upload
rights to the project.

### Which tag a wheel comes from

lsdtt3 has one version number. `python_bindings/pyproject.toml` carries the
same `X.Y.Z` as `project(lsdtt3 VERSION ...)` in `CMakeLists.txt` (checked by
`unit_version_consistency`; set both with `tools/set_version.py`), and a wheel
is built from the release tag `vX.Y.Z`, never from an untagged commit.
Versions up to 0.9.3 are on GitHub releases only, with one wheel per Python
version; PyPI starts at the first release with the single `cp312-abi3` wheel.

### A rebuilt wheel of the same source: `.postN`

PyPI never accepts the same file twice and never frees a version, so a broken
wheel (built against the wrong GDAL, missing a data file) cannot be replaced
in place. When the **source is unchanged** and only the build was wrong,
publish a post-release:

```bash
git switch -c release/0.7 v0.7.1      # or stay on main if it has not moved on
# ... commit the packaging fix ...
tools/set_version.py 0.7.1 --post 1   # pyproject.toml -> 0.7.1.post1, nothing else
git commit -am "Python wheel 0.7.1.post1"
git tag -a v0.7.1.post1 -m "lsdtt3 Python wheel 0.7.1.post1"
git push origin v0.7.1.post1          # release.yml ignores .post tags
tools/make_all_packages.sh --clean --tag v0.7.1.post1 --only wheel
```

`lsdtt3.__version__` and the bundled CLIs then report `0.7.1.post1`. A source
change, even a Python-only fix, is never a `.postN`: it is a new patch release
for everything (`0.7.1`), because PEP 440 reserves post-releases for changes
that do not alter the distributed software.

If the broken wheel is actively harmful, also **yank** `0.7.1` on PyPI (the
project's release page → Options → Yank). pip then skips it unless a user pins
it exactly.

### Release candidates

No pre-release versions are published, to PyPI or in backend bundles (the QGIS
plugin cannot parse a suffixed `backend_version`). Every release goes to
TestPyPI first, and a required reviewer on the `pypi` environment gives the
chance to try it there before PyPI. A version can also be uploaded to TestPyPI
alone with a manual run (`publish_to: testpypi`).

---

## After publishing

1. **Verify from a clean environment**, not the one you built in:

   ```bash
   python -m venv /tmp/verify && /tmp/verify/bin/pip install lsdtt3
   /tmp/verify/bin/python -c "import lsdtt3; print(lsdtt3.__version__)"
   /tmp/verify/bin/lsdtt-channel-extraction -v
   ```

   The second and third lines matter separately: the first proves the
   bindings and their vendored libraries load, the second proves the console
   scripts and the bundled CLI binaries work.

2. **Check the project page** renders the README.

---

## Things that bite

**A version can never be reused.** Deleting a release on PyPI does not free
the version number. Rehearse on TestPyPI.

**The Python package shares the lsdtt3 version number.** `lsdtt3` 0.7.1 on
PyPI is built from tag `v0.7.1`, and `lsdtt3._core.lsdtt3_version()` reports
the version the bundled C++ was built as. A user pinning a version that has no
wheel gets "no matching distribution"; the README's "Versioning" section is
there to explain that.

**`pip install` pulls numpy and geopandas.** They are real dependencies
(`[project] dependencies`), and geopandas is a heavy solve. A user in a
constrained environment will notice.

**Windows wheels are single-threaded.** The Windows leg builds with
`-DLSDTT_PYTHON_OPENMP=OFF`, because MSVC's OpenMP does not cover the
constructs lsdtt3 uses. Windows users get correct results at single-thread
speed; do not quote a Windows wheel timing as an lsdtt3 benchmark.

**There is no source distribution.** Only wheels are published, so a platform
without a wheel cannot `pip install` at all. If you want an sdist, it has to be
built and uploaded deliberately, and it will need the C++ toolchain and GDAL
headers on the user's machine to build.
