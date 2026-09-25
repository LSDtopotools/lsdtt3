# Publishing lsdtt3 to PyPI

How `pip install lsdtt3` gets to be true. Two routes are described: publishing
by hand with `twine` from your own machine, and publishing from GitHub Actions.
Both upload the same artefacts; they differ only in who holds the token and
where the wheels are built.

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
| The bindings extension | `lsdtt3/_core.cpython-*.so` |
| All 39 command-line tools | `lsdtt3/bin/lsdtt-*` |
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

## Before the first upload

These are one-time, and the first three are genuinely blocking.

### 1. Claim the name

Check whether `lsdtt3` is free:

```bash
python -m pip index versions lsdtt3     # or just open https://pypi.org/project/lsdtt3/
```

If the name is taken, decide on an alternative **before** anything else: the
distribution name is baked into `[project] name`, the install instruction in
`docs/python/index.md`, and every downstream mention.

### 2. Decide about publishing binaries from a private repository

The repository is private; PyPI is not. Uploading publishes compiled code and
the CLI executables to anyone who wants them. That is compatible with
AGPL-3.0-or-later, and it is the point of shipping a wheel — but it is a
deliberate disclosure decision, so make it deliberately. The source stays
private; the binaries do not.

### 3. Fill in the PyPI metadata

`python_bindings/pyproject.toml` currently carries only `name`, `version`,
`description`, `license`, `requires-python` and `dependencies`. That is enough
to build, but it produces a bare PyPI page. Add at least:

```toml
[project]
readme = "README.md"                      # python_bindings/README.md exists already
authors = [{ name = "Simon M. Mudd" }]
keywords = ["geomorphology", "topography", "DEM", "LiDAR", "terrain analysis"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: GIS",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows",
]

[project.urls]
Homepage = "https://lsdtopotools.github.io/lsdtt3/"
Documentation = "https://lsdtopotools.github.io/lsdtt3/python/"
Repository = "https://github.com/LSDtopotools/lsdtt3"
```

Without `readme` the project page is empty, which for a scientific package is
a real cost: it is the page people judge it from.

### 4. Get an API token

Create a **project-scoped** token at <https://pypi.org/manage/account/token/>
once the project exists, or an account-scoped token for the very first upload
and replace it with a project-scoped one afterwards. The username is literally
`__token__` and the password is the whole `pypi-...` string.

### 5. Practise on TestPyPI first

Do this for the first release of a new project. TestPyPI is a separate
instance with separate accounts and tokens.

```bash
python -m twine upload --repository testpypi wheelhouse/*.whl
```

then check the result installs cleanly, in a throwaway environment:

```bash
python -m venv /tmp/tp && /tmp/tp/bin/pip install \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    lsdtt3
/tmp/tp/bin/python -c "import lsdtt3; print(lsdtt3.__version__)"
/tmp/tp/bin/lsdtt-raster-preprocessing -v
```

The `--extra-index-url` matters: numpy and geopandas are not on TestPyPI, so
without it the install fails for reasons that have nothing to do with lsdtt3.

**A version number can never be reused on PyPI, even after deleting a
release.** That is why the rehearsal happens on TestPyPI.

---

## Route A: publishing by hand with twine

Use this when you want to control the upload yourself, or when GitHub Actions
minutes are the thing you are trying to save.

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
cross-built; if you want those you need the machines, or route B.

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

### Keep the GitHub release in step

Route A does not create a git tag or a GitHub release; route B does both. If
you publish by hand, tag and release by hand too, so the PyPI version and the
repository history agree:

```bash
# Use the SAME tag convention CI uses: python-v<version>
# (python-wheels.yml builds it as f"python-v{version}")
git tag -a python-v0.2.0 -m "lsdtt3 Python 0.2.0"
git push origin python-v0.2.0
gh release create python-v0.2.0 wheelhouse/*.whl \
    --title "lsdtt3 Python 0.2.0" --generate-notes
```

Note this is a *different* tag series from the `v<version>` tags that mark
lsdtt3 releases. `python-v0.2.0` and `v0.5.2` are independent.

---

## Route B: publishing from GitHub Actions

`.github/workflows/python-wheels.yml` already builds and tests wheels on four
platform legs and creates a GitHub release. **It does not upload to PyPI** —
there is no `twine` step and no token. Adding one is a small change.

### The preferred mechanism: trusted publishing

PyPI's trusted publishing (OIDC) means **no token is stored anywhere**. GitHub
proves the workflow's identity to PyPI directly. Configure it once at
<https://pypi.org/manage/project/lsdtt3/settings/publishing/>, naming the
repository, the workflow file (`python-wheels.yml`) and the environment
(`pypi`).

Then add a job after the existing release job:

```yaml
  publish-pypi:
    name: Publish to PyPI
    needs: [build-wheels, publish-release]
    runs-on: ubuntu-24.04
    environment: pypi          # gate approvals here if you want them
    permissions:
      id-token: write          # required for trusted publishing
    steps:
      - uses: actions/download-artifact@v5
        with:
          pattern: wheels-*
          path: dist
          merge-multiple: true
      - uses: pypa/gh-action-pypi-publish@release/v1
```

The `environment: pypi` line is worth keeping: a GitHub environment can
require manual approval, which makes the upload a deliberate act rather than a
consequence of pressing "Run workflow".

### The alternative: an API token in a secret

If trusted publishing is not an option, store the token as the repository
secret `PYPI_API_TOKEN` and pass it:

```yaml
      - uses: pypa/gh-action-pypi-publish@release/v1
        with:
          password: ${{ secrets.PYPI_API_TOKEN }}
```

Prefer trusted publishing. A long-lived token in a secret is a credential that
can leak; OIDC is scoped to one workflow in one repository.

### Cost

The wheel workflow is the expensive one — four platform legs, five Python
versions each, every wheel compiled from scratch. The publish job itself costs
almost nothing: it downloads artefacts already built and uploads them. So
adding PyPI publishing does **not** meaningfully change the Actions bill;
building the wheels does.

---

## Choosing between the two

| | Route A (twine) | Route B (Actions) |
|---|---|---|
| Actions minutes | none | the full wheel matrix |
| Platforms | Linux only, unless you own Macs and Windows boxes | Linux, macOS ARM, macOS Intel, Windows |
| Python versions | whatever you tell cibuildwheel to build | cp310–cp314 |
| Credentials | a token on your machine | OIDC, nothing stored |
| Tag and GitHub release | you do it | automatic |
| Pre-upload testing | you run it | `test-command` runs against every wheel |

**Route B for real releases**, because Linux-only wheels leave most users
building from source, and because the per-wheel test suite runs there. **Route
A when you want a Linux release without spending minutes**, or to rehearse on
TestPyPI.

---

## After publishing

1. **Verify from a clean environment**, not the one you built in:

   ```bash
   python -m venv /tmp/verify && /tmp/verify/bin/pip install lsdtt3
   /tmp/verify/bin/python -c "import lsdtt3; print(lsdtt3.__version__)"
   /tmp/verify/bin/lsdtt-channel-extraction -v
   ```

   The second and third lines matter separately: the first proves the bindings
   and their vendored libraries load, the second proves the console scripts and
   the bundled CLI binaries work.

2. **Check the project page** renders the README.

3. **Update `docs/python/index.md`** if the install instruction changed.

---

## Things that bite

**A version can never be reused.** Deleting a release on PyPI does not free
the version number. Rehearse on TestPyPI.

**The Python package version is independent of the lsdtt3 version.**
`python_bindings/pyproject.toml` carries its own, separate from
`project(lsdtt3 VERSION ...)` in the top-level `CMakeLists.txt`. Publishing
Python 0.2.0 says nothing about which lsdtt3 version it wraps — worth stating
in the README so users can tell.

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
