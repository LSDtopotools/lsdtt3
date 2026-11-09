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

### 2. Check the licence obligations

Uploading publishes compiled code and the CLI executables to anyone who wants
them. That is compatible with AGPL-3.0-or-later, which also requires the
corresponding source to be available; the public repository and its release
tags provide it. Each wheel bundles GDAL, PROJ and their dependencies, whose
licences travel with them.

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

Use this when you want to control the upload yourself, or to rehearse on
TestPyPI without touching CI.

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

### Which tag a wheel comes from

lsdtt3 has one version number. `python_bindings/pyproject.toml` carries the
same `X.Y.Z` as `project(lsdtt3 VERSION ...)` in `CMakeLists.txt` (checked by
`unit_version_consistency`; set both with `tools/set_version.py`), and a wheel
is built from the release tag `vX.Y.Z`, never from an untagged commit:

```bash
tools/make_all_packages.sh --clean --tag v0.7.1 --only wheel   # needs podman for manylinux
gh release upload v0.7.1 dist/v0.7.1/*.whl                     # attach to the lsdtt3 release
```

There is no separate Python tag series. Before the version numbers were
unified the wheel workflow would have tagged `python-vX.Y.Z`; none was ever
pushed, and the workflow no longer creates tags.

**Not every release gets a wheel.** Wheels are built only for releases where
the Python API is supported; a version missing from PyPI is expected, and
`python_bindings/README.md` (the PyPI project page) says so in its
"Versioning" section. When a release has no wheel, its GitHub release notes say
"Python wheel: not built for this release".

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

### Release candidates: TestPyPI only

No pre-release versions are published to PyPI or in backend bundles (the QGIS
plugin cannot parse a suffixed `backend_version`). To try a wheel before the
real upload, edit `pyproject.toml` to `0.7.1rc1` locally, **do not commit or
tag it**, build, and upload to TestPyPI (step 5 above). That catches most
"uploaded a broken wheel, now need `.post1`" mistakes.

---

## Route B: publishing from GitHub Actions

`.github/workflows/python-wheels.yml` builds and tests wheels for one existing
tag (its `tag` input: `vX.Y.Z` or `vX.Y.Z.postN`), on Linux only or on all four
platform legs (its `platforms` input; see
[A Linux-only wheel build](building_python_wheels.md#a-linux-only-wheel-build)), and
can upload them to the `vX.Y.Z` GitHub release (its `upload` input).
`release.yml` calls it on every `vX.Y.Z` tag with all platforms and upload on;
it can also be run by hand. It never
bumps a version, commits or tags. **It does not upload to PyPI** — there is no
`twine` step and no token. Adding one is a small change.

### The preferred mechanism: trusted publishing

PyPI's trusted publishing (OIDC) means **no token is stored anywhere**. GitHub
proves the workflow's identity to PyPI directly. Configure it once at
<https://pypi.org/manage/project/lsdtt3/settings/publishing/>, naming the
repository, the workflow file and the environment (`pypi`). PyPI matches the
**top-level** workflow of the run, and `python-wheels.yml` runs both as a
reusable workflow called by `release.yml` (tag releases) and on its own
(manual `.postN` runs), so register both workflow files as trusted publishers.

Then add a job after the existing release job:

```yaml
  publish-pypi:
    name: Publish to PyPI
    needs: [build-wheels, upload-release]
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

The wheel workflow is the slow one — four platform legs, five Python
versions each, every wheel compiled from scratch — and it already runs on
every release tag. A publish job costs almost nothing on top: it downloads
artefacts already built and uploads them. Standard GitHub-hosted runners are
free for public repositories, macOS and Windows included.

---

## Choosing between the two

| | Route A (twine) | Route B (Actions) |
|---|---|---|
| CI | none | the wheels `release.yml` already built on the tag |
| Platforms | Linux only, unless you own Macs and Windows boxes | Linux only, or Linux, macOS ARM, macOS Intel, Windows |
| Python versions | whatever you tell cibuildwheel to build | cp310–cp314 |
| Credentials | a token on your machine | OIDC, nothing stored |
| Tag and GitHub release | you tag; `gh release upload` | you tag; the workflow can upload |
| Pre-upload testing | you run it | `test-command` runs against every wheel |

**Route B for releases**: `release.yml` has already built and tested all
twenty wheels on the tag, so publishing them is one more job. With no source
distribution, publishing only Linux wheels would leave macOS and Windows users
with nothing to install. **Route A** to rehearse on TestPyPI, or for a
Linux-only `.postN` fix.

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
