# Building and publishing the documentation site

The user-facing documentation in `docs/` is also the source for a static HTML
site, built with [MkDocs](https://www.mkdocs.org) and the
[Material](https://squidfunk.github.io/mkdocs-material/) theme.

| Thing | Where |
|---|---|
| Configuration, theme, navigation | `mkdocs.yml` (repo root) |
| Page sources | `docs/` |
| Extra styling | `docs/stylesheets/extra.css` |
| Build output | `site/` (gitignored — never commit it) |

---

## Preview it locally

```bash
pixi run docs-serve
```

Then open <http://127.0.0.1:8000>. The server watches `docs/` and `mkdocs.yml`
and reloads the page as you save, so this is the right tool while writing.
`Ctrl-C` stops it.

The first run solves and downloads the `docs` pixi environment (MkDocs plus the
Material theme), which takes a few minutes. After that it starts instantly.

## Build the HTML

```bash
pixi run docs-build
```

This writes the complete static site to `site/` — plain HTML, CSS and
JavaScript, no server-side anything. You can open `site/index.html` directly in
a browser to check it.

The Python API pages are generated from the committed runtime manifest before
release. MkDocs reads those Markdown pages directly, so serving or building the
site does not require the native extension.

`--strict` **turns warnings into errors**.
That is deliberate: without it, a link to a page you renamed or deleted still
"builds successfully" and simply 404s for the reader. If the build fails,
the message names the file and the bad link.

> MkDocs is not needed to build or run lsdtt3, so it is not in the default pixi
> environment. It lives in the `docs` feature, and the two `docs-*` tasks are
> attached to that feature so they select the right environment on their own.

---

## What gets published, and what does not

`mkdocs.yml` has an `exclude_docs` block:

```yaml
exclude_docs: |
  plan_docs/
  tex_technical/
  test_scripts/
  cli_helper_function_audit.md
  cli_parameter_consistency_audit.md
```

Those directories live under `docs/` for convenience but are **not** part of the
site. `plan_docs/` in particular holds uncommitted working notes that are not
intended for readers. If you add another such directory, add it here too —
otherwise MkDocs will happily build every stray Markdown file into the public
site as an orphan page.

`use_directory_urls: false` is set, so pages are emitted as `overview.html`
rather than `overview/index.html`. That is what makes the output work when it is
served from a plain directory, or opened as local files, rather than requiring a
web server that resolves directory indexes.

---

## Publishing

The build output is a plain static site, so it can be served from anywhere
that serves files. Two routes are set up.

### GitHub Pages

`.github/workflows/docs.yml` builds the site with `--strict` and deploys it to
<https://lsdtopotools.github.io/lsdtt3/> **on every release tag** (`vX.Y.Z`,
not the wheel-only `.postN` tags), and on manual dispatch from the Actions tab.
Deploying on releases keeps the site describing a version users can install,
rather than unreleased work on `main`. To publish a documentation-only fix
between releases, run it by hand on `main` (Actions → *Documentation site* →
**Run workflow**).

Both jobs are guarded by `if: !github.event.repository.private`, so nothing
deploys while the repository is private (Pages from a private repository needs
a paid plan).

When the repository goes public, do this once:

1. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
2. **Settings → Environments → `github-pages` → Deployment branches and
   tags:** add a **tag** rule `v*`. The environment allows only the default
   branch out of the box, and a deploy triggered by a tag is rejected without
   this rule.
3. Actions → *Documentation site* → **Run workflow**, to publish without waiting
   for the next release.
4. Check <https://lsdtopotools.github.io/lsdtt3/>.

Why the site works on Pages without changes:

- `site_url` in `mkdocs.yml` is set to the Pages URL, so canonical links and
  the sitemap are correct even though the site lives under `/lsdtt3/`.
- `use_directory_urls: false` emits `page.html` rather than `page/index.html`.
  All links are relative, so the site works under a sub-path, from the
  Edinburgh web space, or opened as local files.
- The deploy uses the Pages artifact actions, so GitHub does not run Jekyll
  over the output and no `.nojekyll` file is needed.
- The QGIS page is built from `qgis_plugin/README.md` with
  `pymdownx.snippets`, which reads files from the checkout at build time.

A deployed Pages site is **public**, including the install pages and the
download URLs they contain.

Do not also use `mkdocs gh-deploy`. It force-pushes a `gh-pages` branch, which
is a different Pages source from the one the workflow uses.

### Edinburgh web space (works today)

This route matches how the release packages are distributed, and works while
the repository is private.

```bash
pixi run docs-build
rsync -avz --delete site/ smudd@<host>:~/public_html/lsdtt3_docs/
```

`--delete` removes files on the server that no longer exist locally, so renamed
or deleted pages do not linger as stale URLs. Check the target path before the
first run — `--delete` against the wrong directory is destructive.

The site is then at roughly `https://www.geos.ed.ac.uk/~smudd/lsdtt3_docs/`.

---

## Checklist when adding a page

1. Put it in the right directory: `docs/what_is_lsdtt3/`,
   `docs/design_and_algorithms/` and `docs/cli_reference/` are for users;
   `docs/developer_notes/` is for contributors. Both are published.
2. Add it to the `nav:` block in `mkdocs.yml`. A page absent from the nav is
   reachable only by search or a direct link, and `--strict` will warn.
3. Run `pixi run docs-build` and fix anything it reports.

---

## See also

- [`building_from_source.md`](building_from_source.md) — building the software itself
- [`making_release_packages.md`](making_release_packages.md) — how releases are published
