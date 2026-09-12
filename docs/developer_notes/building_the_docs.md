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

Both documentation tasks put `python_bindings/src` on `PYTHONPATH`, then serve
or build MkDocs. This lets the Python API reference use the source-only facade
without building the native extension or installing an editable wheel.

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

### Option 1 — copy `site/` to the Edinburgh web space (works today)

This is the route that matches how the release packages are already
distributed, and it does not care that the repository is private.

```bash
pixi run docs-build
rsync -avz --delete site/ smudd@<host>:~/public_html/lsdtt3_docs/
```

`--delete` removes files on the server that no longer exist locally, so renamed
or deleted pages do not linger as stale URLs. Check the target path before the
first run — `--delete` against the wrong directory is destructive.

The site is then at roughly
`https://www.geos.ed.ac.uk/~smudd/lsdtt3_docs/`, alongside
`~smudd/lsdtt_packages/` which serves the tarballs.

### Option 2 — GitHub Pages

```bash
pixi run -e docs mkdocs gh-deploy
```

This builds the site and **force-pushes** it to a `gh-pages` branch, which
GitHub Pages then serves at `https://lsdtopotools.github.io/lsdtt3/`.

Two things to check before relying on it:

1. **The repository is private.** GitHub Pages from a private repository
   requires a paid plan (Pro, Team, or Enterprise). On a Free plan, Pages only
   publishes from public repositories, so `gh-deploy` will push the branch and
   nothing will be served.
2. **A published Pages site is public**, even when the code is not (private
   Pages sites need Enterprise). That is probably what you want during a
   private beta — documentation visible, source not — but it does mean the
   install pages and their download URLs become publicly indexable. Decide
   that deliberately rather than discovering it.

`gh-deploy` force-pushes, so never keep anything by hand on `gh-pages`.

### Should it be automated?

Not yet, and not for free. CI minutes are metered on this private repository and
the release workflow is deliberately restricted to version tags for that reason
(see [`making_release_packages.md`](making_release_packages.md)). A docs build
on every push to `main` would spend more minutes than the releases do.

If it becomes worth automating, the cheapest shape is a job on the **existing
release workflow** — the docs are rebuilt when a version ships, which is when
they change meaningfully — rather than a second workflow with its own triggers.

---

## Checklist when adding a page

1. Put it in the right directory: `docs/` top level and `docs/cli_reference/`
   are user-facing; `docs/developer_notes/` is not.
2. Add it to the `nav:` block in `mkdocs.yml`. A page absent from the nav is
   reachable only by search or a direct link, and `--strict` will warn.
3. Run `pixi run docs-build` and fix anything it reports.

---

## See also

- [`building_from_source.md`](building_from_source.md) — building the software itself
- [`making_release_packages.md`](making_release_packages.md) — how releases are published
