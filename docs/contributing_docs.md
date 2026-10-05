# Edit, preview and publish documentation

The site uses Sphinx, the Read the Docs theme and MyST Markdown in the existing
`docs/` folder. The original logo is reused. Package APIs are extracted with
Sphinx AutoAPI; the README is included from its source, and notebooks are offered
as downloads/links without execution.

## 1. Create a documentation environment

From the repository root:

```bash
python3 -m venv .venv-docs
source .venv-docs/bin/activate
python -m pip install -r docs/requirements.txt
```

On Windows, activation is `.venv-docs\Scripts\activate`. Python 3.11 or newer is
required for this documentation toolchain. This environment is independent of
your MD environment: it does not install or import the simulation package.

## 2. Edit a page

| What you want to change | File |
| --- | --- |
| Homepage/navigation | `docs/index.md` |
| Installation or first example | `docs/installation.md`, `docs/quickstart.md` |
| Workflow tutorial | `docs/workflows/*.md` |
| Current support / limits | `docs/capabilities.md` |
| Notebook downloads | `docs/notebooks.md` |
| API descriptions | Docstrings in `src/iphasimulator/`; generated pages are not edited |
| Theme, extensions or API settings | `docs/conf.py` |
| Small visual adjustments | `docs/_static/custom.css` |
| Build/deployment checks | `.github/workflows/docs.yml` |

Write ordinary Markdown. Add a new page to the appropriate `toctree` in
`index.md`, and use relative links such as `[installation](installation.md)`.
For source-file downloads, use MyST's download role:

```text
{download}`Notebook <../notebooks/01_build_pha_oligomer.ipynb>`
```

For API links, use a Python-domain role, for example
`` {py:func}`iphasimulator.naming.oligomer_name` ``.
Use a normal code fence for examples; no code fence is executed by this build.

## 3. Build and check

```bash
python -m sphinx -b html -W --keep-going docs docs/_build/html
python docs/check_links.py docs/_build/html
```

`-W` makes Sphinx warnings fail the build; `--keep-going` reports all issues in
one pass. The second command validates local links, fragments, images, scripts,
stylesheets and notebook downloads. It is also run on pull requests.

For a completely fresh build after removing/renaming pages:

```bash
python -m sphinx -M clean docs docs/_build
python -m sphinx -b html -E -a -W --keep-going docs docs/_build/html
python docs/check_links.py docs/_build/html
```

Generated `docs/_build/`, `docs/autoapi/` and `.venv-docs/` are ignored by Git.
If you remove/rename a Python module, delete its corresponding generated AutoAPI
page as well before rebuilding (or remove the ignored `docs/autoapi/` directory).
There is no notebook execution extension or simulation hook in this configuration.
Two existing docstrings have formatting defects; a narrowly scoped hook in
`docs/conf.py` formats them for Sphinx without changing the Python sources.
Generated API pages offer their local page source instead of an invalid GitHub
edit link. Edit the package docstrings when maintaining the API documentation.

To inspect remote links when network access is available:

```bash
python -m sphinx -b linkcheck docs docs/_build/linkcheck
```

Remote services may reject automated requests. Review that report separately;
pull requests use deterministic local-link checks and do not depend on remote
availability. Neither kind of documentation check validates an MD protocol.

## 4. Preview locally

```bash
python -m http.server 8000 --bind 127.0.0.1 --directory docs/_build/html
```

Open **http://127.0.0.1:8000/** in your browser. Rebuild after editing, then
refresh. Stop the server with `Ctrl+C`. This serves only on your computer and
does not publish the documentation.

## 5. Publish through the existing GitHub Pages workflow

Pages already uses GitHub Actions for this repository. Respect `main` protection:
submit a pull request and obtain its required approving review before merging.
Do not use the account’s bypass to push directly to `main`. After review:

1. Commit to a topic branch and merge the approved pull request to `main` in
   `MMLabCodes/iPHAsimulatorV2`.
2. In the repository, open **Settings → Pages → Build and deployment** and set
   **Source: GitHub Actions**.
3. Under **Actions**, run the **Documentation** workflow on `main` (or allow a
   subsequent push to trigger it).
4. Check the `build` and `deploy` jobs and use the URL shown by the `github-pages`
   deployment environment. The expected project URL is
   `https://mmlabcodes.github.io/iPHAsimulatorV2/`.

The workflow builds on pull requests with read-only repository permissions and
uploads a preview artifact. **Pull requests never deploy**, including fork PRs.
Deployment is a separate job allowed only on `main` for push/manual events,
after a successful build. Only that job receives `pages: write` and
`id-token: write`. You can protect the `github-pages` environment with a required
reviewer if desired. Adjust branch/repository settings if this is used in a fork.

The Pages setup follows GitHub's [custom workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).
AutoAPI's [static-source configuration](https://sphinx-autoapi.readthedocs.io/en/latest/reference/config.html)
explains API generation options.
