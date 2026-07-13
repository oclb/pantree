---
name: pantree-setup
description: Install, validate, repair, or explain the Pantree development and runtime environment. Use for uv setup, Python compatibility, dependency problems, package installation, or initial command and test validation.
---

# Pantree setup

Use this skill for environment work. Use `pantree-usage` for GFA-to-VCF, simplification, consolidation, or format questions.

Install the locked environment from the repository root:

```bash
uv sync
```

Validate the package and command surface:

```bash
uv run python -c "import pantree"
uv run pantree --help
uv run pytest
```

Pantree supports Python 3.10 or newer. Prefer the checked-in `uv.lock`; do not replace the project environment with ad hoc `pip` installation unless requested.

When a fixture or dependency under Dropbox is zero bytes, treat it as potentially unsynced. Tell the user which file needs syncing rather than regenerating or overwriting it.

For a focused smoke check, run `uv run pytest tests/test_imports.py tests/test_cli.py`. Run the full suite before declaring an environment repair complete.
