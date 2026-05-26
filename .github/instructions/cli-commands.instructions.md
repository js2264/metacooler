---
description: "Use when creating or modifying CLI commands for metacooler. Covers rich-click patterns, command registration, option conventions."
applyTo: "src/metacooler/cli/**"
---
# CLI Command Conventions

## File structure
Each command lives in `cli/cmd_<name>.py`. It imports `cli` from `main.py`:

```python
from metacooler.cli.main import cli

@cli.command(name="<name>", help="Short description.")
@click.argument(...)
@click.option("-v", "--verbose", is_flag=True, help="Increase verbosity (debug logging).")
def cmd_<name>(..., verbose: bool) -> None:
    from metacooler.cli.logging_utils import setup_logging
    setup_logging(verbose, command_name="metacooler <name>")
    # ... implementation (import heavy deps inside the function)
```

## Registration
Add the import to `cli/main.py` in the alphabetically sorted import block:
```python
from metacooler.cli import (  # noqa: E402, F401
    cmd_<name>,
    ...
)
```
And add the command name to the appropriate group in `COMMAND_GROUPS`.

## Option patterns
- **Multi-value comma-separated**: use `multiple=True` + `callback=_parse_comma_list`
- **File paths**: use `click.Path(exists=True, path_type=Path)` for inputs
- **Output paths**: use `click.Path(path_type=Path)` (no exists check)
- **Every command gets `-v/--verbose`**

## Imports
Heavy imports (`MetaCooler`, I/O modules) go **inside** the command function body, not at module level. This keeps `--help` fast.
