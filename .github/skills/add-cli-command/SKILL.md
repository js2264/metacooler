---
name: add-cli-command
description: "Add a new CLI command to metacooler. Use when: creating a new subcommand, adding a command to the CLI, implementing a new feature as a CLI command. Covers file creation, registration, command groups, and test scaffolding."
---
# Add a CLI Command to metacooler

## When to Use
- Adding a new `metacooler <command>` subcommand
- Implementing a feature that needs CLI exposure

## Procedure

### 1. Determine the category
Commands belong to one of four groups (defined in `cli/main.py`):
- **Create and edit**: `create`, `merge`, `binnify`, `zoomify`
- **Extract and view**: `extract`, `flatten`, `dump`, `network`, `aggregate`, `show`
- **Inspect**: `info`, `tree`, `metadata`, `attrs`
- **Manage collections**: `ls`, `cp`, `mv`, `rm`

### 2. Create the command file
Create `src/metacooler/cli/cmd_<name>.py`:

```python
"""metacooler <name> command."""
from __future__ import annotations
from pathlib import Path
import rich_click as click
from metacooler.cli.main import cli

@cli.command(name="<name>", help="Short help text.")
@click.argument("metacool_path", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path), required=True, help="Output path.")
@click.option("-v", "--verbose", is_flag=True, help="Increase verbosity (debug logging).")
def cmd_<name>(metacool_path: Path, output: Path, verbose: bool) -> None:
    """Detailed help text."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.<module> import <function>  # heavy imports inside
    setup_logging(verbose, command_name="metacooler <name>")
    <function>(uri=metacool_path, output=output)
```

### 3. Register in main.py
Add `cmd_<name>` to the import block and the appropriate `COMMAND_GROUPS` entry.

### 4. Create the implementation
- I/O operations → `src/metacooler/io/<name>.py`
- Computational operations → `src/metacooler/ops/<name>.py`

### 5. Write tests
Create `tests/test_<name>.py` using fixtures from `conftest.py`:
- `tiny_metacooler` for single-file tests
- `tiny_mmetacool` for collection tests
- Use `click.testing.CliRunner` for CLI integration tests

### 6. Verify
```bash
uv run pytest tests/test_<name>.py -x -v
uv run ruff check src/metacooler/cli/cmd_<name>.py src/metacooler/io/<name>.py tests/test_<name>.py
uv run ruff format src/metacooler/cli/cmd_<name>.py src/metacooler/io/<name>.py tests/test_<name>.py
uv run pytest tests/ -x -q  # full regression
```
