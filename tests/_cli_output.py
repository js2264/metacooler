"""Helpers for asserting on rich-click's rendered CLI output.

rich-click does three things to an error message that a test never means to
assert on:

  - it draws a Unicode box around it,
  - it hard-wraps the text at the terminal width (mid-token, for a long
    ``tmp_path``), and
  - it interleaves ANSI escape sequences whenever rich decides the stream is
    colour-capable.

The last one is environment-dependent: rich honours ``FORCE_COLOR``, which
GitHub Actions sets, so the same assertion that passes on a developer's
machine sees ``Invalid value for '\x1b[1;36m--level\x1b[0m'`` on CI.  Strip
all three here so message assertions depend only on the message.
"""

from __future__ import annotations

import re

# CSI sequences (colours, and the cursor moves rich's progress bar emits) plus
# OSC sequences (hyperlinks), which are terminated by BEL or ST.
_ANSI_RE = re.compile(
    r"\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07\x1b]*(?:\x07|\x1b\\)|[@-Z\\-_])"
)

_BOX_CHARS = "│╭╮╰╯─"


def flat(output: str) -> str:
    """Collapse rich-click's boxed, wrapped, possibly coloured output.

    Returns the message on a single line with single spaces, so an assertion
    does not depend on where the wrap landed.
    """
    plain = _ANSI_RE.sub("", output)
    stripped = plain.translate({ord(c): " " for c in _BOX_CHARS})
    return " ".join(stripped.split())


def squash(output: str) -> str:
    """Like :func:`flat` but with *all* whitespace removed.

    rich-click wraps mid-token when a long path does not fit the line, so
    assertions that embed a path must ignore whitespace entirely.
    """
    return "".join(flat(output).split())
