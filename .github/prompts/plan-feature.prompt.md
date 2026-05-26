---
description: "Plan a new metacooler feature. Explores the codebase, identifies affected files, and produces a step-by-step implementation plan."
---
# Plan a new feature for metacooler

Given the feature description below, explore the existing codebase to understand:
1. Which existing modules are relevant
2. What patterns to follow (check similar implementations)
3. Which files need to be created or modified
4. What tests are needed

Then produce a detailed implementation plan with:
- New files to create (with their purpose)
- Existing files to modify (with specific changes)
- Test files to create
- Verification steps (`uv run pytest`, `uv run ruff check`)

## Feature: {{input}}

Read the project's architecture from `.github/copilot-instructions.md`, then explore the source tree at `src/metacooler/` and `tests/` before planning.
