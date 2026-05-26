fix:
	uv run ty check src/ tests/
	uv run ruff format .
	uv run ruff check . --fix
	source .venv/bin/activate && quarto render README.qmd --to gfm --output README.md

regenerate:
	source .venv/bin/activate && \
	cd tests/data && \
	./generate-raw-files.py && \
	./raw-to-cool.sh && \
	cd ../../

cli-tests:
	cd tests/data && ./cli_tests.sh && cd ../../
