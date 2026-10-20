PYTHON ?= .venv/bin/python
export PYTHONPATH := src
export MPLCONFIGDIR := $(CURDIR)/.cache/matplotlib

.PHONY: verify reproduce test
verify: test reproduce
reproduce:
	$(PYTHON) -m precovery.reproduce
test:
	$(PYTHON) -m unittest discover -s tests -v
