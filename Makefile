.PHONY: black

default:

black:
	black -l 100 clinvar_tsv tests

test:
	py.test --cov=clinvar_tsv tests
	coverage report
