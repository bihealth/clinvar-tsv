.PHONY: black

default:

black:
	black -l 100 clinvar_tsv

test:
	py.test
