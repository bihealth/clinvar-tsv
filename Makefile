
default:

.PHONY: black
black:
	black -l 100 clinvar_tsv

.PHONY: test
test:
	pytest .

.PHONY: isort
isort:
	isort --force-sort-within-sections -profile .
