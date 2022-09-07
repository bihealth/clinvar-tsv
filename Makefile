
default:

.PHONY: black
black:
	black -l 100 .

.PHONY: test
test:
	pytest .

.PHONY: isort
isort:
	isort --force-sort-within-sections --profile black --line-length 100 .

.PHONY: flake8
flake8:
	flake8
