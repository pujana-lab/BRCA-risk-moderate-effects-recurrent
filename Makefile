.PHONY: test tree

test:
	pytest -q

tree:
	find . -maxdepth 3 -type f | sort
