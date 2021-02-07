clean:
	rm -rf build/* && rm -rf dist/* && rm -rf pdbeccdutils.egg-info && rm -rf htmlcov && rm .coverage

test-report:
	pytest --cov=pdbeccdutils --cov-report=html

package:
	pytest;\
	python setup.py bdist_wheel;\
	pip uninstall pdbeccdutils;\
	pip install -e .;\
	find dist -name "pdbeccdutils-*-py3*" -exec pip install {} \;

make upload:
	twine upload dist/* --verbose

make pypi: package upload clean
