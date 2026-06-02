import re
from pathlib import Path
from setuptools import setup, find_packages

# Read version from surfcluster/__init__.py
init_text = Path('surfcluster/__init__.py').read_text()
match = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", init_text, re.M)
if not match:
    raise RuntimeError("Cannot find __version__ in surfcluster/__init__.py")
version = match.group(1)

setup(
    name='surfcluster',
    version=version,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'surfcluster=surfcluster.cli:main',
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'networkx',
        'pyyaml',
        'pandas',
    ],
    python_requires='>=3.9',
)
