try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'alignSeq',
    'author': 'wchen1994',
    'version': '0.1',
    'install_requires': ['matplotlib'],
    'packages': ['alignSeq'],
}

setup(**config)
