from setuptools import setup

setup(
    name='pdbpython',
    version='0.1',
    author='Christian Feldmann',
    license="GNU GPLv3",
    packages=['pdbpython', ],
    author_email='christian.w.feldmann@gmail.com',
    description='Package to handle PDB files',
    install_requires=['setuptools', 'urllib3', 'numpy']
)
