from distutils.core import setup

from setuptools import find_packages

from prince import __version__

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: MIT License
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3.6
Operating System :: POSIX :: Linux
""".strip().split('\n')

setup(name='prince',
      version=__version__,
      description='PRINCE estimates Variable Number Tandem Repeats (VNTR) copy number from raw next generation sequencing (NGS) data.',
      author='Julius Booth, Margaryta Vityaz, Merhdad Mansouri, Leonid Chindelevitch',
      author_email='',
      url='https://github.com/WGS-TB/PythonPRINCE',
      license='MIT',
      classifiers=classifiers,
      install_requires=[
          'biopython',
          'scipy',
          'numpy'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      packages=find_packages(),
      include_package_data=True,
      scripts=['bin/prince']
)
