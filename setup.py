from setuptools import setup

setup(
    name='py-microwave',
    version='1.0.1',
    packages=['py-microwave'],
    url='https://github.com/sfpeik/py-microwave',
    license='MIT',
    author='speik',
    author_email='speik@hs-bremen.de',
    description='Microwave Toolbox',
    install_requires=[
          'schemdraw',
          'numpy',
          'matplotlib',
          'scipy'
      ],
)
