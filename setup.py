from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='polharmonic',
      version='0.1',
      description='Polarized microscope simulation tool.',
      long_description=readme(),
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      url='https://github.com/talonchandler/polharmonic',
      author='Talon Chandler',
      author_email='talonchandler@talonchandler.com',
      license='MIT',
      packages=['polharmonic'],
      install_requires=[
          'matplotlib',
          'numpy',
      ],
      zip_safe=False,
      test_suite='tests',
      )
