from setuptools import setup

def _LoadFile(filename):
    with open(filename) as infile:
        return infile.read()

setup(
    name='astronomy-engine',
    version='2.1.12',
    description='Astronomy calculation for Sun, Moon, and planets.',
    long_description=_LoadFile('README.md'),
    long_description_content_type='text/markdown',
    author='Donald Cross',
    author_email='cosinekitty@gmail.com',
    url='https://github.com/cosinekitty/astronomy',
    license=_LoadFile('LICENSE'),
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7"
    ],
    packages = ['astronomy']
)
