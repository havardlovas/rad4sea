import pathlib
import setuptools

with open('requirements.txt') as f:
    required = f.read().splitlines()

setuptools.setup(
    name='rad4sea',
    version='0.0.4',    
    description='A Python package for for radiometric transforms for hyperspectral seafloor/water column mapping',
    long_description=pathlib.Path("README.md").read_text(),
    long_description_content_type = "text/markdown",
    url='https://github.com/havardlovas/gref4hsi',
    author='Haavard Snefjellaa Loevaas',
    author_email='havard.s.lovas@ntnu.no',
    license='EUPL-1.2',
    install_requires=required,

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',  
        'Operating System :: Microsoft :: Windows :: Windows 10',        
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    python_requires='>3.7',
    packages=setuptools.find_packages(),
    include_package_data=True
)
