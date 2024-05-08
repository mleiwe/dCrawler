from setuptools import setup,find_packages
setup(
    name='dCrawler',
    version='0.1',
    author='Marcus Lewie,Biswanath Saha',
    author_email='mlewie@gmail.com',
    description='dCrawler is a fully independent clustering algorithm that only requires a distance threshold (Th(d)) to perform clustering.',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'pandas',
        'scipy',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.9',
)
