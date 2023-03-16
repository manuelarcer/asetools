from setuptools import setup

setup(
    name="asetools",
    version='1.0',
    zip_safe=True,
    py_modules=['analysis', 'doscar_analysis', 'databases'],
    scripts=['bin/foldersummary.py']
)
