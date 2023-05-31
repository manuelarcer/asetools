from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        name='asetools',
        version='0.1.0',
        description='Complement to ASE to perform analysis and manipulate atoms',
        author='Juan M Arce-Ramos',
        author_email='manuel.arceram@gmail.com',
        url='https://github.com/manuelarcer/asetools.git',
        packages=find_packages(),
        install_requires=[
        # List of packages your project depends on
        entry_points={
            'console_scripts': [
                'executable_script1=asetools.bin.getenergy:main',
                'executable_script2=asetools.bin.summaryfolders:main',
                # etc.
            ],
        }
        ...
    )
