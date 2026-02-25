from setuptools import setup

setup(
    name="tps_eval",
    packages=[
        "src",
        "vendor.cif_to_pdb",
        "vendor.pymol_scripts",
    ],
)
