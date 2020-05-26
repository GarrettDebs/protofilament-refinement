"""
Setup file for AIdiet.
"""
from setuptools import setup

setup(
    name="pf_refinement",
    version="0.1",
    description="A cryo-EM reconstruction technique to more precisely refine microtubule structures",
    packages=["pf_refinement"],
    scripts=[
        "command/cf_init_project",
        "commands/cf_patch_masks"
        ]
)
