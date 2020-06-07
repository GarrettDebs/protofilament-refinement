"""
Setup file for AIdiet.
"""
import setuptools

setuptools.setup(
    name="pf_refinement",
    version="0.1",
    author="Garrett Debs",
    author_email='gedebs37@gmail.com',
    description="A cryo-EM reconstruction technique to more precisely refine microtubule structures",
    packages=setuptools.find_packages(),
    scripts=[
        "commands/cf_init_project",
        "commands/cf_patch_masks",
        "commands/cf_relion_subtract",
        "commands/cf_protofilament_subtract",
        "commands/cf_focused_classification"
        ],
    package_data={'pf_refinement': ['data/mt50.jpg']}
)
