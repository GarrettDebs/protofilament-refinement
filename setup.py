"""
Setup file for AIdiet.
"""
import setuptools

setuptools.setup(
    name="pf_refinement",
    version="0.1",
    author="Garrett Debs",
    author_email='gedebs37@gmail.com',
    description="A cryo-EM reconstruction technique to more precisely refine "\
    "microtubule structures",
    packages=setuptools.find_packages(),
    scripts=[
        "commands/pf_init_project",
        "commands/pf_wedge_masks",
        "commands/pf_microtubule_subtract",
        "commands/pf_protofilament_subtract",
        "commands/pf_focused_classification",
        "commands/pf_preprocess",
        "commands/pf_plot_distortions"
        ],
    package_data={'pf_refinement': ['data/mt50.jpg']}
)
