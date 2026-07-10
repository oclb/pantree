import importlib
import pkgutil

import pantree


def test_package_imports():
    importlib.import_module("pantree")


def test_all_package_modules_import():
    module_names = [
        module.name
        for module in pkgutil.iter_modules(pantree.__path__, f"{pantree.__name__}.")
    ]

    for module_name in module_names:
        importlib.import_module(module_name)
