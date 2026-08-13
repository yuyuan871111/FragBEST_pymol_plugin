# Development and installation checks

This guide verifies that the plugin and its Python dependencies are available
to the same PyMOL installation a user will run.

## 1. Use PyMOL's Python environment

Run commands with the PyMOL executable, not a potentially different system
Python. On this machine the executable is:

```bash
pymol
```

Confirm the version and Python package location:

```bash
pymol -cq -d "import pymol; print(pymol.__file__); print(cmd.get_version())"
```

## 2. Verify the plugin imports

From the repository root, run:

```bash
pymol -cq -d "import sys; sys.path.insert(0, '.'); import FragBEST_pymol_plugin; from FragBEST_pymol_plugin.loadPLY import load_ply; print('Plugin import passed')"
```

`Plugin import passed` confirms that PyMOL can find the plugin and its required
Python packages, including NumPy.

## 3. Run the automated reader tests

The tests cover both ASCII and binary PLY input and use Python's built-in
`unittest`, so no separate test package needs to be installed. From the
repository root, run:

```bash
pymol -cq -d "import sys, unittest; sys.path.insert(0, '.'); suite = unittest.defaultTestLoader.discover('tests'); result = unittest.TextTestRunner(verbosity=2).run(suite); sys.stdout.write('All tests passed.\n' if result.wasSuccessful() else 'Tests failed.\n'); sys.stdout.flush(); cmd.quit(0 if result.wasSuccessful() else 1)"
```

## 4. Test the development checkout in an open PyMOL session

You do not need to build or install a ZIP plugin while developing. Open PyMOL,
switch its console to Python mode, and run the following. Replace
`/path/to/FragBEST_pymol_plugin` with the repository directory (the directory
that contains the `FragBEST_pymol_plugin` folder).

```python
import sys

sys.path.insert(0, "/path/to/FragBEST_pymol_plugin")
from FragBEST_pymol_plugin import (
    load_grids,
    load_pdb_with_frags,
    load_ply,
    super_ply,
)
```

This imports the development versions of the plugin functions into the current
PyMOL session. Confirm the loader works with a known small PLY file:

```python
load_ply("sample.ply", interest_pt=1, ignore_surface=1)
```

This should create the feature objects and a `mesh_sample.ply` object.
`ignore_surface=1` makes the check faster and isolates vertex loading. Restart
PyMOL after changing code, then repeat the import steps to test the updated
checkout.

If you specifically want to test the PyMOL console command aliases, register
them after the imports:

```python
from FragBEST_pymol_plugin import __init_plugin__

__init_plugin__(None)  # Enables loadply, loadfrag, loadpredict, and superply.
```

## 5. Zip as a pymol plugin
```bash
mkdir pkgs   
cp -r FragBEST_pymol_plugin pkgs/.   
cp LICENSE pkgs/FragBEST_pymol_plugin/.  
cp NOTICE pkgs/FragBEST_pymol_plugin/.  
cp README.md pkgs/FragBEST_pymol_plugin/.
rm -rf pkgs/FragBEST_pymol_plugin/__pycache__
zip -r FragBEST_pymol_plugin.zip pkgs/FragBEST_pymol_plugin
rm -rf pkgs
```


## Common failures

- **`ModuleNotFoundError: FragBEST_pymol_plugin`**: run from the repository
  root for development, or reinstall the plugin and restart PyMOL.
- **PLY parsing error**: confirm the file is a valid ASCII or standard binary
  PLY mesh with `x`, `y`, `z` vertex properties and a face index list.
