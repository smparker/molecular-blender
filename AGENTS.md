# Local Blender installation

After changing the add-on or building its extensions, copy the tested package
into the local Blender installation, including for Python-only changes. Rebuild
changed native extensions for Blender's bundled Python before installing.

For the current Blender 4.3 installation, run this from the repository root:

```sh
python3 blender_setup.py install --mode=copy --link-target=4.3
```

Verify the installed copy in a fresh Blender process started outside the
checkout, using the user's normal preferences.
