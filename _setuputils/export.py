import os
import tempfile
import zipfile
from . generic import *
from . addon_files import iterRelativeAddonFiles, iterRelativeExportCFiles

currentDirectory = os.path.dirname(os.path.abspath(__file__))

def execute_Export(addonDirectory, exportPath, addonName):
    # Replace an existing ZIP only after the new central directory is written.
    destination = os.path.abspath(exportPath)
    files = [name for name in iterRelativeAddonFiles(addonDirectory)
             if os.path.abspath(os.path.join(addonDirectory, name)) != destination]
    with tempfile.TemporaryDirectory(dir=os.path.dirname(destination)) as staging:
        archive = os.path.join(staging, 'addon.zip')
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipFile:
            for relativePath in files:
                absolutePath = os.path.join(addonDirectory, relativePath)
                zipFile.write(absolutePath, os.path.join(addonName, relativePath))
        os.replace(archive, exportPath)

    print(f"Exported add-on to {exportPath}")

def execute_ExportC(addonDirectory, exportCPath, exportCSetupPath, addonName):
    removeFile(exportCPath)
    topdir = "animation_nodes_c"

    with zipfile.ZipFile(exportCPath, "w", zipfile.ZIP_DEFLATED) as zipFile:
        for relativePath in iterRelativeExportCFiles(addonDirectory):
            absolutePath = os.path.join(addonDirectory, relativePath)
            zipFile.write(absolutePath, os.path.join(topdir, addonName, relativePath))

        for path in iterPathsWithExtension(currentDirectory, ".py"):
            relativePath = os.path.relpath(path, currentDirectory)
            zipFile.write(path, os.path.join(topdir, "_setuputils", relativePath))

        setupInfo = getPlatformSummary()
        tmpFileName = "tmp.json"
        writeJsonFile(tmpFileName, setupInfo)
        zipFile.write(tmpFileName, os.path.join(topdir, "setup_info.json"))
        removeFile(tmpFileName)

        zipFile.write(exportCSetupPath, os.path.join(topdir, "setup.py"))

    print("Exported C Build:")
    print("    " + exportCPath)
