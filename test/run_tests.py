#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test runner for molecular blender — discovers and runs all tests via unittest."""

import sys
import os
import unittest
from time import perf_counter

try:
    import bpy
except ModuleNotFoundError as exc:
    if exc.name != "bpy":
        raise
    bpy = None


def blender_extra_args(argv):
    """Return args after '--' when called through Blender, else args after script name."""
    if "--python" in argv:
        return argv[argv.index("--") + 1:] if "--" in argv else []
    return argv[1:]


class GroupedTestResult(unittest.TestResult):
    """Prints results grouped by file, like pytest's default output."""

    def __init__(self, stream):
        super().__init__()
        self.stream = stream
        self.started = perf_counter()
        self._current_module = None

    def _module_of(self, test):
        return test.__class__.__module__

    def _start_module(self, test):
        mod = self._module_of(test)
        if mod != self._current_module:
            if self._current_module is not None:
                self.stream.write("\n")
            self._current_module = mod
            label = f"test/{mod}.py"
            self.stream.write(f"{label} ")
            self.stream.flush()

    def startTest(self, test):
        super().startTest(test)
        self._start_module(test)
        if bpy is not None:
            # A fresh scene isolates frame, selection, render and molecule settings.
            self.original_scene = bpy.context.window.scene
            self.original_ids = set(bpy.data.user_map())
            bpy.context.window.scene = bpy.data.scenes.new("Test scene")

    def stopTest(self, test):
        try:
            if bpy is not None:
                bpy.context.window.scene = self.original_scene
                # Include unlinked meshes, actions, materials and linked libraries.
                added = set(bpy.data.user_map()) - self.original_ids
                bpy.data.batch_remove(added)
                if set(bpy.data.user_map()) != self.original_ids:
                    raise AssertionError("Blender datablocks changed across " + test.id())
        except Exception:
            self.addError(test, sys.exc_info())
        finally:
            super().stopTest(test)

    def addSuccess(self, test):
        self.stream.write(".")
        self.stream.flush()

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.stream.write("F")
        self.stream.flush()

    def addError(self, test, err):
        super().addError(test, err)
        self.stream.write("E")
        self.stream.flush()

    def addSkip(self, test, reason):
        super().addSkip(test, reason)
        self.stream.write("s")
        self.stream.flush()

    def printErrors(self):
        self.stream.write("\n")
        self._print_error_list("ERROR", self.errors)
        self._print_error_list("FAIL", self.failures)

    def _print_error_list(self, flavour, errors):
        for test, err in errors:
            self.stream.write(f"\n{'='*70}\n{flavour}: {test}\n{'-'*70}\n{err}\n")

    def printSummary(self):
        run = self.testsRun
        skipped = len(self.skipped)
        failures = len(self.failures)
        errors = len(self.errors)

        self.stream.write(f"{'-'*70}\n")
        self.stream.write(f"Ran {run} tests in {perf_counter() - self.started:.2f} s\n")
        if sys.platform in ("darwin", "linux"):
            import resource
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            scale = 1024 ** (2 if sys.platform == "darwin" else 1)
            self.stream.write(f"Peak process RSS {rss / scale:.1f} MiB\n")
        self.stream.write("\n")

        if self.wasSuccessful():
            status = "OK"
            if skipped:
                status += f" (skipped={skipped})"
        else:
            counts = []
            if failures:
                counts.append(f"failures={failures}")
            if errors:
                counts.append(f"errors={errors}")
            if skipped:
                counts.append(f"skipped={skipped}")
            status = f"FAILED ({', '.join(counts)})"
        self.stream.write(f"{status}\n")


if __name__ == "__main__":
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(repo_root)

    test_dir = os.path.join(repo_root, "test")
    # Blender does not add the script directory to sys.path. Prefer this
    # checkout over installed addons, and support named tests as well as discovery.
    sys.path[:0] = [repo_root, test_dir]
    extra = blender_extra_args(sys.argv)

    loader = unittest.TestLoader()

    if extra:
        suite = loader.loadTestsFromNames(extra)
    else:
        suite = loader.discover(start_dir=test_dir, pattern="test_*.py", top_level_dir=test_dir)

    stream = sys.stderr
    result = GroupedTestResult(stream)
    suite.run(result)
    result.printErrors()
    result.printSummary()
    sys.exit(0 if result.wasSuccessful() else 1)
