"""Tests for molecular_blender.util"""

import time
import pytest

from molecular_blender.util import Timer, stopwatch, unique_name, round_sigfigs, darken_color


class TestTimer:
    def test_tick_returns_elapsed(self):
        t = Timer()
        time.sleep(0.05)
        elapsed = t.tick()
        assert elapsed >= 0.04

    def test_tick_print(self, capsys):
        import os
        os.environ["MB_DEBUG"] = "1"
        try:
            from molecular_blender import util
            old_debug = util.DEBUG
            util.DEBUG = True
            t = Timer()
            time.sleep(0.01)
            elapsed = t.tick_print("test label")
            assert elapsed >= 0.005
            captured = capsys.readouterr()
            assert "test label" in captured.out
        finally:
            util.DEBUG = old_debug
            os.environ.pop("MB_DEBUG", None)

    def test_tick_print_silent(self, capsys):
        from molecular_blender import util
        old_debug = util.DEBUG
        util.DEBUG = False
        try:
            t = Timer()
            t.tick_print("silent")
            captured = capsys.readouterr()
            assert captured.out == ""
        finally:
            util.DEBUG = old_debug


class TestStopwatch:
    def test_decorator_returns_value(self):
        @stopwatch("test_routine", verbose=False)
        def add(a, b):
            return a + b
        assert add(2, 3) == 5

    def test_decorator_verbose(self, capsys):
        @stopwatch("verbose_routine", verbose=True)
        def noop():
            pass
        noop()
        captured = capsys.readouterr()
        assert "verbose_routine" in captured.out


class TestUniqueName:
    def test_no_conflict(self):
        assert unique_name("foo", ["bar", "baz"]) == "foo"

    def test_conflict(self):
        assert unique_name("foo", ["foo", "foo0"]) == "foo1"

    def test_conflict_chain(self):
        assert unique_name("a", ["a", "a0", "a1", "a2"]) == "a3"

    def test_starting_suffix(self):
        result = unique_name("x", ["x0", "x1"], starting_suffix=0)
        assert result == "x2"

    def test_starting_suffix_no_conflict(self):
        result = unique_name("x", ["y"], starting_suffix=0)
        assert result == "x0"


class TestRoundSigfigs:
    def test_zero(self):
        assert round_sigfigs(0, 4) == 0

    def test_positive(self):
        assert int(round_sigfigs(12345, 2)) == 12000

    def test_negative(self):
        assert int(round_sigfigs(-12345, 2)) == -12000

    def test_small_positive(self):
        assert abs(round_sigfigs(0.00098765, 2) - 0.00099) < 1e-6

    def test_decimal(self):
        assert abs(round_sigfigs(3.1415, 2) - 3.1) < 0.05


class TestDarkenColor:
    def test_rgb(self):
        result = darken_color((1.0, 0.5, 0.0), 0.1)
        assert len(result) == 3
        assert abs(result[0] - 0.9) < 1e-6
        assert abs(result[1] - 0.45) < 1e-6
        assert abs(result[2] - 0.0) < 1e-6

    def test_rgba_preserves_alpha(self):
        result = darken_color((1.0, 1.0, 1.0, 0.5), 0.2)
        assert len(result) == 4
        assert abs(result[0] - 0.8) < 1e-6
        assert abs(result[3] - 0.5) < 1e-6

    def test_zero_amount(self):
        result = darken_color((0.5, 0.5, 0.5), 0.0)
        assert abs(result[0] - 0.5) < 1e-6
