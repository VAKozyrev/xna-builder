import os
import json
import tempfile
import contextlib
from pathlib import Path
from typing import Iterator


def read_config(config_path: Path) -> list[dict]:
    """Reads a JSON config file into a list of dictionaries.

    Args:
        input_path: The path to the input JSON file.

    Returns:
        A list of configuration dictionaries.
    """

    with config_path.open("r") as f:
        config = json.load(f)

    if isinstance(config, dict):
        config_list = [config]
    elif isinstance(config, list):
        config_list = config

    return config_list


@contextlib.contextmanager
def change_dir(new_dir: Path | str) -> Iterator[None]:
    """
    This is a context manager that changes the CWD to a new directory upon
    entering the `with` block and returns to the original upon exiting.

    Args:
        new_dir: The path to the directory to change into.

    Yields:
        None: This context manager does not yield a value.

    Example:
        >>> print(Path.cwd())
        /home/user/project
        >>> with change_dir(".."):
        ...     print(Path.cwd())
        /home/user
        >>> print(Path.cwd())
        /home/user/project
    """
    original_dir = Path.cwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(original_dir)


@contextlib.contextmanager
def in_temp_dir() -> Iterator[Path]:
    """
    A context manager that creates a temporary directory, changes the
    current working directory into it, and cleans up upon exit.

    Yields:
        Path: The path to the newly created temporary directory.
    """
    with tempfile.TemporaryDirectory() as temp_dir_str:
        temp_dir_path = Path(temp_dir_str)
        original_dir = Path.cwd()
        try:
            os.chdir(temp_dir_path)
            yield temp_dir_path
        finally:
            os.chdir(original_dir)