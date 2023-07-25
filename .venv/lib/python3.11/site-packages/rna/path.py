"""
Input output utilities
TODO: look at https://github.com/marian-code/ssh-utilities !!!
"""
import os
import shutil
import functools
import subprocess
import contextlib
import logging

try:
    from pathlib import Path as PathlibPath

    PathlibPath().expanduser()  # not existing in python2 version of pathlib
except (ImportError, AttributeError):
    from pathlib2 import Path as PathlibPath
import glob


class Path(type(PathlibPath())):  # pylint: disable=too-few-public-methods
    """
    pathlib inherited Path, improved by some methods
    """

    def resolve(self):
        """
        TODO: look at this https://bugs.python.org/issue19776
            There are good reasons that Pathlib doesnt do that
        """
        new_path = self.expanduser()
        return super(Path, new_path).resolve()


def resolve(path, *joins):
    """
    Returns:
        path with resolved ~, symbolic links and relative paths.
        *joins: parts to join with path

    Examples:
        Resolves relative paths
        >>> import os
        >>> from rna.path import resolve
        >>> resolve("../rna") == resolve(".")
        True

        Resolves user variables
        >>> resolve("~") == os.path.expanduser("~")
        True

        Also resolves symlincs which i will not test.
    """
    if joins:
        for part in joins:
            path = os.path.join(path, part)
    # py27 pathlib.Path has no expanduser
    # return str(Path(path).expanduser().resolve())
    # resolve environment variable like ~
    path = os.path.expanduser(str(path))
    # make it absolute
    path = os.path.abspath(path)
    # resolve symlinks but preserve '//' in the beginning
    is_nas_path = path.startswith("//")
    path = os.path.realpath(path)
    if is_nas_path:
        path = "/" + path

    return path


def extension(path: str) -> str:
    """
    Returns:
        Extension of the path
    """
    return Path(path).suffix.lstrip(".")


def cp(source, dest, overwrite=True):  # pylint: disable=invalid-name
    """
    copy with shutil
    """
    source = resolve(source)
    dest = resolve(dest)
    if os.path.isfile(dest):
        if not overwrite:
            raise ValueError(
                "Attempt to overwrite destination path" " {0}.".format(dest)
            )
    if os.path.isfile(source):
        shutil.copy(source, dest)
    elif os.path.isdir(source):
        shutil.copytree(source, dest)
    else:
        raise TypeError("Source {source} is not file or dir".format(**locals()))


def scp(source, dest):
    """
    ssh copy
    """
    # subprocess.check_call does not execute in bash, so it does not know ~
    # Carefull when giving remote host : ~ replacement is not tested
    subprocess.check_call(["scp", source, dest])


def mv(source, dest, overwrite=True):  # pylint: disable=invalid-name
    """
    Move files or whole folders
    """
    source = resolve(source)
    dest = resolve(dest)
    if os.path.isfile(dest) and not overwrite:
        raise ValueError(
            "Attempting to move to already existing destination"
            " path {0}.".format(dest)
        )
    logging.debug("Moving file or dir %s to %s", source, dest)
    if os.path.isfile(source):
        shutil.move(source, dest)
    elif os.path.isdir(source):
        shutil.move(source, dest)
        # shutil.copytree(source, dest)
    else:
        raise TypeError("source is not file or dir")


def rm(path, *joins, recursive=False):  # pylint: disable=invalid-name
    """
    Remove path or (empty) directory
    """
    path = resolve(path, *joins)
    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        if recursive:
            shutil.rmtree(path)
        else:
            if len(ls(path)) > 0:
                raise ValueError("Directory {0} is not empty".format(path))
            os.rmdir(path)
    else:
        raise ValueError("Path {0} is not directory and not file".format(path))


def cd(directory):  # pylint: disable=invalid-name
    """
    Change directory
    """
    os.chdir(resolve(directory))


@contextlib.contextmanager
def cd_tmp(tmp_path):
    """
    Temporarily change directory. change back afterwards
    """
    cwd = os.getcwd()
    cd(resolve(tmp_path))
    try:
        yield
    finally:
        cd(cwd)


def ls(directory, *joins, **kwargs):  # pylint: disable=invalid-name
    """
    a bash ls imitation.
    allows wildcards
    """
    directory = resolve(directory, *joins)
    if os.path.isdir(directory):
        directory = os.path.join(directory, "*")
    return [os.path.basename(p) for p in glob.iglob(directory, **kwargs)]


def mkdir(path, is_dir=False):
    """
    Like mkdir -p

    Notes:
        if you want to create a dir, add '/' behind path or use is_dir=True

    Args:
        path (str): directory or path to create.
        is_dir (str): tell mkdir that you are explicitly giving a directory.
    """
    path = resolve(path)
    if is_dir:
        directory = path
    else:
        directory = os.path.dirname(path)
    if not os.path.exists(directory) and not directory == "":
        logging.debug("Create directory %s", directory)
        os.makedirs(directory)
    else:
        logging.debug("Directory %s already exists", directory)


def tree(directory):
    """
    Creates a nested dictionary that represents the folder structure of
    directory

    Examples:
        >>> import rna
        >>> import pathlib
        >>> this = pathlib.Path(__file__)
        >>> directory = this.parent
        >>> this.name in rna.path.tree(directory)
        True
    """
    dir_structure = {}
    directory = resolve(directory)
    directory = directory.rstrip(os.sep)
    start = directory.rfind(os.sep) + 1
    for path, dirs, files in os.walk(directory):  # pylint: disable=unused-variable
        folders = path[start:].split(os.sep)
        subdir = dict.fromkeys(files)
        parent = functools.reduce(dict.get, folders[:-1], dir_structure)
        parent[folders[-1]] = subdir
    if len(dir_structure) == 0:
        return {}
    return dir_structure[
        list(dir_structure.keys())[0]
    ]  # first entry is always directory


def exists(path, **kwargs):
    """
    Forwarding the resolved path to os.path.exists
    """
    return os.path.exists(resolve(path), **kwargs)
