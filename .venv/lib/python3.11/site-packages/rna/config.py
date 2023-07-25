"""
A generic interface for configuration setups using different
 dialects (e.g. cfg, ini, toml, json, etc.).
It allows to provide multiple configuration file locations wich are searched
 in the order they are given for a configuration option.
Also there is the option for falling back a hirarchy level above
 the section requested for "config inheritance".
Dynamic options and links are possible as well.
"""
from modulefinder import Module
import types
import typing
from abc import abstractmethod
import os
import re
import pathlib
import importlib
import logging
import hydra
import omegaconf
import functools


import rna.path
from rna.pattern.strategy import Strategy


OPTION_NOT_FOUND = object()
logger = logging.getLogger(__name__)


class Dialect(Strategy):
    """
    Abstract Class to be derived by specific dialects of a config file (.ini, .cfg, ...)
    """

    def __init__(self):
        super().__init__()
        self._config_cache = {}

    def _get_config(self, config_path) -> "configparser.ConfigParser":
        """
        Retrieve a config object by path (cached)
        """
        if config_path not in self._config_cache:
            if os.path.exists(config_path):
                config = self._load_config(config_path)
                self._config_cache[config_path] = config
            else:
                self._config_cache[config_path] = None
        return self._config_cache[config_path]

    @abstractmethod
    def _load_config(self, config_path: str) -> dict:
        """
        Return the config object (allows itemization).
        """

    def get(self, config_path, *keys: str) -> any:
        """
        Retrieve the value corresponding to the keys from the config file referred by config_path.

        Returns:
            Any or OPTION_NOT_FOUND flag if option or section not found.
        """
        config = self._get_config(config_path)
        if config is None:
            return OPTION_NOT_FOUND
        depth = len(keys)
        for i in range(depth):
            key = keys[i]
            if key not in config:
                return OPTION_NOT_FOUND
            config = config[key]
        return config

    @staticmethod
    def prepare_keys(*keys: str) -> typing.Tuple[str]:
        """
        Prepare keys for get(). Usually splitting along "." and removing empty strings.

        Examples:
            >>> Dialect.prepare_keys("section.subsection", "option")
            ('section', 'subsection', 'option')
            >>> Dialect.prepare_keys("section", "subsection", "option.suboption")
            ('section', 'subsection', 'option', 'suboption')
        """
        return tuple(filter(None, re.split(r"\.", ".".join(keys))))


class TomlParser(Dialect):
    """
    Parser flavour for the '.toml' dialect.
    """

    @staticmethod
    def _load_config(config_path: str) -> dict:
        with open(config_path) as _file:
            import toml

            return toml.load(_file)


class ConfigParser(Dialect):
    """
    Parser flavour for the '.config' dialect.
    """

    def __init__(self, *args, **kwargs):
        # import only when needed
        global configparser
        import configparser

        super().__init__(*args, **kwargs)

    def _load_config(self, config_path: str) -> "configparser.ConfigParser":
        config = configparser.ConfigParser()  # pylint:disable=undefined-variable
        config.read(config_path)
        return config

    def get(self, config_path, *keys: str):
        """
        Abstract method, derived from :meth:`Dialect.get`:
        """
        config_path = self._get_config(config_path)
        try:
            val = self.get(config_path, *keys)
            return val
        # pylint:disable=undefined-variable
        except (configparser.NoOptionError, configparser.NoSectionError):
            return OPTION_NOT_FOUND


class YamlParser(Dialect):
    """
    Parser flavour for the '.yaml' dialect.
    """

    @staticmethod
    def _load_config(config_path: str) -> dict:
        with open(config_path) as _file:
            import yaml

            return yaml.safe_load(_file)


def api(
    module: typing.Union[Module, str],
    *merge_configs: str,
    config_dir: typing.Optional[str] = None,
    config_name: str = "config",
    version_base: typing.Optional[str] = "1.1",
    **compose_kwargs
) -> omegaconf.OmegaConf:
    """
    Provide an OmegaConfig object read from :param:`config_dir`/:param:`config_name`.yaml

    Args:
        module: Module or path to the __init__.py of the module
          (intended use is to call this from the module's __init__.py
          and pass __file__ as argument).
        *merge_configs: paths to configs to be merged in the order given
        config_dir: main directory of hydra config
        config_name: name of main config to read (yaml extended)
        version_base: base version of hydra config (e.g. 1.1)
        **compose_kwargs: forwarded to hydra.compose

    # Environment variables

    This function sets the following additional environment variables:

    * package.name - module name of package
    * package.path - directory of __init__.py

    * os.home - user home directory ("~")
    * os.user - user name used for login
    * os.sysname - operating system name (e.g. 'Linux')
    * os.nodename - name of machine on network (implementation-defined, e.g. 'jdoe_laptop')
    * os.release - operating system release (e.g. '5.15.0-67-generic')
    * os.version - operating system version (e.g. '#74~20.04.1-Ubuntu SMP Wed Feb 22 14:52:34 UTC 2023')
    * os.machine - hardware identifier (e.g. 'x86_64')
    * os.platform - The name of the operating system dependent module imported.

    # Resolvers

    omegaconf provides [builtin resolvers](https://omegaconf.readthedocs.io/en/latest/custom_resolvers.html)

    * oc.env
    * oc.create
    * oc.deprecated
    * oc.decode
    * oc.select
    * oc.dict

    This function registers a new OmegaConf resolver to be used in interpolation as

    * "${path:path/to/stuff}" -> rna.path.resolve

    Notes:
        * Schema valdiation is not out of box with hydra.ConfigStore: [see here](https://suneeta-mall.github.io/2022/03/15/hydra-pydantic-config-management-for-training-application.html#where-hydra-falls-short-)
            Consider using pydantic on top of hydra

            ..code-block:: python

                import rna.config
                from pydantic.dataclasses import dataclass
                from omegaconf import OmegaConf

                @dataclass
                class Employee:
                    name: str
                    age: int


                @dataclass
                class Company:
                    name: str
                    employees: List[Employee]


                cfg = rna.config.api(...)
                company = Company(**OmegaConf.to_container(cfg.company))

            I found that this practice is working nice but careful: It is discarding additional keywords from the config!

        * Have a look at [hydra_zen](_modules/hydra_zen/structured_configs/_implementations.html)

    """
    if isinstance(module, types.ModuleType):
        module_name = module.__name__
        module_path = module.__path__[0]
    elif os.path.basename(module) == "__init__.py":
        module_name = module.split(os.path.sep)[-2]
        module_path = os.path.dirname(module)
    elif isinstance(module, str):
        module_name = module
        module_path = os.path.dirname(importlib.util.find_spec(module_name).origin)
    else:
        raise NotImplementedError

    if config_dir is None:
        config_dir = os.path.join(module_path, "conf")

    omegaconf.OmegaConf.register_new_resolver("path", rna.path.resolve, replace=True)

    logger.info("Reading main config in %s with name %s", config_dir, config_name)
    hydra.initialize_config_dir(config_dir, version_base=version_base)
    cfg = hydra.compose(config_name=config_name, **compose_kwargs)
    for path in merge_configs:
        if os.path.exists(path):
            logger.info("Merging with config from %s", path)
            cfg = omegaconf.OmegaConf.merge(cfg, omegaconf.OmegaConf.load(path))

    # Add dynamic variables
    with omegaconf.open_dict(cfg):
        if not hasattr(cfg, "package"):
            cfg.package = omegaconf.OmegaConf.create()
        cfg.package.name = module_name
        cfg.package.path = module_path

        if not hasattr(cfg, "os"):
            cfg.os = omegaconf.OmegaConf.create()
        cfg.os.home = str(pathlib.Path.home())
        try:
            cfg.os.user = os.getlogin()
        except OSError:
            cfg.os.user = "<not found>"
        uname = os.uname()
        cfg.os.sysname = uname.sysname
        cfg.os.nodename = uname.nodename
        cfg.os.release = uname.release
        cfg.os.version = uname.version
        cfg.os.machine = uname.machine
        cfg.os.platform = os.name

    # free the globalHydra because you might want to include another api that uses this function
    # You would get 'ValueError: GlobalHydra is already initialized' in that case
    hydra.core.global_hydra.GlobalHydra.instance().clear()

    return cfg


def fallback(obj, *attrs):
    """
    Get attributes from object, throwing out attrs before the last from right to left
    in case no omegaconf.errors.ConfigAttributeError was thrown"""
    attrs = list(attrs)
    key = attrs.pop(-1)
    while True:
        try:
            return functools.reduce(getattr, [obj] + attrs + [key])
        except omegaconf.errors.ConfigAttributeError as e:
            if not attrs or e.key != key:
                raise e
            attrs.pop(-1)
