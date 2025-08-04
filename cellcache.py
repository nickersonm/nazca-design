"""Simplify Hashing of Function Signatures with the inspect Module

cashing in memory is still done using a cellnames dict
Add alternative chaching to sqlite and/or redis
s
Time cases:
1- cell create + latency
3- sqlite + cell deserialize

"""
import time
start = time.time()
import os
import inspect
from typing import Callable
import functools
import json
import hashlib  # to guarantee the hash is consistent
import sqlite3  # to store the cache in a database
import nazca as nd
print(f"Import time: {time.time() - start:0.3f}")


# mapping hashme parameters to cellcache parameters
#cfg.hashme -> len(cfg.cellcache) > 0
#cfg.hash_id ->
#cfg.hash_cellnameparams = tuple(parrepr["parameter"] for parrepr in parreprs) -> probably drop, only used in properties dict in Cell for info.
#cfg.hash_params -> paramdict
#cfg.hash_basename -> cellname_base
#cfg.hash_paramsname -> cellname_base_params
#cfg.hash_name -> cellname_full
#cfg.hash_func_id -> func_id
# ? -> funcname
# ? -> paramnames : probably drop, can get from paramdict.keys()
# ? -> hash (long hash) -> add
# ? -> cellname_hash (short hash) -> not needed in Cell?

nd.cfg.cellcachemodel = []  # list of cellcachemodels (to replace hashme data in cfg)

sqlite = True  # Use sqlite for caching if turned on
if sqlite:
    # connect to the sqlite database or create it if it doesn't exist.
    conn  = sqlite3.connect('cellcache.db')
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS my_table
        (key TEXT PRIMARY KEY, value TEXT, timestamp INTEGER)''')

DEFAULT_CELLNAME = "NONAME"
HASH_LENGTH = 4  # number of hash characters to add to the cellname
HASH_IF_NOARGS = False  # add hash to cellname even if no arguments are given, it True
HASH_PREFIX = "$"  # prefix to add to the hash in the cellname

def _create_cellname(
    arguments: dict,
    basename: str,
    parreprs: list,
    hash_short: str,
    hash_prefix: str = HASH_PREFIX,
    addhash: bool = True,
    suffix: str = "",
) -> dict:
    """construct the nazca cellname using the function name and parreprs and hash

    Args:
        arguments (dict): dictionary with the function arguments
        basename (str): base name of the cell
        parreprs (list): list of parameter representations
        hash_short (str): short hash of the function signature
        hash_prefix (str): prefix to add to the hash, default is '$'

    Returns:
        dict: dictionary with the cellname parts
    """
    name_long = basename
    for parrepr in parreprs:
        argkey = parrepr["parameter"]
        parfunc = parrepr["func"]
        if argkey in arguments.keys():
            parvalue = arguments[argkey]
            if not callable(parfunc):
                try:
                    parstr = format(parvalue, parrepr["format"])
                except ValueError:
                    parstr = parvalue
            else:
                try:
                    parstr = parfunc(parvalue)  # call formatting function
                except ValueError:
                    parstr = "*FUNC_ERROR*"
            name_long += f"{parrepr['prefix']}{parstr}"
        else:
            raise ValueError(
                f"Parameter '{argkey}' not found in function parameters {list(arguments.keys())}."
            )

    if not HASH_IF_NOARGS and len(arguments) == 0:
        cellname_full = nd.cfg.gds_cellname_cleanup(f"{basename}{suffix}")  # base + suffix
    elif addhash:
        cellname_full = nd.cfg.gds_cellname_cleanup(f"{name_long}_{hash_prefix}{hash_short}{suffix}")  # base + params + hash
    else:
        cellname_full = nd.cfg.gds_cellname_cleanup(f"{name_long}{suffix}")  # base + params

    nameinfo = {
        "cellname_base": nd.cfg.gds_cellname_cleanup(basename),
        "cellname_base_params": nd.cfg.gds_cellname_cleanup(name_long),
        "cellname_full": cellname_full
    }
    return nameinfo


# alternative implementation of the hashme decorator
def cellcache(
    name: str = "",
    parreprs: list = [],
    suffix: str = None,
    addhash: str = None,
    config:dict = None,
) -> Callable:
#def hashme(*outer_args, **outer_kwargs):
    """Function decorator returning a decorator that captures a signature and caches the function results.
    """
    if config is None:
        config = {}
    if name == "":
        name = DEFAULT_CELLNAME
        nd.main_logger(
            f"No cellname given to @cellcache, calling it '{name}'.",
            "warning"
        )
    # get parameter representations in the right format
    if not isinstance(parreprs, list):
        nd.main_logger(
            f"'parreprs' should be a list. Current value: {parreprs}.",
            "error"
        )
    else:
        for i, parrepr in enumerate(parreprs):
            if isinstance(parrepr, str):
                parreprs[i] = {"parameter": parrepr}
            elif isinstance(parrepr, dict):
                parreprs[i] = parrepr
            else:
                nd.main_logger(
                    "Hashme input variable 'parreprs' (parameter representations) should be a list-of-str or list-of-dict. "
                    f"Current value: {parreprs}.",
                    "error"
                )
                parreprs[i] = {"parameter": "ERROR"}

        # Fill all values for parameter representation in the cellname via 'config' or the default config in cfg:
        for parrepr in parreprs:
            for key in ["prefix", "format", "func"]:
                if key not in parrepr.keys():
                    parrepr[key] = config.get(key, nd.cfg.cellcache_config[key])
    if addhash is None:
        addhash = config.get("addhash", nd.cfg.cellcache_config["addhash"])
    if suffix is None:
        suffix = config.get("suffix", nd.cfg.cellcache_config["suffix"])
    assert isinstance(addhash, bool), f"Hashme parameter 'addhash' should be a boolean. Current value: {addhash}."

    def decorator(cellfunc) -> Callable:
        """Decorator that captures the function signature for a cell and caches the result"""
        # Get the parameter names from the function signature
        sig = inspect.signature(cellfunc)

        @functools.wraps(cellfunc)
        def wrapper(*args, **kwargs):
            """Wrapper function that captures the function signature and caches the result

            Make sure the state of the function is fully captured in the signature.
            Add global parameters to the signature if needed.
            """
            nonlocal name
            basename = name

            # Combine args and kwargs with the function's default parameters
            bound_args = sig.bind_partial(*args, **kwargs)
            bound_args.apply_defaults()
            arguments = bound_args.arguments

            # Handle the parameter formatting
            formatted_params = ', '.join(f"{k}={v}" for k, v in arguments.items())

            # Concatenate the function name with the formatted parameters and hash this string
            full_signature = f"{cellfunc.__name__}({formatted_params})"

            #signature_hash = hashlib.sha256(full_signature.encode()).hexdigest()
            signature_hash = hashlib.md5(full_signature.encode()).hexdigest()
            signature_hash_short = signature_hash[:HASH_LENGTH]

            # TODO: Note that two cellnames can become the same due to formatting, while the signature are actually different
            # How to handle this, as a fab may not allow a certain filename like "cell_L0.50003"
            # nor the hash or counter at the end to make the name unique per signature.

                #for p in sorted(arguments.keys()):
                #    hashstr += "{}_{}".format(p, arguments[p])
                #hash = hashlib.md5(hashstr.encode()).hexdigest()[:HASH_LENGTH]

            # check if the cell (signature) is in the memory or sqlite cache. If not, compute it.
            cache_found = False  # covers both memory and sqlite cache
            cache_stale = False
            if signature_hash in nd.cfg.cellnames:
                cache_found = True
                #print(f"using cellcache for call {func.__name__}({formatted_params})")
            elif sqlite:
                # Get the current timestamp of the source file
                source_file = os.path.realpath(cellfunc.__code__.co_filename)
                current_timestamp = int(os.path.getmtime(source_file))
                cursor.execute("SELECT value, timestamp FROM my_table where key = ?", (signature_hash,))
                result = cursor.fetchone()
                if result is not None:
                    cache_found = True
                    cached_value, cached_timestamp = result
                    if int(cached_timestamp) < current_timestamp:
                        # print(f"{int(cached_timestamp)} < {current_timestamp}"
                        cache_stale = True
                    else:
                        # The database cache is not stale, so deserialize the cached value and add to in-memory cache.
                        cell = nd.Cell.json_init(json.loads(cached_value))
                        print(f"using sqlite cache for {signature_hash}, {name}")
                        #nd.cfg.cellnames[signature_hash] = cell
                        nd.cfg.cellnames[signature_hash] = cell

            if not cache_found or cache_stale:
                # Assign bound_args and other cellhash properties to the decorated function
                cellcachemodel = {
                    "funcname": cellfunc.__name__,
                    "func_id": id(cellfunc),
                    "paramdict": bound_args.arguments,
                    "paramnames": list(bound_args.arguments.keys()),
                    "hash": signature_hash,
                    "cellname_hash": signature_hash_short
                }
                nameinfo = _create_cellname(arguments, basename, parreprs, signature_hash_short, addhash=addhash, suffix=suffix)
                cellcachemodel.update(nameinfo)
                stacklen = len(nd.cfg.cellcachemodel)
                nd.cfg.cellcachemodel.append(cellcachemodel)
                cell = cellfunc(*args, **kwargs)
                if stacklen != len(nd.cfg.cellcachemodel):
                    # cellcachemodel should be popped in the function as (and if!) it call the Cell initialization
                    nd.cfg.cellcachemodel = nd.cfg.cellcachemodel[:stacklen]
                    if len(nd.cfg.cellcachemodel) + 1 != stacklen:
                        unused_cache = ", ".join([f'"{name}"' for i in range(stacklen, len(nd.cfg.cellcachemodel) + 1)])
                        nd.main_logger(
                            f"Cellcache stack length mismatch; Found={len(nd.cfg.cellcachemodel) + 1}, expected={stacklen}. "\
                            "Probably no Cell was created within the @cellcache scope. Removing unused cache element(s): "\
                            f"{unused_cache}",
                            "error"
                        )
                if not isinstance(cell, nd.Cell):
                    nd.main_logger(f'Function {cellfunc.__name__} is expected to return a Cell object with name="{name}".', "warning")
                    return None
                nd.cfg.cellnames[signature_hash] = cell
                # nd.cfg.cellnames[] = cell

                if sqlite:
                       # use cell serialization for database storage
                    serialized_value = json.dumps(nd.cell_tree_serialize(cell))
                    if cache_stale:
                        cursor.execute("UPDATE my_table SET value = ?, timestamp = ? WHERE key = ?", (serialized_value, current_timestamp, signature_hash))
                        print(f"updating sqlite cache for {signature_hash}, {name}")
                    else:
                        cursor.execute("INSERT INTO my_table VALUES (?, ?, ?)", (signature_hash, serialized_value, current_timestamp))
                        print(f"adding to sqlite cache for {signature_hash}, {name}")
                    conn.commit()
            return nd.cfg.cellnames[signature_hash]
        return wrapper
    return decorator


if __name__ == '__main__':
    # test the hahme decorator replacement: cellcache.
    # Use the mock cell class MockCell to emulate the Cell class for testing.


    start = time.time()
    @cellcache()
    def my_function(a, b, c=1):
        """Representative state function to test the cellcache decorator"""
        time.sleep(0.01)
        with nd.Cell() as C:
            pass
        return C

    N = 100
    for i in range(N):
        my_function(5, b=i%4)
        #my_function(5, b=6)  # Prints: Calling my_function(a=5, b=6, c=1)
        #print(my_function.bound_args)

    #my_function(5, b=6)  # repeat call
    #print(my_function.bound_args)

    BB = my_function(a=4, b=3, c=2)  # new call

    from pprint import pprint
    pprint(BB.cellcache)

    print(f"Time for {N} instances: {time.time() - start:0.3f}")