#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Client-side script to import APIs.

TODO: put this script behind an api.
Instead of directly calling functions on the server it should send a request.

TODO: this script just gets the APIs, but over time these files can become outdated.
So there should be some automated check on every use that verifies the API is still in sync with the server version.
One way to do this is by storing a version number in the api itself, then have a function similar to check_server_health()
that returns this version: check_api_version(). If it is not the same on the server and client a message will be printed on the client side to hint the user to update.
The update procedure is simple, you just run get_apis() again.
One downside is that we then also need to maintain a version number in the code.
Another way is that check_api_version() hashes the api.py file itself, and sends the hash to the server that then checks if it is up to date.

@author: Mike Machielsen, Ronald Broeke
Bright Photonics (c) 2024
"""

import os
from api_get_api import get_apis

# Placeholder, this will be the path to where the PDK APIs are stored on the client.
CLIENT_DIR = os.path.dirname(__file__)
LOCAL_PDK_DIR = os.path.join(CLIENT_DIR, "api")
TARGET_DIR = "api"

def save_api(api: str, pdk_name: str) -> None:
    """Store API into a client side file.

    Args:
        api (str): python file containig the API.
        target_dir (str, optional): Directory in which to place the new api. Defaults to "api".
        pdk_name (str, optional): Name of the PDK. Defaults to "demo.py".
    """
    try:
        filedir = os.path.dirname(__file__)

        destination_dir = os.path.join(filedir, TARGET_DIR)
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)

        with open(os.path.join(destination_dir, pdk_name), "w") as file:
            file.write(api)

    except Exception as err:
        print(err)
        print(f"ERROR: could not save file '{pdk_name}'")


def get_apis(payload: dict):
    """Download the apis from a given token from Nazca server.

    Download the client API files of the modules that are subscribed to as specified in the token.

    Args:
        token (dict): dictionary with payload information, including PDK subscriptions.
    """
    apis = get_apis(payload=payload)
    for pdk_name, api in apis.items():
        save_api(api, f"{pdk_name}.py")
        print(f"API for {pdk_name} saved.")


if __name__ == "__main__":
    # create fake token information as test
    subscriptions = {
        "pdk_hhi": "",
        "pdk_smart_hs28pc": "",
        "pdk_demofab": "",
        "pdk_ligentec": "",
    }
    payload = {"subscribe": subscriptions}
    get_apis(payload=payload)

    # now test if it can be downloaded

    # import imported.pdk_smart_hs28pc as sp
    # import imported.pdk_hhi as hhi

    # print("The following PDKs are loaded:")
    # sp.get_pdk_name()
    # hhi.get_pdk_name()
