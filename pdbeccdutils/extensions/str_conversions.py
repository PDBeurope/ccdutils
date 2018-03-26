#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2018 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.


def str_is_int(i):
    """Check if str can be converted to int

    Args:
        i (str): string representation of an integer

    Returns:
        bool: [description]
    """

    try:
        int(i)
        return True
    except ValueError:
        return False


def str_to_int(i):
    """Converts the string into integer

    Args:
        i (str): string representation of an integer

    Returns:
        int: conversion of the string
    """
    try:
        return int(i)
    except ValueError:
        return 0


def str_to_float(f):
    """Converts the string into float

    Args:
        i (str): string representation of a float

    Returns:
        float: conversion of the string
    """
    try:
        return float(f)
    except ValueError:
        return 0.0
