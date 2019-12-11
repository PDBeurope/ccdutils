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


def str_to_int(i):
    """
    Converts a string into integer. Returns 0 if a string cannot be
    converted.

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
    """
    Converts a string into float. Returns 0.0 if a string cannot be
    converted.

    Args:
        i (str): string representation of a float

    Returns:
        float: conversion of the string
    """
    try:
        return float(f)
    except ValueError:
        return 0.0


def listit(t):
    """Format deep tuples into deep list

    Args:
        t (tuple of tuples): deep tuples

    Returns:
        list[list]: deep list
    """
    return list(map(listit, t)) if isinstance(t, (list, tuple)) else t
