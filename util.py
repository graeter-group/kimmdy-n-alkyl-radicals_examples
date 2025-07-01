"""Various utility functions"""
from collections import defaultdict
import numpy as np
from ase.units import kcal, kJ, kB

from tqdm.autonotebook import tqdm
from kimmdy.recipe import Break, RecipeCollection, Recipe
import pandas as pd

h_planck = 4.135667696e-15  # eV/s
R = 1.987e-3  # Gas constant in kcal/(mol*K)


def get_ids(r):
    return tuple(
        sorted((r.recipe.recipe_steps[1].atom_id_1, r.recipe.recipe_steps[1].atom_id_2))
    )


def round_to_stride(times, stride):
    return np.round(times / stride) * stride


def load_results(rates_csv, se_dir, T, stride=100, start_t=0):
    recipe_collection = RecipeCollection.from_csv(rates_csv)[0]

    recipes_d = defaultdict(list)
    recipes_t_d = defaultdict(list)
    for recipe in tqdm(recipe_collection.recipes):
        recipes_d[recipe.recipe_steps[0]].append(recipe)
        i = int(recipe.recipe_steps[1].atom_id_1)  # idx 1 for break/bind, idx2 for break,move,bind
        j = int(recipe.recipe_steps[1].atom_id_2)
        t = round_to_stride(recipe.timespans[0][0] * 1000, stride)
        recipes_t_d[tuple(sorted((i, j)) + [t])].append(recipe)


    for k in tqdm(recipes_d.keys()):
        reaction = recipes_d[k]
        recipes_d[k] = sorted(reaction, key=lambda r: r.timespans[0][0])
    for k in tqdm(recipes_t_d.keys()):
        reaction = recipes_t_d[k]
        recipes_t_d[k] = sorted(reaction, key=lambda r: r.timespans[0][0])

    df = pd.DataFrame(list(recipes_t_d.values()))
    df = df.rename({0: "recipe", 1: "translation", 2: "hash"}, axis=1)

    df["rates (1/s)"] = df.apply(lambda r: arrhenius_to_eyring_rate(r.recipe.rates[0] * 1e12, T=T), axis=1)  # 1/ps to 1/s

    df["ids"] = df.apply(get_ids, axis=1)
    df["time (fs)"] = df.apply(lambda r: round_to_stride(r.recipe.timespans[0][0] * 1000, stride), axis=1)  # ps to fs
    df["time (s)"] = df["time (fs)"].copy() * 1e-15
    df['barriers (kcal/mol)'] = df['rates (1/s)'].apply(
        lambda k: calculate_activation_energy(k, T=T, method="eyring"))

    return df


def arrhenius_to_eyring_rate(k_arrhenius, T, A=0.288e12, kappa=1.0):
    factor = kB / h_planck
    return (k_arrhenius / A) * (kappa * factor * T)


def calculate_activation_energy(rate, T, kappa=1.0, method="eyring", A=None):
    """
    Calculate the activation energy from the rate using the Eyring equation.
    """
    if method == "eyring":
        return -R * T * np.log(rate*h_planck/(kappa*kB*T))
    elif method == "arr":
        A = 0.288e12 if not A else A
        return -R * T * np.log(rate / A)
    else:
        raise ValueError("method must be eyring|arr")


def calculate_rate(E, T, kappa=1.0, method="eyring", A=None):
    """k_*kB*T/h * exp(-E/(RT))"""
    if method == "eyring":
        factor = kB / h_planck
        return kappa * factor * T * np.exp(-E / (R * T))
    elif method == "arr":
        A = 0.288e12 if not A else A
        return A * np.exp(-E / (R * T))
    else:
        raise ValueError("method must be eyring|arr")


def brier_score_divergence(ps, qs):
    return np.sum((ps-qs)**2)

