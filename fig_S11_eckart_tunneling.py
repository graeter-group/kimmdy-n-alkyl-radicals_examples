"""Script for plotting eckart tunneling factor for figure S11. Note that rmgpy must be installed for this to work."""
from rmgpy import kinetics, quantity
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import kimmdy_paper_theme

plot_config = kimmdy_paper_theme.default_plot_config
plot_config['axes.titlesize'] = 10
plot_config['xtick.labelsize'] = 10
plot_config['ytick.labelsize'] = 10
plot_config['axes.labelsize'] = 10
plot_config['lines.linewidth'] = 2
pprint(plot_config)
kimmdy_paper_theme.apply_plot_config(plot_config=plot_config)
kimmdy_paper_theme.init_roboto_font()
plot_colors = kimmdy_paper_theme.plot_colors

Ts = np.linspace(273.15, 1000, 100)

# DFT results taken from SI Tableau 3 of

# Tunneling in Hydrogen-Transfer Isomerization of n-Alkyl Radicals
# Baptiste Sirjean, Enoch Dames, Hai Wang, and Wing Tsang
# The Journal of Physical Chemistry A 2012 116 (1), 319-332
# DOI: 10.1021/jp209360u

dft_results = {
    "1-6":{"freq": 1621.0, "barrier_forward": 14.7, "barrier_backward": 17.6, "color": "HITS_YELLOW"},
    "1-5":{"freq": 1678.2, "barrier_forward": 15.0, "barrier_backward": 17.7, "color": "HITS_GREEN"},
    "1-4":{"freq": 1833.4, "barrier_forward": 22.0, "barrier_backward": 24.7, "color": "HITS_CYAN"},
}


fig = plt.figure()
for shift, values in dft_results.items():


    E_reactant = 0  # set to zero  # kcal/mol
    E_TS = values["barrier_forward"]  # kcal/mol
    E_product = (values["barrier_forward"] - values["barrier_backward"])  # kcal/mol

    assert np.isclose(E_TS-E_product, values["barrier_backward"]) and np.isclose(E_TS-E_reactant, values["barrier_forward"])

    rmgpy_E0_reac = quantity.Energy(E_reactant, 'kcal/mol')
    rmgpy_E0_TS = quantity.Energy(E_TS, 'kcal/mol')
    rmgpy_E0_prod = quantity.Energy(E_product, 'kcal/mol')
    rmgpy_frequency = quantity.Frequency(values["freq"], '1/m')

    K = kinetics.Eckart(rmgpy_frequency, rmgpy_E0_reac, rmgpy_E0_TS, rmgpy_E0_prod)

    tunneling_factors = np.array([
        K.calculate_tunneling_factor(temperature_K)
        for temperature_K in Ts
    ])

    pprint(list(zip(Ts, tunneling_factors)))

    plt.plot(Ts, tunneling_factors, color=plot_colors[values["color"]], label=shift)



plt.xlabel('Temperature [K]')
plt.ylabel(r'Multiplicative Factor $\kappa(T)$')
plt.title('Eckart Tunneling in 1-heptyl radical')
plt.yscale("log")
plt.ylim(1.0, plt.ylim()[1])
plt.legend()
plt.tight_layout()
plt.savefig("tunneling_SI.png", bbox_inches='tight')
plt.show()
